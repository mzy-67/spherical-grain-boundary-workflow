#!/usr/bin/env python3
"""Validate and submit spherical grain-boundary workflows from a YAML config."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from jobflow_remote import submit_flow
from pymatgen.core import Structure
from qtoolkit.core.data_objects import QResources

from spherical_gb.config import merged_simulation, validate_config
from spherical_gb.workflow import SphericalGBWorkflowMaker, species_symbol

REQUIRED_COLUMNS = {
    "id", "sigma", "misorientation_angle_deg",
    "axis_h", "axis_k", "axis_l",
    "plane_h", "plane_k", "plane_l", "gb_type",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--start-row", type=int, default=0, help="Zero-based inclusive row index")
    parser.add_argument("--stop-row", type=int, default=None, help="Zero-based exclusive row index")
    parser.add_argument("--dry-run", action="store_true", help="Validate and build flows without submission")
    return parser.parse_args()


def resolve(repo_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (repo_root / path).resolve()


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    repo_root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load(config_path.read_text())
    validate_config(cfg)

    orientation_file = resolve(repo_root, cfg["inputs"]["orientation_file"])
    structure_file = resolve(repo_root, cfg["inputs"]["structure_file"])
    bulk_energy_file = resolve(repo_root, cfg["inputs"]["bulk_energy_trajectory"])
    potential = cfg["potential"]
    checkpoint_value = potential.get("checkpoint_file")
    checkpoint_file = resolve(repo_root, checkpoint_value) if checkpoint_value else None

    for required in (orientation_file, structure_file, bulk_energy_file):
        if not required.is_file():
            raise FileNotFoundError(required)
    if checkpoint_file is not None and not checkpoint_file.is_file():
        raise FileNotFoundError(checkpoint_file)

    table = pd.read_csv(orientation_file)
    missing = REQUIRED_COLUMNS.difference(table.columns)
    if missing:
        raise ValueError(f"Missing orientation columns: {sorted(missing)}")
    stop = len(table) if args.stop_row is None else args.stop_row
    if not 0 <= args.start_row < stop <= len(table):
        raise ValueError(f"Invalid row interval [{args.start_row}, {stop}) for {len(table)} rows")

    structure = Structure.from_file(structure_file)
    material = cfg["material"]
    transport = cfg.get("transport", {"enabled": True})
    annealing = cfg.get("annealing", {})
    configured_species = tuple(str(value) for value in material["species_order"])
    structure_species = {species_symbol(sp) for sp in structure.species}
    if structure_species != set(configured_species):
        raise ValueError(
            "Configured species do not match the input structure: "
            f"structure={sorted(structure_species)}, config={list(configured_species)}"
        )

    opt = cfg["optimization"]
    gb = cfg["grain_boundary"]
    sim = merged_simulation(cfg)
    sub = cfg["submission"]
    resources = QResources(
        nodes=int(sub["nodes"]),
        processes_per_node=int(sub["processes_per_node"]),
        gpus_per_job=int(sub["gpus_per_job"]),
        scheduler_kwargs={
            "partition": sub["partition"],
            "qverbatim": "\n".join([
                f"#SBATCH --cpus-per-gpu={int(sub['cpus_per_gpu'])}",
                f"#SBATCH --mem={sub['memory']}",
            ]),
        },
    )

    for row_index in range(args.start_row, stop):
        row = table.iloc[row_index]
        metadata = (
            f"gb{int(row.id):04d}_sigma{row.sigma}_"
            f"axis{row.axis_h:g}-{row.axis_k:g}-{row.axis_l:g}_"
            f"ang{row.misorientation_angle_deg:.4f}_"
            f"plane{row.plane_h:g}-{row.plane_k:g}-{row.plane_l:g}_"
            f"{row.gb_type}"
        )
        maker = SphericalGBWorkflowMaker(
            name=str(cfg.get("workflow", {}).get("name", "spherical-grain-boundary")),
            trials=int(opt["trials"]),
            base_estimator=str(opt.get("base_estimator", "GP")),
            acq_func=str(opt.get("acq_func", "EI")),
            acq_optimizer=str(opt.get("acq_optimizer", "lbfgs")),
            random_state=int(opt["random_seed"]),
            crystal_structure=structure,
            check_point_file=str(checkpoint_file) if checkpoint_file else None,
            bulk_energy_traj_file=str(bulk_energy_file),
            species_order=configured_species,
            species_masses={str(k): float(v) for k, v in material.get("masses", {}).items()},
            mobile_species=(str(material["mobile_species"]) if material.get("mobile_species") else None),
            charge_number=float(material.get("charge_number", 1.0)),
            run_transport=bool(transport.get("enabled", True)),
            anneal_temperature_K=float(annealing.get("temperature_K", 1200.0)),
            anneal_equilibration_steps=int(annealing.get("equilibration_steps", 1000)),
            anneal_quench_steps=int(annealing.get("quench_steps", 20000)),
            origin_anchor_species=material.get("origin_anchor_species"),
            pair_style=str(potential["pair_style"]),
            pair_style_args=str(potential.get("pair_style_args", "")),
            pair_coeff=str(potential.get("pair_coeff", "* *")),
            sphere_R=float(gb["sphere_radius"]),
            vaccum_thickness=float(gb["vacuum_thickness"]),
            gb_r=float(gb["energy_radius"]),
            rot_axis=[float(row.axis_h), float(row.axis_k), float(row.axis_l)],
            rot_angle=float(np.deg2rad(row.misorientation_angle_deg)),
            normal=[float(row.plane_h), float(row.plane_k), float(row.plane_l)],
            metadata=metadata,
            lammps_executable=str(sim["lammps_executable"]),
            mobile_radius_A=float(sim["mobile_radius_A"]),
            msd_analysis_radius_A=float(sim["msd_analysis_radius_A"]),
            msd_slab_half_thickness_A=float(sim["msd_slab_half_thickness_A"]),
            msd_timestep_ps=float(sim["timestep_ps"]),
            msd_equilibration_steps=int(sim["equilibration_steps"]),
            msd_production_steps=int(sim["production_steps"]),
            msd_dump_interval=int(sim["trajectory_dump_interval"]),
            msd_fit_min_ps=float(sim["fit_min_ps"]),
            msd_fit_max_ps=float(sim["fit_max_ps"]),
            msd_temperatures_K=tuple(float(value) for value in sim["temperatures_K"]),
        )
        flow = maker.make()
        if args.dry_run:
            print(f"validated row={row_index}, id={int(row.id)}, metadata={metadata}")
        else:
            submit_flow(flow, worker=sub["worker"], resources=resources, project=sub["project"])
            print(f"submitted row={row_index}, id={int(row.id)}, metadata={metadata}")


if __name__ == "__main__":
    main()
