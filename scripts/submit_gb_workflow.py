#!/usr/bin/env python3
"""Submit a selected row range from the LLZO GB orientation table."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from jobflow_remote import submit_flow
from pymatgen.core import Structure
from qtoolkit.core.data_objects import QResources

from jbfuncs.gbmaker2 import SpheregbBOMaker

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
    parser.add_argument("--dry-run", action="store_true", help="Validate and print jobs without submission")
    return parser.parse_args()


def resolve(repo_root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (repo_root / path).resolve()


def main() -> None:
    args = parse_args()
    config_path = args.config.expanduser().resolve()
    repo_root = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load(config_path.read_text())

    orientation_file = resolve(repo_root, cfg["inputs"]["orientation_file"])
    structure_file = resolve(repo_root, cfg["inputs"]["structure_file"])
    bulk_energy_file = resolve(repo_root, cfg["inputs"]["bulk_energy_trajectory"])
    checkpoint_file = resolve(repo_root, cfg["potential"]["checkpoint_file"])

    for required in (orientation_file, structure_file, bulk_energy_file, checkpoint_file):
        if not required.is_file():
            raise FileNotFoundError(required)

    table = pd.read_csv(orientation_file)
    missing = REQUIRED_COLUMNS.difference(table.columns)
    if missing:
        raise ValueError(f"Missing orientation columns: {sorted(missing)}")

    stop = len(table) if args.stop_row is None else args.stop_row
    if not 0 <= args.start_row < stop <= len(table):
        raise ValueError(f"Invalid row interval [{args.start_row}, {stop}) for {len(table)} rows")

    llzo = Structure.from_file(structure_file)
    opt = cfg["optimization"]
    gb = cfg["grain_boundary"]
    sub = cfg["submission"]

    qverbatim = "\n".join([
        f"#SBATCH --cpus-per-gpu={int(sub['cpus_per_gpu'])}",
        f"#SBATCH --mem={sub['memory']}",
    ])
    resources = QResources(
        nodes=int(sub["nodes"]),
        processes_per_node=int(sub["processes_per_node"]),
        gpus_per_job=int(sub["gpus_per_job"]),
        scheduler_kwargs={"partition": sub["partition"], "qverbatim": qverbatim},
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
        maker = SpheregbBOMaker(
            name="llzo-spherical-gb",
            trials=int(opt["trials"]),
            random_state=int(opt["random_seed"]),
            crystal_structure=llzo,
            check_point_file=str(checkpoint_file),
            bulk_energy_traj_file=str(bulk_energy_file),
            sphere_R=float(gb["sphere_radius"]),
            vaccum_thickness=float(gb["vacuum_thickness"]),
            gb_r=float(gb["energy_radius"]),
            rot_axis=[float(row.axis_h), float(row.axis_k), float(row.axis_l)],
            rot_angle=float(np.deg2rad(row.misorientation_angle_deg)),
            normal=[float(row.plane_h), float(row.plane_k), float(row.plane_l)],
            metadata=metadata,
        )
        flow = maker.make()
        if args.dry_run:
            print(f"validated row={row_index}, id={int(row.id)}, metadata={metadata}")
            continue
        submit_flow(flow, worker=sub["worker"], resources=resources, project=sub["project"])
        print(f"submitted row={row_index}, id={int(row.id)}, metadata={metadata}")


if __name__ == "__main__":
    main()
