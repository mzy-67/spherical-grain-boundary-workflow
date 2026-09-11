#!/usr/bin/env python3
"""Perform lightweight checks that do not launch LAMMPS or submit jobs."""
from pathlib import Path
import hashlib
import json
import math
import pandas as pd
import yaml
from pymatgen.core import Structure
from jbfuncs.gbmaker2 import SpheregbBOMaker
from jbfuncs.runtime import spherical_slab_intersection_volume

root = Path(__file__).resolve().parents[1]
table = pd.read_csv(root / "inputs/gb_orientations.csv")
assert len(table) == 1291, len(table)
assert table["id"].is_unique
assert set(table["gb_type"].dropna()) <= {"tilt", "asymmetric_tilt", "twist", "mixed", "unclassified"}
structure = Structure.from_file(root / "inputs/sorted_optimized_bulk_POSCAR")
assert len(structure) > 0
assert SpheregbBOMaker is not None

zenodo = json.loads((root / ".zenodo.json").read_text())
assert zenodo["license"] == "mit"

config = yaml.safe_load((root / "configs/workflow.example.yaml").read_text())
simulation = config["simulation"]
volume = spherical_slab_intersection_volume(
    simulation["msd_analysis_radius_A"],
    simulation["msd_slab_half_thickness_A"],
)
manuscript_volume = math.pi * 30.0 * (35.0**2 - 30.0**2 / 12.0)
assert math.isclose(volume, manuscript_volume, rel_tol=1e-15)
assert simulation["production_steps"] * simulation["timestep_ps"] == 50.0
assert simulation["fit_max_ps"] <= 50.0

for rel in ["inputs/gb_orientations.csv", "inputs/sorted_optimized_bulk_POSCAR", "inputs/bulk_energy_reference.lammpstrj"]:
    path = root / rel
    print(hashlib.sha256(path.read_bytes()).hexdigest(), rel)
print(f"OK: {len(table)} GB entries; bulk atoms={len(structure)}")
