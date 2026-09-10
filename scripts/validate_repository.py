#!/usr/bin/env python3
"""Perform lightweight checks that do not launch LAMMPS or submit jobs."""
from pathlib import Path
import hashlib
import pandas as pd
from pymatgen.core import Structure
from jbfuncs.gbmaker2 import SpheregbBOMaker

root = Path(__file__).resolve().parents[1]
table = pd.read_csv(root / "inputs/gb_orientations.csv")
assert len(table) == 1291, len(table)
assert table["id"].is_unique
assert set(table["gb_type"].dropna()) <= {"tilt", "asymmetric_tilt", "twist", "mixed", "unclassified"}
structure = Structure.from_file(root / "inputs/sorted_optimized_bulk_POSCAR")
assert len(structure) > 0
assert SpheregbBOMaker is not None
for rel in ["inputs/gb_orientations.csv", "inputs/sorted_optimized_bulk_POSCAR", "inputs/bulk_energy_reference.lammpstrj"]:
    path = root / rel
    print(hashlib.sha256(path.read_bytes()).hexdigest(), rel)
print(f"OK: {len(table)} GB entries; bulk atoms={len(structure)}")
