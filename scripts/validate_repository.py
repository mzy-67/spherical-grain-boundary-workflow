#!/usr/bin/env python3
"""Perform repository checks without launching LAMMPS or submitting jobs."""
from pathlib import Path
import hashlib
import json
import math
import pandas as pd
import yaml
from pymatgen.core import Structure

from spherical_gb.config import validate_config
from spherical_gb.geometry import rotation_to_z
from spherical_gb.runtime import spherical_slab_intersection_volume

root = Path(__file__).resolve().parents[1]
table = pd.read_csv(root / "inputs/gb_orientations.csv")
assert len(table) == 1291
assert table["id"].is_unique
structure = Structure.from_file(root / "inputs/sorted_optimized_bulk_POSCAR")
assert len(structure) > 0

config = yaml.safe_load((root / "examples/llzo/config.yaml").read_text())
validate_config(config)
assert set(config["material"]["species_order"]) == {site.specie.symbol for site in structure}
assert config["material"]["mobile_species"] == "Li"
assert (rotation_to_z([0, 0, -1]) @ [0, 0, -1]).tolist() == [0.0, 0.0, 1.0]

volume = spherical_slab_intersection_volume(35.0, 15.0)
expected = math.pi * 30.0 * (35.0**2 - 30.0**2 / 12.0)
assert math.isclose(volume, expected, rel_tol=1e-15)
assert json.loads((root / ".zenodo.json").read_text())["license"] == "mit"

for rel in ["inputs/gb_orientations.csv", "inputs/sorted_optimized_bulk_POSCAR", "inputs/bulk_energy_reference.lammpstrj"]:
    path = root / rel
    print(hashlib.sha256(path.read_bytes()).hexdigest(), rel)
print(f"OK: generic schema validated; LLZO example entries={len(table)}; bulk atoms={len(structure)}")
