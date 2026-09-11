import json
import math
import re
import tomllib
from pathlib import Path

from spherical_gb.runtime import spherical_slab_intersection_volume

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY = "https://github.com/mzy-67/spherical-grain-boundary-workflow"


def test_release_metadata_agrees():
    zenodo = json.loads((ROOT / ".zenodo.json").read_text())
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    citation = (ROOT / "CITATION.cff").read_text()
    cff_version = re.search(r"^version:\s*([^\s]+)$", citation, re.MULTILINE)
    assert zenodo["license"] == "mit"
    assert zenodo["access_right"] == "open"
    assert cff_version is not None
    assert zenodo["version"] == project["version"] == cff_version.group(1)
    assert project["urls"]["Repository"] == REPOSITORY
    assert REPOSITORY in citation


def test_llzo_dataset_counts_remain_documented():
    readme = (ROOT / "README.md").read_text()
    for count in ("1,291", "1,199", "92"):
        assert count in readme


def test_sphere_slab_volume():
    actual = spherical_slab_intersection_volume(35.0, 15.0)
    expected = math.pi * 30.0 * (35.0**2 - 30.0**2 / 12.0)
    assert math.isclose(actual, expected, rel_tol=1e-15)
