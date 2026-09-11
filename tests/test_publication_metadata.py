import json
import math
import re
import tomllib
from pathlib import Path

from jbfuncs.runtime import spherical_slab_intersection_volume


ROOT = Path(__file__).resolve().parents[1]


def test_zenodo_metadata_uses_recognized_mit_id():
    metadata = json.loads((ROOT / ".zenodo.json").read_text())
    assert metadata["license"] == "mit"
    assert metadata["version"] == "0.1.0"
    assert metadata["access_right"] == "open"


def test_release_versions_agree():
    zenodo = json.loads((ROOT / ".zenodo.json").read_text())
    project = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
    citation = (ROOT / "CITATION.cff").read_text()
    cff_version = re.search(r"^version:\s*([^\s]+)$", citation, re.MULTILINE)
    assert cff_version is not None
    assert zenodo["version"] == project["version"] == cff_version.group(1)


def test_citation_points_to_the_repository_without_placeholders():
    citation = (ROOT / "CITATION.cff").read_text()
    assert "https://github.com/mzy-67/LLZO-spherical-GB" in citation
    assert "OWNER" not in citation


def test_dataset_counts_are_explicitly_documented():
    readme = (ROOT / "README.md").read_text()
    for count in ("1,291", "1,230", "1,199"):
        assert count in readme


def test_manuscript_sphere_slab_volume():
    actual = spherical_slab_intersection_volume(35.0, 15.0)
    reported_formula = math.pi * 30.0 * (35.0**2 - 30.0**2 / 12.0)
    assert math.isclose(actual, reported_formula, rel_tol=1e-15)


def test_sphere_slab_volume_rejects_nonpositive_dimensions():
    for radius, half_thickness in ((0, 15), (35, 0), (-1, 15)):
        try:
            spherical_slab_intersection_volume(radius, half_thickness)
        except ValueError:
            continue
        raise AssertionError(f"accepted invalid dimensions: {radius}, {half_thickness}")
