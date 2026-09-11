from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = (ROOT / "src/spherical_gb/workflow.py").read_text()


def test_core_has_no_fixed_llzo_site_count_or_element_dump_order():
    assert "range(192)" not in SOURCE
    assert "element Li La Zr O" not in SOURCE
    assert "type_map = {\"Li\"" not in SOURCE


def test_generic_package_and_llzo_example_are_both_present():
    assert (ROOT / "src/spherical_gb/__init__.py").is_file()
    assert (ROOT / "examples/llzo/config.yaml").is_file()
