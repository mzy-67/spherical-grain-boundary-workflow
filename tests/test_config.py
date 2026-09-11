import pytest

from spherical_gb.config import merged_simulation, validate_config


def config():
    return {
        "inputs": {},
        "potential": {},
        "material": {"species_order": ["Na", "Cl"], "mobile_species": "Na", "charge_number": 1},
        "grain_boundary": {"sphere_radius": 30, "vacuum_thickness": 10, "energy_radius": 25},
        "optimization": {},
        "simulation": {"production_steps": 60000},
        "submission": {},
    }


def test_generic_binary_material_config_is_valid():
    cfg = config()
    validate_config(cfg)
    assert merged_simulation(cfg)["production_steps"] == 60000


def test_mobile_species_must_be_configured():
    cfg = config()
    cfg["material"]["mobile_species"] = "Li"
    with pytest.raises(ValueError, match="mobile_species"):
        validate_config(cfg)
