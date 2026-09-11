"""Configuration helpers for the material-agnostic workflow."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping

DEFAULT_SIMULATION = {
    "lammps_executable": "lmp",
    "mobile_radius_A": 39.0,
    "msd_analysis_radius_A": 35.0,
    "msd_slab_half_thickness_A": 15.0,
    "timestep_ps": 0.001,
    "equilibration_steps": 4000,
    "production_steps": 50000,
    "trajectory_dump_interval": 500,
    "fit_min_ps": 10.0,
    "fit_max_ps": 50.0,
    "temperatures_K": [973.15, 1073.15, 1173.15, 1273.15],
}


def merged_simulation(config: Mapping[str, Any]) -> dict[str, Any]:
    result = deepcopy(DEFAULT_SIMULATION)
    result.update(config.get("simulation", {}))
    return result


def validate_config(config: Mapping[str, Any]) -> None:
    """Validate portable fields before any cluster job is created."""
    for section in ("inputs", "potential", "material", "grain_boundary", "optimization", "submission"):
        if section not in config:
            raise ValueError(f"Missing configuration section: {section}")

    material = config["material"]
    transport_enabled = bool(config.get("transport", {}).get("enabled", True))
    order = list(material.get("species_order", []))
    if not order or len(order) != len(set(order)):
        raise ValueError("material.species_order must be a non-empty list of unique symbols")
    if transport_enabled and material.get("mobile_species") not in order:
        raise ValueError("material.mobile_species must occur in material.species_order when transport is enabled")
    if transport_enabled and float(material.get("charge_number", 0)) == 0:
        raise ValueError("material.charge_number must be non-zero when transport is enabled")
    masses = material.get("masses", {})
    bad_masses = [key for key, value in masses.items() if float(value) <= 0]
    if bad_masses:
        raise ValueError(f"Atomic masses must be positive: {bad_masses}")

    gb = config["grain_boundary"]
    for key in ("sphere_radius", "vacuum_thickness", "energy_radius"):
        if float(gb[key]) <= 0:
            raise ValueError(f"grain_boundary.{key} must be positive")
    if float(gb["energy_radius"]) > float(gb["sphere_radius"]):
        raise ValueError("grain_boundary.energy_radius cannot exceed sphere_radius")

    sim = merged_simulation(config)
    positive = (
        "mobile_radius_A", "msd_analysis_radius_A", "msd_slab_half_thickness_A",
        "timestep_ps", "equilibration_steps", "production_steps", "trajectory_dump_interval",
    )
    invalid = [key for key in positive if float(sim[key]) <= 0]
    if invalid:
        raise ValueError(f"Simulation values must be positive: {invalid}")
    if float(sim["mobile_radius_A"]) < float(sim["msd_analysis_radius_A"]):
        raise ValueError("simulation.mobile_radius_A must enclose msd_analysis_radius_A")
    if not 0 <= float(sim["fit_min_ps"]) < float(sim["fit_max_ps"]):
        raise ValueError("simulation fit interval must satisfy 0 <= fit_min_ps < fit_max_ps")
    if float(sim["fit_max_ps"]) > float(sim["production_steps"]) * float(sim["timestep_ps"]):
        raise ValueError("simulation.fit_max_ps exceeds the production trajectory length")
    temperatures = [float(value) for value in sim["temperatures_K"]]
    if len(temperatures) < 2 or any(value <= 0 for value in temperatures):
        raise ValueError("simulation.temperatures_K must contain at least two positive temperatures")
