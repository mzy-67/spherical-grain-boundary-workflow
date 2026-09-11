"""Small runtime helpers that do not depend on the scientific Python stack."""

from __future__ import annotations

import math
import shlex
import subprocess


def spherical_slab_intersection_volume(radius: float, half_thickness: float) -> float:
    """Return the volume in A^3 of a sphere intersected by a centered slab.

    ``half_thickness`` is the distance from the slab mid-plane to either face.
    When the slab is thicker than the sphere, the full sphere volume is used.
    """
    radius = float(radius)
    half_thickness = float(half_thickness)
    if radius <= 0:
        raise ValueError("radius must be positive")
    if half_thickness <= 0:
        raise ValueError("half_thickness must be positive")
    h = min(radius, half_thickness)
    return 2.0 * math.pi * (radius**2 * h - h**3 / 3.0)


def run_lammps(executable: str, input_file: str, log_file: str) -> None:
    """Run LAMMPS without a shell and fail immediately on a non-zero exit."""
    executable_parts = shlex.split(executable)
    if not executable_parts:
        raise ValueError("LAMMPS executable must not be empty")
    command = [*executable_parts, "-i", input_file, "-log", log_file]
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if completed.returncode != 0:
        output_tail = "\n".join((completed.stdout or "").splitlines()[-40:])
        raise RuntimeError(
            f"LAMMPS failed with exit code {completed.returncode}: "
            f"{shlex.join(command)}\nLast output lines:\n{output_tail}"
        )
