"""Small geometry helpers with no heavy simulation dependencies."""

from __future__ import annotations

import numpy as np


def rotation_to_z(vector) -> np.ndarray:
    """Return a proper rotation matrix that maps ``vector`` onto +z."""
    vector = np.asarray(vector, dtype=float)
    if vector.shape != (3,):
        raise ValueError("vector must contain exactly three components")
    norm = np.linalg.norm(vector)
    if norm < 1e-12:
        raise ValueError("vector cannot be zero")
    source = vector / norm
    target = np.array([0.0, 0.0, 1.0])
    cosine = float(np.clip(np.dot(source, target), -1.0, 1.0))
    if np.isclose(cosine, 1.0):
        return np.eye(3)
    if np.isclose(cosine, -1.0):
        return np.diag([1.0, -1.0, -1.0])
    axis = np.cross(source, target)
    sine = np.linalg.norm(axis)
    skew = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / sine**2)


def inclusive_radius_grid(max_radius: float, min_radius: float = 3.0, step: float = 1.0) -> np.ndarray:
    """Return an increasing radial grid that always contains ``max_radius``."""
    max_radius, min_radius, step = map(float, (max_radius, min_radius, step))
    if min_radius <= 0 or max_radius < min_radius or step <= 0:
        raise ValueError("Require 0 < min_radius <= max_radius and step > 0")
    radii = np.arange(min_radius, max_radius + 1e-12, step)
    if not np.isclose(radii[-1], max_radius):
        radii = np.append(radii, max_radius)
    return radii
