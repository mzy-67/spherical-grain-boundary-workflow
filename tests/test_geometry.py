import numpy as np
import pytest

from spherical_gb.geometry import rotation_to_z


@pytest.mark.parametrize("vector", ([1.0, 1.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [2.0, -3.0, 4.0]))
def test_rotation_to_z_aligns_vector(vector):
    matrix = rotation_to_z(vector)
    unit = np.asarray(vector, dtype=float) / np.linalg.norm(vector)
    assert np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-10)
    assert np.isclose(np.linalg.det(matrix), 1.0)
    assert np.allclose(matrix @ unit, [0.0, 0.0, 1.0], atol=1e-10)


def test_rotation_to_z_rejects_zero_vector():
    with pytest.raises(ValueError):
        rotation_to_z([0.0, 0.0, 0.0])

from spherical_gb.geometry import inclusive_radius_grid


def test_radius_grid_includes_noninteger_endpoint():
    assert np.allclose(inclusive_radius_grid(5.5), [3.0, 4.0, 5.0, 5.5])
