import numpy as np
from jbfuncs.gbmaker2 import rotation_to_z


def test_rotation_to_z_is_orthogonal():
    matrix = rotation_to_z([1.0, 1.0, 1.0])
    assert np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-10)
