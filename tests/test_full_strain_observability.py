import numpy as np
from src.full_strain_observability import strain_observability


def test_basic_tracefree_strains_are_observed():
    strains = [
        np.diag([1.0, -1.0, 0.0]),
        np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]),
        np.diag([1.0, 1.0, -2.0]),
    ]
    for S in strains:
        Q, N = strain_observability(S)
        assert Q >= 13.0 / 20.0 * N - 1e-12
