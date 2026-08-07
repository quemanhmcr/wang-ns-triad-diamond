import numpy as np

from src.extremal_helicity_symplectic import (
    child_energy_helicity_tensor,
    factorized_tensor,
    parent_wedge,
    transfer_relevant_strain_observability,
)


def test_isosceles_tensor_factorizes():
    for r in (0.55, 0.610904101586766, 0.9):
        assert np.linalg.norm(child_energy_helicity_tensor(r) - factorized_tensor(r)) < 2e-12


def test_common_sl2_parent_deformation_is_exactly_neutral():
    u = np.array([1.0 + 0.2j, -0.3 + 0.4j])
    v = np.array([-0.2 + 0.1j, 0.7 - 0.5j])
    M = np.array([[1.2, 0.4], [0.5, 1.0]])
    M /= np.sqrt(np.linalg.det(M))
    assert abs(parent_wedge(M @ u, M @ v) - parent_wedge(u, v)) < 1e-12


def test_transfer_relevant_strain_half_bound_and_equality_mode():
    S = np.diag([1.0, 0.0, -1.0])
    Q, N = transfer_relevant_strain_observability(S)
    assert abs(Q - 0.5 * N) < 2e-12
    S2 = np.array([[0.3, 0.2, 0.4], [0.2, -0.1, -0.2], [0.4, -0.2, -0.2]])
    Q2, N2 = transfer_relevant_strain_observability(S2)
    assert Q2 >= 0.5 * N2 - 2e-12
