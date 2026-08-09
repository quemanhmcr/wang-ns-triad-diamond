import numpy as np

from src.nonaffine_role_interface_work import (
    binary_interface_decomposition,
    interface_work,
    partition_interface_balance,
)


def test_skew_interface_is_conservative_between_binary_roles():
    K = np.array([[0.0, 2.0], [-2.0, 0.0]], dtype=complex)
    Q = np.diag([1.0, 0.0]).astype(complex)
    u = np.array([1.0 + 0.5j, -0.7 + 0.2j])
    out = binary_interface_decomposition(Q, K, u)
    assert abs(float(out["skew_conservation_residual"])) < 1e-12
    assert abs(float(out["Q_symmetric_work"])) < 1e-12


def test_symmetric_interface_is_off_diagonal_strain_work():
    S = np.array([[1.0, 0.4], [0.4, -0.5]], dtype=complex)
    Q = np.diag([1.0, 0.0]).astype(complex)
    u = np.array([1.2, -0.8], dtype=complex)
    out = binary_interface_decomposition(Q, S, u)
    assert abs(float(out["Q_symmetric_formula_residual"])) < 1e-12
    assert abs(float(out["symmetric_pair_equality_residual"])) < 1e-12
    assert abs(float(out["Q_skew_work"])) < 1e-12


def test_complete_partition_reconstructs_global_strain_balance():
    Ps = [
        np.diag([1.0, 0.0, 0.0]).astype(complex),
        np.diag([0.0, 1.0, 0.0]).astype(complex),
        np.diag([0.0, 0.0, 1.0]).astype(complex),
    ]
    L = np.array(
        [[0.2, 1.0, -0.3], [-0.4, -0.1, 0.8], [0.5, -0.2, 0.7]],
        dtype=complex,
    )
    u = np.array([1.0, -0.7, 0.4], dtype=complex)
    out = partition_interface_balance(Ps, L, u)
    assert abs(float(out["total_skew_interface"])) < 1e-12
    assert abs(float(out["symmetric_global_balance_residual"])) < 1e-12
    assert float(out["skew_antisymmetry_residual"]) < 1e-12
    assert float(out["symmetric_symmetry_residual"]) < 1e-12


def test_interface_work_vanishes_when_operator_preserves_role():
    Q = np.diag([1.0, 0.0]).astype(complex)
    L = np.diag([2.0, -1.0]).astype(complex)
    u = np.array([1.0, 3.0], dtype=complex)
    assert abs(interface_work(Q, L, u)) < 1e-12
