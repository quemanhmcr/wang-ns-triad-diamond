import numpy as np
import pytest

from src.event_anchored_role_registration import envelope_registration_residual
from src.smooth_quadratic_carrier_interface import (
    RELINK_OWNER,
    STATUS,
    STRAIN_OWNER,
    coefficient_obstruction_energy_reentry,
    hard_linear_complement_skew_defect,
    positive_smooth_interface_split,
    quadratic_complement,
    quadratic_partition_diagnostics,
    smooth_carrier_energy_identity,
    smooth_quadratic_interface_balance,
    theorem_certificate,
)


def test_I_minus_Q_is_not_the_energy_complement_of_a_smooth_envelope():
    K = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=complex)
    Q = np.diag([0.5, 0.0]).astype(complex)
    u = np.array([1.0, 1.0], dtype=complex)
    out = hard_linear_complement_skew_defect(Q, K, u)
    assert out["linear_complement_skew_defect"] == pytest.approx(-1.0)
    assert out["overlap_defect_formula"] == pytest.approx(-1.0)
    assert out["overlap_formula_residual"] == pytest.approx(0.0)
    assert out["quadratic_complement_skew_residual"] == pytest.approx(0.0)


def test_smooth_angle_pair_is_an_exact_moving_quadratic_energy_partition():
    theta = np.array([0.0, 0.4, 1.1])
    dtheta = np.array([0.3, -0.2, 0.5])
    Q = np.diag(np.cos(theta)).astype(complex)
    R = np.diag(np.sin(theta)).astype(complex)
    dQ = np.diag(-np.sin(theta) * dtheta).astype(complex)
    dR = np.diag(np.cos(theta) * dtheta).astype(complex)
    out = quadratic_partition_diagnostics((Q, R), (dQ, dR))
    assert out["quadratic_partition_residual"] < 1e-12
    assert out["quadratic_partition_derivative_residual"] < 1e-12


def test_smooth_interface_recombines_outer_commutator_with_diagonal_role_work():
    theta = np.array([0.2, 0.7, 1.2])
    dtheta = np.array([0.3, -0.4, 0.1])
    Q = np.diag(np.cos(theta)).astype(complex)
    R = np.diag(np.sin(theta)).astype(complex)
    dQ = np.diag(-np.sin(theta) * dtheta).astype(complex)
    dR = np.diag(np.cos(theta) * dtheta).astype(complex)
    L = np.array(
        [[0.2, 1.0, -0.3], [-0.4, -0.1, 0.8], [0.5, -0.2, 0.7]],
        dtype=complex,
    )
    u = np.array([1.0 + 0.2j, -0.7, 0.4 - 0.3j], dtype=complex)
    out = smooth_quadratic_interface_balance((Q, R), (dQ, dR), L, u)
    assert out["native_outer_recombination_residual"] < 1e-12
    assert out["moving_pair_antisymmetry_residual"] < 1e-12
    assert out["moving_pair_row_sum_residual"] < 1e-12
    assert out["skew_pair_antisymmetry_residual"] < 1e-12
    assert out["relink_pair_antisymmetry_residual"] < 1e-12
    assert out["strain_pair_symmetry_residual"] < 1e-12
    assert out["skew_pair_row_sum_residual"] < 1e-12
    assert out["relink_pair_row_sum_residual"] < 1e-12
    assert out["strain_pair_row_sum_residual"] < 1e-12
    assert np.linalg.norm(
        out["conservative_relink_pair_matrix"]
        + out["conservative_relink_pair_matrix"].T
    ) < 1e-12
    assert np.linalg.norm(
        out["conservative_relink_pair_matrix"].sum(axis=1)
        - out["signed_conservative_relink_work"]
    ) < 1e-12
    assert abs(out["total_conservative_relink_work"]) < 1e-12
    assert abs(out["global_strain_reconstruction_residual"]) < 1e-12


def test_three_role_partition_has_intrinsic_antisymmetric_moving_relink_flux():
    alpha = np.array([0.35, 0.8])
    beta = np.array([0.25, 1.05])
    dalpha = np.array([0.2, -0.3])
    dbeta = np.array([-0.4, 0.15])
    a1 = np.cos(alpha)
    a2 = np.sin(alpha) * np.cos(beta)
    a3 = np.sin(alpha) * np.sin(beta)
    da1 = -np.sin(alpha) * dalpha
    da2 = np.cos(alpha) * dalpha * np.cos(beta) - np.sin(alpha) * np.sin(beta) * dbeta
    da3 = np.cos(alpha) * dalpha * np.sin(beta) + np.sin(alpha) * np.cos(beta) * dbeta
    As = tuple(np.diag(a).astype(complex) for a in (a1, a2, a3))
    dAs = tuple(np.diag(a).astype(complex) for a in (da1, da2, da3))
    L = np.array([[0.3, 1.1], [-0.7, -0.2]], dtype=complex)
    u = np.array([1.0 - 0.4j, -0.6 + 0.2j])
    out = smooth_quadratic_interface_balance(As, dAs, L, u)
    C = out["conservative_relink_pair_matrix"]
    assert np.linalg.norm(C + C.T) < 1e-12
    assert np.linalg.norm(C.sum(axis=1) - out["signed_conservative_relink_work"]) < 1e-12
    assert abs(C[0, 1]) + abs(C[0, 2]) + abs(C[1, 2]) > 1e-4


def test_direct_Q_squared_NS_energy_law_has_exact_low_low_HH_interface_repartition():
    T = np.zeros((2, 2, 2), dtype=complex)
    T[0, 0, 0] = 1.0
    T[1, 0, 1] = 0.4
    T[1, 1, 0] = -0.2
    T[0, 1, 1] = 0.7
    u = np.array([1.2 + 0.1j, -0.8 + 0.3j])
    V = np.array([1.0, 0.0], dtype=complex)
    Q = np.diag([0.0, 0.6]).astype(complex)
    dQ = np.diag([0.0, -0.15]).astype(complex)
    D = np.diag([-1.0, -3.0]).astype(complex)
    out = smooth_carrier_energy_identity(
        tensor=T,
        state=u,
        resolved_state=V,
        smooth_envelope=Q,
        envelope_rate=dQ,
        viscosity_operator=D,
        viscosity=0.2,
    )
    assert out["low_low_work"] == pytest.approx(0.0)
    assert abs(out["direct_energy_identity_residual"]) < 1e-12
    assert abs(out["resolved_repartition_residual"]) < 1e-12
    assert abs(out["outer_to_native_interface_residual"]) < 1e-12
    assert out["carrier_viscous_dissipation"] >= 0.0


def test_actual_nonidempotent_event_envelope_registers_and_completes_quadratically():
    P = np.diag([1.0, 0.0, 0.0]).astype(complex)
    Q = np.diag([1.0, 0.4, 0.2]).astype(complex)
    R = quadratic_complement(Q)
    u = np.array([1.0 + 2.0j, 3.0 - 1.0j, -2.0j])
    phi = np.array([2.0 - 1.0j, 1.0, 4.0j])
    assert np.linalg.norm(Q @ Q - Q) > 0.1
    assert abs(envelope_registration_residual(P, Q, u, phi)) < 1e-12
    out = quadratic_partition_diagnostics((Q, R))
    assert out["quadratic_partition_residual"] < 1e-12


def test_positive_native_interface_work_has_only_relink_or_existing_strain_owners():
    relink = np.array([3.0, -2.0, 0.5, -0.2])
    strain = np.array([-1.0, 4.0, -0.1, -0.3])
    total = relink + strain
    out = positive_smooth_interface_split(total, relink, strain)
    assert out["positive_native_interface_work"] == pytest.approx(4.4)
    assert out["positive_conservative_relink_work"] == pytest.approx(3.5)
    assert out["positive_existing_strain_work"] == pytest.approx(4.0)
    assert set(out["joint_physical_owners"]) == {RELINK_OWNER, STRAIN_OWNER}
    assert out["new_interface_currency_created"] is False
    assert out["primary_selected"] is False


def test_coefficient_obstruction_only_locates_physical_energy_reentry():
    common = dict(
        terminal_coefficient=0.8 + 0.1j,
        terminal_probe_l2=1.0,
        terminal_carrier_energy=1.0,
        initial_carrier_energy=0.1,
        strain_action=0.01,
        signed_interface_atoms=(0.06, 0.04),
        signed_relink_atoms=(0.04, 0.02),
        signed_strain_atoms=(0.02, 0.02),
    )
    small = coefficient_obstruction_energy_reentry(
        coefficient_obstruction_impulse=0.25,
        **common,
    )
    huge = coefficient_obstruction_energy_reentry(
        coefficient_obstruction_impulse=1.0e12j,
        **common,
    )
    assert small["branch"] == "physical_high_high_transfer_generation"
    assert huge["branch"] == small["branch"]
    assert huge["energy_gate"] == small["energy_gate"]
    assert huge["coefficient_impulse_used_as_physical_work"] is False
    assert huge["coefficient_obstruction_is_interval_locator"] is True
    assert huge["causal_weight_source"] == "actual_smooth_carrier_energy_work"
    assert huge["energy_gate"]["physical_hh_work_lower"] >= 8.0 / 15.0 - 1e-14


def test_energy_selected_interface_branch_routes_only_after_actual_work_gate():
    out = coefficient_obstruction_energy_reentry(
        terminal_coefficient=0.7,
        terminal_probe_l2=1.0,
        terminal_carrier_energy=1.0,
        initial_carrier_energy=0.05,
        strain_action=0.01,
        coefficient_obstruction_impulse=1000.0,
        signed_interface_atoms=(0.25,),
        signed_relink_atoms=(0.15,),
        signed_strain_atoms=(0.10,),
    )
    assert out["branch"] == "smooth_interface_physical_work"
    assert out["energy_gate_branch"] == "classified_residual_physical_work"
    assert out["joint_interface_owners"] == (RELINK_OWNER,)
    assert out["coefficient_impulse_used_as_physical_work"] is False


def test_terminal_carrier_energy_must_dominate_registered_coefficient_energy():
    with pytest.raises(ValueError, match="coefficient lower bound"):
        coefficient_obstruction_energy_reentry(
            terminal_coefficient=2.0,
            terminal_probe_l2=1.0,
            terminal_carrier_energy=1.0,
            initial_carrier_energy=0.1,
            strain_action=0.01,
            coefficient_obstruction_impulse=1.0,
            signed_interface_atoms=(0.1,),
            signed_relink_atoms=(0.1,),
            signed_strain_atoms=(0.0,),
        )


def test_certificate_preserves_PDE_ontology_and_global_scope():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "Q^2" in cert["native_object"]
    assert "hard orthogonal P" in cert["hard_smooth_separation"]
    assert "must not be interpreted alone" in cert["outer_recombination"]
    assert "never used as work" in cert["coefficient_reentry"]
    assert "complementary" in cert["relation_to_donor_quotient"]
    assert "does not prove" in cert["scope"]
