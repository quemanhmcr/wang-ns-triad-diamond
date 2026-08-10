import numpy as np
import pytest

from src.event_anchored_role_registration import envelope_registration_residual
from src.smooth_quadratic_carrier_interface import (
    GaugeQuotientedInterfaceWork,
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


def _physical_two_role_balance():
    Q = np.diag([1.0, 0.0]).astype(complex)
    R = np.diag([0.0, 1.0]).astype(complex)
    dQ = np.zeros((2, 2), dtype=complex)
    dR = np.zeros((2, 2), dtype=complex)
    G = np.zeros((2, 2), dtype=complex)
    Kphys = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=complex)
    S = np.zeros((2, 2), dtype=complex)
    u = np.array([1.0, 1.0], dtype=complex)
    return smooth_quadratic_interface_balance(
        (Q, R),
        (dQ, dR),
        Kphys + S,
        u,
        common_gauge_skew_operator=G,
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


def test_square_partition_is_exact_but_not_yet_a_physical_owner_law():
    theta = np.array([0.0, 0.4, 1.1])
    dtheta = np.array([0.3, -0.2, 0.5])
    Q = np.diag(np.cos(theta)).astype(complex)
    R = np.diag(np.sin(theta)).astype(complex)
    dQ = np.diag(-np.sin(theta) * dtheta).astype(complex)
    dR = np.diag(np.cos(theta) * dtheta).astype(complex)
    out = quadratic_partition_diagnostics((Q, R), (dQ, dR))
    assert out["quadratic_partition_residual"] < 1e-12
    assert out["quadratic_partition_derivative_residual"] < 1e-12
    with pytest.raises(ValueError, match="quotient observer motion"):
        smooth_quadratic_interface_balance(
            (Q, R),
            (dQ, dR),
            np.zeros((3, 3), dtype=complex),
            np.array([1.0, -0.4, 0.7], dtype=complex),
            common_gauge_skew_operator=np.zeros((3, 3), dtype=complex),
        )


def test_common_transport_gauge_cancels_before_physical_relink():
    theta = np.array([0.2, 0.7, 1.2])
    Q = np.diag(np.cos(theta)).astype(complex)
    R = np.diag(np.sin(theta)).astype(complex)
    G = np.array(
        [[0.0, -0.3, 0.2], [0.3, 0.0, -0.1], [-0.2, 0.1, 0.0]],
        dtype=complex,
    )
    dQ = -(G @ Q - Q @ G)
    dR = -(G @ R - R @ G)
    Kphys = np.array(
        [[0.0, -0.4, 0.1], [0.4, 0.0, -0.2], [-0.1, 0.2, 0.0]],
        dtype=complex,
    )
    S = np.array(
        [[0.2, 0.1, -0.05], [0.1, -0.1, 0.08], [-0.05, 0.08, 0.3]],
        dtype=complex,
    )
    u = np.array([1.0 + 0.2j, -0.7, 0.4 - 0.3j], dtype=complex)
    out = smooth_quadratic_interface_balance(
        (Q, R),
        (dQ, dR),
        G + Kphys + S,
        u,
        common_gauge_skew_operator=G,
    )
    assert out["gauge_transport_operator_residual"] < 1e-12
    assert out["gauge_work_cancellation_residual"] < 1e-12
    assert out["native_gauge_quotient_residual"] < 1e-12
    assert out["native_outer_recombination_residual"] < 1e-12
    assert out["relink_pair_antisymmetry_residual"] < 1e-12
    assert out["strain_pair_symmetry_residual"] < 1e-12
    assert out["relink_pair_row_sum_residual"] < 1e-12
    assert out["strain_pair_row_sum_residual"] < 1e-12
    assert abs(out["total_physical_relink_work"]) < 1e-12
    assert out["observer_partition_motion_charged_as_physics"] is False
    assert isinstance(out["work_certificate"], GaugeQuotientedInterfaceWork)


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


def test_positive_interface_split_requires_gauge_quotiented_certificate():
    work = GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=(2.0, 2.0, 0.4, -0.5),
        signed_physical_relink_atoms=(3.0, -2.0, 0.5, -0.2),
        signed_existing_strain_atoms=(-1.0, 4.0, -0.1, -0.3),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
    )
    out = positive_smooth_interface_split(work)
    assert out["positive_native_interface_work"] == pytest.approx(4.4)
    assert out["positive_conservative_relink_work"] == pytest.approx(3.5)
    assert out["positive_existing_strain_work"] == pytest.approx(4.0)
    assert set(out["joint_physical_owners"]) == {RELINK_OWNER, STRAIN_OWNER}
    assert out["new_interface_currency_created"] is False
    assert out["observer_partition_motion_charged_as_physics"] is False


def test_coefficient_obstruction_magnitude_never_changes_physical_reentry():
    work = _physical_two_role_balance()["work_certificate"]
    common = dict(
        terminal_coefficient=1.0,
        terminal_probe_l2=1.0,
        terminal_carrier_energy=20.0,
        initial_carrier_energy=1.0,
        strain_action=0.01,
        interface_work=work,
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
    assert huge["observer_partition_motion_charged_as_physics"] is False
    assert huge["causal_weight_source"] == "actual_gauge_quotiented_smooth_carrier_energy_work"


def test_energy_selected_interface_branch_uses_actual_gauge_quotiented_work():
    work = _physical_two_role_balance()["work_certificate"]
    out = coefficient_obstruction_energy_reentry(
        terminal_coefficient=1.0,
        terminal_probe_l2=1.0,
        terminal_carrier_energy=5.0,
        initial_carrier_energy=0.1,
        strain_action=0.01,
        coefficient_obstruction_impulse=1000.0,
        interface_work=work,
    )
    assert out["branch"] == "smooth_interface_physical_work"
    assert out["energy_gate_branch"] == "classified_residual_physical_work"
    assert out["joint_interface_owners"] == (RELINK_OWNER,)
    assert out["coefficient_impulse_used_as_physical_work"] is False


def test_terminal_carrier_energy_must_dominate_registered_coefficient_energy():
    work = _physical_two_role_balance()["work_certificate"]
    with pytest.raises(ValueError, match="coefficient lower bound"):
        coefficient_obstruction_energy_reentry(
            terminal_coefficient=2.0,
            terminal_probe_l2=1.0,
            terminal_carrier_energy=1.0,
            initial_carrier_energy=0.1,
            strain_action=0.01,
            coefficient_obstruction_impulse=1.0,
            interface_work=work,
        )


def test_certificate_preserves_PDE_ontology_and_global_scope():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "Q^2" in cert["native_object"]
    assert "hard orthogonal P" in cert["hard_smooth_separation"]
    assert "must not be interpreted alone" in cert["outer_recombination"]
    assert "dot A_a+[G,A_a]=0" in cert["observer_gauge_quotient"]
    assert "cannot be Hahn-routed" in cert["forbidden_observer_motion"]
    assert "never used as work" in cert["coefficient_reentry"]
    assert "separately typed measures" in cert["relation_to_donor_quotient"]
    assert "does not prove" in cert["scope"]
