import numpy as np
import pytest

from src.resolved_interface_donor_quotient import (
    SKEW_OWNER,
    STATUS,
    SYMMETRIC_OWNER,
    high_tail_interface_component_route,
    positive_interface_component_split,
    resolved_role_work_decomposition,
    skew_donor_closure,
    skew_flux_balance,
    skew_subset_balance,
    theorem_certificate,
)



def test_full_resolved_low_high_work_splits_into_conservative_K_flux_and_existing_S_strain():
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
    out = resolved_role_work_decomposition(Ps, L, u)
    assert out["work_split_residual"] < 1e-12
    assert out["skew_antisymmetry_residual"] < 1e-12
    assert out["symmetric_symmetry_residual"] < 1e-12
    assert out["total_skew_work"] == pytest.approx(0.0)
    assert out["symmetric_is_existing_strain_provenance"] is True
    assert out["new_interface_source_created"] is False

def test_positive_interface_is_covered_by_skew_or_existing_strain_without_new_currency():
    skew = np.array([3.0, -2.0, 0.5, -0.2])
    sym = np.array([-1.0, 4.0, -0.1, -0.3])
    total = skew + sym
    out = positive_interface_component_split(total, skew, sym)
    assert out["positive_interface_work"] == pytest.approx(4.4)
    assert out["positive_skew_redistribution_work"] == pytest.approx(3.5)
    assert out["positive_symmetric_strain_work"] == pytest.approx(4.0)
    assert out["positive_cover_margin"] == pytest.approx(3.1)
    assert set(out["joint_owners"]) == {SKEW_OWNER, SYMMETRIC_OWNER}
    assert out["new_source_currency_created"] is False
    assert out["primary_selected"] is False


def test_skew_pair_work_is_exact_directed_flux_divergence():
    T = np.array(
        [
            [0.0, 3.0, -1.0],
            [-3.0, 0.0, 2.0],
            [1.0, -2.0, 0.0],
        ]
    )
    out = skew_flux_balance(T)
    assert out["net_skew_role_work"] == pytest.approx(np.array([2.0, -1.0, -1.0]))
    assert out["incoming_positive_flux"] == pytest.approx(np.array([3.0, 2.0, 1.0]))
    assert out["outgoing_positive_flux"] == pytest.approx(np.array([1.0, 3.0, 2.0]))
    assert out["worst_role_divergence_residual"] == pytest.approx(0.0)
    assert out["total_net_skew_work"] == pytest.approx(0.0)


def test_internal_skew_circulation_cancels_from_every_subset_balance():
    T = np.array(
        [
            [0.0, 4.0, -2.0, 0.0],
            [-4.0, 0.0, 3.0, 1.0],
            [2.0, -3.0, 0.0, -1.0],
            [0.0, -1.0, 1.0, 0.0],
        ]
    )
    out = skew_subset_balance(T, (0, 1, 2))
    assert out["subset_net_skew_work"] == pytest.approx(out["boundary_signed_flux"])
    assert out["subset_net_skew_work"] == pytest.approx(
        out["boundary_positive_inflow"] - out["boundary_positive_outflow"]
    )
    assert out["internal_circulation_contribution"] == 0.0


def test_positive_skew_gain_has_finite_same_event_negative_net_donor_set():
    T = np.array(
        [
            [0.0, 3.0, -1.0],
            [-3.0, 0.0, 2.0],
            [1.0, -2.0, 0.0],
        ]
    )
    out = skew_donor_closure(T, (0,))
    assert out["recipient_net_gain"] == pytest.approx(2.0)
    assert out["recipient_positive_incoming_flux"] >= out["recipient_net_gain"]
    assert set(out["terminal_negative_net_donor_roles"]) == {1, 2}
    assert out["maximum_shortest_donor_path_length"] <= 2
    assert out["same_physical_event"] is True
    assert out["same_physical_time"] is True
    assert out["new_causal_charge_created"] is False
    assert out["recursive_generation_created"] is False
    assert out["scale_progress_created"] is False
    assert out["primary_selected"] is False


def test_high_tail_interface_owner_gives_clean_quarter_tail_component_in_same_unit():
    out = high_tail_interface_component_route(
        2.0,
        0.5,
        positive_interface_common_work=0.6,
        positive_skew_common_work=0.31,
        positive_symmetric_common_work=0.29,
    )
    assert out["clean_interface_owner_lower"] == pytest.approx(0.5)
    assert out["clean_component_owner_lower"] == pytest.approx(0.25)
    assert out["owner_threshold"] == pytest.approx(0.3)
    assert out["joint_owners"] == (SKEW_OWNER,)
    assert out["shell_scale_reweighting_used"] is False


def test_high_tail_interface_exact_component_tie_remains_joint():
    out = high_tail_interface_component_route(1.0, 1.0, 0.6, 0.3, 0.3)
    assert set(out["joint_owners"]) == {SKEW_OWNER, SYMMETRIC_OWNER}
    assert out["primary_selected"] is False


def test_certificate_quotients_circulation_without_claiming_global_termination():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "not an identified measure" in cert["operator_bridge"]
    assert "same physical event time" in cert["finite_same_event"]
    assert "not a new source currency" in cert["symmetric_semantics"]
    assert "creates neither energy nor a second causal charge" in cert["skew_semantics"]
    assert "no scale progress" in cert["scale"]
    assert "does not prove" in cert["scope"]
