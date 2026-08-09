import math

import pytest

from src.affine_sgs_boundary_ledger import sgs_increment_cubic_upper
from src.coherent_increment_service import cubic_to_square_threshold
from src.objective_source_routing_compiler import (
    OWNER_NAMES,
    STATUS,
    compile_objective_source_owners,
    local_dv_reentry,
    objective_owner_weight_threshold,
    objective_sgs_aggregate_route,
    objective_sgs_episode_thresholds,
    objective_sgs_high_frequency_physical_reentry,
    objective_sgs_square_service_per_source,
    pressure_source_alternatives,
    theorem_certificate,
    viscous_dv_reentry,
)
from src.resolved_objective_strain_collision import sgs_gradient_stress_lower


def test_exact_sgs_power_cancellation_closed_form():
    g1 = 1.7
    clp = 1.3
    cb = 1.2
    rho = 0.021
    r = sgs_gradient_stress_lower(rho)
    q = r**1.5 / sgs_increment_cubic_upper(g1, 1.0)
    composed = cubic_to_square_threshold(q, g1, clp, cb)
    closed = objective_sgs_square_service_per_source(g1, clp, cb) * rho
    assert composed == pytest.approx(closed, rel=2e-14, abs=1e-14)


def test_joint_owner_tie_is_not_lexicographic():
    c = 0.8
    A = 0.32
    sigma = objective_owner_weight_threshold(A, c)
    weights = {k: sigma for k in OWNER_NAMES}
    out = compile_objective_source_owners(
        A,
        c,
        weights,
        viscosity=1.0,
        filter_l1=1.2,
        lp_constant=1.1,
        bernstein_constant=1.05,
    )
    assert set(out["joint_owners"]) == set(OWNER_NAMES)
    assert out["additive_reset_created"] is False
    assert out["packet_synchronization_created"] is False


def test_local_and_viscous_feed_generic_shell_but_remain_recursive():
    loc = local_dv_reentry(0.03, 0.5, 1.0)
    vis = viscous_dv_reentry(0.01, 0.5, 1.0)
    for out in (loc, vis):
        assert out["resolved_DV_lower"] > 0
        assert out["critical_shell_mass_lower"] > 0
        assert out["own_scale_service_lower_on_full_survivor"] > 0
        assert out["master_semantics"] == "RECURSE_CRITICAL"


def test_pressure_mass_is_not_silently_promoted_to_shell():
    sigma = 0.04
    out = pressure_source_alternatives(sigma, 0.5, 1.2, 1.1, 1.05)
    assert out["resolved_lowpass_mass_occupation_lower"] / 5700.0 == pytest.approx(sigma / 2.0)
    assert out["stress_l32_occupation_lower"] / 380.0 == pytest.approx(sigma / 2.0)
    assert out["effective_sgs_source_weight_if_stress_branch"] == pytest.approx(sigma / 2.0)
    assert out["resolved_mass_is_generic_critical_shell"] == "NO"
    assert out["fixed_pair_service_ratio_upper"] < 1.0 / 3.0
    assert out["fixed_pair_total_future_multiplier_upper"] < 1.5


def test_integrated_sgs_high_frequency_owner_is_not_DV():
    th = objective_sgs_episode_thresholds(0.02, 0.7, 1.2, 1.1, 1.05)
    Y = th["integrated_forced_square_service"]
    out = objective_sgs_aggregate_route(
        0.02,
        0.7,
        1.2,
        1.1,
        1.05,
        high_frequency_dissipation=Y / 4.0,
        old_pool_integrated_capacity=0.0,
        old_old_integrated_service=0.0,
        selected_interface_integrated_service=0.0,
        new_edge_integrated_services=[],
    )
    assert out["branch"] == "high_frequency_dissipation"
    assert out["resolved_DV_supplier"] == "NO"


def test_integrated_sgs_dominant_fresh_edge_supplies_generic_shell():
    sigma = 0.02
    c = 0.7
    th = objective_sgs_episode_thresholds(sigma, c, 1.2, 1.1, 1.05)
    Y = th["integrated_forced_square_service"]
    # d=0 -> low service lower Y. Keep old and Xi below Y/8; fresh realizes the rest.
    w = [0.30 * Y, 0.20 * Y, 0.20 * Y, 0.20 * Y]
    out = objective_sgs_aggregate_route(
        sigma,
        c,
        1.2,
        1.1,
        1.05,
        high_frequency_dissipation=0.0,
        old_pool_integrated_capacity=0.05 * Y,
        old_old_integrated_service=0.05 * Y,
        selected_interface_integrated_service=0.05 * Y,
        new_edge_integrated_services=w,
    )
    assert out["branch"] == "dominant_fresh_critical_shell"
    assert out["clean_peak_whole_shell_mass_lower"] == pytest.approx(Y / (64.0 * c))
    assert out["master_semantics"] == "RECURSE_CRITICAL_VIA_GENERIC_SHELL"


def test_integrated_sgs_entropy_or_cycle_when_no_dominant_edge():
    sigma = 0.02
    c = 0.7
    th = objective_sgs_episode_thresholds(sigma, c, 1.2, 1.1, 1.05)
    Y = th["integrated_forced_square_service"]
    w = [0.18 * Y] * 5
    entropy = objective_sgs_aggregate_route(
        sigma,
        c,
        1.2,
        1.1,
        1.05,
        high_frequency_dissipation=0.0,
        old_pool_integrated_capacity=0.05 * Y,
        old_old_integrated_service=0.05 * Y,
        selected_interface_integrated_service=0.05 * Y,
        new_edge_integrated_services=w,
    )
    assert entropy["branch"] == "fresh_service_collision_entropy"
    assert entropy["entropy_lower"] == pytest.approx(math.log(4.0))

    cyc = objective_sgs_aggregate_route(
        sigma,
        c,
        1.2,
        1.1,
        1.05,
        high_frequency_dissipation=0.0,
        old_pool_integrated_capacity=0.05 * Y,
        old_old_integrated_service=0.05 * Y,
        selected_interface_integrated_service=0.05 * Y,
        new_edge_integrated_services=w,
        ancestry_labels=["A", "A", "A", "A", "A"],
    )
    assert cyc["branch"] == "fresh_service_same_ancestry_cycle"
    assert cyc["hidden_pair_lower"] == pytest.approx(0.25)


def test_certificate_records_two_forbidden_identifications():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    text = cert["forbidden_identifications"]
    assert "pressure low-pass mass is not generic critical-shell mass" in text
    assert "high-frequency SGS dissipation is not resolved D_V" in text


def test_high_frequency_owner_uses_tail_energy_not_resolved_DV():
    D = 0.4
    nu = 1.0
    c_lp = 0.25
    threshold = nu * c_lp * D
    out = objective_sgs_high_frequency_physical_reentry(
        D,
        nu,
        inherited_scaled_tail_energy=threshold,
        positive_scaled_tail_work=threshold,
    )
    assert set(out["energy_gate"]["joint_owners"]) == {
        "inherited_tail_energy",
        "positive_nonlinear_regeneration",
    }
    assert out["resolved_DV_supplier"] is False
    assert out["lp_to_physical_tail_lower"] == pytest.approx(c_lp)
    assert out["physical_tail_dissipation_lower"] == pytest.approx(c_lp * D)
    assert out["clean_thresholds"]["inherited_critical_shell_mass"] == pytest.approx(threshold)
