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
    objective_sgs_aggregate_scale_route,
    objective_sgs_episode_thresholds,
    objective_sgs_high_frequency_physical_reentry,
    objective_sgs_square_service_per_source,
    pressure_canonical_interface,
    pressure_source_alternatives,
    realized_pressure_source_route,
    theorem_certificate,
    viscous_dv_reentry,
)
from src.resolved_objective_strain_collision import sgs_gradient_stress_lower
from src.fresh_service_scale_reentry import pushforward_fresh_edges_to_bands


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
    assert out["routes"]["pressure"]["aggregate_muV_canonical_route"] is False


def test_local_and_viscous_feed_generic_shell_but_remain_recursive():
    loc = local_dv_reentry(0.03, 0.5, 1.0)
    vis = viscous_dv_reentry(0.01, 0.5, 1.0)
    for out in (loc, vis):
        assert out["resolved_DV_lower"] > 0
        assert out["critical_shell_mass_lower"] > 0
        assert out["own_scale_service_lower_on_full_survivor"] > 0
        assert out["master_semantics"] == "RECURSE_CRITICAL"


def test_pressure_muV_split_survives_only_as_diagnostic():
    sigma = 0.04
    out = pressure_source_alternatives(sigma, 0.5, 1.2, 1.1, 1.05)
    assert out["resolved_lowpass_mass_occupation_lower"] / 5700.0 == pytest.approx(sigma / 2.0)
    assert out["stress_l32_occupation_lower"] / 380.0 == pytest.approx(sigma / 2.0)
    assert out["resolved_mass_is_generic_critical_shell"] == "NO"
    assert out["canonical_pressure_route"] is False
    assert out["master_semantics"] == "DIAGNOSTIC_ONLY"
    assert out["fixed_pair_service_ratio_upper"] < 1.0 / 5.0


def test_symbolic_pressure_compiler_has_no_aggregate_muV_state():
    sigma = 0.04
    out = pressure_canonical_interface(sigma, 0.5, 1.2, 1.1, 1.05)
    assert out["positive_source_owner_threshold"] == pytest.approx(sigma / 2.0)
    assert out["sgs_stress_occupation_lower_if_sgs_owner"] == pytest.approx(190.0 * sigma)
    assert out["pair_quarter_shell_corollary_lower"] == pytest.approx(80.0 * sigma / 0.5)
    assert out["pair_quarter_entropy_corollary_lower"] == pytest.approx(math.log(4.0))
    assert out["aggregate_muV_canonical_route"] is False
    assert "resolved_lowpass_mass_occupation_lower" not in out
    assert out["realized_positive_source_law_required"] is True


def test_realized_diffuse_pressure_pair_law_still_enters_critical_shell():
    sigma = 1.0
    out = realized_pressure_source_route(
        sigma,
        1.0,
        1.0,
        1.2,
        1.1,
        1.05,
        block_frequency=4.0,
        sgs_positive_source_weight=0.0,
        pair_positive_weights=[0.2] * 5,
        pair_shell_indices=[(i, i) for i in range(5)],
        pair_frequencies=[(1.0 / (2**i), 1.0 / (2**i)) for i in range(5)],
    )
    assert out["joint_primary_owners"] == ("resolved_pressure_pair_law",)
    rr = out["routes"]["resolved_pressure_pair_law"]
    assert rr["critical_shell_mass_lower"] > 0
    assert rr["full_survivor_own_scale_service_lower"] > 0
    assert rr["full_survivor_integrated_service_lower"] > 0
    assert rr["entropy_weighted_critical_shell_mass_lower"] >= 320.0 - 1e-10
    assert rr["entropy_weighted_full_survivor_service_lower"] >= rr["clean_entropy_weighted_full_survivor_service_lower"] - 1e-12
    assert rr["entropy_weighted_full_survivor_integrated_service_lower"] >= rr["clean_entropy_weighted_full_survivor_integrated_service_lower"] - 1e-12
    assert rr["next_owner"] == "generic_critical_shell_first_stop"
    assert out["aggregate_muV_used"] is False


def test_realized_pressure_sgs_owner_uses_actual_positive_weight():
    sigma = 1.0
    out = realized_pressure_source_route(
        sigma,
        1.0,
        1.0,
        1.2,
        1.1,
        1.05,
        block_frequency=4.0,
        sgs_positive_source_weight=1.0,
        pair_positive_weights=[],
        pair_shell_indices=[],
        pair_frequencies=[],
    )
    assert out["joint_primary_owners"] == ("sgs_pressure_source",)
    rr = out["routes"]["sgs_pressure_source"]
    assert rr["stress_l32_occupation_lower"] == pytest.approx(380.0)
    assert rr["integrated_forced_square_service_lower"] > 0
    assert rr["next_owner"] == "coherent_service"


def test_realized_pressure_exact_half_tie_keeps_both_owners():
    sigma = 1.0
    out = realized_pressure_source_route(
        sigma,
        1.0,
        1.0,
        1.2,
        1.1,
        1.05,
        block_frequency=4.0,
        sgs_positive_source_weight=0.5,
        pair_positive_weights=[0.5],
        pair_shell_indices=[(0, 0)],
        pair_frequencies=[(1.0, 1.0)],
    )
    assert set(out["joint_primary_owners"]) == {"sgs_pressure_source", "resolved_pressure_pair_law"}
    assert set(out["routes"]) == {"sgs_pressure_source", "resolved_pressure_pair_law"}


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


def test_legacy_integrated_sgs_dominant_fresh_edge_is_sideledger_only():
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
    assert out["canonical_renewal_fate"] is False
    assert out["master_semantics"] == "SIDELEDGER_ONLY__LEGACY_CELL_CLUSTER"


def test_legacy_integrated_sgs_entropy_or_cycle_is_sideledger_only():
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
    assert entropy["canonical_renewal_fate"] is False
    assert entropy["master_semantics"] == "SIDELEDGER_ONLY__LEGACY_CELL_ENTROPY"

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
    assert cyc["canonical_renewal_fate"] is False
    assert cyc["master_semantics"] == "SIDELEDGER_ONLY__LEGACY_CELL_CYCLE"



def test_canonical_sgs_fresh_route_uses_band_pushforward_not_cell_argmax():
    sigma = 0.02
    c = 0.7
    g1, clp, cb = 1.2, 1.1, 1.05
    th = objective_sgs_episode_thresholds(sigma, c, g1, clp, cb)
    Y = th["integrated_forced_square_service"]
    law = {0: 0.30 * Y, -1: 0.30 * Y, -2: 0.30 * Y}
    out = objective_sgs_aggregate_scale_route(
        sigma, c, 1.0, g1, clp, cb,
        block_frequency=8.0,
        high_frequency_dissipation=0.0,
        old_pool_integrated_capacity=0.05 * Y,
        old_old_integrated_service=0.05 * Y,
        selected_interface_integrated_service=0.05 * Y,
        fresh_band_integrated_services=law,
    )
    assert out["joint_primary_owners"] == ("fresh_scale_critical_shell",)
    assert out["coherent_cell_priority_used"] is False
    fresh = out["routes"]["fresh_scale_critical_shell"]
    assert fresh["coherent_cell_argmax_used"] is False
    assert fresh["cell_ancestry_sideledger_optional"] is True
    assert fresh["critical_shell_mass_lower"] > 0
    assert fresh["master_semantics"] == "RECURSE_CRITICAL_VIA_REFINEMENT_INVARIANT_SCALE_SHELL"
    scale = fresh["scale_route"]
    assert scale["H_inf_weighted_hard_shell_mass_lower"] >= Y / (24.0 * c) - 1e-13


def test_canonical_sgs_fresh_route_is_invariant_under_coherent_cell_refinement():
    sigma = 0.02
    c = 0.7
    g1, clp, cb = 1.2, 1.1, 1.05
    Y = objective_sgs_episode_thresholds(sigma, c, g1, clp, cb)["integrated_forced_square_service"]
    coarse = pushforward_fresh_edges_to_bands(
        [0.30 * Y, 0.30 * Y, 0.30 * Y],
        [0, -1, -2],
        [False, False, False],
        [False, False, False],
    )
    fine = pushforward_fresh_edges_to_bands(
        [0.10 * Y, 0.20 * Y, 0.12 * Y, 0.18 * Y, 0.05 * Y, 0.25 * Y],
        [0, 0, -1, -1, -2, -2],
        [False] * 6,
        [False] * 6,
    )
    assert coarse == pytest.approx(fine)
    common = dict(
        source_weight=sigma, scaled_lifetime=c, viscosity=1.0,
        filter_l1=g1, lp_constant=clp, bernstein_constant=cb,
        block_frequency=8.0, high_frequency_dissipation=0.0,
        old_pool_integrated_capacity=0.05 * Y, old_old_integrated_service=0.05 * Y,
        selected_interface_integrated_service=0.05 * Y,
    )
    a = objective_sgs_aggregate_scale_route(**common, fresh_band_integrated_services=coarse)
    b = objective_sgs_aggregate_scale_route(**common, fresh_band_integrated_services=fine)
    ar = a["routes"]["fresh_scale_critical_shell"]
    br = b["routes"]["fresh_scale_critical_shell"]
    assert ar["critical_shell_mass_lower"] == pytest.approx(br["critical_shell_mass_lower"])
    assert ar["H_inf_scale"] == pytest.approx(br["H_inf_scale"])
    assert ar["H2_scale"] == pytest.approx(br["H2_scale"])



def test_canonical_sgs_realized_ties_keep_all_physical_owners():
    sigma = 0.02
    c = 0.7
    g1, clp, cb = 1.2, 1.1, 1.05
    Y = objective_sgs_episode_thresholds(sigma, c, g1, clp, cb)["integrated_forced_square_service"]
    out = objective_sgs_aggregate_scale_route(
        sigma, c, 1.0, g1, clp, cb,
        block_frequency=8.0,
        high_frequency_dissipation=0.25 * Y,
        old_pool_integrated_capacity=0.13 * Y,
        old_old_integrated_service=0.13 * Y,
        selected_interface_integrated_service=0.125 * Y,
        fresh_band_integrated_services={0: 0.10 * Y, -1: 0.08 * Y, -2: 0.07 * Y},
    )
    assert set(out["joint_primary_owners"]) == {
        "high_frequency_dissipation",
        "old_pool_not_yet_eroded",
        "selected_interface_Xi",
        "fresh_scale_critical_shell",
    }
    assert out["coherent_cell_priority_used"] is False
    assert out["master_semantics"] == "JOINT_NATIVE_OWNERS__NO_LEXICOGRAPHIC_PRIORITY"


def test_symbolic_sgs_interface_marks_cell_dominance_noncanonical():
    th = objective_sgs_episode_thresholds(0.02, 0.7, 1.2, 1.1, 1.05)
    assert th["fresh_cell_dominance_is_canonical_renewal"] == "NO"
    assert "FRESH_MATERIAL_SERVICE_TO_REFINEMENT_INVARIANT_SCALE_LAW" in th["canonical_fresh_route"]


def test_certificate_records_pressure_pair_native_route_and_forbidden_identifications():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "mu_child exp(H2_pair)>=320 Sigma_P/c" in cert["pressure_pair_route"]
    assert "diagnostic only" in cert["pressure_legacy_muV"]
    assert "full no-hit shell corridor" in cert["pressure_service_conjugacy"]
    text = cert["forbidden_identifications"]
    assert "aggregate pressure mu_V is not a canonical renewal state" in text
    assert "pressure H2 is not a causal child-energy probability" in text
    assert "high-frequency SGS dissipation is not resolved D_V" in text
    assert "coherent-cell entropy is not a canonical fresh renewal fate" in text
    assert "mu_hard exp(H_inf_scale)>=Y/(24c)" in cert["sgs_fresh_scale_route"]
    assert "retained jointly" in cert["sgs_joint_owner_rule"]


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
