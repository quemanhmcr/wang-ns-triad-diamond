from dataclasses import replace

import pytest

from src.canonical_positive_edge_work_routing import (
    HardCellWork,
    HardProductCell,
    _near_extremal_positive_fiber,
    _nonforward_positive_fiber,
    route_canonical_positive_edge_work,
    single_hard_role_map,
)
from src.continuum_helical_edge_measure_registration import continuum_edge_measure_ledger
from src.mixed_fate_reserved_young_handoff import (
    STATUS,
    certify_low_deficit_work_contamination,
    certify_reserved_young_handoff,
    coherent_fresh_hahn_kernel_counterexample,
    good_work_capacity_density_lower,
    inherited_negative_work,
    low_deficit_contamination_ratio_upper,
    reserved_failure_good_fraction_upper,
    reserved_failure_good_work_upper,
    phase_fate_role_refinement_counterexample,
    positive_subtraction_not_trilinear_counterexample,
    theorem_certificate,
)


def _cell(*, g: float, b: float, n: float) -> HardCellWork:
    T = g + b - n
    P = g + b
    fresh = max(T, 0.0)
    return HardCellWork(
        cell=HardProductCell(parent_roles=("p0", "p1"), child_role="c"),
        signed_work=T,
        inherited_positive_work=P,
        inherited_negative_work=n,
        inherited_good_positive_work=g,
        inherited_bad_positive_work=b,
        fresh_cell_hahn_positive=fresh,
        cancellation_gap=P - fresh,
    )


def test_mixed_fate_reserved_gate_uses_full_signed_work_but_needs_no_bad_assistance():
    cell = _cell(g=0.99, b=0.005, n=0.002)
    out = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.001,
        christ_modulus=0.02,
    )
    assert out.inherited_bad_positive_work > 0.0
    assert out.inherited_negative_work == pytest.approx(0.002)
    assert out.signed_full_cell_work == pytest.approx(0.993)
    assert out.reservation_value == pytest.approx(0.988)
    assert out.full_signed_young_deficit == pytest.approx(0.007)
    assert out.terminal_bad_assistance_fraction == pytest.approx(0.005)
    assert out.reserved_young_deficit == pytest.approx(0.012)
    assert out.reserved_young_deficit == pytest.approx(
        out.full_signed_young_deficit + out.terminal_bad_assistance_fraction
    )
    assert out.full_signed_christ_gate
    assert out.reserved_christ_gate
    assert out.young_handoff_certified
    assert out.full_signed_trilinear_work_used
    assert not out.reservation_value_used_as_trilinear_work
    assert not out.inherited_good_causal_mass_changed
    assert not out.later_hahn_used_as_causal_law


def test_bad_positive_work_cannot_make_an_assisted_full_young_cell_pass_reserved_gate():
    # Full T_C exactly saturates the supplied Young upper only because b_C fills
    # the last four percent.  The theorem must refuse to bind the good cause.
    cell = _cell(g=0.96, b=0.04, n=0.0)
    out = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.02,
    )
    assert out.full_signed_young_deficit == pytest.approx(0.0)
    assert out.full_signed_christ_gate
    assert out.reserved_young_deficit == pytest.approx(0.04)
    assert not out.reserved_christ_gate
    assert not out.young_handoff_certified


def test_negative_work_stays_inside_the_signed_young_input_before_any_gate():
    cell = _cell(g=1.10, b=0.01, n=0.15)
    out = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.03,
    )
    assert inherited_negative_work(cell) == pytest.approx(0.15)
    assert out.signed_full_cell_work == pytest.approx(0.96)
    assert out.reservation_value == pytest.approx(0.95)
    assert out.full_signed_young_deficit == pytest.approx(0.04)
    assert out.reserved_young_deficit == pytest.approx(0.05)
    assert not out.young_handoff_certified


def test_fate_pure_case_is_exact_special_case_of_reserved_handoff():
    cell = _cell(g=0.99, b=0.0, n=0.005)
    out = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.001,
        christ_modulus=0.02,
    )
    assert out.terminal_bad_assistance_fraction == 0.0
    assert out.reserved_young_deficit == pytest.approx(out.full_signed_young_deficit)
    assert out.reserved_christ_gate == out.full_signed_christ_gate


def test_reserved_gate_exact_tie_is_deterministic_and_allowed():
    # T=7/8, b=1/16 gives reservation 13/16 and deficit 3/16.  With
    # xi=1/16 the reserved Christ gate lands exactly at 1/4.
    cell = _cell(g=0.875, b=0.0625, n=0.0625)
    out = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0625,
        christ_modulus=0.25,
    )
    assert out.reserved_young_deficit == 0.1875
    assert out.reserved_young_deficit + out.normalized_symbol_freezing_error == 0.25
    assert out.young_handoff_certified


def test_christ_modulus_is_capped_below_one_for_failure_alternative_semantics():
    cell = _cell(g=0.99, b=0.005, n=0.002)
    with pytest.raises(ValueError, match="Christ modulus"):
        certify_reserved_young_handoff(
            cell,
            certified_full_cell_young_upper=1.0,
            normalized_symbol_freezing_error=0.0,
            christ_modulus=1.0,
        )


def test_positive_rescaling_changes_no_dimensionless_handoff_deficit():
    base = _cell(g=0.99, b=0.005, n=0.002)
    scaled = _cell(g=99.0, b=0.5, n=0.2)
    a = certify_reserved_young_handoff(
        base,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.001,
        christ_modulus=0.02,
    )
    b = certify_reserved_young_handoff(
        scaled,
        certified_full_cell_young_upper=100.0,
        normalized_symbol_freezing_error=0.001,
        christ_modulus=0.02,
    )
    assert a.full_signed_young_deficit == pytest.approx(b.full_signed_young_deficit, abs=2e-15)
    assert a.terminal_bad_assistance_fraction == pytest.approx(b.terminal_bad_assistance_fraction, abs=2e-15)
    assert a.reserved_young_deficit == pytest.approx(b.reserved_young_deficit, abs=2e-15)


def test_low_native_deficit_converts_reference_capacity_to_actual_work_only_through_pointwise_rn_lower():
    kappa = good_work_capacity_density_lower()
    assert 0.0 < kappa < 1.0
    assert kappa == pytest.approx(0.19)
    assert low_deficit_contamination_ratio_upper(0.0) == 0.0
    assert low_deficit_contamination_ratio_upper(1.0e-8) > 0.0
    assert low_deficit_contamination_ratio_upper(5.0e-8) > low_deficit_contamination_ratio_upper(1.0e-8)
    assert low_deficit_contamination_ratio_upper(1.0e-4) == float("inf")


def test_physical_low_deficit_edge_law_bounds_actual_mixed_and_negative_work_without_changing_cause():
    # Use registered NS helical edges, with tiny bad/backscatter quadrature mass,
    # then maximally coarsen only to force the mixed-cell seam.  The low-deficit
    # theorem must compare actual dW masses, not normalized capacity probabilities.
    ledger = continuum_edge_measure_ledger((
        _near_extremal_positive_fiber(1.0),
        _nonforward_positive_fiber(1.0e-10),
        _nonforward_positive_fiber(1.0e-10, phase_sign=-1.0),
    ))
    assert 0.0 <= ledger.block_transfer_deficit < 1.0e-4
    roles = single_hard_role_map(ledger)
    routing = route_canonical_positive_edge_work(ledger, tau=0.1, mode_roles=roles)
    out = certify_low_deficit_work_contamination(ledger, routing)
    assert out.good_work == pytest.approx(routing.good_positive_work, rel=2e-10)
    assert out.bad_positive_work == pytest.approx(routing.bad_positive_work, rel=2e-10)
    assert out.negative_work == pytest.approx(ledger.negative_edge_work, rel=2e-10)
    assert out.actual_contamination_to_good_ratio <= out.universal_contamination_to_good_ratio_upper * (1.0 + 2e-8)
    assert not out.canonical_good_mass_changed
    assert not out.capacity_used_as_causal_law


def test_low_deficit_contamination_certificate_rejects_forged_summary():
    ledger = continuum_edge_measure_ledger((_near_extremal_positive_fiber(1.0),))
    routing = route_canonical_positive_edge_work(
        ledger,
        tau=0.1,
        mode_roles=single_hard_role_map(ledger),
    )
    forged = replace(ledger, block_transfer_deficit=ledger.block_transfer_deficit + 1.0e-6)
    with pytest.raises(AssertionError, match="replay field block_transfer_deficit"):
        certify_low_deficit_work_contamination(forged, routing)


def test_phase_only_fate_change_defeats_geometry_only_hard_role_purification():
    out = phase_fate_role_refinement_counterexample()
    assert out["good_is_good"]
    assert out["bad_is_bad"]
    assert out["same_normalized_geometric_multiplier"] > 1.0 - 1.0e-4
    assert out["good_phase"] != out["bad_phase"]
    assert not out["mode_role_geometry_changed"]


def test_subtracting_hahn_positive_bad_work_is_not_a_new_trilinear_form():
    out = positive_subtraction_not_trilinear_counterexample()
    assert out["b_of_sum"] == 0.0
    assert out["sum_of_b"] == 1.0
    assert not out["positive_part_is_additive"]
    assert not out["reservation_is_trilinear_work"]


def test_fresh_non_diagonal_coherent_hahn_cannot_be_the_canonical_kernel_pushforward():
    out = coherent_fresh_hahn_kernel_counterexample()
    assert out["canonical_positive_mass"] == 0.0
    assert out["fresh_localized_hahn_positive"] == pytest.approx(0.5)
    assert not out["fresh_hahn_can_be_positive_kernel_pushforward_of_zero_canonical_mass"]


def test_theorem_certificate_keeps_scope_local_and_causal_law_fixed():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "T_C-b_C" in cert["young_semantics"]
    assert "capacity is reference only" in cert["low_deficit_actual_work_bound"]
    assert not cert["claims_global_regularity"]


def test_reserved_failure_with_full_christ_margin_is_dominated_by_existing_bad_and_negative_work():
    cell = _cell(g=0.979, b=0.021, n=0.0)
    out = certify_reserved_young_handoff(cell, certified_full_cell_young_upper=1.0, normalized_symbol_freezing_error=0.0, christ_modulus=0.02)
    assert out.full_signed_christ_gate and not out.young_handoff_certified
    upper = reserved_failure_good_work_upper(out, full_signed_christ_margin_floor=0.02)
    assert out.inherited_good_positive_work <= upper * (1.0 + 2e-9)
    assert out.terminal_bad_assistance_fraction > 0.02


def test_low_native_deficit_bounds_fraction_of_full_margin_cells_that_can_fail_reservation():
    assert reserved_failure_good_fraction_upper(0.0, full_signed_christ_margin_floor=0.02) == 0.0
    bound = reserved_failure_good_fraction_upper(1.0e-8, full_signed_christ_margin_floor=0.02)
    contamination = low_deficit_contamination_ratio_upper(1.0e-8)
    assert 0.0 < contamination < bound < 1.0
    assert bound == pytest.approx(49.0 * contamination)
