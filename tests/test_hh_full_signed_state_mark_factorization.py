import pytest

from src.canonical_positive_edge_work_routing import HardCellWork, HardProductCell
from src.hh_full_signed_state_mark_factorization import (
    STATUS,
    certify_full_signed_role_state_mark,
    factor_role_state_mark_onto_good_causal_work,
    theorem_certificate,
)
from src.mixed_fate_reserved_young_handoff import certify_reserved_young_handoff


def _cell(g: float, b: float, n: float, label: str) -> HardCellWork:
    T = g + b - n
    P = g + b
    fresh = max(T, 0.0)
    return HardCellWork(
        cell=HardProductCell(parent_roles=(f"{label}-p1", f"{label}-p2"), child_role=f"{label}-c"),
        signed_work=T,
        inherited_positive_work=P,
        inherited_negative_work=n,
        inherited_good_positive_work=g,
        inherited_bad_positive_work=b,
        fresh_cell_hahn_positive=fresh,
        cancellation_gap=P - fresh,
    )


def test_certificate_is_state_causal_type_separation():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "shared full hard-role state" in cert["type_separation"]
    assert "projection" in cert["factorization"]
    assert "never reopened" in cert["bad_semantics"]
    assert "never a payment" in cert["negative_semantics"]
    assert not cert["claims_global_regularity"]


def test_state_metadata_lift_projects_exactly_to_good_dW_plus():
    cell = _cell(0.82, 0.0, 0.02, "pure")
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=0.84,
        normalized_symbol_freezing_error=0.005,
        christ_modulus=0.06,
        event_state_key="event-pure",
    )
    assert mark.role_state_mark_available
    out = factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key="event-pure")
    assert out.marked_good_positive_work == pytest.approx(0.82)
    assert out.projected_good_positive_work == pytest.approx(0.82)
    assert not out.state_mark_charged_as_currency
    assert not out.creates_recursion_depth


def test_reservation_can_fail_while_the_actual_full_role_state_is_christ_marked():
    # T/Y=0.95 passes the ordinary full-state Christ gate at modulus 0.10.
    # But removing b=0.25 counterfactually leaves T-b=0.70, so the predecessor
    # reserved gate fails.  The state theorem still holds because it describes
    # the actual roles that existed; the continuing causal mass remains g=0.70.
    cell = _cell(0.70, 0.25, 0.0, "mixed")
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.10,
        event_state_key="event-mixed",
    )
    assert mark.role_state_mark_available

    reserved = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.10,
    )
    assert reserved.full_signed_christ_gate
    assert not reserved.young_handoff_certified

    out = factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key="event-mixed")
    assert out.canonical_good_positive_work == pytest.approx(0.70)
    assert out.canonical_bad_positive_work == pytest.approx(0.25)
    assert out.marked_good_positive_work == pytest.approx(0.70)
    assert out.bad_work_remains_terminal
    assert not out.canonical_bad_causal_law_reopened


def test_negative_work_remains_signed_evidence_not_payment():
    cell = _cell(0.90, 0.08, 0.04, "negative")
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.01,
        christ_modulus=0.08,
        event_state_key="event-negative",
    )
    assert mark.role_state_mark_available
    out = factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key="event-negative")
    assert out.canonical_negative_work == pytest.approx(0.04)
    assert not out.negative_work_used_as_payment
    assert out.projected_good_positive_work == pytest.approx(0.90)


def test_full_christ_failure_fails_closed_without_state_mark():
    cell = _cell(0.60, 0.05, 0.05, "fail")
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.01,
        christ_modulus=0.10,
        event_state_key="event-fail",
    )
    assert not mark.role_state_mark_available
    with pytest.raises(ValueError, match="did not pass Christ"):
        factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key="event-fail")


def test_state_mark_cannot_cross_to_another_hard_role_state():
    cell1 = _cell(0.95, 0.0, 0.0, "one")
    cell2 = _cell(0.95, 0.0, 0.0, "two")
    mark = certify_full_signed_role_state_mark(
        cell1,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.10,
        event_state_key="event-one",
    )
    with pytest.raises(ValueError, match="across hard-role cells"):
        factor_role_state_mark_onto_good_causal_work(cell2, mark, event_state_key="event-one")



def test_state_mark_cannot_cross_time_with_same_hard_role_labels():
    cell = _cell(0.95, 0.0, 0.0, "same-role")
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.10,
        event_state_key="event-t0",
    )
    with pytest.raises(ValueError, match="across physical event states"):
        factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key="event-t1")
