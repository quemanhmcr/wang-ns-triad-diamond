import math

import pytest

from src.canonical_positive_edge_work_routing import HardCellWork, HardProductCell
from src.hh_full_signed_state_mark_factorization import (
    certify_full_signed_role_state_mark,
    factor_role_state_mark_onto_good_causal_work,
)
from src.hh_state_mark_generated_event_bridge import (
    STATUS,
    bridge_state_mark_to_generated_event,
    certify_physical_hh_generation_gate,
    coalesce_state_marked_parent_slots,
    compose_state_marked_event_with_master,
    theorem_certificate,
)
from src.mixed_fate_reserved_young_handoff import certify_reserved_young_handoff
from src.physical_branch_compiler import CauseHit, MasterDisposition, PhysicalCause


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


def _factor(cell: HardCellWork, *, key: str, Y: float, xi: float, christ: float):
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=Y,
        normalized_symbol_freezing_error=xi,
        christ_modulus=christ,
        event_state_key=key,
    )
    assert mark.role_state_mark_available
    return factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key=key)


def _energy_gate(key: str, actual_hh_work: float):
    return certify_physical_hh_generation_gate(
        event_state_key=key,
        terminal_energy=0.5,
        initial_energy=0.05,
        residual_positive_work=0.05,
        strain_action=0.0,
        actual_positive_hh_work=actual_hh_work,
    )


def test_certificate_states_exact_composition_not_new_compiler():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "existing recursive_physical_witness_constructor" in cert["compiler"]
    assert "canonical good dW+" in cert["mass_identity"]
    assert "coalesce" in cert["root_reuse"]


def test_mixed_cell_full_christ_mark_promotes_good_event_even_when_reservation_fails():
    cell = _cell(0.70, 0.25, 0.0, "mixed")
    factor = _factor(cell, key="event-mixed", Y=1.0, xi=0.0, christ=0.10)
    reserved = certify_reserved_young_handoff(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.0,
        christ_modulus=0.10,
    )
    assert reserved.full_signed_christ_gate
    assert not reserved.young_handoff_certified

    bridge = bridge_state_mark_to_generated_event(factor, pair_cell=0, registration_good=True, energy_generation=_energy_gate("event-mixed", 0.95))
    assert bridge.generated_event.marking_good
    assert bridge.generated_event.mass == pytest.approx(0.70)
    assert bridge.canonical_bad_positive_work == pytest.approx(0.25)
    assert not bridge.reservation_gate_required
    assert not bridge.bad_work_enters_generated_law

    out = compose_state_marked_event_with_master(bridge, pair_cells_upper=1)
    assert out.generated_master.continuation_mass == pytest.approx(0.70)
    assert out.terminal_bad_positive_work == pytest.approx(0.25)
    assert out.total_accounted_positive_work == pytest.approx(0.95)
    assert out.total_canonical_positive_work == pytest.approx(0.95)


def test_negative_work_stays_outside_generated_mass_but_remains_provenance():
    cell = _cell(0.90, 0.08, 0.04, "negative")
    factor = _factor(cell, key="event-neg", Y=1.0, xi=0.01, christ=0.08)
    bridge = bridge_state_mark_to_generated_event(factor, pair_cell=1, registration_good=True, energy_generation=_energy_gate("event-neg", 0.98))
    assert bridge.generated_event.mass == pytest.approx(0.90)
    assert bridge.canonical_negative_work == pytest.approx(0.04)
    assert not bridge.negative_work_enters_generated_law
    out = compose_state_marked_event_with_master(bridge, pair_cells_upper=2)
    assert out.total_accounted_positive_work == pytest.approx(0.98)


def test_state_mark_does_not_override_common_slice_failure():
    cell = _cell(0.95, 0.0, 0.0, "registration")
    factor = _factor(cell, key="event-reg", Y=1.0, xi=0.0, christ=0.10)
    with pytest.raises(ValueError, match="first-stop provenance"):
        bridge_state_mark_to_generated_event(factor, pair_cell=0, registration_good=False, energy_generation=_energy_gate("event-reg", 0.95))

    hit = CauseHit(0.2, PhysicalCause.RESOLVED_SOURCE, 1.0, "existing source stop")
    bridge = bridge_state_mark_to_generated_event(
        factor,
        pair_cell=0,
        registration_good=False,
        energy_generation=_energy_gate("event-reg", 0.95),
        physical_hits=(hit,),
    )
    out = compose_state_marked_event_with_master(bridge, pair_cells_upper=1)
    assert out.generated_master.continuation_mass == 0.0
    assert out.generated_master.master_mass[MasterDisposition.RECURSE_CRITICAL.value] == pytest.approx(0.95)


def test_shared_state_anchor_coalesces_in_existing_reuse_pushforward():
    slots = {"e0-p0": 0.2, "e0-p1": 0.3, "e1-p0": 0.4, "e1-p1": 0.1}
    anchors = {
        "e0-p0": (4, 2, 1),
        "e1-p0": (4, 2, 1),
        "e0-p1": (8, 0, 3),
        "e1-p1": (8, 0, 3),
    }
    pushed = coalesce_state_marked_parent_slots(slots, anchors)
    assert pushed[(4, 2, 1)] == pytest.approx(0.6)
    assert pushed[(8, 0, 3)] == pytest.approx(0.4)
    assert math.isclose(sum(pushed.values()), sum(slots.values()))
    assert len(pushed) == 2


def test_full_state_christ_failure_cannot_enter_adapter():
    cell = _cell(0.60, 0.05, 0.05, "fail")
    mark = certify_full_signed_role_state_mark(
        cell,
        certified_full_cell_young_upper=1.0,
        normalized_symbol_freezing_error=0.01,
        christ_modulus=0.10,
        event_state_key="event-fail",
    )
    assert not mark.role_state_mark_available
    with pytest.raises(ValueError):
        factor_role_state_mark_onto_good_causal_work(cell, mark, event_state_key="event-fail")


def test_state_mark_cannot_mint_generated_event_without_physical_energy_generation_branch():
    cell = _cell(0.80, 0.05, 0.0, "energy-guard")
    factor = _factor(cell, key="event-energy-guard", Y=1.0, xi=0.0, christ=0.20)
    with pytest.raises(ValueError, match="not on the certified physical HH-generation branch"):
        certify_physical_hh_generation_gate(
            event_state_key="event-energy-guard",
            terminal_energy=0.5,
            initial_energy=0.2,
            residual_positive_work=0.0,
            strain_action=0.0,
            actual_positive_hh_work=0.85,
        )


def test_energy_gate_must_belong_to_same_event_state():
    cell = _cell(0.80, 0.05, 0.0, "energy-key")
    factor = _factor(cell, key="event-energy-key", Y=1.0, xi=0.0, christ=0.20)
    wrong = _energy_gate("different-event", 0.85)
    with pytest.raises(ValueError, match="different event state"):
        bridge_state_mark_to_generated_event(
            factor, pair_cell=0, registration_good=True, energy_generation=wrong
        )
