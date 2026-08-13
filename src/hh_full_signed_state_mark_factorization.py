from __future__ import annotations

import math
from dataclasses import dataclass

from src.canonical_positive_edge_work_routing import HardCellWork, HardProductCell
from src.complex_young_parent_marking import christ_complex_parent_mark_available

STATUS = (
    "DRAFT_HH_FULL_SIGNED_STATE_MARK_FACTORIZATION__"
    "CHRIST_MARK_IS_SHARED_EVENT_ROLE_WITNESS_NOT_CAUSAL_CURRENCY__"
    "CANONICAL_GOOD_DW_PLUS_INHERITS_METADATA_WITH_EXACT_MASS_PROJECTION__"
    "TERMINAL_BAD_DW_PLUS_NOT_REOPENED"
)


def _close(a: float, b: float, *, factor: float = 8.0e-10) -> bool:
    return abs(float(a) - float(b)) <= factor * max(abs(float(a)), abs(float(b)), 1.0e-300)


@dataclass(frozen=True)
class FullSignedRoleStateMark:
    cell: HardProductCell
    event_state_key: str
    signed_full_cell_work: float
    certified_full_cell_young_upper: float
    full_signed_state_deficit: float
    normalized_symbol_freezing_error: float
    christ_modulus: float
    role_state_mark_available: bool
    mark_reads_full_physical_roles: bool = True
    mark_is_causal_mass: bool = False
    mark_is_new_owner: bool = False
    mark_creates_event: bool = False
    mark_creates_recursion_depth: bool = False
    mark_creates_scale_progress: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        if not str(self.event_state_key):
            raise ValueError("nonempty physical event-state key required")
        vals = (
            self.signed_full_cell_work,
            self.certified_full_cell_young_upper,
            self.full_signed_state_deficit,
            self.normalized_symbol_freezing_error,
            self.christ_modulus,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite full-signed state-mark data required")
        if self.certified_full_cell_young_upper <= 0.0:
            raise ValueError("positive full-cell Young upper required")
        if self.normalized_symbol_freezing_error < 0.0:
            raise ValueError("nonnegative symbol-freezing error required")
        if not (0.0 < self.christ_modulus < 1.0):
            raise ValueError("external Christ modulus must lie in (0,1)")
        if not self.mark_reads_full_physical_roles:
            raise ValueError("Christ mark must read the unchanged full physical hard-role state")
        if (
            self.mark_is_causal_mass
            or self.mark_is_new_owner
            or self.mark_creates_event
            or self.mark_creates_recursion_depth
            or self.mark_creates_scale_progress
            or self.later_hahn_used
        ):
            raise ValueError("state mark was promoted into causal/event semantics")
        expected = 1.0 - self.signed_full_cell_work / self.certified_full_cell_young_upper
        if not _close(self.full_signed_state_deficit, expected, factor=5.0e-10):
            raise AssertionError("full-signed state deficit changed")
        expected_gate = (
            self.signed_full_cell_work > 0.0
            and self.full_signed_state_deficit + self.normalized_symbol_freezing_error <= self.christ_modulus
        )
        if self.role_state_mark_available != expected_gate:
            raise AssertionError("state-mark availability is not the ordinary full-signed Christ gate")


@dataclass(frozen=True)
class GoodCausalStateMarkFactorization:
    cell: HardProductCell
    event_state_key: str
    canonical_good_positive_work: float
    canonical_bad_positive_work: float
    canonical_negative_work: float
    marked_good_positive_work: float
    projected_good_positive_work: float
    state_mark: FullSignedRoleStateMark
    bad_work_remains_terminal: bool = True
    canonical_good_causal_law_reweighted: bool = False
    canonical_bad_causal_law_reopened: bool = False
    negative_work_used_as_payment: bool = False
    state_mark_charged_as_currency: bool = False
    causal_fate_restriction_applied_to_velocity_field: bool = False
    good_only_velocity_field_synthesized: bool = False
    hard_role_projectors_changed_by_routing: bool = False
    creates_new_owner: bool = False
    creates_new_event: bool = False
    creates_recursion_depth: bool = False
    creates_scale_progress: bool = False

    def __post_init__(self) -> None:
        if not str(self.event_state_key):
            raise ValueError("nonempty physical event-state key required")
        numeric = (
            self.canonical_good_positive_work,
            self.canonical_bad_positive_work,
            self.canonical_negative_work,
            self.marked_good_positive_work,
            self.projected_good_positive_work,
        )
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in numeric):
            raise ValueError("finite nonnegative causal work masses required")
        if self.state_mark.cell != self.cell:
            raise ValueError("state mark may only annotate the same physical hard cell")
        if self.state_mark.event_state_key != self.event_state_key:
            raise ValueError("state mark may only annotate the same physical event state")
        if not self.state_mark.role_state_mark_available:
            raise ValueError("cannot attach an unavailable full-signed role-state mark")
        if not self.canonical_good_positive_work > 0.0:
            raise ValueError("state-mark factorization requires positive canonical good dW+ mass")
        if not _close(self.marked_good_positive_work, self.canonical_good_positive_work):
            raise AssertionError("metadata lift changed canonical good dW+ mass")
        if not _close(self.projected_good_positive_work, self.canonical_good_positive_work):
            raise AssertionError("forgetting state metadata did not recover canonical good dW+")
        if not self.bad_work_remains_terminal:
            raise ValueError("state marking may not reopen terminal bad dW+")
        if (
            self.canonical_good_causal_law_reweighted
            or self.canonical_bad_causal_law_reopened
            or self.negative_work_used_as_payment
            or self.state_mark_charged_as_currency
            or self.causal_fate_restriction_applied_to_velocity_field
            or self.good_only_velocity_field_synthesized
            or self.hard_role_projectors_changed_by_routing
            or self.creates_new_owner
            or self.creates_new_event
            or self.creates_recursion_depth
            or self.creates_scale_progress
        ):
            raise ValueError("state-mark factorization changed physical causal semantics")


def certify_full_signed_role_state_mark(
    cell: HardCellWork,
    *,
    certified_full_cell_young_upper: float,
    normalized_symbol_freezing_error: float,
    christ_modulus: float,
    event_state_key: str,
) -> FullSignedRoleStateMark:
    """Read ordinary Young/Christ state information from the unchanged hard roles.

    No reserved deficit and no Hahn-positive sublaw is used to manufacture the
    mark.  The causal routing of dW+ remains a separate ledger.
    """
    Y = float(certified_full_cell_young_upper)
    xi = float(normalized_symbol_freezing_error)
    modulus = float(christ_modulus)
    if not all(math.isfinite(v) for v in (Y, xi, modulus)):
        raise ValueError("finite Young upper, symbol error and Christ modulus required")
    if Y <= 0.0 or xi < 0.0 or not (0.0 < modulus < 1.0):
        raise ValueError("invalid full-signed Christ data")
    T = float(cell.signed_work)
    if abs(T) > Y + 5.0e-10 * max(abs(T), Y, 1.0e-300):
        raise AssertionError("certified full-cell Young upper does not bound signed work")
    deficit = 1.0 - T / Y
    available = T > 0.0 and christ_complex_parent_mark_available(
        weighted_deficit=deficit,
        normalized_symbol_freezing_error=xi,
        christ_modulus_for_target_distance=modulus,
    )
    return FullSignedRoleStateMark(
        cell=cell.cell,
        event_state_key=str(event_state_key),
        signed_full_cell_work=T,
        certified_full_cell_young_upper=Y,
        full_signed_state_deficit=deficit,
        normalized_symbol_freezing_error=xi,
        christ_modulus=modulus,
        role_state_mark_available=available,
    )


def factor_role_state_mark_onto_good_causal_work(
    cell: HardCellWork,
    state_mark: FullSignedRoleStateMark,
    *,
    event_state_key: str,
) -> GoodCausalStateMarkFactorization:
    """Lift metadata onto canonical good dW+ and preserve its exact projection."""
    if state_mark.cell != cell.cell:
        raise ValueError("cannot move a state mark across hard-role cells")
    if state_mark.event_state_key != str(event_state_key):
        raise ValueError("cannot move a state mark across physical event states")
    if not state_mark.role_state_mark_available:
        raise ValueError("full signed hard-role state did not pass Christ")
    good = float(cell.inherited_good_positive_work)
    if not good > 0.0:
        raise ValueError("no canonical good positive work is present to annotate")
    return GoodCausalStateMarkFactorization(
        cell=cell.cell,
        event_state_key=str(event_state_key),
        canonical_good_positive_work=good,
        canonical_bad_positive_work=float(cell.inherited_bad_positive_work),
        canonical_negative_work=float(cell.inherited_negative_work),
        marked_good_positive_work=good,
        projected_good_positive_work=good,
        state_mark=state_mark,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "type_separation": "Young/Christ marks the shared full hard-role state; canonical dW+ is the positive causal work law",
        "ns_restriction_invariance": "dW_G+=1_G dW+ restricts the interaction/work measure only; it does not act on u(t), synthesize a good-only field, or change the event projectors P_a",
        "state_gate": "uses christ_complex_parent_mark_available from the existing complex-Young theorem on delta_state=1-T_C/Y_C and the same Xi; PR7 does not weaken the analytic premise",
        "factorization": "on one fixed event-state fiber, mu_G,C lifts to (id,M_C)_# mu_G,C and projection forgetting M_C returns exactly mu_G,C",
        "bad_semantics": "canonical bad dW+ remains on existing TRANSFER_WORK_LOSS and is never reopened by sharing the same event state",
        "negative_semantics": "canonical dW- remains signed evidence/donor provenance and is never a payment for good work",
        "reservation": "T_C-b_C is a stronger counterfactual robustness certificate but is not required for state marking",
        "christ_margin": "the spare-margin variable disappears from state marking: full Christ pass gives a mark; full Christ fail takes the existing Young/symbol/transfer route",
        "downstream": "dual Gaussian is an analysis probe of the role; Bargmann is an energy identity anchor; physical productivity stays weighted by actual retained dW+",
        "coherent_scope": "no fresh coherent Hahn law or general coherent positive kernel is inferred",
        "claims_generic_hh_termination": False,
        "claims_mixed_owner_termination": False,
        "claims_global_regularity": False,
    }
