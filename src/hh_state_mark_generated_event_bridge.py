from __future__ import annotations

import math
from dataclasses import dataclass

from src.bargmann_root_cell_registration import pushforward_parent_slot_weights
from src.hh_full_signed_state_mark_factorization import GoodCausalStateMarkFactorization
from src.physical_branch_compiler import CauseHit
from src.physical_energy_causal_bridge import route_physical_energy_causality
from src.recursive_physical_witness_constructor import (
    GeneratedMasterPartition,
    GeneratedPairEvent,
    RegenerationHit,
    compile_generated_pair_master_measure,
)

STATUS = (
    "DRAFT_HH_STATE_MARK_TO_GENERATED_EVENT_BRIDGE__"
    "PHYSICAL_ENERGY_GENERATION_GATE_REQUIRED_BEFORE_GENERATED_EVENT__"
    "FULL_SIGNED_ROLE_WITNESS_PROMOTES_MARKING_GOOD_ON_CANONICAL_GOOD_DW_PLUS_ONLY__"
    "EXISTING_MASTER_COMPILER_UNCHANGED__BAD_AND_NEGATIVE_WORK_EXCLUDED_FROM_GENERATED_MASS"
)


def _close(a: float, b: float, *, factor: float = 8.0e-10) -> bool:
    return abs(float(a) - float(b)) <= factor * max(abs(float(a)), abs(float(b)), 1.0e-300)


@dataclass(frozen=True)
class PhysicalHHGenerationGateWitness:
    event_state_key: str
    terminal_energy: float
    initial_energy: float
    residual_positive_work: float
    strain_action: float
    actual_positive_hh_work: float
    certified_physical_hh_work_lower: float
    clean_generation_threshold: float
    branch: str = "physical_high_high_transfer_generation"

    def __post_init__(self) -> None:
        if not str(self.event_state_key):
            raise ValueError("nonempty event-state key required for physical HH generation")
        vals = (
            self.terminal_energy, self.initial_energy, self.residual_positive_work,
            self.strain_action, self.actual_positive_hh_work,
            self.certified_physical_hh_work_lower, self.clean_generation_threshold,
        )
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in vals):
            raise ValueError("finite nonnegative physical-energy generation data required")
        if self.terminal_energy <= 0.0 or self.actual_positive_hh_work <= 0.0:
            raise ValueError("positive terminal energy and actual HH work required")
        if self.branch != "physical_high_high_transfer_generation":
            raise ValueError("only the physical-energy HH-generation branch may mint a generated-event witness")
        tol = 8.0e-10 * max(self.actual_positive_hh_work, self.certified_physical_hh_work_lower, 1.0e-300)
        if self.actual_positive_hh_work + tol < self.certified_physical_hh_work_lower:
            raise AssertionError("actual positive HH work is below the physical-energy generation lower bound")
        tol2 = 8.0e-10 * max(self.certified_physical_hh_work_lower, self.clean_generation_threshold, 1.0e-300)
        if self.certified_physical_hh_work_lower + tol2 < self.clean_generation_threshold:
            raise AssertionError("physical-energy gate lost the clean 8/15 generation threshold")


def certify_physical_hh_generation_gate(
    *,
    event_state_key: str,
    terminal_energy: float,
    initial_energy: float,
    residual_positive_work: float,
    strain_action: float,
    actual_positive_hh_work: float,
) -> PhysicalHHGenerationGateWitness:
    """Replay the certified physical-energy branch before any GeneratedPairEvent exists."""
    route = route_physical_energy_causality(
        terminal_energy=float(terminal_energy),
        initial_energy=float(initial_energy),
        residual_positive_work=float(residual_positive_work),
        strain_action=float(strain_action),
    )
    if route.get("branch") != "physical_high_high_transfer_generation":
        raise ValueError("selected event state is not on the certified physical HH-generation branch")
    return PhysicalHHGenerationGateWitness(
        event_state_key=str(event_state_key),
        terminal_energy=float(terminal_energy),
        initial_energy=float(initial_energy),
        residual_positive_work=float(residual_positive_work),
        strain_action=float(strain_action),
        actual_positive_hh_work=float(actual_positive_hh_work),
        certified_physical_hh_work_lower=float(route["physical_hh_work_lower"]),
        clean_generation_threshold=float(route["clean_threshold"]),
    )


@dataclass(frozen=True)
class StateMarkedGeneratedEventBridge:
    event_state_key: str
    pair_cell: int
    canonical_good_positive_work: float
    canonical_bad_positive_work: float
    canonical_negative_work: float
    generated_event: GeneratedPairEvent
    energy_generation: PhysicalHHGenerationGateWitness
    state_mark_promotes_marking_good: bool = True
    reservation_gate_required: bool = False
    bad_work_enters_generated_law: bool = False
    negative_work_enters_generated_law: bool = False
    causal_mass_reweighted: bool = False
    creates_new_owner: bool = False
    creates_new_event: bool = False
    creates_recursion_depth: bool = False
    creates_scale_progress: bool = False

    def __post_init__(self) -> None:
        if not str(self.event_state_key):
            raise ValueError("nonempty event-state key required")
        if self.pair_cell < 0:
            raise ValueError("nonnegative pair-cell index required")
        vals = (self.canonical_good_positive_work, self.canonical_bad_positive_work, self.canonical_negative_work)
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in vals):
            raise ValueError("finite nonnegative canonical work masses required")
        if not self.canonical_good_positive_work > 0.0:
            raise ValueError("generated bridge requires positive canonical good dW+ mass")
        if self.energy_generation.event_state_key != self.event_state_key:
            raise ValueError("physical-energy generation witness belongs to a different event state")
        if self.canonical_good_positive_work > self.energy_generation.actual_positive_hh_work + 8.0e-10*max(self.canonical_good_positive_work, self.energy_generation.actual_positive_hh_work, 1.0e-300):
            raise AssertionError("canonical good submeasure exceeds the generated physical HH work law")
        if self.generated_event.pair_cell != self.pair_cell:
            raise AssertionError("generated event changed the selected pair-cell label")
        if not self.generated_event.marking_good:
            raise AssertionError("full-signed role-state witness did not promote marking_good")
        if not _close(self.generated_event.mass, self.canonical_good_positive_work):
            raise AssertionError("generated event mass is not the canonical good dW+ restriction")
        if (
            not self.state_mark_promotes_marking_good
            or self.reservation_gate_required
            or self.bad_work_enters_generated_law
            or self.negative_work_enters_generated_law
            or self.causal_mass_reweighted
            or self.creates_new_owner
            or self.creates_new_event
            or self.creates_recursion_depth
            or self.creates_scale_progress
        ):
            raise ValueError("state-mark bridge changed causal/event semantics")


@dataclass(frozen=True)
class StateMarkedMasterComposition:
    bridge: StateMarkedGeneratedEventBridge
    generated_master: GeneratedMasterPartition
    terminal_bad_positive_work: float
    total_canonical_positive_work: float
    total_accounted_positive_work: float
    good_work_charged_once: bool = True
    bad_work_charged_once: bool = True
    state_mark_is_master_currency: bool = False

    def __post_init__(self) -> None:
        good = self.bridge.canonical_good_positive_work
        bad = self.bridge.canonical_bad_positive_work
        if not _close(self.generated_master.total_mass, good):
            raise AssertionError("existing generated master did not receive exactly canonical good dW+ mass")
        if not _close(self.terminal_bad_positive_work, bad):
            raise AssertionError("terminal bad work changed during state-mark/master composition")
        if not _close(self.total_canonical_positive_work, good + bad):
            raise AssertionError("canonical positive-work decomposition changed")
        if not _close(self.total_accounted_positive_work, self.generated_master.total_mass + bad):
            raise AssertionError("state-mark/master composition lost or duplicated positive work")
        if not _close(self.total_accounted_positive_work, self.total_canonical_positive_work):
            raise AssertionError("good generated law plus terminal bad law failed canonical dW+ reconstruction")
        if not self.good_work_charged_once or not self.bad_work_charged_once or self.state_mark_is_master_currency:
            raise ValueError("state mark was promoted into a second master charge")


def bridge_state_mark_to_generated_event(
    factorization: GoodCausalStateMarkFactorization,
    *,
    pair_cell: int,
    registration_good: bool,
    energy_generation: PhysicalHHGenerationGateWitness,
    physical_hits: tuple[CauseHit, ...] = (),
    regeneration_hits: tuple[RegenerationHit, ...] = (),
) -> StateMarkedGeneratedEventBridge:
    """Promote state evidence to marking_good without changing causal mass."""
    if pair_cell < 0:
        raise ValueError("nonnegative pair-cell index required")
    if not isinstance(energy_generation, PhysicalHHGenerationGateWitness):
        raise TypeError("typed physical HH-generation gate witness required")
    if energy_generation.event_state_key != factorization.event_state_key:
        raise ValueError("physical-energy generation witness belongs to a different event state")
    if not factorization.state_mark.role_state_mark_available:
        raise ValueError("full-signed event role state did not pass the ordinary Christ gate")
    if factorization.state_mark.event_state_key != factorization.event_state_key:
        raise ValueError("factorization lost same-event state identity")
    if not factorization.bad_work_remains_terminal:
        raise ValueError("state-mark bridge may not reopen terminal bad positive work")
    if not registration_good and not physical_hits and not regeneration_hits:
        raise ValueError("failed common-slice registration must retain its existing first-stop provenance")

    event = GeneratedPairEvent(
        mass=float(factorization.canonical_good_positive_work),
        pair_cell=int(pair_cell),
        marking_good=True,
        registration_good=bool(registration_good),
        physical_hits=tuple(physical_hits),
        regeneration_hits=tuple(regeneration_hits),
    )
    return StateMarkedGeneratedEventBridge(
        event_state_key=factorization.event_state_key,
        pair_cell=int(pair_cell),
        canonical_good_positive_work=float(factorization.canonical_good_positive_work),
        canonical_bad_positive_work=float(factorization.canonical_bad_positive_work),
        canonical_negative_work=float(factorization.canonical_negative_work),
        generated_event=event,
        energy_generation=energy_generation,
    )


def compose_state_marked_event_with_master(
    bridge: StateMarkedGeneratedEventBridge,
    *,
    pair_cells_upper: int | None = None,
    scaled_lifetime: float = 1.0,
) -> StateMarkedMasterComposition:
    master = compile_generated_pair_master_measure(
        events=(bridge.generated_event,),
        pair_cells_upper=pair_cells_upper,
        scaled_lifetime=scaled_lifetime,
    )
    good = bridge.canonical_good_positive_work
    bad = bridge.canonical_bad_positive_work
    return StateMarkedMasterComposition(
        bridge=bridge,
        generated_master=master,
        terminal_bad_positive_work=bad,
        total_canonical_positive_work=good + bad,
        total_accounted_positive_work=master.total_mass + bad,
    )


def coalesce_state_marked_parent_slots(
    slot_weights: dict[str, float],
    parent_to_anchor: dict[str, tuple[int, ...]],
) -> dict[tuple[int, ...], float]:
    return pushforward_parent_slot_weights(slot_weights, parent_to_anchor)


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "input": "same-event GoodCausalStateMarkFactorization plus a replayed certified physical-energy HH-generation gate on the same event state",
        "generation_guard": "GeneratedPairEvent is forbidden unless route_physical_energy_causality returns physical_high_high_transfer_generation and the actual positive HH law meets its certified lower bound",
        "promotion": "role-state mark availability sets marking_good=True only on the already-canonical geometry-good dW+ restriction",
        "mass_identity": "GeneratedPairEvent.mass equals canonical good dW+ exactly; bad dW+ and dW- are excluded from generated causal mass",
        "compiler": "the existing recursive_physical_witness_constructor is used unchanged; registration and first-stop provenance remain independent downstream facts",
        "master_recombination": "existing generated-master output on good work plus already-terminal bad work reconstructs total canonical positive dW+ exactly once",
        "reservation": "no reserved Young/Christ gate is required for marking_good; the stronger counterfactual reservation certificate remains optional analysis",
        "root_reuse": "state-derived parent anchors use the existing Bargmann positive pushforward; coincident anchors coalesce instead of cloning parent energy identity",
        "negative": "canonical dW- never enters GeneratedPairEvent.mass and is never a payment or recursive owner created by this bridge",
        "scope": "does not infer generation from Young-good work alone; it also does not prove common-slice registration, generic HH termination, mixed-owner termination, singular-time closure, or global regularity",
    }
