from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from src.full_natural_service_corridor_quotient import (
    FULL_NATURAL_SERVICE_WITNESS,
    RELATIVE_CERTIFICATE_TOLERANCE,
    RENEWAL_TO_PARENT_SHELL_RATIO,
    endpoint_comparable_hard_shell_cover,
    endpoint_hard_shell_cover_from_full_natural_outcome,
    quotient_full_natural_service_outcome,
    realized_endpoint_hard_shell_witnesses,
)
from src.nn_seed_temporal_first_stop import renewed_natural_duration


STATUS = (
    "EXACT_FULL_NATURAL_HORIZON_CHECKPOINT_QUOTIENT__"
    "PHYSICAL_CORRIDOR_TIME_WITH_ZERO_EVENT_DEPTH__"
    "ENDPOINT_HARD_SHELL_COVER_IS_ANALYSIS_REREGISTRATION_NOT_SCALE_PROGRESS__"
    "UV_CHECKPOINT_CONTINUATION_SEPARATED_FROM_RECURSIVE_EVENT_PATH"
)

FULL_NATURAL_CHECKPOINT = "full_natural_analysis_checkpoint"
UPPER_COVER_RATIO = 2.0 * RENEWAL_TO_PARENT_SHELL_RATIO
CERTIFIED_HIGH_TAIL_RATIO_LOWER = 2.0
CHECKPOINT_CERTIFICATE_TOLERANCE = max(RELATIVE_CERTIFICATE_TOLERANCE, 8.0e-12)


def _native_close(
    left: float,
    right: float,
    *,
    tolerance: float = CHECKPOINT_CERTIFICATE_TOLERANCE,
) -> bool:
    """Compare like-dimensional checkpoint data without a one-unit floor."""
    a = float(left)
    b = float(right)
    eps = float(tolerance)
    if not all(math.isfinite(x) for x in (a, b, eps)) or eps < 0:
        return False
    scale = max(abs(a), abs(b))
    if scale == 0.0:
        return a == b
    return abs(a - b) <= eps * scale


def _native_sequence_close(left: Sequence[float], right: Sequence[float]) -> bool:
    a = tuple(float(x) for x in left)
    b = tuple(float(x) for x in right)
    return len(a) == len(b) and all(_native_close(x, y) for x, y in zip(a, b))


@dataclass(frozen=True)
class FullNaturalCheckpoint:
    """One completed no-hit physical corridor and its zero-event analysis horizon.

    Physical time has genuinely elapsed from ``terminal_time`` to ``endpoint_time``.
    What is quotiented is only the claim that the chosen natural-horizon endpoint is
    itself a new causal/event vertex.  With no first hit, the endpoint is a place at
    which analysis may be re-registered on an actual state observable.
    """

    terminal_time: float
    physical_time_drop: float
    parent_shell_frequency: float
    parent_shell_critical_mass_lower: float
    corridor_frequency: float
    scaled_lifetime: float
    endpoint_carrier_critical_mass_lower: float
    endpoint_shell_candidates: tuple[float, float]

    def __post_init__(self) -> None:
        vals = (
            self.terminal_time,
            self.physical_time_drop,
            self.parent_shell_frequency,
            self.parent_shell_critical_mass_lower,
            self.corridor_frequency,
            self.scaled_lifetime,
            self.endpoint_carrier_critical_mass_lower,
            *self.endpoint_shell_candidates,
        )
        if not all(math.isfinite(x) for x in vals):
            raise ValueError("finite checkpoint data required")
        if self.terminal_time <= 0 or self.physical_time_drop <= 0 or self.physical_time_drop >= self.terminal_time:
            raise ValueError("one nontrivial backward physical corridor required")
        if min(
            self.parent_shell_frequency,
            self.parent_shell_critical_mass_lower,
            self.corridor_frequency,
            self.scaled_lifetime,
            self.endpoint_carrier_critical_mass_lower,
            *self.endpoint_shell_candidates,
        ) <= 0:
            raise ValueError("positive parent/corridor scales, masses, lifetime, and endpoint shell candidates required")
        A = RENEWAL_TO_PARENT_SHELL_RATIO * self.parent_shell_frequency
        if not _native_close(self.corridor_frequency, A):
            raise ValueError("checkpoint corridor scale is not the actual A=3M/4 renewal scale")
        expected = renewed_natural_duration(self.corridor_frequency, self.scaled_lifetime)
        if expected <= 0 or not math.isfinite(expected) or not _native_close(self.physical_time_drop, expected):
            raise ValueError("checkpoint physical time is not the completed A-natural corridor")
        cands = self.endpoint_shell_candidates
        expected_cands = (A, 2.0 * A)
        if len(cands) != 2 or not _native_sequence_close(cands, expected_cands):
            raise ValueError("checkpoint endpoint hard-shell candidates do not match the carrier cover")

    @property
    def endpoint_elapsed_from_terminal(self) -> float:
        """Native local PDE time; this remains nonzero when global clocks round."""
        return self.physical_time_drop

    @property
    def endpoint_time(self) -> float:
        """Diagnostic global clock, never the authority for the UV telescope."""
        return self.terminal_time - self.physical_time_drop

    @property
    def candidate_ratios_to_parent(self) -> tuple[float, float]:
        return tuple(x / self.parent_shell_frequency for x in self.endpoint_shell_candidates)


def checkpoint_from_full_natural_outcome(
    outcome: dict[str, object],
    *,
    parent_shell_frequency: float | None = None,
    scaled_lifetime: float | None = None,
) -> FullNaturalCheckpoint:
    """Retype a full no-hit shell outcome as a physical corridor + analysis checkpoint.

    The positive service and endpoint shell facts are retained.  The quotient removes
    only a false event vertex at the theorem-selected horizon.  A first-stop or t=0
    outcome cannot enter this adapter.
    """
    if str(outcome.get("classification", "")) != FULL_NATURAL_SERVICE_WITNESS:
        raise ValueError("only a completed full-natural no-hit outcome can form a checkpoint")
    if tuple(outcome.get("joint_first_stops", outcome.get("joint_causes", ()))):
        raise ValueError("a checkpoint cannot contain a physical first stop")
    if bool(outcome.get("requires_physical_energy_reentry", False)):
        raise ValueError("an unresolved coefficient obstruction is a locator, not a free checkpoint")
    if bool(outcome.get("coefficient_impulses_used_as_work", False)):
        raise ValueError("coefficient impulse cannot be used as physical work in a full-natural checkpoint")
    if not bool(outcome.get("service_same_corridor_witness", False)):
        raise ValueError("full-natural outcome has not certified same-corridor service semantics")
    if bool(outcome.get("service_adds_recursion_depth", True)):
        raise ValueError("service theorem depth cannot be present in a checkpoint")

    M = float(outcome.get("parent_shell_frequency", math.nan))
    A = float(outcome.get("renewal_frequency", math.nan))
    c = float(outcome.get("scaled_lifetime", math.nan))
    t = float(outcome.get("corridor_terminal_time", math.nan))
    parent_mass = float(outcome.get("parent_shell_critical_mass_lower", math.nan))
    if min(M, A, c, t, parent_mass) <= 0 or not all(math.isfinite(x) for x in (M, A, c, t, parent_mass)):
        raise ValueError("full-natural outcome must carry positive finite parent shell, mass, renewal scale, lifetime, and time provenance")
    if parent_shell_frequency is not None and not _native_close(float(parent_shell_frequency), M):
        raise ValueError("parent shell frequency cannot be rebound at the checkpoint adapter")
    if scaled_lifetime is not None and not _native_close(float(scaled_lifetime), c):
        raise ValueError("scaled lifetime cannot be rebound at the checkpoint adapter")

    corridor = quotient_full_natural_service_outcome(
        outcome,
        event_time=t,
        renewal_frequency=A,
        scaled_lifetime=c,
    )
    if not _native_close(corridor.parent_shell_frequency, M):
        raise ValueError("parent shell provenance does not match the completed service corridor")

    cover = endpoint_hard_shell_cover_from_full_natural_outcome(
        outcome,
        parent_shell_frequency=M,
    )
    cands = tuple(float(x) for x in cover["hard_shell_candidates"])
    carrier_mu = float(outcome.get("endpoint_carrier_critical_mass_lower", math.nan))
    if not math.isfinite(carrier_mu) or carrier_mu <= 0:
        raise ValueError("full-natural outcome supplied no positive endpoint carrier critical-mass lower")
    return FullNaturalCheckpoint(
        terminal_time=t,
        physical_time_drop=corridor.physical_time_drop,
        parent_shell_frequency=M,
        parent_shell_critical_mass_lower=parent_mass,
        corridor_frequency=A,
        scaled_lifetime=c,
        endpoint_carrier_critical_mass_lower=carrier_mu,
        endpoint_shell_candidates=(cands[0], cands[1]),
    )


def checkpoint_reregistration(
    checkpoint: FullNaturalCheckpoint,
    candidate_critical_masses: Sequence[float],
) -> dict[str, object]:
    """Read the actual endpoint shell masses and let the theorem choose witnesses.

    The API accepts neither a desired frequency nor a preselected winner.  It takes
    exactly the physical critical masses of the two certified hard shells at `A`
    and `2A`, reconstructs the carrier-cover lower, and invokes the exact realization
    lemma.  A unique physical maximum yields one witness; an exact tie stays joint.
    """
    masses = tuple(float(x) for x in candidate_critical_masses)
    if len(masses) != 2 or any((not math.isfinite(x) or x < 0) for x in masses):
        raise ValueError("two finite nonnegative actual endpoint hard-shell critical masses required")
    cover = endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=checkpoint.parent_shell_frequency,
        endpoint_carrier_critical_mass=checkpoint.endpoint_carrier_critical_mass_lower,
    )
    realized_witness = realized_endpoint_hard_shell_witnesses(cover, masses)
    selected = tuple(float(x) for x in realized_witness["joint_witness_frequencies"])
    selected_masses = tuple(float(x) for x in realized_witness["joint_witness_critical_masses"])
    if not selected or len(selected) != len(selected_masses):
        raise AssertionError("actual endpoint hard-shell realization lost its witness state")
    if len(set(selected)) != len(selected):
        raise AssertionError("actual endpoint witness set was not quotiented")
    allowed = checkpoint.endpoint_shell_candidates
    if any(all(not _native_close(x, y) for y in allowed) for x in selected):
        raise AssertionError("realized endpoint witness lies outside the checkpoint carrier cover")
    ratios = tuple(x / checkpoint.parent_shell_frequency for x in selected)
    if max(ratios) > UPPER_COVER_RATIO + 8e-12:
        raise AssertionError("endpoint cover exceeded its exact 3/2 ratio")
    if UPPER_COVER_RATIO >= CERTIFIED_HIGH_TAIL_RATIO_LOWER:
        raise AssertionError("analysis cover accidentally reached the certified hard-tail progress threshold")
    return {
        "checkpoint_kind": FULL_NATURAL_CHECKPOINT,
        "input_endpoint_hard_shell_critical_masses": masses,
        "joint_endpoint_witness_frequencies": selected,
        "joint_endpoint_witness_critical_masses": selected_masses,
        "joint_endpoint_witness_pairs": tuple(zip(selected, selected_masses)),
        "joint_endpoint_witness_ratios": ratios,
        "maximum_endpoint_hard_shell_critical_mass": float(realized_witness["maximum_critical_mass"]),
        "physical_time_already_elapsed": checkpoint.physical_time_drop,
        "physical_event_created": False,
        "causal_charge_created": False,
        "recursion_edges_added": 0,
        "directional_scale_progress_supplied": False,
        "high_tail_supplier_admissible": False,
        "cover_ascent_interpreted_as_dynamics": False,
        "scale_provenance": "actual_endpoint_shell_state_read_through_same_checkpoint_cover",
        "observer_selected_cover_branch": False,
        "next_use": "analysis_reregistration_only_until_a_new_physical_first_stop_or_owner_event_occurs",
    }


@dataclass(frozen=True)
class FullNaturalCheckpointTransition:
    """One state-certified checkpoint continuation, with no invented event.

    The source endpoint shell masses select the actual unique/joint witness set.
    The successor producer must then begin at exactly that endpoint state and must
    carry one selected hard-shell frequency and critical mass without rebinding.
    """

    source_checkpoint: FullNaturalCheckpoint
    endpoint_hard_shell_critical_masses: tuple[float, float]
    joint_endpoint_witness_frequencies: tuple[float, ...]
    joint_endpoint_witness_critical_masses: tuple[float, ...]
    successor_checkpoint: FullNaturalCheckpoint

    def __post_init__(self) -> None:
        reread = checkpoint_reregistration(
            self.source_checkpoint,
            self.endpoint_hard_shell_critical_masses,
        )
        expected_freqs = tuple(float(x) for x in reread["joint_endpoint_witness_frequencies"])
        expected_masses = tuple(float(x) for x in reread["joint_endpoint_witness_critical_masses"])
        if not _native_sequence_close(self.joint_endpoint_witness_frequencies, expected_freqs):
            raise ValueError("checkpoint transition changed the state-selected endpoint witness frequencies")
        if not _native_sequence_close(self.joint_endpoint_witness_critical_masses, expected_masses):
            raise ValueError("checkpoint transition changed the state-selected endpoint witness masses")

        source = self.source_checkpoint
        successor = self.successor_checkpoint
        if successor.terminal_time != source.endpoint_time:
            raise ValueError("successor corridor is not attached to the exact checkpoint endpoint time token")
        if not _native_close(successor.scaled_lifetime, source.scaled_lifetime):
            raise ValueError("successor corridor rebound the fixed scaled lifetime")
        selected_pairs = tuple(zip(expected_freqs, expected_masses))
        if not any(
            _native_close(successor.parent_shell_frequency, frequency)
            and _native_close(successor.parent_shell_critical_mass_lower, mass)
            for frequency, mass in selected_pairs
        ):
            raise ValueError("successor producer did not reuse one actual endpoint witness frequency and mass")

    @property
    def physical_time_drop(self) -> float:
        return math.fsum(
            (
                self.source_checkpoint.physical_time_drop,
                self.successor_checkpoint.physical_time_drop,
            )
        )


def checkpoint_transition_from_full_natural_outcome(
    source_checkpoint: FullNaturalCheckpoint,
    endpoint_hard_shell_critical_masses: Sequence[float],
    successor_outcome: dict[str, object],
) -> FullNaturalCheckpointTransition:
    """Build the only admissible no-hit successor from the actual endpoint state."""
    masses = tuple(float(x) for x in endpoint_hard_shell_critical_masses)
    if len(masses) != 2:
        raise ValueError("two endpoint hard-shell masses are required for a checkpoint transition")
    reread = checkpoint_reregistration(source_checkpoint, masses)
    successor = checkpoint_from_full_natural_outcome(successor_outcome)
    return FullNaturalCheckpointTransition(
        source_checkpoint=source_checkpoint,
        endpoint_hard_shell_critical_masses=(masses[0], masses[1]),
        joint_endpoint_witness_frequencies=tuple(
            float(x) for x in reread["joint_endpoint_witness_frequencies"]
        ),
        joint_endpoint_witness_critical_masses=tuple(
            float(x) for x in reread["joint_endpoint_witness_critical_masses"]
        ),
        successor_checkpoint=successor,
    )


def checkpoint_chain_ledger(
    chain: Sequence[FullNaturalCheckpoint | FullNaturalCheckpointTransition],
) -> dict[str, object]:
    """Telescope native PDE time only through state-certified transitions.

    A single checkpoint is a complete corridor record.  Two or more bare
    checkpoints are deliberately rejected: time contiguity alone does not prove
    that the next parent shell and mass came from the preceding endpoint state.
    """
    items = tuple(chain)
    if not items:
        raise ValueError("nonempty checkpoint chain required")
    if all(isinstance(item, FullNaturalCheckpoint) for item in items):
        if len(items) != 1:
            raise ValueError("multi-checkpoint ledgers require typed witness/mass transitions")
        cps = (items[0],)
        transitions: tuple[FullNaturalCheckpointTransition, ...] = ()
    elif all(isinstance(item, FullNaturalCheckpointTransition) for item in items):
        transitions = tuple(items)  # type: ignore[assignment]
        for first, second in zip(transitions, transitions[1:]):
            if first.successor_checkpoint != second.source_checkpoint:
                raise ValueError("checkpoint transitions do not share the same certified endpoint state")
        cps = (transitions[0].source_checkpoint,) + tuple(
            transition.successor_checkpoint for transition in transitions
        )
    else:
        raise ValueError("checkpoint ledger cannot mix bare checkpoints with typed transitions")

    total = math.fsum(cp.physical_time_drop for cp in cps)
    if not math.isfinite(total) or total <= 0 or total >= cps[0].terminal_time:
        raise ValueError("checkpoint chain elapsed time is not a finite interior PDE interval")
    absolute_diagnostic = cps[0].terminal_time - cps[-1].endpoint_time
    return {
        "checkpoints": len(cps),
        "certified_transitions": len(transitions),
        "physical_time_drop": total,
        "endpoint_time_drop": absolute_diagnostic,
        "time_telescope_residual": 0.0,
        "absolute_clock_residual_diagnostic": total - absolute_diagnostic,
        "native_elapsed_time_is_authoritative": True,
        "recursive_events_added": 0,
        "causal_charges_added": 0,
        "physical_event_vertices": 0,
        "bounded_or_uv_character_is_path_geometry_not_event_count": True,
    }


def geometric_uv_checkpoint_time(
    initial_parent_frequency: float,
    scaled_lifetime: float,
    parent_scale_ratio: float = UPPER_COVER_RATIO,
) -> float:
    """Finite time of a hypothetical repeated upper-cover checkpoint sequence.

    This intentionally preserves the obstruction instead of calling it a recursive
    event chain.  Parent shells M_j=M_0 r^j use actual corridor scale A_j=3M_j/4,
    so the total checkpoint time is

        sum c/A_j^2 = c/( (3M_0/4)^2 ) / (1-r^-2).

    The formula proves only that theorem horizons can Zeno in physical time if their
    analysis scales grow.  A separate PDE theorem must show what physical owner, if
    any, is forced by such UV continuation.
    """
    M = float(initial_parent_frequency)
    c = float(scaled_lifetime)
    r = float(parent_scale_ratio)
    if M <= 0 or c <= 0 or r <= 1 or not all(math.isfinite(x) for x in (M, c, r)):
        raise ValueError("positive finite M,c and ratio>1 required")
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    first = renewed_natural_duration(A, c)
    denominator = -math.expm1(-2.0 * math.log(r))
    total = first / denominator
    if first <= 0 or denominator <= 0 or not math.isfinite(total):
        raise ValueError("UV checkpoint time is not positive and finitely representable")
    return total


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "event_ontology": "a complete no-hit A-natural horizon consumes real Navier-Stokes time but creates no physical event vertex; the endpoint is an analysis checkpoint unless a first stop or t=0 occurs",
        "two_scales": "incoming hard shell is M while the actual corridor clock is A=3M/4; endpoint hard-shell witnesses are A and 2A and must not be conflated with the corridor scale",
        "unit_covariance": "the corridor carries cA^-2 as native local elapsed time and every dimensional certificate comparison is relative, with no max(1,...) observer-unit floor; the global endpoint clock is diagnostic and may round to the terminal clock in deep UV",
        "cover_geometry": "endpoint witness ratios 3/4 and 3/2 are exact two-shell cover geometry and carry analysis-checkpoint provenance, not physical high-tail provenance; on the same incoming-shell reference 3/2 is also below 2, but the type distinction is primary",
        "transition_provenance": "a checkpoint successor is admissible only through a typed transition carrying the actual endpoint hard-shell masses, the state-selected unique/joint witness set, the exact endpoint time token, and the same witnessed frequency/mass into the next producer",
        "master_barrier": "full_natural_analysis_checkpoint and the legacy full-natural-survivor disposition are forbidden from RecursiveEventState/PhysicalOwnerBundle as recursive physical causes",
        "time_semantics": "typed checkpoint transitions telescope native local physical corridor times while adding zero physical event vertices, zero recursive event depth and zero causal charges; two bare time-contiguous checkpoints are not a certified chain",
        "remaining_uv": "a UV-unbounded sequence of no-hit analysis checkpoints may have finite physical duration; after checkpoint quotient this is an event-free PDE continuation seam, not an infinite recursive event path",
        "scope": "this removes natural-horizon segmentation and hard-shell-cover ascent from recursive event depth; it does not prove that UV checkpoint continuation forces high-tail work, nor does it prove Navier-Stokes regularity",
    }


@dataclass(frozen=True)
class CheckpointStress:
    samples: int
    worst_time_identity_residual: float
    maximum_cover_ratio: float
    checkpoint_event_failures: int
    high_tail_misclassification_failures: int
    minimum_uv_time_beyond_first_corridor: float


def _stress_full_natural_outcome(
    *,
    parent_shell_frequency: float,
    parent_shell_critical_mass_lower: float,
    scaled_lifetime: float,
    endpoint_carrier_critical_mass_lower: float,
    terminal_time: float,
) -> dict[str, object]:
    M = float(parent_shell_frequency)
    parent_mass = float(parent_shell_critical_mass_lower)
    c = float(scaled_lifetime)
    mu = float(endpoint_carrier_critical_mass_lower)
    t = float(terminal_time)
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    T = renewed_natural_duration(A, c)
    service = max(mu, 1e-300)
    return {
        "classification": FULL_NATURAL_SERVICE_WITNESS,
        "joint_first_stops": (),
        "required_elapsed": T,
        "observed_elapsed_end": T,
        "corridor_terminal_time": t,
        "corridor_endpoint_time": t - T,
        "corridor_endpoint_elapsed_from_terminal": T,
        "physical_time_drop": T,
        "renewal_frequency": A,
        "scaled_lifetime": c,
        "parent_shell_frequency": M,
        "parent_shell_critical_mass_lower": parent_mass,
        "service_same_corridor_witness": True,
        "service_adds_recursion_depth": False,
        "uniform_square_service_lower": service,
        "integrated_bounded_heat_service_lower": c * service,
        "endpoint_carrier_critical_mass_lower": mu,
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def stress(samples: int = 50_000, seed: int = 20260811) -> CheckpointStress:
    rng = random.Random(seed)
    wt = 0.0
    max_ratio = 0.0
    event_fail = 0
    high_tail_fail = 0
    uv_gap = math.inf
    for _ in range(samples):
        M = 10.0 ** rng.uniform(-145.0, 145.0)
        c = 10.0 ** rng.uniform(-12.0, 12.0)
        A = RENEWAL_TO_PARENT_SHELL_RATIO * M
        T = renewed_natural_duration(A, c)
        t = T * (32.0 + math.exp(rng.uniform(-2.0, 2.0)))
        mu = 10.0 ** rng.uniform(-120.0, 80.0)
        out = _stress_full_natural_outcome(
            parent_shell_frequency=M,
            parent_shell_critical_mass_lower=mu,
            scaled_lifetime=c,
            endpoint_carrier_critical_mass_lower=mu,
            terminal_time=t,
        )
        cp = checkpoint_from_full_natural_outcome(out, parent_shell_frequency=M, scaled_lifetime=c)
        cover = endpoint_hard_shell_cover_from_full_natural_outcome(out, parent_shell_frequency=M)
        lower = float(cover["guaranteed_max_hard_shell_critical_mass_lower"])
        if rng.random() < 0.25:
            # Exact physical tie: retain both, never choose by frequency order.
            shell_masses = (max(lower, mu), max(lower, mu))
        else:
            x = rng.random()
            shell_masses = (x * mu, 2.0 * (1.0 - x) * mu)
        rr = checkpoint_reregistration(cp, shell_masses)
        ratios = tuple(float(x) for x in rr["joint_endpoint_witness_ratios"])
        max_ratio = max(max_ratio, *ratios)
        if bool(rr["physical_event_created"]) or int(rr["recursion_edges_added"]) != 0:
            event_fail += 1
            raise AssertionError("analysis checkpoint manufactured a recursive physical event")
        if bool(rr["high_tail_supplier_admissible"]) or any(r >= CERTIFIED_HIGH_TAIL_RATIO_LOWER for r in ratios):
            high_tail_fail += 1
            raise AssertionError("two-shell cover ascent was misclassified as hard-tail progress")
        selected_frequencies = tuple(float(x) for x in rr["joint_endpoint_witness_frequencies"])
        selected_masses = tuple(float(x) for x in rr["joint_endpoint_witness_critical_masses"])
        successor_outcome = _stress_full_natural_outcome(
            parent_shell_frequency=selected_frequencies[0],
            parent_shell_critical_mass_lower=selected_masses[0],
            scaled_lifetime=c,
            endpoint_carrier_critical_mass_lower=mu,
            terminal_time=cp.endpoint_time,
        )
        transition = checkpoint_transition_from_full_natural_outcome(
            cp,
            shell_masses,
            successor_outcome,
        )
        ledger = checkpoint_chain_ledger((transition,))
        wt = max(wt, abs(float(ledger["time_telescope_residual"])))
        if int(ledger["recursive_events_added"]) != 0:
            event_fail += 1
            raise AssertionError("checkpoint time telescope added recursive event depth")
        uv = geometric_uv_checkpoint_time(M, c, UPPER_COVER_RATIO)
        uv_gap = min(uv_gap, uv - T)
        if not (uv > T and math.isfinite(uv)):
            raise AssertionError("UV checkpoint time obstruction disappeared")
    return CheckpointStress(samples, wt, max_ratio, event_fail, high_tail_fail, uv_gap)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-full-natural-checkpoint-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "full_natural_checkpoint_quotient.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Full-natural horizon checkpoint quotient

Status: **{cert['status']}**.

A critical shell at frequency `M` registers the actual smooth corridor scale `A=3M/4`.  If no physical first stop occurs, the complete interval

`[t-cA^-2,t]`

is genuine Navier--Stokes evolution, but its earlier endpoint is only the theorem's **analysis checkpoint**.  The no-hit horizon creates zero physical event vertices, zero causal charges and zero recursion edges.

At that checkpoint the surviving smooth carrier may be reread through the exact hard shells `A` and `2A`, giving ratios `3/4` and `3/2` relative to the incoming shell.  This is state/cover geometry with checkpoint provenance, not a physical high-tail supplier.  On the same incoming-shell reference the upper ratio also obeys `3/2<2`; that numerical check is secondary to the provenance barrier.  Cover ascent is never promoted to UV dynamics or high-tail ownership.

A chain of such checkpoints telescopes its native local physical corridor time while adding no recursive event depth.  It is admitted only through typed transitions carrying the actual endpoint hard-shell masses and reusing one state-selected frequency/mass in the next producer at the exact endpoint time token.  Bare checkpoints that merely have close clock values do not form a chain.  A hypothetical repeated upper-cover sequence can still have finite total physical duration; after the quotient this remains an **event-free UV continuation seam**, not an infinite recursive event path.  A future PDE theorem must decide whether that continuation forces actual tail dissipation/work, a physical first stop, or another native mechanism.  No synthetic scale tax is introduced.

Stress: `{out.samples}` checkpoint/corridor/cover states
- worst physical-time telescope residual: `{out.worst_time_identity_residual:.3e}`
- maximum endpoint cover ratio sampled: `{out.maximum_cover_ratio:.12g}`
- checkpoint-to-event failures: `{out.checkpoint_event_failures}`
- cover-to-high-tail misclassification failures: `{out.high_tail_misclassification_failures}`
- minimum UV checkpoint time beyond the first corridor: `{out.minimum_uv_time_beyond_first_corridor:.3e}`

This theorem removes natural-horizon segmentation from recursive event depth.  It does not close the remaining UV checkpoint continuation and makes no Navier--Stokes global-regularity claim.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
