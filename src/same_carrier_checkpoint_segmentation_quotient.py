from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    RESIDUAL_FRACTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)
from src.full_natural_checkpoint_quotient import FULL_NATURAL_CHECKPOINT
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION
from src.smooth_sgs_first_hit_extraction import (
    PhysicalPathMonitor,
    ThresholdTopology,
    first_physical_corridor_exit,
)


STATUS = (
    "EXACT_SAME_CARRIER_CHECKPOINT_SEGMENTATION_QUOTIENT__"
    "NATURAL_HORIZONS_DO_NOT_RESET_FIRST_HIT__"
    "CUMULATIVE_NATIVE_MONITORS_FROM_ONE_PHYSICAL_EVENT__"
    "INTERIOR_CHECKPOINT_ZENO_IS_STOP_OR_CONTINUATION__"
    "HARDEN_ONLY_AT_A_NEW_PHYSICAL_EVENT"
)

SAME_CARRIER_CONTINUATION = "same_event_anchored_smooth_carrier_continuation"


@dataclass(frozen=True)
class SameCarrierMonitorSegment:
    """One analysis segment of cumulative monitors for one fixed smooth carrier.

    The segment boundaries may be natural-horizon checkpoints, plotting times, or
    any other observer-chosen subdivision.  The monitor values are *not* increments
    local to this segment.  They are the cumulative observables from the same
    terminal physical event and the same terminal carrier/dual:

      K(s,t), |I_R(s,t)|, |I_HH(s,t)|.

    Only K is monotone by construction.  The two impulse magnitudes may decrease
    because the underlying complex cumulative impulses can cancel.  They must not
    be added segment-by-segment or converted into work.
    """

    carrier_id: str
    terminal_amplitude: float
    elapsed_times: tuple[float, ...]
    strain_action: tuple[float, ...]
    residual_impulse_abs: tuple[float, ...]
    hh_impulse_abs: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.carrier_id:
            raise ValueError("same-carrier segment requires a nonempty carrier id")
        amp = float(self.terminal_amplitude)
        if amp <= 0 or not math.isfinite(amp):
            raise ValueError("positive finite fixed terminal amplitude required")
        rows = (
            tuple(float(x) for x in self.elapsed_times),
            tuple(float(x) for x in self.strain_action),
            tuple(float(x) for x in self.residual_impulse_abs),
            tuple(float(x) for x in self.hh_impulse_abs),
        )
        n = len(rows[0])
        if n < 2 or any(len(x) != n for x in rows):
            raise ValueError("matching monitor paths of length at least two required")
        if any(not math.isfinite(x) for row in rows for x in row):
            raise ValueError("finite same-carrier monitor data required")
        if any(rows[0][j + 1] <= rows[0][j] for j in range(n - 1)):
            raise ValueError("elapsed times must increase strictly inside each segment")
        if any(x < 0 for row in rows[1:] for x in row):
            raise ValueError("native cumulative monitor magnitudes/actions are nonnegative")
        if any(rows[1][j + 1] + 2e-14 < rows[1][j] for j in range(n - 1)):
            raise ValueError("cumulative strain action cannot decrease")


def _scale_tol(*values: float) -> float:
    return 2e-11 * max(1.0, *(abs(float(x)) for x in values))


def join_same_carrier_segments(segments: Sequence[SameCarrierMonitorSegment]) -> dict[str, object]:
    """Remove checkpoint segmentation without changing the cumulative path.

    Adjacent segments must meet with the same carrier id, same terminal amplitude,
    same elapsed time, and the same cumulative monitor values.  A reset to zero at
    a checkpoint is rejected.  The coefficient-impulse magnitudes are only required
    to be continuous across a boundary; they are deliberately not required to be
    monotone or additive.
    """
    segs = tuple(segments)
    if not segs:
        raise ValueError("nonempty same-carrier segment family required")
    first = segs[0]
    amp = float(first.terminal_amplitude)
    if abs(first.elapsed_times[0]) > _scale_tol(first.elapsed_times[0]):
        raise ValueError("the first segment must start at the terminal event elapsed time zero")
    for name, value in (
        ("strain", first.strain_action[0]),
        ("role-interface impulse", first.residual_impulse_abs[0]),
        ("HH impulse", first.hh_impulse_abs[0]),
    ):
        if abs(value) > _scale_tol(value, amp):
            raise ValueError(f"cumulative {name} must start at zero at the terminal event")

    elapsed = list(first.elapsed_times)
    strain = list(first.strain_action)
    residual = list(first.residual_impulse_abs)
    hh = list(first.hh_impulse_abs)

    for left, right in zip(segs, segs[1:]):
        if right.carrier_id != first.carrier_id or left.carrier_id != first.carrier_id:
            raise TypeError("a no-event checkpoint cannot replace the event-anchored smooth carrier")
        if abs(float(right.terminal_amplitude) - amp) > _scale_tol(right.terminal_amplitude, amp):
            raise TypeError("a no-event checkpoint cannot reset the fixed terminal coefficient amplitude")
        boundary_values = (
            (left.elapsed_times[-1], right.elapsed_times[0], "elapsed time"),
            (left.strain_action[-1], right.strain_action[0], "cumulative strain action"),
            (left.residual_impulse_abs[-1], right.residual_impulse_abs[0], "cumulative role-interface impulse magnitude"),
            (left.hh_impulse_abs[-1], right.hh_impulse_abs[0], "cumulative HH impulse magnitude"),
        )
        for a, b, name in boundary_values:
            if abs(float(a) - float(b)) > _scale_tol(a, b, amp):
                raise TypeError(f"checkpoint reset/discontinuity detected in {name}")
        elapsed.extend(right.elapsed_times[1:])
        strain.extend(right.strain_action[1:])
        residual.extend(right.residual_impulse_abs[1:])
        hh.extend(right.hh_impulse_abs[1:])

    if any(strain[j + 1] + _scale_tol(strain[j], strain[j + 1]) < strain[j] for j in range(len(strain) - 1)):
        raise ValueError("joined cumulative strain action cannot decrease")
    return {
        "carrier_id": first.carrier_id,
        "terminal_amplitude": amp,
        "elapsed_times": tuple(elapsed),
        "strain_action": tuple(strain),
        "residual_impulse_abs": tuple(residual),
        "hh_impulse_abs": tuple(hh),
        "analysis_segments": len(segs),
        "inserted_checkpoint_boundaries": len(segs) - 1,
        "carrier_restarts": 0,
        "monitor_resets": 0,
    }


def same_carrier_first_exit(
    segments: Sequence[SameCarrierMonitorSegment],
    *,
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """First stop of one fixed carrier, independent of inserted checkpoint cuts."""
    path = join_same_carrier_segments(segments)
    amp = float(path["terminal_amplitude"])
    monitors = (
        PhysicalPathMonitor(
            "high_strain_critical_dissipation",
            LOW_STRAIN_ACTION,
            tuple(float(x) for x in path["strain_action"]),
            ThresholdTopology.CLOSED,
        ),
        PhysicalPathMonitor(
            ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
            RESIDUAL_FRACTION * amp,
            tuple(float(x) for x in path["residual_impulse_abs"]),
            ThresholdTopology.CLOSED,
        ),
        PhysicalPathMonitor(
            HH_COEFFICIENT_OBSTRUCTION,
            GENERATED_FRACTION * amp,
            tuple(float(x) for x in path["hh_impulse_abs"]),
            ThresholdTopology.CLOSED,
        ),
    )
    out = first_physical_corridor_exit(
        tuple(float(x) for x in path["elapsed_times"]),
        monitors,
        tie_tolerance=tie_tolerance,
    )
    stops = out.joint_first_stops
    needs_reentry = any(
        x in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for x in stops
    )
    return {
        "classification": "same_carrier_named_first_stop" if out.first_time is not None else "same_carrier_no_hit_continuation",
        "carrier_id": path["carrier_id"],
        "terminal_amplitude": amp,
        "first_elapsed": out.first_time,
        "joint_first_stops": stops,
        "individual_debuts": out.individual_debuts,
        "observed_elapsed_end": float(path["elapsed_times"][-1]),
        "analysis_segments": int(path["analysis_segments"]),
        "inserted_checkpoint_boundaries": int(path["inserted_checkpoint_boundaries"]),
        "carrier_restarts": 0,
        "monitor_resets": 0,
        "requires_physical_energy_reentry": needs_reentry,
        "coefficient_impulses_used_as_work": False,
        "checkpoint_segmentation_used_as_causal_order": False,
    }


def partition_same_carrier_path(
    *,
    carrier_id: str,
    terminal_amplitude: float,
    elapsed_times: Sequence[float],
    strain_action: Sequence[float],
    residual_impulse_abs: Sequence[float],
    hh_impulse_abs: Sequence[float],
    checkpoint_indices: Sequence[int],
) -> tuple[SameCarrierMonitorSegment, ...]:
    """Insert observer-chosen checkpoint cuts into one already-defined cumulative path."""
    rows = tuple(tuple(float(x) for x in row) for row in (elapsed_times, strain_action, residual_impulse_abs, hh_impulse_abs))
    n = len(rows[0])
    if n < 2 or any(len(row) != n for row in rows):
        raise ValueError("matching global cumulative monitor paths required")
    cuts = tuple(int(x) for x in checkpoint_indices)
    if tuple(sorted(set(cuts))) != cuts or any(x <= 0 or x >= n - 1 for x in cuts):
        raise ValueError("checkpoint indices must be unique increasing interior sample indices")
    starts = (0,) + cuts
    ends = cuts + (n - 1,)
    return tuple(
        SameCarrierMonitorSegment(
            carrier_id=carrier_id,
            terminal_amplitude=float(terminal_amplitude),
            elapsed_times=rows[0][a : b + 1],
            strain_action=rows[1][a : b + 1],
            residual_impulse_abs=rows[2][a : b + 1],
            hh_impulse_abs=rows[3][a : b + 1],
        )
        for a, b in zip(starts, ends)
    )


def segmentation_invariance(
    *,
    carrier_id: str,
    terminal_amplitude: float,
    elapsed_times: Sequence[float],
    strain_action: Sequence[float],
    residual_impulse_abs: Sequence[float],
    hh_impulse_abs: Sequence[float],
    checkpoint_indices: Sequence[int],
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Compare no checkpoint against any finite checkpoint segmentation of the same path."""
    whole = partition_same_carrier_path(
        carrier_id=carrier_id,
        terminal_amplitude=terminal_amplitude,
        elapsed_times=elapsed_times,
        strain_action=strain_action,
        residual_impulse_abs=residual_impulse_abs,
        hh_impulse_abs=hh_impulse_abs,
        checkpoint_indices=(),
    )
    segmented = partition_same_carrier_path(
        carrier_id=carrier_id,
        terminal_amplitude=terminal_amplitude,
        elapsed_times=elapsed_times,
        strain_action=strain_action,
        residual_impulse_abs=residual_impulse_abs,
        hh_impulse_abs=hh_impulse_abs,
        checkpoint_indices=checkpoint_indices,
    )
    a = same_carrier_first_exit(whole, tie_tolerance=tie_tolerance)
    b = same_carrier_first_exit(segmented, tie_tolerance=tie_tolerance)
    ta = a["first_elapsed"]
    tb = b["first_elapsed"]
    if (ta is None) != (tb is None):
        raise AssertionError("checkpoint insertion changed whether the fixed carrier has a first stop")
    residual = 0.0 if ta is None else abs(float(ta) - float(tb))
    if residual > _scale_tol(float(ta or 0.0), float(tb or 0.0), float(elapsed_times[-1])):
        raise AssertionError("checkpoint insertion moved the fixed-carrier first-stop time")
    if tuple(a["joint_first_stops"]) != tuple(b["joint_first_stops"]):
        raise AssertionError("checkpoint insertion changed the fixed-carrier joint first-stop set")
    return {
        "first_elapsed": ta,
        "joint_first_stops": tuple(a["joint_first_stops"]),
        "first_time_residual": residual,
        "checkpoint_count": len(tuple(checkpoint_indices)),
        "carrier_restarts": 0,
        "monitor_resets": 0,
        "segmentation_changed_first_hit": False,
    }


def checkpoint_continuation_policy(
    checkpoint_record: dict[str, object],
    *,
    request_carrier_replacement: bool = False,
    request_terminal_amplitude_reset: bool = False,
    request_monitor_reset: bool = False,
) -> dict[str, object]:
    """Fail closed if an analysis checkpoint is used to restart the causal carrier.

    Hard-shell observations at the checkpoint remain legitimate state sidecars.
    They can become event-anchored hard roles only if a later physical interaction
    or another named physical first stop supplies that event semantics.
    """
    if str(checkpoint_record.get("checkpoint_kind", "")) != FULL_NATURAL_CHECKPOINT:
        raise ValueError("same-carrier checkpoint policy requires a full-natural analysis checkpoint")
    if bool(checkpoint_record.get("physical_event_created", True)):
        raise ValueError("checkpoint record unexpectedly claims a physical event")
    if bool(checkpoint_record.get("causal_charge_created", True)):
        raise ValueError("checkpoint record unexpectedly claims a causal charge")
    if int(checkpoint_record.get("recursion_edges_added", 1)) != 0:
        raise ValueError("checkpoint record unexpectedly adds recursive event depth")
    if request_carrier_replacement:
        raise TypeError("a no-event checkpoint cannot replace the event-anchored smooth carrier")
    if request_terminal_amplitude_reset:
        raise TypeError("a no-event checkpoint cannot reset the terminal coefficient baseline")
    if request_monitor_reset:
        raise TypeError("a no-event checkpoint cannot reset cumulative native first-hit monitors")
    return {
        "canonical_continuation": SAME_CARRIER_CONTINUATION,
        "hard_shell_checkpoint_witnesses": "state_sidecars_only",
        "carrier_replacement_authorized": False,
        "terminal_amplitude_reset_authorized": False,
        "monitor_reset_authorized": False,
        "checkpoint_scale_path_is_physical_lineage": False,
        "hardening_requires_new_physical_event": True,
    }


def maximal_same_carrier_outcome(
    segments: Sequence[SameCarrierMonitorSegment],
    *,
    event_time: float,
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Continue one carrier through arbitrary horizons until a true stop or t=0.

    Natural-horizon service estimates may be read on subintervals, but the exact
    carrier equation and cumulative first-hit filtration do not expire there.
    """
    t = float(event_time)
    if t <= 0 or not math.isfinite(t):
        raise ValueError("positive finite terminal event time required")
    out = same_carrier_first_exit(segments, tie_tolerance=tie_tolerance)
    observed = float(out["observed_elapsed_end"])
    if observed > t + _scale_tol(observed, t):
        raise ValueError("same-carrier observation extends before the initial boundary")
    if out["first_elapsed"] is not None:
        return {
            **out,
            "endpoint_time": t - float(out["first_elapsed"]),
            "hits_initial_boundary": False,
            "physical_event_created_by_checkpoint": False,
        }
    if observed >= t - _scale_tol(observed, t):
        return {
            **out,
            "classification": "absorbing_initial_boundary",
            "joint_first_stops": ("t=0",),
            "first_elapsed": t,
            "endpoint_time": 0.0,
            "hits_initial_boundary": True,
            "requires_physical_energy_reentry": False,
            "physical_event_created_by_checkpoint": False,
        }
    return {
        **out,
        "classification": "same_carrier_event_free_continuation",
        "endpoint_time": t - observed,
        "hits_initial_boundary": False,
        "remaining_backward_time": t - observed,
        "physical_event_created_by_checkpoint": False,
        "next_action": "continue_the_same_event_anchored_carrier_and_cumulative_monitors",
    }


def interior_checkpoint_accumulation_outcome(
    *,
    event_time: float,
    accumulation_time: float,
    terminal_amplitude: float,
    strain_action_limit: float,
    residual_impulse_abs_limit: float,
    hh_impulse_abs_limit: float,
) -> dict[str, object]:
    """Classify a limit of infinitely many inserted horizons for one fixed carrier.

    On a pre-singular smooth interval the cumulative strain action is continuous
    (indeed AC/monotone) and the complex coefficient impulses are AC, hence their
    magnitudes are continuous.  If a closed threshold is attained at an interior
    accumulation time, that time is the first-stop face (coefficient faces still
    require Q^2 physical-energy reentry).  If all faces remain strictly below,
    continuity gives an open continuation past the checkpoint accumulation.  If
    the accumulation is t=0, the initial boundary absorbs.
    """
    t = float(event_time)
    s = float(accumulation_time)
    amp = float(terminal_amplitude)
    K = float(strain_action_limit)
    IR = float(residual_impulse_abs_limit)
    IH = float(hh_impulse_abs_limit)
    if t <= 0 or s < 0 or s >= t or amp <= 0 or min(K, IR, IH) < 0:
        raise ValueError("valid positive event/amplitude, interior-or-zero accumulation time, and nonnegative limits required")
    if not all(math.isfinite(x) for x in (t, s, amp, K, IR, IH)):
        raise ValueError("finite accumulation data required")
    if s == 0.0:
        return {
            "classification": "absorbing_initial_boundary",
            "joint_first_stops": ("t=0",),
            "requires_physical_energy_reentry": False,
            "checkpoint_accumulation_is_obstruction": False,
            "same_carrier_extends_past_accumulation": False,
        }
    tol = _scale_tol(amp, K, IR, IH)
    hits: list[str] = []
    if K >= LOW_STRAIN_ACTION - tol:
        hits.append("high_strain_critical_dissipation")
    if IR >= RESIDUAL_FRACTION * amp - tol:
        hits.append(ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION)
    if IH >= GENERATED_FRACTION * amp - tol:
        hits.append(HH_COEFFICIENT_OBSTRUCTION)
    if hits:
        stops = tuple(sorted(hits))
        return {
            "classification": "first_stop_at_interior_checkpoint_accumulation",
            "joint_first_stops": stops,
            "requires_physical_energy_reentry": any(
                x in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
                for x in stops
            ),
            "coefficient_impulses_used_as_work": False,
            "checkpoint_accumulation_is_obstruction": False,
            "same_carrier_extends_past_accumulation": False,
        }
    return {
        "classification": "same_carrier_extends_across_interior_checkpoint_accumulation",
        "joint_first_stops": (),
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
        "checkpoint_accumulation_is_obstruction": False,
        "same_carrier_extends_past_accumulation": True,
        "continuation_reason": "all cumulative native monitors retain strict margin and are continuous/AC on the smooth pre-singular interval",
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "fixed_carrier": "between genuine physical events the canonical object is the same smooth event-anchored carrier Q and terminal dual; a theorem horizon does not create a replacement carrier",
        "cumulative_monitors": "K_A[s,t], |I_role-interface[s,t]| and |I_HH[s,t]| are always read cumulatively from the same terminal event; K is monotone while impulse magnitudes may cancel and are therefore never summed segment-by-segment",
        "segmentation": "inserting or deleting any finite set of natural-horizon checkpoints leaves the first physical stop time and complete joint first-stop set unchanged",
        "checkpoint_policy": "hard-shell witnesses exposed at a no-event checkpoint are legitimate state sidecars only; they cannot reset terminal amplitude, monitor baselines, or harden into a new causal carrier without a new physical event",
        "long_horizon": "the exact outer-role/adjoint equations do not expire after one natural service window; the same carrier continues while cumulative support/strain and coefficient obstruction faces remain unhit",
        "interior_accumulation": "on a smooth pre-singular interval an interior accumulation of analysis checkpoints is either a closed physical first-stop/energy-reentry face at the limit or is crossed by the same carrier; checkpoint Zeno count and checkpoint scale path are not PDE obstructions",
        "boundary": "if the maximal same-carrier no-hit continuation reaches t=0, the initial boundary absorbs",
        "high_tail_separation": "this does not deny physical UV dynamics: actual high-tail dissipation/work remains a separate event theorem; only UV scale motion manufactured by re-hardening at no-event checkpoints is removed from canonical lineage",
        "scope": "this closes natural-horizon segmentation as a continuation obstruction for a fixed event-anchored carrier; it does not telescope infinitely recurring genuine physical owner events and does not prove Navier-Stokes global regularity",
    }


@dataclass(frozen=True)
class SameCarrierStress:
    samples: int
    worst_segmentation_first_time_residual: float
    segmentation_failures: int
    reset_barrier_failures: int
    nonmonotone_impulse_paths: int
    accumulation_stop_cases: int
    accumulation_continue_cases: int
    maximum_checkpoint_count: int


def _engineered_path(rng: random.Random) -> tuple[float, tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...], int]:
    n = rng.randint(8, 30)
    times = tuple(j / (n - 1) for j in range(n))
    amp = math.exp(rng.uniform(-3.0, 3.0))
    mode = rng.randrange(5)

    # K is a cumulative integral and therefore monotone.
    if mode == 1:
        tau = times[rng.randint(2, n - 3)]
        K = tuple((LOW_STRAIN_ACTION / tau) * x for x in times)
    else:
        kend = rng.uniform(0.05, 0.85) * LOW_STRAIN_ACTION
        exponent = rng.uniform(0.7, 1.8)
        K = tuple(kend * (x ** exponent) for x in times)

    ir_th = RESIDUAL_FRACTION * amp
    hh_th = GENERATED_FRACTION * amp

    def safe_wave(th: float, phase: float) -> tuple[float, ...]:
        out = []
        for x in times:
            base = 0.42 + 0.25 * math.sin(7.0 * x + phase) + 0.12 * math.sin(17.0 * x + 0.3 * phase)
            out.append(max(0.0, min(0.82, base)) * th)
        out[0] = 0.0
        return tuple(out)

    IR = safe_wave(ir_th, rng.uniform(-math.pi, math.pi))
    IH = safe_wave(hh_th, rng.uniform(-math.pi, math.pi))
    if mode == 2:
        j = rng.randint(2, n - 3)
        tau = times[j]
        IR = tuple((ir_th / tau) * x for x in times)
    elif mode == 3:
        j = rng.randint(2, n - 3)
        tau = times[j]
        IH = tuple((hh_th / tau) * x for x in times)
    elif mode == 4:
        j = rng.randint(2, n - 3)
        tau = times[j]
        IR = tuple((ir_th / tau) * x for x in times)
        IH = tuple((hh_th / tau) * x for x in times)
    return amp, times, K, IR, IH, mode


def stress(samples: int = 50_000, seed: int = 20260811) -> SameCarrierStress:
    rng = random.Random(seed)
    worst = 0.0
    segmentation_fail = 0
    reset_fail = 0
    nonmono = 0
    stop_cases = 0
    continue_cases = 0
    max_checkpoints = 0

    for i in range(samples):
        amp, times, K, IR, IH, _mode = _engineered_path(rng)
        n = len(times)
        possible = list(range(1, n - 1))
        rng.shuffle(possible)
        count = rng.randint(0, min(len(possible), 12))
        cuts = tuple(sorted(possible[:count]))
        max_checkpoints = max(max_checkpoints, len(cuts))
        inv = segmentation_invariance(
            carrier_id=f"carrier-{i}",
            terminal_amplitude=amp,
            elapsed_times=times,
            strain_action=K,
            residual_impulse_abs=IR,
            hh_impulse_abs=IH,
            checkpoint_indices=cuts,
            tie_tolerance=5e-11,
        )
        residual = float(inv["first_time_residual"])
        worst = max(worst, residual)
        if residual > 5e-11 or bool(inv["segmentation_changed_first_hit"]):
            segmentation_fail += 1
            raise AssertionError("analysis checkpoint segmentation changed the same-carrier first hit")

        if any(IR[j + 1] < IR[j] for j in range(n - 1)) or any(IH[j + 1] < IH[j] for j in range(n - 1)):
            nonmono += 1

        # Deliberately manufacture a baseline reset at one inserted boundary.
        if cuts:
            segs = list(
                partition_same_carrier_path(
                    carrier_id=f"carrier-{i}",
                    terminal_amplitude=amp,
                    elapsed_times=times,
                    strain_action=K,
                    residual_impulse_abs=IR,
                    hh_impulse_abs=IH,
                    checkpoint_indices=cuts,
                )
            )
            b = 1
            bad = segs[b]
            segs[b] = SameCarrierMonitorSegment(
                carrier_id=bad.carrier_id,
                terminal_amplitude=bad.terminal_amplitude,
                elapsed_times=bad.elapsed_times,
                strain_action=(0.0,) + bad.strain_action[1:],
                residual_impulse_abs=(0.0,) + bad.residual_impulse_abs[1:],
                hh_impulse_abs=(0.0,) + bad.hh_impulse_abs[1:],
            )
            try:
                join_same_carrier_segments(segs)
            except TypeError:
                pass
            else:
                reset_fail += 1
                raise AssertionError("checkpoint baseline reset crossed the same-carrier type barrier")

        # Interior accumulation: either a face is attained or the fixed carrier crosses it.
        if i % 2 == 0:
            acc = interior_checkpoint_accumulation_outcome(
                event_time=2.0,
                accumulation_time=0.7,
                terminal_amplitude=amp,
                strain_action_limit=0.7 * LOW_STRAIN_ACTION,
                residual_impulse_abs_limit=0.7 * RESIDUAL_FRACTION * amp,
                hh_impulse_abs_limit=0.7 * GENERATED_FRACTION * amp,
            )
            if not bool(acc["same_carrier_extends_past_accumulation"]):
                raise AssertionError("strictly subthreshold interior checkpoint accumulation blocked the same carrier")
            continue_cases += 1
        else:
            acc = interior_checkpoint_accumulation_outcome(
                event_time=2.0,
                accumulation_time=0.7,
                terminal_amplitude=amp,
                strain_action_limit=LOW_STRAIN_ACTION,
                residual_impulse_abs_limit=0.4 * RESIDUAL_FRACTION * amp,
                hh_impulse_abs_limit=0.4 * GENERATED_FRACTION * amp,
            )
            if "high_strain_critical_dissipation" not in tuple(acc["joint_first_stops"]):
                raise AssertionError("closed strain face disappeared at checkpoint accumulation")
            stop_cases += 1

    return SameCarrierStress(samples, worst, segmentation_fail, reset_fail, nonmono, stop_cases, continue_cases, max_checkpoints)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-same-carrier-checkpoint-segmentation-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "same_carrier_checkpoint_segmentation_quotient.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    md = f"""# Same-carrier checkpoint segmentation quotient

Status: **{cert['status']}**.

One physical event anchors one smooth carrier and one terminal dual.  Natural service horizons may be inserted while following that carrier, but they do not reset the event-to-endpoint cumulative observables

`K_A[s,t]`, `|I_role-interface[s,t]|`, `|I_HH[s,t]|`.

`K_A` is a monotone physical strain action.  The coefficient-impulse magnitudes need not be monotone because the underlying complex cumulative impulses may cancel; they are therefore **not** summed per segment and never used as work.

Deleting or inserting checkpoint cuts leaves the same native first stop and exact joint tie set.  A no-event checkpoint cannot replace the carrier, reset the fixed terminal amplitude, or reset the cumulative monitor baselines.  Its hard-shell readings remain state sidecars until a new physical event actually hardens a role.

An interior accumulation of arbitrarily many checkpoint cuts is likewise not a PDE obstruction.  On a smooth pre-singular interval the native cumulative monitors are continuous/AC: if a closed face is attained at the limit, that is the existing first stop (with coefficient faces routed through physical-energy reentry); if every face retains strict margin, the same carrier continues across the accumulation time.  If the continuation reaches `t=0`, the initial boundary absorbs.

Stress: `{out.samples}` cumulative-path/segmentation states
- worst segmentation first-time residual: `{out.worst_segmentation_first_time_residual:.3e}`
- segmentation failures: `{out.segmentation_failures}`
- checkpoint reset-barrier failures: `{out.reset_barrier_failures}`
- sampled nonmonotone coefficient-impulse paths: `{out.nonmonotone_impulse_paths}`
- interior accumulation stop cases: `{out.accumulation_stop_cases}`
- interior accumulation continuation cases: `{out.accumulation_continue_cases}`
- maximum inserted checkpoint count sampled: `{out.maximum_checkpoint_count}`

This theorem removes natural-horizon segmentation and no-event re-hardening from the continuation topology of a fixed event-anchored carrier.  It does not remove genuine high-tail physics, does not telescope infinitely recurring physical owner events, and makes no Navier--Stokes global-regularity claim.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
