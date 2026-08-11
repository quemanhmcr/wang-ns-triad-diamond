from __future__ import annotations

import argparse
import cmath
import json
import math
import random
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    RESIDUAL_FRACTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)
from src.full_natural_checkpoint_quotient import FullNaturalCheckpoint
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION
from src.smooth_sgs_first_hit_extraction import superlevel_debut_piecewise_linear


STATUS = (
    "EXACT_SAME_CARRIER_CHECKPOINT_SEGMENTATION_QUOTIENT__"
    "ONE_EVENT_CARRIER_DUAL_AND_PDE_PATH_PROVENANCE__"
    "CUMULATIVE_COMPLEX_NATIVE_MONITORS__"
    "OBSERVER_CUTS_ARE_NOT_NATURAL_WINDOWS__"
    "HARDEN_ONLY_AT_A_NEW_PHYSICAL_EVENT"
)

SAME_CARRIER_CONTINUATION = "same_event_anchored_smooth_carrier_continuation"
STRAIN_STOP = "high_strain_critical_dissipation"


def _finite_complex(value: complex) -> bool:
    z = complex(value)
    return math.isfinite(z.real) and math.isfinite(z.imag) and math.isfinite(abs(z))


def _native_duration(carrier_frequency: float, scaled_lifetime: float) -> float:
    """Evaluate c A^-2 without first forming A^2, which may overflow."""
    A = float(carrier_frequency)
    c = float(scaled_lifetime)
    if A <= 0 or c <= 0 or not math.isfinite(A) or not math.isfinite(c):
        raise ValueError("positive finite carrier frequency and scaled lifetime required")
    duration = (math.sqrt(c) / A) ** 2
    if duration <= 0 or not math.isfinite(duration):
        raise ValueError("fixed-carrier natural duration must be positive and representable in native time")
    return duration


@dataclass(frozen=True)
class SameCarrierProvenance:
    """Identity of one event-anchored carrier on one actual PDE trajectory.

    Equality is deliberately exact.  A carrier, terminal dual, trajectory, terminal
    state, scale, lifetime, or terminal coefficient cannot be rebound merely because
    the corresponding floating values are numerically close.
    """

    event_id: str
    carrier_id: str
    terminal_dual_id: str
    trajectory_id: str
    terminal_state_token: str
    terminal_time: float
    carrier_frequency: float
    scaled_lifetime: float
    terminal_coefficient: complex

    def __post_init__(self) -> None:
        labels = (
            self.event_id,
            self.carrier_id,
            self.terminal_dual_id,
            self.trajectory_id,
            self.terminal_state_token,
        )
        if any(not isinstance(x, str) or not x for x in labels):
            raise ValueError("nonempty physical event/carrier/dual/trajectory/state tokens required")
        t = float(self.terminal_time)
        A = float(self.carrier_frequency)
        c = float(self.scaled_lifetime)
        z = complex(self.terminal_coefficient)
        if min(t, A, c) <= 0 or not all(math.isfinite(x) for x in (t, A, c)):
            raise ValueError("positive finite terminal time, carrier frequency, and lifetime required")
        if not _finite_complex(z) or abs(z) <= 0:
            raise ValueError("nonzero finite terminal coefficient required")
        _native_duration(A, c)
        object.__setattr__(self, "terminal_time", t)
        object.__setattr__(self, "carrier_frequency", A)
        object.__setattr__(self, "scaled_lifetime", c)
        object.__setattr__(self, "terminal_coefficient", z)

    @property
    def terminal_amplitude(self) -> float:
        return abs(self.terminal_coefficient)

    @property
    def native_window_duration(self) -> float:
        return _native_duration(self.carrier_frequency, self.scaled_lifetime)


@dataclass(frozen=True)
class SameCarrierMonitorSegment:
    """A restriction of one cumulative monitor path to one observer interval.

    The complex impulses themselves are retained.  Their magnitudes are derived
    observables, not independently supplied paths.  Shared state tokens and exact
    cumulative boundary values are the gluing authority.
    """

    provenance: SameCarrierProvenance
    state_tokens: tuple[str, ...]
    elapsed_times: tuple[float, ...]
    strain_action: tuple[float, ...]
    residual_impulse: tuple[complex, ...]
    hh_impulse: tuple[complex, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.provenance, SameCarrierProvenance):
            raise TypeError("typed same-carrier provenance required")
        tokens = tuple(self.state_tokens)
        elapsed = tuple(float(x) for x in self.elapsed_times)
        strain = tuple(float(x) for x in self.strain_action)
        residual = tuple(complex(x) for x in self.residual_impulse)
        hh = tuple(complex(x) for x in self.hh_impulse)
        n = len(elapsed)
        if n < 2 or any(len(row) != n for row in (tokens, strain, residual, hh)):
            raise ValueError("matching same-carrier paths of length at least two required")
        if any(not isinstance(token, str) or not token for token in tokens):
            raise ValueError("nonempty PDE state token required at every path sample")
        if any(tokens[j + 1] == tokens[j] for j in range(n - 1)):
            raise ValueError("distinct elapsed samples require distinct PDE state tokens")
        if any(not math.isfinite(x) or x < 0 for x in elapsed):
            raise ValueError("finite nonnegative native elapsed times required")
        if any(elapsed[j + 1] <= elapsed[j] for j in range(n - 1)):
            raise ValueError("native elapsed times must increase strictly")
        if any(not math.isfinite(x) or x < 0 for x in strain):
            raise ValueError("finite nonnegative cumulative strain action required")
        if any(strain[j + 1] < strain[j] for j in range(n - 1)):
            raise ValueError("cumulative strain action cannot decrease")
        if any(not _finite_complex(x) for row in (residual, hh) for x in row):
            raise ValueError("finite cumulative complex coefficient impulses required")
        object.__setattr__(self, "state_tokens", tokens)
        object.__setattr__(self, "elapsed_times", elapsed)
        object.__setattr__(self, "strain_action", strain)
        object.__setattr__(self, "residual_impulse", residual)
        object.__setattr__(self, "hh_impulse", hh)

    @property
    def carrier_id(self) -> str:
        return self.provenance.carrier_id

    @property
    def terminal_amplitude(self) -> float:
        return self.provenance.terminal_amplitude

    @property
    def residual_impulse_abs(self) -> tuple[float, ...]:
        return tuple(abs(x) for x in self.residual_impulse)

    @property
    def hh_impulse_abs(self) -> tuple[float, ...]:
        return tuple(abs(x) for x in self.hh_impulse)


def join_same_carrier_segments(segments: Sequence[SameCarrierMonitorSegment]) -> dict[str, object]:
    """Glue exact restrictions of one event/carrier/dual/PDE path.

    Numerical closeness has no causal authority.  Every shared boundary must carry
    the same state token, native elapsed token, strain value, and complex impulses.
    """
    segs = tuple(segments)
    if not segs:
        raise ValueError("nonempty same-carrier segment family required")
    if any(not isinstance(segment, SameCarrierMonitorSegment) for segment in segs):
        raise TypeError("typed same-carrier monitor segments required")
    provenance = segs[0].provenance
    if any(segment.provenance != provenance for segment in segs):
        raise TypeError("checkpoint cannot replace or rebind event, carrier, terminal dual, trajectory, scale, or terminal coefficient")

    first = segs[0]
    if first.elapsed_times[0] != 0.0:
        raise ValueError("the cumulative path must begin at terminal event elapsed time zero")
    if first.state_tokens[0] != provenance.terminal_state_token:
        raise TypeError("the first PDE state token is not the event terminal state")
    if first.strain_action[0] != 0.0:
        raise ValueError("cumulative strain action must start at zero at the terminal event")
    if first.residual_impulse[0] != 0.0j or first.hh_impulse[0] != 0.0j:
        raise ValueError("cumulative coefficient impulses must start at zero at the terminal event")

    tokens = list(first.state_tokens)
    elapsed = list(first.elapsed_times)
    strain = list(first.strain_action)
    residual = list(first.residual_impulse)
    hh = list(first.hh_impulse)
    for left, right in zip(segs, segs[1:]):
        boundary_values = (
            (left.state_tokens[-1], right.state_tokens[0], "PDE state token"),
            (left.elapsed_times[-1], right.elapsed_times[0], "native elapsed time token"),
            (left.strain_action[-1], right.strain_action[0], "cumulative strain action"),
            (left.residual_impulse[-1], right.residual_impulse[0], "cumulative complex role-interface impulse"),
            (left.hh_impulse[-1], right.hh_impulse[0], "cumulative complex HH impulse"),
        )
        for left_value, right_value, name in boundary_values:
            if left_value != right_value:
                raise TypeError(f"checkpoint reset/discontinuity detected in {name}")
        tokens.extend(right.state_tokens[1:])
        elapsed.extend(right.elapsed_times[1:])
        strain.extend(right.strain_action[1:])
        residual.extend(right.residual_impulse[1:])
        hh.extend(right.hh_impulse[1:])

    if any(strain[j + 1] < strain[j] for j in range(len(strain) - 1)):
        raise ValueError("joined cumulative strain action cannot decrease")
    return {
        "provenance": provenance,
        "event_id": provenance.event_id,
        "carrier_id": provenance.carrier_id,
        "terminal_dual_id": provenance.terminal_dual_id,
        "trajectory_id": provenance.trajectory_id,
        "terminal_amplitude": provenance.terminal_amplitude,
        "state_tokens": tuple(tokens),
        "elapsed_times": tuple(elapsed),
        "strain_action": tuple(strain),
        "residual_impulse": tuple(residual),
        "hh_impulse": tuple(hh),
        "residual_impulse_abs": tuple(abs(x) for x in residual),
        "hh_impulse_abs": tuple(abs(x) for x in hh),
        "analysis_segments": len(segs),
        "inserted_checkpoint_boundaries": len(segs) - 1,
        "carrier_restarts": 0,
        "monitor_resets": 0,
        "complex_phase_preserved": True,
        "same_pde_trajectory_certified": True,
    }


def complex_modulus_debut_piecewise_linear(
    times: Sequence[float],
    values: Sequence[complex],
    threshold: float,
) -> float | None:
    """First contact |z|>=threshold for a piecewise-linear complex path.

    Interpolating endpoint magnitudes is wrong when phase changes.  This routine
    solves the quadratic circle intersection of the actual complex chord.
    """
    t = tuple(float(x) for x in times)
    z = tuple(complex(x) for x in values)
    level = float(threshold)
    if len(t) < 2 or len(t) != len(z):
        raise ValueError("matching complex path samples of length at least two required")
    if level <= 0 or not math.isfinite(level):
        raise ValueError("positive finite modulus threshold required")
    if any(not math.isfinite(x) for x in t) or any(not _finite_complex(x) for x in z):
        raise ValueError("finite complex path and times required")
    if any(t[j + 1] <= t[j] for j in range(len(t) - 1)):
        raise ValueError("strictly increasing native times required")
    if abs(z[0]) >= level:
        return t[0]

    for j in range(1, len(t)):
        left = z[j - 1]
        right = z[j]
        if abs(right) < level:
            continue
        direction = right - left
        scale = max(level, abs(left), abs(right))
        l = left / scale
        d = direction / scale
        radius = level / scale
        qa = d.real * d.real + d.imag * d.imag
        qb = 2.0 * (l.real * d.real + l.imag * d.imag)
        qc = l.real * l.real + l.imag * l.imag - radius * radius
        if qa == 0.0:
            return t[j]
        discriminant = max(0.0, qb * qb - 4.0 * qa * qc)
        root_disc = math.sqrt(discriminant)
        if qb >= 0.0 and qb + root_disc != 0.0:
            fraction = (-2.0 * qc) / (qb + root_disc)
        else:
            fraction = (-qb + root_disc) / (2.0 * qa)
        fraction = min(1.0, max(0.0, fraction))
        return t[j - 1] + fraction * (t[j] - t[j - 1])
    return None


def same_carrier_first_exit(
    segments: Sequence[SameCarrierMonitorSegment],
    *,
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Native first stop of one fixed carrier, independent of observer cuts."""
    path = join_same_carrier_segments(segments)
    elapsed = tuple(float(x) for x in path["elapsed_times"])
    amplitude = float(path["terminal_amplitude"])
    debuts = {
        STRAIN_STOP: superlevel_debut_piecewise_linear(
            elapsed,
            tuple(float(x) for x in path["strain_action"]),
            LOW_STRAIN_ACTION,
        ),
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION: complex_modulus_debut_piecewise_linear(
            elapsed,
            tuple(complex(x) for x in path["residual_impulse"]),
            RESIDUAL_FRACTION * amplitude,
        ),
        HH_COEFFICIENT_OBSTRUCTION: complex_modulus_debut_piecewise_linear(
            elapsed,
            tuple(complex(x) for x in path["hh_impulse"]),
            GENERATED_FRACTION * amplitude,
        ),
    }
    finite = tuple(value for value in debuts.values() if value is not None)
    first = min(finite) if finite else None
    tolerance = 0.0 if tie_tolerance is None else float(tie_tolerance)
    if tolerance < 0 or not math.isfinite(tolerance):
        raise ValueError("finite nonnegative regression tie tolerance required")
    stops = () if first is None else tuple(
        sorted(label for label, value in debuts.items() if value is not None and abs(value - first) <= tolerance)
    )
    needs_reentry = any(
        label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for label in stops
    )
    return {
        "classification": "same_carrier_named_first_stop" if first is not None else "same_carrier_no_hit_continuation",
        "event_id": path["event_id"],
        "carrier_id": path["carrier_id"],
        "terminal_dual_id": path["terminal_dual_id"],
        "trajectory_id": path["trajectory_id"],
        "terminal_amplitude": amplitude,
        "first_elapsed": first,
        "joint_first_stops": stops,
        "individual_debuts": debuts,
        "observed_elapsed_end": elapsed[-1],
        "endpoint_state_token": tuple(path["state_tokens"])[-1],
        "analysis_segments": int(path["analysis_segments"]),
        "inserted_checkpoint_boundaries": int(path["inserted_checkpoint_boundaries"]),
        "carrier_restarts": 0,
        "monitor_resets": 0,
        "requires_physical_energy_reentry": needs_reentry,
        "coefficient_impulses_used_as_work": False,
        "checkpoint_segmentation_used_as_causal_order": False,
        "complex_phase_preserved": True,
        "same_pde_trajectory_certified": True,
    }


def partition_same_carrier_path(
    *,
    provenance: SameCarrierProvenance,
    state_tokens: Sequence[str],
    elapsed_times: Sequence[float],
    strain_action: Sequence[float],
    residual_impulse: Sequence[complex],
    hh_impulse: Sequence[complex],
    checkpoint_indices: Sequence[int],
) -> tuple[SameCarrierMonitorSegment, ...]:
    """Restrict one already-certified cumulative PDE path at observer cuts."""
    tokens = tuple(state_tokens)
    elapsed = tuple(float(x) for x in elapsed_times)
    strain = tuple(float(x) for x in strain_action)
    residual = tuple(complex(x) for x in residual_impulse)
    hh = tuple(complex(x) for x in hh_impulse)
    n = len(elapsed)
    if n < 2 or any(len(row) != n for row in (tokens, strain, residual, hh)):
        raise ValueError("matching global cumulative same-carrier paths required")
    cuts = tuple(int(x) for x in checkpoint_indices)
    if tuple(sorted(set(cuts))) != cuts or any(x <= 0 or x >= n - 1 for x in cuts):
        raise ValueError("checkpoint indices must be unique increasing interior sample indices")
    starts = (0,) + cuts
    ends = cuts + (n - 1,)
    return tuple(
        SameCarrierMonitorSegment(
            provenance=provenance,
            state_tokens=tokens[a : b + 1],
            elapsed_times=elapsed[a : b + 1],
            strain_action=strain[a : b + 1],
            residual_impulse=residual[a : b + 1],
            hh_impulse=hh[a : b + 1],
        )
        for a, b in zip(starts, ends)
    )


def segmentation_invariance(
    *,
    provenance: SameCarrierProvenance,
    state_tokens: Sequence[str],
    elapsed_times: Sequence[float],
    strain_action: Sequence[float],
    residual_impulse: Sequence[complex],
    hh_impulse: Sequence[complex],
    checkpoint_indices: Sequence[int],
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Compare no cut with finite restrictions of the same certified PDE path."""
    common = {
        "provenance": provenance,
        "state_tokens": state_tokens,
        "elapsed_times": elapsed_times,
        "strain_action": strain_action,
        "residual_impulse": residual_impulse,
        "hh_impulse": hh_impulse,
    }
    whole = partition_same_carrier_path(**common, checkpoint_indices=())
    segmented = partition_same_carrier_path(**common, checkpoint_indices=checkpoint_indices)
    left = same_carrier_first_exit(whole, tie_tolerance=tie_tolerance)
    right = same_carrier_first_exit(segmented, tie_tolerance=tie_tolerance)
    if left["first_elapsed"] != right["first_elapsed"]:
        raise AssertionError("observer cuts moved the same-path first-stop time")
    if tuple(left["joint_first_stops"]) != tuple(right["joint_first_stops"]):
        raise AssertionError("observer cuts changed the same-path exact first-stop set")
    return {
        "first_elapsed": left["first_elapsed"],
        "joint_first_stops": tuple(left["joint_first_stops"]),
        "first_time_residual": 0.0,
        "checkpoint_count": len(tuple(checkpoint_indices)),
        "carrier_restarts": 0,
        "monitor_resets": 0,
        "segmentation_changed_first_hit": False,
        "complex_phase_preserved": True,
        "same_pde_trajectory_certified": True,
    }


def checkpoint_continuation_policy(
    checkpoint_record: FullNaturalCheckpoint | object,
    *,
    provenance: SameCarrierProvenance | None = None,
    request_carrier_replacement: bool = False,
    request_terminal_amplitude_reset: bool = False,
    request_monitor_reset: bool = False,
) -> dict[str, object]:
    """Authorize continuation only from a typed checkpoint bound to this carrier."""
    if not isinstance(checkpoint_record, FullNaturalCheckpoint):
        raise TypeError("typed FullNaturalCheckpoint record required for same-carrier continuation")
    if not isinstance(provenance, SameCarrierProvenance):
        raise TypeError("typed same-carrier provenance required at the checkpoint")
    if checkpoint_record.terminal_time != provenance.terminal_time:
        raise TypeError("checkpoint terminal time is not the event terminal time token")
    if checkpoint_record.corridor_frequency != provenance.carrier_frequency:
        raise TypeError("checkpoint corridor scale is not the fixed carrier frequency token")
    if checkpoint_record.scaled_lifetime != provenance.scaled_lifetime:
        raise TypeError("checkpoint rebound the fixed carrier lifetime token")
    if request_carrier_replacement:
        raise TypeError("a no-event checkpoint cannot replace the event-anchored smooth carrier")
    if request_terminal_amplitude_reset:
        raise TypeError("a no-event checkpoint cannot reset the terminal coefficient baseline")
    if request_monitor_reset:
        raise TypeError("a no-event checkpoint cannot reset cumulative native first-hit monitors")
    return {
        "canonical_continuation": SAME_CARRIER_CONTINUATION,
        "event_id": provenance.event_id,
        "carrier_id": provenance.carrier_id,
        "terminal_dual_id": provenance.terminal_dual_id,
        "trajectory_id": provenance.trajectory_id,
        "hard_shell_checkpoint_witnesses": "state_sidecars_only",
        "carrier_replacement_authorized": False,
        "terminal_amplitude_reset_authorized": False,
        "monitor_reset_authorized": False,
        "checkpoint_scale_path_is_physical_lineage": False,
        "hardening_requires_new_physical_event": True,
        "typed_checkpoint_and_carrier_provenance_verified": True,
    }


def fixed_carrier_natural_window_capacity(
    *,
    event_time: float,
    carrier_frequency: float,
    scaled_lifetime: float,
) -> dict[str, object]:
    """Finite capacity of genuine c A^-2 windows for one fixed carrier.

    Arbitrary observer cuts may accumulate, but they are not full-natural service
    windows.  Genuine fixed-A, fixed-c windows have one positive native duration and
    therefore cannot form an interior Zeno sequence before t=0.
    """
    t = float(event_time)
    if t <= 0 or not math.isfinite(t):
        raise ValueError("positive finite event time required")
    duration = _native_duration(carrier_frequency, scaled_lifetime)
    time_fraction = Fraction.from_float(t)
    duration_fraction = Fraction.from_float(duration)
    maximum = int(time_fraction // duration_fraction)
    remainder = float(time_fraction - maximum * duration_fraction)
    return {
        "native_window_duration": duration,
        "maximum_complete_windows_before_t0": maximum,
        "remaining_native_time_after_complete_windows": remainder,
        "interior_zeno_possible": False,
        "arbitrary_observer_cuts_are_natural_windows": False,
        "carrier_frequency_fixed": True,
        "scaled_lifetime_fixed": True,
    }


def maximal_same_carrier_outcome(
    segments: Sequence[SameCarrierMonitorSegment],
    *,
    event_time: float,
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Continue the certified path to its true first stop, t=0, or a later sample."""
    path = join_same_carrier_segments(segments)
    provenance = path["provenance"]
    if not isinstance(provenance, SameCarrierProvenance):
        raise AssertionError("joined path lost same-carrier provenance")
    t = float(event_time)
    if t != provenance.terminal_time:
        raise TypeError("event time was rebound outside same-carrier provenance")
    out = same_carrier_first_exit(segments, tie_tolerance=tie_tolerance)
    observed = float(out["observed_elapsed_end"])
    if observed > t:
        raise ValueError("same-carrier observation extends before the initial boundary")
    if out["first_elapsed"] is not None:
        first = float(out["first_elapsed"])
        return {
            **out,
            "endpoint_time": t - first,
            "hits_initial_boundary": False,
            "physical_event_created_by_checkpoint": False,
        }
    if observed == t:
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
    remaining = t - observed
    if remaining <= 0:
        raise AssertionError("positive native remaining time was lost")
    return {
        **out,
        "classification": "same_carrier_event_free_continuation",
        "endpoint_time": remaining,
        "hits_initial_boundary": False,
        "remaining_backward_time": remaining,
        "physical_event_created_by_checkpoint": False,
        "next_action": "continue_the_same_event_anchored_carrier_and_cumulative_monitors",
    }


@dataclass(frozen=True)
class SameCarrierPrelimitCertificate:
    """Actual cumulative path restrictions approaching one proposed cut limit."""

    segments: tuple[SameCarrierMonitorSegment, ...]

    def __post_init__(self) -> None:
        segments = tuple(self.segments)
        join_same_carrier_segments(segments)
        object.__setattr__(self, "segments", segments)


@dataclass(frozen=True)
class SmoothPDEExtensionToken:
    """Typed statement that the same trajectory is smooth on an open interval."""

    trajectory_id: str
    state_token: str
    physical_time: float
    open_interval: tuple[float, float]

    def __post_init__(self) -> None:
        if not self.trajectory_id or not self.state_token:
            raise ValueError("nonempty trajectory and state tokens required")
        s = float(self.physical_time)
        interval = tuple(float(x) for x in self.open_interval)
        if len(interval) != 2 or not all(math.isfinite(x) for x in (s, *interval)):
            raise ValueError("finite physical time and open smooth interval required")
        if not interval[0] < s < interval[1]:
            raise ValueError("smooth extension interval must be open around the accumulation time")
        object.__setattr__(self, "physical_time", s)
        object.__setattr__(self, "open_interval", (interval[0], interval[1]))


def interior_checkpoint_accumulation_outcome(
    *,
    event_time: float,
    accumulation_time: float,
    prelimit_certificate: SameCarrierPrelimitCertificate,
    smooth_extension_token: SmoothPDEExtensionToken | None,
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Classify a cut accumulation from the actual no-earlier-hit PDE path.

    Endpoint scalars alone cannot certify a first hit.  The cumulative complex path
    must show no earlier debut, and strict-margin continuation at an interior time
    additionally requires an open smooth-extension token for the same trajectory.
    """
    if not isinstance(prelimit_certificate, SameCarrierPrelimitCertificate):
        raise TypeError("typed actual prelimit path certificate required")
    path = join_same_carrier_segments(prelimit_certificate.segments)
    provenance = path["provenance"]
    if not isinstance(provenance, SameCarrierProvenance):
        raise AssertionError("prelimit path lost same-carrier provenance")
    t = float(event_time)
    s = float(accumulation_time)
    if t != provenance.terminal_time:
        raise TypeError("accumulation event time was rebound outside carrier provenance")
    if s < 0 or s >= t or not math.isfinite(s):
        raise ValueError("interior-or-zero finite accumulation time required")
    elapsed_limit = t - s
    observed = float(tuple(path["elapsed_times"])[-1])
    if observed != elapsed_limit:
        raise TypeError("prelimit path endpoint is not the native accumulation-time token")
    out = same_carrier_first_exit(prelimit_certificate.segments, tie_tolerance=tie_tolerance)
    first = out["first_elapsed"]
    if first is not None and float(first) < elapsed_limit:
        raise ValueError("an earlier physical first hit contradicts a no-hit checkpoint accumulation")

    if s == 0.0:
        return {
            "classification": "absorbing_initial_boundary",
            "joint_first_stops": ("t=0",),
            "requires_physical_energy_reentry": False,
            "checkpoint_accumulation_is_obstruction": False,
            "same_carrier_extends_past_accumulation": False,
            "prelimit_path_certified": True,
        }
    if first is not None:
        if float(first) != elapsed_limit:
            raise ValueError("closed face did not occur at the accumulation endpoint")
        stops = tuple(out["joint_first_stops"])
        return {
            "classification": "first_stop_at_interior_checkpoint_accumulation",
            "joint_first_stops": stops,
            "requires_physical_energy_reentry": any(
                label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
                for label in stops
            ),
            "coefficient_impulses_used_as_work": False,
            "checkpoint_accumulation_is_obstruction": False,
            "same_carrier_extends_past_accumulation": False,
            "prelimit_path_certified": True,
        }

    if not isinstance(smooth_extension_token, SmoothPDEExtensionToken):
        raise TypeError("strict-margin continuation requires a typed smooth PDE extension token")
    if smooth_extension_token.trajectory_id != provenance.trajectory_id:
        raise TypeError("smooth extension belongs to a different PDE trajectory")
    if smooth_extension_token.state_token != tuple(path["state_tokens"])[-1]:
        raise TypeError("smooth extension is not attached to the accumulation endpoint state")
    if smooth_extension_token.physical_time != s:
        raise TypeError("smooth extension is attached to a different physical-time token")
    return {
        "classification": "same_carrier_extends_across_interior_checkpoint_accumulation",
        "joint_first_stops": (),
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
        "checkpoint_accumulation_is_obstruction": False,
        "same_carrier_extends_past_accumulation": True,
        "prelimit_path_certified": True,
        "smooth_pde_extension_certified": True,
        "continuation_reason": "the actual cumulative path has strict native margins and the same PDE trajectory is smooth on an open interval around the cut accumulation",
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "fixed_carrier": "between genuine physical events the canonical object carries one exact event id, smooth carrier, terminal dual, terminal coefficient, scale/lifetime, terminal state token and actual PDE trajectory id",
        "cumulative_monitors": "K_A[s,t] is monotone; the actual complex I_role-interface[s,t] and I_HH[s,t] paths are retained, their magnitudes are derived, and segment magnitudes are never added or used as work",
        "segmentation": "finite checkpoint insertion is only restriction and exact gluing of the same cumulative PDE path; exact shared state/time/complex-boundary tokens leave the first stop unchanged",
        "checkpoint_policy": "only a typed FullNaturalCheckpoint bound to the fixed carrier scale, lifetime and event time can authorize sidecar rereading; a dictionary or close floating data cannot reset or replace the carrier",
        "natural_windows": "arbitrary observer cuts are not full-natural service windows; for fixed A and c every genuine window has one fixed positive native duration c A^-2, so fixed-carrier natural windows cannot have an interior Zeno accumulation",
        "interior_accumulation": "an arbitrary-cut accumulation is classified only from the actual cumulative prelimit path: an exact endpoint face is the first stop, an earlier crossing invalidates the no-hit premise, and strict-margin continuation requires a matching open smooth-PDE extension token",
        "boundary": "t=0 is absorbing only when the native cumulative path reaches exactly the full event time; positive remaining native time is never rounded away",
        "high_tail_separation": "actual high-tail dissipation/work remains a separate certified physical UV route; diagnostic shell readings carry neither fictitious c A_j^-2 durations nor causal scale lineage",
        "scope": "this closes observer segmentation for one certified event-anchored PDE path; it does not telescope infinitely recurring genuine physical owner events and does not prove Navier-Stokes global regularity",
    }


@dataclass(frozen=True)
class SameCarrierStress:
    samples: int
    worst_segmentation_first_time_residual: float
    segmentation_failures: int
    reset_barrier_failures: int
    nonmonotone_impulse_magnitudes: int
    accumulation_stop_cases: int
    accumulation_continue_cases: int
    maximum_checkpoint_count: int
    fixed_window_zeno_failures: int


def _engineered_path(
    rng: random.Random,
) -> tuple[float, tuple[float, ...], tuple[float, ...], tuple[complex, ...], tuple[complex, ...]]:
    n = rng.randint(8, 30)
    times = tuple(j / (n - 1) for j in range(n))
    amplitude = math.exp(rng.uniform(-12.0, 12.0))
    mode = rng.randrange(5)
    if mode == 1:
        tau = times[rng.randint(2, n - 3)]
        strain = tuple((LOW_STRAIN_ACTION / tau) * x for x in times)
    else:
        end = rng.uniform(0.05, 0.85) * LOW_STRAIN_ACTION
        exponent = rng.uniform(0.7, 1.8)
        strain = tuple(end * x**exponent for x in times)

    def wave(threshold: float, phase: float) -> tuple[complex, ...]:
        values: list[complex] = []
        for x in times:
            radius = max(0.0, min(0.82, 0.42 + 0.25 * math.sin(7.0 * x + phase))) * threshold
            values.append(radius * cmath.exp(1j * (phase + 9.0 * x)))
        values[0] = 0.0j
        return tuple(values)

    residual = wave(RESIDUAL_FRACTION * amplitude, rng.uniform(-math.pi, math.pi))
    hh = wave(GENERATED_FRACTION * amplitude, rng.uniform(-math.pi, math.pi))
    if mode in (2, 4):
        tau = times[rng.randint(2, n - 3)]
        residual = tuple(complex((RESIDUAL_FRACTION * amplitude / tau) * x) for x in times)
    if mode in (3, 4):
        tau = times[rng.randint(2, n - 3)]
        hh = tuple(complex((GENERATED_FRACTION * amplitude / tau) * x) for x in times)
    return amplitude, times, strain, residual, hh


def _provenance(index: int, amplitude: float, terminal_time: float = 2.0) -> SameCarrierProvenance:
    return SameCarrierProvenance(
        event_id=f"event-{index}",
        carrier_id=f"carrier-{index}",
        terminal_dual_id=f"dual-{index}",
        trajectory_id=f"trajectory-{index}",
        terminal_state_token=f"state-{index}-0",
        terminal_time=terminal_time,
        carrier_frequency=2.0,
        scaled_lifetime=1.0,
        terminal_coefficient=complex(amplitude),
    )


def stress(samples: int = 50_000, seed: int = 20260811) -> SameCarrierStress:
    rng = random.Random(seed)
    worst = 0.0
    segmentation_failures = 0
    reset_failures = 0
    nonmonotone = 0
    stop_cases = 0
    continue_cases = 0
    maximum_cuts = 0
    zeno_failures = 0
    capacity = fixed_carrier_natural_window_capacity(
        event_time=1.0,
        carrier_frequency=2.0,
        scaled_lifetime=1.0,
    )
    if capacity["interior_zeno_possible"] or capacity["maximum_complete_windows_before_t0"] != 4:
        zeno_failures += 1
        raise AssertionError("positive fixed-carrier natural duration admitted interior Zeno")

    for index in range(samples):
        amplitude, times, strain, residual, hh = _engineered_path(rng)
        provenance = _provenance(index, amplitude)
        tokens = tuple(f"state-{index}-{j}" for j in range(len(times)))
        possible = list(range(1, len(times) - 1))
        rng.shuffle(possible)
        cuts = tuple(sorted(possible[: rng.randint(0, min(len(possible), 12))]))
        maximum_cuts = max(maximum_cuts, len(cuts))
        invariant = segmentation_invariance(
            provenance=provenance,
            state_tokens=tokens,
            elapsed_times=times,
            strain_action=strain,
            residual_impulse=residual,
            hh_impulse=hh,
            checkpoint_indices=cuts,
            tie_tolerance=0.0,
        )
        current = float(invariant["first_time_residual"])
        worst = max(worst, current)
        if current != 0.0 or invariant["segmentation_changed_first_hit"]:
            segmentation_failures += 1
            raise AssertionError("observer segmentation changed a same-PDE-path first hit")
        residual_abs = tuple(abs(x) for x in residual)
        hh_abs = tuple(abs(x) for x in hh)
        if any(residual_abs[j + 1] < residual_abs[j] for j in range(len(times) - 1)) or any(
            hh_abs[j + 1] < hh_abs[j] for j in range(len(times) - 1)
        ):
            nonmonotone += 1

        if cuts:
            segments = list(
                partition_same_carrier_path(
                    provenance=provenance,
                    state_tokens=tokens,
                    elapsed_times=times,
                    strain_action=strain,
                    residual_impulse=residual,
                    hh_impulse=hh,
                    checkpoint_indices=cuts,
                )
            )
            bad = segments[1]
            segments[1] = SameCarrierMonitorSegment(
                provenance=bad.provenance,
                state_tokens=bad.state_tokens,
                elapsed_times=bad.elapsed_times,
                strain_action=(0.0,) + bad.strain_action[1:],
                residual_impulse=(0.0j,) + bad.residual_impulse[1:],
                hh_impulse=(0.0j,) + bad.hh_impulse[1:],
            )
            try:
                join_same_carrier_segments(segments)
            except TypeError:
                pass
            else:
                reset_failures += 1
                raise AssertionError("checkpoint reset crossed exact same-carrier boundary tokens")

        elapsed_limit = 1.3
        acc_provenance = _provenance(index + samples, amplitude, terminal_time=2.0)
        acc_tokens = (acc_provenance.terminal_state_token, f"acc-state-{index}")
        if index % 2 == 0:
            acc_segment = SameCarrierMonitorSegment(
                provenance=acc_provenance,
                state_tokens=acc_tokens,
                elapsed_times=(0.0, elapsed_limit),
                strain_action=(0.0, 0.7 * LOW_STRAIN_ACTION),
                residual_impulse=(0.0j, complex(0.7 * RESIDUAL_FRACTION * amplitude)),
                hh_impulse=(0.0j, complex(0.7 * GENERATED_FRACTION * amplitude)),
            )
            extension = SmoothPDEExtensionToken(
                trajectory_id=acc_provenance.trajectory_id,
                state_token=acc_tokens[-1],
                physical_time=0.7,
                open_interval=(0.6, 0.8),
            )
            outcome = interior_checkpoint_accumulation_outcome(
                event_time=2.0,
                accumulation_time=0.7,
                prelimit_certificate=SameCarrierPrelimitCertificate((acc_segment,)),
                smooth_extension_token=extension,
            )
            if not outcome["same_carrier_extends_past_accumulation"]:
                raise AssertionError("strict actual margins failed to cross an observer-cut accumulation")
            continue_cases += 1
        else:
            acc_segment = SameCarrierMonitorSegment(
                provenance=acc_provenance,
                state_tokens=acc_tokens,
                elapsed_times=(0.0, elapsed_limit),
                strain_action=(0.0, LOW_STRAIN_ACTION),
                residual_impulse=(0.0j, complex(0.4 * RESIDUAL_FRACTION * amplitude)),
                hh_impulse=(0.0j, complex(0.4 * GENERATED_FRACTION * amplitude)),
            )
            outcome = interior_checkpoint_accumulation_outcome(
                event_time=2.0,
                accumulation_time=0.7,
                prelimit_certificate=SameCarrierPrelimitCertificate((acc_segment,)),
                smooth_extension_token=None,
            )
            if STRAIN_STOP not in tuple(outcome["joint_first_stops"]):
                raise AssertionError("exact endpoint strain face disappeared")
            stop_cases += 1

    return SameCarrierStress(
        samples,
        worst,
        segmentation_failures,
        reset_failures,
        nonmonotone,
        stop_cases,
        continue_cases,
        maximum_cuts,
        zeno_failures,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=50_000)
    parser.add_argument("--outdir", type=Path, default=Path("results-same-carrier-checkpoint-segmentation-quotient"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    certificate = theorem_certificate()
    result = stress(args.samples)
    payload = {"certificate": certificate, "stress": asdict(result)}
    (args.outdir / "same_carrier_checkpoint_segmentation_quotient.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    summary = f"""# Same-carrier checkpoint segmentation quotient

Status: **{certificate['status']}**.

Finite observer cuts are quotiented only when exact event/carrier/terminal-dual,
terminal-coefficient, scale/lifetime, PDE-trajectory, shared-state, native-time and
cumulative-complex-impulse provenance reconstructs one physical path.  Complex
phase is retained when locating coefficient faces; magnitudes are never summed as
segment work.

Arbitrary observer cuts are not full-natural service windows.  For one fixed
carrier, every genuine natural window has the same positive native duration
`c A^-2`, so such windows cannot form an interior Zeno sequence.  An accumulation
of arbitrary cuts is classified only from the actual prelimit path and a matching
smooth-PDE extension token.

Stress: `{result.samples}` same-PDE-path/segmentation states
- worst segmentation first-time residual: `{result.worst_segmentation_first_time_residual:.3e}`
- segmentation failures: `{result.segmentation_failures}`
- checkpoint reset-barrier failures: `{result.reset_barrier_failures}`
- sampled nonmonotone impulse-magnitude paths: `{result.nonmonotone_impulse_magnitudes}`
- interior accumulation stop cases: `{result.accumulation_stop_cases}`
- interior accumulation continuation cases: `{result.accumulation_continue_cases}`
- maximum inserted observer-cut count: `{result.maximum_checkpoint_count}`
- fixed-window Zeno failures: `{result.fixed_window_zeno_failures}`

This closes segmentation only for one certified event-anchored PDE path.  It does
not telescope recurring genuine physical owners and does not prove Navier--Stokes
global regularity.
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
