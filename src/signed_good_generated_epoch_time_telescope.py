from __future__ import annotations

import argparse
import json
import math
import random
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from src.asynchronous_duhamel_sync import (
    BACKWARD_FRACTION,
    INITIAL_HALF_SPAN,
    LIFETIME_GROWTH_MIN,
    MIN_REFERENCE_BACKSTEP,
    PARENT_CHILD_HI,
    PARENT_CHILD_LO,
    SHARP_SYNC_BOUND,
    interior_depth_upper,
    minimum_backward_displacement,
)
from src.common_slice_coefficient_registration import HH_COEFFICIENT_OBSTRUCTION

STATUS = (
    "EXACT_SIGNED_GOOD_GENERATED_HH_PARABOLIC_EPOCH_TELESCOPE__"
    "ACTUAL_PHYSICAL_HH_WORK_AFTER_ENERGY_REENTRY_ONLY__"
    "ASYNCHRONOUS_COMMON_SLICE_BACKSHIFT_TO_T0__"
    "NO_DUHAMEL_WEIGHT_OR_EVENT_COUNT_BUDGET"
)

ACTUAL_HH_GENERATION_BRANCH = "physical_high_high_transfer_generation"
ACTUAL_HH_CAUSAL_MEASURE = "actual_positive_hh_child_energy_work"
NATIVE_RELATIVE_TOLERANCE = 8.0e-12


def _finite_positive_ratio(numerator: float, denominator: float) -> float:
    x = float(numerator)
    y = float(denominator)
    if not (math.isfinite(x) and x > 0.0 and math.isfinite(y) and y > 0.0):
        raise ValueError("positive finite native quantities required")
    try:
        ratio = math.exp(math.log(x) - math.log(y))
    except OverflowError as exc:
        raise ValueError("native ratio left the positive finite floating range") from exc
    if not math.isfinite(ratio) or ratio <= 0.0:
        raise ValueError("native ratio left the positive finite floating range")
    return ratio


def _native_lower_holds(value: float, lower: float) -> bool:
    """Compare positive same-unit data without an absolute observer unit."""
    log_ratio = math.log(float(value)) - math.log(float(lower))
    return log_ratio >= math.log1p(-NATIVE_RELATIVE_TOLERANCE)


def _native_upper_holds(value: float, upper: float) -> bool:
    log_ratio = math.log(float(value)) - math.log(float(upper))
    return log_ratio <= math.log1p(NATIVE_RELATIVE_TOLERANCE)


def _parabolic_lifetime(scaled_lifetime: float, frequency: float) -> float:
    log_lifetime = math.log(float(scaled_lifetime)) - 2.0 * math.log(float(frequency))
    try:
        lifetime = math.exp(log_lifetime)
    except OverflowError as exc:
        raise ValueError("parabolic lifetime left the finite native range") from exc
    if not math.isfinite(lifetime) or lifetime <= 0.0:
        raise ValueError("parabolic lifetime left the positive finite native range")
    return lifetime


def _require_native_time_resolution(time: float, lifetime: float) -> None:
    """Fail closed when a float clock cannot resolve one native lifetime."""
    t = float(time)
    T = float(lifetime)
    if t > 0.0 and math.ulp(t) > NATIVE_RELATIVE_TOLERANCE * T:
        raise ValueError("native backward-time resolution is too coarse for this lifetime")


def _native_common_reference_slice(
    interval_start: float,
    interval_end: float,
    parent_lifetime: float,
) -> float:
    """Dimensionless form of the common-slice law at the PDE's native clock."""
    a = float(interval_start)
    b = float(interval_end)
    T = float(parent_lifetime)
    _require_native_time_resolution(a, T)
    _require_native_time_resolution(b, T)
    width_coordinate = (b - a) / T
    if width_coordinate < -NATIVE_RELATIVE_TOLERANCE:
        raise ValueError("invalid physical support interval")
    if width_coordinate > float(SHARP_SYNC_BOUND) * (1.0 + NATIVE_RELATIVE_TOLERANCE):
        raise ValueError("event slab lies outside the sharp synchronization cone")
    if width_coordinate + float(BACKWARD_FRACTION) >= 1.0 - NATIVE_RELATIVE_TOLERANCE:
        raise AssertionError("common reference slice left a natural backward window")
    return math.fsum((a, -float(BACKWARD_FRACTION) * T))


@dataclass(frozen=True)
class SignedGoodGeneratedWorkProvenance:
    """Identity of one actual positive HH-work law on one PDE history."""

    event_id: str
    trajectory_id: str
    child_carrier_id: str
    generated_parent_carrier_id: str
    work_law_id: str
    child_frequency: float
    parent_frequency: float
    scaled_lifetime: float
    slab_start: float
    slab_end: float

    def __post_init__(self) -> None:
        identities = (
            self.event_id,
            self.trajectory_id,
            self.child_carrier_id,
            self.generated_parent_carrier_id,
            self.work_law_id,
        )
        if not all(isinstance(x, str) and x for x in identities):
            raise ValueError("nonempty generated-work provenance identities required")
        positive = (
            self.child_frequency,
            self.parent_frequency,
            self.scaled_lifetime,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in positive):
            raise ValueError("positive finite generated-work provenance scales required")
        times = (self.slab_start, self.slab_end)
        if not all(math.isfinite(x) and x >= 0.0 for x in times):
            raise ValueError("finite nonnegative physical provenance times required")
        child_lifetime = _parabolic_lifetime(self.scaled_lifetime, self.child_frequency)
        _require_native_time_resolution(self.slab_start, child_lifetime)
        _require_native_time_resolution(self.slab_end, child_lifetime)
        if self.slab_end <= self.slab_start:
            raise ValueError("positive ordered physical provenance slab required")
        if not _native_upper_holds(self.slab_end - self.slab_start, child_lifetime):
            raise ValueError("generated-work provenance exceeds one child natural slab")


@dataclass(frozen=True)
class SignedGoodGeneratedHHStep:
    """One actual signed-good HH generation and its physical parent-work support.

    The step is created canonically by :func:`signed_good_step_from_energy_reentry`.
    ``work_support_start/end`` are the times of the selected heavy half of the
    **actual positive child-energy HH work law**, viewed simultaneously as the
    support of the hard parent-pair events.  The parent role has its own natural
    lifetime ``c N_parent^-2`` and its common registration surface is

        s = work_support_start - (2/5) T_parent.

    This surface is not a recursive event.  It is the physical backward surface
    on which the same generated parent carrier is registered.
    """

    child_frequency: float
    parent_frequency: float
    scaled_lifetime: float
    work_support_start: float
    work_support_end: float
    physical_hh_work_mass: float
    physical_hh_work_total: float
    physical_hh_work_lower: float
    provenance: SignedGoodGeneratedWorkProvenance
    energy_reentry_branch: str = ACTUAL_HH_GENERATION_BRANCH
    causal_measure: str = ACTUAL_HH_CAUSAL_MEASURE
    coefficient_impulse_used_as_work: bool = False
    observer_partition_motion_charged_as_physics: bool = False

    def __post_init__(self) -> None:
        positive = (
            self.child_frequency,
            self.parent_frequency,
            self.scaled_lifetime,
            self.physical_hh_work_mass,
            self.physical_hh_work_total,
            self.physical_hh_work_lower,
        )
        if not all(math.isfinite(x) and x > 0 for x in positive):
            raise ValueError("positive finite signed-good generated-HH data required")
        times = (self.work_support_start, self.work_support_end)
        if not all(math.isfinite(x) and x >= 0 for x in times) or self.work_support_end < self.work_support_start:
            raise ValueError("finite nonnegative ordered physical work support required")
        if self.energy_reentry_branch != ACTUAL_HH_GENERATION_BRANCH:
            raise TypeError("only actual physical HH generation after energy reentry may enter the signed-good epoch")
        if self.causal_measure != ACTUAL_HH_CAUSAL_MEASURE:
            raise TypeError("generated epoch must use actual positive HH child-energy work as its causal law")
        if self.coefficient_impulse_used_as_work:
            raise TypeError("raw HH coefficient impulse cannot be used as physical work in a generated epoch")
        if self.observer_partition_motion_charged_as_physics:
            raise TypeError("observer partition motion cannot be used as physical generated work")

        if not isinstance(self.provenance, SignedGoodGeneratedWorkProvenance):
            raise TypeError("typed actual-PDE generated-work provenance required")
        if not math.isclose(
            _finite_positive_ratio(self.child_frequency, self.provenance.child_frequency),
            1.0,
            rel_tol=NATIVE_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        ):
            raise ValueError("step child frequency disagrees with its physical work provenance")
        if not math.isclose(
            _finite_positive_ratio(self.parent_frequency, self.provenance.parent_frequency),
            1.0,
            rel_tol=NATIVE_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        ):
            raise ValueError("step parent frequency disagrees with its physical work provenance")
        if not math.isclose(
            _finite_positive_ratio(self.scaled_lifetime, self.provenance.scaled_lifetime),
            1.0,
            rel_tol=NATIVE_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        ):
            raise ValueError("step lifetime constant disagrees with its physical work provenance")
        slab_duration = self.provenance.slab_end - self.provenance.slab_start
        slab_tol = NATIVE_RELATIVE_TOLERANCE * slab_duration
        if (
            self.work_support_start < self.provenance.slab_start - slab_tol
            or self.work_support_end > self.provenance.slab_end + slab_tol
        ):
            raise ValueError("selected HH-work support left its provenance slab")

        ratio = _finite_positive_ratio(self.parent_frequency, self.child_frequency)
        ratio_tol = 8.0e-13
        if not ratio > float(PARENT_CHILD_LO) * (1.0 + ratio_tol):
            raise ValueError("signed-good generated parent violates the strict 3/5 lower frequency ratio")
        if not ratio < float(PARENT_CHILD_HI) * (1.0 - ratio_tol):
            raise ValueError("signed-good generated parent violates the strict 5/8 upper frequency ratio")

        if not _native_upper_holds(self.physical_hh_work_mass, self.physical_hh_work_total):
            raise ValueError("selected positive HH-work mass cannot exceed its total")
        if not _native_lower_holds(self.physical_hh_work_mass, 0.5 * self.physical_hh_work_total):
            raise ValueError("selected generated work support does not carry a physical heavy half")
        if not _native_lower_holds(self.physical_hh_work_total, self.physical_hh_work_lower):
            raise ValueError("actual positive HH work law does not realize the energy-gate lower bound")

        width = self.work_support_end - self.work_support_start
        child_half = 0.5 * self.child_natural_lifetime
        if width > 0.0 and not _native_upper_holds(width, child_half):
            raise ValueError("selected generated support is wider than one half child natural slab")
        parent_span = float(INITIAL_HALF_SPAN) * self.parent_natural_lifetime
        if width > 0.0 and not _native_upper_holds(width, parent_span):
            raise AssertionError("signed-good parent support exceeded the certified 25/128 natural-lifetime span")

    @property
    def child_natural_lifetime(self) -> float:
        return _parabolic_lifetime(self.scaled_lifetime, self.child_frequency)

    @property
    def parent_natural_lifetime(self) -> float:
        return _parabolic_lifetime(self.scaled_lifetime, self.parent_frequency)

    @property
    def parent_child_ratio(self) -> float:
        return _finite_positive_ratio(self.parent_frequency, self.child_frequency)

    @property
    def normalized_parent_span(self) -> float:
        return (self.work_support_end - self.work_support_start) / self.parent_natural_lifetime

    @property
    def common_reference_time(self) -> float:
        return math.fsum(
            (
                self.work_support_start,
                -float(BACKWARD_FRACTION) * self.parent_natural_lifetime,
            )
        )


def _physical_energy_gate(reentry: Mapping[str, object]) -> Mapping[str, object]:
    nested = reentry.get("energy_gate")
    if isinstance(nested, Mapping):
        return nested
    return reentry


def signed_good_step_from_energy_reentry(
    *,
    reentry: Mapping[str, object],
    selected_physical_half_slab: Mapping[str, object],
    child_frequency: float,
    parent_frequency: float,
    scaled_lifetime: float,
) -> SignedGoodGeneratedHHStep:
    """Canonical type boundary: coefficient locator -> Q^2 energy -> actual HH work.

    A raw ``HH_COEFFICIENT_OBSTRUCTION`` is deliberately insufficient.  The
    reentry must contain an actual physical-energy gate whose branch is
    ``physical_high_high_transfer_generation`` and whose positive HH-work lower is
    realized by the physical work law.  The selected half-slab must be a positive
    sublaw carrying at least half of that same work.
    """
    top_branch = str(reentry.get("branch", ""))
    if top_branch == HH_COEFFICIENT_OBSTRUCTION:
        raise TypeError("HH coefficient obstruction is only an interval locator; actual Q^2 energy reentry is required")
    if bool(reentry.get("requires_physical_energy_reentry", False)) and "energy_gate" not in reentry:
        raise TypeError("coefficient obstruction has not yet passed through actual physical-energy reentry")
    if bool(reentry.get("coefficient_impulse_used_as_physical_work", False)):
        raise TypeError("coefficient impulse magnitude cannot be promoted to physical HH work")
    if bool(reentry.get("observer_partition_motion_charged_as_physics", False)):
        raise TypeError("observer partition motion cannot be promoted to physical HH work")

    gate = _physical_energy_gate(reentry)
    if bool(gate.get("coefficient_impulse_used_as_physical_work", False)):
        raise TypeError("coefficient impulse magnitude cannot be promoted to physical HH work")
    if bool(gate.get("observer_partition_motion_charged_as_physics", False)):
        raise TypeError("observer partition motion cannot be promoted to physical HH work")
    branch = str(gate.get("branch", top_branch))
    if branch != ACTUAL_HH_GENERATION_BRANCH:
        raise TypeError("energy reentry did not select actual physical high-high generation")
    try:
        work_lower = float(gate["physical_hh_work_lower"])
    except (KeyError, TypeError, ValueError) as exc:
        raise TypeError("generated epoch requires the physical HH-work lower from the energy gate") from exc
    if not math.isfinite(work_lower) or work_lower <= 0:
        raise TypeError("energy gate must certify positive actual HH work")

    reentry_provenance = reentry.get("provenance")
    gate_provenance = gate.get("provenance")
    half_provenance = selected_physical_half_slab.get("provenance")
    if not isinstance(reentry_provenance, SignedGoodGeneratedWorkProvenance):
        raise TypeError("energy reentry requires typed actual-PDE work provenance")
    if not isinstance(gate_provenance, SignedGoodGeneratedWorkProvenance):
        raise TypeError("physical energy gate requires typed actual-PDE work provenance")
    if not isinstance(half_provenance, SignedGoodGeneratedWorkProvenance):
        raise TypeError("selected physical half-slab requires typed actual-PDE work provenance")
    if not (reentry_provenance == gate_provenance == half_provenance):
        raise TypeError("energy gate and half-slab provenance identify foreign physical work laws")
    if bool(selected_physical_half_slab.get("coefficient_impulse_used_as_physical_work", False)):
        raise TypeError("coefficient impulse magnitude cannot be promoted to physical HH work")
    if bool(selected_physical_half_slab.get("observer_partition_motion_charged_as_physics", False)):
        raise TypeError("observer partition motion cannot be promoted to physical HH work")

    try:
        start = float(selected_physical_half_slab["start"])
        end = float(selected_physical_half_slab["end"])
        mass = float(selected_physical_half_slab["mass"])
        total = float(selected_physical_half_slab["total"])
        span_bound = float(selected_physical_half_slab["normalized_parent_span_upper"])
    except (KeyError, TypeError, ValueError) as exc:
        raise TypeError("selected physical HH half-slab certificate is incomplete") from exc
    if not math.isfinite(span_bound) or span_bound < 0:
        raise TypeError("selected HH work support has an invalid normalized parent-span certificate")
    if span_bound > float(INITIAL_HALF_SPAN) * (1.0 + 8.0e-13):
        raise TypeError("selected HH work support lacks the certified signed-good parent-span bound")

    return SignedGoodGeneratedHHStep(
        child_frequency=float(child_frequency),
        parent_frequency=float(parent_frequency),
        scaled_lifetime=float(scaled_lifetime),
        work_support_start=start,
        work_support_end=end,
        physical_hh_work_mass=mass,
        physical_hh_work_total=total,
        physical_hh_work_lower=work_lower,
        provenance=reentry_provenance,
        energy_reentry_branch=branch,
        causal_measure=ACTUAL_HH_CAUSAL_MEASURE,
        coefficient_impulse_used_as_work=False,
        observer_partition_motion_charged_as_physics=False,
    )


@dataclass(frozen=True)
class SignedGoodGeneratedEpochCertificate:
    layer_count: int
    transition_count: int
    initial_parent_frequency: float
    final_parent_frequency: float
    initial_parent_lifetime: float
    first_common_reference_time: float
    last_common_reference_time: float
    cumulative_reference_backshift: float
    minimum_cumulative_backshift: float
    maximum_observed_parent_ratio: float
    minimum_observed_lifetime_growth: float
    total_selected_physical_hh_work: float
    interior_transition_upper: int
    total_layer_upper_before_or_at_boundary: int
    hits_initial_boundary: bool
    common_slices_are_recursive_events: bool = False
    duhamel_weights_used_as_causal_law: bool = False
    event_count_budget_used: bool = False
    generic_hh_claimed: bool = False

    def __post_init__(self) -> None:
        if self.layer_count < 1 or self.transition_count != self.layer_count - 1:
            raise ValueError("valid nonempty generated epoch layer counts required")
        positive = (
            self.initial_parent_frequency,
            self.final_parent_frequency,
            self.initial_parent_lifetime,
            self.total_selected_physical_hh_work,
        )
        if not all(math.isfinite(x) and x > 0 for x in positive):
            raise ValueError("positive finite generated epoch certificate data required")
        finite = (
            self.first_common_reference_time,
            self.last_common_reference_time,
            self.cumulative_reference_backshift,
            self.minimum_cumulative_backshift,
            self.maximum_observed_parent_ratio,
            self.minimum_observed_lifetime_growth,
        )
        if not all(math.isfinite(x) for x in finite):
            raise ValueError("finite generated epoch geometry required")
        native_time_tol = NATIVE_RELATIVE_TOLERANCE * self.initial_parent_lifetime
        if self.cumulative_reference_backshift < -native_time_tol or self.minimum_cumulative_backshift < -native_time_tol:
            raise ValueError("backward displacement cannot be negative")
        if self.interior_transition_upper < 0 or self.total_layer_upper_before_or_at_boundary < self.layer_count:
            raise ValueError("finite generated epoch depth bound must cover the observed layers")
        if self.common_slices_are_recursive_events or self.duhamel_weights_used_as_causal_law or self.event_count_budget_used:
            raise ValueError("generated epoch theorem may not manufacture event depth, Duhamel causality, or an event-count currency")
        if self.generic_hh_claimed:
            raise ValueError("the theorem is signed-good generated HH only, not a generic HH recursion theorem")


def signed_good_generated_epoch_telescope(
    steps: Sequence[SignedGoodGeneratedHHStep],
) -> SignedGoodGeneratedEpochCertificate:
    """Finite physical-time depth of one consecutive signed-good generated epoch.

    For the parent-work support ``H_j=[a_j,b_j]`` let

        T_j = c N_j^-2,   s_j = a_j - (2/5) T_j.

    Every support has normalized width at most ``25/128 < 10/39``.  Consecutive
    generated support must lie inside the previous common registration interval
    ``[s_j,b_j]``.  Signed-good scale geometry gives

        T_(j+1)/T_j > 64/25.

    The existing asynchronous algebra therefore yields

        s_j-s_(j+1) >= (1792/4875) T_j,
        s_0-s_L >= (1792/7605) T_0[(64/25)^L-1].

    These ``s_j`` are registration surfaces, not event vertices.  Once the next
    required surface reaches ``t=0`` the adjoint/physical carrier gate is truncated
    there and the initial boundary is absorbing.  Hence an interior signed-good
    generated ancestry has finite depth at finite physical time.
    """
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty signed-good generated epoch required")

    c0 = rows[0].scaled_lifetime
    for row in rows[1:]:
        if not math.isclose(
            _finite_positive_ratio(row.scaled_lifetime, c0),
            1.0,
            rel_tol=NATIVE_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        ):
            raise ValueError("one generated epoch must use one supplied scaled natural-lifetime constant")

    common = tuple(row.common_reference_time for row in rows)
    lifetimes = tuple(row.parent_natural_lifetime for row in rows)
    max_ratio = max(row.parent_child_ratio for row in rows)
    min_growth = math.inf if len(rows) == 1 else min(lifetimes[j + 1] / lifetimes[j] for j in range(len(rows) - 1))
    if len(rows) == 1:
        min_growth = float(LIFETIME_GROWTH_MIN)

    for j, row in enumerate(rows):
        _require_native_time_resolution(row.work_support_start, row.parent_natural_lifetime)
        _require_native_time_resolution(row.work_support_end, row.parent_natural_lifetime)
        if row.normalized_parent_span > float(SHARP_SYNC_BOUND) + 8e-13:
            raise ValueError("generated support left the sharp asynchronous synchronization cone")
        # Recompute through the existing common-slice helper as a regression against
        # silently changing the registration geometry.
        s = _native_common_reference_slice(
            row.work_support_start,
            row.work_support_end,
            row.parent_natural_lifetime,
        )
        if abs(s - common[j]) > NATIVE_RELATIVE_TOLERANCE * row.parent_natural_lifetime:
            raise AssertionError("generated common reference surface disagrees with the certified asynchronous theorem")

        if common[j] <= 0.0 and j + 1 < len(rows):
            raise ValueError("generated ancestry continued after its required registration surface had reached t=0")

        if j + 1 == len(rows):
            continue
        nxt = rows[j + 1]
        if row.provenance.trajectory_id != nxt.provenance.trajectory_id:
            raise ValueError("generated epoch spliced foreign PDE trajectories")
        if not math.isclose(
            _finite_positive_ratio(nxt.child_frequency, row.parent_frequency),
            1.0,
            rel_tol=NATIVE_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        ):
            raise ValueError("generated epoch must continue through the actual signed-good parent carrier scale")
        if row.provenance.generated_parent_carrier_id != nxt.provenance.child_carrier_id:
            raise ValueError("generated epoch spliced a foreign physical parent carrier")

        _require_native_time_resolution(nxt.work_support_start, row.parent_natural_lifetime)
        _require_native_time_resolution(nxt.work_support_end, row.parent_natural_lifetime)
        lower_coordinate = (nxt.work_support_start - common[j]) / row.parent_natural_lifetime
        upper_coordinate = (nxt.work_support_end - row.work_support_end) / row.parent_natural_lifetime
        if lower_coordinate < -NATIVE_RELATIVE_TOLERANCE or upper_coordinate > NATIVE_RELATIVE_TOLERANCE:
            raise ValueError("next generated physical-work support is not contained in the previous common registration interval")

        growth = _finite_positive_ratio(nxt.parent_natural_lifetime, row.parent_natural_lifetime)
        if growth < float(LIFETIME_GROWTH_MIN) * (1.0 - NATIVE_RELATIVE_TOLERANCE):
            raise AssertionError("signed-good generated parent natural lifetime failed the 64/25 growth law")
        backshift = common[j] - common[j + 1]
        required = float(MIN_REFERENCE_BACKSTEP) * row.parent_natural_lifetime
        if backshift / row.parent_natural_lifetime < float(MIN_REFERENCE_BACKSTEP) * (
            1.0 - NATIVE_RELATIVE_TOLERANCE
        ):
            raise AssertionError("asynchronous generated common surfaces did not move backward by the certified physical amount")

    transitions = len(rows) - 1
    cumulative = common[0] - common[-1]
    minimum = minimum_backward_displacement(lifetimes[0], transitions)
    if transitions:
        if cumulative <= 0.0 or minimum <= 0.0 or not _native_lower_holds(cumulative, minimum):
            raise AssertionError("generated epoch lost the geometric physical-time backshift telescope")

    first_s = common[0]
    if first_s <= 0:
        interior_upper = 0
        total_upper = 1
    else:
        interior_upper = interior_depth_upper(first_s, lifetimes[0])
        # ``interior_upper`` counts further reference-surface transitions which may
        # remain interior.  One additional generated layer may be the layer whose
        # required registration window first reaches the absorbing t=0 surface.
        total_upper = interior_upper + 2
    if len(rows) > total_upper:
        raise AssertionError("signed-good generated epoch exceeded its finite physical-time depth bound")

    hits_boundary = common[-1] <= 0.0
    total_work = math.fsum(row.physical_hh_work_mass for row in rows)
    if not math.isfinite(total_work) or total_work <= 0.0:
        raise ValueError("generated epoch total physical HH work left the finite native range")
    return SignedGoodGeneratedEpochCertificate(
        layer_count=len(rows),
        transition_count=transitions,
        initial_parent_frequency=rows[0].parent_frequency,
        final_parent_frequency=rows[-1].parent_frequency,
        initial_parent_lifetime=lifetimes[0],
        first_common_reference_time=common[0],
        last_common_reference_time=common[-1],
        cumulative_reference_backshift=cumulative,
        minimum_cumulative_backshift=minimum,
        maximum_observed_parent_ratio=max_ratio,
        minimum_observed_lifetime_growth=min_growth,
        total_selected_physical_hh_work=total_work,
        interior_transition_upper=interior_upper,
        total_layer_upper_before_or_at_boundary=total_upper,
        hits_initial_boundary=hits_boundary,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "physical_input": "only an energy reentry which selects actual positive HH child-energy work may enter; the raw |I_HH| coefficient obstruction is an interval locator and is rejected",
        "signed_good_geometry": "each actual generated parent satisfies strict 3/5<N_parent/N_child<5/8, hence its natural lifetime grows by more than 64/25 backward",
        "physical_support": "the positive HH work law is restricted to an actual heavy half-slab; same-time parent-pair support then occupies at most 25/128 of the parent natural lifetime",
        "common_surface": "for H_j=[a_j,b_j], the s_j=a_j-(2/5)T_j are physical registration surfaces, not event vertices, recursive events, or checkpoint currency",
        "asynchronous_backshift": "H_(j+1) subset [s_j,b_j] and alpha_j<=10/39 imply s_j-s_(j+1)>=(1792/4875)T_j even when the next generated support begins anywhere in the previous common interval",
        "time_telescope": "Delta s_L>=(1792/7605)T_0[(64/25)^L-1], so at finite physical time a consecutive signed-good generated ancestry reaches the absorbing t=0 registration boundary after finite depth",
        "causal_weights": "all generation masses are actual positive child-energy HH work; raw Duhamel/adjoint amplitude weights are not used as a causal probability law",
        "reuse_relation": "Shannon/Renyi remains the native breadth/collision law of the generated ancestry; it is not needed to manufacture a scalar time cost for single-lineage depth",
        "master_consequence": "an infinite recursive event path cannot eventually consist only of consecutive signed-good generated HH renewals; generic HH/high-tail events not certified signed-good remain outside this theorem",
        "anti_count": "the finite-depth conclusion comes from physical parabolic lifetimes and the absorbing initial surface, not from charging one unit per event or from a normalized global reset",
        "scope": "this closes eventually-pure signed-good generated-HH recurrence only; mixed recurrence and generic non-signed-good HH/source/material/high-tail routes remain open, and no Navier-Stokes global-regularity claim is made",
    }


@dataclass(frozen=True)
class SignedGoodGeneratedEpochStress:
    samples: int
    minimum_scale_window_margin: float
    minimum_heavy_half_work_margin: float
    minimum_reference_backshift_margin: float
    minimum_cumulative_backshift_margin: float
    maximum_certified_layer_upper: int
    boundary_forced_samples: int
    raw_coefficient_reentry_rejections: int
    nonsigned_good_rejections: int
    invalid_support_restart_rejections: int
    minimum_sampled_frequency: float
    maximum_sampled_frequency: float
    minimum_sampled_natural_lifetime: float
    maximum_sampled_natural_lifetime: float
    minimum_sampled_physical_work: float
    maximum_sampled_physical_work: float


def _synthetic_reentry(
    work_lower: float,
    provenance: SignedGoodGeneratedWorkProvenance,
) -> dict[str, object]:
    return {
        "branch": ACTUAL_HH_GENERATION_BRANCH,
        "energy_gate": {
            "branch": ACTUAL_HH_GENERATION_BRANCH,
            "physical_hh_work_lower": work_lower,
            "provenance": provenance,
        },
        "provenance": provenance,
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
    }


def stress(samples: int = 50_000, seed: int = 20260811) -> SignedGoodGeneratedEpochStress:
    rng = random.Random(seed)
    mscale = mhalf = mback = mcum = math.inf
    maxupper = boundary = rawreject = ratio_reject = support_reject = 0
    min_frequency = min_lifetime = min_work = math.inf
    max_frequency = max_lifetime = max_work = 0.0

    for k in range(samples):
        # Sample in dimensionless logarithmic coordinates.  This exercises the
        # same theorem over more than 150 decades without inserting a unit-size
        # observer floor into frequency, time, or physical work.
        log_child = rng.uniform(-180.0, 180.0)
        log_child_lifetime = rng.uniform(-180.0, 180.0)
        child = math.exp(log_child)
        c = math.exp(log_child_lifetime + 2.0 * log_child)
        target = rng.randint(1, 7)
        first_ratio = rng.uniform(0.603, 0.622)
        first_parent = first_ratio * child
        first_parent_T = _parabolic_lifetime(c, first_parent)
        # Keep enough physical room for the requested interior layers.  The
        # theorem itself still checks the exact support/backshift geometry.
        end = 4.0 * first_parent_T * (2.7 ** target)
        rows: list[SignedGoodGeneratedHHStep] = []

        trajectory_id = f"stress-trajectory-{k}"
        for _j in range(target):
            ratio = rng.uniform(0.603, 0.622)
            parent = ratio * child
            Tchild = _parabolic_lifetime(c, child)
            width = rng.uniform(0.05, 0.45) * Tchild
            if rows:
                prev = rows[-1]
                end = min(end, prev.work_support_end)
                lower = prev.common_reference_time
                if end - width < lower:
                    end = lower + width
                if end > prev.work_support_end:
                    break
            start = max(0.0, end - width)
            total = math.exp(rng.uniform(-180.0, 180.0))
            work_lower = rng.uniform(0.3, 0.9) * total
            mass = rng.uniform(0.5, 0.9) * total
            provenance = SignedGoodGeneratedWorkProvenance(
                event_id=f"event-{k}-{_j}",
                trajectory_id=trajectory_id,
                child_carrier_id=f"carrier-{k}-{_j}",
                generated_parent_carrier_id=f"carrier-{k}-{_j + 1}",
                work_law_id=f"positive-hh-work-{k}-{_j}",
                child_frequency=child,
                parent_frequency=parent,
                scaled_lifetime=c,
                slab_start=end - Tchild,
                slab_end=end,
            )
            half = {
                "start": start,
                "end": end,
                "mass": mass,
                "total": total,
                "normalized_parent_span_upper": float(INITIAL_HALF_SPAN),
                "provenance": provenance,
            }
            step = signed_good_step_from_energy_reentry(
                reentry=_synthetic_reentry(work_lower, provenance),
                selected_physical_half_slab=half,
                child_frequency=child,
                parent_frequency=parent,
                scaled_lifetime=c,
            )
            rows.append(step)
            min_frequency = min(min_frequency, child, parent)
            max_frequency = max(max_frequency, child, parent)
            min_lifetime = min(
                min_lifetime,
                step.child_natural_lifetime,
                step.parent_natural_lifetime,
            )
            max_lifetime = max(
                max_lifetime,
                step.child_natural_lifetime,
                step.parent_natural_lifetime,
            )
            min_work = min(min_work, work_lower, mass, total)
            max_work = max(max_work, work_lower, mass, total)
            mscale = min(
                mscale,
                step.parent_child_ratio - float(PARENT_CHILD_LO),
                float(PARENT_CHILD_HI) - step.parent_child_ratio,
            )
            mhalf = min(mhalf, step.physical_hh_work_mass - 0.5 * step.physical_hh_work_total)
            child = parent
            end = step.work_support_end

        if not rows:
            continue
        out = signed_good_generated_epoch_telescope(rows)
        maxupper = max(maxupper, out.total_layer_upper_before_or_at_boundary)
        if out.hits_initial_boundary:
            boundary += 1
        if len(rows) > 1:
            for j in range(len(rows) - 1):
                actual = rows[j].common_reference_time - rows[j + 1].common_reference_time
                required = float(MIN_REFERENCE_BACKSTEP) * rows[j].parent_natural_lifetime
                mback = min(mback, actual - required)
        mcum = min(mcum, out.cumulative_reference_backshift - out.minimum_cumulative_backshift)

        # A raw coefficient obstruction cannot be smuggled directly into the epoch.
        if k % 7 == 0:
            try:
                signed_good_step_from_energy_reentry(
                    reentry={
                        "branch": HH_COEFFICIENT_OBSTRUCTION,
                        "requires_physical_energy_reentry": True,
                    },
                    selected_physical_half_slab={
                        "start": 0.0,
                        "end": 0.0,
                        "mass": 1.0,
                        "total": 1.0,
                        "normalized_parent_span_upper": float(INITIAL_HALF_SPAN),
                    },
                    child_frequency=2.0,
                    parent_frequency=1.22,
                    scaled_lifetime=1.0,
                )
            except TypeError:
                rawreject += 1
            else:
                raise AssertionError("raw HH coefficient obstruction crossed the actual-work epoch boundary")

        # Generic/non-signed-good HH generation remains outside this theorem.
        if k % 11 == 0:
            Tchild = 1.0 / 4.0
            provenance = SignedGoodGeneratedWorkProvenance(
                event_id=f"bad-ratio-event-{k}",
                trajectory_id=f"bad-ratio-trajectory-{k}",
                child_carrier_id="bad-ratio-child",
                generated_parent_carrier_id="bad-ratio-parent",
                work_law_id=f"bad-ratio-work-{k}",
                child_frequency=2.0,
                parent_frequency=1.4,
                scaled_lifetime=1.0,
                slab_start=2.0,
                slab_end=2.0 + Tchild,
            )
            try:
                signed_good_step_from_energy_reentry(
                    reentry=_synthetic_reentry(0.5, provenance),
                    selected_physical_half_slab={
                        "start": 2.0,
                        "end": 2.0 + 0.1 * Tchild,
                        "mass": 0.6,
                        "total": 1.0,
                        "normalized_parent_span_upper": float(INITIAL_HALF_SPAN),
                        "provenance": provenance,
                    },
                    child_frequency=2.0,
                    parent_frequency=1.4,
                    scaled_lifetime=1.0,
                )
            except ValueError:
                ratio_reject += 1
            else:
                raise AssertionError("non-signed-good HH event crossed the signed-good epoch theorem")

        # A next support outside the previous common interval is not the same
        # asynchronous generated lineage even if its scale ratio is signed-good.
        if len(rows) >= 2 and k % 13 == 0:
            bad = list(rows)
            row = bad[-1]
            shifted_time = bad[-2].work_support_end + 2.0 * bad[-2].parent_natural_lifetime
            shifted_provenance = SignedGoodGeneratedWorkProvenance(
                event_id=f"shifted-{row.provenance.event_id}",
                trajectory_id=row.provenance.trajectory_id,
                child_carrier_id=row.provenance.child_carrier_id,
                generated_parent_carrier_id=row.provenance.generated_parent_carrier_id,
                work_law_id=f"shifted-{row.provenance.work_law_id}",
                child_frequency=row.child_frequency,
                parent_frequency=row.parent_frequency,
                scaled_lifetime=row.scaled_lifetime,
                slab_start=shifted_time,
                slab_end=shifted_time + row.child_natural_lifetime,
            )
            shifted = SignedGoodGeneratedHHStep(
                child_frequency=row.child_frequency,
                parent_frequency=row.parent_frequency,
                scaled_lifetime=row.scaled_lifetime,
                work_support_start=shifted_time,
                work_support_end=shifted_time,
                physical_hh_work_mass=row.physical_hh_work_mass,
                physical_hh_work_total=row.physical_hh_work_total,
                physical_hh_work_lower=row.physical_hh_work_lower,
                provenance=shifted_provenance,
            )
            bad[-1] = shifted
            try:
                signed_good_generated_epoch_telescope(bad)
            except ValueError:
                support_reject += 1
            else:
                raise AssertionError("nonconsecutive generated support restart crossed the epoch telescope")

    if not math.isfinite(mback):
        mback = 0.0
    return SignedGoodGeneratedEpochStress(
        samples=samples,
        minimum_scale_window_margin=mscale,
        minimum_heavy_half_work_margin=mhalf,
        minimum_reference_backshift_margin=mback,
        minimum_cumulative_backshift_margin=mcum,
        maximum_certified_layer_upper=maxupper,
        boundary_forced_samples=boundary,
        raw_coefficient_reentry_rejections=rawreject,
        nonsigned_good_rejections=ratio_reject,
        invalid_support_restart_rejections=support_reject,
        minimum_sampled_frequency=min_frequency,
        maximum_sampled_frequency=max_frequency,
        minimum_sampled_natural_lifetime=min_lifetime,
        maximum_sampled_natural_lifetime=max_lifetime,
        minimum_sampled_physical_work=min_work,
        maximum_sampled_physical_work=max_work,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-signed-good-generated-epoch-time-telescope"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "signed_good_generated_epoch_time_telescope.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    md = f"""# Signed-good generated-HH epoch: parabolic physical-time telescope

Status: **{cert['status']}**.

A raw HH Duhamel coefficient obstruction is only an interval locator.  This theorem begins only after the same carrier has reentered its actual `Q^2` energy law and the physical-energy gate has selected positive high--high child-energy work.  Restrict that same positive work to a heavy half-slab and read the same-time hard parent pair there.

For every signed-good parent,

`3/5 < N_parent/N_child < 5/8`,

so its natural lifetime satisfies

`64/25 < T_parent/T_child < 25/9`.

If `H_j=[a_j,b_j]` is the physical generated parent-work support and `T_j` its parent natural lifetime, put

`s_j=a_j-(2/5)T_j`.

The half-slab geometry gives `|H_j|/T_j<=25/128<10/39`.  The next generated support may begin anywhere in the previous common interval; requiring only `H_(j+1) subset [s_j,b_j]`, the exact asynchronous theorem still gives

`s_j-s_(j+1) >= (1792/4875)T_j`.

Hence after `L` interior transitions,

`s_0-s_L >= (1792/7605)T_0[(64/25)^L-1]`.

The `s_j` are **registration surfaces**, not recursive events.  Because the right-hand side grows geometrically, at finite physical time the required backward registration surface reaches `t=0` after finite signed-good generated depth.  The adjoint/material-carrier gate is then truncated at the initial surface, which is absorbing.

No Duhamel amplitude law is used as causal probability; every mass here is actual positive HH child-energy work.  No unit event-count cost is charged.  Shannon/Renyi remains the separate native law for ancestry breadth/collision/reuse.

Stress: `{out.samples}` signed-good generated physical-work epochs
- minimum strict scale-window margin: `{out.minimum_scale_window_margin:.3e}`
- minimum heavy-half physical-work margin: `{out.minimum_heavy_half_work_margin:.3e}`
- minimum one-step reference-backshift margin: `{out.minimum_reference_backshift_margin:.3e}`
- minimum cumulative-backshift margin: `{out.minimum_cumulative_backshift_margin:.3e}`
- maximum certified layer upper bound sampled: `{out.maximum_certified_layer_upper}`
- epochs whose sampled terminal registration reached `t=0`: `{out.boundary_forced_samples}`
- raw coefficient-obstruction rejections: `{out.raw_coefficient_reentry_rejections}`
- non-signed-good HH rejections: `{out.nonsigned_good_rejections}`
- nonconsecutive support-restart rejections: `{out.invalid_support_restart_rejections}`
- native frequency range: `[{out.minimum_sampled_frequency:.3e}, {out.maximum_sampled_frequency:.3e}]`
- native lifetime range: `[{out.minimum_sampled_natural_lifetime:.3e}, {out.maximum_sampled_natural_lifetime:.3e}]`
- physical-work range: `[{out.minimum_sampled_physical_work:.3e}, {out.maximum_sampled_physical_work:.3e}]`

Master consequence: a recursive path cannot eventually consist only of consecutive **signed-good** physical HH generation.  Generic HH/high-tail generation without this physical signed-good geometry remains a genuine separate route.  Mixed-owner recurrence remains open, and no Navier--Stokes global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
