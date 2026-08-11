from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, Sequence

from src.high_strain_resolved_ancestor import TRANSPORTER_RADIUS
from src.full_natural_checkpoint_quotient import FULL_NATURAL_CHECKPOINT, FullNaturalCheckpoint
from src.full_natural_service_corridor_quotient import FULL_NATURAL_SERVICE_WITNESS
from src.same_carrier_checkpoint_segmentation_quotient import (
    SameCarrierCheckpointPathCertificate,
    SameCarrierMonitorSegment,
    SameCarrierProvenance,
    checkpoint_continuation_policy,
)
from src.common_slice_coefficient_registration import (
    HH_COEFFICIENT_OBSTRUCTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)


STATUS = (
    "EXACT_CONTINUUM_MASTER_EVENT_QUOTIENT__ZERO_CHARGE_RELAYS_COLLAPSED__"
    "NATIVE_PHYSICAL_TIME_RECURSION__SUPPLIER_SPECIFIC_SCALE_PROGRESS__"
    "NATURAL_HORIZON_CHECKPOINTS_ZERO_EVENT_DEPTH__"
    "SAME_CARRIER_FIRST_HIT_NOT_RESET_BY_CHECKPOINTS__"
    "CHECKPOINT_ZENO_NOT_CANONICAL_LINEAGE__NO_COMMON_CLOCK_OR_CAUSAL_REWEIGHTING"
)


COEFFICIENT_OBSTRUCTION_LABELS = frozenset(
    {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
)

NON_EVENT_CORRIDOR_WITNESS_LABELS = frozenset({FULL_NATURAL_SERVICE_WITNESS})
NON_EVENT_CHECKPOINT_LABELS = frozenset({FULL_NATURAL_CHECKPOINT})


def require_routed_physical_owner_labels(labels: Iterable[str]) -> tuple[str, ...]:
    """Reject Duhamel coefficient locators at the canonical physical-owner boundary.

    Coefficient obstructions are measurable first-stop *locators*.  They become
    causal owners only after actual Q^2 carrier energy/work has passed through a
    physical-energy reentry theorem.  Compatibility fields named ``causes`` are
    therefore not sufficient evidence for admission to the master state.
    """
    out = tuple(str(x) for x in labels if str(x))
    bad = tuple(sorted(set(out).intersection(COEFFICIENT_OBSTRUCTION_LABELS)))
    if bad:
        raise TypeError(
            "unrouted coefficient obstruction cannot enter the canonical physical owner state; resolve actual Q^2 energy/work first: "
            + ", ".join(bad)
        )
    witness_bad = tuple(sorted(set(out).intersection(NON_EVENT_CORRIDOR_WITNESS_LABELS)))
    if witness_bad:
        raise TypeError(
            "full-natural own-scale service is a same-corridor physical witness, not a separate recursive owner event: "
            + ", ".join(witness_bad)
        )
    checkpoint_bad = tuple(sorted(set(out).intersection(NON_EVENT_CHECKPOINT_LABELS)))
    if checkpoint_bad:
        raise TypeError(
            "full-natural horizon checkpoint is analysis re-registration after real corridor time, not a recursive physical owner event: "
            + ", ".join(checkpoint_bad)
        )
    return out


def owner_bundle_from_energy_reentry(
    physical_measure: str,
    mass: float,
    reentry: dict[str, object],
) -> "PhysicalOwnerBundle":
    """Typed adapter from an actual energy/work reentry result to master ownership.

    The adapter refuses raw obstruction impulses.  For the smooth-interface branch
    it uses only the gauge-quotiented physical relink/strain owners returned after
    the energy gate.  Other gate branches use their physical branch name directly.
    """
    if bool(reentry.get("coefficient_impulse_used_as_physical_work", True)):
        raise TypeError("energy reentry must certify that coefficient impulse magnitude was not used as physical work")
    if bool(reentry.get("observer_partition_motion_charged_as_physics", True)):
        raise TypeError("energy reentry must quotient observer partition motion before physical ownership")
    branch = str(reentry.get("branch", ""))
    if branch == "smooth_interface_physical_work":
        owners = tuple(str(x) for x in reentry.get("joint_interface_owners", ()))
        if not owners:
            raise TypeError("smooth-interface energy reentry supplied no physical relink/strain owner")
    elif branch in {
        "material_energy_inheritance",
        "high_strain_critical_dissipation",
        "physical_high_high_transfer_generation",
    }:
        owners = (branch,)
    else:
        raise TypeError("unrecognized or unresolved physical-energy reentry branch")
    owners = require_routed_physical_owner_labels(owners)
    return canonical_owner_bundle(physical_measure, mass, owners)


class SupplierKind(str, Enum):
    """Physical supplier geometry, never theorem priority."""

    GENERATED_SIGNED_GOOD_HH = "generated_signed_good_hh_parent"
    RESOLVED_DISSIPATION = "resolved_dissipation_ancestor"
    PRESSURE_PAIR = "resolved_pressure_pair_shell"
    FRESH_SGS_SCALE = "fresh_sgs_scale_shell"
    HIGH_TAIL = "high_tail_hard_shell"
    GENERIC_CRITICAL_SHELL = "generic_critical_shell_corridor"
    MATERIAL_REUSE = "material_reuse_or_service"
    HH_REGENERATION = "hh_regeneration_owner"


class EventDisposition(str, Enum):
    NAMED_RECURSIVE_PHYSICAL_EVENT = "named_recursive_physical_event"
    FULL_NATURAL_SURVIVOR = "full_natural_survivor"  # legacy compatibility only; rejected by RecursiveEventState
    ABSORBING_INITIAL_BOUNDARY = "absorbing_initial_boundary"
    TERMINAL_COST = "terminal_multiplicative_or_global_resource_cost"


@dataclass(frozen=True)
class RecursiveEventState:
    """Canonical master state after representation and relay quotient.

    `time` is physical time. `frequency` is an actual hard-shell/carrier scale
    already supplied by the PDE-facing theorem.  There is deliberately no step
    index, normalized common clock, coherent-cell label, or theorem priority.
    Material/diagnostic information may be retained as sidecars, but it cannot
    change the carrier state by itself.
    """

    time: float
    frequency: float
    physical_measure: str
    joint_causes: tuple[str, ...] = ()
    disposition: EventDisposition = EventDisposition.NAMED_RECURSIVE_PHYSICAL_EVENT
    sidecars: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not math.isfinite(self.time) or self.time < 0:
            raise ValueError("physical event time must be finite and nonnegative")
        if not math.isfinite(self.frequency) or self.frequency <= 0:
            raise ValueError("physical event frequency must be positive and finite")
        if not self.physical_measure:
            raise ValueError("a named physical observable/measure is required")
        require_routed_physical_owner_labels(self.joint_causes)
        if len(set(self.joint_causes)) != len(self.joint_causes):
            raise ValueError("joint physical cause set must already be quotiented")
        if self.disposition is EventDisposition.FULL_NATURAL_SURVIVOR:
            raise TypeError("full-natural no-hit horizon is an analysis checkpoint, not a RecursiveEventState")
        if self.disposition is EventDisposition.ABSORBING_INITIAL_BOUNDARY and self.time != 0.0:
            raise ValueError("the absorbing initial boundary is exactly t=0")

    @property
    def absorbing(self) -> bool:
        return self.disposition is EventDisposition.ABSORBING_INITIAL_BOUNDARY


@dataclass(frozen=True)
class PhysicalOwnerBundle:
    """One positive physical measure with a set-valued owner/provenance mark.

    Exact ties and downstream theorem manifestations may add owner names, but the
    mass of the underlying measure occurs once.  No owner fractions are created.
    """

    physical_measure: str
    mass: float
    owners: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.physical_measure:
            raise ValueError("named physical measure required")
        if not math.isfinite(self.mass) or self.mass < 0:
            raise ValueError("finite nonnegative physical mass required")
        if not self.owners:
            raise ValueError("at least one physical owner required")
        require_routed_physical_owner_labels(self.owners)
        if tuple(sorted(set(self.owners))) != self.owners:
            raise ValueError("owners must be a sorted quotiented set")


def canonical_owner_bundle(
    physical_measure: str,
    mass: float,
    owners: Iterable[str],
) -> PhysicalOwnerBundle:
    """Quotient duplicate owner manifestations without splitting physical mass."""
    canonical = tuple(sorted({str(x) for x in owners if str(x)}))
    return PhysicalOwnerBundle(str(physical_measure), float(mass), canonical)


def zero_charge_owner_relay(bundle: PhysicalOwnerBundle, downstream_owners: Iterable[str]) -> PhysicalOwnerBundle:
    """Attach downstream supplier provenance to the *same* physical measure.

    This is legal only as a same-measure relay: it changes no mass and creates no
    causal entropy.  A theorem that creates a genuinely new positive physical law
    must construct a new bundle instead of using this function.
    """
    return canonical_owner_bundle(bundle.physical_measure, bundle.mass, (*bundle.owners, *tuple(downstream_owners)))


@dataclass(frozen=True)
class WitnessRelay:
    """A certified consequence map which creates state, not a second causal charge.

    Unlike ``zero_charge_owner_relay``, the downstream observable may have
    different units (for example pressure-pair work -> critical shell mass).
    The relay therefore asserts no equality of measures or masses.  It records
    only that the downstream state is certified from the upstream physical law
    and must not be charged again as an independent cause.
    """

    upstream_physical_measure: str
    downstream_state_observable: str
    certification_relation: str
    causal_charge_created: bool = False
    diagnostic_probability_created: bool = False

    def __post_init__(self) -> None:
        if not self.upstream_physical_measure or not self.downstream_state_observable or not self.certification_relation:
            raise ValueError("named upstream law, downstream state and certification relation required")
        if self.causal_charge_created or self.diagnostic_probability_created:
            raise ValueError("a witness relay may create neither a second causal charge nor a diagnostic probability law")


def zero_charge_witness_relay(upstream_physical_measure: str, downstream_state_observable: str, certification_relation: str) -> WitnessRelay:
    return WitnessRelay(upstream_physical_measure, downstream_state_observable, certification_relation)


@dataclass(frozen=True)
class SupplierScaleCertificate:
    supplier: SupplierKind
    ratio_lower: float | None
    ratio_upper: float | None
    lower_strict: bool
    upper_strict: bool
    directional_progress: str
    physical_meaning: str


def supplier_scale_certificate(supplier: SupplierKind) -> SupplierScaleCertificate:
    """Certified next-shell/current-scale geometry, and nothing stronger."""
    if supplier is SupplierKind.GENERATED_SIGNED_GOOD_HH:
        return SupplierScaleCertificate(
            supplier, 3.0 / 5.0, 5.0 / 8.0, True, True, "strictly_to_lower_frequency",
            "actual signed-good HH parent role; T_parent/T_child is in (64/25,25/9)",
        )
    if supplier is SupplierKind.RESOLVED_DISSIPATION:
        return SupplierScaleCertificate(
            supplier, None, TRANSPORTER_RADIUS, False, False, "to_lower_frequency",
            "actual resolved D_V shell satisfies M<=N/4; natural lifetime is at least 16 times longer",
        )
    if supplier is SupplierKind.PRESSURE_PAIR:
        return SupplierScaleCertificate(
            supplier, None, 1.0 / 4.0, False, False, "to_lower_frequency",
            "actual resolved pressure pair lies inside V=S_(N/4)u and supplies a u hard shell M<=N/4",
        )
    if supplier is SupplierKind.FRESH_SGS_SCALE:
        return SupplierScaleCertificate(
            supplier, None, 2.0, False, False, "no_directional_progress_supplied",
            "refinement-invariant low/base LP pushforward gives one of two hard shells with M<=2N; no lower ratio is certified",
        )
    if supplier is SupplierKind.HIGH_TAIL:
        return SupplierScaleCertificate(
            supplier, 2.0, None, False, False, "to_higher_frequency",
            "hard-tail support gives M/N>=2 and T_M/T_N<=1/4",
        )
    if supplier is SupplierKind.GENERIC_CRITICAL_SHELL:
        return SupplierScaleCertificate(
            supplier, None, None, False, False, "no_supplier_relative_scale_progress",
            "generic shell theorem gives named first stop / t=0 / full own-scale service; A=3M/4 is analysis registration, not master scale progress",
        )
    if supplier is SupplierKind.MATERIAL_REUSE:
        return SupplierScaleCertificate(
            supplier, None, None, False, False, "no_scale_progress_without_external_epoch_geometry",
            "full-natural own-scale service is a same-corridor witness; genuinely independent material reuse or source-service ownership is physical ancestry routing. No scale progress is inferred without a separate PDE theorem",
        )
    if supplier is SupplierKind.HH_REGENERATION:
        return SupplierScaleCertificate(
            supplier, None, None, False, False, "no_scale_progress_until_physical_HH_owner_is_resolved",
            "coefficient regeneration is only an earlier-generation obstruction until re-entered through actual positive HH work",
        )
    raise ValueError(f"unknown supplier {supplier}")


def validate_supplier_scale(supplier: SupplierKind, current_frequency: float, next_frequency: float) -> dict[str, object]:
    N = float(current_frequency)
    M = float(next_frequency)
    if min(N, M) <= 0 or not math.isfinite(N + M):
        raise ValueError("positive finite current/next frequencies required")
    cert = supplier_scale_certificate(supplier)
    ratio = M / N
    tol = 8e-13 * max(1.0, abs(ratio))
    if cert.ratio_lower is not None:
        if cert.lower_strict:
            if not ratio > cert.ratio_lower + tol:
                raise ValueError("supplier violates strict certified lower scale ratio")
        elif ratio + tol < cert.ratio_lower:
            raise ValueError("supplier violates certified lower scale ratio")
    if cert.ratio_upper is not None:
        if cert.upper_strict:
            if not ratio < cert.ratio_upper - tol:
                raise ValueError("supplier violates strict certified upper scale ratio")
        elif ratio > cert.ratio_upper + tol:
            raise ValueError("supplier violates certified upper scale ratio")
    return {
        "supplier": supplier.value,
        "actual_frequency_ratio": ratio,
        "directional_progress": cert.directional_progress,
        "physical_meaning": cert.physical_meaning,
    }


def natural_duration(frequency: float, scaled_lifetime: float) -> float:
    M = float(frequency)
    c = float(scaled_lifetime)
    if min(M, c) <= 0 or not math.isfinite(M + c):
        raise ValueError("positive finite frequency and scaled lifetime required")
    return c / (M * M)


def full_natural_survivor_endpoint(event_time: float, frequency: float, scaled_lifetime: float) -> dict[str, float | str | bool]:
    """One no-hit corridor uses its own physical natural time; an interior horizon returns an analysis checkpoint, not an event."""
    t = float(event_time)
    if t < 0 or not math.isfinite(t):
        raise ValueError("finite nonnegative event time required")
    T = natural_duration(frequency, scaled_lifetime)
    if T >= t:
        return {
            "start_time": t,
            "requested_duration": T,
            "end_time": 0.0,
            "actual_time_drop": t,
            "hits_initial_boundary": True,
            "disposition": EventDisposition.ABSORBING_INITIAL_BOUNDARY.value,
        }
    return {
        "start_time": t,
        "requested_duration": T,
        "end_time": t - T,
        "actual_time_drop": T,
        "hits_initial_boundary": False,
        "disposition": FULL_NATURAL_CHECKPOINT,
        "checkpoint_kind": FULL_NATURAL_CHECKPOINT,
        "physical_event_created": False,
        "recursion_edges_added": 0,
    }


def physical_time_telescope(times: Sequence[float]) -> dict[str, float]:
    """Exact physical-time identity for any ordered event or checkpoint times."""
    t = tuple(float(x) for x in times)
    if len(t) < 1 or any((not math.isfinite(x) or x < 0) for x in t):
        raise ValueError("nonempty finite nonnegative physical times required")
    if any(t[j + 1] > t[j] for j in range(len(t) - 1)):
        raise ValueError("backward recursive event times must be nonincreasing")
    drops = sum(t[j] - t[j + 1] for j in range(len(t) - 1))
    endpoint = t[0] - t[-1]
    return {
        "sum_physical_time_drops": drops,
        "endpoint_time_drop": endpoint,
        "residual": drops - endpoint,
    }


def log_scale_telescope(frequencies: Sequence[float]) -> dict[str, float]:
    """Exact coordinate identity; it does not assert monotone scale progress."""
    f = tuple(float(x) for x in frequencies)
    if len(f) < 1 or any((not math.isfinite(x) or x <= 0) for x in f):
        raise ValueError("nonempty positive finite frequencies required")
    increments = sum(math.log(f[j + 1] / f[j]) for j in range(len(f) - 1))
    endpoint = math.log(f[-1] / f[0])
    return {
        "sum_log_scale_changes": increments,
        "endpoint_log_scale_change": endpoint,
        "residual": increments - endpoint,
    }


def bounded_scale_full_survivor_steps_to_boundary(
    initial_time: float,
    frequency_upper: float,
    scaled_lifetime: float,
) -> int:
    """Number of bounded-scale no-hit checkpoint corridors sufficient to force t=0.

    If every survivor scale M_j<=Mbar then every requested duration is at least
    c/Mbar^2. Therefore K=ceil(t0 Mbar^2/c) such checkpoint corridors cannot all remain interior.
    """
    t0 = float(initial_time)
    M = float(frequency_upper)
    c = float(scaled_lifetime)
    if t0 < 0 or M <= 0 or c <= 0 or not all(math.isfinite(x) for x in (t0, M, c)):
        raise ValueError("finite t0>=0 and positive scale bound/lifetime required")
    if t0 == 0:
        return 0
    return int(math.ceil(t0 * M * M / c))


def trace_full_natural_survivors(
    initial_time: float,
    survivor_frequencies: Sequence[float],
    scaled_lifetime: float,
) -> dict[str, object]:
    """Compatibility helper for a finite prefix of no-hit physical corridors; interior endpoints are checkpoints."""
    t = float(initial_time)
    if t < 0 or not math.isfinite(t):
        raise ValueError("finite nonnegative initial time required")
    times = [t]
    used = 0
    boundary = t == 0.0
    for frequency in survivor_frequencies:
        if boundary:
            break
        out = full_natural_survivor_endpoint(t, float(frequency), scaled_lifetime)
        t = float(out["end_time"])
        times.append(t)
        used += 1
        boundary = bool(out["hits_initial_boundary"])
    tel = physical_time_telescope(times)
    return {
        "times": tuple(times),
        "survivors_used": used,
        "checkpoints_used": used,
        "recursive_events_added": 0,
        "hits_initial_boundary": boundary,
        "final_time": t,
        "physical_time_telescope_residual": float(tel["residual"]),
    }


def geometric_uv_natural_time_sum(initial_frequency: float, scaled_lifetime: float, scale_ratio: float = 2.0) -> float:
    """Total time of independently certified scales M_j=M0*r^j, r>1.

    The finite sum warns that time alone cannot terminate a genuinely
    changing-scale physical path. It is not a duration assignment for arbitrary
    observer cuts on one fixed carrier.
    """
    M = float(initial_frequency)
    c = float(scaled_lifetime)
    r = float(scale_ratio)
    if M <= 0 or c <= 0 or r <= 1 or not all(math.isfinite(x) for x in (M, c, r)):
        raise ValueError("positive finite M,c and r>1 required")
    return (c / (M * M)) / (1.0 - r ** -2)


@dataclass(frozen=True)
class NativePathLedger:
    """Typed product ledger: no exchange rate between unlike physical objects."""

    physical_time_drop: float = 0.0
    log_scale_change: float = 0.0
    multiplicative_transfer_cost: float = 0.0
    causal_reuse_action: float = 0.0
    xi_penalty: float = 0.0
    global_resources: tuple[tuple[str, float], ...] = ()

    def __post_init__(self) -> None:
        nonnegative = (
            self.physical_time_drop,
            self.multiplicative_transfer_cost,
            self.causal_reuse_action,
            self.xi_penalty,
        )
        if any((not math.isfinite(x) or x < 0) for x in nonnegative):
            raise ValueError("finite nonnegative physical ledger actions required")
        if not math.isfinite(self.log_scale_change):
            raise ValueError("finite signed log-scale change required")
        names = [k for k, _ in self.global_resources]
        if names != sorted(set(names)):
            raise ValueError("global resources must be a sorted unique tuple")
        if any((not k or not math.isfinite(v) or v < 0) for k, v in self.global_resources):
            raise ValueError("valid nonnegative global resource consumptions required")


def compose_native_ledgers(increments: Sequence[NativePathLedger]) -> NativePathLedger:
    resources: dict[str, float] = {}
    for inc in increments:
        for name, value in inc.global_resources:
            resources[name] = resources.get(name, 0.0) + value
    return NativePathLedger(
        physical_time_drop=sum(x.physical_time_drop for x in increments),
        log_scale_change=sum(x.log_scale_change for x in increments),
        multiplicative_transfer_cost=sum(x.multiplicative_transfer_cost for x in increments),
        causal_reuse_action=sum(x.causal_reuse_action for x in increments),
        xi_penalty=sum(x.xi_penalty for x in increments),
        global_resources=tuple(sorted(resources.items())),
    )


@dataclass(frozen=True)
class DiagnosticConcentration:
    """Noncausal coordinate attached to a positive law; never a master action."""

    name: str
    value: float
    underlying_measure: str

    def __post_init__(self) -> None:
        if not self.name or not self.underlying_measure or not math.isfinite(self.value):
            raise ValueError("finite named diagnostic coordinate and underlying measure required")


def forbid_diagnostic_as_causal_action(_: DiagnosticConcentration) -> None:
    raise TypeError(
        "diagnostic/concentration coordinates (fresh-scale H_inf/H2, high-tail H_inf, etc.) cannot be charged as causal entropy or a global reset"
    )


def master_escape_dichotomy() -> dict[str, str]:
    """Analytic infinite-event consequence after checkpoint segmentation is quotiented."""
    return {
        "statement": (
            "after zero-charge relays, same-corridor service witnesses, natural-horizon checkpoints, and same-carrier checkpoint segmentation are quotiented, any infinite recursive EVENT path avoiding t=0 must contain infinitely many genuine named non-free physical owner events"
        ),
        "proof": (
            "a no-event natural horizon cannot replace the event-anchored smooth carrier or reset its cumulative strain/coefficient first-hit monitors. Genuine fixed-A,c natural windows have one positive duration and cannot accumulate before t=0. Arbitrary observer cuts carry no window duration; a checkpoint must be bound to an actual no-hit restriction of the exact event/carrier/dual/PDE path, whose gluing leaves the same first stop, and at an interior cut accumulation the actual no-earlier-hit path plus a matching smooth-extension token gives either an endpoint face or continuation of the same carrier"
        ),
        "remaining_physics": (
            "the remaining infinite-path problem is recurrence of genuine physical owner events and must telescope through their native work/source/service/reuse laws. Actual high-tail dissipation/work remains a genuine event route when its independent physical hypotheses are met; checkpoint re-hardening cannot manufacture that provenance"
        ),
    }


def theorem_certificate() -> dict[str, object]:
    scale = {kind.value: asdict(supplier_scale_certificate(kind)) for kind in SupplierKind}
    return {
        "status": STATUS,
        "canonical_state": (
            "one physical event time, one supplied physical shell/carrier frequency, one named physical measure, an unsplit joint physical-owner set after any required coefficient-energy reentry, and optional sidecars; no theorem-depth counter or common normalized clock"
        ),
        "relay_quotient": (
            "same-law owner relays preserve one unsplit physical mass; certified witness relays may change units/observable (for example source work to shell mass) but create no second causal charge. Inserting/removing either relay cannot manufacture causal entropy"
        ),
        "joint_owners": (
            "exact ties and multiple certified downstream owner names remain a set-valued provenance mark on one unsplit physical mass; no lexicographic priority and no heterogeneous RN tie normalization"
        ),
        "coefficient_obstruction_barrier": (
            "Duhamel HH/interface coefficient threshold hits are first-stop locators, not physical owners; the canonical event state and owner bundle reject those labels until actual Q^2 energy/work reentry returns a physical inheritance, high-strain, HH-work, relink or strain owner"
        ),
        "full_natural_service_barrier": (
            "full_natural_own_scale_service is a positive witness carried by the already-completed natural corridor, not a second recursive owner event; canonical owner states reject this classification label as an owner"
        ),
        "full_natural_checkpoint_barrier": (
            "a complete no-hit natural horizon consumes actual physical corridor time but is only an analysis checkpoint; RecursiveEventState rejects both the checkpoint label and the legacy full-natural-survivor disposition"
        ),
        "same_carrier_checkpoint_segmentation": (
            "a no-event checkpoint cannot replace the event-anchored smooth carrier or reset its terminal coefficient, cumulative strain action, or cumulative complex coefficient-impulse monitors; a bare checkpoint is insufficient, because an actual no-hit cumulative restriction ending at its exact native duration and matching the expected event/carrier/terminal-dual/PDE trajectory is required; hard-shell readings remain sidecars until a new physical event hardens a role"
        ),
        "universal_time_identity": "sum_j (t_j-t_(j+1)) = t_0-t_L on any ordered physical event or checkpoint times; event counting is a separate ontology",
        "natural_survivor": "a no-hit full-natural corridor consumes its theorem-supplied physical duration, but its horizon endpoint is an analysis checkpoint unless a first stop or t=0 occurs",
        "compact_scale_no_escape": (
            "for the same fixed carrier every genuine cA^-2 window has positive fixed duration and reaches t=0 after finitely many complete windows; arbitrary checkpoint cuts are duration-free restrictions and checkpoint scale readings are noncanonical"
        ),
        "infinite_escape_dichotomy": master_escape_dichotomy()["statement"],
        "uv_obstruction": (
            "time alone still allows an independently certified geometrically changing-scale physical path to have finite parabolic duration; that anti-theorem cannot attach natural-window durations to checkpoint state readings, which are not a canonical physical lineage, and genuine UV progression requires physical tail work/dissipation or another event"
        ),
        "scale_progress": scale,
        "bellman_coordinate": (
            "there is no canonical scalar exchange rate between log scale, physical work, service/reuse and global resources; the natural master ledger is typed/direct-product. Physical time and actual log shell scale telescope kinematically, while transfer cost, causal reuse and each genuinely global resource telescope only on their own physical laws"
        ),
        "service_semantics": (
            "own-scale service produced by a completed full-natural shell corridor is a same-interval physical witness and adds zero recursion depth; material rereading of that same law is also a witness relay. Independent source/service/reuse events remain physical owners. No service observable is promoted to an additive globally bounded currency"
        ),
        "diagnostic_separation": (
            "fresh-scale H_inf/H2 and high-tail scale/time concentration coordinates remain conjugate lower-bound diagnostics; they are not causal Shannon/Renyi action and cannot be inserted into the causal ledger"
        ),
        "boundary": "t=0 is absorbing",
        "scope": (
            "this is a continuum master assembly/quotient theorem, not a global no-escape or Navier-Stokes regularity proof; checkpoint segmentation is no longer an independent escape seam, while termination of infinitely recurring genuine physical owners remains open"
        ),
    }


@dataclass(frozen=True)
class QuotientStress:
    samples: int
    worst_owner_mass_residual: float
    worst_time_telescope_residual: float
    worst_scale_telescope_residual: float
    bounded_scale_boundary_failures: int
    supplier_scale_failures: int
    coefficient_obstruction_barrier_failures: int
    service_witness_barrier_failures: int
    checkpoint_barrier_failures: int
    checkpoint_segmentation_barrier_failures: int
    maximum_relay_owner_count: int
    minimum_uv_time_gap_to_naive_infinite_sum: float


def stress(samples: int = 50_000, seed: int = 20260810) -> QuotientStress:
    rng = random.Random(seed)
    wom = wtime = wscale = 0.0
    boundary_fail = supplier_fail = obstruction_fail = service_witness_fail = checkpoint_fail = segmentation_fail = 0
    maxowners = 0
    uv_gap = math.inf

    suppliers_with_sampled_geometry = (
        SupplierKind.GENERATED_SIGNED_GOOD_HH,
        SupplierKind.RESOLVED_DISSIPATION,
        SupplierKind.PRESSURE_PAIR,
        SupplierKind.FRESH_SGS_SCALE,
        SupplierKind.HIGH_TAIL,
    )

    try:
        canonical_owner_bundle(
            "raw coefficient locator",
            1.0,
            (ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,),
        )
    except TypeError:
        pass
    else:
        obstruction_fail += 1
        raise AssertionError("raw coefficient obstruction crossed the canonical physical-owner boundary")

    try:
        canonical_owner_bundle(
            "completed natural-corridor service witness",
            1.0,
            (FULL_NATURAL_SERVICE_WITNESS,),
        )
    except TypeError:
        pass
    else:
        service_witness_fail += 1
        raise AssertionError("full-natural service witness crossed the canonical recursive-owner boundary")

    try:
        canonical_owner_bundle(
            "full-natural checkpoint rereading",
            1.0,
            (FULL_NATURAL_CHECKPOINT,),
        )
    except TypeError:
        pass
    else:
        checkpoint_fail += 1
        raise AssertionError("analysis checkpoint crossed the canonical physical-owner boundary")

    try:
        RecursiveEventState(
            0.5,
            4.0,
            "no-hit natural horizon",
            (),
            EventDisposition.FULL_NATURAL_SURVIVOR,
        )
    except TypeError:
        pass
    else:
        checkpoint_fail += 1
        raise AssertionError("legacy full-natural survivor disposition crossed into RecursiveEventState")

    forged_checkpoint_record = {
        "checkpoint_kind": FULL_NATURAL_CHECKPOINT,
        "physical_event_created": False,
        "causal_charge_created": False,
        "recursion_edges_added": 0,
    }
    try:
        checkpoint_continuation_policy(forged_checkpoint_record)
    except TypeError:
        pass
    else:
        segmentation_fail += 1
        raise AssertionError("dictionary flags forged same-carrier checkpoint authority")

    checkpoint_record = FullNaturalCheckpoint(
        terminal_time=2.0,
        physical_time_drop=0.25,
        parent_shell_frequency=8.0 / 3.0,
        parent_shell_critical_mass_lower=2.0,
        corridor_frequency=2.0,
        scaled_lifetime=1.0,
        endpoint_carrier_critical_mass_lower=2.0,
        endpoint_shell_candidates=(2.0, 4.0),
    )
    provenance = SameCarrierProvenance(
        event_id="master-stress-event",
        carrier_id="master-stress-Q",
        terminal_dual_id="master-stress-dual",
        trajectory_id="master-stress-NS-trajectory",
        terminal_state_token="master-stress-terminal-state",
        terminal_time=2.0,
        carrier_frequency=2.0,
        scaled_lifetime=1.0,
        terminal_coefficient=1.0 + 0.0j,
    )
    checkpoint_path = SameCarrierCheckpointPathCertificate(
        checkpoint_record,
        (
            SameCarrierMonitorSegment(
                provenance=provenance,
                state_tokens=(provenance.terminal_state_token, "master-stress-checkpoint-state"),
                elapsed_times=(0.0, checkpoint_record.physical_time_drop),
                strain_action=(0.0, 0.0),
                residual_impulse=(0.0j, 0.0j),
                hh_impulse=(0.0j, 0.0j),
            ),
        ),
    )
    policy = checkpoint_continuation_policy(checkpoint_path, provenance=provenance)
    if bool(policy["carrier_replacement_authorized"]) or bool(policy["monitor_reset_authorized"]):
        segmentation_fail += 1
        raise AssertionError("full-natural checkpoint authorized a same-carrier restart/reset")
    for request in (
        {"request_carrier_replacement": True},
        {"request_terminal_amplitude_reset": True},
        {"request_monitor_reset": True},
    ):
        try:
            checkpoint_continuation_policy(checkpoint_path, provenance=provenance, **request)
        except TypeError:
            pass
        else:
            segmentation_fail += 1
            raise AssertionError("checkpoint restart/reset request crossed the same-carrier master barrier")

    routed = owner_bundle_from_energy_reentry(
        "actual positive q^2-weighted HH work",
        1.0,
        {
            "branch": "physical_high_high_transfer_generation",
            "coefficient_impulse_used_as_physical_work": False,
            "observer_partition_motion_charged_as_physics": False,
        },
    )
    if routed.owners != ("physical_high_high_transfer_generation",):
        obstruction_fail += 1
        raise AssertionError("actual physical-energy reentry failed to enter canonical ownership")

    for _ in range(samples):
        mass = math.exp(rng.uniform(-14.0, 5.0))
        bundle = canonical_owner_bundle("actual_positive_physical_law", mass, ("source", "source"))
        nrelay = rng.randint(1, 7)
        for j in range(nrelay):
            bundle = zero_charge_owner_relay(bundle, (f"relay_owner_{j}", "source"))
        maxowners = max(maxowners, len(bundle.owners))
        wom = max(wom, abs(bundle.mass - mass))
        if bundle.mass != mass:
            raise AssertionError("zero-charge relay changed physical mass")

        # Physical time is the only universal recursion clock and telescopes exactly.
        t0 = math.exp(rng.uniform(-5.0, 2.0))
        drops = [rng.random() for _j in range(rng.randint(1, 8))]
        s = sum(drops)
        drops = [d * (0.8 * t0 / max(s, 1e-300)) for d in drops]
        times = [t0]
        for d in drops:
            times.append(times[-1] - d)
        tout = physical_time_telescope(times)
        wtime = max(wtime, abs(float(tout["residual"])))
        if abs(float(tout["residual"])) > 2e-12 * max(1.0, t0):
            raise AssertionError("physical event-time telescope failed")

        freqs = [math.exp(rng.uniform(-4.0, 4.0))]
        for _j in range(rng.randint(1, 8)):
            freqs.append(freqs[-1] * math.exp(rng.uniform(-1.0, 1.0)))
        sout = log_scale_telescope(freqs)
        wscale = max(wscale, abs(float(sout["residual"])))
        if abs(float(sout["residual"])) > 3e-12 * max(1.0, abs(float(sout["endpoint_log_scale_change"]))):
            raise AssertionError("physical log-scale telescope failed")

        # A bounded-scale tail cannot contain arbitrarily many free full corridors.
        c = math.exp(rng.uniform(-3.0, 1.0))
        Mbar = math.exp(rng.uniform(-2.0, 3.0))
        K = bounded_scale_full_survivor_steps_to_boundary(t0, Mbar, c)
        if K > 0:
            ff = [Mbar * (0.2 + 0.8 * rng.random()) for _j in range(K)]
            traced = trace_full_natural_survivors(t0, ff, c)
            if not bool(traced["hits_initial_boundary"]):
                boundary_fail += 1
                raise AssertionError("bounded-scale full-survivor checkpoint continuation escaped t=0")
            if int(traced["recursive_events_added"]) != 0:
                checkpoint_fail += 1
                raise AssertionError("full-natural checkpoints manufactured recursive event depth")

        kind = rng.choice(suppliers_with_sampled_geometry)
        N = math.exp(rng.uniform(-3.0, 3.0))
        if kind is SupplierKind.GENERATED_SIGNED_GOOD_HH:
            ratio = rng.uniform(0.602, 0.623)
        elif kind in (SupplierKind.RESOLVED_DISSIPATION, SupplierKind.PRESSURE_PAIR):
            ratio = rng.uniform(0.01, 0.25)
        elif kind is SupplierKind.FRESH_SGS_SCALE:
            ratio = math.exp(rng.uniform(-5.0, math.log(2.0)))
        else:
            ratio = math.exp(rng.uniform(math.log(2.0), math.log(64.0)))
        try:
            validate_supplier_scale(kind, N, ratio * N)
        except ValueError as exc:
            supplier_fail += 1
            raise AssertionError("certified supplier scale sample was rejected") from exc

        uv = geometric_uv_natural_time_sum(N, c, 2.0)
        first = natural_duration(N, c)
        uv_gap = min(uv_gap, uv - first)
        if not (uv > first and math.isfinite(uv)):
            raise AssertionError("UV geometric natural-time obstruction disappeared")

    return QuotientStress(samples, wom, wtime, wscale, boundary_fail, supplier_fail, obstruction_fail, service_witness_fail, checkpoint_fail, segmentation_fail, maxowners, uv_gap)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-continuum-master-event-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "continuum_master_event_quotient.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Continuum master event quotient

Status: **{cert['status']}**.

The canonical recursive state contains physical event vertices only.  Raw HH/interface coefficient thresholds are interval locators until actual `Q^2` energy/work reentry; `full_natural_own_scale_service` is a same-corridor witness; and a complete no-hit natural horizon is an **analysis checkpoint**, not a recursive event.

Physical time remains exact. If checkpoint times are `t_0>=...>=t_L`, then `sum_j(t_j-t_(j+1))=t_0-t_L`. More strongly, a no-event horizon does not restart the event-anchored smooth carrier: its exact event/carrier/terminal-dual/PDE-path provenance, terminal coefficient and cumulative complex first-hit monitors continue. Genuine fixed-carrier natural windows have one positive `cA^-2` duration; arbitrary observer cuts have no service duration.

Endpoint hard-shell rereading at a full-natural checkpoint is likewise witness geometry.  The companion checkpoint theorem keeps the incoming hard shell `M`, the actual corridor scale `A=3M/4`, and endpoint hard-shell candidates `A,2A` distinct.  Their ratios `3/4,3/2` do not supply directional progress and, without a new physical event, those hard-shell readings do not replace the smooth carrier or define a causal scale lineage.

Thus, after zero-charge relays, observer gauges, coefficient locators, same-event donor circulation, same-corridor service layers, natural-horizon checkpoints and exact same-path segmentation are quotiented, an infinite recursive **event** path avoiding `t=0` must contain infinitely many genuine physical owner events. A geometric finite-time sum still warns about genuinely changing-scale physical producers, but it cannot be manufactured from duration-free checkpoint readings. Genuine UV progression enters only through certified tail work/dissipation or another physical event.

Stress: `{out.samples}` quotient/path states
- worst zero-charge owner-mass residual: `{out.worst_owner_mass_residual:.3e}`
- worst physical-time telescope residual: `{out.worst_time_telescope_residual:.3e}`
- worst log-scale telescope residual: `{out.worst_scale_telescope_residual:.3e}`
- bounded-scale boundary failures: `{out.bounded_scale_boundary_failures}`
- supplier-scale failures: `{out.supplier_scale_failures}`
- coefficient-obstruction barrier failures: `{out.coefficient_obstruction_barrier_failures}`
- full-natural service-witness barrier failures: `{out.service_witness_barrier_failures}`
- full-natural checkpoint barrier failures: `{out.checkpoint_barrier_failures}`
- same-carrier checkpoint-segmentation barrier failures: `{out.checkpoint_segmentation_barrier_failures}`
- largest relayed joint-owner set sampled: `{out.maximum_relay_owner_count}`
- minimum sampled changing-scale geometric time beyond its first window: `{out.minimum_uv_time_gap_to_naive_infinite_sum:.3e}`

This theorem does **not** prove global no-escape or Navier--Stokes regularity.  It removes no-event checkpoint segmentation from physical lineage; termination of infinitely recurring genuine physical owner events remains the open master problem.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
