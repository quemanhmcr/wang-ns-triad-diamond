from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from src.objective_pressure_pair_atomization import (
    canonical_all_pair_absolute_capacity_upper,
    objective_pressure_pair_route,
)
from src.objective_source_routing_compiler import objective_owner_weight_threshold


STATUS = (
    "DRAFT_ALL_RESOLVED_PRESSURE_PAIR_DESCENDING_EPOCH_TELESCOPE__"
    "ACTUAL_POSITIVE_PAIR_SOURCE_VS_POINTWISE_ABSOLUTE_PAIR_CAPACITY__"
    "PARENT_FREQUENCY_FLOOR_PLUS_CERTIFIED_ONE_QUARTER_CHILD_SCALE__"
    "DIFFUSE_AND_DOMINANT_PAIR_LAWS_INCLUDED__NO_ENTROPY_COST_NO_CAPACITY_OWNER"
)

PRESSURE_PAIR_SCALE_RATIO_UPPER = 0.25
PAIR_CAPACITY_DENOMINATOR = 2560.0


@dataclass(frozen=True)
class ResolvedPressurePairOwnerCertificate:
    """Typed certificate for one actual resolved positive pressure-pair owner.

    The certificate is built by the already-certified objective pressure route.
    ``maximum_pair_frequencies`` is the actual unordered pair attaining q_max.
    At least one of those two hard u-shells carries the stored positive critical
    mass lower.  Pair entropy is retained only as diagnostic provenance.
    """

    parent_frequency: float
    pressure_source_weight: float
    pair_positive_source_weight: float
    maximum_pair_mass: float
    maximum_pair_frequencies: tuple[float, float]
    selected_pair_shell_mass_lower: float
    scaled_lifetime: float
    pair_entropy: float
    resolved_pair_owner: bool = True
    pair_entropy_used_as_cost: bool = False
    aggregate_resolved_energy_used_as_causal_owner: bool = False

    def __post_init__(self) -> None:
        positive = (
            self.parent_frequency,
            self.pressure_source_weight,
            self.pair_positive_source_weight,
            self.maximum_pair_mass,
            self.selected_pair_shell_mass_lower,
            self.scaled_lifetime,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in positive):
            raise ValueError("positive finite resolved pressure-pair owner data required")
        if not math.isfinite(self.pair_entropy) or self.pair_entropy < 0.0:
            raise ValueError("finite nonnegative pair entropy required")
        if not (0.0 < self.maximum_pair_mass <= 1.0):
            raise ValueError("maximum normalized pair mass must lie in (0,1]")
        fa, fb = self.maximum_pair_frequencies
        if not all(math.isfinite(x) and x > 0.0 for x in (fa, fb)):
            raise ValueError("positive finite maximum-pair frequencies required")
        tol = 8.0e-13 * max(1.0, self.parent_frequency, fa, fb)
        if max(fa, fb) > self.parent_frequency / 4.0 + tol:
            raise ValueError("resolved pressure pair escaped V=S_(N/4)u")
        if self.pair_positive_source_weight + 8.0e-13 * max(1.0, self.pressure_source_weight) < self.pressure_source_weight / 2.0:
            raise ValueError("resolved pair law is not a certified pressure half-owner")
        if not self.resolved_pair_owner:
            raise ValueError("typed certificate must represent an actual resolved pair owner")
        if self.pair_entropy_used_as_cost or self.aggregate_resolved_energy_used_as_causal_owner:
            raise ValueError("pressure entropy/capacity cannot be promoted to causal owner currency")


@dataclass(frozen=True)
class ResolvedPressurePairRenewalStep:
    owner: ResolvedPressurePairOwnerCertificate
    child_frequency: float
    child_critical_mass: float

    def __post_init__(self) -> None:
        if not isinstance(self.owner, ResolvedPressurePairOwnerCertificate):
            raise TypeError("typed resolved pressure-pair owner certificate required")
        M = float(self.child_frequency)
        mu = float(self.child_critical_mass)
        if not all(math.isfinite(x) and x > 0.0 for x in (M, mu)):
            raise ValueError("positive finite selected pressure child shell required")
        fa, fb = self.owner.maximum_pair_frequencies
        freq_tol = 8.0e-13 * max(1.0, M, fa, fb)
        if min(abs(M - fa), abs(M - fb)) > freq_tol:
            raise ValueError("selected child must be one of the actual maximal-pair hard shells")
        mass_tol = 8.0e-13 * max(1.0, mu, self.owner.selected_pair_shell_mass_lower)
        if mu + mass_tol < self.owner.selected_pair_shell_mass_lower:
            raise ValueError("selected maximal-pair child lost its certified critical-shell lower")
        if M > self.owner.parent_frequency / 4.0 + freq_tol:
            raise AssertionError("pressure-pair child lost the certified N_next<=N/4 support geometry")


@dataclass(frozen=True)
class ResolvedPressurePairEpochCertificate:
    event_count: int
    root_frequency: float
    last_parent_frequency: float
    final_child_frequency: float
    global_energy_upper: float
    pressure_source_weight_floor: float
    scaled_lifetime: float
    parent_frequency_floor: float
    maximum_event_count: int
    maximum_observed_scale_ratio: float
    minimum_observed_pair_mass: float
    maximum_observed_pair_entropy: float
    dominant_cut_used: bool = False
    pressure_entropy_used_as_cost: bool = False
    pair_capacity_used_as_causal_owner: bool = False
    pair_capacity_summed_across_events: bool = False
    generic_shell_registration_used_as_progress: bool = False

    def __post_init__(self) -> None:
        if self.event_count < 1:
            raise ValueError("nonempty resolved pressure-pair epoch required")
        positive = (
            self.root_frequency,
            self.last_parent_frequency,
            self.final_child_frequency,
            self.global_energy_upper,
            self.pressure_source_weight_floor,
            self.scaled_lifetime,
            self.parent_frequency_floor,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in positive):
            raise ValueError("positive finite pressure epoch certificate values required")
        if self.maximum_event_count < self.event_count:
            raise ValueError("parent-frequency floor does not cover observed pressure-pair epoch")
        if self.maximum_observed_scale_ratio > PRESSURE_PAIR_SCALE_RATIO_UPPER + 8.0e-13:
            raise ValueError("epoch contains non-pressure-pair scale progress")
        if not (0.0 < self.minimum_observed_pair_mass <= 1.0):
            raise ValueError("invalid observed maximum-pair mass diagnostic")
        if not math.isfinite(self.maximum_observed_pair_entropy) or self.maximum_observed_pair_entropy < 0.0:
            raise ValueError("invalid pair entropy diagnostic")
        if (
            self.dominant_cut_used
            or self.pressure_entropy_used_as_cost
            or self.pair_capacity_used_as_causal_owner
            or self.pair_capacity_summed_across_events
            or self.generic_shell_registration_used_as_progress
        ):
            raise ValueError("resolved pressure epoch used a forbidden shortcut")


def certify_resolved_pressure_pair_owner(
    pressure_source_weight: float,
    scaled_lifetime: float,
    parent_frequency: float,
    *,
    sgs_positive_source_weight: float,
    pair_positive_weights: Sequence[float],
    pair_shell_indices: Sequence[tuple[int, int]],
    pair_frequencies: Sequence[tuple[float, float]],
) -> ResolvedPressurePairOwnerCertificate:
    """Bind the epoch step to the certified actual positive pressure-pair route."""
    route = objective_pressure_pair_route(
        pressure_source_weight,
        scaled_lifetime,
        parent_frequency,
        sgs_positive_source_weight=sgs_positive_source_weight,
        pair_positive_weights=pair_positive_weights,
        pair_shell_indices=pair_shell_indices,
        pair_frequencies=pair_frequencies,
    )
    owners = tuple(route["joint_primary_owners"])
    if "resolved_pressure_pair_law" not in owners:
        raise ValueError("actual pressure source did not certify a resolved positive pair owner")
    imax = int(route["max_pair_witness_index"])
    freqs = tuple((float(a), float(b)) for a, b in pair_frequencies)
    if not (0 <= imax < len(freqs)):
        raise AssertionError("pressure route lost its maximal physical pair index")
    entropy = route["pair_source_entropy"]
    if not isinstance(entropy, dict):
        raise TypeError("pressure route lost its physical pair entropy diagnostic")
    return ResolvedPressurePairOwnerCertificate(
        parent_frequency=float(parent_frequency),
        pressure_source_weight=float(pressure_source_weight),
        pair_positive_source_weight=float(route["pair_positive_source_total"]),
        maximum_pair_mass=float(entropy["maximum_atom_mass"]),
        maximum_pair_frequencies=freqs[imax],
        selected_pair_shell_mass_lower=float(route["max_pair_u_shell_mass_lower"]),
        scaled_lifetime=float(scaled_lifetime),
        pair_entropy=float(entropy["H2_pair_source"]),
    )


def resolved_pressure_pair_source_floor_from_objective_stop(
    objective_variation_action_floor: float,
    scaled_lifetime: float,
) -> float:
    return objective_owner_weight_threshold(
        objective_variation_action_floor,
        scaled_lifetime,
    )


def resolved_pressure_pair_parent_frequency_floor(
    global_energy_upper: float,
    pressure_source_weight_floor: float,
    scaled_lifetime: float,
) -> float:
    """Local source-capacity inequality converted to a parent-frequency lower.

    At each scaled time the certified absolute pair envelope is

        sum cap_ab <= N ||V||_2^2 / 2560.

    Over a source interval of length at most c and with ||V||_2^2<=E_global,

        R_pair <= c N E_global / 2560.

    A resolved pressure-pair owner has R_pair>=Sigma_P/2>=sigma_*/2.  Hence

        N >= 1280 sigma_*/(c E_global).

    Capacity is used only as a pointwise upper envelope on an already-existing
    positive source law. It is not normalized, charged, or summed across events.
    """
    E = float(global_energy_upper)
    sigma = float(pressure_source_weight_floor)
    c = float(scaled_lifetime)
    if not all(math.isfinite(x) and x > 0.0 for x in (E, sigma, c)):
        raise ValueError("positive finite energy/source-floor/lifetime required")
    # Bind the numerical constant to the certified pointwise helper.
    unit_cap = canonical_all_pair_absolute_capacity_upper(E, 1.0)
    expected = E / PAIR_CAPACITY_DENOMINATOR
    if abs(unit_cap - expected) > 1.0e-14 * max(1.0, expected):
        raise AssertionError("canonical pressure pair absolute-capacity constant changed")
    return 0.5 * sigma / (c * unit_cap)


def _event_count_upper(root_frequency: float, parent_floor: float) -> int:
    N0 = float(root_frequency)
    Nmin = float(parent_floor)
    if not all(math.isfinite(x) and x > 0.0 for x in (N0, Nmin)):
        raise ValueError("positive finite root/floor frequencies required")
    if N0 + 8.0e-13 * max(1.0, N0, Nmin) < Nmin:
        return 0
    # Event j has parent N_j<=N_0*4^-j and every resolved-pair owner requires
    # N_j>=Nmin.  Thus j<=floor(log_4(N0/Nmin)); event count is j_max+1.
    raw = math.log(N0 / Nmin) / math.log(4.0)
    jmax = max(0, int(math.floor(raw + 8.0e-13)))
    return jmax + 1


def resolved_pressure_pair_epoch_telescope(
    steps: Sequence[ResolvedPressurePairRenewalStep],
    *,
    global_energy_upper: float,
    pressure_source_weight_floor: float,
    scaled_lifetime: float,
) -> ResolvedPressurePairEpochCertificate:
    """Close a consecutive epoch of arbitrary resolved pressure-pair owners.

    Dominant and diffuse pair laws are treated identically.  At every event the
    *actual positive* pair source supplies the lower R_pair>=Sigma_P/2, while the
    canonical pair capacity supplies only the local upper

        R_pair <= c N E_global / 2560.

    Hence every pressure-pair parent has one physical frequency floor.  The
    independent support theorem gives N_(j+1)<=N_j/4.  A geometrically descending
    sequence of parent frequencies cannot cross the fixed floor infinitely often.
    """
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty resolved pressure-pair epoch required")
    E = float(global_energy_upper)
    sigma_floor = float(pressure_source_weight_floor)
    c = float(scaled_lifetime)
    if not all(math.isfinite(x) and x > 0.0 for x in (E, sigma_floor, c)):
        raise ValueError("positive finite energy/source-floor/lifetime required")

    parent_floor = resolved_pressure_pair_parent_frequency_floor(E, sigma_floor, c)
    root = rows[0].owner.parent_frequency
    max_ratio = 0.0
    min_pmax = 1.0
    max_h2 = 0.0

    for i, row in enumerate(rows):
        owner = row.owner
        ctol = 8.0e-13 * max(1.0, c, owner.scaled_lifetime)
        if abs(owner.scaled_lifetime - c) > ctol:
            raise ValueError("one pressure epoch must use the same registered scaled lifetime c")
        if owner.pressure_source_weight + 8.0e-13 * max(1.0, sigma_floor) < sigma_floor:
            raise ValueError("resolved pressure owner fell below the uniform source floor")
        if i > 0:
            prev = rows[i - 1].child_frequency
            ftol = 8.0e-13 * max(1.0, prev, owner.parent_frequency)
            if abs(owner.parent_frequency - prev) > ftol:
                raise ValueError("consecutive pressure-pair events must recurse through the selected physical child shell")

        # Verify the local positive-source/capacity inequality using only the
        # supplied global kinetic-energy upper and the certified pair envelope.
        cap_upper = c * canonical_all_pair_absolute_capacity_upper(E, owner.parent_frequency)
        cap_tol = 8.0e-13 * max(1.0, cap_upper, owner.pair_positive_source_weight)
        if owner.pair_positive_source_weight > cap_upper + cap_tol:
            raise ValueError("resolved positive pair source exceeds the certified energy/capacity envelope")
        if owner.parent_frequency + 8.0e-13 * max(1.0, owner.parent_frequency, parent_floor) < parent_floor:
            raise ValueError("resolved pressure-pair parent fell below its source/energy-imposed frequency floor")

        max_ratio = max(max_ratio, row.child_frequency / owner.parent_frequency)
        min_pmax = min(min_pmax, owner.maximum_pair_mass)
        max_h2 = max(max_h2, owner.pair_entropy)

    count_upper = _event_count_upper(root, parent_floor)
    if len(rows) > count_upper:
        raise ValueError("observed resolved pressure-pair epoch exceeds the physical parent-floor/scale telescope")

    return ResolvedPressurePairEpochCertificate(
        event_count=len(rows),
        root_frequency=root,
        last_parent_frequency=rows[-1].owner.parent_frequency,
        final_child_frequency=rows[-1].child_frequency,
        global_energy_upper=E,
        pressure_source_weight_floor=sigma_floor,
        scaled_lifetime=c,
        parent_frequency_floor=parent_floor,
        maximum_event_count=count_upper,
        maximum_observed_scale_ratio=max_ratio,
        minimum_observed_pair_mass=min_pmax,
        maximum_observed_pair_entropy=max_h2,
    )


def theorem_certificate() -> dict[str, object]:
    sigma = resolved_pressure_pair_source_floor_from_objective_stop(1.0, 1.0)
    floor = resolved_pressure_pair_parent_frequency_floor(1.0, sigma, 1.0)
    if abs(sigma - 0.25) > 1.0e-14:
        raise AssertionError("objective pressure owner floor lost A/(4c)")
    if abs(floor - 320.0) > 1.0e-11:
        raise AssertionError("resolved pressure parent floor lost 320 A/(c^2 E_global)")
    return {
        "status": STATUS,
        "actual_positive_owner": "resolved pair source is actual positive pressure service with R_pair>=Sigma_P/2; capacity never defines the owner",
        "local_capacity": "pointwise sum_(a<=b) cap_ab<=N||V||_2^2/2560 and ||V||_2<=||u||_2 give R_pair<=c N E_global/2560 on one source episode",
        "parent_floor": "if Sigma_P>=sigma_*, every resolved-pair event has N>=1280 sigma_*/(c E_global); under sigma_*=A_*/(4c), N>=320 A_*/(c^2 E_global)",
        "scale_progress": "the same pressure-pair supplier exposes an actual hard child with N_next<=N/4; generic shell A=3M/4 is not used",
        "epoch_bound": "L consecutive resolved-pair events satisfy L<=1+floor(log_4(N_0/N_min)), N_min=1280 sigma_*/(c E_global)",
        "entropy": "no q_max cut and no H2 cost is used; dominant and diffuse resolved pressure-pair laws are both included",
        "capacity_semantics": "the absolute pair capacity is used only pointwise inside each source event to derive a parent-frequency impossibility; it is never normalized, charged, summed across events, or promoted to recurrence currency",
        "scope": "draft eventually-pure resolved pressure-pair epoch telescope; pressure-SGS/fresh-SGS service and arbitrary cross-family source/strain/HH alternation remain open; no global-regularity claim",
    }
