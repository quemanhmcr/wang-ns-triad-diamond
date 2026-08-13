from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping, Sequence

from src.fresh_service_scale_reentry import (
    canonical_annular_frame_registration,
    fresh_service_scale_route,
)
from src.objective_source_routing_compiler import (
    objective_owner_weight_threshold,
    objective_sgs_integrated_square_service_lower,
)


STATUS = (
    "DRAFT_DESCENDING_FRESH_SGS_SCALE_EPOCH_TELESCOPE__"
    "SQUARE_FRAME_LOCAL_SERVICE_ENVELOPE_GIVES_PARENT_FREQUENCY_FLOOR__"
    "ONLY_SELECTED_PHYSICAL_CHILD_AT_MOST_ONE_HALF_RECURSION__"
    "TOP_BAND_J0_JMINUS1_RENEWAL_REMAINS_OPEN"
)

DESCENDING_FRESH_SCALE_RATIO_UPPER = 0.5
INCREMENT_TO_BAND_ENERGY = 4.0
FRESH_OWNER_FRACTION = 0.25


@dataclass(frozen=True)
class FreshSGSScaleOwnerCertificate:
    """One typed fresh NN scale owner produced by the certified SGS route."""

    parent_frequency: float
    forced_square_service_threshold: float
    fresh_service: float
    selected_band_index: int
    selected_band_service: float
    hard_shell_candidates: tuple[float, float]
    selected_hard_shell_mass_lower: float
    scaled_lifetime: float
    scale_entropy_inf: float
    scale_entropy2: float
    scale_entropy_used_as_cost: bool = False
    material_freshness_promoted_to_shell_owner: bool = False

    def __post_init__(self) -> None:
        positive = (
            self.parent_frequency,
            self.forced_square_service_threshold,
            self.fresh_service,
            self.selected_band_service,
            self.selected_hard_shell_mass_lower,
            self.scaled_lifetime,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in positive):
            raise ValueError("positive finite fresh-SGS owner data required")
        if self.selected_band_index > 0:
            raise ValueError("fresh low/base route requires j<=0")
        a, b = self.hard_shell_candidates
        if not all(math.isfinite(x) and x > 0.0 for x in (a, b)):
            raise ValueError("positive finite hard-shell candidates required")
        if max(a, b) > 2.0 * self.parent_frequency + 8.0e-13 * max(1.0, self.parent_frequency, a, b):
            raise ValueError("fresh scale route escaped the certified <=2N shell range")
        if self.fresh_service + 8.0e-13 * max(1.0, self.forced_square_service_threshold) < FRESH_OWNER_FRACTION * self.forced_square_service_threshold:
            raise ValueError("fresh service lost the certified Y/4 owner lower")
        if not all(math.isfinite(x) and x >= 0.0 for x in (self.scale_entropy_inf, self.scale_entropy2)):
            raise ValueError("finite nonnegative fresh scale diagnostics required")
        if self.scale_entropy_used_as_cost or self.material_freshness_promoted_to_shell_owner:
            raise ValueError("fresh scale diagnostics/material provenance cannot be promoted to causal cost or shell ownership")


@dataclass(frozen=True)
class DescendingFreshSGSRenewalStep:
    owner: FreshSGSScaleOwnerCertificate
    child_frequency: float
    child_critical_mass: float

    def __post_init__(self) -> None:
        if not isinstance(self.owner, FreshSGSScaleOwnerCertificate):
            raise TypeError("typed fresh SGS scale owner required")
        M = float(self.child_frequency)
        mu = float(self.child_critical_mass)
        if not all(math.isfinite(x) and x > 0.0 for x in (M, mu)):
            raise ValueError("positive finite selected fresh child shell required")
        a, b = self.owner.hard_shell_candidates
        tol = 8.0e-13 * max(1.0, M, a, b)
        if min(abs(M - a), abs(M - b)) > tol:
            raise ValueError("selected fresh child must be one of the certified two hard shells")
        if mu + 8.0e-13 * max(1.0, mu, self.owner.selected_hard_shell_mass_lower) < self.owner.selected_hard_shell_mass_lower:
            raise ValueError("selected fresh child lost its certified hard-shell lower")
        if M > DESCENDING_FRESH_SCALE_RATIO_UPPER * self.owner.parent_frequency + tol:
            raise ValueError("step is not in the strictly descending fresh-scale subepoch N_next<=N/2")


@dataclass(frozen=True)
class DescendingFreshSGSEpochCertificate:
    event_count: int
    root_frequency: float
    last_parent_frequency: float
    final_child_frequency: float
    global_energy_upper: float
    forced_square_service_floor: float
    scaled_lifetime: float
    parent_frequency_floor: float
    maximum_event_count: int
    maximum_observed_scale_ratio: float
    maximum_selected_band_index: int
    scale_entropy_used_as_cost: bool = False
    fresh_service_summed_as_global_reset: bool = False
    generic_shell_registration_used_as_progress: bool = False
    fixed_event_gap_used: bool = False

    def __post_init__(self) -> None:
        if self.event_count < 1:
            raise ValueError("nonempty descending fresh SGS epoch required")
        positive = (
            self.root_frequency,
            self.last_parent_frequency,
            self.final_child_frequency,
            self.global_energy_upper,
            self.forced_square_service_floor,
            self.scaled_lifetime,
            self.parent_frequency_floor,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in positive):
            raise ValueError("positive finite fresh epoch certificate values required")
        if self.maximum_event_count < self.event_count:
            raise ValueError("fresh parent-floor/half-scale telescope does not cover observed epoch")
        if self.maximum_observed_scale_ratio > DESCENDING_FRESH_SCALE_RATIO_UPPER + 8.0e-13:
            raise ValueError("non-descending fresh event entered descending epoch")
        if self.maximum_selected_band_index > 0:
            raise ValueError("invalid selected fresh band index")
        if (
            self.scale_entropy_used_as_cost
            or self.fresh_service_summed_as_global_reset
            or self.generic_shell_registration_used_as_progress
            or self.fixed_event_gap_used
        ):
            raise ValueError("descending fresh epoch used a forbidden shortcut")


def certify_fresh_sgs_scale_owner(
    forced_square_service_threshold: float,
    scaled_lifetime: float,
    parent_frequency: float,
    fresh_band_weights: Mapping[int, float],
    *,
    viscosity: float = 1.0,
) -> FreshSGSScaleOwnerCertificate:
    route = fresh_service_scale_route(
        forced_square_service_threshold,
        scaled_lifetime,
        parent_frequency,
        fresh_band_weights,
        viscosity=viscosity,
    )
    candidates = tuple(float(x) for x in route["hard_shell_candidates"])
    if len(candidates) != 2:
        raise AssertionError("fresh scale theorem lost its two-shell physical cover")
    return FreshSGSScaleOwnerCertificate(
        parent_frequency=float(parent_frequency),
        forced_square_service_threshold=float(forced_square_service_threshold),
        fresh_service=float(route["fresh_service"]),
        selected_band_index=int(route["selected_band_index"]),
        selected_band_service=float(route["selected_band_service"]),
        hard_shell_candidates=(candidates[0], candidates[1]),
        selected_hard_shell_mass_lower=float(route["hard_shell_mass_lower"]),
        scaled_lifetime=float(scaled_lifetime),
        scale_entropy_inf=float(route["H_inf_scale"]),
        scale_entropy2=float(route["H2_scale"]),
    )


def fresh_sgs_forced_service_floor_from_objective_stop(
    objective_variation_action_floor: float,
    scaled_lifetime: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> float:
    sigma = objective_owner_weight_threshold(
        objective_variation_action_floor,
        scaled_lifetime,
    )
    return objective_sgs_integrated_square_service_lower(
        sigma,
        filter_l1,
        lp_constant,
        bernstein_constant,
    )


def fresh_sgs_parent_frequency_floor(
    global_energy_upper: float,
    forced_square_service_floor: float,
    scaled_lifetime: float,
) -> float:
    """Fresh F>=Y/4 and the square-frame envelope imply N>=Y/(16cE).

    For the canonical square-normalized annular frame,

      F = sum_(j<=0) F_j
        <= 4 int sum_(j<=0) M_j ||u_j||_2^2 d tau
        <= 4 N int sum_j ||u_j||_2^2 d tau
        = 4 N int ||u||_2^2 d tau
        <= 4 c N E_global.

    The fresh owner has F>=Y/4.  No service is debited from a global account.
    """
    E = float(global_energy_upper)
    Y = float(forced_square_service_floor)
    c = float(scaled_lifetime)
    if not all(math.isfinite(x) and x > 0.0 for x in (E, Y, c)):
        raise ValueError("positive finite energy/service-floor/lifetime required")
    frame = canonical_annular_frame_registration()
    if "sum_j phi_j(xi)^2=1" not in str(frame["square_partition"]):
        raise AssertionError("canonical fresh annular frame lost square normalization")
    return Y / (16.0 * c * E)


def _event_count_upper(root_frequency: float, parent_floor: float) -> int:
    N0 = float(root_frequency)
    Nmin = float(parent_floor)
    if not all(math.isfinite(x) and x > 0.0 for x in (N0, Nmin)):
        raise ValueError("positive finite root/floor frequencies required")
    if N0 + 8.0e-13 * max(1.0, N0, Nmin) < Nmin:
        return 0
    raw = math.log(N0 / Nmin) / math.log(2.0)
    jmax = max(0, int(math.floor(raw + 8.0e-13)))
    return jmax + 1


def descending_fresh_sgs_epoch_telescope(
    steps: Sequence[DescendingFreshSGSRenewalStep],
    *,
    global_energy_upper: float,
    forced_square_service_floor: float,
    scaled_lifetime: float,
) -> DescendingFreshSGSEpochCertificate:
    """Close a consecutive fresh SGS epoch only while N_next<=N/2.

    The theorem intentionally stops at the first top-band/comparable-scale fresh
    renewal.  It does not use j=0 or j=-1 material vocabulary as a fake cost.
    """
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty descending fresh SGS epoch required")
    E = float(global_energy_upper)
    Yfloor = float(forced_square_service_floor)
    c = float(scaled_lifetime)
    if not all(math.isfinite(x) and x > 0.0 for x in (E, Yfloor, c)):
        raise ValueError("positive finite energy/service-floor/lifetime required")

    parent_floor = fresh_sgs_parent_frequency_floor(E, Yfloor, c)
    root = rows[0].owner.parent_frequency
    max_ratio = 0.0
    max_j = -10**9

    for i, row in enumerate(rows):
        owner = row.owner
        ctol = 8.0e-13 * max(1.0, c, owner.scaled_lifetime)
        if abs(owner.scaled_lifetime - c) > ctol:
            raise ValueError("one fresh epoch must use the same scaled lifetime c")
        if owner.forced_square_service_threshold + 8.0e-13 * max(1.0, Yfloor) < Yfloor:
            raise ValueError("fresh SGS owner fell below the uniform forced-service floor")
        if i > 0:
            prev = rows[i - 1].child_frequency
            ftol = 8.0e-13 * max(1.0, prev, owner.parent_frequency)
            if abs(owner.parent_frequency - prev) > ftol:
                raise ValueError("consecutive fresh SGS events must recurse through the selected physical child shell")

        # Local square-frame envelope only; no cross-event service accounting.
        local_upper = 4.0 * c * owner.parent_frequency * E
        stol = 8.0e-13 * max(1.0, local_upper, owner.fresh_service)
        if owner.fresh_service > local_upper + stol:
            raise ValueError("fresh service exceeds the canonical square-frame energy envelope")
        if owner.parent_frequency + 8.0e-13 * max(1.0, owner.parent_frequency, parent_floor) < parent_floor:
            raise ValueError("fresh SGS parent fell below its local service/energy-imposed frequency floor")

        max_ratio = max(max_ratio, row.child_frequency / owner.parent_frequency)
        max_j = max(max_j, owner.selected_band_index)

    count_upper = _event_count_upper(root, parent_floor)
    if len(rows) > count_upper:
        raise ValueError("observed descending fresh SGS epoch exceeds the physical parent-floor/half-scale telescope")

    return DescendingFreshSGSEpochCertificate(
        event_count=len(rows),
        root_frequency=root,
        last_parent_frequency=rows[-1].owner.parent_frequency,
        final_child_frequency=rows[-1].child_frequency,
        global_energy_upper=E,
        forced_square_service_floor=Yfloor,
        scaled_lifetime=c,
        parent_frequency_floor=parent_floor,
        maximum_event_count=count_upper,
        maximum_observed_scale_ratio=max_ratio,
        maximum_selected_band_index=max_j,
    )


def theorem_certificate() -> dict[str, object]:
    floor = fresh_sgs_parent_frequency_floor(1.0, 1.0, 1.0)
    if abs(floor - 1.0 / 16.0) > 1.0e-14:
        raise AssertionError("fresh parent-frequency floor lost Y/(16cE)")
    return {
        "status": STATUS,
        "local_envelope": "canonical square frame gives F_fresh<=4c N E_global, while the actual fresh owner has F_fresh>=Y/4; hence N>=Y/(16cE_global)",
        "objective_floor": "if SGS is a qualifying objective owner, Y>=C_Y A_*/(4c), so N>=C_Y A_*/(64 c^2 E_global)",
        "descent": "only physical selected hard-shell renewals with N_next<=N/2 enter this telescope; no generic A=3M/4 progress is used",
        "epoch_bound": "L consecutive strictly descending fresh-SGS events satisfy L<=1+floor(log_2(N_0/N_min)), N_min=Y_*/(16cE_global)",
        "top_band_survivor": "non-descending/comparable fresh renewals necessarily come from the top canonical LP geometry (j=0 or j=-1, depending on which of the two hard shells is selected); they are deliberately excluded and remain the hard source-service frontier",
        "semantics": "fresh scale H_inf/H2 remain deterministic diagnostics; NN provenance is not promoted to whole-shell freshness; service is never summed as an additive reset; no fixed event-gap clock is assumed",
        "scope": "draft strictly-descending fresh-SGS epoch telescope only; top-band fresh renewal, cross-family source/strain/HH alternation, and global regularity remain open",
    }
