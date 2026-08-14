from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Mapping, Sequence

from src.asynchronous_duhamel_sync import BACKWARD_FRACTION, common_reference_slice

STATUS = (
    "DRAFT_INTRINSIC_GENERATED_ENERGY_DONOR_TIME_TELESCOPE__"
    "ACTUAL_POSITIVE_CHILD_WORK_PLUS_DONOR_RATIO_LT_5_OVER_8__"
    "HEAVY_HALF_SUPPORT_FORCES_FINITE_BACKWARD_PHYSICAL_DEPTH__NO_J"
)

DONOR_CHILD_RATIO_UPPER = Fraction(5, 8)
LIFETIME_GROWTH_MIN = Fraction(64, 25)
HALF_CHILD_TO_DONOR_SPAN_MAX = Fraction(25, 128)
BACKWARD_FRACTION_LOCAL = Fraction(2, 5)
ONE_STEP_BACKSHIFT_MIN = Fraction(6859, 16000)
BACKWARD_SUM_COEFF = Fraction(6859, 24960)
ACTUAL_HH_GENERATION_BRANCH = "physical_high_high_transfer_generation"


@dataclass(frozen=True)
class GeneratedEnergyDonorStep:
    child_frequency: float
    donor_frequency: float
    scaled_lifetime: float
    work_support_start: float
    work_support_end: float
    physical_hh_work_mass: float
    physical_hh_work_total: float
    physical_hh_work_lower: float
    energy_reentry_branch: str = ACTUAL_HH_GENERATION_BRANCH
    donor_is_physical_energy_donor: bool = True
    uses_log_progress_J: bool = False
    coefficient_impulse_used_as_work: bool = False

    def __post_init__(self) -> None:
        positive = (
            self.child_frequency,
            self.donor_frequency,
            self.scaled_lifetime,
            self.physical_hh_work_mass,
            self.physical_hh_work_total,
            self.physical_hh_work_lower,
        )
        if not all(math.isfinite(v) and v > 0.0 for v in positive):
            raise ValueError("positive finite generated energy-donor data required")
        if not (math.isfinite(self.work_support_start) and math.isfinite(self.work_support_end)):
            raise ValueError("finite physical work support required")
        if self.work_support_start < 0.0 or self.work_support_end < self.work_support_start:
            raise ValueError("ordered nonnegative work support required")
        if self.energy_reentry_branch != ACTUAL_HH_GENERATION_BRANCH:
            raise TypeError("only actual physical HH generation may enter the donor telescope")
        if not self.donor_is_physical_energy_donor:
            raise TypeError("interaction-parent labels alone do not define energy ancestry")
        if self.uses_log_progress_J:
            raise ValueError("intrinsic donor telescope does not use the log-progress J classifier")
        if self.coefficient_impulse_used_as_work:
            raise TypeError("coefficient impulse cannot replace physical child-energy work")
        ratio = self.donor_frequency / self.child_frequency
        if not ratio < float(DONOR_CHILD_RATIO_UPPER):
            raise ValueError("generated physical energy donor must lie strictly below 5/8 child scale")
        tol = 8.0e-12 * max(1.0, self.physical_hh_work_mass, self.physical_hh_work_total, self.physical_hh_work_lower)
        if self.physical_hh_work_mass + tol < 0.5 * self.physical_hh_work_total:
            raise ValueError("selected work support does not carry a physical heavy half")
        if self.physical_hh_work_total + tol < self.physical_hh_work_lower:
            raise ValueError("actual HH work law does not realize the energy-gate lower")
        width = self.work_support_end - self.work_support_start
        time_tol = 8.0e-12 * max(1.0, self.child_natural_lifetime, self.donor_natural_lifetime)
        if width > 0.5 * self.child_natural_lifetime + time_tol:
            raise ValueError("selected work support is wider than one half child natural slab")
        if self.normalized_donor_span > float(HALF_CHILD_TO_DONOR_SPAN_MAX) + 8.0e-12:
            raise AssertionError("donor support exceeded the intrinsic 25/128 span")

    @property
    def child_natural_lifetime(self) -> float:
        return self.scaled_lifetime / (self.child_frequency * self.child_frequency)

    @property
    def donor_natural_lifetime(self) -> float:
        return self.scaled_lifetime / (self.donor_frequency * self.donor_frequency)

    @property
    def donor_child_ratio(self) -> float:
        return self.donor_frequency / self.child_frequency

    @property
    def normalized_donor_span(self) -> float:
        return (self.work_support_end - self.work_support_start) / self.donor_natural_lifetime

    @property
    def common_reference_time(self) -> float:
        return self.work_support_start - float(BACKWARD_FRACTION_LOCAL) * self.donor_natural_lifetime


def donor_step_from_energy_reentry(
    *,
    reentry: Mapping[str, object],
    selected_physical_half_slab: Mapping[str, object],
    child_frequency: float,
    donor_frequency: float,
    scaled_lifetime: float,
) -> GeneratedEnergyDonorStep:
    gate = reentry.get("energy_gate") if isinstance(reentry.get("energy_gate"), Mapping) else reentry
    if not isinstance(gate, Mapping) or str(gate.get("branch", "")) != ACTUAL_HH_GENERATION_BRANCH:
        raise TypeError("actual physical-energy HH generation gate required")
    try:
        lower = float(gate["physical_hh_work_lower"])
        start = float(selected_physical_half_slab["start"])
        end = float(selected_physical_half_slab["end"])
        mass = float(selected_physical_half_slab["mass"])
        total = float(selected_physical_half_slab["total"])
    except (KeyError, TypeError, ValueError) as exc:
        raise TypeError("incomplete physical energy/work half-slab certificate") from exc
    return GeneratedEnergyDonorStep(
        child_frequency=float(child_frequency),
        donor_frequency=float(donor_frequency),
        scaled_lifetime=float(scaled_lifetime),
        work_support_start=start,
        work_support_end=end,
        physical_hh_work_mass=mass,
        physical_hh_work_total=total,
        physical_hh_work_lower=lower,
    )


@dataclass(frozen=True)
class GeneratedEnergyDonorEpoch:
    layer_count: int
    cumulative_reference_backshift: float
    minimum_cumulative_backshift: float
    minimum_lifetime_growth: float
    maximum_donor_child_ratio: float
    first_common_reference_time: float
    last_common_reference_time: float
    hits_initial_boundary: bool
    uses_event_count_budget: bool = False
    uses_log_progress_J: bool = False

    def __post_init__(self) -> None:
        if self.layer_count < 1:
            raise ValueError("nonempty generated donor epoch required")
        if self.cumulative_reference_backshift < -1e-12 or self.minimum_cumulative_backshift < -1e-12:
            raise ValueError("backward displacement cannot be negative")
        if self.minimum_lifetime_growth + 1e-12 < float(LIFETIME_GROWTH_MIN):
            raise AssertionError("donor lineage lost the 64/25 lifetime growth")
        if self.maximum_donor_child_ratio >= float(DONOR_CHILD_RATIO_UPPER) + 1e-12:
            raise AssertionError("donor lineage escaped the strict 5/8 scale contraction")
        if self.uses_event_count_budget or self.uses_log_progress_J:
            raise ValueError("physical-time donor telescope may not use event counting or J")


def minimum_backward_displacement(initial_donor_lifetime: float, transitions: int) -> float:
    if initial_donor_lifetime <= 0.0 or transitions < 0:
        raise ValueError("positive initial lifetime and nonnegative transition count required")
    if transitions == 0:
        return 0.0
    g = float(LIFETIME_GROWTH_MIN)
    return float(BACKWARD_SUM_COEFF) * initial_donor_lifetime * (g**transitions - 1.0)


def generated_energy_donor_epoch_telescope(
    steps: Sequence[GeneratedEnergyDonorStep],
) -> GeneratedEnergyDonorEpoch:
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty generated energy-donor epoch required")
    c0 = rows[0].scaled_lifetime
    if any(abs(row.scaled_lifetime-c0) > 8e-13*max(1.0,abs(c0)) for row in rows[1:]):
        raise ValueError("one donor epoch must use one scaled natural-lifetime constant")

    common = tuple(row.common_reference_time for row in rows)
    lifetimes = tuple(row.donor_natural_lifetime for row in rows)
    growths: list[float] = []
    ratios = [row.donor_child_ratio for row in rows]

    for j, row in enumerate(rows):
        # Same common-slice theorem, now with the stronger intrinsic 25/128 span.
        s = common_reference_slice(row.work_support_start, row.work_support_end, row.donor_natural_lifetime)
        if abs(s-common[j]) > 8e-12*max(1.0,abs(s),abs(common[j])):
            raise AssertionError("donor common reference surface changed")
        if j+1 == len(rows):
            continue
        nxt = rows[j+1]
        freq_tol = 8e-12*max(1.0,row.donor_frequency,nxt.child_frequency)
        if abs(nxt.child_frequency-row.donor_frequency) > freq_tol:
            raise ValueError("next generated child is not the previous physical energy donor scale")
        time_tol = 8e-12*max(1.0,abs(common[j]),row.work_support_end,nxt.work_support_start,nxt.work_support_end,lifetimes[j])
        if nxt.work_support_start < common[j]-time_tol or nxt.work_support_end > row.work_support_end+time_tol:
            raise ValueError("next generated donor support is not inside the previous common physical interval")
        growth = lifetimes[j+1]/lifetimes[j]
        growths.append(growth)
        if growth < float(LIFETIME_GROWTH_MIN)-8e-12:
            raise AssertionError("generated energy-donor lifetime failed 64/25 growth")
        backshift = common[j]-common[j+1]
        required = float(ONE_STEP_BACKSHIFT_MIN)*lifetimes[j]
        if backshift + time_tol < required:
            raise AssertionError("generated energy-donor common surfaces did not move backward fast enough")

    transitions=len(rows)-1
    cumulative=common[0]-common[-1]
    minimum=minimum_backward_displacement(lifetimes[0],transitions)
    if cumulative + 1e-11*max(1.0,abs(cumulative),abs(minimum),lifetimes[0]) < minimum:
        raise AssertionError("energy-donor epoch lost geometric backward-time telescope")
    min_growth=min(growths) if growths else float(LIFETIME_GROWTH_MIN)
    return GeneratedEnergyDonorEpoch(
        layer_count=len(rows),
        cumulative_reference_backshift=cumulative,
        minimum_cumulative_backshift=minimum,
        minimum_lifetime_growth=min_growth,
        maximum_donor_child_ratio=max(ratios),
        first_common_reference_time=common[0],
        last_common_reference_time=common[-1],
        hits_initial_boundary=common[-1] <= 0.0,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "primitive_input": "actual positive HH child-energy work after physical-energy reentry, plus a proved physical energy-donor ratio N_d/N_c<5/8",
        "span": "one physical heavy half has width <=T_child/2 <(25/128)T_donor",
        "common_surface": "s=a-(2/5)T_donor is a registration surface, not a new event",
        "growth": "successive physical donor lifetimes grow by >64/25 backward",
        "backstep": "with alpha<=25/128, s_j-s_(j+1) >= (6859/16000)T_j",
        "telescope": "Delta s_L >= (6859/24960)T_0[(64/25)^L-1]",
        "scope": "no lower donor/child ratio, no log-progress J, no Duhamel probability, no event-count budget, and no global regularity claim",
    }
