from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

from src.curl_spectral_curvature_balance import critical_extremizer_parent_ratio

STATUS = (
    "DRAFT_INTRINSIC_CRITICAL_PARABOLIC_DONOR_THRESHOLD__"
    "C_GT_ONE_OVER_8_SQRT2_FORCES_MEDIAN_DONOR_RATIO_LT_5_OVER_8__"
    "BACKWARD_NATURAL_LIFETIME_GT_64_OVER_25__NO_LOG_PROGRESS_J"
)

CRITICAL_THRESHOLD_SQUARED = Fraction(1, 128)
DONOR_RATIO_UPPER = Fraction(5, 8)
LIFETIME_RATIO_LOWER = Fraction(64, 25)
SLICE_SAFE_UPPER_SQUARED = Fraction(3025, 393216)
SLICE_GAP_SQUARED = CRITICAL_THRESHOLD_SQUARED - SLICE_SAFE_UPPER_SQUARED


def critical_threshold() -> float:
    return 1.0 / (8.0 * math.sqrt(2.0))


def barycentric_majority_critical_squared(D: float, t: float) -> float:
    """Exact C_A^2 in barycentric coordinates.

    D is the same-helicity median-donor/child radius ratio and
    t=(1+D-S)/(2D) lies in (1/2,1) for a strict-UV majority-child triad.
    """
    d, q = float(D), float(t)
    if not (0.0 < d < 1.0 and 0.5 < q < 1.0):
        raise ValueError("strict-UV barycentric domain requires 0<D<1 and 1/2<t<1")
    return 2.0*d*d*(1.0-d)**2*q**3*(1.0-q)*(1.0-d*q)*(1.0+d-d*q)


def donor_monotonicity_polynomial(D: float, u: float) -> float:
    """Positive factor controlling -partial_D C_A^2 for D>=5/8.

    Here u=1-t belongs to (0,1/2).
    """
    d, x = float(D), float(u)
    return (
        6*d**3*x*x - 6*d**3*x - 4*d*d*x*x + 14*d*d*x
        - 5*d*d - 6*d*x + 7*d - 2
    )


@dataclass(frozen=True)
class CriticalParabolicDonorThreshold:
    critical_efficiency_threshold: float
    donor_ratio_upper: float
    backward_natural_lifetime_ratio_lower: float
    sharp_extremizer_donor_ratio: float
    slice_safe_upper_squared: Fraction
    slice_gap_squared: Fraction
    uses_log_progress_J: bool = False
    creates_event_budget: bool = False

    def __post_init__(self) -> None:
        if self.critical_efficiency_threshold != critical_threshold():
            raise AssertionError("critical threshold changed")
        if self.donor_ratio_upper != float(DONOR_RATIO_UPPER):
            raise AssertionError("donor ratio threshold changed from 5/8")
        if self.backward_natural_lifetime_ratio_lower != float(LIFETIME_RATIO_LOWER):
            raise AssertionError("parabolic lifetime ratio changed from 64/25")
        if not self.sharp_extremizer_donor_ratio < self.donor_ratio_upper:
            raise AssertionError("sharp critical extremizer left D<5/8")
        if self.slice_safe_upper_squared != SLICE_SAFE_UPPER_SQUARED:
            raise AssertionError("D=5/8 safe upper changed")
        if self.slice_gap_squared != SLICE_GAP_SQUARED or self.slice_gap_squared <= 0:
            raise AssertionError("critical threshold lost its exact positive gap")
        if self.uses_log_progress_J or self.creates_event_budget:
            raise ValueError("intrinsic critical threshold may not import J or create a synthetic budget")


def theorem_certificate() -> CriticalParabolicDonorThreshold:
    return CriticalParabolicDonorThreshold(
        critical_efficiency_threshold=critical_threshold(),
        donor_ratio_upper=float(DONOR_RATIO_UPPER),
        backward_natural_lifetime_ratio_lower=float(LIFETIME_RATIO_LOWER),
        sharp_extremizer_donor_ratio=critical_extremizer_parent_ratio(),
        slice_safe_upper_squared=SLICE_SAFE_UPPER_SQUARED,
        slice_gap_squared=SLICE_GAP_SQUARED,
    )


def strong_critical_event_backward_lifetime_ratio(
    *,
    median_donor_ratio: float,
    normalized_critical_production: float,
) -> float:
    """Fail-closed corollary after the analytic threshold theorem is supplied.

    A strict-UV event with normalized critical production above 1/(8sqrt2)
    must be the majority-child heterochiral sector and have D<5/8.  Its median
    donor therefore lives for D^{-2}>64/25 child natural lifetimes backward.
    """
    d = float(median_donor_ratio)
    c = float(normalized_critical_production)
    if not math.isfinite(d) or not (0.0 < d < 1.0):
        raise ValueError("physical median donor ratio in (0,1) required")
    if not math.isfinite(c) or not c > critical_threshold():
        raise ValueError("strictly super-threshold critical production required")
    if not d < float(DONOR_RATIO_UPPER):
        raise ValueError("analytic critical threshold theorem requires D<5/8")
    ratio = 1.0/(d*d)
    if not ratio > float(LIFETIME_RATIO_LOWER):
        raise AssertionError("critical donor lost strict backward parabolic expansion")
    return ratio
