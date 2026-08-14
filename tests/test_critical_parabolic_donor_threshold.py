from fractions import Fraction
import math

import pytest

from src.critical_parabolic_donor_threshold import (
    CRITICAL_THRESHOLD_SQUARED,
    DONOR_RATIO_UPPER,
    LIFETIME_RATIO_LOWER,
    SLICE_GAP_SQUARED,
    SLICE_SAFE_UPPER_SQUARED,
    barycentric_majority_critical_squared,
    critical_threshold,
    donor_monotonicity_polynomial,
    strong_critical_event_backward_lifetime_ratio,
    theorem_certificate,
)


def test_exact_threshold_has_positive_elementary_slice_gap():
    cert = theorem_certificate()
    assert CRITICAL_THRESHOLD_SQUARED == Fraction(1,128)
    assert SLICE_SAFE_UPPER_SQUARED == Fraction(3025,393216)
    assert SLICE_GAP_SQUARED == Fraction(47,393216)
    assert cert.critical_efficiency_threshold == pytest.approx(1/(8*math.sqrt(2)))
    assert not cert.uses_log_progress_J
    assert not cert.creates_event_budget


def test_barycentric_D_monotonicity_factor_is_positive_in_threshold_region():
    for D in (5/8, 0.64, 2/3, 0.75, 0.9, 0.99):
        for u in (1e-6, 0.1, 0.25, 0.499999):
            assert donor_monotonicity_polynomial(D, u) > 0.0


def test_D_five_eighth_slice_stays_below_intrinsic_threshold():
    D = 5/8
    maximum = 0.0
    for j in range(1, 5000):
        t = 0.5 + 0.5*j/5000
        if t >= 1.0:
            continue
        maximum = max(maximum, barycentric_majority_critical_squared(D, t))
    assert maximum < float(CRITICAL_THRESHOLD_SQUARED)
    assert maximum < critical_threshold()**2


def test_superthreshold_donor_has_strict_backward_parabolic_expansion():
    ratio = strong_critical_event_backward_lifetime_ratio(
        median_donor_ratio=0.60,
        normalized_critical_production=critical_threshold()*1.01,
    )
    assert ratio > float(LIFETIME_RATIO_LOWER)
    assert ratio > 64/25
    with pytest.raises(ValueError, match="D<5/8"):
        strong_critical_event_backward_lifetime_ratio(
            median_donor_ratio=float(DONOR_RATIO_UPPER),
            normalized_critical_production=critical_threshold()*1.01,
        )
