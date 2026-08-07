import math

from src.log_scale_flux import (
    adverse_lower_to_upper_ratio,
    maximizing_orbit_mellin_coefficients,
    maximizing_orbit_rate_coefficients,
    mellin_flux_segments,
    sharp_cutoff_triad_flux,
)
from src.triad_extremizer import symmetric_rstar


def test_sharp_cutoff_piecewise_flux_and_mellin_identity():
    k, p, q = 0.6, 0.75, 1.0
    rates = (0.2, -0.9, 0.7)
    assert sharp_cutoff_triad_flux(0.5, k, p, q, *rates) == 0.0
    assert sharp_cutoff_triad_flux(0.7, k, p, q, *rates) == -0.2
    assert sharp_cutoff_triad_flux(0.9, k, p, q, *rates) == 0.7
    assert sharp_cutoff_triad_flux(1.1, k, p, q, *rates) == 0.0
    row = mellin_flux_segments(k, p, q, *rates)
    expected = -0.2 * math.log(p / k) + 0.7 * math.log(q / p)
    assert abs(row.total - expected) < 1e-15


def test_equal_parents_kill_lower_mellin_segment_exactly():
    r = symmetric_rstar()
    rates = maximizing_orbit_rate_coefficients(r, r, -1)
    row = mellin_flux_segments(r, r, 1.0, *rates)
    assert row.lower == 0.0
    assert abs(row.upper - (2.0 * r) * math.log(1.0 / r)) < 2e-15


def test_adverse_maximizing_orbit_has_backscatter_below_parent_split():
    x, y = 0.59, 0.63
    row = maximizing_orbit_mellin_coefficients(x, y, -1)
    assert row.lower < 0.0
    assert row.upper > 0.0
    assert row.total > 0.0
    assert abs(row.upper - (x + y) * math.log(1.0 / y)) < 2e-15
    assert abs(row.lower + (1.0 - y) * math.log(y / x)) < 2e-15


def test_favorable_child_sign_adds_lower_forward_segment():
    x, y = 0.59, 0.63
    row = maximizing_orbit_mellin_coefficients(x, y, 1)
    assert row.lower > 0.0
    assert row.upper > 0.0


def test_local_grid_is_well_inside_certified_ninety_percent_retention():
    r = symmetric_rstar()
    worst = 0.0
    for i in range(21):
        u = 0.08 * i / 20.0
        for j in range(21):
            v = -0.08 + 0.16 * j / 20.0
            R = r * math.exp(-v)
            x = R * math.exp(-u / 2.0)
            y = R * math.exp(u / 2.0)
            worst = max(worst, adverse_lower_to_upper_ratio(x, y))
    assert worst < 0.09
