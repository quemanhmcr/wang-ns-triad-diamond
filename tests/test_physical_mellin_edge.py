import math

from src.log_scale_flux import maximizing_orbit_rate_coefficients, mellin_flux_segments
from src.physical_mellin_edge import forward_mellin_coefficient
from src.triad_extremizer import symmetric_jstar, symmetric_rstar
from src.helical import coupling_magnitude_closed


def test_full_mellin_matches_flux_segment_identity():
    x, y = 0.57, 0.66
    for sq in (-1, 1):
        rates = maximizing_orbit_rate_coefficients(x, y, sq)
        seg = mellin_flux_segments(x, y, 1.0, *rates)
        g = coupling_magnitude_closed(x, y, 1.0, 1, -1, sq)
        # maximizing orbit has positive child coefficient x+y, so phase sign is +.
        assert math.isclose(forward_mellin_coefficient(x, y, 1, -1, sq), g * seg.total, rel_tol=1e-12, abs_tol=1e-12)


def test_equal_parent_full_mellin_reduces_to_old_J():
    r = symmetric_rstar()
    j = symmetric_jstar(r)
    for sq in (-1, 1):
        val = forward_mellin_coefficient(r, r, 1, -1, sq)
        assert abs(val - j) < 2e-12
