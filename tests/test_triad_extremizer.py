import math

from src.triad_extremizer import (
    RSTAR_HI, RSTAR_LO, critical_equation, symmetric_gamma, symmetric_jstar, symmetric_rstar,
)


def test_certified_bracket_straddles_unique_critical_root():
    assert critical_equation(RSTAR_LO) > 0.0
    assert critical_equation(RSTAR_HI) < 0.0
    r = symmetric_rstar()
    assert RSTAR_LO < r < RSTAR_HI
    assert abs(critical_equation(r)) < 2e-15


def test_corrected_extremizer_constants():
    r = symmetric_rstar()
    assert abs(r - 0.6109041015867660) < 2e-15
    assert abs(symmetric_gamma(r) - 0.4928152853421352) < 2e-15
    assert abs(symmetric_jstar(r) - 0.1001101758561887) < 2e-15
    assert math.isclose(1.0 / r, 1.6369181306895697, rel_tol=0.0, abs_tol=3e-15)
