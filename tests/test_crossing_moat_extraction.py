import math

from src.crossing_moat_extraction import (
    BINS,
    CROSSING_HODGE_COEFF,
    CrossingEdge,
    extract_four_bin_core,
    theoretical_moat_margin,
    theoretical_shell_halfwidth,
)
from src.triad_extremizer import symmetric_gamma, symmetric_rstar


def test_four_bin_crossing_pigeonhole_and_shell():
    gamma = symmetric_gamma(symmetric_rstar())
    tau = 1.3
    rows = []
    for j in range(20):
        gap = gamma + (-1) ** j * 0.01
        p = tau - gap * (j + 1) / 21
        rows.append(CrossingEdge(p - 0.004, p, p + gap, 1.0 + j, 0.1 + 0.03 * j))
    ex = extract_four_bin_core(rows, tau, gamma)
    assert ex.transfer_core.transfer >= ex.total_transfer / BINS
    assert ex.hodge_core.hodge_numerator >= ex.total_hodge_numerator / BINS
    assert max(ex.transfer_core.max_parent_deviation, ex.transfer_core.max_child_deviation) <= ex.theoretical_shell_halfwidth + 1e-12


def test_derived_constants_fit_old_common_moat():
    gamma = symmetric_gamma(symmetric_rstar())
    assert theoretical_shell_halfwidth(gamma) < 2 / 25
    assert theoretical_moat_margin(gamma) > 1 / 25
    assert math.isclose(float(CROSSING_HODGE_COEFF), 25 / 424)
