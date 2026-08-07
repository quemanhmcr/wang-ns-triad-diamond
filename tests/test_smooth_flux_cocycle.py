import math
import numpy as np

from src.smooth_flux_cocycle import (
    ConservativeInteraction,
    ProgressEdge,
    bad_capacity_mass_bound,
    certify_midgap_block,
    compact_cdf,
    mellin_moment,
    near_extremal_gap_radius,
    child_transfer_density_condition_number,
    positive_core_mass_lower_bound,
    polarization_certificate,
    ramp_potential,
    smooth_tail_flux,
)


def test_even_kernel_cdf_symmetry():
    delta = 0.07
    for z in np.linspace(-0.1, 0.1, 101):
        assert abs(compact_cdf(float(z), delta) + compact_cdf(float(-z), delta) - 1.0) < 2e-15


def test_ramp_is_exactly_flat_and_linear_outside_transition():
    tau, delta = 0.4, 0.05
    assert ramp_potential(tau - 0.08, tau, delta) == 0.0
    assert abs(ramp_potential(tau + 0.09, tau, delta) - 0.18) < 1e-15


def test_filter_independent_mellin_moment_matches_sharp_triad_formula():
    logs = (math.log(0.57), math.log(0.66), 0.0)
    rates = (0.34, -1.12, 0.78)
    I = ConservativeInteraction(logs, rates)
    expected = -rates[0] * math.log(0.66 / 0.57) + rates[2] * math.log(1.0 / 0.66)
    assert abs(mellin_moment(I) - expected) < 1e-14


def test_single_edge_midgap_tail_recovers_upper_progress_exactly():
    p, q = math.log(0.64), 0.0
    tau = 0.5 * (p + q)
    delta = 0.03
    # k is below p; its rate and the parent redistribution do not enter once the smooth transition lies in the p-q gap.
    I = ConservativeInteraction((math.log(0.59), p, q), (0.2, -0.9, 0.7))
    assert p <= tau - delta and q >= tau + delta
    actual = smooth_tail_flux(I, tau, delta)
    expected = 0.7 * (q - p)
    assert abs(actual - expected) < 2e-14


def test_transfer_weighted_common_midgap_cancels_centering_error():
    edges = [
        ProgressEdge(-0.04, 0.46, 1.0),
        ProgressEdge(0.03, 0.52, 2.0),
        ProgressEdge(-0.01, 0.50, 4.0),
    ]
    cert = certify_midgap_block(edges, 0.04)
    assert cert.moat_margin > 0.0
    assert abs(cert.equality_residual) < 2e-15


def test_polarization_identity_and_bad_mass_bound():
    cap = [1.0, 2.0, 3.0]
    mult = [0.95, 0.8, 1.0]
    phase = [1.0, 0.5, -0.2]
    cert = polarization_certificate(cap, mult, phase)
    assert abs(cert.exact_residual) < 2e-15
    assert cert.total_deficit + 1e-15 >= cert.multiplier_deficit
    assert abs(bad_capacity_mass_bound(0.02, 0.1) - 0.2) < 1e-15


def test_near_extremal_good_core_forces_gap_and_positive_weight_comparability():
    eta = 1e-4
    a = near_extremal_gap_radius(eta)
    assert abs(a - (0.01 + 0.0025)) < 1e-15
    cond = child_transfer_density_condition_number(eta)
    assert 1.0 < cond < 1.06
    assert abs(positive_core_mass_lower_bound(2e-6, eta) - 0.98) < 1e-15
