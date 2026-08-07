import math

from src.curvature_moat_balance import (
    affine_window_curvature_coeff,
    balanced_schedule_partial,
    critical_mass_lower_from_curvature,
    critical_mass_lower_from_unavoidable_error,
    localization_error,
    old_quadratic_schedule_partial,
    optimal_localization_error,
    optimal_moat_width,
    scalar_hessian_bernstein_constant,
)


def test_exact_moat_optimization():
    a, b, kappa = 2.0, 3.0, 0.04
    M = optimal_moat_width(a, b, kappa)
    assert math.isclose(localization_error(a, b, kappa, M), optimal_localization_error(a, b, kappa), rel_tol=1e-14)


def test_affine_window_curvature_coefficient():
    assert math.isclose(affine_window_curvature_coeff(2.0, 3.0, 4.0), 12.0)


def test_curvature_failure_forces_critical_mass():
    a, b, C, eta = 1.3, 0.7, 2.1, 0.05
    kappa = eta * eta / (4 * a * b)
    direct = critical_mass_lower_from_curvature(kappa, C)
    packaged = critical_mass_lower_from_unavoidable_error(eta, a, b, C)
    assert math.isclose(direct, packaged, rel_tol=1e-14)


def test_old_quadratic_schedule_has_curvature_harmonic_countermodel():
    c1, h1 = old_quadratic_schedule_partial(1000)
    c2, h2 = old_quadratic_schedule_partial(100000)
    assert c2 - c1 < 0.002
    assert h2 - h1 > 4.0


def test_balanced_schedule_makes_both_terms_summable():
    c, h, tail = balanced_schedule_partial(100000)
    assert math.isclose(c, h, rel_tol=1e-14)
    assert tail < 0.007
    # Total infinite sum is bounded by the partial sum plus the integral tail.
    assert c + tail < 1.3


def test_scalar_hessian_bernstein_constant_positive():
    C = scalar_hessian_bernstein_constant(1.0)
    assert 0.08 < C < 0.09
