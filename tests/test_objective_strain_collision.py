import math

from src.objective_strain_collision import (
    derivative_bernstein_constant,
    far_pressure_hessian_no_fresh_coefficient,
    objective_source_channel_level,
    pressure_hessian_fresh_mass_threshold,
    pressure_hessian_kernel_component_bound,
    quadratic_objective_source_mass_lower,
    viscous_objective_source_mass_lower,
)


def test_pressure_hessian_kernel_constant():
    assert pressure_hessian_kernel_component_bound() == 150.0


def test_far_pressure_hessian_has_5_minus_3_gain():
    c3 = far_pressure_hessian_no_fresh_coefficient(1.0, first_shell=3)
    c4 = far_pressure_hessian_no_fresh_coefficient(1.0, first_shell=4)
    assert math.isclose(c4, c3 / 4.0)


def test_objective_source_pigeonhole_level():
    assert math.isclose(objective_source_channel_level(0.12, 2.0), 0.001)


def test_quadratic_and_viscous_mass_thresholds_invert_bounds():
    rho = 0.03
    CA = 0.4
    muq = quadratic_objective_source_mass_lower(rho, CA)
    assert math.isclose(4 * CA * CA * muq, rho)
    nu, C3 = 0.7, 0.3
    muv = viscous_objective_source_mass_lower(rho, nu, C3)
    assert math.isclose(nu * C3 * math.sqrt(muv), rho)


def test_pressure_fresh_threshold_excludes_both_half_budgets():
    rho, cn, cf = 0.2, 3.0, 5.0
    mu = pressure_hessian_fresh_mass_threshold(rho, cn, cf)
    assert cn * mu <= rho / 2 + 1e-14
    assert cf * mu <= rho / 2 + 1e-14


def test_derivative_bernstein_constants_scale_with_lambda():
    c1 = derivative_bernstein_constant(1, 1.0)
    c2 = derivative_bernstein_constant(1, 2.0)
    assert math.isclose(c2 / c1, 2 ** 2.5)
