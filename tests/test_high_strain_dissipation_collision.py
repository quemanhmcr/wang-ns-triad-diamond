import math

from src.high_strain_dissipation_collision import (
    clean_high_strain_dissipation_lower,
    geometric_fresh_energy_sum,
    geometric_physical_cost_sum,
    gradient_linf_bernstein_constant,
    normalized_dissipation_lower,
)


def test_clean_bernstein_constant():
    assert math.isclose(gradient_linf_bernstein_constant(.25),1/(8*math.sqrt(6)*math.pi),rel_tol=1e-14)


def test_high_strain_clean_constant():
    assert math.isclose(clean_high_strain_dissipation_lower(1.0),32*math.pi**2/75,rel_tol=1e-14)


def test_collision_scaling_is_dimensionless():
    K=.07; c=.4
    assert math.isclose(normalized_dissipation_lower(K,c),384*math.pi**2*K*K/c,rel_tol=1e-14)


def test_geometric_critical_costs_are_summable():
    assert math.isclose(geometric_physical_cost_sum(2.,4.,2.,.5),.5)
    assert math.isclose(geometric_fresh_energy_sum(3.,6.,2.),1.)
