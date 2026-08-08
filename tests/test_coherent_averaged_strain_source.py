import math

import numpy as np

from src.coherent_averaged_strain_source import (
    coherent_average_corotational_strain_rhs,
    coherent_average_gradient_rhs,
    coherent_bar_quadratic_scaled_weight_upper,
    coherent_local_source_weight_upper,
    coherent_reynolds_scaled_source_weight_upper,
    coherent_reynolds_source_bounds,
    coherent_variance_scaled_weight_constant,
    inherited_filtered_source_routes,
    skew,
    sym,
    transport_integration_by_parts_bound,
)


def test_coherent_averaged_corotational_identity_is_exact_algebraically():
    rng = np.random.default_rng(17)
    A = rng.normal(size=(3, 3))
    A -= np.trace(A) / 3 * np.eye(3)
    F2 = rng.normal(size=(3, 3))
    HP = rng.normal(size=(3, 3))
    GR = rng.normal(size=(3, 3))
    DA = rng.normal(size=(3, 3))
    TR = rng.normal(size=(3, 3))
    nu = 0.37
    dA = coherent_average_gradient_rhs(A, F2, HP, GR, DA, TR, nu)
    direct = sym(dA) + sym(A) @ skew(A) - skew(A) @ sym(A)
    formula = coherent_average_corotational_strain_rhs(A, F2, HP, GR, DA, TR, nu)
    assert np.linalg.norm(direct - formula) < 1e-13


def test_reynolds_source_terms_are_quadratic_in_coherent_deformation():
    K = 0.2
    kap = 1.1
    out = coherent_reynolds_source_bounds(K, kap)
    assert math.isclose(out["fluctuation_square_upper"], kap * kap * K * K)
    assert math.isclose(out["residual_transport_upper"], math.sqrt(7.0) * kap * K * K)


def test_transport_integration_by_parts_cauchy_bound_matches_clean_factors():
    K = 0.3
    kap = 1.12
    b = transport_integration_by_parts_bound(7 * K * K, kap * kap * K * K)
    assert math.isclose(b, math.sqrt(7.0) * kap * K * K)


def test_scaled_averaging_sources_route_linearly_to_critical_DV():
    D = 2.5
    assert coherent_variance_scaled_weight_constant() > 0
    assert coherent_reynolds_scaled_source_weight_upper(D) > 0
    assert coherent_bar_quadratic_scaled_weight_upper(D) > 0
    assert math.isclose(
        coherent_local_source_weight_upper(D),
        coherent_reynolds_scaled_source_weight_upper(D) + coherent_bar_quadratic_scaled_weight_upper(D),
    )


def test_filtered_pressure_sgs_viscous_routes_are_inherited_without_new_average_factor():
    out = inherited_filtered_source_routes(0.01, 0.5)
    assert out["pressure"]["resolved_critical_mass"] == 28.5
    assert out["pressure"]["stress_l32"] == 1.9
    assert out["sgs_stress_l32"] == 3.8
    assert out["viscous_DV"] > 0
