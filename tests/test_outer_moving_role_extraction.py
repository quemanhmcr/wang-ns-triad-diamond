import math

import numpy as np

from src.outer_moving_role_extraction import (
    LOW_STRAIN_ACTION,
    ROLE_LOWER,
    TRANSPORT_CUT,
    bilinear_apply,
    outer_role_identity_residual,
    persistent_low_low_gap,
    transported_role_lower_radius,
)


def test_low_strain_affine_transport_keeps_role_above_low_low_output():
    lower = transported_role_lower_radius(ROLE_LOWER, LOW_STRAIN_ACTION)
    assert lower > 0.5
    assert persistent_low_low_gap(ROLE_LOWER, LOW_STRAIN_ACTION, TRANSPORT_CUT) > 0.0


def test_support_bound_has_expected_closed_form():
    assert math.isclose(
        transported_role_lower_radius(3.0 / 5.0, 1.0 / 30.0),
        (3.0 / 5.0) * math.exp(-1.0 / 30.0),
    )


def test_exact_outer_role_algebra_in_commuting_viscosity_model():
    rng = np.random.default_rng(19)
    n = 3
    T = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
    u = rng.normal(size=n) + 1j * rng.normal(size=n)
    V = rng.normal(size=n) + 1j * rng.normal(size=n)
    Q = np.diag([0.2, 0.7, 1.0]).astype(complex)
    dQ = np.diag([0.1, -0.2, 0.05]).astype(complex)
    D = np.diag([-1.0, -2.0, -3.0]).astype(complex)
    res = outer_role_identity_residual(
        tensor=T, u=u, V=V, Q=Q, dQ=dQ,
        viscosity_operator=D, viscosity=0.3,
    )
    assert np.linalg.norm(res) < 1e-11


def test_quadratic_low_high_algebra_is_exact_before_projection():
    rng = np.random.default_rng(23)
    n = 4
    T = rng.normal(size=(n, n, n))
    V = rng.normal(size=n)
    h = rng.normal(size=n)
    u = V + h
    lhs = -bilinear_apply(T, u, u) + bilinear_apply(T, V, u) + bilinear_apply(T, u, V)
    rhs = bilinear_apply(T, V, V) - bilinear_apply(T, h, h)
    assert np.linalg.norm(lhs - rhs) < 1e-11
