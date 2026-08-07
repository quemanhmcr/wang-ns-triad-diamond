import math

import numpy as np

from src.affine_grain_dynamics import (
    AVERAGE_DEFICIT,
    FROZEN_STRAIN_TIME,
    POINTWISE_HODGE,
    affine_window_material_error_bound,
    dual_center_derivative,
    effective_shape_driver,
    extremal_shape_rates,
    frozen_strain_average_deficit_lower,
    frozen_strain_shape_coords,
    fourier_gaussian_rhs,
    kelvin_amplitude_rhs,
    kelvin_energy_rate,
    logdet_precision_rate,
    planar_strain_coercivity_ratio,
    sym,
    tracefree_2x2,
)


def test_extremal_strain_kernel_is_planar_conformal():
    scalar = 1.7 * np.eye(2)
    u, v, h = extremal_shape_rates(scalar)
    assert abs(u) < 1e-14
    assert abs(v) < 1e-14
    assert abs(h) < 1e-14


def test_extremal_planar_strain_coercivity():
    D = np.array([[1.0, 0.3], [0.3, -1.0]])
    assert planar_strain_coercivity_ratio(D) > 0.43


def test_frozen_strain_stays_local_and_has_hodge_cost():
    psi = 0.731
    R = np.array([[math.cos(psi), -math.sin(psi)], [math.sin(psi), math.cos(psi)]])
    D = R @ np.diag([1.0, -1.0]) @ R.T
    t = float(FROZEN_STRAIN_TIME)
    u, v, H = frozen_strain_shape_coords(D, t)
    assert abs(u) <= 2 / 25 + 1e-12
    assert abs(v) <= 2 / 25 + 1e-12
    assert H >= float(POINTWISE_HODGE) * t * t - 1e-12
    assert frozen_strain_average_deficit_lower(1.0, t) == float(AVERAGE_DEFICIT) * t * t


def test_gaussian_dual_center_and_logdet_identities():
    A = np.array([[0.2, 0.4, -0.1], [0.0, -0.3, 0.2], [0.3, -0.1, 0.1]])
    assert abs(np.trace(A)) < 1e-14
    P = np.array([[2.0, 0.2, 0.0], [0.2, 1.5, 0.1], [0.0, 0.1, 1.1]])
    k = np.array([1.0, -2.0, 0.7])
    nu = 0.8
    assert np.allclose(dual_center_derivative(P, k, A, nu), A @ (P @ k), atol=2e-13)
    Pdot, _ = fourier_gaussian_rhs(P, k, A, nu)
    assert math.isclose(np.trace(np.linalg.solve(P, Pdot)), logdet_precision_rate(P, A, nu), rel_tol=1e-13, abs_tol=1e-13)


def test_isotropic_viscosity_is_first_order_shape_neutral():
    A = np.array([[0.4, -0.2, 0.1], [0.3, -0.1, 0.0], [0.0, 0.2, -0.3]])
    P = 2.3 * np.eye(3)
    Mnu = effective_shape_driver(P, A, 1.7)
    assert np.allclose(tracefree_2x2(Mnu[:2, :2]), tracefree_2x2(sym(A)[:2, :2]), atol=2e-13)


def test_kelvin_amplitude_preserves_transversality_and_energy_identity():
    A = np.array([[0.2, 0.4, -0.1], [0.0, -0.3, 0.2], [0.3, -0.1, 0.1]])
    k = np.array([1.0, 2.0, -1.0])
    a = np.array([1.0, 0.0, 1.0])
    a -= k * (k @ a) / (k @ k)
    nu = 0.5
    kdot = -A.T @ k
    adot = kelvin_amplitude_rhs(a, k, A, nu)
    assert abs(kdot @ a + k @ adot) < 1e-12
    assert math.isclose(2 * a @ adot, kelvin_energy_rate(a, k, A, nu), rel_tol=1e-12, abs_tol=1e-12)


def test_affine_window_error_bound_scales_linearly_in_radius():
    a = affine_window_material_error_bound(3.0, 0.1, 1.2, 1.1, 2.0)
    b = affine_window_material_error_bound(3.0, 0.2, 1.2, 1.1, 2.0)
    assert math.isclose(b, 2 * a)
