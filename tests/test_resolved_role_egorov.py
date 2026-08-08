import math

import numpy as np

from src.resolved_role_egorov import (
    affine_subtracted_increment_upper,
    center_jet_shear_countermodel,
    coherent_deformation_variance_from_hermite,
    exact_affine_subtracted_commutator,
    isotropic_gaussian_packet_moments,
    scalar_cell_egorov_l2_upper,
)


def test_exact_affine_subtracted_commutator_identity():
    K = np.array([0.4, -0.2, 0.1, 0.3], complex)
    y = np.array([0.0, 0.2, -0.3, 0.5])
    V = np.array([0.1, 0.7, -0.4, 1.2])
    df = np.array([1 + 2j, -0.3j, 0.4 - 0.1j, -0.9 + 0.2j])
    lhs, rhs = exact_affine_subtracted_commutator(K, y, V, df, 0.35)
    assert np.linalg.norm(lhs - rhs) < 1e-13


def test_integral_taylor_curvature_bound_has_correct_geometry():
    H = 2.0
    assert math.isclose(affine_subtracted_increment_upper(H, 3.0, 4.0), 40.0)


def test_l2_egorov_bound_has_H_over_N_scaling():
    out = scalar_cell_egorov_l2_upper(
        hessian_sup=8.0,
        frequency_scale=4.0,
        kernel_l1=1.0,
        kernel_m1=2.0,
        kernel_m2=3.0,
        packet_x_moment=1.5,
        packet_grad_moment=0.5,
        packet_xgrad_moment=2.0,
    )
    assert math.isclose(out, 2.0 * (1.5 + 4.0 + 2.25))


def test_flat_center_jet_does_not_control_full_nonaffine_lowpass_residual():
    c = center_jet_shear_countermodel()
    assert c["second_mode_over_N"] == 0.25
    assert c["gradient_center"] == 0.0
    assert c["hessian_center"] == 0.0
    assert c["third_derivative_center"] > 0.0


def test_gaussian_poincare_says_coherent_deformation_variance_is_curvature_controlled():
    var, grad = coherent_deformation_variance_from_hermite([1, 2, 5], [3.0, 4.0, 2.0])
    assert math.isclose(var, 9.0)
    assert math.isclose(grad, 21.0)
    assert grad >= var


def test_isotropic_gaussian_moments_are_exact():
    m = isotropic_gaussian_packet_moments(1.0)
    assert math.isclose(m["Mx"], math.sqrt(3.0))
    assert math.isclose(m["Mg"], math.sqrt(7.0 / 4.0))
    assert math.isclose(m["Mxg"], math.sqrt(27.0 / 4.0))
