import math

import numpy as np

from src.dual_gaussian_root_registration import (
    affine_log_distance_upper_from_frobenius,
    covariance_cover_number_upper,
    dual_gaussian_pairing_lower,
    dual_probe_critical_mass_lower,
    dual_probe_l2_norm_sq,
    effective_root_frame_budget,
    frobenius_radius_for_affine_log_radius,
    normalized_covariance_eigenvalue_bounds,
    phase_space_color_count,
    representative_radius_ratio_lower,
    same_color_scale_shells_are_disjoint,
    scale_bin_shell_union,
)


def test_exact_dual_pairing_is_one_without_covariance_quantization():
    assert math.isclose(dual_gaussian_pairing_lower(0.0), 1.0)


def test_default_quantized_dual_probe_beats_clean_one_fifth_root_mass():
    eta = dual_probe_critical_mass_lower(0.01, 0.4)
    assert eta > 0.2


def test_dual_l2_norm_has_inverse_profile_radius_scaling():
    r = 2.5
    assert math.isclose(dual_probe_l2_norm_sq(r), 3 * math.sqrt(math.pi) / (2 * r))


def test_log_covariance_quantization_only_changes_radius_by_determinant_factor():
    d = 0.4
    assert math.isclose(representative_radius_ratio_lower(d), math.exp(-math.sqrt(3) * d / 6))


def test_covariance_bins_and_phase_space_colors_give_finite_scale_independent_budget():
    assert covariance_cover_number_upper() > 0
    assert phase_space_color_count(4, 6) == 5**6
    assert effective_root_frame_budget() > 0


def test_frobenius_covariance_net_implies_affine_log_radius_pointwise():
    m, _ = normalized_covariance_eigenvalue_bounds()
    d = 0.4
    eps = frobenius_radius_for_affine_log_radius(d, m)
    assert math.isclose(affine_log_distance_upper_from_frobenius(eps, m), d)


def test_four_scale_colors_make_outer_shell_projectors_disjoint():
    assert same_color_scale_shells_are_disjoint()
    lo0, hi0 = scale_bin_shell_union(0)
    lo4, hi4 = scale_bin_shell_union(4)
    assert hi0 < lo4
