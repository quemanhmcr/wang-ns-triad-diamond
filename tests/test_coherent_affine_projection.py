import math

from src.coherent_affine_projection import (
    EXTENDED_ASPECT,
    affine_residual_spectral_gap,
    coherent_deformation_to_dissipation_constant,
    gaussian_core_nonaffine_forcing_upper,
    gaussian_position_weight_upper,
    intrinsic_carrier_upper,
    longest_axis_from_radius_aspect,
    normalized_dissipation_from_coherent_deformation,
)


def test_ou_gap_after_removing_constant_and_affine_modes():
    out = affine_residual_spectral_gap([2, 3, 7], [4.0, 1.0, 2.0])
    assert math.isclose(out["residual_velocity_l2_sq"], 7.0)
    assert math.isclose(out["deformation_variance"], 25.0)
    assert math.isclose(out["spectral_gap_margin"], 11.0)


def test_creation_annihilation_plus_spectral_gap_gives_clean_seven():
    k2 = 10.0
    r2 = k2 / 2.0
    assert math.isclose(gaussian_position_weight_upper(r2, k2), 7.0 * k2)


def test_full_nonaffine_gaussian_core_forcing_constant():
    K = 0.2
    q = 1.5
    expected = (1.0 + q / math.sqrt(2.0) + math.sqrt(7.0) / 2.0) * K
    assert math.isclose(gaussian_core_nonaffine_forcing_upper(K, q), expected)


def test_radius_aspect_bounds_long_axis_and_intrinsic_carrier():
    assert math.isclose(longest_axis_from_radius_aspect(2.0, 8.0), 8.0)
    assert math.isclose(intrinsic_carrier_upper(2.0, 8.0, 1.25), 10.0)


def test_coherent_deformation_action_forces_critical_dissipation():
    C = coherent_deformation_to_dissipation_constant()
    assert C > 0
    I = 0.1
    c = 0.5
    lower = normalized_dissipation_from_coherent_deformation(I, c)
    assert math.isclose(lower, I * I / (C * c))
    assert coherent_deformation_to_dissipation_constant(EXTENDED_ASPECT) == C
