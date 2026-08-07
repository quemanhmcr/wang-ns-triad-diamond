import math

from src.transfer_profile_extraction import (
    finite_energy_shell_lp_upper,
    gaussian_profile_transfer_lower,
    one_shot_profile_certificate,
    shell_critical_l2_mass_lower,
    trilinear_replacement_loss,
)


def test_trilinear_replacement_polynomial():
    eps = 0.01
    assert math.isclose(trilinear_replacement_loss(eps), eps + (1 + eps) * eps + (1 + eps) ** 2 * eps)


def test_profile_transfer_is_original_minus_remainder_budget():
    assert math.isclose(gaussian_profile_transfer_lower(0.999, 0.001), 0.999 - trilinear_replacement_loss(0.001))


def test_shell_lp_to_critical_l2_bridge_inverts():
    lp = 0.8
    cvol = 10.0
    mass = shell_critical_l2_mass_lower(lp, cvol)
    assert finite_energy_shell_lp_upper(mass, cvol) >= lp - 1e-12


def test_one_shot_profile_keeps_fixed_critical_mass():
    cert = one_shot_profile_certificate(0.999, 0.01, 8.0)
    assert cert.gaussian_transfer_lower > 0.96
    assert cert.gaussian_lp_mass_inside_shell_lower == 0.99
    assert cert.gaussian_critical_l2_mass_lower > 0.48
