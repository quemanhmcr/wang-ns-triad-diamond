import math
import numpy as np
import pytest

from src.sgs_source_collision import (
    affine_factor_from_radius_upper,
    arb_constants_certificate,
    clean_sgs_source_upper,
    clean_viscous_source_upper,
    cubic_increment_from_sgs_source_lower,
    enstrophy_from_viscous_source_lower,
    fresh_radius_mass_lower,
    stress_l32_from_source_lower,
    stress_support_radius,
)


def test_support_radius():
    assert stress_support_radius(8.0) == 4.0


def test_affine_factor_radius_isotropic():
    N = 5.0
    L = np.eye(3) * (2.0 / N)
    actual, upper, kappa = affine_factor_from_radius_upper(L, N)
    assert abs(actual - 2.0) < 1e-12
    assert abs(upper - 2.0) < 1e-12
    assert abs(kappa - 1.0) < 1e-12


def test_sgs_collision_inverse():
    rho, s0 = 0.002, 1.7
    r = stress_l32_from_source_lower(rho, s0)
    assert clean_sgs_source_upper(s0, r) >= rho - 1e-15
    assert cubic_increment_from_sgs_source_lower(rho, s0, 1.0) > 0


def test_viscous_collision_inverse():
    rho, nu, s0 = 0.001, 1.0, 2.0
    d = enstrophy_from_viscous_source_lower(rho, nu, s0)
    assert clean_viscous_source_upper(nu, s0, d) >= rho - 1e-15


def test_radius_mass():
    assert abs(fresh_radius_mass_lower(3.0) - 0.9) < 1e-14


def test_arb_constants_optional():
    pytest.importorskip("flint")
    c = arb_constants_certificate()
    assert c["sgs_clean_constant"] == "1/800"
    assert c["viscous_clean_constant"] == "1/6000"
