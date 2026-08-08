import math
from src.physical_flat_episode import (
    barycenter_step_error_upper, clean_potential_zeta_upper,
    half_cosine_error_clean, half_cosine_from_uv, tau_flat_kappa0_lower,
    tau_flat_zeta_upper,
)
from src.spherical_erosion import C_STAR


def test_exact_half_cosine_at_extremizer():
    r=1/(2*C_STAR)
    assert abs(half_cosine_from_uv(r,0,0)-C_STAR)<1e-13


def test_clean_zeta_tau_one_percent():
    z=tau_flat_zeta_upper(.01,1.0); assert z<.026
    assert tau_flat_kappa0_lower(.01,1.0)>.17


def test_barycenter_clean_formula():
    H=1e-4; d=.01
    e=barycenter_step_error_upper(H,d); assert e<.03
    assert clean_potential_zeta_upper(H,d)<=2*e+1e-15
