import math
from src.hermite_triad_selection import (
    A3, exact_arithmetic_certificate, parity_selection, sideband_rescue_upper,
    single_role_transfer_deficit_lower, small_odd_sideband_norm_growth,
)


def test_parity_selection():
    assert parity_selection(1) and parity_selection(3) and parity_selection(5)
    assert not parity_selection(0) and not parity_selection(2)


def test_rescue_starts_quadratic():
    rho=0.01
    assert math.isclose(sideband_rescue_upper(rho,rho,0.0),A3*rho*rho)
    assert sideband_rescue_upper(rho,0.0,0.0)==0.0


def test_small_odd_norm_growth_and_deficit():
    s=1/100
    assert small_odd_sideband_norm_growth(s)>1
    assert single_role_transfer_deficit_lower(s)==s*s/16


def test_exact_certificate():
    c=exact_arithmetic_certificate(); assert c['small_sideband_threshold']=='1/80'
