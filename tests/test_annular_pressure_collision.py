import math

from src.annular_pressure_collision import (
    PressureShell,
    collision_certificate,
    critical_l2_mass_from_l3_charge,
    far_multipole_moment,
    no_fresh_far_coefficient,
    pressure_kernel_component_gradient_bound,
)


def test_pressure_kernel_clean_constant():
    assert pressure_kernel_component_gradient_bound() == 10.0


def test_fourth_power_beats_three_dimensional_packing():
    mu = 0.03
    cgeom = 2.0
    n0 = 3
    shells = []
    for n in range(n0, n0 + 6):
        count = int(cgeom * 2 ** (3 * n))
        shells.append(PressureShell(n, tuple([mu] * count)))
    moment = far_multipole_moment(shells, n0)
    infinite_bound = cgeom * mu * 2 ** (1 - n0)
    assert moment < infinite_bound


def test_collision_threshold_is_strict_contrapositive():
    rho = 0.25
    cnear = 4.0
    cfar = 9.0
    cert = collision_certificate(rho, cnear, cfar)
    mu = 0.999 * cert.fresh_mass_threshold
    assert cnear * mu ** 1.5 + cfar * mu < rho


def test_critical_l2_mass_from_l3_is_scale_invariant_formula():
    q = 0.125
    cb = 2.0
    got = critical_l2_mass_from_l3_charge(q, cb)
    assert math.isclose(got, q ** (2 / 3) / 4.0)


def test_no_fresh_far_coefficient_has_spare_power():
    c = no_fresh_far_coefficient(3.0, 0.5, first_shell=4, kernel_constant=400.0)
    assert math.isclose(c, 400.0 * 0.5 * 3.0 * 2 ** (1 - 4))
