from fractions import Fraction

from src.pressure_reservoir_sync import (
    amortized_pressure_pair_capacity_upper,
    exact_pressure_pair_certificate,
    pair_energy_service_ratio_upper,
    pressure_hessian_pair_energy_service_ratio_upper,
    pressure_hessian_pair_service_coefficient,
    total_pressure_hessian_pair_capacity_upper,
    pressure_pair_service_coefficient,
    total_pressure_pair_capacity_upper,
)



def test_objective_pressure_hessian_pair_one_fifth_life_exact():
    r = pressure_hessian_pair_energy_service_ratio_upper()
    assert r == Fraction(4084101, 20971520)
    assert r < Fraction(1, 5)


def test_objective_pressure_hessian_pair_amortized_capacity():
    base = pressure_hessian_pair_service_coefficient(1.0, 0.8, 4.0) * 3.0
    assert total_pressure_hessian_pair_capacity_upper(1.0, 0.8, 4.0, 3.0) < 1.25 * base

def test_pair_one_third_life_exact():
    assert pair_energy_service_ratio_upper() == Fraction(194481, 655360)
    assert pair_energy_service_ratio_upper() < Fraction(1, 3)


def test_pair_amortized_capacity():
    base = pressure_pair_service_coefficient(1.0, 0.8, 4.0) * 3.0
    for q in range(10):
        assert amortized_pressure_pair_capacity_upper(q, 1.0, 0.8, 4.0, 3.0) <= base * (1 / 3) ** q + 1e-14
    assert total_pressure_pair_capacity_upper(1.0, 0.8, 4.0, 3.0) < 1.5 * base


def test_certificate():
    c = exact_pressure_pair_certificate()
    assert c["hessian_pair_service_clean"] == "<1/5"
    assert c["pressure_third_pair_service_clean"] == "<1/3"
