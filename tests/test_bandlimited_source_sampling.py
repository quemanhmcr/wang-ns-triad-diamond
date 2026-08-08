import pytest

from src.bandlimited_source_sampling import (
    exact_scaling_certificate,
    germano_l32_power_upper,
    pressure_l32_power_upper,
    pressure_source_samples_from_mass_and_increments,
    resolvable_cluster_count_upper,
    sgs_sampling_scaling_exponent,
    sgs_source_sample_sum_from_increments,
    viscous_sampling_scaling_exponent,
    viscous_source_sample_square_sum_upper,
)


def test_exact_scale_cancellation():
    assert sgs_sampling_scaling_exponent() == 0
    assert viscous_sampling_scaling_exponent() == 0


def test_sgs_increment_route_positive():
    out = sgs_source_sample_sum_from_increments(0.2, 1.0, 1.1, 2.0, 3.0, 4.0)
    assert out > 0


def test_viscous_replication_is_additive_dissipation():
    out = viscous_source_sample_square_sum_upper(0.3, 1.0, 1.1, 2.0, 3.0, 4.0)
    assert out > 0


def test_pressure_routes_to_mass_and_increment():
    rpow = germano_l32_power_upper(1.0, 0.2)
    ppow = pressure_l32_power_upper(0.4, rpow, 2.0, 1.5)
    out = pressure_source_samples_from_mass_and_increments(0.4, 0.2, 1.0, 1.1, 2.0, 3.0, 4.0, 2.0, 1.5)
    assert ppow > 0 and out > 0


def test_certificate():
    c = exact_scaling_certificate()
    assert c["SGS_pressure_sampling_exponent"] == "0"


def test_resolvable_cluster_count():
    assert resolvable_cluster_count_upper(8.0, 2.0, 1.5) == pytest.approx(8.0 / (2.0**1.5))
