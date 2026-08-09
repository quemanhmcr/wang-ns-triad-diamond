from fractions import Fraction
import math

import pytest

from src.first_hit_heat_reservoir_erosion import (
    CLEAN_HEAT_POOL_RATIO,
    PHYSICAL_HEAT_POOL_RATIO_UPPER,
    first_forced_non_oo_generation,
    first_hit_epoch_certificate,
    first_hit_kelvin_growth_upper,
    forced_non_oo_service_lower,
    heat_band_service_upper,
    old_pool_heat_capacity_upper,
    total_future_old_pool_heat_capacity_upper,
)


def test_clean_heat_pool_ratio_is_441_over_640_and_above_physical_bound():
    assert CLEAN_HEAT_POOL_RATIO == Fraction(441, 640)
    assert PHYSICAL_HEAT_POOL_RATIO_UPPER < float(CLEAN_HEAT_POOL_RATIO) < 0.7


def test_first_hit_boundary_stays_in_kelvin_growth_corridor():
    assert first_hit_kelvin_growth_upper(1 / 30) < 21 / 20
    with pytest.raises(ValueError):
        first_hit_kelvin_growth_upper(1 / 30 + 1e-4)


def test_heat_band_service_capacity_formula():
    out = heat_band_service_upper(
        child_frequency=100.0,
        band_frequency_upper=20.0,
        scaled_lifetime=2.0,
        energy_upper=3.0,
    )
    assert math.isclose(out, 24.0)


def test_old_pool_capacity_has_clean_heat_specific_geometric_ratio():
    c0 = old_pool_heat_capacity_upper(
        generation=0,
        initial_low_cut_ratio=0.2,
        initial_block_frequency=50.0,
        frame_energy_bound=1.5,
        global_energy=4.0,
    )
    c1 = old_pool_heat_capacity_upper(
        generation=1,
        initial_low_cut_ratio=0.2,
        initial_block_frequency=50.0,
        frame_energy_bound=1.5,
        global_energy=4.0,
    )
    assert math.isclose(c1 / c0, 441 / 640)
    total = total_future_old_pool_heat_capacity_upper(
        initial_low_cut_ratio=0.2,
        initial_block_frequency=50.0,
        frame_energy_bound=1.5,
        global_energy=4.0,
    )
    assert math.isclose(total / c0, 640 / 199)


def test_first_forced_non_oo_generation_is_minimal():
    S = 2.0
    C0 = 30.0
    q = first_forced_non_oo_generation(
        high_strain_heat_threshold=S,
        initial_old_capacity=C0,
        forced_non_oo_fraction=0.5,
    )
    r = 441 / 640
    assert C0 * r**q <= 1.0
    assert q == 0 or C0 * r ** (q - 1) > 1.0


def test_non_oo_is_exact_positive_remainder():
    assert math.isclose(forced_non_oo_service_lower(total_heat_service=5.0, old_old_service=1.75), 3.25)


def test_certificate_keeps_universal_renewal_open():
    cert = first_hit_epoch_certificate()
    assert "441/640" in cert["clean_ratio"]
    assert "does not prove universal slab renewal" in cert["scope"]
