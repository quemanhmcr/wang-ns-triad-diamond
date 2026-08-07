import math

from src.source_episode_collision import (
    h1_channel_normalized_integral_lower,
    sgs_source_linear_collision_coefficients,
    source_weight_partition_lower,
    source_weighted_sgs_episode_costs,
    source_weighted_viscous_episode_costs,
)


def test_h1_scaled_source_weight():
    assert abs(h1_channel_normalized_integral_lower(0.2, 0.5) - 0.2 / 66.0) < 1e-14
    assert abs(h1_channel_normalized_integral_lower(0.2, 0.5, 1800) - 0.2 / 900.0) < 1e-14


def test_sgs_linear_coefficients_positive():
    c = sgs_source_linear_collision_coefficients(2.0, 1.0, 1.0, 1.0)
    assert c["low_band_mass_per_source"] > 0
    assert c["high_enstrophy_per_source"] > 0


def test_clean_entropy_source_routing():
    out = source_weighted_sgs_episode_costs(0.1, 0.5, 2.0, 1.0, 1.0, 1.0)
    assert abs(out["atomic_entropy"] - math.log(4)) < 1e-14
    assert abs(out["ancestry_entropy"] - math.log(2)) < 1e-14
    assert abs(out["same_ancestry_pair_mass"] - 0.25) < 1e-14
    assert out["high_frequency_dissipation"] > 0


def test_viscous_no_persistence_cost():
    out = source_weighted_viscous_episode_costs(0.1, 0.5, 2.0, 1.0)
    assert out["resolved_dissipation"] > 0


def test_binary_source_weight_partition():
    assert abs(source_weight_partition_lower(1.0, 4) - 1.0 / 16.0) < 1e-14


def test_extended_source_divisor_routes_positive_costs():
    out = source_weighted_sgs_episode_costs(0.1, 0.5, 2.0, 1.0, 1.0, 1.0, source_divisor=1800)
    assert out["total_source_weight"] > 0
    assert out["high_frequency_dissipation"] > 0
    v = source_weighted_viscous_episode_costs(0.1, 0.5, 2.0, 1.0, source_divisor=1800)
    assert v["resolved_dissipation"] > 0
