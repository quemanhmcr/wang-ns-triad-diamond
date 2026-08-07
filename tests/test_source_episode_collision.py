import math

from src.source_episode_collision import (
    h1_channel_normalized_integral_lower,
    h1_source_level_threshold,
    sgs_episode_thresholds,
    sgs_persistent_episode_costs,
    temporal_concentration_alternative,
    viscous_episode_thresholds,
    viscous_persistent_episode_costs,
)


def test_h1_scaled_source_threshold():
    I1, c = 0.2, 0.5
    sigma = h1_channel_normalized_integral_lower(I1, c)
    rho = h1_source_level_threshold(I1, c)
    assert abs(2 * c * rho - sigma) < 1e-14


def test_clean_entropy_constants():
    th = sgs_episode_thresholds(0.1, 0.5, 2.0, 1.0, 1.0, 1.0)
    assert abs(th["atomic_entropy_if_no_dominant"] - math.log(4)) < 1e-14
    assert abs(th["ancestry_entropy_or_pair_entropy"] - math.log(2)) < 1e-14
    assert abs(th["same_ancestry_pair_mass"] - 0.25) < 1e-14


def test_persistent_sgs_cost():
    th = sgs_episode_thresholds(0.1, 0.5, 2.0, 1.0, 1.0, 1.0)
    out = sgs_persistent_episode_costs(th, 0.2)
    assert abs(out["high_frequency_dissipation"] - 0.05 * th["high_enstrophy"]) < 1e-14


def test_persistent_viscous_cost():
    th = viscous_episode_thresholds(0.1, 0.5, 2.0, 1.0)
    out = viscous_persistent_episode_costs(th, 0.2)
    assert abs(out["resolved_dissipation"] - 0.1 * th["resolved_enstrophy"]) < 1e-14


def test_temporal_concentration():
    out = temporal_concentration_alternative(0.4, 0.5, 0.1)
    assert abs(out["superlevel_threshold"] - 0.4) < 1e-14
    assert abs(out["source_integral_on_superlevel"] - 0.2) < 1e-14
