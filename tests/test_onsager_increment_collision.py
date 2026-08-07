import math

from src.onsager_increment_collision import (
    coarse_dyadic_bound_from_masses,
    exact_dyadic_square_weight,
    increment_collision_thresholds,
    packet_mass_entropy_route,
    persistent_dissipation_lower,
)


def test_dyadic_bound():
    masses = {-2: 0.2, -1: 0.3, 0: 0.4, 1: 0.1, 2: 0.05}
    assert exact_dyadic_square_weight(masses) <= coarse_dyadic_bound_from_masses(masses) + 1e-14


def test_collision_thresholds():
    th = increment_collision_thresholds(1.0, 1.0, 1.0, 1.0)
    assert abs(th["square_mass_threshold"] - 1.0) < 1e-14
    assert abs(th["low_band_critical_mass"] - 3.0 / 8.0) < 1e-14
    assert abs(th["high_normalized_enstrophy"] - 1.0 / 4.0) < 1e-14


def test_dominant_packet_route():
    out = packet_mass_entropy_route([0.7, 0.1, 0.1, 0.1], [0, 1, 2, 3], dominant_fraction=0.25)
    assert out["branch"] == "dominant_packet"


def test_entropy_or_cycle_route():
    masses = [1.0] * 8
    labels = [0, 0, 0, 0, 1, 1, 1, 1]
    out = packet_mass_entropy_route(masses, labels, dominant_fraction=0.25, ancestry_alpha=0.5)
    assert out["branch"] in {"ancestry_Bellman_entropy", "same_ancestry_pair_cycle"}


def test_persistent_dissipation():
    assert abs(persistent_dissipation_lower(0.3, 0.4) - 0.12) < 1e-14
