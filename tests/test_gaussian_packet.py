import math

import numpy as np

from src.gaussian_packet import (
    packet_center,
    scalar_frequency_ratio,
    scalar_position_ratio,
    scalar_width_ratio,
)


def test_packet_center_closes_and_matches_known_optimum():
    p, q, z, signs, j = packet_center()
    assert np.allclose(p + q, z)
    assert signs in ((1, -1, -1), (-1, 1, 1))
    assert abs(np.linalg.norm(p) - 0.6109040) < 2e-6
    assert abs(j - 0.100110175856) < 2e-8


def test_exact_scalar_ratios_are_one_only_at_alignment():
    assert scalar_frequency_ratio(np.zeros(3), 1.0) == 1.0
    assert scalar_frequency_ratio(np.array([1.0, 0.0, 0.0]), 1.0) == math.exp(-1.0 / 12.0)
    xs = np.zeros((3, 3))
    assert scalar_position_ratio(xs, 1.0) == 1.0
    xs[1, 0] = 1.0
    assert scalar_position_ratio(xs, 1.0) == math.exp(-2.0 / 3.0)
    assert abs(scalar_width_ratio(np.ones(3)) - 1.0) < 1e-15
    assert scalar_width_ratio(np.array([0.5, 1.0, 1.5])) < 1.0


def test_sharp_symmetric_young_constant():
    from src.gaussian_packet import sharp_young_constant
    assert abs(sharp_young_constant(3) - (math.sqrt(3.0) / 2.0) ** 3) < 1e-15
