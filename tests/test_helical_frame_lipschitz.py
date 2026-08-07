import math
import numpy as np

from src.helical_frame_lipschitz import (
    helical_normal_gauge_derivative,
    normal_derivative,
    parent_angle_cos_from_uv,
    parent_helical_derivative_bound,
)
from src.helical_spin_transport import triad_normal


def test_good_core_angle_is_far_from_degenerate():
    r = 0.610904101586766
    for u in (0.0, 1 / 200):
        for v in (-1 / 100, 1 / 100):
            c = parent_angle_cos_from_uv(r, u, v)
            assert 1 / 4 < c < 2 / 5
            assert math.sqrt(1 - c * c) > 9 / 10


def test_normal_and_helical_differential_bounds():
    theta = 1.2
    a = np.array([math.cos(theta / 2), math.sin(theta / 2), 0.0])
    b = np.array([math.cos(theta / 2), -math.sin(theta / 2), 0.0])
    da = np.array([-a[1], a[0], 0.1]); da -= np.dot(da, a) * a
    db = np.array([b[1], -b[0], -0.1]); db -= np.dot(db, b) * b
    dn = normal_derivative(a, b, da, db)
    assert np.linalg.norm(dn) <= (10 / 9) * (np.linalg.norm(da) + np.linalg.norm(db))
    n = triad_normal(a, b)
    dh = helical_normal_gauge_derivative(a, n, da, dn, 1)
    assert np.linalg.norm(dh) <= parent_helical_derivative_bound(np.linalg.norm(da), np.linalg.norm(db))
