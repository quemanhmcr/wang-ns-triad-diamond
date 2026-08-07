import math
import numpy as np

from src.helical_phase_holonomy import (
    diamond_edge_data,
    diamond_holonomy_from_edges,
    diamond_phase_residuals,
    residual_holonomy,
    sharp_four_phase_cost,
    weighted_phase_cost_lower,
    wrap_angle,
)


def fixture():
    a = np.array([1.0, 0.2, 0.1])
    b = np.array([-0.1, 0.9, 0.3])
    c = np.array([0.2, -0.15, 1.1])
    signs = (1, -1, 1, 1, -1, 1)
    return a, b, c, signs


def test_edge_normal_gauge_reconstructs_global_phase():
    a, b, c, signs = fixture()
    for e in diamond_edge_data(a, b, c, signs).values():
        assert abs(wrap_angle(e["global_phase"] - e["reconstructed_phase"])) < 2e-12


def test_diamond_residual_holonomy_cancels_modal_phases():
    a, b, c, signs = fixture()
    edges = diamond_edge_data(a, b, c, signs)
    _, _, H = diamond_holonomy_from_edges(edges)
    phases = {"a": 0.2, "b": -1.1, "c": 0.7, "m": 2.2, "n": -0.9, "d": 1.4}
    r = diamond_phase_residuals(a, b, c, signs, phases)
    assert abs(wrap_angle(residual_holonomy(r) - H)) < 2e-12


def test_sharp_phase_cost_at_equal_residuals():
    H = 0.8
    ds = [H / 4, H / 4, -H / 4, -H / 4]
    actual = sum(1 - math.cos(x) for x in ds)
    assert math.isclose(actual, sharp_four_phase_cost(H), rel_tol=1e-12, abs_tol=1e-12)


def test_weighted_phase_cost_lower():
    H = 0.5
    floor = 0.03
    assert math.isclose(weighted_phase_cost_lower(H, floor), floor * sharp_four_phase_cost(H))
