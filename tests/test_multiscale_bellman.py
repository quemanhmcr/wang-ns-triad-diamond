import math

import numpy as np

from src.multiscale_bellman import (
    component_score,
    equal_branch_cascade,
    holonomy_convex_cost,
    refine_masses,
    refinement_certificate,
    single_edge_optimum,
)


def test_deterministic_refinement_has_zero_cost():
    x = np.array([0.3, 0.7])
    y = np.array([0.4, 0.6])
    z = np.array([0.5, 0.5])
    a = np.array([[1.0, 0.0], [0.0, 1.0]])
    c = refinement_certificate(x, y, z, a, a, a)
    assert abs(c.exact_ratio - 1.0) < 1e-12
    assert abs(c.total_cost) < 1e-12


def test_equal_binary_split_costs_log_two():
    x = y = z = np.array([1.0])
    a = np.array([[0.5, 0.5]])
    c = refinement_certificate(x, y, z, a, a, a)
    assert abs(c.exact_ratio - 0.5) < 1e-12
    assert abs(c.entropy_cost - math.log(2.0)) < 1e-12


def test_reuse_factor_is_additive_log_cost():
    x = y = z = np.array([1.0])
    a = np.array([[1.0]])
    rho = np.array([0.9])
    c = refinement_certificate(x, y, z, a, a, a, rho)
    assert abs(c.exact_ratio - 0.9) < 1e-12
    assert abs(c.reuse_cost + math.log(0.9)) < 1e-12


def test_refinement_bound_random():
    rng = np.random.default_rng(4)
    for _ in range(500):
        m, k = 4, 5
        x, y, z = (rng.dirichlet(np.ones(m)) for _ in range(3))
        ax, ay, az = (rng.dirichlet(np.ones(k), size=m) for _ in range(3))
        rho = rng.uniform(0.75, 1.0, size=m)
        c = refinement_certificate(x, y, z, ax, ay, az, rho)
        assert c.exact_ratio <= c.combined_bound + 1e-11
        assert -math.log(max(c.exact_ratio, 1e-300)) + 1e-11 >= c.total_cost


def test_pure_branching_telescopes():
    c = equal_branch_cascade(depth=7, branches=2, reuse_factor=1.0)
    assert abs(c["product_exact"] - 2.0 ** -7) < 1e-12
    assert abs(c["final_component_score"] - 2.0 ** -7) < 1e-12
    assert abs(c["product_exact"] - c["predicted_product"]) < 1e-12


def test_holonomy_soft_threshold_formula():
    gamma, a, b = 0.5, 0.2, 4.0
    cost = holonomy_convex_cost(gamma, a, b)
    expected = a * gamma - a * a / (4.0 * b)
    assert abs(cost - expected) < 1e-14


def test_single_edge_optimum_equation():
    r, _, gamma = single_edge_optimum()
    assert abs(-math.log(r) - (4.0 * r * r - 1.0)) < 2e-6
    assert abs(gamma + math.log(r)) < 1e-14
