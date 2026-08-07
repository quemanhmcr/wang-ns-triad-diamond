import math

import numpy as np

from src.spherical_erosion import (
    C_STAR,
    KAPPA_STAR,
    TETRAHEDRAL_DIAMETER,
    THETA_STAR,
    angular_diameter,
    convex_hull_distance,
    equal_marginal_barycenter_test,
    hemisphere_barrier,
    lineage_barrier_certificate,
    make_companion,
    midpoint_barrier_margin,
    no_fresh_lifespan_bound,
    normalize,
    regular_tetrahedron,
    spherical_angle,
    spherical_midpoint,
)


def test_midpoint_geometry_and_support_identity():
    p = normalize(np.array([0.3, -0.2, 1.0]))
    q = make_companion(p, np.array([1.0, 0.5, 0.2]))
    m = spherical_midpoint(p, q)
    assert abs(spherical_angle(p, q) - THETA_STAR) < 1e-12
    assert abs(spherical_angle(p, m) - THETA_STAR / 2.0) < 1e-12
    pole = normalize(np.array([0.1, 0.2, 1.0]))
    lhs = float(pole @ m)
    rhs = float(pole @ (p + q)) / (2.0 * C_STAR)
    assert abs(lhs - rhs) < 1e-12


def test_exact_barrier_inequality():
    p = normalize(np.array([0.2, 0.1, 1.0]))
    q = make_companion(p, np.array([1.0, -0.3, 0.5]))
    pole = normalize(spherical_midpoint(p, q) + np.array([0.0, 0.0, 0.3]))
    assert min(pole @ p, pole @ q) > 0
    assert midpoint_barrier_margin(p, q, pole) >= -1e-12


def test_lineage_companion_barrier_telescope():
    pole = np.array([0.0, 0.0, 1.0])
    x0 = normalize(np.array([0.5, 0.0, 1.0]))
    states = [x0]
    companions = []
    for tangent in [
        np.array([1.0, 1.0, 2.0]),
        np.array([-1.0, 1.0, 2.0]),
        np.array([1.0, -1.0, 2.0]),
    ]:
        q = make_companion(states[-1], tangent)
        assert pole @ q > 0
        companions.append(q)
        states.append(spherical_midpoint(states[-1], q))
    cert = lineage_barrier_certificate(x0, companions, pole)
    assert cert["exact_margin"] >= -1e-11
    assert cert["simple_margin"] >= -1e-11
    assert cert["total_companion_barrier"] >= 2 * KAPPA_STAR * len(companions) - hemisphere_barrier(x0, pole) - 1e-11


def test_cap_lifespan_formula():
    assert no_fresh_lifespan_bound(math.radians(30)) == 0
    assert no_fresh_lifespan_bound(math.radians(60)) == 3
    assert no_fresh_lifespan_bound(math.radians(85)) > 10


def test_regular_tetrahedron_is_sharp_nonhemisphere_certificate():
    tetra = regular_tetrahedron()
    hull = convex_hull_distance(tetra)
    assert hull["distance"] < 1e-9
    assert abs(angular_diameter(tetra) - TETRAHEDRAL_DIAMETER) < 1e-12


def test_equal_marginal_barycenter_amplification():
    rng = np.random.default_rng(8)
    pairs = []
    for _ in range(20):
        p = normalize(rng.normal(size=3))
        q = make_companion(p, rng.normal(size=3))
        pairs.append((p, q))
    cert = equal_marginal_barycenter_test(pairs, rng.random(len(pairs)))
    assert cert["error"] < 1e-12
    assert abs(cert["amplification"] - 1.0 / C_STAR) < 1e-12


def test_barycenter_collision_entropy_bound():
    from src.spherical_erosion import barycenter_collision_certificate
    rng = np.random.default_rng(22)
    points = np.array([normalize(rng.normal(size=3)) for _ in range(30)])
    weights = rng.random(30)
    cert = barycenter_collision_certificate(points, weights)
    assert cert["margin"] >= -1e-12
    assert cert["entropy"] >= cert["entropy_lower_bound"] - 1e-12


def test_balanced_chain_entropy_is_linear_up_to_bounded_loss():
    from src.spherical_erosion import balanced_chain_entropy_bound
    a = balanced_chain_entropy_bound(12)
    b = balanced_chain_entropy_bound(48)
    assert a["total_entropy_lower_bound"] > 0
    assert b["bounded_deficit"] < 3.77
    assert b["total_entropy_lower_bound"] > 48 * math.log(2) - 3.77
