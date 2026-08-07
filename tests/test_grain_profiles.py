import numpy as np

from src.grain_profiles import (
    component_score,
    dominant_component_certificate,
    holder_component_bound,
    two_branch_kernel,
)


def test_component_bellman_bound_random():
    rng = np.random.default_rng(123)
    for n in [2, 3, 5, 10]:
        for _ in range(300):
            x, y, z = [rng.dirichlet(np.ones(n)) for _ in range(3)]
            s = component_score(x, y, z)
            assert s <= holder_component_bound(x, y, z) + 1e-12
            assert s <= 1.0 + 1e-12


def test_two_equal_components_cost_half():
    x = np.array([0.5, 0.5])
    assert abs(component_score(x, x, x) - 0.5) < 1e-14


def test_dominant_certificate_is_valid():
    rng = np.random.default_rng(321)
    for _ in range(1000):
        x, y, z = [rng.dirichlet(np.ones(4) * 0.4) for _ in range(3)]
        cert = dominant_component_certificate(x, y, z)
        assert cert.x_mass + 1e-12 >= cert.x_lower
        assert cert.y_mass + 1e-12 >= cert.yz_lower_each
        assert cert.z_mass + 1e-12 >= cert.yz_lower_each


def test_gaussian_matched_branches_are_normalized():
    k = two_branch_kernel(0.02, 0.16, 2.0)
    assert abs(k[0,0,0] - 1.0) < 1e-9
    assert abs(k[1,1,1] - 1.0) < 1e-8


def test_cross_interaction_decays_with_spatial_separation():
    k0 = two_branch_kernel(0.02, 0.12, 0.0)
    k2 = two_branch_kernel(0.02, 0.12, 2.0)
    cross0 = max(k0[i,j,l] for i in range(2) for j in range(2) for l in range(2) if not (i==j==l))
    cross2 = max(k2[i,j,l] for i in range(2) for j in range(2) for l in range(2) if not (i==j==l))
    assert cross2 < cross0

def test_sharp_shared_component_certificate_near_one():
    x = np.array([0.995, 0.005])
    y = np.array([0.994, 0.006])
    z = np.array([0.996, 0.004])
    cert = dominant_component_certificate(x, y, z)
    assert cert.index == 0
    assert cert.y_mass + 1e-12 >= cert.yz_lower_each
    assert cert.z_mass + 1e-12 >= cert.yz_lower_each
