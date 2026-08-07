import numpy as np

from src.relative_polarization_transport import (
    J2,
    child_factor,
    common_hyperbolic_countermodel,
    parent_wedge,
    pointwise_capacity_bound,
    polarization_numerator,
    polarization_rhs,
    sl2_symmetric_step,
    tracefree_symmetric,
    wedge_rhs,
)


def test_tracefree_generators_are_symplectic_lie_algebra():
    for delta, beta in [(1.2, -0.3), (-0.7, 0.9), (0.0, 2.0)]:
        D = tracefree_symmetric(delta, beta)
        assert np.linalg.norm(D.T @ J2 + J2 @ D) < 1e-14


def test_common_piecewise_step_preserves_parent_wedge():
    U = np.array([1 + 2j, -0.3 + 0.4j])
    V = np.array([-0.2 + 0.7j, 0.9 - 0.1j])
    D = tracefree_symmetric(1.3, -0.8)
    M = sl2_symmetric_step(D, 0.37)
    assert abs(np.linalg.det(M) - 1.0) < 1e-12
    assert abs(parent_wedge(M @ U, M @ V) - parent_wedge(U, V)) < 1e-11


def test_wedge_and_polarization_rhs_centered_difference():
    U = np.array([0.7 + 0.2j, -0.4 + 0.1j])
    V = np.array([-0.1 + 0.8j, 0.6 - 0.3j])
    Z = np.array([0.5 - 0.4j, -0.2 + 0.7j])
    D1 = tracefree_symmetric(0.3, 0.2)
    D2 = tracefree_symmetric(-0.1, 0.5)
    D3 = tracefree_symmetric(0.4, -0.2)
    eps = 1e-7
    def step(D, X, t):
        return sl2_symmetric_step(D, t) @ X
    numw = (parent_wedge(step(D1,U,eps), step(D2,V,eps)) - parent_wedge(step(D1,U,-eps), step(D2,V,-eps))) / (2*eps)
    assert abs(numw - wedge_rhs(U,V,D1,D2)) < 1e-8
    nump = (polarization_numerator(step(D1,U,eps), step(D2,V,eps), step(D3,Z,eps)) - polarization_numerator(step(D1,U,-eps), step(D2,V,-eps), step(D3,Z,-eps))) / (2*eps)
    assert abs(nump - polarization_rhs(U,V,Z,D1,D2,D3)) < 1e-8


def test_capacity_bound():
    rng = np.random.default_rng(7)
    for _ in range(200):
        U = rng.normal(size=2)+1j*rng.normal(size=2)
        V = rng.normal(size=2)+1j*rng.normal(size=2)
        Z = rng.normal(size=2)+1j*rng.normal(size=2)
        D1 = tracefree_symmetric(*rng.normal(size=2))
        D2 = tracefree_symmetric(*rng.normal(size=2))
        D3 = tracefree_symmetric(*rng.normal(size=2))
        assert abs(polarization_rhs(U,V,Z,D1,D2,D3)) <= pointwise_capacity_bound(U,V,Z,D1,D2,D3) + 1e-12


def test_hyperbolic_distance_is_not_a_defect():
    cm = common_hyperbolic_countermodel(8.0)
    assert cm["propagator_distance"] > 1000
    assert cm["condition_number"] > 1e6
    assert cm["wedge_relative_residual"] < 1e-9


def test_child_factor_is_signed_helicity_linear_form():
    Z = np.array([3+2j, 1-4j])
    assert child_factor(Z) == (2+6j)
