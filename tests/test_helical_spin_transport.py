import math
import numpy as np

from src.helical_spin_transport import (
    berry_chern_number,
    berry_connection_phi,
    berry_curvature_theta_phi,
    forward_normal_coupling,
    frozen_pure_strain_conversion,
    helical_with_normal,
    normal_transition_phase,
    rotation_matrix,
    spherical_helical,
    transverse_strain_helicity_matrix,
    triad_normal,
)


def test_chern_number_is_spin_one_charge():
    assert berry_chern_number(1) == -2
    assert berry_chern_number(-1) == 2
    # direct analytic flux integral: int F dtheta dphi = -4 pi s
    for s in (-1, 1):
        assert math.isclose((-4 * math.pi * s) / (2 * math.pi), berry_chern_number(s))


def test_normal_gauge_is_helical_and_rotation_covariant():
    k = np.array([1.2, -0.4, 0.8])
    n0 = np.cross(k, np.array([0.3, 0.9, -0.1]))
    n0 /= np.linalg.norm(n0)
    for s in (-1, 1):
        h = helical_with_normal(k, s, n0)
        assert np.linalg.norm(1j * np.cross(k, h) - s * np.linalg.norm(k) * h) < 1e-12
        R = rotation_matrix(np.array([0.2, 0.5, 1.0]), 0.7)
        assert np.linalg.norm(helical_with_normal(R @ k, s, R @ n0) - R @ h) < 1e-12


def test_spin_transition_is_dihedral_phase():
    k = np.array([0.0, 0.0, 2.0])
    n1 = np.array([0.0, 1.0, 0.0])
    psi = 0.73
    n2 = rotation_matrix(k, psi) @ n1
    for s in (-1, 1):
        h1 = helical_with_normal(k, s, n1)
        h2 = helical_with_normal(k, s, n2)
        assert np.linalg.norm(h2 - np.exp(1j * normal_transition_phase(k, s, n1, n2)) * h1) < 1e-12


def test_forward_normal_coupling_is_pure_quadrature():
    x = np.array([0.7, 0.1, 0.2])
    y = np.array([-0.1, 0.9, 0.3])
    z = x + y
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                g = forward_normal_coupling(x, y, z, sx, sy, sz)
                assert abs(g.real) < 1e-12


def test_berry_local_formula():
    theta, phi = 1.1, -0.4
    for s in (-1, 1):
        assert math.isclose(berry_connection_phi(theta, s), s * math.cos(theta))
        assert math.isclose(berry_curvature_theta_phi(theta, s), -s * math.sin(theta))
        h = spherical_helical(theta, phi, s)
        assert math.isclose(float(np.vdot(h, h).real), 1.0, rel_tol=1e-12)


def test_tracefree_transverse_strain_is_offdiagonal_in_helicity():
    delta, beta = 0.37, -0.22
    M = transverse_strain_helicity_matrix(delta, beta)
    assert abs(M[0, 0]) < 1e-15 and abs(M[1, 1]) < 1e-15
    assert math.isclose(abs(M[0, 1]), math.hypot(delta, beta), rel_tol=1e-12)


def test_frozen_strain_converts_helicity_exactly():
    d, t = 0.4, 0.7
    plus, minus, frac = frozen_pure_strain_conversion(d, t)
    assert math.isclose(plus, math.cosh(d * t), rel_tol=1e-12)
    assert math.isclose(minus, math.sinh(d * t), rel_tol=1e-12)
    assert 0 < frac < 0.5
