import numpy as np

from src.divfree_coherent_parseval import (
    parseval_projected_energy,
    pressure_pairing,
    projected_coefficient,
    unprojected_coefficient,
)


def test_projected_coefficients_equal_on_projected_subspace():
    P=np.diag([1.,1.,0.]).astype(complex)
    u=np.array([1+2j,-.5j,0.])
    phi=np.array([.2j,1-.3j,4+2j])
    assert abs(projected_coefficient(P,phi,u)-unprojected_coefficient(phi,u))<1e-14


def test_projected_orthonormal_basis_is_parseval_on_subspace():
    P=np.diag([1.,1.,0.]).astype(complex)
    F=np.eye(3,dtype=complex)
    u=np.array([1+1j,2-.5j,0.])
    assert abs(parseval_projected_energy(P,F,u)-np.vdot(u,u).real)<1e-14


def test_gradient_pressure_is_orthogonal():
    P=np.diag([1.,1.,0.]).astype(complex)
    phi=np.array([1.,2.,3.],complex)
    grad=np.array([0.,0.,4+2j])
    assert abs(pressure_pairing(P,phi,grad))<1e-14
