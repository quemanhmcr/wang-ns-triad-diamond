import math
import numpy as np

from src.affine_gaussian_forcing import (
    full_quadratic_advection_residual_sq,
    gaussian_laplacian_multiplier,
    normalized_carrier,
    osculating_transverse_bound,
    osculating_transverse_residual_sq,
    quadratic_phase_tensor,
    symmetrize_rank3,
    transform_hessian,
    whitened_velocity_hessian,
)


def test_pure_quadratic_phase_is_tangent_after_projection():
    # B with zero fully symmetric part can still generate q.B chirp but no
    # third-Hermite transverse residual.
    B=np.zeros((3,3,3))
    B[0,1,2]=B[0,2,1]=1.0
    B[1,0,2]=B[1,2,0]=-0.5
    B[2,0,1]=B[2,1,0]=-0.5
    assert np.linalg.norm(symmetrize_rank3(B)) < 1e-12
    q=np.array([1.0,.3,-.2])
    assert np.linalg.norm(quadratic_phase_tensor(B,q)) > 0
    assert osculating_transverse_residual_sq(B) < 1e-24
    assert full_quadratic_advection_residual_sq(B,q) > 0


def test_transverse_bound():
    rng=np.random.default_rng(3)
    B=rng.normal(size=(3,3,3)); B=.5*(B+B.swapaxes(1,2))
    lhs=math.sqrt(osculating_transverse_residual_sq(B))
    assert lhs <= osculating_transverse_bound(B)+1e-12


def test_affine_coordinate_covariance_exact():
    rng=np.random.default_rng(4)
    L=np.array([[1.2,.1,0.],[0.,.8,.2],[.1,0.,1.1]])
    H=rng.normal(size=(3,3,3)); H=.5*(H+H.swapaxes(1,2))
    k=np.array([2.,-.4,.7])
    S=np.array([[2.,.3,0.],[0.,.7,.2],[.1,0.,1.4]])
    B=whitened_velocity_hessian(L,H); q=normalized_carrier(L,k)
    Hp=transform_hessian(S,H); Lp=S@L; kp=np.linalg.solve(S.T,k)
    assert np.linalg.norm(whitened_velocity_hessian(Lp,Hp)-B) < 1e-11
    assert np.linalg.norm(normalized_carrier(Lp,kp)-q) < 1e-12


def test_bulk_viscosity_is_gaussian_tangent_polynomial():
    G=np.array([[1.2+.1j,.05,0.],[.05,.8-.03j,.02j],[0.,.02j,1.1]],dtype=complex)
    G=.5*(G+G.T)
    k=np.array([1.1,-.4,.3]); y=np.array([.2,-.7,.5])
    # Finite-difference Laplacian check on psi itself.
    def psi(x): return np.exp(-.5*x@(G@x)+1j*k@x)
    h=2e-4
    lap=sum((psi(y+h*np.eye(3)[j])-2*psi(y)+psi(y-h*np.eye(3)[j]))/h**2 for j in range(3))
    exact=gaussian_laplacian_multiplier(G,k,y)*psi(y)
    assert abs(lap-exact)<2e-7
