import numpy as np

from src.affine_gaussian_forcing import symmetrize_rank3
from src.quadratic_swirl_kernel import (
    divergence_gradient,
    gaussian_envelope_advection,
    quadratic_velocity,
    reconstruct_M_from_kernel,
    swirl_tensor,
    swirl_velocity,
)


def test_swirl_kernel_exact_properties():
    M=np.array([[1.,.2,.3],[.2,-.4,.1],[.3,.1,-.6]])
    B=swirl_tensor(M)
    assert np.linalg.norm(symmetrize_rank3(B))<1e-12
    assert np.linalg.norm(divergence_gradient(B))<1e-12
    assert np.linalg.norm(reconstruct_M_from_kernel(B)-M)<1e-12


def test_swirl_velocity_is_radial_tangent():
    M=np.diag([1.,-.2,-.8]); z=np.array([.7,-1.1,.4])
    B=swirl_tensor(M)
    assert np.linalg.norm(quadratic_velocity(B,z)-swirl_velocity(M,z))<1e-12
    assert abs(gaussian_envelope_advection(M,z))<1e-12
