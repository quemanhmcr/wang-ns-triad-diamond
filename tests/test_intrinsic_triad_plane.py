import math
import numpy as np

from src.intrinsic_triad_plane import (
    gaussian_effective_carrier_driver,
    gram_rhs,
    intrinsic_coercivity_ratio,
    oriented_extremal_parents,
    plane_normal_rhs,
    restricted_shape_driver,
)


def test_gram_law_under_general_3d_driver():
    E=np.array([[1.,0.],[0.,1.],[0.,0.]])
    K=oriented_extremal_parents(E)
    B=np.array([[0.2,0.3,-0.4],[-0.1,-0.3,0.2],[0.4,0.1,0.1]])
    kd=-B@K
    assert np.allclose(gram_rhs(K,B),kd.T@K+K.T@kd,atol=2e-13)


def test_full3d_coercivity_ignores_extrinsic_tilt():
    th=0.71
    E=np.array([[math.cos(th),0.], [0.,1.], [math.sin(th),0.]])
    B=np.array([[1.0,0.4,2.0],[0.4,-1.0,-3.0],[5.0,4.0,0.2]])
    assert intrinsic_coercivity_ratio(B,E)>0.43


def test_plane_normal_evolution_preserves_orthogonality():
    E=np.array([[1.,0.],[0.,1.],[0.,0.]])
    n=np.array([0.,0.,1.])
    K=oriented_extremal_parents(E)
    B=np.array([[0.2,0.3,-0.4],[-0.1,-0.3,0.2],[0.4,0.1,0.1]])
    nd=plane_normal_rhs(n,B); kd=-B@K
    assert np.linalg.norm(nd@K+n@kd)<2e-13


def test_isotropic_viscosity_is_intrinsic_shape_gauge_in_any_plane():
    q,_=np.linalg.qr(np.array([[1.,2.,3.],[2.,-1.,1.],[1.,1.,-1.]]))
    E=q[:,:2]
    A=np.array([[0.3,0.2,-0.1],[0.0,-0.1,0.4],[0.2,-0.3,-0.2]])
    P=1.7*np.eye(3)
    B=gaussian_effective_carrier_driver(A,P,0.9)
    assert np.allclose(restricted_shape_driver(B,E),restricted_shape_driver(A.T,E),atol=2e-13)
