import numpy as np
from src.hermite_helicity_ledger import sym3, third_hermite_norm_sq, affine_envelope_sideband_norm_sq, first_hermite_sideband_norm_sq, h3_projection_onto_degree2, h1_projection_onto_base


def test_h3_exact_norm_coefficient():
    T=np.zeros((3,3,3)); T[0,0,0]=1
    assert abs(third_hermite_norm_sq(T)-6)<1e-12
    assert abs(affine_envelope_sideband_norm_sq(T)-3/8)<1e-12


def test_h3_orthogonal_to_quadratic_tangent():
    T=sym3(np.arange(27,dtype=float).reshape(3,3,3)/10)
    C2=np.array([[1.,.2,.1],[.2,-.3,.4],[.1,.4,.7]])
    p=h3_projection_onto_degree2(T,(.7,np.array([.2,-.1,.4]),C2),order=6)
    assert abs(p)<1e-10


def test_h1_has_zero_base_projection_and_exact_norm():
    C=np.array([[1.,2.,-1.],[.5,-.3,.2]])
    assert np.linalg.norm(h1_projection_onto_base(C,np.array([1.,0.]),order=4))<1e-12
    assert abs(first_hermite_sideband_norm_sq(C)-np.sum(C*C))<1e-12
