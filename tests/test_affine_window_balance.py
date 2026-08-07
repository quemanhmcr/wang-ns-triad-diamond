import numpy as np
from src.affine_window_balance import affine_balance_optimum, clean_shell_gradient_bound, material_window_leakage_bound, normalized_material_coordinate_rhs, quadratic_velocity


def test_exact_material_coordinate_remainder_for_quadratic_flow():
    rng=np.random.default_rng(7);L=np.array([[2.,.1,0.],[.1,.8,0.],[0.,0.,1.2]]);U0=rng.normal(size=3);A=rng.normal(size=(3,3));H=rng.normal(size=(3,3,3));H=.5*(H+H.swapaxes(1,2));z=np.array([.4,-.2,.7]);y=L@z;Ux=quadratic_velocity(U0,A,H,y)
    rhs=normalized_material_coordinate_rhs(L,Ux,U0,A,y)
    B=np.einsum('ai,ijk,jb,kc->abc',np.linalg.inv(L),H,L,L)
    expected=.5*np.einsum('abc,b,c->a',B,z,z)
    assert np.linalg.norm(rhs-expected)<1e-11


def test_clean_shell_gradient_scale():
    assert clean_shell_gradient_bound(100.,10.,1.)==.15


def test_affine_balance_optimizer():
    M,E=affine_balance_optimum(2.,3.,.04)
    assert abs(E-(2/M+3*.04*M))<1e-12
    assert material_window_leakage_bound(.1,4.,1.,2.)==.8
