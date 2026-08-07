import numpy as np
from src.affine_curvature_connection import curvature_connection_residual, pressure_third_far_shell_exponent, far_pressure_third_geometric_sum


def test_curvature_connection_exact_algebra():
    A=np.array([[.2,.1,0.],[-.3,.1,.2],[.1,-.2,-.3]]); A-=np.trace(A)/3*np.eye(3)
    H=np.arange(27,dtype=float).reshape(3,3,3)/20; H=.5*(H+H.swapaxes(1,2))
    F2=np.arange(27,dtype=float)[::-1].reshape(3,3,3)/30; F2=.5*(F2+F2.swapaxes(1,2))
    L=np.array([[1.2,.1,0.],[0.,.8,.2],[.1,0.,1.4]])
    assert np.linalg.norm(curvature_connection_residual(A,H,F2,L))<1e-11


def test_pressure_third_far_exponent():
    assert pressure_third_far_shell_exponent()==3
    assert far_pressure_third_geometric_sum(3)>0
