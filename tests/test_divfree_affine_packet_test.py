import numpy as np
from src.divfree_affine_packet_test import leray_matrix, pressure_projection_residual, affine_shell_commutator_coefficient, discrete_convolution_zero, discrete_commutator, discrete_first_moment, discrete_lipschitz


def test_leray_kills_pressure_gradient():
    k=np.array([1.,2.,-.5]); assert pressure_projection_residual(k,2+.3j)<1e-12
    P=leray_matrix(k); assert np.linalg.norm(P@k)<1e-12


def test_affine_shell_commutator_coefficient():
    assert abs(affine_shell_commutator_coefficient(2.,3.,6.)-1.5)<1e-12


def test_discrete_commutator_young_bound():
    ker={-1:.2,0:.7,1:-.1}; f=np.arange(10,dtype=float)+1j*np.arange(10,dtype=float)[::-1]; chi=.03*np.arange(10)
    lhs=np.linalg.norm(discrete_commutator(ker,chi,f)); rhs=discrete_first_moment(ker)*discrete_lipschitz(chi)*np.linalg.norm(f)
    assert lhs<=rhs+1e-12
