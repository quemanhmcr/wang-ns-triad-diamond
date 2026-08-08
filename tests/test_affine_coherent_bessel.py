import numpy as np
import pytest
from src.affine_coherent_bessel import (
    arb_clean_bessel_certificate, coherent_overlap_magnitude, gram_matrix,
    intrinsic_phase_point, overlap_from_intrinsic_points, synthesis_coefficient_energy_upper,
)

def test_overlap_intrinsic_formula():
    L=np.diag([2.,1.,.5]);X=np.array([1.,2.,3.]);Y=np.array([-.3,.4,1.]);k=np.array([.2,-.1,.7]);l=np.array([-.4,.3,.1])
    z=intrinsic_phase_point(X,k,L);e=intrinsic_phase_point(Y,l,L)
    assert abs(coherent_overlap_magnitude(X,k,Y,l,L)-overlap_from_intrinsic_points(z,e))<1e-14

def test_common_affine_invariance():
    L=np.array([[1.2,.2,0],[0,.8,.1],[.1,0,1.1]]);S=np.array([[1.1,.2,0],[-.1,.9,.3],[0,.1,1.]])
    X=np.array([1.,.2,-.3]);k=np.array([.4,-.2,.1])
    assert np.linalg.norm(intrinsic_phase_point(S@X,np.linalg.solve(S.T,k),S@L)-intrinsic_phase_point(X,k,L))<1e-13

def test_simple_gram_positive_and_bounded():
    Z=np.array([[0,0,0,0,0,0],[4,0,0,0,0,0],[0,4,0,0,0,0]],float);G=gram_matrix(Z)
    assert np.linalg.eigvalsh(G)[0]>0
    assert np.linalg.eigvalsh(G)[-1]<25/4

def test_arb_optional():
    pytest.importorskip('flint');c=arb_clean_bessel_certificate();assert c['clean_Bessel_constant']=='25/4'


def test_clean_synthesis_budget():
    assert synthesis_coefficient_energy_upper(47.0) == pytest.approx(50.0)
