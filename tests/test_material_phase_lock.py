import numpy as np

from src.material_phase_lock import hessian_phase_rhs, kelvin_rhs


def test_common_kelvin_resonance_is_exact():
    A=np.array([[.2,.1,0.],[-.3,.4,.2],[0.,-.1,-.6]])
    k1=np.array([1.,2.,-.3]);k2=np.array([-.2,.4,1.]);k3=k1+k2
    assert np.linalg.norm(kelvin_rhs(k1,A)+kelvin_rhs(k2,A)-kelvin_rhs(k3,A))<1e-12


def test_common_hessian_chirp_cancels_at_resonance():
    rng=np.random.default_rng(4)
    A=rng.normal(size=(3,3));H=rng.normal(size=(3,3,3));H=.5*(H+H.swapaxes(1,2))
    k1=rng.normal(size=3);k2=rng.normal(size=3);k3=k1+k2
    K1=rng.normal(size=(3,3));K1=.5*(K1+K1.T)
    K2=rng.normal(size=(3,3));K2=.5*(K2+K2.T);K3=K1+K2
    lock=hessian_phase_rhs(K1,k1,A,H)+hessian_phase_rhs(K2,k2,A,H)-hessian_phase_rhs(K3,k3,A,H)
    assert np.linalg.norm(lock)<1e-11
