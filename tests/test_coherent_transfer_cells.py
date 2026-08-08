import numpy as np
from src.affine_coherent_moyal import periodic_discrete_stft
from src.coherent_transfer_cells import (
    affine_phase_covariance_residual, cell_energy, cell_work, discrete_cross_moyal,
    selection_jump, service_no_escape, symmetric_difference_energy,
)


def test_polarized_moyal_discrete():
    rng=np.random.default_rng(1); n=17
    f=rng.normal(size=n)+1j*rng.normal(size=n); F=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g)
    assert abs(discrete_cross_moyal(f,F,g)-np.vdot(F,f)) < 1e-11


def test_cell_work_and_energy_partition():
    rng=np.random.default_rng(2); n=19
    f=rng.normal(size=n)+1j*rng.normal(size=n); F=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g)
    Vf=periodic_discrete_stft(f,g); VF=periodic_discrete_stft(F,g); labels=np.indices(Vf.shape).sum(axis=0)%5
    assert abs(cell_work(Vf,VF,labels).sum()-2*np.vdot(F,f).real)<1e-10
    assert abs(cell_energy(Vf,labels).sum()-np.vdot(f,f).real)<1e-10


def test_affine_intrinsic_phase_is_invariant():
    rng=np.random.default_rng(3)
    M=rng.normal(size=(3,3)); M += 3*np.eye(3)
    L=rng.normal(size=(3,3)); L += 2*np.eye(3)
    assert affine_phase_covariance_residual(M,L,rng.normal(size=3),rng.normal(size=3)) < 1e-12


def test_selection_jump_is_paid_by_symmetric_difference():
    e=np.array([1.,2.,3.,4.]); a=[0,1]; b=[1,2]
    assert abs(selection_jump(e,a,b)) <= symmetric_difference_energy(e,a,b)


def test_service_no_escape_thirds():
    r=service_no_escape(3.0, .4, .5, 2.2)
    assert r['branch']=='relink_symmetric_difference'
    assert r['branch_value'] >= 1.0
