import math

import numpy as np

from src.smooth_material_carrier_relay import (
    coefficient_energy_margin,
    hahn_weighted_generation_chain,
    orthogonal_hard_event_resolution,
    registered_coefficient_energy_lower,
    relay_energy_route,
    two_step_affine_reanchor_zeta_residual,
)


def test_registered_coefficient_is_energy_of_same_carrier_not_new_hard_role():
    w=np.array([1+2j,-3+1j,2-1j],complex)
    psi=np.array([2-1j,1+0j,-1+2j],complex)
    z=np.vdot(psi,w)
    assert math.isclose(registered_coefficient_energy_lower(z,np.linalg.norm(psi)),abs(z)**2/np.vdot(psi,psi).real)
    assert coefficient_energy_margin(w,psi)>=-1e-13


def test_two_step_common_affine_reanchor_preserves_intrinsic_material_label():
    M1=np.array([[1.1,.2,0],[0,.9,.1],[.1,0,1.02]])
    M2=np.array([[.95,0,.1],[.1,1.05,0],[0,.05,1.]])
    L=np.array([[1.2,.1,0],[0,.8,.2],[.1,0,1.1]])
    X=np.array([.3,-.7,1.2]); k=np.array([2.,-1.,.5])
    assert two_step_affine_reanchor_zeta_residual(M1,M2,L,X,k)<1e-12


def test_smooth_carrier_positive_generation_is_dominated_by_physical_hahn_work():
    r=np.array([4.,-10.,3.,-1.,2.])
    q=np.array([1.,.2,.6,.9,.5])
    out=hahn_weighted_generation_chain(r,q)
    assert out['hahn_margin']>=-1e-13
    assert out['contraction_margin']>=-1e-13
    assert out['carrier_positive_generation']<=out['physical_hahn_positive']+1e-13


def test_hard_resolution_occurs_exactly_at_work_event():
    w=np.array([1+1j,2-1j,-1+.5j,3+0j])
    Fs=[np.array([1,2j,0,-1],complex),np.array([-.5,1,2,-2j],complex)]
    P1=np.diag([1,1,0,0]).astype(complex)
    P2=np.diag([0,0,1,1]).astype(complex)
    out=orthogonal_hard_event_resolution(w,Fs,(P1,P2))
    assert abs(out['reconstruction_residual'])<1e-13
    assert out['hard_positive_mass']+1e-13>=out['aggregate_positive_work']


def test_previous_slab_uses_existing_energy_gate_on_same_carrier():
    out=relay_energy_route(terminal_carrier_energy=1.0,initial_carrier_energy=.1,residual_positive_work=.1,strain_action=.01)
    assert out['branch']=='physical_high_high_transfer_generation'
    assert out['physical_hh_work_lower']>=8/15-1e-14
