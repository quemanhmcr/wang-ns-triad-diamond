import numpy as np
from src.affine_sgs_boundary_ledger import (
    cubic_boundary_pointwise_bound, affine_cubic_charge_lower_bound,
    affine_pressure_charge_lower_bound, viscous_boundary_lifetime_bound,
    clean_viscous_boundary_lifetime_bound, partition_flux_residual,
    sgs_increment_identity, sgs_increment_cubic_upper,
)


def test_cubic_boundary_young_bound():
    U=np.array([1.,-2.,.5]); W=np.array([.3,.1,-.4]); R=np.arange(9,dtype=float).reshape(3,3)/10
    lhs,rhs=cubic_boundary_pointwise_bound(U,W,R)
    assert lhs<=rhs+1e-12


def test_affine_charge_scales_with_moat_width():
    assert affine_cubic_charge_lower_bound(2.,10.,1.,8.) > affine_cubic_charge_lower_bound(2.,10.,1.,4.)
    assert affine_pressure_charge_lower_bound(2.,10.,1.,8.) > 0


def test_clean_viscous_bound_dominates_delta_one_twentieth():
    exact=viscous_boundary_lifetime_bound(.3,.2,.05,1.4,10.,3.)
    clean=clean_viscous_boundary_lifetime_bound(.3,.2,1.4,10.,3.)
    assert exact<clean


def test_partition_flux_cancels_when_weights_sum_to_one():
    dt=np.array([.2,-.1,-.1]); G=np.array([[1.,2.,3.],[-.5,.1,-1.],[ -.5,-2.1,-2.]])
    assert np.linalg.norm(G.sum(axis=0))<1e-12
    assert abs(partition_flux_residual(dt,G,np.array([2.]),np.array([.3,-.7,1.1])))<1e-12


def test_sgs_increment_identity_and_cubic_bound():
    w=np.array([.8,.4,-.2]); du=np.array([[1.,0.,.2],[-.3,.7,.1],[.2,-.4,1.1]])
    R=sgs_increment_identity(w,du)
    g1=np.sum(np.abs(w)); cubic=np.sum(np.abs(w)*np.linalg.norm(du,axis=1)**3)
    assert np.linalg.norm(R,'fro')**1.5 <= sgs_increment_cubic_upper(g1,cubic)+1e-12
