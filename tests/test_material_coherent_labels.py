import math
import numpy as np

from src.material_coherent_labels import (
    affine_transport_zeta_residual,
    covariance_representation_xi_upper,
    dyadic_address,
    frequency_representation_xi_upper,
    geometric_schedule_sum_upper,
    intrinsic_zeta,
    parent_address,
    refinement_energy_residual,
)
from src.smooth_symbol_freezing import sharp_young_constant_3d


def test_common_affine_transport_preserves_intrinsic_label():
    L=np.diag([1.2,.8,1.1]); X=np.array([1.,-2.,.5]); k=np.array([3.,1.,-1.])
    M=np.array([[1.,.2,0.],[0.,1.,.1],[.1,0.,1.]])
    assert affine_transport_zeta_residual(M,L,X,k)<1e-13


def test_dyadic_addresses_are_nested():
    z=np.array([-.1,.2,1.7,-2.1,.8,.03])
    for j in range(8):
        assert parent_address(dyadic_address(z,j+1,.7))==dyadic_address(z,j,.7)


def test_refinement_has_zero_moyal_switch_charge():
    child=[.1,.4,.2,.3]
    assert abs(refinement_energy_residual(sum(child),child))<1e-15


def test_geometric_representation_budgets():
    A3=sharp_young_constant_3d()
    assert math.isclose(geometric_schedule_sum_upper(.01),.02)
    assert math.isclose(frequency_representation_xi_upper(.01,2.,3.),2*A3*2*3*.01)
    assert math.isclose(covariance_representation_xi_upper(.01,4.),math.sqrt(2)*.01*4.)


def test_intrinsic_zeta_shape():
    z=intrinsic_zeta(np.eye(3),np.ones(3),2*np.ones(3))
    assert z.shape==(6,)
