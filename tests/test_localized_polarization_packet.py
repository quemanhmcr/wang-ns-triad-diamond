import math
import numpy as np

from src.localized_polarization_packet import (
    combined_localization_optimum,
    generator_freezing_bound,
    kelvin_direction_lipschitz_bound,
    kelvin_direction_rhs,
    kelvin_packet_direction_bound,
    low_strain_packet_error,
    simplified_low_strain_bound,
    unit,
)


def test_kelvin_direction_rhs_is_tangent():
    k=unit(np.array([1.0,2.0,-0.5])); A=np.array([[1.,2.,0.],[-1.,0.,3.],[.2,-.1,-1.]])
    rhs=kelvin_direction_rhs(k,A)
    assert abs(float(np.dot(k,rhs)))<1e-12


def test_kelvin_lipschitz_bound():
    A=np.array([[.2,.1,0.],[0.,-.3,.4],[.2,0.,.1]])
    B=A+0.03*np.ones((3,3))
    a=unit(np.array([1.,.2,.1])); b=unit(np.array([1.,.21,.09]))
    lhs=np.linalg.norm(kelvin_direction_rhs(a,A)-kelvin_direction_rhs(b,B))
    assert lhs<=kelvin_direction_lipschitz_bound(A,B,a,b)+1e-12


def test_scale_free_direction_bound():
    c=.4; sigma=.05; kappa=2e-4; M=10.; h=1e-3
    assert kelvin_packet_direction_bound(c,sigma,kappa,M,h) > h


def test_simplified_low_strain_bound_at_threshold():
    for c in [.1,.5,1.0]:
        sigma=1/(30*c)
        for kappa,M,h in [(1e-4,8.,1e-3),(3e-5,20.,2e-4)]:
            assert low_strain_packet_error(c,sigma,kappa,M,h) <= simplified_low_strain_bound(c,kappa,M,h)


def test_combined_optimum_exact():
    a,b,c,kappa,h=1.3,.7,.4,2e-4,1e-3
    M,E=combined_localization_optimum(a,b,c,kappa,h)
    direct=a/M+(b+7.5*c)*kappa*M+3*h
    assert math.isclose(E,direct,rel_tol=1e-13,abs_tol=1e-13)
    assert E <= a/(.9*M)+(b+7.5*c)*kappa*(.9*M)+3*h+1e-12
    assert E <= a/(1.1*M)+(b+7.5*c)*kappa*(1.1*M)+3*h+1e-12


def test_generator_bound_is_positive_homogeneous():
    x=generator_freezing_bound(.2,.3,.01)
    y=generator_freezing_bound(.4,.6,.01)
    assert y>x
