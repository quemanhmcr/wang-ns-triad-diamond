import numpy as np
from src.affine_critical_grain import (affine_critical_mass, aspect_upper_from_geometric_scale, fresh_radius_budget, geometric_radius, incompressible_log_radius_rate, log_geometric_radius_rate, physical_covariance_rhs)


def test_geometric_radius():
    S=np.diag([4.,9.,16.])
    assert abs(geometric_radius(S)-(2*3*4)**(1/3))<1e-12


def test_affine_mass_is_radius_normalized():
    assert affine_critical_mass(.3,1.)==.3
    assert affine_critical_mass(.6,2.)==.3


def test_fresh_energy_radius_budget():
    r=np.array([.2,.3,.4]);eta=.25;E=.3*r
    sr,b=fresh_radius_budget(E,r,eta,1.)
    assert sr<=b+1e-12


def test_aspect_geometric_scale_relation():
    s=2.
    assert aspect_upper_from_geometric_scale(s)==18.


def test_incompressible_affine_strain_preserves_geometric_radius_inviscid():
    Sigma=np.array([[2.,.1,0.],[.1,1.3,.2],[0.,.2,.8]])
    A=np.array([[.4,.2,0.],[-.1,-.3,.5],[.2,0.,-.1]])
    assert abs(np.trace(A))<1e-14
    assert abs(log_geometric_radius_rate(Sigma,A,0.0))<1e-12


def test_viscosity_geometric_radius_square_lower_rate():
    Sigma=np.diag([4.,1.,.25]); nu=.7
    rg=geometric_radius(Sigma)
    rate_log=incompressible_log_radius_rate(Sigma,nu)
    rate_sq=2*rg*rg*rate_log
    assert rate_sq>=nu-1e-12
