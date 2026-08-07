import numpy as np
from src.affine_critical_grain import affine_critical_mass, aspect_upper_from_geometric_scale, fresh_radius_budget, geometric_radius


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
