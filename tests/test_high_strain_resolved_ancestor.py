import math

import numpy as np

from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower
from src.high_strain_resolved_ancestor import (
    ancestor_lifetime_ratio,
    dyadic_resolved_shell_upper,
    high_strain_ancestor_mass_threshold,
    infinite_resolved_shell_upper_sum,
    low_mass_dissipation_upper,
    retained_fraction_lower,
    shell_dissipation_law,
)


def test_resolved_dyadic_upper_radii_have_exact_geometric_sum():
    N=16.0
    partial=sum(dyadic_resolved_shell_upper(N,j) for j in range(60))
    assert math.isclose(infinite_resolved_shell_upper_sum(N),N/2)
    assert abs(partial-N/2)<1e-15*N


def test_clean_high_strain_threshold_leaves_at_least_half_dissipation_law():
    c=1.7
    Dstar=clean_high_strain_dissipation_lower(c)
    mustar=high_strain_ancestor_mass_threshold(c)
    assert math.isclose(low_mass_dissipation_upper(10.0,c,mustar),Dstar/2)
    assert retained_fraction_lower(Dstar,c,mustar)>=0.5-1e-14
    assert retained_fraction_lower(3*Dstar,c,mustar)>=5/6-1e-14


def test_finite_shell_physical_law_obeys_geometric_bad_mass_bound():
    N=32.0; c=1.0
    mustar=high_strain_ancestor_mass_threshold(c)
    # Two shells, three time cells; V is an arbitrary contraction of u.
    E=np.array([[8.0,0.2,12.0],[20.0,0.1,0.3]])
    a=np.array([[1.0,.7,.4],[.8,.2,.9]])
    rho=np.array([[1.0,.6,.4],[.7,.3,.9]])
    out=shell_dissipation_law(N,c,E,a,rho,critical_mass_threshold=mustar)
    assert out['bad_upper_margin']>=-1e-13
    assert 0<=out['retained_fraction']<=1


def test_resolved_ancestor_natural_lifetime_is_at_least_sixteen_child_lifetimes():
    N=100.0
    for j in range(12):
        M=dyadic_resolved_shell_upper(N,j)
        assert ancestor_lifetime_ratio(N,M)>=16.0
