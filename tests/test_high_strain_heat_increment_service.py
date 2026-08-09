import math

import numpy as np

from src.high_strain_heat_increment_service import (
    HEAT_DEFECT_LOWER,
    critical_ancestor_heat_service_fraction_lower,
    heat_increment_multiplier,
    heat_to_gradient_ratio,
    high_strain_heat_service_lower,
    retained_heat_service_bounds,
    spectral_heat_service,
)
from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower


def test_heat_increment_multiplier_is_exact_heat_defect_symbol():
    N=7.0; k=1.25
    expected=2*N*N*(1-math.exp(-k*k/(2*N*N)))
    assert math.isclose(heat_increment_multiplier(k,N),expected,rel_tol=2e-15)


def test_resolved_heat_defect_is_uniformly_comparable_to_gradient():
    N=13.0
    for k in np.linspace(0,N/4,101):
        ratio=heat_to_gradient_ratio(float(k),N)
        assert HEAT_DEFECT_LOWER-2e-15<=ratio<=1+2e-15


def test_spectral_heat_service_bounds_arbitrary_resolved_energy():
    N=20.0
    k=np.array([0.1,1.0,2.5,4.9])
    e=np.array([4.0,.2,10.0,2.0])
    out=spectral_heat_service(N,k,e)
    assert out['lower_margin']>=-1e-13
    assert out['upper_margin']>=-1e-13


def test_high_strain_forces_fixed_heat_increment_service():
    c=.8
    D=clean_high_strain_dissipation_lower(c)
    assert math.isclose(high_strain_heat_service_lower(c),HEAT_DEFECT_LOWER*D)


def test_half_dissipation_ancestor_retains_nearly_half_heat_service():
    D=3.0
    out=retained_heat_service_bounds(D,D/2)
    assert math.isclose(out['fraction_lower'],critical_ancestor_heat_service_fraction_lower())
    assert out['fraction_lower']>.48
