import math
import numpy as np

from src.affine_shell_aspect import (
    affine_young_scaling_factors,
    aspect_mass_lower,
    local_ellipsoid_mass_coefficient,
    physical_axis_lower_constant,
    uncertainty_matrix_residual,
)


def test_uncertainty_matrix_identity():
    A=np.array([[2.,.2,0.],[.2,1.3,.1],[0.,.1,.7]])
    assert uncertainty_matrix_residual(A)<1e-12


def test_clean_axis_constant_numerically():
    assert physical_axis_lower_constant()>2/3


def test_clean_local_mass_constants_numerically():
    assert local_ellipsoid_mass_coefficient()>0.3
    assert aspect_mass_lower(1.0)>0.2


def test_affine_young_symmetry_scaling():
    for d in [1e-12, .01, 1., 20., 1e12]:
        a,b=affine_young_scaling_factors(d)
        assert math.isclose(a,1.,rel_tol=1e-13,abs_tol=1e-13)
        assert math.isclose(b,1.,rel_tol=1e-13,abs_tol=1e-13)
