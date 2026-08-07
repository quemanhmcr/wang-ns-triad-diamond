import math

import numpy as np

from src.affine_grain_dynamics import sym
from src.strain_coherence import (
    COHERENT_AVG_DEFICIT,
    COHERENT_STRAIN_TIME,
    COHERENCE_FRACTION,
    coherence_failure_action_lower,
    coherent_strain_average_deficit_lower,
    corotational_strain_rhs,
    velocity_gradient_material_rhs,
)


def test_coherent_average_cost_constant():
    T = float(COHERENT_STRAIN_TIME)
    assert math.isclose(coherent_strain_average_deficit_lower(1.0, T), float(COHERENT_AVG_DEFICIT) * T * T)


def test_coherence_failure_has_total_variation_cost():
    assert math.isclose(coherence_failure_action_lower(3.0), 3.0 * float(COHERENCE_FRACTION))


def test_velocity_gradient_and_corotational_strain_algebra():
    A = np.array([[0.3, 0.4, -0.2], [-0.1, -0.2, 0.5], [0.2, 0.1, -0.1]])
    assert abs(np.trace(A)) < 1e-14
    Hp = np.array([[0.2, 0.1, 0.0], [0.1, -0.3, 0.2], [0.0, 0.2, 0.1]])
    LA = np.array([[0.5, -0.2, 0.1], [0.3, -0.1, 0.0], [0.2, 0.4, -0.4]])
    nu = 0.7
    Arhs = velocity_gradient_material_rhs(A, Hp, LA, nu)
    S = sym(A); O = 0.5 * (A - A.T)
    materialS = sym(Arhs)
    expected = -(S @ S) - (O @ O) - Hp + nu * sym(LA)
    assert np.allclose(materialS, expected, atol=2e-13)
    objective = materialS + S @ O - O @ S
    assert np.allclose(objective, corotational_strain_rhs(A, Hp, LA, nu), atol=2e-13)
