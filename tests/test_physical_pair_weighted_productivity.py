import math

import numpy as np

from src.physical_pair_weighted_productivity import (
    generated_energy_sup_factor,
    physical_log_productivity_constant,
    tempered_pair_cell_penalty_upper,
    variable_productivity_root_log_lower,
)


def test_generated_energy_sup_is_below_two_positive_work_units():
    assert generated_energy_sup_factor() < 2.0


def test_physical_productivity_decays_only_as_inverse_pair_count():
    a=physical_log_productivity_constant(1)
    b=physical_log_productivity_constant(7)
    assert a>0
    assert math.isclose(b,a/7.0)


def test_polynomial_pair_refinement_has_finite_binary_discounted_penalty():
    M0=5.0; p=12.0; L=80
    j=np.arange(L,dtype=float)
    M=M0*(j+3.0)**p
    weights=2.0**(-(j+1.0))
    actual=float(np.sum(weights*np.log(M)))
    assert actual <= tempered_pair_cell_penalty_upper(M0,p)+1e-12


def test_variable_productivity_recursion_matches_closed_form():
    lam=np.array([0.1,0.05,0.02,0.01])
    ellT=-2.0
    got=variable_productivity_root_log_lower(lam,ellT)
    weights=2.0**(-(np.arange(len(lam))+1.0))
    want=float(np.sum(weights*np.log(lam))+(2.0**(-len(lam)))*ellT)
    assert math.isclose(got,want,rel_tol=1e-13,abs_tol=1e-13)
