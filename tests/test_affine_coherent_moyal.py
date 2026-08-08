import numpy as np
from src.affine_coherent_moyal import discrete_moyal_residual, exact_moyal_certificate, moyal_energy, phase_cell_budget

def test_moyal_formula():
    assert moyal_energy(3.0, 1.0) == 3.0

def test_discrete_periodic_moyal():
    f=np.array([1+1j,2,-.3j,1.2],complex); g=np.array([1.,2.,-1.,.5],complex); g/=np.linalg.norm(g)
    assert abs(discrete_moyal_residual(f,g)) < 1e-12

def test_positive_cells_sum():
    assert phase_cell_budget([.2,.3,.5]) == 1.0

def test_certificate():
    assert exact_moyal_certificate()['normalized_window_budget'] == 'P=1'
