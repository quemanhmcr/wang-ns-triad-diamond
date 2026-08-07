import math
import numpy as np

from src.objective_helical_polarization import (
    coherence_second_magnus_bound,
    objective_transverse_generator,
    strain_area_commutator,
    strain_commutator,
)


def test_rigid_rotation_is_removed_by_objective_frame():
    omega = 0.7
    A = np.array([[0.0, -omega], [omega, 0.0]])
    assert np.linalg.norm(objective_transverse_generator(A)) < 1e-14


def test_strain_commutator_is_area_rotation():
    D1 = np.array([[0.3, 0.2], [0.2, -0.3]])
    D2 = np.array([[-0.1, 0.4], [0.4, 0.1]])
    assert np.linalg.norm(strain_commutator(D1, D2) - strain_area_commutator(0.3, 0.2, -0.1, 0.4)) < 1e-14


def test_repo_coherence_second_magnus_is_tiny():
    b = coherence_second_magnus_bound(1 / 20, 1 / 30)
    assert b < 1 / 10000
    assert b > 0
