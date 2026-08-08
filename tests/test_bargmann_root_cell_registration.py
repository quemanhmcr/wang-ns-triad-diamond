import math

from src.bargmann_root_cell_registration import (
    bargmann_ball_energy_fraction,
    canonical_cell_critical_mass_lower,
    canonical_cell_frame_budget,
    actual_canonical_cell_quantum,
    normalized_canonical_cell_quantum,
    deterministic_energy_anchor,
    optimal_bargmann_fraction,
    optimal_bargmann_radius,
    pushforward_parent_slot_weights,
    unit_grid_cells_intersecting_ball_upper,
)
from src.dual_gaussian_root_registration import normalized_dual_probe_critical_mass_lower


def test_bargmann_optimal_ball_fraction_is_exact():
    R = math.sqrt(3.0)
    assert math.isclose(optimal_bargmann_radius(), R)
    assert math.isclose(optimal_bargmann_fraction(), 9.0 / (2.0 * math.e**3))
    assert math.isclose(bargmann_ball_energy_fraction(R), optimal_bargmann_fraction())


def test_sqrt3_ball_meets_at_most_five_cells_per_phase_coordinate():
    assert unit_grid_cells_intersecting_ball_upper(math.sqrt(3.0), 6) == 5**6


def test_dual_probe_quantum_forces_positive_canonical_material_cell_quantum():
    eta_probe = normalized_dual_probe_critical_mass_lower()
    eta_cell = canonical_cell_critical_mass_lower(eta_probe)
    assert eta_cell == normalized_canonical_cell_quantum()
    assert eta_cell > 0


def test_energy_anchor_is_deterministic_and_actual():
    cells = {(2,): 0.5, (0,): 0.7, (1,): 0.7}
    assert deterministic_energy_anchor(cells) == (0,)


def test_positive_parent_law_pushforward_preserves_mass_and_merges_reuse():
    weights = {"a": 0.2, "b": 0.3, "c": 0.5}
    labels = {"a": (1,), "b": (1,), "c": (2,)}
    out = pushforward_parent_slot_weights(weights, labels)
    assert math.isclose(sum(out.values()), 1.0)
    assert math.isclose(out[(1,)], 0.5)
    assert math.isclose(out[(2,)], 0.5)


def test_cell_frame_budget_is_positive_and_depth_independent_by_construction():
    assert canonical_cell_frame_budget() > 0


def test_actual_material_cell_quantum_scales_with_squared_parent_amplitude():
    eta = normalized_canonical_cell_quantum()
    assert math.isclose(actual_canonical_cell_quantum(2.5), 6.25 * eta)
    assert actual_canonical_cell_quantum(0.0) == 0.0
