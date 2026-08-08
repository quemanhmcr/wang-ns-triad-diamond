import math

import numpy as np

from src.physical_energy_causal_bridge import (
    PHYSICAL_HH_WORK_FRACTION,
    adjoint_projection_work_decomposition,
    flat_scalar_measure_countermodel,
    heavy_half_physical_transfer,
    physical_hh_work_lower_bound,
    route_physical_energy_causality,
)


def test_flat_scalar_raw_duhamel_and_physical_transfer_laws_are_not_equal():
    out = flat_scalar_measure_countermodel(0.5)
    assert math.isclose(out["duhamel_cdf"], 0.5)
    assert math.isclose(out["physical_transfer_cdf"], 0.25)
    assert math.isclose(out["cdf_gap"], 0.25)


def test_adjoint_projection_identity_splits_actual_physical_work_exactly():
    c = np.array([1 + 2j, -0.5 + 0.25j])
    psi = np.array([0.75 - 0.2j, 1.1 + 0.4j])
    F = np.array([-0.3 + 0.7j, 0.9 - 0.1j])
    out = adjoint_projection_work_decomposition(c, psi, F)
    assert abs(out["residual"]) < 1e-13
    assert math.isclose(
        out["physical_work"],
        out["response_energy_work"] + out["orthogonal_child_work"],
        rel_tol=1e-13,
        abs_tol=1e-13,
    )


def test_low_strain_energy_gate_forces_clean_physical_high_high_work():
    out = route_physical_energy_causality(
        terminal_energy=1.0,
        initial_energy=0.19,
        residual_positive_work=0.19,
        strain_action=1 / 30,
    )
    assert out["branch"] == "physical_high_high_transfer_generation"
    assert out["physical_hh_work_lower"] >= float(PHYSICAL_HH_WORK_FRACTION) - 1e-14


def test_energy_gate_routes_inheritance_residual_and_high_strain_before_generation():
    assert route_physical_energy_causality(
        terminal_energy=1.0, initial_energy=0.2, residual_positive_work=0.0, strain_action=0.0
    )["branch"] == "material_energy_inheritance"
    assert route_physical_energy_causality(
        terminal_energy=1.0, initial_energy=0.0, residual_positive_work=0.2, strain_action=0.0
    )["branch"] == "classified_residual_physical_work"
    assert route_physical_energy_causality(
        terminal_energy=1.0, initial_energy=0.0, residual_positive_work=0.0, strain_action=0.034
    )["branch"] == "high_strain_critical_dissipation"


def test_linear_gronwall_bound_has_exact_clean_constant_at_threshold():
    low = physical_hh_work_lower_bound(
        terminal_energy=15.0,
        initial_energy=3.0,
        residual_positive_work=3.0,
        strain_action=1 / 30,
    )
    assert math.isclose(low, 8.0)


def test_synchronization_uses_arbitrary_positive_physical_transfer_weights():
    out = heavy_half_physical_transfer(
        times=[0.05, 0.1, 0.8, 0.9],
        positive_work_weights=[1.0, 1.0, 10.0, 2.0],
        slab_start=0.0,
        slab_end=1.0,
    )
    assert out["half"] == 1
    assert out["mass"] >= 0.5 * out["total"]
    assert math.isclose(out["normalized_parent_span_upper"], 25 / 128)
