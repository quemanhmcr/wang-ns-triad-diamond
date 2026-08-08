import math

from src.complex_young_parent_marking import (
    christ_complex_parent_mark_available,
    complex_parent_marking_budget,
    complex_young_deficit_upper,
    complex_young_efficiency_lower,
    convolution_pair_efficiency_lower,
)


def test_weighted_physical_transfer_reduces_to_complex_unweighted_young():
    assert math.isclose(complex_young_efficiency_lower(0.99, 0.01), 0.98)
    assert math.isclose(convolution_pair_efficiency_lower(0.99, 0.01), 0.98)
    assert math.isclose(complex_young_deficit_upper(0.01, 0.01), 0.02)


def test_no_lower_bound_on_frozen_symbol_modulus_is_needed_in_reduction():
    # The proof only uses |m0|<=m_*; weighted near saturation itself rules out
    # a tiny m0 unless the Xi error is already large.
    assert complex_young_efficiency_lower(1.0, 0.001) == 0.999


def test_external_christ_modulus_is_used_as_a_gate_not_invented_numerically():
    assert christ_complex_parent_mark_available(
        weighted_deficit=0.002,
        normalized_symbol_freezing_error=0.001,
        christ_modulus_for_target_distance=0.004,
    )
    assert not christ_complex_parent_mark_available(
        weighted_deficit=0.003,
        normalized_symbol_freezing_error=0.002,
        christ_modulus_for_target_distance=0.004,
    )


def test_symbol_lipschitz_schedule_enters_once_as_normalized_xi():
    out = complex_parent_marking_budget(
        weighted_deficit=0.001,
        symbol_lipschitz_constant=2.0,
        relative_cell_diameter=0.0005,
        christ_modulus_for_target_distance=0.003,
    )
    assert math.isclose(out["normalized_symbol_freezing_error"], 0.001)
    assert math.isclose(out["complex_young_deficit_upper"], 0.002)
    assert out["complex_parent_mark_available"] is True
