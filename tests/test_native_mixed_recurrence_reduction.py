import math

import pytest

from src.native_closed_triad_no_free_power_branching import POWER_LOG_COST_LOWER
from src.native_mixed_recurrence_reduction import (
    conditional_master_power_bound,
    theorem_certificate,
)


def test_conditional_master_formula_uses_only_free_power_vertices():
    out = conditional_master_power_bound(17, 5, finite_summable_prefactor=2.25)
    assert out.guaranteed_free_power_vertices == 12
    assert out.remainder_upper == pytest.approx(2.25 * (10.0 / 13.0) ** 12)
    assert out.cyclic_log_cost_lower == pytest.approx(12 * POWER_LOG_COST_LOWER)
    assert out.exact_action_power_ties_joint
    assert not out.power_cocharge_on_action_vertices
    assert not out.master_composition_certified


def test_all_action_containing_vertices_remove_power_cocharge():
    out = conditional_master_power_bound(9, 9)
    assert out.guaranteed_free_power_vertices == 0
    assert out.remainder_upper == pytest.approx(1.0)
    assert out.cyclic_log_cost_lower == pytest.approx(0.0)


def test_action_upper_bound_larger_than_depth_gives_no_forced_power_charge():
    out = conditional_master_power_bound(4, 5)
    assert out.guaranteed_free_power_vertices == 0
    assert out.remainder_upper == pytest.approx(1.0)


def test_negative_action_upper_bound_is_rejected():
    with pytest.raises(ValueError):
        conditional_master_power_bound(4, -1)


def test_certificate_places_stock_power_action_above_smaller_intrinsic_law_set():
    cert = theorem_certificate()
    assert set(cert["native_dependencies"]) == {
        "closed_triad_current",
        "mode_stock_continuity",
        "local_action_speed_lock",
    }
    laws = cert["primitive_law_set"]
    assert set(laws) == {
        "I_closed_triad_current",
        "II_mode_stock_continuity",
        "III_local_action_speed_lock",
    }
    assert "recurrence ontology" in cert["ontology_status"]
    assert "10/13" in cert["candidate_master_bound"]
    assert "joint causes" in cert["tie_policy"]
    assert "conditional" in cert["master_composition_status"]
    assert "no FIFO/LIFO/proportional temporal matching" in cert["temporal_guard"]
    assert cert["global_regularity_claimed"] is False
