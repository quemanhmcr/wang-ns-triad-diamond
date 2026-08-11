import math

import pytest

from src.full_natural_service_corridor_quotient import (
    FULL_NATURAL_SERVICE_WITNESS,
    STATUS,
    endpoint_comparable_hard_shell_cover,
    material_partition_is_same_corridor_measure,
    quotient_full_natural_service_outcome,
    realized_endpoint_hard_shell_witnesses,
    theorem_certificate,
)


def _outcome(A: float, c: float, y: float = 0.3, *, event_time: float | None = None):
    T = c / A**2
    t = 4.0 * T if event_time is None else float(event_time)
    return {
        "classification": FULL_NATURAL_SERVICE_WITNESS,
        "joint_first_stops": (),
        "required_elapsed": T,
        "observed_elapsed_end": T,
        "uniform_square_service_lower": y,
        "integrated_bounded_heat_service_lower": c * y,
        "endpoint_carrier_critical_mass_lower": 2.0 * y,
        "corridor_terminal_time": t,
        "corridor_endpoint_time": t - T,
        "corridor_endpoint_elapsed_from_terminal": T,
        "physical_time_drop": T,
        "renewal_frequency": A,
        "scaled_lifetime": c,
        "parent_shell_frequency": A / 0.75,
        "service_same_corridor_witness": True,
        "service_adds_recursion_depth": False,
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def test_full_natural_service_is_attached_to_the_already_completed_corridor():
    A, c = 3.0, 1.2
    T = c / A**2
    corridor = quotient_full_natural_service_outcome(
        _outcome(A, c),
        event_time=4.0 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
    )
    assert corridor.physical_time_drop == pytest.approx(T)
    assert corridor.endpoint_time == pytest.approx(3.0 * T)


def test_t0_root_cannot_be_retyped_as_full_natural_service():
    A, c = 2.0, 1.0
    T = c / A**2
    with pytest.raises(ValueError, match="absorbing boundary"):
        quotient_full_natural_service_outcome(
            _outcome(A, c, event_time=T),
            event_time=T,
            renewal_frequency=A,
            scaled_lifetime=c,
        )


def test_material_rereading_partitions_same_service_without_recursion_depth():
    A, c = 4.0, 0.8
    T = c / A**2
    corridor = quotient_full_natural_service_outcome(
        _outcome(A, c, event_time=3.0 * T),
        event_time=3.0 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
    )
    out = material_partition_is_same_corridor_measure(
        corridor,
        [1.0, 2.0, 3.0, 4.0],
        [True, True, False, False],
        [True, False, True, False],
    )
    assert out["old_old"] == pytest.approx(1.0)
    assert out["old_new_interface"] == pytest.approx(5.0)
    assert out["new_new"] == pytest.approx(4.0)
    assert out["partition_residual"] == pytest.approx(0.0)
    assert out["recursion_edges_added"] == 0
    assert out["causal_charge_created"] is False
    assert out["service_mass_duplicated"] is False


def test_endpoint_smooth_carrier_has_two_comparable_hard_shell_witnesses():
    cover = endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=8.0,
        endpoint_carrier_critical_mass=3.0,
    )
    assert cover["renewal_frequency"] == pytest.approx(6.0)
    assert cover["hard_shell_candidates"] == pytest.approx((6.0, 12.0))
    assert cover["next_corridor_renewal_candidates"] == pytest.approx((4.5, 9.0))
    assert cover["candidate_ratios_to_parent"] == pytest.approx((0.75, 1.5))
    assert cover["guaranteed_max_hard_shell_critical_mass_lower"] == pytest.approx(2.0)
    assert cover["new_causal_charge_created"] is False
    assert cover["new_physical_time_edge_created"] is False


def test_endpoint_hard_shell_exact_tie_stays_joint():
    cover = endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=8.0,
        endpoint_carrier_critical_mass=3.0,
    )
    out = realized_endpoint_hard_shell_witnesses(cover, (2.2, 2.2))
    assert out["joint_witness_frequencies"] == pytest.approx((6.0, 12.0))
    assert out["joint_next_corridor_renewal_frequencies"] == pytest.approx((4.5, 9.0))
    assert out["causal_primary_selected"] is False
    assert out["recursion_edges_added"] == 0
    assert out["physical_time_drop_added"] == 0.0


def test_endpoint_cover_rejects_insufficient_actual_shell_mass():
    cover = endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=8.0,
        endpoint_carrier_critical_mass=3.0,
    )
    with pytest.raises(ValueError, match="do not realize"):
        realized_endpoint_hard_shell_witnesses(cover, (0.2, 0.3))


def test_certificate_removes_service_theorem_depth_without_uv_overclaim():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "second event time or recursion edge" in cert["time_ontology"]
    assert "zero causal charge and zero recursion depth" in cert["material_ontology"]
    assert "3/4 or 3/2" in cert["scale_geometry"]
    assert "does not terminate UV-unbounded" in cert["scope"]
