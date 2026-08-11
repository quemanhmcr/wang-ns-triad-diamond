import math

import pytest

from src.continuum_master_event_quotient import SupplierKind, validate_supplier_scale
from src.full_natural_checkpoint_quotient import (
    CERTIFIED_HIGH_TAIL_RATIO_LOWER,
    FULL_NATURAL_CHECKPOINT,
    STATUS,
    UPPER_COVER_RATIO,
    checkpoint_chain_ledger,
    checkpoint_from_full_natural_outcome,
    checkpoint_reregistration,
    geometric_uv_checkpoint_time,
    theorem_certificate,
)
from src.full_natural_service_corridor_quotient import (
    FULL_NATURAL_SERVICE_WITNESS,
)


def _full_outcome(M: float, c: float, mu: float = 2.0, t_factor: float = 4.0):
    A = 0.75 * M
    T = c / A**2
    t = t_factor * T
    return {
        "classification": FULL_NATURAL_SERVICE_WITNESS,
        "joint_first_stops": (),
        "required_elapsed": T,
        "observed_elapsed_end": T,
        "corridor_terminal_time": t,
        "corridor_endpoint_time": t - T,
        "physical_time_drop": T,
        "service_same_corridor_witness": True,
        "service_adds_recursion_depth": False,
        "uniform_square_service_lower": 0.2,
        "integrated_bounded_heat_service_lower": 0.2 * c,
        "endpoint_carrier_critical_mass_lower": mu,
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def test_checkpoint_keeps_parent_shell_and_actual_corridor_scale_distinct():
    M, c = 8.0, 1.2
    cp = checkpoint_from_full_natural_outcome(_full_outcome(M, c), parent_shell_frequency=M, scaled_lifetime=c)
    assert cp.parent_shell_frequency == pytest.approx(M)
    assert cp.corridor_frequency == pytest.approx(0.75 * M)
    assert cp.endpoint_shell_candidates == pytest.approx((0.75 * M, 1.5 * M))
    assert cp.candidate_ratios_to_parent == pytest.approx((0.75, 1.5))
    assert cp.physical_time_drop == pytest.approx(c / (0.75 * M) ** 2)


def test_endpoint_cover_reregistration_adds_no_event_or_causal_charge():
    M, c = 10.0, 0.8
    cp = checkpoint_from_full_natural_outcome(_full_outcome(M, c), parent_shell_frequency=M, scaled_lifetime=c)
    out = checkpoint_reregistration(cp, (1.4, 2.0))
    assert out["checkpoint_kind"] == FULL_NATURAL_CHECKPOINT
    assert out["joint_endpoint_witness_ratios"] == pytest.approx((1.5,))
    assert out["physical_event_created"] is False
    assert out["causal_charge_created"] is False
    assert out["recursion_edges_added"] == 0
    assert out["directional_scale_progress_supplied"] is False
    assert out["high_tail_supplier_admissible"] is False
    assert out["cover_ascent_interpreted_as_dynamics"] is False
    assert out["observer_selected_cover_branch"] is False
    assert out["joint_endpoint_witness_frequencies"] == pytest.approx((1.5 * M,))


def test_upper_three_halves_cover_witness_is_not_certified_high_tail_progress():
    assert UPPER_COVER_RATIO == pytest.approx(1.5)
    assert CERTIFIED_HIGH_TAIL_RATIO_LOWER == pytest.approx(2.0)
    assert UPPER_COVER_RATIO < CERTIFIED_HIGH_TAIL_RATIO_LOWER
    with pytest.raises(ValueError, match="lower scale ratio"):
        validate_supplier_scale(SupplierKind.HIGH_TAIL, 8.0, 12.0)


def test_checkpoint_chain_telescopes_real_time_with_zero_recursive_events():
    M1, c = 8.0, 1.0
    out1 = _full_outcome(M1, c, t_factor=6.0)
    cp1 = checkpoint_from_full_natural_outcome(out1, parent_shell_frequency=M1, scaled_lifetime=c)

    reread1 = checkpoint_reregistration(cp1, (0.8, 2.0))
    M2 = float(reread1["joint_endpoint_witness_frequencies"][0])
    assert M2 == pytest.approx(1.5 * M1)
    A2 = 0.75 * M2
    T2 = c / A2**2
    out2 = _full_outcome(M2, c)
    # Make the second corridor begin exactly where the first ended.
    out2["corridor_terminal_time"] = cp1.endpoint_time
    out2["corridor_endpoint_time"] = cp1.endpoint_time - T2
    out2["physical_time_drop"] = T2
    cp2 = checkpoint_from_full_natural_outcome(out2, parent_shell_frequency=M2, scaled_lifetime=c)

    led = checkpoint_chain_ledger((cp1, cp2))
    assert led["time_telescope_residual"] == pytest.approx(0.0, abs=1e-15)
    assert led["recursive_events_added"] == 0
    assert led["causal_charges_added"] == 0
    assert led["physical_event_vertices"] == 0


def test_uv_checkpoint_zano_time_is_preserved_as_event_free_obstruction():
    M, c = 2.0, 1.0
    first = c / (0.75 * M) ** 2
    total = geometric_uv_checkpoint_time(M, c, 1.5)
    assert total == pytest.approx(first / (1.0 - 1.0 / 1.5**2))
    assert total > first
    assert math.isfinite(total)


def test_checkpoint_adapter_rejects_first_stop_and_t0_like_endpoint():
    M, c = 4.0, 1.0
    bad = _full_outcome(M, c)
    bad["joint_first_stops"] = ("high_strain",)
    with pytest.raises(ValueError, match="physical first stop"):
        checkpoint_from_full_natural_outcome(bad, parent_shell_frequency=M, scaled_lifetime=c)

    root = _full_outcome(M, c)
    root["corridor_endpoint_time"] = 0.0
    root["corridor_terminal_time"] = root["physical_time_drop"]
    with pytest.raises(ValueError, match="absorbing boundary"):
        checkpoint_from_full_natural_outcome(root, parent_shell_frequency=M, scaled_lifetime=c)



def test_checkpoint_rejects_any_use_of_coefficient_impulse_as_physical_work():
    M, c = 4.0, 1.0
    bad = _full_outcome(M, c)
    bad["coefficient_impulses_used_as_work"] = True
    with pytest.raises(ValueError, match="coefficient impulse"):
        checkpoint_from_full_natural_outcome(bad, parent_shell_frequency=M, scaled_lifetime=c)


def test_observer_cannot_choose_a_cover_frequency_instead_of_actual_shell_masses():
    M, c = 8.0, 1.0
    cp = checkpoint_from_full_natural_outcome(_full_outcome(M, c), parent_shell_frequency=M, scaled_lifetime=c)
    with pytest.raises(ValueError, match="two finite nonnegative actual endpoint hard-shell"):
        checkpoint_reregistration(cp, (1.5 * M,))


def test_actual_lower_shell_winner_cannot_be_replaced_by_desired_upper_cover_branch():
    M, c = 8.0, 1.0
    cp = checkpoint_from_full_natural_outcome(_full_outcome(M, c), parent_shell_frequency=M, scaled_lifetime=c)
    out = checkpoint_reregistration(cp, (2.0, 0.8))
    assert out["joint_endpoint_witness_frequencies"] == pytest.approx((0.75 * M,))
    assert out["joint_endpoint_witness_ratios"] == pytest.approx((0.75,))
    assert out["observer_selected_cover_branch"] is False


def test_exact_endpoint_shell_mass_tie_remains_joint_at_checkpoint():
    M, c = 8.0, 1.0
    cp = checkpoint_from_full_natural_outcome(_full_outcome(M, c), parent_shell_frequency=M, scaled_lifetime=c)
    out = checkpoint_reregistration(cp, (2.0, 2.0))
    assert out["joint_endpoint_witness_frequencies"] == pytest.approx((0.75 * M, 1.5 * M))
    assert out["joint_endpoint_witness_ratios"] == pytest.approx((0.75, 1.5))


def test_certificate_separates_event_recursion_from_uv_checkpoint_continuation():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "zero physical event vertices" in cert["time_semantics"]
    assert "3/4 and 3/2" in cert["cover_geometry"]
    assert "below" in cert["cover_geometry"]
    assert "event-free PDE continuation seam" in cert["remaining_uv"]
    assert "does not prove" in cert["scope"]
