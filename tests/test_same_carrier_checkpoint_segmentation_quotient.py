import pytest

from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    RESIDUAL_FRACTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)
from src.full_natural_checkpoint_quotient import FULL_NATURAL_CHECKPOINT
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION
from src.same_carrier_checkpoint_segmentation_quotient import (
    SAME_CARRIER_CONTINUATION,
    STATUS,
    SameCarrierMonitorSegment,
    checkpoint_continuation_policy,
    interior_checkpoint_accumulation_outcome,
    join_same_carrier_segments,
    maximal_same_carrier_outcome,
    partition_same_carrier_path,
    same_carrier_first_exit,
    segmentation_invariance,
    theorem_certificate,
)


def _paths():
    amp = 4.0
    t = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    K = (0.0, 0.004, 0.009, 0.014, 0.020, 0.025)
    # Magnitudes are cumulative from the same terminal event but may decrease by complex phase cancellation.
    IR = (0.0, 0.30, 0.55, 0.40, 0.62, 0.50)
    IH = (0.0, 0.40, 0.72, 0.60, 0.85, 0.70)
    return amp, t, K, IR, IH


def test_checkpoint_segmentation_does_not_move_same_carrier_first_hit():
    amp, t, K, IR, IH = _paths()
    out = segmentation_invariance(
        carrier_id="fixed-Q",
        terminal_amplitude=amp,
        elapsed_times=t,
        strain_action=K,
        residual_impulse_abs=IR,
        hh_impulse_abs=IH,
        checkpoint_indices=(2, 4),
    )
    assert out["first_elapsed"] is None
    assert out["joint_first_stops"] == ()
    assert out["first_time_residual"] == pytest.approx(0.0)
    assert out["carrier_restarts"] == 0
    assert out["monitor_resets"] == 0


def test_first_stop_is_global_cumulative_not_reset_per_natural_window():
    amp = 4.0
    t = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    K = (0.0, 0.004, 0.008, 0.012, 0.016, 0.020)
    IR = (0.0, 0.2, 0.5, 0.75, 1.0, 1.1)  # threshold amp/4 = 1 at elapsed .8
    IH = (0.0, 0.4, 0.6, 0.8, 0.7, 0.9)
    segs = partition_same_carrier_path(
        carrier_id="fixed-Q",
        terminal_amplitude=amp,
        elapsed_times=t,
        strain_action=K,
        residual_impulse_abs=IR,
        hh_impulse_abs=IH,
        checkpoint_indices=(2, 3),
    )
    out = same_carrier_first_exit(segs)
    assert out["first_elapsed"] == pytest.approx(0.8)
    assert out["joint_first_stops"] == (ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,)
    assert out["requires_physical_energy_reentry"] is True
    assert out["coefficient_impulses_used_as_work"] is False


def test_coefficient_impulse_magnitudes_need_not_be_monotone_or_additive():
    amp, t, K, IR, IH = _paths()
    assert any(IR[j + 1] < IR[j] for j in range(len(IR) - 1))
    assert any(IH[j + 1] < IH[j] for j in range(len(IH) - 1))
    segs = partition_same_carrier_path(
        carrier_id="fixed-Q",
        terminal_amplitude=amp,
        elapsed_times=t,
        strain_action=K,
        residual_impulse_abs=IR,
        hh_impulse_abs=IH,
        checkpoint_indices=(1, 3, 4),
    )
    joined = join_same_carrier_segments(segs)
    assert joined["residual_impulse_abs"] == pytest.approx(IR)
    assert joined["hh_impulse_abs"] == pytest.approx(IH)


def test_checkpoint_baseline_reset_is_rejected():
    amp, t, K, IR, IH = _paths()
    segs = list(
        partition_same_carrier_path(
            carrier_id="fixed-Q",
            terminal_amplitude=amp,
            elapsed_times=t,
            strain_action=K,
            residual_impulse_abs=IR,
            hh_impulse_abs=IH,
            checkpoint_indices=(2,),
        )
    )
    b = segs[1]
    segs[1] = SameCarrierMonitorSegment(
        carrier_id=b.carrier_id,
        terminal_amplitude=b.terminal_amplitude,
        elapsed_times=b.elapsed_times,
        strain_action=(0.0,) + b.strain_action[1:],
        residual_impulse_abs=(0.0,) + b.residual_impulse_abs[1:],
        hh_impulse_abs=(0.0,) + b.hh_impulse_abs[1:],
    )
    with pytest.raises(TypeError, match="reset/discontinuity"):
        join_same_carrier_segments(segs)


def test_checkpoint_cannot_reset_terminal_amplitude_or_replace_carrier():
    amp, t, K, IR, IH = _paths()
    segs = list(
        partition_same_carrier_path(
            carrier_id="fixed-Q",
            terminal_amplitude=amp,
            elapsed_times=t,
            strain_action=K,
            residual_impulse_abs=IR,
            hh_impulse_abs=IH,
            checkpoint_indices=(2,),
        )
    )
    b = segs[1]
    segs[1] = SameCarrierMonitorSegment(
        carrier_id="new-Q",
        terminal_amplitude=b.terminal_amplitude,
        elapsed_times=b.elapsed_times,
        strain_action=b.strain_action,
        residual_impulse_abs=b.residual_impulse_abs,
        hh_impulse_abs=b.hh_impulse_abs,
    )
    with pytest.raises(TypeError, match="cannot replace"):
        join_same_carrier_segments(segs)


def test_checkpoint_policy_keeps_hard_shell_readings_as_sidecars_only():
    checkpoint = {
        "checkpoint_kind": FULL_NATURAL_CHECKPOINT,
        "physical_event_created": False,
        "causal_charge_created": False,
        "recursion_edges_added": 0,
    }
    out = checkpoint_continuation_policy(checkpoint)
    assert out["canonical_continuation"] == SAME_CARRIER_CONTINUATION
    assert out["hard_shell_checkpoint_witnesses"] == "state_sidecars_only"
    assert out["carrier_replacement_authorized"] is False
    assert out["terminal_amplitude_reset_authorized"] is False
    assert out["monitor_reset_authorized"] is False
    assert out["checkpoint_scale_path_is_physical_lineage"] is False
    with pytest.raises(TypeError, match="cannot replace"):
        checkpoint_continuation_policy(checkpoint, request_carrier_replacement=True)
    with pytest.raises(TypeError, match="cannot reset the terminal"):
        checkpoint_continuation_policy(checkpoint, request_terminal_amplitude_reset=True)
    with pytest.raises(TypeError, match="cannot reset cumulative"):
        checkpoint_continuation_policy(checkpoint, request_monitor_reset=True)


def test_interior_checkpoint_accumulation_with_strict_margin_is_crossed_by_same_carrier():
    amp = 5.0
    out = interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.7,
        terminal_amplitude=amp,
        strain_action_limit=0.8 * LOW_STRAIN_ACTION,
        residual_impulse_abs_limit=0.8 * RESIDUAL_FRACTION * amp,
        hh_impulse_abs_limit=0.8 * GENERATED_FRACTION * amp,
    )
    assert out["classification"] == "same_carrier_extends_across_interior_checkpoint_accumulation"
    assert out["same_carrier_extends_past_accumulation"] is True
    assert out["checkpoint_accumulation_is_obstruction"] is False


def test_closed_face_at_interior_checkpoint_accumulation_is_the_existing_first_stop():
    amp = 5.0
    out = interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.7,
        terminal_amplitude=amp,
        strain_action_limit=LOW_STRAIN_ACTION,
        residual_impulse_abs_limit=RESIDUAL_FRACTION * amp,
        hh_impulse_abs_limit=0.3 * GENERATED_FRACTION * amp,
    )
    assert set(out["joint_first_stops"]) == {
        "high_strain_critical_dissipation",
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
    }
    assert out["requires_physical_energy_reentry"] is True
    assert out["coefficient_impulses_used_as_work"] is False


def test_t0_absorbs_same_carrier_checkpoint_accumulation():
    out = interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.0,
        terminal_amplitude=5.0,
        strain_action_limit=0.0,
        residual_impulse_abs_limit=0.0,
        hh_impulse_abs_limit=0.0,
    )
    assert out["classification"] == "absorbing_initial_boundary"
    assert out["joint_first_stops"] == ("t=0",)


def test_maximal_no_hit_same_carrier_reaches_t0_without_hardening_at_checkpoints():
    amp = 4.0
    t = (0.0, 0.5, 1.0, 1.5, 2.0)
    K = (0.0, 0.004, 0.008, 0.012, 0.016)
    IR = (0.0, 0.2, 0.4, 0.35, 0.5)
    IH = (0.0, 0.3, 0.5, 0.45, 0.6)
    segs = partition_same_carrier_path(
        carrier_id="fixed-Q",
        terminal_amplitude=amp,
        elapsed_times=t,
        strain_action=K,
        residual_impulse_abs=IR,
        hh_impulse_abs=IH,
        checkpoint_indices=(1, 2, 3),
    )
    out = maximal_same_carrier_outcome(segs, event_time=2.0)
    assert out["classification"] == "absorbing_initial_boundary"
    assert out["joint_first_stops"] == ("t=0",)
    assert out["carrier_restarts"] == 0
    assert out["monitor_resets"] == 0


def test_exact_joint_coefficient_faces_survive_checkpoint_segmentation():
    amp = 4.0
    t = (0.0, 0.25, 0.5, 0.75, 1.0)
    K = (0.0, 0.003, 0.006, 0.009, 0.012)
    IR = tuple(RESIDUAL_FRACTION * amp * x / 0.75 for x in t)
    IH = tuple(GENERATED_FRACTION * amp * x / 0.75 for x in t)
    out = segmentation_invariance(
        carrier_id="fixed-Q",
        terminal_amplitude=amp,
        elapsed_times=t,
        strain_action=K,
        residual_impulse_abs=IR,
        hh_impulse_abs=IH,
        checkpoint_indices=(1, 2),
        tie_tolerance=1e-12,
    )
    assert out["first_elapsed"] == pytest.approx(0.75)
    assert set(out["joint_first_stops"]) == {
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
        HH_COEFFICIENT_OBSTRUCTION,
    }


def test_certificate_closes_checkpoint_segmentation_not_genuine_owner_recurrence():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "same smooth event-anchored carrier" in cert["fixed_carrier"]
    assert "never summed segment-by-segment" in cert["cumulative_monitors"]
    assert "leaves the first physical stop time" in cert["segmentation"]
    assert "cannot reset terminal amplitude" in cert["checkpoint_policy"]
    assert "stop or is crossed" in cert["interior_accumulation"]
    assert "does not telescope infinitely recurring genuine physical owner events" in cert["scope"]
