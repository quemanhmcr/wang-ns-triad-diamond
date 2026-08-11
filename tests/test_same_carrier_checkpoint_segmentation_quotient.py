import cmath

import pytest

from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    RESIDUAL_FRACTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)
from src.full_natural_checkpoint_quotient import FullNaturalCheckpoint
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION
from src.same_carrier_checkpoint_segmentation_quotient import (
    SAME_CARRIER_CONTINUATION,
    STATUS,
    SameCarrierCheckpointPathCertificate,
    SameCarrierMonitorSegment,
    SameCarrierPrelimitCertificate,
    SameCarrierProvenance,
    SmoothPDEExtensionToken,
    checkpoint_continuation_policy,
    interior_checkpoint_accumulation_outcome,
    join_same_carrier_segments,
    maximal_same_carrier_outcome,
    partition_same_carrier_path,
    same_carrier_first_exit,
    segmentation_invariance,
    theorem_certificate,
)


def _provenance(amplitude: float = 4.0, terminal_time: float = 2.0) -> SameCarrierProvenance:
    return SameCarrierProvenance(
        event_id="physical-event",
        carrier_id="fixed-Q",
        terminal_dual_id="fixed-terminal-dual",
        trajectory_id="actual-NS-path",
        terminal_state_token="state-0",
        terminal_time=terminal_time,
        carrier_frequency=2.0,
        scaled_lifetime=1.0,
        terminal_coefficient=complex(amplitude),
    )


def _paths():
    amplitude = 4.0
    elapsed = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    strain = (0.0, 0.004, 0.009, 0.014, 0.020, 0.025)
    residual_radii = (0.0, 0.30, 0.55, 0.40, 0.62, 0.50)
    hh_radii = (0.0, 0.40, 0.72, 0.60, 0.85, 0.70)
    residual = tuple(radius * cmath.exp(1j * 2.0 * j) for j, radius in enumerate(residual_radii))
    hh = tuple(radius * cmath.exp(-1j * 1.7 * j) for j, radius in enumerate(hh_radii))
    tokens = tuple(f"state-{j}" for j in range(len(elapsed)))
    return _provenance(amplitude), tokens, elapsed, strain, residual, hh


def _partition(cuts=()):
    provenance, tokens, elapsed, strain, residual, hh = _paths()
    return partition_same_carrier_path(
        provenance=provenance,
        state_tokens=tokens,
        elapsed_times=elapsed,
        strain_action=strain,
        residual_impulse=residual,
        hh_impulse=hh,
        checkpoint_indices=cuts,
    )


def _checkpoint_path_certificate(
    provenance: SameCarrierProvenance | None = None,
) -> SameCarrierCheckpointPathCertificate:
    bound = provenance or _provenance(4.0)
    checkpoint = FullNaturalCheckpoint(
        terminal_time=bound.terminal_time,
        physical_time_drop=0.25,
        parent_shell_frequency=8.0 / 3.0,
        parent_shell_critical_mass_lower=2.0,
        corridor_frequency=2.0,
        scaled_lifetime=1.0,
        endpoint_carrier_critical_mass_lower=2.0,
        endpoint_shell_candidates=(2.0, 4.0),
    )
    segment = SameCarrierMonitorSegment(
        provenance=bound,
        state_tokens=(bound.terminal_state_token, "checkpoint-endpoint-state"),
        elapsed_times=(0.0, 0.25),
        strain_action=(0.0, 0.0),
        residual_impulse=(0.0j, 0.0j),
        hh_impulse=(0.0j, 0.0j),
    )
    return SameCarrierCheckpointPathCertificate(checkpoint, (segment,))


def _accumulation_segment(
    *,
    amplitude: float = 5.0,
    strain_end: float,
    residual_end: complex,
    hh_end: complex,
    accumulation_time: float = 0.7,
) -> SameCarrierMonitorSegment:
    provenance = _provenance(amplitude)
    elapsed_end = provenance.terminal_time - accumulation_time
    return SameCarrierMonitorSegment(
        provenance=provenance,
        state_tokens=(provenance.terminal_state_token, "accumulation-state"),
        elapsed_times=(0.0, elapsed_end),
        strain_action=(0.0, strain_end),
        residual_impulse=(0.0j, residual_end),
        hh_impulse=(0.0j, hh_end),
    )


def test_checkpoint_segmentation_does_not_move_same_carrier_first_hit():
    provenance, tokens, elapsed, strain, residual, hh = _paths()
    out = segmentation_invariance(
        provenance=provenance,
        state_tokens=tokens,
        elapsed_times=elapsed,
        strain_action=strain,
        residual_impulse=residual,
        hh_impulse=hh,
        checkpoint_indices=(2, 4),
    )
    assert out["first_elapsed"] is None
    assert out["joint_first_stops"] == ()
    assert out["first_time_residual"] == 0.0
    assert out["same_pde_trajectory_certified"] is True


def test_first_stop_is_global_cumulative_not_reset_per_natural_window():
    provenance = _provenance(4.0)
    elapsed = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)
    segments = partition_same_carrier_path(
        provenance=provenance,
        state_tokens=tuple(f"state-{j}" for j in range(6)),
        elapsed_times=elapsed,
        strain_action=(0.0, 0.004, 0.008, 0.012, 0.016, 0.020),
        residual_impulse=tuple(complex(x) for x in (0.0, 0.2, 0.5, 0.75, 1.0, 1.1)),
        hh_impulse=tuple(complex(x) for x in (0.0, 0.4, 0.6, 0.8, 0.7, 0.9)),
        checkpoint_indices=(2, 3),
    )
    out = same_carrier_first_exit(segments)
    assert out["first_elapsed"] == pytest.approx(0.8)
    assert out["joint_first_stops"] == (ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,)
    assert out["requires_physical_energy_reentry"] is True
    assert out["coefficient_impulses_used_as_work"] is False


def test_complex_impulse_magnitudes_need_not_be_monotone_or_additive():
    _, _, _, _, residual, hh = _paths()
    assert any(abs(residual[j + 1]) < abs(residual[j]) for j in range(len(residual) - 1))
    assert any(abs(hh[j + 1]) < abs(hh[j]) for j in range(len(hh) - 1))
    joined = join_same_carrier_segments(_partition((1, 3, 4)))
    assert joined["residual_impulse"] == residual
    assert joined["hh_impulse"] == hh
    assert joined["complex_phase_preserved"] is True


def test_checkpoint_baseline_reset_is_rejected():
    segments = list(_partition((2,)))
    segment = segments[1]
    segments[1] = SameCarrierMonitorSegment(
        provenance=segment.provenance,
        state_tokens=segment.state_tokens,
        elapsed_times=segment.elapsed_times,
        strain_action=(0.0,) + segment.strain_action[1:],
        residual_impulse=(0.0j,) + segment.residual_impulse[1:],
        hh_impulse=(0.0j,) + segment.hh_impulse[1:],
    )
    with pytest.raises(TypeError, match="reset/discontinuity"):
        join_same_carrier_segments(segments)


def test_checkpoint_cannot_reset_terminal_coefficient_or_replace_carrier():
    segments = list(_partition((2,)))
    segment = segments[1]
    foreign = SameCarrierProvenance(
        event_id=segment.provenance.event_id,
        carrier_id="new-Q",
        terminal_dual_id=segment.provenance.terminal_dual_id,
        trajectory_id=segment.provenance.trajectory_id,
        terminal_state_token=segment.provenance.terminal_state_token,
        terminal_time=segment.provenance.terminal_time,
        carrier_frequency=segment.provenance.carrier_frequency,
        scaled_lifetime=segment.provenance.scaled_lifetime,
        terminal_coefficient=2.0 * segment.provenance.terminal_coefficient,
    )
    segments[1] = SameCarrierMonitorSegment(
        provenance=foreign,
        state_tokens=segment.state_tokens,
        elapsed_times=segment.elapsed_times,
        strain_action=segment.strain_action,
        residual_impulse=segment.residual_impulse,
        hh_impulse=segment.hh_impulse,
    )
    with pytest.raises(TypeError, match="replace or rebind"):
        join_same_carrier_segments(segments)


def test_checkpoint_policy_keeps_hard_shell_readings_as_sidecars_only():
    provenance = _provenance(4.0)
    checkpoint = _checkpoint_path_certificate(provenance)
    out = checkpoint_continuation_policy(checkpoint, provenance=provenance)
    assert out["canonical_continuation"] == SAME_CARRIER_CONTINUATION
    assert out["hard_shell_checkpoint_witnesses"] == "state_sidecars_only"
    assert out["typed_checkpoint_and_carrier_provenance_verified"] is True
    assert out["actual_no_hit_pde_restriction_verified"] is True
    assert out["checkpoint_endpoint_state_token"] == "checkpoint-endpoint-state"
    for request, message in (
        ({"request_carrier_replacement": True}, "cannot replace"),
        ({"request_terminal_amplitude_reset": True}, "cannot reset the terminal"),
        ({"request_monitor_reset": True}, "cannot reset cumulative"),
    ):
        with pytest.raises(TypeError, match=message):
            checkpoint_continuation_policy(checkpoint, provenance=provenance, **request)


def test_checkpoint_path_certificate_rejects_an_earlier_physical_first_stop():
    provenance = _provenance(4.0)
    base = _checkpoint_path_certificate(provenance)
    hit = SameCarrierMonitorSegment(
        provenance=provenance,
        state_tokens=(provenance.terminal_state_token, "hit-before-checkpoint"),
        elapsed_times=(0.0, 0.25),
        strain_action=(0.0, 0.0),
        residual_impulse=(0.0j, RESIDUAL_FRACTION * provenance.terminal_amplitude),
        hh_impulse=(0.0j, 0.0j),
    )
    with pytest.raises(ValueError, match="earlier named first stop|no-hit"):
        SameCarrierCheckpointPathCertificate(base.checkpoint, (hit,))


def test_interior_checkpoint_accumulation_with_strict_margin_is_crossed_by_same_carrier():
    amplitude = 5.0
    segment = _accumulation_segment(
        amplitude=amplitude,
        strain_end=0.8 * LOW_STRAIN_ACTION,
        residual_end=0.8 * RESIDUAL_FRACTION * amplitude,
        hh_end=0.8 * GENERATED_FRACTION * amplitude,
    )
    out = interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.7,
        prelimit_certificate=SameCarrierPrelimitCertificate((segment,)),
        smooth_extension_token=SmoothPDEExtensionToken(
            trajectory_id=segment.provenance.trajectory_id,
            state_token=segment.state_tokens[-1],
            physical_time=0.7,
            open_interval=(0.6, 0.8),
        ),
    )
    assert out["classification"] == "same_carrier_extends_across_interior_checkpoint_accumulation"
    assert out["same_carrier_extends_past_accumulation"] is True


def test_closed_faces_at_accumulation_are_existing_first_stops():
    amplitude = 5.0
    segment = _accumulation_segment(
        amplitude=amplitude,
        strain_end=LOW_STRAIN_ACTION,
        residual_end=RESIDUAL_FRACTION * amplitude,
        hh_end=0.3 * GENERATED_FRACTION * amplitude,
    )
    out = interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.7,
        prelimit_certificate=SameCarrierPrelimitCertificate((segment,)),
        smooth_extension_token=None,
        tie_tolerance=1.0e-12,
    )
    assert set(out["joint_first_stops"]) == {
        "high_strain_critical_dissipation",
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
    }
    assert out["requires_physical_energy_reentry"] is True


def test_t0_absorbs_same_carrier_checkpoint_accumulation():
    provenance = _provenance(5.0)
    segment = SameCarrierMonitorSegment(
        provenance=provenance,
        state_tokens=(provenance.terminal_state_token, "initial-state"),
        elapsed_times=(0.0, 2.0),
        strain_action=(0.0, 0.5 * LOW_STRAIN_ACTION),
        residual_impulse=(0.0j, 0.0j),
        hh_impulse=(0.0j, 0.0j),
    )
    out = interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.0,
        prelimit_certificate=SameCarrierPrelimitCertificate((segment,)),
        smooth_extension_token=None,
    )
    assert out["classification"] == "absorbing_initial_boundary"
    assert out["joint_first_stops"] == ("t=0",)


def test_maximal_no_hit_same_carrier_reaches_t0_without_hardening_at_cuts():
    provenance = _provenance(4.0)
    elapsed = (0.0, 0.5, 1.0, 1.5, 2.0)
    segments = partition_same_carrier_path(
        provenance=provenance,
        state_tokens=tuple(f"state-{j}" for j in range(5)),
        elapsed_times=elapsed,
        strain_action=(0.0, 0.004, 0.008, 0.012, 0.016),
        residual_impulse=tuple(complex(x) for x in (0.0, 0.2, 0.4, 0.35, 0.5)),
        hh_impulse=tuple(complex(x) for x in (0.0, 0.3, 0.5, 0.45, 0.6)),
        checkpoint_indices=(1, 2, 3),
    )
    out = maximal_same_carrier_outcome(segments, event_time=2.0)
    assert out["classification"] == "absorbing_initial_boundary"
    assert out["joint_first_stops"] == ("t=0",)
    assert out["carrier_restarts"] == 0


def test_exact_joint_coefficient_faces_survive_checkpoint_segmentation():
    provenance = _provenance(4.0)
    elapsed = (0.0, 0.25, 0.5, 0.75, 1.0)
    out = segmentation_invariance(
        provenance=provenance,
        state_tokens=tuple(f"state-{j}" for j in range(5)),
        elapsed_times=elapsed,
        strain_action=(0.0, 0.003, 0.006, 0.009, 0.012),
        residual_impulse=tuple(complex(RESIDUAL_FRACTION * 4.0 * x / 0.75) for x in elapsed),
        hh_impulse=tuple(complex(GENERATED_FRACTION * 4.0 * x / 0.75) for x in elapsed),
        checkpoint_indices=(1, 2),
        tie_tolerance=1.0e-12,
    )
    assert out["first_elapsed"] == pytest.approx(0.75)
    assert set(out["joint_first_stops"]) == {
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
        HH_COEFFICIENT_OBSTRUCTION,
    }


def test_certificate_closes_only_certified_same_path_segmentation():
    certificate = theorem_certificate()
    assert certificate["status"] == STATUS
    assert "actual PDE trajectory id" in certificate["fixed_carrier"]
    assert "actual complex" in certificate["cumulative_monitors"]
    assert "exact gluing" in certificate["segmentation"]
    assert "arbitrary observer cuts are not full-natural service windows" in certificate["natural_windows"]
    assert "does not telescope infinitely recurring genuine physical owner events" in certificate["scope"]
