"""Adversarial audit of the same-carrier checkpoint quotient.

The quotient is legitimate only for restrictions of one event-anchored PDE path.
Observer-scale tolerances, matching impulse magnitudes, or a repeated string label
cannot certify a common carrier, terminal dual, trajectory, or native time token.
"""

import math

import pytest

import src.same_carrier_checkpoint_segmentation_quotient as same_carrier
from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    RESIDUAL_FRACTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)
from src.full_natural_checkpoint_quotient import FullNaturalCheckpoint
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION


def _segment(
    *,
    amplitude: float = 4.0,
    terminal_time: float = 2.0,
    elapsed: tuple[float, ...] = (0.0, 1.0),
    strain: tuple[float, ...] = (0.0, 0.0),
    residual: tuple[float, ...] = (0.0, 0.0),
    hh: tuple[float, ...] = (0.0, 0.0),
    state_tokens: tuple[str, ...] | None = None,
):
    provenance = same_carrier.SameCarrierProvenance(
        event_id="event",
        carrier_id="event-Q",
        terminal_dual_id="terminal-dual",
        trajectory_id="actual-NS-trajectory",
        terminal_state_token="terminal-state",
        terminal_time=terminal_time,
        carrier_frequency=2.0,
        scaled_lifetime=1.0,
        terminal_coefficient=complex(amplitude),
    )
    tokens = state_tokens or tuple(
        "terminal-state" if j == 0 and elapsed[0] == 0.0 else f"state-{j}"
        for j in range(len(elapsed))
    )
    return same_carrier.SameCarrierMonitorSegment(
        provenance=provenance,
        state_tokens=tokens,
        elapsed_times=elapsed,
        strain_action=strain,
        residual_impulse=tuple(complex(x) for x in residual),
        hh_impulse=tuple(complex(x) for x in hh),
    )


def test_positive_native_elapsed_origin_is_not_terminal_event_zero():
    segment = _segment(elapsed=(1.0e-30, 2.0e-30))

    with pytest.raises(ValueError, match="zero|terminal event|origin"):
        same_carrier.join_same_carrier_segments((segment,))


def test_nonzero_cumulative_impulse_at_event_is_not_erased_by_unit_floor():
    amplitude = 1.0e-30
    segment = _segment(
        amplitude=amplitude,
        elapsed=(0.0, 1.0e-30),
        residual=(0.5 * RESIDUAL_FRACTION * amplitude, 0.5 * RESIDUAL_FRACTION * amplitude),
    )

    with pytest.raises(ValueError, match="start at zero|terminal event"):
        same_carrier.join_same_carrier_segments((segment,))


def test_native_time_gap_cannot_be_hidden_when_checkpoint_rows_are_joined():
    left = _segment(
        elapsed=(0.0, 1.0e-20),
        state_tokens=("terminal-state", "shared-boundary"),
    )
    right = _segment(
        elapsed=(2.0e-20, 3.0e-20),
        state_tokens=("shared-boundary", "later-state"),
    )

    with pytest.raises(TypeError, match="elapsed time|time token|discontinuity"):
        same_carrier.join_same_carrier_segments((left, right))


def test_tiny_terminal_coefficient_cannot_be_rebound_by_order_one_factor():
    amplitude = 1.0e-30
    left = _segment(
        amplitude=amplitude,
        elapsed=(0.0, 1.0),
        state_tokens=("terminal-state", "shared-boundary"),
    )
    right = _segment(
        amplitude=2.0 * amplitude,
        elapsed=(1.0, 2.0),
        state_tokens=("shared-boundary", "later-state"),
    )

    with pytest.raises(TypeError, match="terminal coefficient|amplitude|reset"):
        same_carrier.join_same_carrier_segments((left, right))


def test_distinct_native_debuts_are_not_fused_into_an_exact_tie():
    amplitude = 4.0
    times = (0.0, 1.0e-16, 2.0e-16, 3.0e-16)
    segment = _segment(
        amplitude=amplitude,
        elapsed=times,
        strain=(0.0,) * 4,
        residual=(0.0, RESIDUAL_FRACTION * amplitude, RESIDUAL_FRACTION * amplitude, RESIDUAL_FRACTION * amplitude),
        hh=(0.0, 0.0, GENERATED_FRACTION * amplitude, GENERATED_FRACTION * amplitude),
    )

    out = same_carrier.same_carrier_first_exit((segment,))
    assert out["first_elapsed"] == times[1]
    assert out["joint_first_stops"] == (ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,)


def test_tiny_coefficient_scale_retains_strict_accumulation_margin():
    amplitude = 1.0e-30
    segment = _segment(
        amplitude=amplitude,
        elapsed=(0.0, 1.3),
        strain=(0.0, 0.5 * LOW_STRAIN_ACTION),
        residual=(0.0, 0.5 * RESIDUAL_FRACTION * amplitude),
        hh=(0.0, 0.5 * GENERATED_FRACTION * amplitude),
        state_tokens=("terminal-state", "accumulation-state"),
    )
    out = same_carrier.interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.7,
        prelimit_certificate=same_carrier.SameCarrierPrelimitCertificate((segment,)),
        smooth_extension_token=same_carrier.SmoothPDEExtensionToken(
            trajectory_id="actual-NS-trajectory",
            state_token="accumulation-state",
            physical_time=0.7,
            open_interval=(0.6, 0.8),
        ),
    )

    assert out["classification"] == "same_carrier_extends_across_interior_checkpoint_accumulation"
    assert out["joint_first_stops"] == ()


def test_above_threshold_limit_cannot_be_declared_first_hit_at_accumulation():
    segment = _segment(
        elapsed=(0.0, 1.3),
        strain=(0.0, 1.1 * LOW_STRAIN_ACTION),
        state_tokens=("terminal-state", "accumulation-state"),
    )
    with pytest.raises(ValueError, match="earlier|first|no-hit|continuity|overshoot"):
        same_carrier.interior_checkpoint_accumulation_outcome(
            event_time=2.0,
            accumulation_time=0.7,
            prelimit_certificate=same_carrier.SameCarrierPrelimitCertificate((segment,)),
            smooth_extension_token=None,
        )


def test_positive_remaining_native_time_is_not_rounded_to_t0():
    segment = _segment(
        terminal_time=1.0e-20,
        elapsed=(0.0, 5.0e-21),
        state_tokens=("terminal-state", "interior-state"),
    )
    out = same_carrier.maximal_same_carrier_outcome((segment,), event_time=1.0e-20)

    assert out["classification"] == "same_carrier_event_free_continuation"
    assert out["endpoint_time"] > 0.0
    assert out["remaining_backward_time"] == 5.0e-21


def test_plain_dictionary_cannot_forge_checkpoint_continuation_authority():
    forged = {
        "checkpoint_kind": "full_natural_analysis_checkpoint",
        "physical_event_created": False,
        "causal_charge_created": False,
        "recursion_edges_added": 0,
    }

    with pytest.raises(TypeError, match="FullNaturalCheckpoint|typed|checkpoint record"):
        same_carrier.checkpoint_continuation_policy(forged)


def test_typed_checkpoint_matching_only_time_scale_and_lifetime_cannot_claim_same_pde_path():
    checkpoint = FullNaturalCheckpoint(
        terminal_time=2.0,
        physical_time_drop=0.25,
        parent_shell_frequency=8.0 / 3.0,
        parent_shell_critical_mass_lower=2.0,
        corridor_frequency=2.0,
        scaled_lifetime=1.0,
        endpoint_carrier_critical_mass_lower=2.0,
        endpoint_shell_candidates=(2.0, 4.0),
    )
    foreign = same_carrier.SameCarrierProvenance(
        event_id="foreign-event",
        carrier_id="foreign-Q",
        terminal_dual_id="foreign-dual",
        trajectory_id="foreign-NS-trajectory",
        terminal_state_token="foreign-terminal-state",
        terminal_time=2.0,
        carrier_frequency=2.0,
        scaled_lifetime=1.0,
        terminal_coefficient=4.0 + 0.0j,
    )

    with pytest.raises(TypeError, match="path|restriction|bound|provenance"):
        same_carrier.checkpoint_continuation_policy(checkpoint, provenance=foreign)


def test_checkpoint_path_certificate_rejects_a_foreign_expected_trajectory():
    checkpoint = FullNaturalCheckpoint(
        terminal_time=2.0,
        physical_time_drop=0.25,
        parent_shell_frequency=8.0 / 3.0,
        parent_shell_critical_mass_lower=2.0,
        corridor_frequency=2.0,
        scaled_lifetime=1.0,
        endpoint_carrier_critical_mass_lower=2.0,
        endpoint_shell_candidates=(2.0, 4.0),
    )
    actual_segment = _segment(
        elapsed=(0.0, 0.25),
        state_tokens=("terminal-state", "actual-checkpoint-state"),
    )
    certificate = same_carrier.SameCarrierCheckpointPathCertificate(
        checkpoint,
        (actual_segment,),
    )
    actual = actual_segment.provenance
    foreign = same_carrier.SameCarrierProvenance(
        event_id=actual.event_id,
        carrier_id=actual.carrier_id,
        terminal_dual_id=actual.terminal_dual_id,
        trajectory_id="foreign-NS-trajectory",
        terminal_state_token=actual.terminal_state_token,
        terminal_time=actual.terminal_time,
        carrier_frequency=actual.carrier_frequency,
        scaled_lifetime=actual.scaled_lifetime,
        terminal_coefficient=actual.terminal_coefficient,
    )

    with pytest.raises(TypeError, match="different.*trajectory|event/carrier/dual/PDE"):
        same_carrier.checkpoint_continuation_policy(certificate, provenance=foreign)


def test_segment_carries_complex_path_and_physical_provenance_not_only_magnitudes():
    fields = set(same_carrier.SameCarrierMonitorSegment.__dataclass_fields__)
    assert {"provenance", "state_tokens", "residual_impulse", "hh_impulse"} <= fields


def test_piecewise_linear_complex_phase_sets_the_modulus_debut():
    debut = getattr(same_carrier, "complex_modulus_debut_piecewise_linear", None)
    assert debut is not None, "the audit requires the complex impulse path, not linear interpolation of magnitudes"
    actual = debut((0.0, 1.0), (0.5 + 0.0j, -2.0 + 0.0j), 1.0)
    assert actual == pytest.approx(0.6, rel=0.0, abs=8.0 * math.ulp(0.6))


def test_fixed_carrier_natural_windows_have_positive_native_duration_and_no_interior_zeno():
    capacity = getattr(same_carrier, "fixed_carrier_natural_window_capacity", None)
    assert capacity is not None, "natural service windows must be typed separately from arbitrary observer cuts"
    out = capacity(event_time=1.0, carrier_frequency=2.0, scaled_lifetime=1.0)
    assert out["native_window_duration"] == 0.25
    assert out["maximum_complete_windows_before_t0"] == 4
    assert out["interior_zeno_possible"] is False


def test_accumulation_classifier_requires_actual_no_hit_path_and_smooth_extension_tokens():
    signature = __import__("inspect").signature(
        same_carrier.interior_checkpoint_accumulation_outcome
    )
    assert "prelimit_certificate" in signature.parameters
    assert "smooth_extension_token" in signature.parameters


def test_certificate_does_not_attach_diagnostic_uv_scales_to_fictitious_natural_durations():
    cert = same_carrier.theorem_certificate()
    text = " ".join(str(value) for value in cert.values())
    assert "arbitrary observer cuts are not full-natural service windows" in text
    assert "fixed positive native duration" in text
    assert "does not prove Navier-Stokes global regularity" in text
