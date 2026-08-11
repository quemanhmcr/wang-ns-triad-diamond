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
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION


def _segment(
    *,
    amplitude: float = 4.0,
    elapsed: tuple[float, ...] = (0.0, 1.0),
    strain: tuple[float, ...] = (0.0, 0.0),
    residual: tuple[float, ...] = (0.0, 0.0),
    hh: tuple[float, ...] = (0.0, 0.0),
):
    return same_carrier.SameCarrierMonitorSegment(
        carrier_id="event-Q",
        terminal_amplitude=amplitude,
        elapsed_times=elapsed,
        strain_action=strain,
        residual_impulse_abs=residual,
        hh_impulse_abs=hh,
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
    left = _segment(elapsed=(0.0, 1.0e-20))
    right = _segment(elapsed=(2.0e-20, 3.0e-20))

    with pytest.raises(TypeError, match="elapsed time|time token|discontinuity"):
        same_carrier.join_same_carrier_segments((left, right))


def test_tiny_terminal_coefficient_cannot_be_rebound_by_order_one_factor():
    amplitude = 1.0e-30
    left = _segment(amplitude=amplitude, elapsed=(0.0, 1.0))
    right = _segment(amplitude=2.0 * amplitude, elapsed=(1.0, 2.0))

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
    out = same_carrier.interior_checkpoint_accumulation_outcome(
        event_time=2.0,
        accumulation_time=0.7,
        terminal_amplitude=amplitude,
        strain_action_limit=0.5 * LOW_STRAIN_ACTION,
        residual_impulse_abs_limit=0.5 * RESIDUAL_FRACTION * amplitude,
        hh_impulse_abs_limit=0.5 * GENERATED_FRACTION * amplitude,
    )

    assert out["classification"] == "same_carrier_extends_across_interior_checkpoint_accumulation"
    assert out["joint_first_stops"] == ()


def test_above_threshold_limit_cannot_be_declared_first_hit_at_accumulation():
    with pytest.raises(ValueError, match="earlier|first|no-hit|continuity|overshoot"):
        same_carrier.interior_checkpoint_accumulation_outcome(
            event_time=2.0,
            accumulation_time=0.7,
            terminal_amplitude=4.0,
            strain_action_limit=1.1 * LOW_STRAIN_ACTION,
            residual_impulse_abs_limit=0.0,
            hh_impulse_abs_limit=0.0,
        )


def test_positive_remaining_native_time_is_not_rounded_to_t0():
    segment = _segment(elapsed=(0.0, 5.0e-21))
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
