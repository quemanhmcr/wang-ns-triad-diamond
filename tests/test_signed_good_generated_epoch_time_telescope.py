import pytest

from src.asynchronous_duhamel_sync import (
    BACKWARD_SUM_COEFF,
    INITIAL_HALF_SPAN,
    LIFETIME_GROWTH_MIN,
    MIN_REFERENCE_BACKSTEP,
)
from src.common_slice_coefficient_registration import HH_COEFFICIENT_OBSTRUCTION
from src.physical_energy_causal_bridge import heavy_half_physical_transfer, route_physical_energy_causality
from src.signed_good_generated_epoch_time_telescope import (
    ACTUAL_HH_GENERATION_BRANCH,
    STATUS,
    signed_good_generated_epoch_telescope,
    signed_good_step_from_energy_reentry,
    theorem_certificate,
)


def _reentry(lower: float = 1.0):
    return {
        "branch": ACTUAL_HH_GENERATION_BRANCH,
        "energy_gate": {
            "branch": ACTUAL_HH_GENERATION_BRANCH,
            "physical_hh_work_lower": lower,
        },
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
    }


def _step(child: float, parent: float, c: float, start: float, end: float, total: float = 2.0):
    return signed_good_step_from_energy_reentry(
        reentry=_reentry(0.8),
        selected_physical_half_slab={
            "start": start,
            "end": end,
            "mass": 1.1,
            "total": total,
            "normalized_parent_span_upper": float(INITIAL_HALF_SPAN),
        },
        child_frequency=child,
        parent_frequency=parent,
        scaled_lifetime=c,
    )


def test_raw_hh_coefficient_obstruction_is_not_a_generated_physical_event():
    with pytest.raises(TypeError, match="interval locator"):
        signed_good_step_from_energy_reentry(
            reentry={"branch": HH_COEFFICIENT_OBSTRUCTION, "requires_physical_energy_reentry": True},
            selected_physical_half_slab={
                "start": 0.0,
                "end": 0.0,
                "mass": 1.0,
                "total": 1.0,
                "normalized_parent_span_upper": float(INITIAL_HALF_SPAN),
            },
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_energy_reentry_must_supply_actual_positive_hh_work_lower():
    with pytest.raises(TypeError, match="HH-work lower"):
        signed_good_step_from_energy_reentry(
            reentry={
                "branch": ACTUAL_HH_GENERATION_BRANCH,
                "coefficient_impulse_used_as_physical_work": False,
                "observer_partition_motion_charged_as_physics": False,
            },
            selected_physical_half_slab={
                "start": 1.0,
                "end": 1.01,
                "mass": 1.0,
                "total": 1.5,
                "normalized_parent_span_upper": float(INITIAL_HALF_SPAN),
            },
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_generic_hh_outside_signed_good_scale_window_is_rejected():
    with pytest.raises(ValueError, match="5/8 upper"):
        _step(10.0, 7.0, 1.0, 1.0, 1.001)


def test_actual_physical_half_slab_binds_same_work_law_and_parent_span():
    child = 10.0
    parent = 6.1
    c = 1.0
    Tchild = c / child**2
    row = _step(child, parent, c, 2.0, 2.0 + 0.4 * Tchild)
    assert row.physical_hh_work_mass >= 0.5 * row.physical_hh_work_total
    assert row.physical_hh_work_total >= row.physical_hh_work_lower
    assert row.normalized_parent_span <= float(INITIAL_HALF_SPAN)
    assert row.parent_natural_lifetime / row.child_natural_lifetime > float(LIFETIME_GROWTH_MIN)


def test_consecutive_generated_common_surfaces_have_certified_backward_shift():
    c = 1.0
    r = 0.61
    child0 = 10.0
    parent0 = r * child0
    Tchild0 = c / child0**2
    row0 = _step(child0, parent0, c, 5.0, 5.0 + 0.25 * Tchild0)

    child1 = parent0
    parent1 = r * child1
    Tchild1 = c / child1**2
    end1 = row0.work_support_end
    width1 = 0.2 * Tchild1
    row1 = _step(child1, parent1, c, end1 - width1, end1)
    assert row1.work_support_start >= row0.common_reference_time

    out = signed_good_generated_epoch_telescope((row0, row1))
    actual = row0.common_reference_time - row1.common_reference_time
    required = float(MIN_REFERENCE_BACKSTEP) * row0.parent_natural_lifetime
    assert actual >= required
    assert out.cumulative_reference_backshift >= out.minimum_cumulative_backshift
    assert out.minimum_cumulative_backshift == pytest.approx(
        float(BACKWARD_SUM_COEFF)
        * row0.parent_natural_lifetime
        * (float(LIFETIME_GROWTH_MIN) ** 1 - 1.0)
    )


def test_nonconsecutive_support_restart_cannot_hide_inside_one_generated_epoch():
    c = 1.0
    r = 0.61
    row0 = _step(10.0, 6.1, c, 5.0, 5.001)
    row1 = _step(6.1, r * 6.1, c, row0.work_support_end + 1.0, row0.work_support_end + 1.0)
    with pytest.raises(ValueError, match="not contained"):
        signed_good_generated_epoch_telescope((row0, row1))


def test_parent_carrier_scale_must_continue_exactly_between_generated_layers():
    c = 1.0
    row0 = _step(10.0, 6.1, c, 5.0, 5.001)
    row1 = _step(6.0, 3.66, c, 4.99, 4.991)
    with pytest.raises(ValueError, match="actual signed-good parent carrier scale"):
        signed_good_generated_epoch_telescope((row0, row1))


def test_required_common_surface_reaching_t0_terminates_interior_generated_depth():
    c = 1.0
    # Choose an actual generated event whose parent natural window already reaches
    # the initial surface.  This event may exist, but no further interior generated
    # layer can be registered behind it.
    child = 1.0
    parent = 0.61
    row = _step(child, parent, c, 0.2, 0.2)
    assert row.common_reference_time < 0.0
    out = signed_good_generated_epoch_telescope((row,))
    assert out.hits_initial_boundary
    assert out.total_layer_upper_before_or_at_boundary == 1

    nxt = _step(parent, 0.61 * parent, c, 0.1, 0.1)
    with pytest.raises(ValueError, match="continued after.*t=0"):
        signed_good_generated_epoch_telescope((row, nxt))


def test_common_slices_are_not_recursive_events_or_event_count_currency():
    row = _step(10.0, 6.1, 1.0, 5.0, 5.001)
    out = signed_good_generated_epoch_telescope((row,))
    assert out.common_slices_are_recursive_events is False
    assert out.event_count_budget_used is False
    assert out.duhamel_weights_used_as_causal_law is False
    assert out.generic_hh_claimed is False


def test_certificate_closes_only_signed_good_generated_recurrence():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "raw |I_HH| coefficient obstruction" in cert["physical_input"]
    assert "3/5<N_parent/N_child<5/8" in cert["signed_good_geometry"]
    assert "1792/4875" in cert["asynchronous_backshift"]
    assert "1792/7605" in cert["time_telescope"]
    assert "registration surfaces, not event vertices" in cert["common_surface"]
    assert "actual positive child-energy HH work" in cert["causal_weights"]
    assert "generic HH/high-tail" in cert["master_consequence"]
    assert "mixed recurrence" in cert["scope"]
    assert "no Navier-Stokes global-regularity claim" in cert["scope"]


def test_canonical_adapter_accepts_the_real_physical_energy_and_heavy_half_outputs():
    child = 10.0
    parent = 6.1
    c = 1.0
    Tchild = c / child**2
    gate = route_physical_energy_causality(
        terminal_energy=1.0,
        initial_energy=0.1,
        residual_positive_work=0.1,
        strain_action=0.0,
    )
    assert gate["branch"] == ACTUAL_HH_GENERATION_BRANCH
    half = heavy_half_physical_transfer(
        times=[0.1 * Tchild, 0.3 * Tchild, 0.7 * Tchild, 0.9 * Tchild],
        positive_work_weights=[0.2, 0.2, 0.8, 0.4],
        slab_start=0.0,
        slab_end=Tchild,
    )
    row = signed_good_step_from_energy_reentry(
        reentry=gate,
        selected_physical_half_slab=half,
        child_frequency=child,
        parent_frequency=parent,
        scaled_lifetime=c,
    )
    assert row.physical_hh_work_total == pytest.approx(sum([0.2, 0.2, 0.8, 0.4]))
    assert row.physical_hh_work_total >= row.physical_hh_work_lower
    assert row.normalized_parent_span <= float(INITIAL_HALF_SPAN)
