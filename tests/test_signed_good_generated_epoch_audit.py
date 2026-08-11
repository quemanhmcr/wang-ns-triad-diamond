"""Adversarial boundary tests for the signed-good generated-HH epoch.

The continuum telescope is homogeneous and follows one actual positive HH-work
law on one PDE history.  These tests make sure the executable cannot replace
that native law by absolute observer units or splice unrelated certificates.
"""

import math

import pytest

from src.asynchronous_duhamel_sync import INITIAL_HALF_SPAN
from src.signed_good_generated_epoch_time_telescope import (
    ACTUAL_HH_GENERATION_BRANCH,
    SignedGoodGeneratedWorkProvenance,
    signed_good_generated_epoch_telescope,
    signed_good_step_from_energy_reentry,
)


def _provenance(
    tag: str = "a",
    *,
    child: float = 10.0,
    parent: float = 6.1,
    c: float = 1.0,
    slab_start: float = 1.0,
    slab_end: float = 1.01,
) -> SignedGoodGeneratedWorkProvenance:
    return SignedGoodGeneratedWorkProvenance(
        event_id=f"event-{tag}",
        trajectory_id="actual-NS-history",
        child_carrier_id=f"carrier-{child:.17g}",
        generated_parent_carrier_id=f"carrier-{parent:.17g}",
        work_law_id=f"positive-HH-work-{tag}",
        child_frequency=child,
        parent_frequency=parent,
        scaled_lifetime=c,
        slab_start=slab_start,
        slab_end=slab_end,
    )


def _reentry(
    lower: float,
    *,
    provenance: object | None,
    nested_coefficient_flag: bool = False,
) -> dict[str, object]:
    out: dict[str, object] = {
        "branch": ACTUAL_HH_GENERATION_BRANCH,
        "energy_gate": {
            "branch": ACTUAL_HH_GENERATION_BRANCH,
            "physical_hh_work_lower": lower,
            "coefficient_impulse_used_as_physical_work": nested_coefficient_flag,
            "provenance": provenance,
        },
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
    }
    if provenance is not None:
        out["provenance"] = provenance
    return out


def _half(
    *,
    start: float,
    end: float,
    mass: float,
    total: float,
    provenance: object | None,
    span_upper: float = float(INITIAL_HALF_SPAN),
) -> dict[str, object]:
    out: dict[str, object] = {
        "start": start,
        "end": end,
        "mass": mass,
        "total": total,
        "normalized_parent_span_upper": span_upper,
    }
    if provenance is not None:
        out["provenance"] = provenance
    return out


def _step(
    child: float,
    parent: float,
    c: float,
    start: float,
    end: float,
    *,
    provenance: object | None = None,
):
    if provenance is None:
        Tchild = math.exp(math.log(c) - 2.0 * math.log(child))
        slab_end = end if end > 0.0 else Tchild
        slab_start = max(0.0, slab_end - Tchild)
        token = _provenance(
            f"{child:.17g}-{start:.17g}",
            child=child,
            parent=parent,
            c=c,
            slab_start=slab_start,
            slab_end=slab_end,
        )
    else:
        token = provenance
    return signed_good_step_from_energy_reentry(
        reentry=_reentry(0.8, provenance=token),
        selected_physical_half_slab=_half(
            start=start,
            end=end,
            mass=1.1,
            total=2.0,
            provenance=token,
        ),
        child_frequency=child,
        parent_frequency=parent,
        scaled_lifetime=c,
    )


def test_untyped_scalar_dictionaries_cannot_claim_actual_pde_work_provenance():
    with pytest.raises(TypeError, match="provenance"):
        signed_good_step_from_energy_reentry(
            reentry=_reentry(0.8, provenance=None),
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.001,
                mass=1.1,
                total=2.0,
                provenance=None,
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_energy_gate_and_half_slab_cannot_splice_foreign_work_laws():
    with pytest.raises(TypeError, match="provenance|work law"):
        signed_good_step_from_energy_reentry(
            reentry=_reentry(0.8, provenance=_provenance("gate")),
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.001,
                mass=1.1,
                total=2.0,
                provenance=_provenance("half"),
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_nested_energy_gate_cannot_splice_a_foreign_typed_work_law():
    gate_token = _provenance("gate")
    foreign_token = _provenance("nested")
    reentry = _reentry(0.8, provenance=gate_token)
    reentry["energy_gate"]["provenance"] = foreign_token
    with pytest.raises(TypeError, match="provenance|foreign physical work laws"):
        signed_good_step_from_energy_reentry(
            reentry=reentry,
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.001,
                mass=1.1,
                total=2.0,
                provenance=gate_token,
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_nested_energy_gate_cannot_hide_forbidden_coefficient_work():
    token = _provenance()
    with pytest.raises(TypeError, match="coefficient impulse"):
        signed_good_step_from_energy_reentry(
            reentry=_reentry(
                0.8,
                provenance=token,
                nested_coefficient_flag=True,
            ),
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.001,
                mass=1.1,
                total=2.0,
                provenance=token,
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_tiny_native_work_cannot_fake_heavy_half_or_energy_lower_bound():
    token = _provenance()
    with pytest.raises(ValueError, match="heavy half|lower bound|work"):
        signed_good_step_from_energy_reentry(
            reentry=_reentry(2.0e-120, provenance=token),
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.0,
                mass=1.0e-130,
                total=1.0e-120,
                provenance=token,
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_selected_positive_sublaw_mass_cannot_exceed_its_total():
    token = _provenance()
    with pytest.raises(ValueError, match="mass|total"):
        signed_good_step_from_energy_reentry(
            reentry=_reentry(0.8, provenance=token),
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.0,
                mass=3.0,
                total=2.0,
                provenance=token,
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_nonfinite_parent_span_certificate_fails_closed():
    token = _provenance()
    with pytest.raises((TypeError, ValueError), match="span|finite"):
        signed_good_step_from_energy_reentry(
            reentry=_reentry(0.8, provenance=token),
            selected_physical_half_slab=_half(
                start=1.0,
                end=1.0,
                mass=1.1,
                total=2.0,
                provenance=token,
                span_upper=math.nan,
            ),
            child_frequency=10.0,
            parent_frequency=6.1,
            scaled_lifetime=1.0,
        )


def test_typed_work_provenance_cannot_exceed_the_child_natural_slab():
    with pytest.raises(ValueError, match="child natural slab"):
        _provenance("long-slab", slab_start=1.0, slab_end=1.02)


def test_tiny_native_frequency_cannot_accept_a_foreign_next_carrier():
    c = 1.0e-240
    row0 = _step(1.0e-120, 0.61e-120, c, 10.0, 10.0)
    foreign_child = 1.01 * row0.parent_frequency
    row1 = _step(
        foreign_child,
        0.61 * foreign_child,
        c,
        10.0,
        10.0,
    )
    with pytest.raises(ValueError, match="foreign physical parent carrier|actual signed-good parent carrier scale"):
        signed_good_generated_epoch_telescope((row0, row1))


def test_tiny_scaled_lifetime_cannot_change_inside_one_epoch():
    c = 1.0e-240
    row0 = _step(1.0e-120, 0.61e-120, c, 10.0, 10.0)
    row1 = _step(
        row0.parent_frequency,
        0.61 * row0.parent_frequency,
        1.01 * c,
        10.0,
        10.0,
    )
    with pytest.raises(ValueError, match="scaled natural-lifetime constant"):
        signed_good_generated_epoch_telescope((row0, row1))


def test_tiny_native_time_cannot_move_common_surfaces_forward():
    c = 1.0e-120
    with pytest.raises((ValueError, AssertionError), match="native backward-time"):
        _step(10.0, 6.1, c, 1.0e-13, 1.0e-13)


def test_native_unit_covariance_accepts_the_same_dimensionless_epoch():
    c = 1.0e-180
    child = 1.0e-60
    parent = 0.61 * child
    Tchild = math.exp(math.log(c) - 2.0 * math.log(child))
    row0 = _step(child, parent, c, 100.0 * Tchild, 100.25 * Tchild)
    Tchild1 = math.exp(math.log(c) - 2.0 * math.log(parent))
    row1 = _step(
        parent,
        0.61 * parent,
        c,
        row0.work_support_end - 0.2 * Tchild1,
        row0.work_support_end,
    )
    out = signed_good_generated_epoch_telescope((row0, row1))
    assert out.cumulative_reference_backshift >= out.minimum_cumulative_backshift
