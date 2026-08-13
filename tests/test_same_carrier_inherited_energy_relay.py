from dataclasses import replace

import pytest

from src.continuum_master_event_quotient import energy_reentry_master_route
from src.material_label_carrier_quotient import carrier_registration_with_material_sidecars
from src.physical_energy_causal_bridge import route_physical_energy_causality
from src.same_carrier_checkpoint_segmentation_quotient import partition_same_carrier_path
from src.same_carrier_inherited_energy_relay import (
    SAME_CARRIER_INHERITED_STOCK_RELAY,
    same_carrier_inheritance_master_projection,
    same_carrier_inherited_energy_relay,
    theorem_certificate,
)


def _segments(*, strain_end=0.02, residual_end=0.1, hh_end=0.2, cuts=(2,)):
    return partition_same_carrier_path(
        carrier_id="fixed-Q-psi",
        terminal_amplitude=1.0,
        elapsed_times=(0.0, 0.1, 0.2, 0.3, 0.4),
        strain_action=(0.0, 0.004, 0.009, 0.014, strain_end),
        residual_impulse_abs=(0.0, 0.03, 0.08, 0.07, residual_end),
        hh_impulse_abs=(0.0, 0.08, 0.16, 0.14, hh_end),
        checkpoint_indices=cuts,
    )


def _material(*, same_role=True, same_probe=True):
    z = 1.0 + 0.0j
    ih = 0.08 + 0.02j
    ii = 0.03 - 0.01j
    return carrier_registration_with_material_sidecars(
        z,
        z - ih - ii,
        ih,
        ii,
        intrinsic_material_membership_change=True,
        selected_family_change=True,
        selected_family_switch_energy=0.01,
        same_smooth_role=same_role,
        same_analysis_probe=same_probe,
    )


def _relay():
    return same_carrier_inherited_energy_relay(
        _segments(),
        initial_time=1.6,
        terminal_time=2.0,
        initial_energy=0.31,
        terminal_energy=1.0,
        residual_positive_work=0.12,
        strain_action=0.02,
        material_registration=_material(),
        initial_endpoint_is_non_event_carrier_slice=True,
    )


def test_certificate_states_stock_not_generation_and_no_temporal_matching():
    cert = theorem_certificate()
    assert "earlier physical carrier stock" in cert["physical_gate"]
    assert "zero generation depth" in cert["master_ontology"]
    assert "sidecars remain separately routed" in cert["master_ontology"]
    assert cert["later_hahn_used"] is False
    assert cert["claims_global_regularity"] is False


def test_same_carrier_inheritance_survives_checkpoints_and_material_sidecars_without_new_event():
    cert = _relay()
    assert cert.carrier_id == "fixed-Q-psi"
    assert cert.initial_energy == pytest.approx(0.31)
    assert cert.terminal_energy == pytest.approx(1.0)
    assert cert.inherited_fraction == pytest.approx(0.31)
    assert cert.residual_positive_work == pytest.approx(0.12)
    assert cert.residual_owner_threshold == pytest.approx(0.2)
    assert cert.analysis_segments == 2
    assert cert.inserted_checkpoint_boundaries == 1
    assert cert.material_sidecars
    assert not cert.recursive_generation_created
    assert not cert.new_causal_charge_created
    assert not cert.between_time_deposit_matching_used


def test_checkpoint_partition_does_not_change_same_carrier_stock_relay():
    a = same_carrier_inherited_energy_relay(
        _segments(cuts=()),
        initial_time=1.6,
        terminal_time=2.0,
        initial_energy=0.31,
        terminal_energy=1.0,
        residual_positive_work=0.12,
        strain_action=0.02,
        initial_endpoint_is_non_event_carrier_slice=True,
    )
    b = same_carrier_inherited_energy_relay(
        _segments(cuts=(1, 3)),
        initial_time=1.6,
        terminal_time=2.0,
        initial_energy=0.31,
        terminal_energy=1.0,
        residual_positive_work=0.12,
        strain_action=0.02,
        initial_endpoint_is_non_event_carrier_slice=True,
    )
    assert a.carrier_id == b.carrier_id
    assert a.inherited_fraction == pytest.approx(b.inherited_fraction)
    assert a.observed_elapsed == pytest.approx(b.observed_elapsed)
    assert a.inserted_checkpoint_boundaries == 0
    assert b.inserted_checkpoint_boundaries == 2


def test_named_first_stop_role_change_and_noninheritance_all_fail_closed():
    with pytest.raises(TypeError):
        same_carrier_inherited_energy_relay(
            _segments(strain_end=1.0 / 30.0),
            initial_time=1.6,
            terminal_time=2.0,
            initial_energy=0.31,
            terminal_energy=1.0,
            residual_positive_work=0.12,
            strain_action=1.0 / 30.0,
            initial_endpoint_is_non_event_carrier_slice=True,
        )
    with pytest.raises(TypeError):
        same_carrier_inherited_energy_relay(
            _segments(),
            initial_time=1.6,
            terminal_time=2.0,
            initial_energy=0.31,
            terminal_energy=1.0,
            residual_positive_work=0.12,
            strain_action=0.02,
            material_registration=_material(same_role=False),
            initial_endpoint_is_non_event_carrier_slice=True,
        )
    with pytest.raises(TypeError):
        same_carrier_inherited_energy_relay(
            _segments(),
            initial_time=1.6,
            terminal_time=2.0,
            initial_energy=0.10,
            terminal_energy=1.0,
            residual_positive_work=0.10,
            strain_action=0.02,
            initial_endpoint_is_non_event_carrier_slice=True,
        )


def test_typed_master_projection_makes_stock_zero_depth_while_legacy_master_stays_fail_closed():
    cert = _relay()
    gate = route_physical_energy_causality(
        terminal_energy=1.0,
        initial_energy=0.31,
        residual_positive_work=0.12,
        strain_action=0.02,
    )
    assert gate["branch"] == "material_energy_inheritance"
    projection = same_carrier_inheritance_master_projection(
        "same smooth carrier energy stock",
        0.31,
        {
            **gate,
            "coefficient_impulse_used_as_physical_work": False,
            "observer_partition_motion_charged_as_physics": False,
            "classified_residual_positive_work": 0.12,
        },
        cert,
    )
    assert projection.between_time_stock_relays == (SAME_CARRIER_INHERITED_STOCK_RELAY,)
    assert projection.sidecar_events == cert.material_sidecars
    assert not projection.stock_owner_bundle_created
    assert not projection.stock_recursive_event_created
    assert not projection.sidecars_quotiented_as_stock

    # Until the second-phase central-master wiring is certified, the old API is
    # deliberately fail-closed and retains untyped inheritance as an owner.
    legacy = energy_reentry_master_route(
        "untyped inherited carrier energy",
        0.31,
        {
            **gate,
            "coefficient_impulse_used_as_physical_work": False,
            "observer_partition_motion_charged_as_physics": False,
        },
    )
    assert legacy.recursive_event_created
    assert legacy.owner_bundle is not None
    assert legacy.owner_bundle.owners == ("material_energy_inheritance",)

def test_genuine_physical_event_at_earlier_endpoint_cannot_be_erased_into_stock_relay():
    with pytest.raises(TypeError):
        same_carrier_inherited_energy_relay(
            _segments(),
            initial_time=1.6,
            terminal_time=2.0,
            initial_energy=0.31,
            terminal_energy=1.0,
            residual_positive_work=0.12,
            strain_action=0.02,
            initial_endpoint_is_non_event_carrier_slice=False,
        )


def test_simultaneous_residual_owner_prevents_stock_only_quotient_even_when_E0_is_large():
    with pytest.raises(TypeError):
        same_carrier_inherited_energy_relay(
            _segments(),
            initial_time=1.6,
            terminal_time=2.0,
            initial_energy=0.80,
            terminal_energy=1.0,
            residual_positive_work=0.21,
            strain_action=0.02,
            initial_endpoint_is_non_event_carrier_slice=True,
        )


def test_stock_projection_rejects_same_endpoint_certificate_from_a_different_residual_work_law():
    cert = _relay()
    gate = route_physical_energy_causality(
        terminal_energy=1.0,
        initial_energy=0.31,
        residual_positive_work=0.05,
        strain_action=0.02,
    )
    with pytest.raises(TypeError):
        same_carrier_inheritance_master_projection(
            "same endpoints, different residual work",
            0.31,
            {
                **gate,
                "coefficient_impulse_used_as_physical_work": False,
                "observer_partition_motion_charged_as_physics": False,
                "classified_residual_positive_work": 0.05,
            },
            cert,
        )


def test_stock_projection_rejects_tampered_or_mismatched_energy_gate():
    cert = _relay()
    with pytest.raises(ValueError):
        replace(cert, later_hahn_used=True)

    gate = route_physical_energy_causality(
        terminal_energy=2.0,
        initial_energy=0.50,
        residual_positive_work=0.1,
        strain_action=0.02,
    )
    with pytest.raises(TypeError):
        same_carrier_inheritance_master_projection(
            "mismatched inherited stock",
            0.5,
            {
                **gate,
                "coefficient_impulse_used_as_physical_work": False,
                "observer_partition_motion_charged_as_physics": False,
                "classified_residual_positive_work": 0.10,
            },
            cert,
        )
