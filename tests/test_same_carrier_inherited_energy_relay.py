from dataclasses import replace

import pytest

from src.continuum_master_event_quotient import canonical_owner_bundle, energy_reentry_master_route
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


def _sidecar_free_relay():
    return same_carrier_inherited_energy_relay(
        _segments(),
        initial_time=1.6,
        terminal_time=2.0,
        initial_energy=0.31,
        terminal_energy=1.0,
        residual_positive_work=0.12,
        strain_action=0.02,
        initial_endpoint_is_non_event_carrier_slice=True,
    )


def _inheritance_reentry(cert=None, *, residual_work=0.12):
    gate = route_physical_energy_causality(
        terminal_energy=1.0,
        initial_energy=0.31,
        residual_positive_work=residual_work,
        strain_action=0.02,
    )
    out = {
        **gate,
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
        "classified_residual_positive_work": residual_work,
    }
    if cert is not None:
        out["same_carrier_inherited_energy_relay_certificate"] = cert
    return out


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
    assert cert.selected_family_switch_energy == pytest.approx(0.01)
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


def test_typed_master_projection_preserves_sidecars_separately_from_stock():
    cert = _relay()
    projection = same_carrier_inheritance_master_projection(
        "same smooth carrier energy stock",
        0.31,
        _inheritance_reentry(residual_work=0.12),
        cert,
    )
    assert projection.between_time_stock_relays == (SAME_CARRIER_INHERITED_STOCK_RELAY,)
    assert projection.sidecar_events == cert.material_sidecars
    assert not projection.stock_owner_bundle_created
    assert not projection.stock_recursive_event_created
    assert not projection.sidecars_quotiented_as_stock


def test_central_master_collapses_only_sidecar_free_typed_inheritance_to_stock_zero_depth():
    cert = _sidecar_free_relay()
    route = energy_reentry_master_route(
        "same smooth carrier inherited stock",
        0.31,
        _inheritance_reentry(cert, residual_work=0.12),
    )
    assert route.owner_bundle is None
    assert route.recursive_event_created is False
    assert route.same_event_relays == ()
    assert route.between_time_stock_relays == (SAME_CARRIER_INHERITED_STOCK_RELAY,)
    assert route.same_carrier_inherited_energy_certificate is cert


def test_central_master_keeps_sidecar_bearing_typed_inheritance_event_facing_fail_closed():
    cert = _relay()
    route = energy_reentry_master_route(
        "sidecar-bearing inherited carrier branch",
        0.31,
        _inheritance_reentry(cert, residual_work=0.12),
    )
    assert route.owner_bundle is not None
    assert route.owner_bundle.owners == ("material_energy_inheritance",)
    assert route.recursive_event_created is True
    assert route.between_time_stock_relays == ()
    assert route.same_carrier_inherited_energy_certificate is None
    assert route.material_sidecar_stock_relay_certificate is None


def test_central_master_keeps_untyped_inheritance_event_facing():
    route = energy_reentry_master_route(
        "untyped inherited carrier energy",
        0.31,
        _inheritance_reentry(),
    )
    assert route.owner_bundle is not None
    assert route.owner_bundle.owners == ("material_energy_inheritance",)
    assert route.recursive_event_created is True
    assert route.between_time_stock_relays == ()


def test_stock_relay_label_cannot_reenter_recursive_owner_bundle():
    with pytest.raises(TypeError, match="between-time physical stock continuation"):
        canonical_owner_bundle("bad stock promotion", 0.31, (SAME_CARRIER_INHERITED_STOCK_RELAY,))

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


def test_central_master_rejects_inherited_stock_certificate_with_a_different_master_mass():
    cert = _sidecar_free_relay()
    with pytest.raises(TypeError, match="mass must equal"):
        energy_reentry_master_route(
            "wrongly reweighted inherited stock",
            0.62,
            _inheritance_reentry(cert, residual_work=0.12),
        )


def test_central_master_rejects_sidecar_free_certificate_transplanted_to_different_residual_work():
    cert = _sidecar_free_relay()
    with pytest.raises(TypeError, match="classified residual physical-work law"):
        energy_reentry_master_route(
            "mismatched stock certificate",
            0.31,
            _inheritance_reentry(cert, residual_work=0.05),
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
