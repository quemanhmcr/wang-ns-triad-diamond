from dataclasses import replace

import pytest

from src.continuum_master_event_quotient import energy_reentry_master_route
from src.material_label_carrier_quotient import (
    MATERIAL_MEMBERSHIP_EVENT,
    SELECTED_FAMILY_EVENT,
    selected_family_switch_sidecar,
)
from src.material_sidecar_stock_central_routing import (
    material_sidecar_joint_stop_projection,
    material_sidecar_stock_central_relay,
)
from src.material_sidecar_stock_owner_decomposition import (
    MEMBERSHIP_PROVENANCE_CURRENCY,
    SELECTED_FAMILY_MOYAL_CURRENCY,
    material_sidecar_stock_decomposition,
)
from src.physical_branch_compiler import CauseHit, MasterDisposition, PhysicalCause
from src.same_carrier_inherited_energy_relay import (
    SAME_CARRIER_INHERITED_STOCK_RELAY,
    SameCarrierInheritedEnergyRelayCertificate,
)


def _inherited(*events: str, boundary: float = 0.0) -> SameCarrierInheritedEnergyRelayCertificate:
    sidecars = tuple(sorted(events))
    return SameCarrierInheritedEnergyRelayCertificate(
        carrier_id="physical-helical-carrier-(7,6,5)+",
        initial_time=1.0,
        terminal_time=1.25,
        initial_energy=0.31,
        terminal_energy=1.0,
        inherited_fraction=0.31,
        residual_positive_work=0.12,
        residual_owner_threshold=0.2,
        observed_elapsed=0.25,
        analysis_segments=1,
        inserted_checkpoint_boundaries=0,
        material_sidecars=sidecars,
        selected_family_switch_energy=boundary,
    )


def _decomposition(cert: SameCarrierInheritedEnergyRelayCertificate):
    switch = None
    if SELECTED_FAMILY_EVENT in cert.material_sidecars:
        # Exact Moyal cell construction: the symmetric difference is cell 0.
        switch = selected_family_switch_sidecar(
            (cert.selected_family_switch_energy, 0.17, 0.09),
            (),
            (0,),
        )
    return material_sidecar_stock_decomposition(cert, selected_family_switch_certificate=switch)


def _reentry(cert=None, decomposition=None):
    out = {
        "branch": "material_energy_inheritance",
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
        "value": 0.31,
        "threshold": 0.2,
        "classified_residual_positive_work": 0.12,
    }
    if cert is not None:
        out["same_carrier_inherited_energy_relay_certificate"] = cert
    if decomposition is not None:
        out["material_sidecar_stock_decomposition"] = decomposition
    return out


def test_membership_sidecar_is_preserved_next_to_one_stock_charge_without_recursive_event():
    cert = _inherited(MATERIAL_MEMBERSHIP_EVENT)
    dec = _decomposition(cert)
    route = energy_reentry_master_route("same carrier inherited energy", 0.31, _reentry(cert, dec))
    relay = route.material_sidecar_stock_relay_certificate
    assert route.owner_bundle is None
    assert not route.recursive_event_created
    assert route.between_time_stock_relays == (SAME_CARRIER_INHERITED_STOCK_RELAY,)
    assert route.same_carrier_inherited_energy_certificate is cert
    assert relay is not None
    assert relay.stock_charge_count == 1
    assert relay.sidecar_events == (MATERIAL_MEMBERSHIP_EVENT,)
    assert relay.sidecar_currencies == (MEMBERSHIP_PROVENANCE_CURRENCY,)
    assert relay.selected_family_boundary_energy == 0.0
    assert not relay.physical_work_created
    assert not relay.sidecar_promoted_to_physical_hit


def test_selected_family_positive_moyal_boundary_is_preserved_but_does_not_create_a_stop():
    cert = _inherited(SELECTED_FAMILY_EVENT, boundary=0.23)
    dec = _decomposition(cert)
    route = energy_reentry_master_route("same carrier inherited energy", 0.31, _reentry(cert, dec))
    relay = route.material_sidecar_stock_relay_certificate
    assert relay is not None
    assert relay.sidecar_currencies == (SELECTED_FAMILY_MOYAL_CURRENCY,)
    assert relay.selected_family_boundary_energy == pytest.approx(0.23)
    assert route.owner_bundle is None

    with pytest.raises(ValueError, match="no causal stop"):
        material_sidecar_joint_stop_projection(relay)


def test_real_source_stop_keeps_moyal_boundary_as_non_event_provenance_only():
    cert = _inherited(SELECTED_FAMILY_EVENT, boundary=0.23)
    route = energy_reentry_master_route("same carrier inherited energy", 0.31, _reentry(cert, _decomposition(cert)))
    relay = route.material_sidecar_stock_relay_certificate
    assert relay is not None
    out = material_sidecar_joint_stop_projection(
        relay,
        physical_hits=(CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 1.0, "actual source service"),),
    )
    assert out.first_time == pytest.approx(0.4)
    assert out.joint_physical_causes == (PhysicalCause.RESOLVED_SOURCE.value,)
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value
    assert out.non_event_material_sidecar_events == (SELECTED_FAMILY_EVENT,)
    assert out.non_event_material_sidecar_currencies == (SELECTED_FAMILY_MOYAL_CURRENCY,)
    assert out.selected_family_boundary_energy == pytest.approx(0.23)


def test_raw_material_state_locator_cannot_be_promoted_by_moyal_boundary_sidecar():
    cert = _inherited(SELECTED_FAMILY_EVENT, boundary=0.23)
    route = energy_reentry_master_route("same carrier inherited energy", 0.31, _reentry(cert, _decomposition(cert)))
    relay = route.material_sidecar_stock_relay_certificate
    assert relay is not None
    with pytest.raises(TypeError, match="carrier/material-state locator"):
        material_sidecar_joint_stop_projection(
            relay,
            physical_hits=(CauseHit(0.4, PhysicalCause.MATERIAL_RELINK, 3.0, "unresolved material-state exit"),),
        )


def test_exact_source_strain_tie_is_not_fractionalized_or_erased_by_material_sidecars():
    cert = _inherited(MATERIAL_MEMBERSHIP_EVENT, SELECTED_FAMILY_EVENT, boundary=0.23)
    route = energy_reentry_master_route("same carrier inherited energy", 0.31, _reentry(cert, _decomposition(cert)))
    relay = route.material_sidecar_stock_relay_certificate
    hits = (
        CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 1e-12),
        CauseHit(0.4, PhysicalCause.HIGH_STRAIN_DISSIPATION, 1e12),
    )
    out = material_sidecar_joint_stop_projection(relay, physical_hits=hits)
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value
    assert set(out.joint_physical_causes) == {
        PhysicalCause.RESOLVED_SOURCE.value,
        PhysicalCause.HIGH_STRAIN_DISSIPATION.value,
    }
    assert not out.fine_rn_split_required
    assert out.non_event_material_sidecar_events == cert.material_sidecars
    assert out.selected_family_boundary_energy == pytest.approx(0.23)


def test_missing_phase_a_decomposition_stays_event_facing_fail_closed():
    cert = _inherited(SELECTED_FAMILY_EVENT, boundary=0.23)
    route = energy_reentry_master_route("unresolved sidecar-bearing inheritance", 0.31, _reentry(cert))
    assert route.owner_bundle is not None
    assert route.owner_bundle.owners == ("material_energy_inheritance",)
    assert route.between_time_stock_relays == ()
    assert route.material_sidecar_stock_relay_certificate is None


def test_transplanted_sidecar_decomposition_is_rejected_instead_of_silently_reclassified():
    cert = _inherited(SELECTED_FAMILY_EVENT, boundary=0.23)
    dec = _decomposition(cert)
    foreign = replace(dec, carrier_id="different-physical-carrier")
    with pytest.raises(TypeError, match="transplanted"):
        material_sidecar_stock_central_relay(cert, foreign)
    with pytest.raises(TypeError, match="transplanted"):
        energy_reentry_master_route("same carrier inherited energy", 0.31, _reentry(cert, foreign))


def test_sidecar_decomposition_without_inherited_certificate_is_not_a_causal_substitute():
    cert = _inherited(MATERIAL_MEMBERSHIP_EVENT)
    with pytest.raises(TypeError, match="without its inherited-stock certificate"):
        energy_reentry_master_route("bad detached sidecar", 0.31, _reentry(None, _decomposition(cert)))
