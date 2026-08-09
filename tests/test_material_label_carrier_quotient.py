import math

import numpy as np

from src.coherent_transfer_cells import symmetric_difference_energy
from src.common_slice_coefficient_registration import registration_first_stop
from src.material_label_carrier_quotient import (
    HH_STOP,
    INTERFACE_STOP,
    MATERIAL_MEMBERSHIP_EVENT,
    ROLE_DELEGATE_EVENT,
    SELECTED_FAMILY_EVENT,
    carrier_registration_with_material_sidecars,
    legacy_relink_refinement_certificate,
    reclassify_positive_service_sidecar,
    selected_family_switch_sidecar,
    theorem_certificate,
)


def test_pure_material_membership_change_does_not_add_carrier_impulse():
    ze = 1.0 + 0j
    ih = 0.1 + 0j
    ir = 0.1j
    zs = ze - ih - ir
    out = carrier_registration_with_material_sidecars(
        ze,
        zs,
        ih,
        ir,
        intrinsic_material_membership_change=True,
    )
    assert out["carrier_continuation_certified"] is True
    assert out["carrier_stop_causes"] == ()
    assert out["sidecar_events"] == (MATERIAL_MEMBERSHIP_EVENT,)
    assert out["slice_amplitude"] >= 0.25


def test_selected_family_switch_keeps_moyal_charge_but_same_carrier_can_continue():
    e = np.array([1.0, 2.0, 3.0, 4.0])
    old = [0, 1]
    new = [1, 2]
    sw = selected_family_switch_sidecar(e, old, new)
    ze = 2.0 + 0j
    ih = 0.2 + 0j
    ir = 0.1j
    out = carrier_registration_with_material_sidecars(
        ze,
        ze - ih - ir,
        ih,
        ir,
        selected_family_change=True,
        selected_family_switch_energy=float(sw["symmetric_difference_energy"]),
    )
    assert out["carrier_continuation_certified"] is True
    assert out["sidecar_events"] == (SELECTED_FAMILY_EVENT,)
    assert out["selected_family_switch_energy"] == symmetric_difference_energy(e, old, new)
    assert out["sidecar_requires_accounting"] is True
    assert out["same_carrier_reusable_after_sidecar"] is True


def test_material_sidecars_do_not_hide_simultaneous_interface_and_hh_carrier_stops():
    ze = 1.0 + 0j
    ir = 0.25 + 0j
    ih = 0.5j
    out = carrier_registration_with_material_sidecars(
        ze,
        ze - ir - ih,
        ih,
        ir,
        intrinsic_material_membership_change=True,
        selected_family_change=True,
        selected_family_switch_energy=3.0,
    )
    assert out["carrier_continuation_certified"] is False
    assert set(out["carrier_stop_causes"]) == {INTERFACE_STOP, HH_STOP}
    assert set(out["sidecar_events"]) == {MATERIAL_MEMBERSHIP_EVENT, SELECTED_FAMILY_EVENT}
    assert len(out["joint_physical_events"]) == 4
    assert out["primary_selected"] is False


def test_true_role_or_probe_change_is_never_declared_transparent():
    out = carrier_registration_with_material_sidecars(
        1.0 + 0j,
        0j,
        0j,
        0j,
        intrinsic_material_membership_change=True,
        same_smooth_role=False,
        same_analysis_probe=True,
    )
    assert out["quotient_applicable"] is False
    assert out["carrier_continuation_certified"] is False
    assert out["carrier_stop_causes"] == (ROLE_DELEGATE_EVENT,)


def test_rereading_ownership_reclassifies_but_never_creates_service():
    w = np.array([1.0, 2.0, 3.0, 4.0])
    out = reclassify_positive_service_sidecar(
        w,
        [True, True, False, False],
        [True, False, True, False],
        [False, True, False, True],
        [False, True, False, False],
    )
    assert out["before"]["total"] == 10.0
    assert out["after"]["total"] == 10.0
    assert out["total_service_residual"] == 0.0
    assert out["category_delta_sum"] == 0.0
    assert out["service_created_by_relabel"] == 0.0


def test_selected_family_energy_jump_is_paid_by_preserved_symmetric_difference_charge():
    e = [1.0, 2.0, 3.0, 4.0]
    sw = selected_family_switch_sidecar(e, [0, 1], [1, 2])
    assert abs(sw["selection_energy_jump"]) <= sw["symmetric_difference_energy"]
    assert sw["jump_bound_margin"] >= 0.0


def test_new_quotient_refines_only_pure_label_subtype_of_legacy_conservative_stop():
    ze = 1.0 + 0j
    ih = 0.1 + 0j
    ir = 0.1j
    zs = ze - ih - ir
    legacy = registration_first_stop(ze, zs, ih, ir, material_relink=True)
    refined = carrier_registration_with_material_sidecars(
        ze,
        zs,
        ih,
        ir,
        intrinsic_material_membership_change=True,
        same_smooth_role=True,
        same_analysis_probe=True,
    )
    assert legacy["continuing"] is False
    assert refined["carrier_continuation_certified"] is True
    cert = legacy_relink_refinement_certificate()
    assert cert["old_theorem_status"] == "retained_as_conservative_superset"


def test_certificate_preserves_ancestry_charge_and_refuses_role_change_quotient():
    cert = theorem_certificate()
    assert "no second independent coefficient impulse" in cert["no_double_count"]
    assert "remains ancestry/service currency" in cert["selected_family"]
    assert "no transparency is claimed" in cert["nonquotient"]
    assert "does not erase material ancestry" in cert["master"]
