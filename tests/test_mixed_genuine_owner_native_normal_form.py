import pytest

from src.material_service_native_owner_factorization import project_material_recurrence_to_native_owners
from src.mixed_genuine_owner_native_normal_form import (
    NativeEventRoot,
    certified_objective_source_normal_form,
    certified_resolved_contact_normal_form,
    native_owner_normal_form,
    theorem_certificate,
)
from src.physical_branch_compiler import PhysicalCause
from src.resolved_contact_native_binding import K_RELAY


def test_material_sidecars_do_not_create_normal_form_event_root():
    material = project_material_recurrence_to_native_owners(
        membership_reread=True,
        selected_family_boundary_energy=2.0,
    )
    out = native_owner_normal_form(material_projection=material)
    assert out.event_roots == ()
    assert out.recursive_core_roots == ()
    assert out.material_sidecars


def test_legacy_material_roots_fail_closed():
    for cause in (PhysicalCause.MATERIAL_RELINK, PhysicalCause.NEW_COHERENT_ANCESTRY):
        with pytest.raises(TypeError, match="cannot survive"):
            native_owner_normal_form(physical_causes=(cause,))


def test_source_strain_and_actual_hh_are_the_only_recursive_core_classes_in_this_projection():
    out = native_owner_normal_form(
        physical_causes=(PhysicalCause.RESOLVED_SOURCE, PhysicalCause.HIGH_STRAIN_DISSIPATION),
        internal_hh_regeneration=True,
    )
    assert set(out.recursive_core_roots) == {
        NativeEventRoot.RESOLVED_SOURCE.value,
        NativeEventRoot.STRAIN_DISSIPATION.value,
        NativeEventRoot.ACTUAL_HH_OR_HARD_TAIL.value,
    }
    assert out.terminal_roots == ()


def test_transfer_reuse_sideband_boundary_and_uniform_resource_are_terminal_classes():
    out = native_owner_normal_form(
        physical_causes=(
            PhysicalCause.TRANSFER_WORK_LOSS,
            PhysicalCause.CAUSAL_REUSE,
            PhysicalCause.INTRINSIC_SIDEBAND,
            PhysicalCause.UNIFORM_GLOBAL_RESOURCE,
            PhysicalCause.INITIAL_BOUNDARY,
        )
    )
    assert out.recursive_core_roots == ()
    assert len(out.terminal_roots) == 5


def test_certified_resolved_contact_k_only_is_zero_depth_relay():
    out = certified_resolved_contact_normal_form(
        physical_tail_dissipation=8.0,
        viscosity=1.0,
        actual_contact_common_work=4.0,
        hh_complement_common_work=0.0,
        mixed_common_work=4.0,
        positive_skew_common_work=4.0,
        positive_strain_common_work=0.0,
    )
    assert out.event_roots == ()
    assert out.zero_depth_relays == (K_RELAY,)


def test_certified_resolved_contact_strain_reduces_to_existing_strain_root():
    out = certified_resolved_contact_normal_form(
        physical_tail_dissipation=8.0,
        viscosity=1.0,
        actual_contact_common_work=4.0,
        hh_complement_common_work=0.0,
        mixed_common_work=4.0,
        positive_skew_common_work=0.0,
        positive_strain_common_work=4.0,
    )
    assert out.event_roots == (NativeEventRoot.STRAIN_DISSIPATION.value,)


def test_certified_resolved_contact_hh_reduces_to_actual_hh_root():
    out = certified_resolved_contact_normal_form(
        physical_tail_dissipation=8.0,
        viscosity=1.0,
        actual_contact_common_work=4.0,
        hh_complement_common_work=4.0,
        mixed_common_work=0.0,
        positive_skew_common_work=0.0,
        positive_strain_common_work=0.0,
    )
    assert out.event_roots == (NativeEventRoot.ACTUAL_HH_OR_HARD_TAIL.value,)


def test_candidate_scope_keeps_surviving_three_root_mixed_problem_open():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert set(cert["recursive_core"]) == {
        NativeEventRoot.RESOLVED_SOURCE.value,
        NativeEventRoot.STRAIN_DISSIPATION.value,
        NativeEventRoot.ACTUAL_HH_OR_HARD_TAIL.value,
    }
    assert "not yet proved finite" in cert["scope"]


def _source_normal_form(weights):
    return certified_objective_source_normal_form(
        objective_variation_action=1.0,
        scaled_lifetime=1.0,
        owner_weights=weights,
        viscosity=1.0,
        filter_l1=1.0,
        lp_constant=1.0,
        bernstein_constant=1.0,
    )


def test_local_objective_source_owner_refines_to_existing_strain_dissipation_root():
    out = _source_normal_form({"local_dv": 1.0, "pressure": 0.0, "sgs": 0.0, "viscous": 0.0})
    assert out.recursive_core_roots == (NativeEventRoot.STRAIN_DISSIPATION.value,)


def test_viscous_objective_source_owner_refines_to_existing_strain_dissipation_root():
    out = _source_normal_form({"local_dv": 0.0, "pressure": 0.0, "sgs": 0.0, "viscous": 1.0})
    assert out.recursive_core_roots == (NativeEventRoot.STRAIN_DISSIPATION.value,)


def test_pressure_and_sgs_stay_source_service_without_realized_positive_law_route():
    for weights in (
        {"local_dv": 0.0, "pressure": 1.0, "sgs": 0.0, "viscous": 0.0},
        {"local_dv": 0.0, "pressure": 0.0, "sgs": 1.0, "viscous": 0.0},
    ):
        out = _source_normal_form(weights)
        assert out.recursive_core_roots == (NativeEventRoot.RESOLVED_SOURCE.value,)


def test_exact_local_sgs_owner_tie_keeps_both_native_classes_without_priority():
    out = _source_normal_form({"local_dv": 0.5, "pressure": 0.0, "sgs": 0.5, "viscous": 0.0})
    assert set(out.recursive_core_roots) == {
        NativeEventRoot.RESOLVED_SOURCE.value,
        NativeEventRoot.STRAIN_DISSIPATION.value,
    }
