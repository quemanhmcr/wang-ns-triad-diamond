import pytest

from src.continuum_master_event_quotient import SupplierKind
from src.fresh_service_scale_reentry import fresh_service_scale_route
from src.joint_causal_stop_projection import InternalHit, InternalRecursiveCause
from src.material_source_service_native_normal_form import (
    FRESH_SERVICE_SHELL_RELAY,
    STATUS,
    pde_native_material_source_service_projection,
    theorem_certificate,
)
from src.physical_branch_compiler import CauseHit, MasterDisposition, PhysicalCause
from src.smooth_quadratic_carrier_interface import RELINK_OWNER
from src.smooth_relink_donor_quotient import (
    SMOOTH_RELINK_SAME_EVENT_RELAY,
    SmoothRelinkDonorCertificate,
)


def _smooth_relink_certificate() -> SmoothRelinkDonorCertificate:
    return SmoothRelinkDonorCertificate(
        relink_owner=RELINK_OWNER,
        recipient_roles=(0,),
        terminal_negative_net_donor_roles=(1,),
        maximum_shortest_donor_path_length=1,
        role_count=2,
        positive_relink_work=0.7,
        recipient_positive_incoming_flux=0.7,
        pair_antisymmetry_residual=0.0,
        row_binding_residual=0.0,
        total_relink_work_residual=0.0,
    )


def test_naked_material_relink_is_rejected_as_unresolved_pde_provenance():
    with pytest.raises(TypeError, match="not a primitive PDE generator"):
        pde_native_material_source_service_projection(
            physical_hits=(CauseHit(0.4, PhysicalCause.MATERIAL_RELINK, 1.0),)
        )


def test_naked_new_coherent_ancestry_is_rejected_as_unresolved_pde_provenance():
    with pytest.raises(TypeError, match="not a primitive PDE generator"):
        pde_native_material_source_service_projection(
            physical_hits=(CauseHit(0.4, PhysicalCause.NEW_COHERENT_ANCESTRY, 1.0),)
        )


def test_source_hit_requires_native_supplier_binding():
    with pytest.raises(TypeError, match="cannot remain an abstract recursive label"):
        pde_native_material_source_service_projection(
            physical_hits=(CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 2.0),)
        )


def test_source_hit_keeps_source_owner_once_native_supplier_is_bound():
    out = pde_native_material_source_service_projection(
        physical_hits=(CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 2.0),),
        source_supplier_kinds=(SupplierKind.RESOLVED_DISSIPATION,),
    )
    assert out.projection is not None
    assert out.recursive_physical_causes == (PhysicalCause.RESOLVED_SOURCE.value,)
    assert out.source_supplier_kinds == (SupplierKind.RESOLVED_DISSIPATION.value,)
    assert out.projection.master_disposition == MasterDisposition.RECURSE_CRITICAL.value


def test_conservative_smooth_relink_is_zero_depth_and_mints_no_stop():
    out = pde_native_material_source_service_projection(
        smooth_relink_certificate=_smooth_relink_certificate(),
    )
    assert out.projection is None
    assert out.recursive_physical_causes == ()
    assert out.recursive_internal_causes == ()
    assert out.zero_depth_relays == (SMOOTH_RELINK_SAME_EVENT_RELAY,)


def test_actual_hh_generation_remains_internal_recursive_owner():
    out = pde_native_material_source_service_projection(
        internal_hits=(InternalHit(0.35, InternalRecursiveCause.HH_REGENERATION, "actual positive HH work"),)
    )
    assert out.projection is not None
    assert out.recursive_internal_causes == (InternalRecursiveCause.HH_REGENERATION.value,)
    assert out.projection.master_disposition == MasterDisposition.RECURSE_CRITICAL.value


def test_exact_source_strain_hh_tie_is_retained_without_material_priority():
    out = pde_native_material_source_service_projection(
        physical_hits=(
            CauseHit(0.5, PhysicalCause.RESOLVED_SOURCE, 3.0, "objective source"),
            CauseHit(0.5, PhysicalCause.HIGH_STRAIN_DISSIPATION, 7.0, "strain first hit"),
        ),
        internal_hits=(InternalHit(0.5, InternalRecursiveCause.HH_REGENERATION, "positive HH work"),),
        source_supplier_kinds=(SupplierKind.PRESSURE_PAIR,),
    )
    assert out.projection is not None
    assert set(out.recursive_physical_causes) == {
        PhysicalCause.RESOLVED_SOURCE.value,
        PhysicalCause.HIGH_STRAIN_DISSIPATION.value,
    }
    assert out.recursive_internal_causes == (InternalRecursiveCause.HH_REGENERATION.value,)
    assert out.projection.master_disposition == MasterDisposition.RECURSE_CRITICAL.value


def test_fresh_nn_service_is_supplier_relay_to_shell_not_new_material_owner():
    route = fresh_service_scale_route(
        integrated_square_service_threshold=1.0,
        scaled_lifetime=1.0,
        block_frequency=16.0,
        fresh_band_weights={0: 0.3},
    )
    out = pde_native_material_source_service_projection(
        fresh_service_scale_route=route,
    )
    assert out.projection is None
    assert out.recursive_physical_causes == ()
    assert out.zero_depth_relays == (FRESH_SERVICE_SHELL_RELAY,)


def test_forged_fresh_service_route_is_rejected_fail_closed():
    with pytest.raises(TypeError, match="generic critical-shell"):
        pde_native_material_source_service_projection(
            fresh_service_scale_route={
                "next_owner": "material_relink",
                "master_semantics": "RECURSE_CRITICAL_VIA_GENERIC_SHELL",
                "hard_shell_mass_lower": 1.0,
            }
        )


def test_draft_certificate_does_not_overclaim_mixed_termination_or_regularity():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert set(cert["primitive_material_causes_rejected"]) == {
        PhysicalCause.MATERIAL_RELINK.value,
        PhysicalCause.NEW_COHERENT_ANCESTRY.value,
    }
    assert cert["claims_mixed_owner_termination"] is False
    assert cert["claims_global_regularity"] is False
