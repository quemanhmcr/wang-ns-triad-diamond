import pytest

from src.material_service_native_owner_factorization import (
    FORBIDDEN_NATIVE_OWNER_ROOTS,
    NativePositiveServiceAtom,
    factor_positive_service_by_material,
    project_material_recurrence_to_native_owners,
    theorem_certificate,
)
from src.material_sidecar_stock_owner_decomposition import (
    MEMBERSHIP_PROVENANCE_CURRENCY,
    SELECTED_FAMILY_MOYAL_CURRENCY,
)
from src.physical_branch_compiler import PhysicalCause
from src.smooth_quadratic_carrier_interface import RELINK_OWNER, STRAIN_OWNER
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
        positive_relink_work=1.0,
        recipient_positive_incoming_flux=1.0,
        pair_antisymmetry_residual=0.0,
        row_binding_residual=0.0,
        total_relink_work_residual=0.0,
    )


def test_material_partition_commutes_with_native_owner_disintegration():
    atoms = (
        NativePositiveServiceAtom(2.0, "resolved_source_or_sgs_service", True, True),
        NativePositiveServiceAtom(3.0, "resolved_source_or_sgs_service", True, False),
        NativePositiveServiceAtom(5.0, "actual_positive_hh_child_energy_work", False, False),
        NativePositiveServiceAtom(7.0, "existing_symmetric_strain_deformation", False, True),
    )
    out = factor_positive_service_by_material(atoms)
    assert out.total_service == 17.0
    assert out.oo_service == 2.0
    assert out.on_service == 10.0
    assert out.nn_service == 5.0
    assert out.owner_total == (
        ("actual_positive_hh_child_energy_work", 5.0),
        ("existing_symmetric_strain_deformation", 7.0),
        ("resolved_source_or_sgs_service", 5.0),
    )
    assert out.material_native_owner_provenance == (
        "actual_positive_hh_child_energy_work",
        "existing_symmetric_strain_deformation",
        "resolved_source_or_sgs_service",
    )
    assert out.ownership_partition_residual == 0.0
    assert out.owner_disintegration_residual == 0.0
    assert not out.recursive_material_owner_created
    assert not out.recursive_new_ancestry_owner_created
    assert not out.new_recursive_vertex_created
    assert not out.later_hahn_used


def test_fresh_nn_inherits_source_root_without_minting_new_ancestry_or_event():
    out = project_material_recurrence_to_native_owners(
        service_atoms=(
            NativePositiveServiceAtom(4.0, "resolved_source_or_sgs_service", False, False),
        )
    )
    assert out.native_owner_provenance == ("resolved_source_or_sgs_service",)
    assert PhysicalCause.NEW_COHERENT_ANCESTRY.value not in out.native_owner_provenance
    assert not out.recursive_new_ancestry_owner_created
    assert not out.new_recursive_vertex_created


def test_interface_on_inherits_existing_root_without_minting_material_relink_or_event():
    out = project_material_recurrence_to_native_owners(
        service_atoms=(
            NativePositiveServiceAtom(4.0, "actual_positive_hh_child_energy_work", True, False),
        )
    )
    assert out.native_owner_provenance == ("actual_positive_hh_child_energy_work",)
    assert PhysicalCause.MATERIAL_RELINK.value not in out.native_owner_provenance
    assert not out.recursive_material_owner_created
    assert not out.new_recursive_vertex_created


def test_membership_and_selected_family_are_sidecars_only():
    out = project_material_recurrence_to_native_owners(
        membership_reread=True,
        selected_family_boundary_energy=3.25,
    )
    assert out.native_owner_provenance == ()
    assert out.same_event_relays == ()
    assert set(out.sidecar_currencies) == {
        MEMBERSHIP_PROVENANCE_CURRENCY,
        SELECTED_FAMILY_MOYAL_CURRENCY,
    }
    assert out.selected_family_boundary_energy == 3.25
    assert not out.new_recursive_vertex_created


def test_smooth_kphys_relink_is_same_event_relay_not_recursive_material_owner():
    out = project_material_recurrence_to_native_owners(smooth_relink=_smooth_relink_certificate())
    assert out.native_owner_provenance == ()
    assert out.same_event_relays == (SMOOTH_RELINK_SAME_EVENT_RELAY,)
    assert not out.recursive_material_owner_created
    assert not out.new_recursive_vertex_created


def test_relink_strain_tie_keeps_only_strain_as_native_owner_provenance():
    out = project_material_recurrence_to_native_owners(
        smooth_relink=_smooth_relink_certificate(),
        additional_native_owner_provenance=(STRAIN_OWNER,),
        require_physical_role_change_owner=True,
    )
    assert out.native_owner_provenance == (STRAIN_OWNER,)
    assert out.same_event_relays == (SMOOTH_RELINK_SAME_EVENT_RELAY,)
    assert not out.new_recursive_vertex_created


def test_unbacked_role_change_fails_closed():
    with pytest.raises(TypeError, match="fail-closed"):
        project_material_recurrence_to_native_owners(require_physical_role_change_owner=True)


@pytest.mark.parametrize("owner", sorted(FORBIDDEN_NATIVE_OWNER_ROOTS))
def test_material_or_relink_bookkeeping_cannot_be_native_service_root(owner: str):
    with pytest.raises(ValueError):
        NativePositiveServiceAtom(1.0, owner, False, False)


def test_negative_service_atom_is_rejected_before_material_classification():
    with pytest.raises(ValueError):
        NativePositiveServiceAtom(-1.0, "resolved_source_or_sgs_service", False, False)


def test_candidate_certificate_states_narrow_scope_and_no_regularity_claim():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "material/new-ancestry labels do not form an independent letter" in cert["master_consequence"]
    assert "no global-regularity claim" in cert["scope"]
