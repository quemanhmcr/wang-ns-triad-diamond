import pytest

from src.joint_causal_stop_projection import InternalHit, InternalRecursiveCause
from src.material_native_joint_stop_projection import material_native_joint_stop_projection, theorem_certificate
from src.material_service_native_owner_factorization import NativePositiveServiceAtom, project_material_recurrence_to_native_owners
from src.physical_branch_compiler import CauseHit, MasterDisposition, PhysicalCause


def _sidecar_only():
    return project_material_recurrence_to_native_owners(
        membership_reread=True,
        selected_family_boundary_energy=2.0,
    )


def test_material_only_observation_creates_no_stop_and_is_not_kelvin_flat():
    out = material_native_joint_stop_projection(_sidecar_only())
    assert out.no_causal_stop
    assert out.physical_projection is None
    assert out.first_time is None
    assert out.joint_physical_causes == ()
    assert out.joint_internal_causes == ()
    assert out.master_disposition is None


def test_material_provenance_does_not_move_native_source_first_time():
    material = project_material_recurrence_to_native_owners(
        service_atoms=(NativePositiveServiceAtom(3.0, "resolved_source_or_sgs_service", False, False),)
    )
    out = material_native_joint_stop_projection(
        material,
        physical_hits=(CauseHit(0.4, PhysicalCause.RESOLVED_SOURCE, 3.0, "native source"),),
    )
    assert not out.no_causal_stop
    assert out.first_time == 0.4
    assert out.joint_physical_causes == (PhysicalCause.RESOLVED_SOURCE.value,)
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value
    assert not out.material_changed_first_time
    assert not out.material_joined_physical_tie


def test_exact_source_high_strain_hh_tie_stays_unsplit_after_material_attachment():
    material = project_material_recurrence_to_native_owners(membership_reread=True)
    out = material_native_joint_stop_projection(
        material,
        physical_hits=(
            CauseHit(0.5, PhysicalCause.RESOLVED_SOURCE, 1.0, "source"),
            CauseHit(0.5, PhysicalCause.HIGH_STRAIN_DISSIPATION, 1.0, "strain"),
        ),
        internal_hits=(InternalHit(0.5, InternalRecursiveCause.HH_REGENERATION, "HH route"),),
    )
    assert out.first_time == 0.5
    assert set(out.joint_physical_causes) == {
        PhysicalCause.RESOLVED_SOURCE.value,
        PhysicalCause.HIGH_STRAIN_DISSIPATION.value,
    }
    assert out.joint_internal_causes == (InternalRecursiveCause.HH_REGENERATION.value,)
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value


@pytest.mark.parametrize("cause", [PhysicalCause.MATERIAL_RELINK, PhysicalCause.NEW_COHERENT_ANCESTRY])
def test_legacy_material_roots_fail_closed(cause: PhysicalCause):
    with pytest.raises(TypeError, match="legacy material/new-ancestry CauseHit rejected"):
        material_native_joint_stop_projection(
            _sidecar_only(),
            physical_hits=(CauseHit(0.3, cause, 2.0, "legacy material root"),),
        )


def test_native_hh_internal_stop_survives_material_sidecars_without_replacement_event():
    out = material_native_joint_stop_projection(
        _sidecar_only(),
        internal_hits=(InternalHit(0.25, InternalRecursiveCause.HH_REGENERATION, "actual HH route"),),
    )
    assert out.first_time == 0.25
    assert out.joint_physical_causes == ()
    assert out.joint_internal_causes == (InternalRecursiveCause.HH_REGENERATION.value,)
    assert out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value


def test_candidate_certificate_does_not_claim_mixed_recurrence_termination():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "no remaining mixed-owner telescope" in cert["scope"]
