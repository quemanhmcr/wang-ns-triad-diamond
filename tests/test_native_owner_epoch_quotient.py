from src.material_native_joint_stop_projection import material_native_joint_stop_projection
from src.material_service_native_owner_factorization import project_material_recurrence_to_native_owners
from src.native_owner_epoch_quotient import (
    NativeEpochKind,
    NativeEpochRecord,
    quotient_material_observations_from_epoch_path,
    theorem_certificate,
)
from src.physical_branch_compiler import CauseHit, PhysicalCause


def _material_only(t: float) -> NativeEpochRecord:
    material = project_material_recurrence_to_native_owners(membership_reread=True)
    projection = material_native_joint_stop_projection(material)
    return NativeEpochRecord(t, NativeEpochKind.MATERIAL_ONLY_OBSERVATION, projection)


def _high_strain(t: float) -> NativeEpochRecord:
    material = project_material_recurrence_to_native_owners()
    projection = material_native_joint_stop_projection(
        material,
        physical_hits=(CauseHit(t, PhysicalCause.HIGH_STRAIN_DISSIPATION, 1.0, "strain"),),
    )
    return NativeEpochRecord(t, NativeEpochKind.HIGH_STRAIN, projection)


def _other_source(t: float) -> NativeEpochRecord:
    material = project_material_recurrence_to_native_owners()
    projection = material_native_joint_stop_projection(
        material,
        physical_hits=(CauseHit(t, PhysicalCause.RESOLVED_SOURCE, 1.0, "source"),),
    )
    return NativeEpochRecord(t, NativeEpochKind.OTHER_NATIVE_EVENT, projection)


def test_material_observations_do_not_split_high_strain_epoch():
    records = (
        _high_strain(1.0),
        _material_only(0.9),
        _material_only(0.8),
        _high_strain(0.7),
        _material_only(0.6),
        _high_strain(0.5),
    )
    out = quotient_material_observations_from_epoch_path(records)
    assert out.event_kinds == (
        NativeEpochKind.HIGH_STRAIN.value,
        NativeEpochKind.HIGH_STRAIN.value,
        NativeEpochKind.HIGH_STRAIN.value,
    )
    assert out.high_strain_runs == ((0, 3),)
    assert out.material_observations_removed == 3
    assert out.material_epoch_breakers_created == 0


def test_real_source_event_still_breaks_high_strain_epoch():
    records = (
        _high_strain(1.0),
        _material_only(0.9),
        _other_source(0.8),
        _material_only(0.7),
        _high_strain(0.6),
    )
    out = quotient_material_observations_from_epoch_path(records)
    assert out.event_kinds == (
        NativeEpochKind.HIGH_STRAIN.value,
        NativeEpochKind.OTHER_NATIVE_EVENT.value,
        NativeEpochKind.HIGH_STRAIN.value,
    )
    assert out.high_strain_runs == ((0, 1), (2, 3))


def test_material_observations_do_not_split_signed_good_hh_epoch():
    material = project_material_recurrence_to_native_owners(membership_reread=True)
    no_stop = material_native_joint_stop_projection(material)
    records = (
        NativeEpochRecord(1.0, NativeEpochKind.SIGNED_GOOD_GENERATED_HH),
        NativeEpochRecord(0.8, NativeEpochKind.MATERIAL_ONLY_OBSERVATION, no_stop),
        NativeEpochRecord(0.6, NativeEpochKind.SIGNED_GOOD_GENERATED_HH),
    )
    out = quotient_material_observations_from_epoch_path(records)
    assert out.signed_good_hh_runs == ((0, 2),)
    assert out.event_times == (1.0, 0.6)


def test_candidate_scope_keeps_mixed_native_word_open():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "does not prove" in cert["scope"]
