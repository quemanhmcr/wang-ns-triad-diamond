import pytest

from src.material_label_carrier_quotient import MATERIAL_MEMBERSHIP_EVENT, SELECTED_FAMILY_EVENT
from src.material_native_first_exit_projection import (
    material_native_first_physical_exit,
    theorem_certificate,
)
from src.material_service_native_owner_factorization import project_material_recurrence_to_native_owners
from src.smooth_sgs_first_hit_extraction import PhysicalPathMonitor, ThresholdTopology


def _material():
    return project_material_recurrence_to_native_owners(
        membership_reread=True,
        selected_family_boundary_energy=1.0,
    )


def _monitor(label: str, values=(0.0, 1.0), threshold=0.5):
    return PhysicalPathMonitor(
        label=label,
        values=values,
        threshold=threshold,
        topology=ThresholdTopology.CLOSED,
    )


def test_material_only_observation_has_no_first_exit():
    out = material_native_first_physical_exit((0.0, 1.0), (), _material())
    assert out.no_physical_exit
    assert out.first_time is None
    assert out.joint_first_stops == ()


@pytest.mark.parametrize(
    "label",
    [
        "material_state_exit",
        "material_relink",
        "new_coherent_ancestry",
        MATERIAL_MEMBERSHIP_EVENT,
        SELECTED_FAMILY_EVENT,
    ],
)
def test_material_bookkeeping_cannot_be_inserted_as_first_exit_face(label: str):
    with pytest.raises(TypeError, match="not a native PhysicalPathMonitor"):
        material_native_first_physical_exit(
            (0.0, 1.0),
            (_monitor(label),),
            _material(),
        )


def test_native_strain_source_tie_is_unchanged_by_material_sidecar():
    times = (0.0, 1.0)
    monitors = (
        _monitor("strain_action"),
        _monitor("objective_source_action"),
    )
    out = material_native_first_physical_exit(times, monitors, _material())
    assert out.first_time == 0.5
    assert set(out.joint_first_stops) == {"strain_action", "objective_source_action"}
    assert not out.material_changed_first_time
    assert not out.material_joined_first_stop_set


def test_physical_role_interface_locator_is_not_hidden():
    out = material_native_first_physical_exit(
        (0.0, 1.0),
        (_monitor("role_interface_coefficient_obstruction"),),
        _material(),
    )
    assert out.first_time == 0.5
    assert out.joint_first_stops == ("role_interface_coefficient_obstruction",)


def test_scope_remains_narrow():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "does not terminate" in cert["scope"]
