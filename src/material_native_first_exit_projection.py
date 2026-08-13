from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from src.material_label_carrier_quotient import (
    MATERIAL_MEMBERSHIP_EVENT,
    SELECTED_FAMILY_EVENT,
)
from src.material_service_native_owner_factorization import MaterialRecurrenceProjection
from src.physical_branch_compiler import PhysicalCause
from src.smooth_sgs_first_hit_extraction import (
    JointFirstExit,
    PhysicalPathMonitor,
    first_physical_corridor_exit,
)


STATUS = (
    "DRAFT_NATIVE_MATERIAL_FIRST_EXIT_FILTER__"
    "MATERIAL_MEMBERSHIP_AND_FAMILY_SIDECARS_CANNOT_ENTER_DEBUT_SET__"
    "PHYSICAL_INTERFACE_OR_WORK_MONITORS_ONLY"
)

LEGACY_MATERIAL_MONITOR_LABELS = frozenset(
    {
        "material_state_exit",
        "material_relink",
        "new_coherent_ancestry",
        PhysicalCause.MATERIAL_RELINK.value,
        PhysicalCause.NEW_COHERENT_ANCESTRY.value,
        MATERIAL_MEMBERSHIP_EVENT,
        SELECTED_FAMILY_EVENT,
    }
)


@dataclass(frozen=True)
class MaterialNativeFirstExit:
    physical_exit: JointFirstExit | None
    material_projection: MaterialRecurrenceProjection
    no_physical_exit: bool
    material_attached_after_debut: bool = True
    material_changed_first_time: bool = False
    material_joined_first_stop_set: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.material_projection, MaterialRecurrenceProjection):
            raise TypeError("typed material native-owner projection required")
        if self.no_physical_exit != (self.physical_exit is None or self.physical_exit.first_time is None):
            raise ValueError("no-physical-exit flag disagrees with physical debut")
        if not self.material_attached_after_debut:
            raise ValueError("material provenance must be attached after the physical debut topology")
        if self.material_changed_first_time or self.material_joined_first_stop_set:
            raise ValueError("material sidecar illegally altered physical first-exit topology")
        if self.material_projection.new_recursive_vertex_created:
            raise ValueError("material projection illegally created a recursive vertex")

    @property
    def first_time(self) -> float | None:
        if self.physical_exit is None:
            return None
        return self.physical_exit.first_time

    @property
    def joint_first_stops(self) -> tuple[str, ...]:
        if self.physical_exit is None:
            return ()
        return self.physical_exit.joint_first_stops


def _validate_native_monitors(monitors: Sequence[PhysicalPathMonitor]) -> tuple[PhysicalPathMonitor, ...]:
    out = tuple(monitors)
    labels = tuple(m.label for m in out)
    bad = tuple(sorted(set(labels) & LEGACY_MATERIAL_MONITOR_LABELS))
    if bad:
        raise TypeError(
            "material membership/family bookkeeping is not a native PhysicalPathMonitor: "
            + ", ".join(bad)
            + "; register an independently physical role-interface/work/source/strain monitor instead"
        )
    return out


def material_native_first_physical_exit(
    times: Sequence[float],
    monitors: Sequence[PhysicalPathMonitor],
    material_projection: MaterialRecurrenceProjection,
    *,
    tie_tolerance: float = 0.0,
) -> MaterialNativeFirstExit:
    """Take the debut only over native physical monitors, then attach material data.

    Pure membership rereading, selected-family switching, and OO/ON/NN service
    restrictions cannot be inserted as corridor faces. A genuine Q/probe/role
    change must expose its already-registered physical interface/coefficient/work
    observable under its own non-material label; this wrapper does not hide such
    a monitor.
    """
    if not isinstance(material_projection, MaterialRecurrenceProjection):
        raise TypeError("typed material native-owner projection required")
    native = _validate_native_monitors(monitors)
    if not native:
        return MaterialNativeFirstExit(
            physical_exit=None,
            material_projection=material_projection,
            no_physical_exit=True,
        )
    out = first_physical_corridor_exit(times, native, tie_tolerance=tie_tolerance)
    return MaterialNativeFirstExit(
        physical_exit=out,
        material_projection=material_projection,
        no_physical_exit=out.first_time is None,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "debut_rule": "the first-exit minimum is taken only over already-registered native physical observables; material membership/family/restriction data are attached afterward",
        "forbidden_faces": tuple(sorted(LEGACY_MATERIAL_MONITOR_LABELS)),
        "genuine_role_change": "a genuine Q/probe/role change remains visible only through its independently registered role-interface/coefficient/work observable, never through material membership alone",
        "material_only": "with no native physical monitor, a material-only observation creates no first exit and no recursive vertex",
        "tie_rule": "material provenance cannot join or alter an exact tie among native physical corridor faces",
        "scope": "draft upstream topology theorem; it removes synthetic material debuts but does not terminate remaining native-owner recurrence or claim Navier-Stokes regularity",
    }
