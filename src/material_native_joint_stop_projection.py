from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from src.joint_causal_stop_projection import InternalHit, JointStopProjection, joint_stop_master_projection
from src.material_service_native_owner_factorization import MaterialRecurrenceProjection
from src.physical_branch_compiler import CauseHit, PhysicalCause, UniformResourceCertificate


STATUS = (
    "DRAFT_MATERIAL_NATIVE_JOINT_STOP_PROJECTION__"
    "PHYSICAL_FIRST_HIT_BEFORE_MATERIAL_PROVENANCE__"
    "LEGACY_MATERIAL_ROOTS_FAIL_CLOSED__"
    "SIDECAR_ONLY_OBSERVATIONS_CREATE_NO_CAUSAL_STOP"
)

LEGACY_MATERIAL_ROOTS = frozenset(
    {PhysicalCause.MATERIAL_RELINK, PhysicalCause.NEW_COHERENT_ANCESTRY}
)


@dataclass(frozen=True)
class MaterialNativeJointStopProjection:
    """A physical/internal joint stop with material provenance attached afterward."""

    physical_projection: JointStopProjection | None
    material_projection: MaterialRecurrenceProjection
    no_causal_stop: bool
    material_attached_after_first_hit: bool = True
    material_changed_first_time: bool = False
    material_joined_physical_tie: bool = False
    material_changed_disposition: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.material_projection, MaterialRecurrenceProjection):
            raise TypeError("typed material native-owner projection required")
        if self.no_causal_stop != (self.physical_projection is None):
            raise ValueError("no-causal-stop flag disagrees with physical projection")
        if not self.material_attached_after_first_hit:
            raise ValueError("material provenance must be attached only after physical first-hit classification")
        if (
            self.material_changed_first_time
            or self.material_joined_physical_tie
            or self.material_changed_disposition
        ):
            raise ValueError("material provenance altered physical first-hit topology")
        if self.material_projection.new_recursive_vertex_created:
            raise ValueError("material projection illegally created a recursive event vertex")

    @property
    def first_time(self) -> float | None:
        return None if self.physical_projection is None else self.physical_projection.first_time

    @property
    def joint_physical_causes(self) -> tuple[str, ...]:
        return () if self.physical_projection is None else self.physical_projection.joint_physical_causes

    @property
    def joint_internal_causes(self) -> tuple[str, ...]:
        return () if self.physical_projection is None else self.physical_projection.joint_internal_causes

    @property
    def master_disposition(self) -> str | None:
        return None if self.physical_projection is None else self.physical_projection.master_disposition


def _reject_legacy_material_hits(physical_hits: tuple[CauseHit, ...]) -> None:
    bad = tuple(hit for hit in physical_hits if hit.cause in LEGACY_MATERIAL_ROOTS)
    if bad:
        names = ", ".join(sorted({hit.cause.value for hit in bad}))
        raise TypeError(
            "legacy material/new-ancestry CauseHit rejected by native material projection: "
            f"{names}; supply the independently witnessed native PDE cause and keep material data as provenance"
        )


def material_native_joint_stop_projection(
    material_projection: MaterialRecurrenceProjection,
    *,
    physical_hits: tuple[CauseHit, ...] = (),
    internal_hits: tuple[InternalHit, ...] = (),
    fixed_transfer_loss: bool = False,
    uniform_certificates: Mapping[PhysicalCause, UniformResourceCertificate] | None = None,
) -> MaterialNativeJointStopProjection:
    """Run physical first-hit classification before attaching material provenance.

    This is deliberately a wrapper around the certified core. It does not alter
    the core enum or legacy API. New code opting into the material factorization
    must provide the native PDE hit itself; MATERIAL_RELINK and
    NEW_COHERENT_ANCESTRY are no longer admissible stand-alone roots here.

    A material-only observation produces no_causal_stop=True. It is not sent
    through kelvin_flat_certified=True because absence of a material cause is not
    a proof of Kelvin flatness.
    """
    if not isinstance(material_projection, MaterialRecurrenceProjection):
        raise TypeError("typed material native-owner projection required")
    _reject_legacy_material_hits(physical_hits)

    if not physical_hits and not internal_hits and not fixed_transfer_loss:
        return MaterialNativeJointStopProjection(
            physical_projection=None,
            material_projection=material_projection,
            no_causal_stop=True,
        )

    physical = joint_stop_master_projection(
        physical_hits=physical_hits,
        internal_hits=internal_hits,
        fixed_transfer_loss=fixed_transfer_loss,
        uniform_certificates=uniform_certificates,
    )
    return MaterialNativeJointStopProjection(
        physical_projection=physical,
        material_projection=material_projection,
        no_causal_stop=False,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "topology": "physical/internal first-hit projection runs before material provenance is attached; material data cannot change first time, exact tie set, or master disposition",
        "legacy_roots": "MATERIAL_RELINK and NEW_COHERENT_ANCESTRY CauseHit inputs fail closed in the new wrapper; callers must supply the independently witnessed native PDE cause",
        "material_only": "membership/family/OO-ON-NN observations with no independently witnessed physical/internal hit produce no causal stop rather than a synthetic Kelvin-flat or recursive event",
        "core_compatibility": "the certified joint_stop_master_projection core is unchanged; the candidate is an opt-in fail-closed wrapper",
        "ties": "exact ties among independently witnessed native physical/internal causes remain unsplit because the same core projection is used",
        "scope": "draft topology/integration theorem only; no remaining mixed-owner telescope or Navier-Stokes global-regularity claim",
    }
