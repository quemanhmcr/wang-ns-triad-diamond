from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

from src.continuum_master_event_quotient import SupplierKind
from src.joint_causal_stop_projection import (
    InternalHit,
    JointStopProjection,
    joint_stop_master_projection,
)
from src.material_sidecar_stock_central_routing import (
    MaterialSidecarStockCentralRelayCertificate,
)
from src.physical_branch_compiler import (
    CauseHit,
    PhysicalCause,
    UniformResourceCertificate,
)
from src.smooth_relink_donor_quotient import (
    SMOOTH_RELINK_SAME_EVENT_RELAY,
    SmoothRelinkDonorCertificate,
)


STATUS = (
    "DRAFT_PDE_NATIVE_MATERIAL_SOURCE_SERVICE_NORMAL_FORM__"
    "NO_PRIMITIVE_MATERIAL_RELINK_GENERATOR__"
    "KPHYS_ZERO_DEPTH__SOURCE_STRAIN_HH_RETAIN_NATIVE_OWNERS__"
    "FRESH_SERVICE_RELAYS_TO_GENERIC_SHELL"
)

NON_PRIMITIVE_MATERIAL_CAUSES = frozenset(
    {
        PhysicalCause.MATERIAL_RELINK,
        PhysicalCause.NEW_COHERENT_ANCESTRY,
    }
)

SOURCE_SUPPLIER_KINDS = frozenset(
    {
        SupplierKind.RESOLVED_DISSIPATION,
        SupplierKind.PRESSURE_PAIR,
        SupplierKind.FRESH_SGS_SCALE,
        SupplierKind.HIGH_TAIL,
        SupplierKind.GENERIC_CRITICAL_SHELL,
        SupplierKind.MATERIAL_REUSE,
    }
)

FRESH_SERVICE_SHELL_RELAY = "fresh_NN_service_to_generic_critical_shell_first_stop"
MATERIAL_BOUNDARY_SIDECAR_RELAY = "material_boundary_sidecar_zero_generation_depth"


@dataclass(frozen=True)
class MaterialSourceServiceNativeNormalForm:
    """Strict master-facing normal form for material/source-service recurrence.

    This object is downstream of existing exact PDE identities.  It manufactures
    neither physical work nor a new first-stop clock; it records only which
    already-proved native causes remain after material bookkeeping and
    conservative same-event relink have been quotiented.
    """

    projection: JointStopProjection | None
    recursive_physical_causes: tuple[str, ...]
    recursive_internal_causes: tuple[str, ...]
    source_supplier_kinds: tuple[str, ...]
    zero_depth_relays: tuple[str, ...]
    sidecar_events: tuple[str, ...]
    primitive_material_generator_created: bool = False
    selected_family_boundary_used_as_work: bool = False
    conservative_relink_promoted_to_recursion: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        if self.primitive_material_generator_created:
            raise ValueError("native material normal form cannot mint a primitive material generator")
        if self.selected_family_boundary_used_as_work:
            raise ValueError("selected-family R_switch is a boundary sidecar, not physical work")
        if self.conservative_relink_promoted_to_recursion:
            raise ValueError("smooth K_phys relink is same-event donor provenance, not recursive generation")
        if self.later_hahn_used:
            raise ValueError("native material normal form cannot re-Hahn an inherited physical cause")
        forbidden = {x.value for x in NON_PRIMITIVE_MATERIAL_CAUSES}
        if forbidden.intersection(self.recursive_physical_causes):
            raise ValueError("material/new-ancestry manifestation survived as a primitive recursive cause")
        if tuple(sorted(set(self.zero_depth_relays))) != self.zero_depth_relays:
            raise ValueError("zero-depth relay set must be sorted and quotiented")
        if tuple(sorted(set(self.sidecar_events))) != self.sidecar_events:
            raise ValueError("material sidecar event set must be sorted and quotiented")
        if tuple(sorted(set(self.source_supplier_kinds))) != self.source_supplier_kinds:
            raise ValueError("source supplier set must be sorted and quotiented")


def _validate_fresh_service_scale_route(route: Mapping[str, object]) -> None:
    if route.get("next_owner") != "generic_critical_shell_first_stop":
        raise TypeError("fresh material service must relay to the certified generic critical-shell first stop")
    if route.get("master_semantics") != "RECURSE_CRITICAL_VIA_GENERIC_SHELL":
        raise TypeError("fresh material service changed the certified shell-reentry master semantics")
    mass = float(route.get("hard_shell_mass_lower", 0.0))
    if not mass > 0.0:
        raise TypeError("fresh material service requires a positive certified hard-shell seed")


def _validate_source_suppliers(
    physical_hits: tuple[CauseHit, ...],
    source_supplier_kinds: tuple[SupplierKind, ...],
) -> tuple[SupplierKind, ...]:
    suppliers = tuple(sorted(set(source_supplier_kinds), key=lambda x: x.value))
    if any(x not in SOURCE_SUPPLIER_KINDS for x in suppliers):
        raise TypeError("source-service recurrence used a supplier outside the certified PDE-facing supplier set")
    has_source = any(hit.cause is PhysicalCause.RESOLVED_SOURCE for hit in physical_hits)
    if has_source and not suppliers:
        raise TypeError(
            "resolved source cannot remain an abstract recursive label: bind it to a native dissipation/pressure/fresh/high-tail/shell/material-reuse supplier"
        )
    if suppliers and not has_source:
        raise TypeError("source supplier provenance was supplied without an independently witnessed resolved-source hit")
    return suppliers


def pde_native_material_source_service_projection(
    *,
    physical_hits: tuple[CauseHit, ...] = (),
    internal_hits: tuple[InternalHit, ...] = (),
    source_supplier_kinds: tuple[SupplierKind, ...] = (),
    material_sidecar_relay: MaterialSidecarStockCentralRelayCertificate | None = None,
    smooth_relink_certificate: SmoothRelinkDonorCertificate | None = None,
    fresh_service_scale_route: Mapping[str, object] | None = None,
    fixed_transfer_loss: bool = False,
    kelvin_flat_certified: bool = False,
    uniform_certificates: Mapping[PhysicalCause, UniformResourceCertificate] | None = None,
) -> MaterialSourceServiceNativeNormalForm:
    """Descend material/source-service recurrence to native PDE causes, fail closed.

    Upstream exact identities imply:
      * same Q/same probe continuation has only HH + interface terms;
      * Q^2 native interface work is K_phys + S after observer-gauge quotient;
      * K_phys has antisymmetric same-event donor closure and zero depth;
      * selected-family switching is a non-event Moyal boundary sidecar;
      * fresh NN service pushes to a hard-shell seed before materiality is reread;
      * source-owned service keeps the original source cause and native supplier;
        actual HH generation keeps the internal HH cause.

    MATERIAL_RELINK and NEW_COHERENT_ANCESTRY are therefore rejected as primitive
    first-stop causes.  A caller holding only such a label has unresolved PDE
    provenance and must descend further rather than manufacture recursion.
    """
    for hit in physical_hits:
        if hit.cause in NON_PRIMITIVE_MATERIAL_CAUSES:
            raise TypeError(
                f"{hit.cause.value} is not a primitive PDE generator in the native material normal form; descend to source, strain/dissipation, HH work, inherited stock, or conservative same-event relink"
            )

    suppliers = _validate_source_suppliers(physical_hits, source_supplier_kinds)

    zero_depth: set[str] = set()
    sidecars: tuple[str, ...] = ()
    if material_sidecar_relay is not None:
        if not isinstance(material_sidecar_relay, MaterialSidecarStockCentralRelayCertificate):
            raise TypeError("typed central material-sidecar relay certificate required")
        sidecars = material_sidecar_relay.sidecar_events
        zero_depth.add(MATERIAL_BOUNDARY_SIDECAR_RELAY)

    if smooth_relink_certificate is not None:
        if not isinstance(smooth_relink_certificate, SmoothRelinkDonorCertificate):
            raise TypeError("typed smooth K_phys donor certificate required")
        if (
            smooth_relink_certificate.recursive_generation_created
            or smooth_relink_certificate.new_causal_charge_created
            or smooth_relink_certificate.scale_progress_created
        ):
            raise TypeError("conservative smooth relink certificate was promoted beyond same-event provenance")
        zero_depth.add(SMOOTH_RELINK_SAME_EVENT_RELAY)

    if fresh_service_scale_route is not None:
        _validate_fresh_service_scale_route(fresh_service_scale_route)
        zero_depth.add(FRESH_SERVICE_SHELL_RELAY)

    projection: JointStopProjection | None
    if physical_hits or internal_hits or fixed_transfer_loss or kelvin_flat_certified:
        projection = joint_stop_master_projection(
            physical_hits=physical_hits,
            internal_hits=internal_hits,
            fixed_transfer_loss=fixed_transfer_loss,
            kelvin_flat_certified=kelvin_flat_certified,
            uniform_certificates=uniform_certificates,
        )
        recursive_physical = projection.joint_physical_causes
        recursive_internal = projection.joint_internal_causes
    else:
        projection = None
        recursive_physical = ()
        recursive_internal = ()

    return MaterialSourceServiceNativeNormalForm(
        projection=projection,
        recursive_physical_causes=recursive_physical,
        recursive_internal_causes=recursive_internal,
        source_supplier_kinds=tuple(x.value for x in suppliers),
        zero_depth_relays=tuple(sorted(zero_depth)),
        sidecar_events=sidecars,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "local_ns_normal_form": (
            "same-carrier material labels do not enter the coefficient equation; after Q^2 energy reentry and common-gauge quotient, native interface work has only conservative K_phys relink plus existing symmetric strain, while actual generation is source or HH work"
        ),
        "material_sidecar": (
            "membership rereading and selected-family R_switch remain non-event sidecars; R_switch is never reused as dW, stock, K_phys, or a first-stop weight"
        ),
        "conservative_relink": (
            "smooth K_phys relink has finite same-event donor closure, zero net generation, zero scale progress, and zero recursive depth"
        ),
        "fresh_service": (
            "fresh NN coherent service is pushed to its refinement-invariant frequency law and hard critical-shell seed; freshness is provenance and the next recursive owner is the generic shell first-stop law"
        ),
        "source_service": (
            "a resolved-source hit is accepted only when bound to an already-certified PDE-facing supplier kind; the material manifestation does not mint a second owner"
        ),
        "primitive_material_causes_rejected": tuple(sorted(x.value for x in NON_PRIMITIVE_MATERIAL_CAUSES)),
        "claims_mixed_owner_termination": False,
        "claims_global_regularity": False,
    }
