from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Sequence

from src.joint_causal_stop_projection import InternalRecursiveCause
from src.material_service_native_owner_factorization import MaterialRecurrenceProjection
from src.objective_source_routing_compiler import compile_objective_source_owners
from src.physical_branch_compiler import PhysicalCause
from src.resolved_contact_native_binding import (
    HH_WINDOW,
    K_RELAY,
    S_STRAIN,
    hard_tail_resolved_contact_route,
)


STATUS = (
    "DRAFT_MIXED_GENUINE_OWNER_NATIVE_NORMAL_FORM__"
    "MATERIAL_AND_RESOLVED_CONTACT_VOCABULARY_QUOTIENTED__"
    "ONLY_NATIVE_EVENT_ROOTS_SURVIVE__NO_TERMINATION_OVERCLAIM"
)


class NativeEventRoot(str, Enum):
    """Event roots which remain after currently proved zero-depth quotients."""

    TRANSFER_WORK_LOSS = "transfer_work_loss"
    RESOLVED_SOURCE = "resolved_source_or_independent_service"
    STRAIN_DISSIPATION = "strain_or_critical_dissipation"
    ACTUAL_HH_OR_HARD_TAIL = "actual_hh_or_hard_tail_work"
    CAUSAL_REUSE = "causal_reuse_cost"
    INTRINSIC_SIDEBAND = "intrinsic_sideband_cost"
    UNIFORM_GLOBAL_RESOURCE = "uniform_global_resource"
    INITIAL_BOUNDARY = "initial_boundary"


TERMINAL_ROOTS = frozenset(
    {
        NativeEventRoot.TRANSFER_WORK_LOSS,
        NativeEventRoot.CAUSAL_REUSE,
        NativeEventRoot.INTRINSIC_SIDEBAND,
        NativeEventRoot.UNIFORM_GLOBAL_RESOURCE,
        NativeEventRoot.INITIAL_BOUNDARY,
    }
)

RECURSIVE_CORE_ROOTS = frozenset(
    {
        NativeEventRoot.RESOLVED_SOURCE,
        NativeEventRoot.STRAIN_DISSIPATION,
        NativeEventRoot.ACTUAL_HH_OR_HARD_TAIL,
    }
)


PHYSICAL_CAUSE_NORMAL_FORM = {
    PhysicalCause.TRANSFER_WORK_LOSS: NativeEventRoot.TRANSFER_WORK_LOSS,
    PhysicalCause.RESOLVED_SOURCE: NativeEventRoot.RESOLVED_SOURCE,
    PhysicalCause.HIGH_STRAIN_DISSIPATION: NativeEventRoot.STRAIN_DISSIPATION,
    PhysicalCause.CAUSAL_REUSE: NativeEventRoot.CAUSAL_REUSE,
    PhysicalCause.INTRINSIC_SIDEBAND: NativeEventRoot.INTRINSIC_SIDEBAND,
    PhysicalCause.UNIFORM_GLOBAL_RESOURCE: NativeEventRoot.UNIFORM_GLOBAL_RESOURCE,
    PhysicalCause.INITIAL_BOUNDARY: NativeEventRoot.INITIAL_BOUNDARY,
}

FORBIDDEN_LEGACY_PRIMITIVES = frozenset(
    {PhysicalCause.MATERIAL_RELINK, PhysicalCause.NEW_COHERENT_ANCESTRY}
)


@dataclass(frozen=True)
class NativeOwnerNormalForm:
    event_roots: tuple[str, ...]
    zero_depth_relays: tuple[str, ...]
    material_sidecars: tuple[str, ...]
    material_native_owner_provenance: tuple[str, ...]
    recursive_core_roots: tuple[str, ...]
    terminal_roots: tuple[str, ...]
    material_primitive_survived: bool = False
    resolved_contact_primitive_survived: bool = False
    descendant_created_event_depth: bool = False

    def __post_init__(self) -> None:
        for values in (
            self.event_roots,
            self.zero_depth_relays,
            self.material_sidecars,
            self.material_native_owner_provenance,
            self.recursive_core_roots,
            self.terminal_roots,
        ):
            if tuple(sorted(set(values))) != values:
                raise ValueError("normal-form labels must be unique sorted tuples")
        if self.material_primitive_survived or self.resolved_contact_primitive_survived:
            raise ValueError("a proved intermediate vocabulary item survived as a primitive root")
        if self.descendant_created_event_depth:
            raise ValueError("a same-event consequence/relay created theorem-depth recursion")
        roots = set(self.event_roots)
        if set(self.recursive_core_roots) != roots & {r.value for r in RECURSIVE_CORE_ROOTS}:
            raise ValueError("recursive core projection disagrees with event roots")
        if set(self.terminal_roots) != roots & {r.value for r in TERMINAL_ROOTS}:
            raise ValueError("terminal projection disagrees with event roots")


def _sorted(values: set[str]) -> tuple[str, ...]:
    return tuple(sorted(values))


def native_owner_normal_form(
    *,
    physical_causes: Sequence[PhysicalCause] = (),
    internal_hh_regeneration: bool = False,
    material_projection: MaterialRecurrenceProjection | None = None,
    resolved_contact_route: dict[str, object] | None = None,
) -> NativeOwnerNormalForm:
    """Collapse only theorem manifestations already proved to be zero-depth.

    This function is intentionally conservative.  It does not turn a source into
    strain/HH merely because downstream source theorems offer such alternatives;
    doing so would require the appropriate typed source certificate at that
    physical event.  Likewise it does not infer signed-good HH geometry from a
    generic HH/hard-tail root.
    """
    roots: set[str] = set()
    relays: set[str] = set()
    sidecars: set[str] = set()
    provenance: set[str] = set()

    for cause in physical_causes:
        if not isinstance(cause, PhysicalCause):
            raise TypeError("typed PhysicalCause inputs required")
        if cause in FORBIDDEN_LEGACY_PRIMITIVES:
            raise TypeError(
                f"legacy {cause.value} cannot survive the native normal form; attach material provenance to its independently witnessed native event"
            )
        roots.add(PHYSICAL_CAUSE_NORMAL_FORM[cause].value)

    if internal_hh_regeneration:
        # The generic internal label remains recursive only after its physical
        # energy/work route has been registered.  This flag means that binding
        # has already happened upstream and the root is actual HH/hard-tail work.
        roots.add(NativeEventRoot.ACTUAL_HH_OR_HARD_TAIL.value)

    if material_projection is not None:
        if not isinstance(material_projection, MaterialRecurrenceProjection):
            raise TypeError("typed MaterialRecurrenceProjection required")
        if material_projection.new_recursive_vertex_created:
            raise ValueError("material projection illegally created event depth")
        relays.update(material_projection.same_event_relays)
        sidecars.update(material_projection.sidecar_currencies)
        provenance.update(material_projection.native_owner_provenance)
        # Native owner provenance is metadata on an already-existing event.  It
        # does not mint a replacement event root here.

    if resolved_contact_route is not None:
        continuations = resolved_contact_route.get("joint_physical_continuations")
        if not isinstance(continuations, tuple) or not continuations:
            raise TypeError("typed resolved-contact route output with joint_physical_continuations required")
        if resolved_contact_route.get("canonical_cause_replaced") is not False:
            raise ValueError("resolved-contact route replaced the canonical cause")
        if resolved_contact_route.get("later_hahn_used") is not False:
            raise ValueError("resolved-contact route used a later Hahn split")
        if resolved_contact_route.get("recipient_shell_reweighting_used") is not False:
            raise ValueError("resolved-contact route changed the common N dW causal unit")
        known = {HH_WINDOW, K_RELAY, S_STRAIN}
        if any(x not in known for x in continuations):
            raise ValueError("unknown resolved-contact continuation")
        if HH_WINDOW in continuations:
            roots.add(NativeEventRoot.ACTUAL_HH_OR_HARD_TAIL.value)
        if S_STRAIN in continuations:
            roots.add(NativeEventRoot.STRAIN_DISSIPATION.value)
        if K_RELAY in continuations:
            relays.add(K_RELAY)

    recursive = roots & {r.value for r in RECURSIVE_CORE_ROOTS}
    terminal = roots & {r.value for r in TERMINAL_ROOTS}
    return NativeOwnerNormalForm(
        event_roots=_sorted(roots),
        zero_depth_relays=_sorted(relays),
        material_sidecars=_sorted(sidecars),
        material_native_owner_provenance=_sorted(provenance),
        recursive_core_roots=_sorted(recursive),
        terminal_roots=_sorted(terminal),
    )


def certified_resolved_contact_normal_form(
    *,
    physical_tail_dissipation: float,
    viscosity: float,
    actual_contact_common_work: float,
    hh_complement_common_work: float,
    mixed_common_work: float,
    positive_skew_common_work: float,
    positive_strain_common_work: float,
    material_projection: MaterialRecurrenceProjection | None = None,
) -> NativeOwnerNormalForm:
    """Call the certified resolved-contact owner cover, then quotient its vocabulary."""
    route = hard_tail_resolved_contact_route(
        physical_tail_dissipation=physical_tail_dissipation,
        viscosity=viscosity,
        actual_contact_common_work=actual_contact_common_work,
        hh_complement_common_work=hh_complement_common_work,
        mixed_common_work=mixed_common_work,
        positive_skew_common_work=positive_skew_common_work,
        positive_strain_common_work=positive_strain_common_work,
    )
    return native_owner_normal_form(
        material_projection=material_projection,
        resolved_contact_route=route,
    )


def certified_objective_source_normal_form(
    *,
    objective_variation_action: float,
    scaled_lifetime: float,
    owner_weights: dict[str, float],
    viscosity: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    material_projection: MaterialRecurrenceProjection | None = None,
) -> NativeOwnerNormalForm:
    """Refine one coarse objective-source stop by its certified physical owners.

    The source compiler proves a non-lexicographic four-owner cover.  Local
    quadratic and viscous owners already pass through resolved D_V and therefore
    enter the existing strain/critical-dissipation class.  Pressure and SGS
    retain source/service ownership unless a later typed realized-route theorem
    resolves them further.  Exact owner ties remain joint.

    This function changes the *owner normal form* of the same source first stop;
    it does not create another event vertex for the compiler or critical shell.
    """
    compiled = compile_objective_source_owners(
        objective_variation_action,
        scaled_lifetime,
        owner_weights,
        viscosity=viscosity,
        filter_l1=filter_l1,
        lp_constant=lp_constant,
        bernstein_constant=bernstein_constant,
    )
    if compiled.get("additive_reset_created") is not False:
        raise ValueError("objective source compiler created an illicit additive reset")
    if compiled.get("packet_synchronization_created") is not False:
        raise ValueError("objective source compiler created an illicit packet synchronization interface")
    owners = compiled.get("joint_owners")
    if not isinstance(owners, tuple) or not owners:
        raise TypeError("certified source compiler returned no joint physical owners")
    known = {"local_dv", "pressure", "sgs", "viscous"}
    if any(owner not in known for owner in owners):
        raise ValueError("unknown certified objective-source owner")

    causes: list[PhysicalCause] = []
    # local_dv/viscous are already bound by the compiler to D_V -> critical shell.
    if any(owner in {"local_dv", "viscous"} for owner in owners):
        causes.append(PhysicalCause.HIGH_STRAIN_DISSIPATION)
    # pressure/SGS remain source/service roots until a realized positive-law route
    # supplies a stronger typed descendant.
    if any(owner in {"pressure", "sgs"} for owner in owners):
        causes.append(PhysicalCause.RESOLVED_SOURCE)

    return native_owner_normal_form(
        physical_causes=tuple(causes),
        material_projection=material_projection,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "proved_quotients": "material membership/family/OO-ON-NN and smooth/fixed-event K relays add provenance but no causal vertex; certified resolved contact reduces to K relay, existing strain, or actual HH natural-window continuation",
        "legacy_material": "MATERIAL_RELINK and NEW_COHERENT_ANCESTRY fail closed as primitive roots in the normal form",
        "resolved_contact": "resolved contact fails closed as a primitive root; the certified u=V+h signed-before-Hahn theorem supplies only existing K/S/HH continuations",
        "recursive_core": tuple(sorted(r.value for r in RECURSIVE_CORE_ROOTS)),
        "terminal_roots": tuple(sorted(r.value for r in TERMINAL_ROOTS)),
        "source_rule": "typed objective-source ownership refines local_DV/viscous owners into the existing strain/dissipation root; pressure/SGS remain source/service until a realized positive-law certificate resolves them further; no untyped compiler-name reduction is used",
        "hh_rule": "generic actual HH/hard-tail remains one native recursive root; signed-good and high-tail subroutes retain their own certified supplier-specific telescopes/continuations and are not conflated",
        "scope": "normal-form/ontology reduction only; the surviving mixed three-root recursive core is not yet proved finite, and no Navier-Stokes global-regularity claim is made",
    }
