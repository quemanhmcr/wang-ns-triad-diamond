from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
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


@dataclass(frozen=True)
class NativeNormalFormStress:
    samples: int
    naked_material_rejections: int
    naked_new_ancestry_rejections: int
    unbound_source_rejections: int
    source_supplier_routes: int
    conservative_relink_zero_depth_routes: int
    hh_routes: int
    source_strain_hh_ties: int
    fresh_service_shell_relays: int
    primitive_material_generators_created: int
    conservative_relink_recursions_created: int
    maximum_fresh_shell_mass_lower: float


def _stress_smooth_relink_certificate(scale: float) -> SmoothRelinkDonorCertificate:
    x = float(scale)
    if not math.isfinite(x) or x <= 0:
        raise ValueError("positive finite relink stress scale required")
    return SmoothRelinkDonorCertificate(
        relink_owner="smooth_physical_conservative_relink",
        recipient_roles=(0,),
        terminal_negative_net_donor_roles=(1,),
        maximum_shortest_donor_path_length=1,
        role_count=2,
        positive_relink_work=x,
        recipient_positive_incoming_flux=x,
        pair_antisymmetry_residual=0.0,
        row_binding_residual=0.0,
        total_relink_work_residual=0.0,
    )


def stress(samples: int = 50_000, seed: int = 2026081307) -> NativeNormalFormStress:
    """Randomized topology/refinement referee for the exact normal-form guards.

    This is regression evidence only.  It deliberately samples the already-proved
    routing identities; it does not replace their analytic PDE proofs.
    """
    from src.fresh_service_scale_reentry import fresh_service_scale_route as make_fresh_route

    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = random.Random(int(seed))
    material_reject = ancestry_reject = source_reject = 0
    source_routes = relink_routes = hh_routes = ties = fresh_routes = 0
    primitive = relink_recursive = 0
    max_fresh_mass = 0.0

    supplier_pool = (
        SupplierKind.RESOLVED_DISSIPATION,
        SupplierKind.PRESSURE_PAIR,
        SupplierKind.FRESH_SGS_SCALE,
        SupplierKind.HIGH_TAIL,
        SupplierKind.GENERIC_CRITICAL_SHELL,
        SupplierKind.MATERIAL_REUSE,
    )

    for j in range(count):
        mode = j % 7
        t = rng.uniform(1.0e-8, 3.0)
        weight = 10.0 ** rng.uniform(-9.0, 6.0)

        if mode == 0:
            try:
                pde_native_material_source_service_projection(
                    physical_hits=(CauseHit(t, PhysicalCause.MATERIAL_RELINK, weight),)
                )
            except TypeError:
                material_reject += 1
            else:
                primitive += 1
                raise AssertionError("naked MATERIAL_RELINK escaped the native provenance barrier")
            continue

        if mode == 1:
            try:
                pde_native_material_source_service_projection(
                    physical_hits=(CauseHit(t, PhysicalCause.NEW_COHERENT_ANCESTRY, weight),)
                )
            except TypeError:
                ancestry_reject += 1
            else:
                primitive += 1
                raise AssertionError("naked NEW_COHERENT_ANCESTRY escaped the native provenance barrier")
            continue

        if mode == 2:
            try:
                pde_native_material_source_service_projection(
                    physical_hits=(CauseHit(t, PhysicalCause.RESOLVED_SOURCE, weight),)
                )
            except TypeError:
                source_reject += 1
            else:
                raise AssertionError("unbound resolved source escaped the native supplier barrier")
            continue

        if mode == 3:
            supplier = supplier_pool[j % len(supplier_pool)]
            out = pde_native_material_source_service_projection(
                physical_hits=(CauseHit(t, PhysicalCause.RESOLVED_SOURCE, weight, "native source"),),
                source_supplier_kinds=(supplier,),
            )
            if out.recursive_physical_causes != (PhysicalCause.RESOLVED_SOURCE.value,):
                raise AssertionError("native source supplier changed the independently witnessed source owner")
            source_routes += 1

        elif mode == 4:
            cert = _stress_smooth_relink_certificate(weight)
            out = pde_native_material_source_service_projection(smooth_relink_certificate=cert)
            if out.projection is not None or out.recursive_physical_causes or out.recursive_internal_causes:
                relink_recursive += 1
                raise AssertionError("conservative K_phys relink manufactured recursive depth")
            if SMOOTH_RELINK_SAME_EVENT_RELAY not in out.zero_depth_relays:
                raise AssertionError("conservative K_phys donor provenance was lost")
            relink_routes += 1

        elif mode == 5:
            out = pde_native_material_source_service_projection(
                internal_hits=(InternalHit(t),)
            )
            if not out.recursive_internal_causes:
                raise AssertionError("actual HH regeneration owner was lost")
            hh_routes += 1

        else:
            supplier = supplier_pool[j % len(supplier_pool)]
            out = pde_native_material_source_service_projection(
                physical_hits=(
                    CauseHit(t, PhysicalCause.RESOLVED_SOURCE, weight, "source"),
                    CauseHit(t, PhysicalCause.HIGH_STRAIN_DISSIPATION, weight * rng.uniform(0.2, 5.0), "strain"),
                ),
                internal_hits=(InternalHit(t),),
                source_supplier_kinds=(supplier,),
            )
            if set(out.recursive_physical_causes) != {
                PhysicalCause.RESOLVED_SOURCE.value,
                PhysicalCause.HIGH_STRAIN_DISSIPATION.value,
            } or not out.recursive_internal_causes:
                raise AssertionError("exact heterogeneous physical tie was split or prioritized")
            ties += 1

        # Independently exercise the fresh-service supplier relay on a sparse
        # deterministic subsequence so the stress covers the real scale theorem,
        # not a forged dictionary fixture.
        if j % 23 == 0:
            Y = 10.0 ** rng.uniform(-6.0, 2.0)
            c = 10.0 ** rng.uniform(-2.0, 0.5)
            N = 10.0 ** rng.uniform(-1.0, 4.0)
            fresh = rng.uniform(0.25, 0.95) * Y
            route = make_fresh_route(Y, c, N, {0: fresh})
            fresh_out = pde_native_material_source_service_projection(fresh_service_scale_route=route)
            if fresh_out.projection is not None or fresh_out.recursive_physical_causes:
                raise AssertionError("fresh service supplier relay minted a material first stop")
            if FRESH_SERVICE_SHELL_RELAY not in fresh_out.zero_depth_relays:
                raise AssertionError("fresh service lost its hard-shell supplier relay")
            max_fresh_mass = max(max_fresh_mass, float(route["hard_shell_mass_lower"]))
            fresh_routes += 1

        primitive += int(out.primitive_material_generator_created)
        relink_recursive += int(out.conservative_relink_promoted_to_recursion)

    if primitive or relink_recursive:
        raise AssertionError("native material normal form created forbidden recursive topology")

    return NativeNormalFormStress(
        samples=count,
        naked_material_rejections=material_reject,
        naked_new_ancestry_rejections=ancestry_reject,
        unbound_source_rejections=source_reject,
        source_supplier_routes=source_routes,
        conservative_relink_zero_depth_routes=relink_routes,
        hh_routes=hh_routes,
        source_strain_hh_ties=ties,
        fresh_service_shell_relays=fresh_routes,
        primitive_material_generators_created=primitive,
        conservative_relink_recursions_created=relink_recursive,
        maximum_fresh_shell_mass_lower=max_fresh_mass,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081307)
    ap.add_argument("--outdir", type=Path, default=Path("results-material-source-service-native-normal-form"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = stress(args.samples, args.seed)
    payload = {"certificate": theorem_certificate(), "stress": asdict(result)}
    (args.outdir / "material_source_service_native_normal_form.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    summary = f"""# PDE-native material/source-service normal form

Status: **{STATUS}**.

This is a draft composition/refinement theorem.  It does not claim mixed-owner termination or Navier--Stokes regularity.

Stress: `{result.samples}` routed/rejected topology states
- naked MATERIAL_RELINK rejections: `{result.naked_material_rejections}`
- naked NEW_COHERENT_ANCESTRY rejections: `{result.naked_new_ancestry_rejections}`
- unbound source rejections: `{result.unbound_source_rejections}`
- source + native supplier routes: `{result.source_supplier_routes}`
- conservative K_phys zero-depth routes: `{result.conservative_relink_zero_depth_routes}`
- actual HH routes: `{result.hh_routes}`
- exact source/strain/HH ties: `{result.source_strain_hh_ties}`
- fresh-service hard-shell relays: `{result.fresh_service_shell_relays}`
- primitive material generators created: `{result.primitive_material_generators_created}`
- conservative relinks promoted to recursion: `{result.conservative_relink_recursions_created}`
- maximum sampled certified fresh-shell mass lower: `{result.maximum_fresh_shell_mass_lower:.12e}`
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
