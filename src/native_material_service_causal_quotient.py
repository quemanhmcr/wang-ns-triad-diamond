from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from src.heat_edge_material_ownership import partition_positive_edge_measure
from src.physical_branch_compiler import CauseHit, PhysicalCause
from src.smooth_quadratic_carrier_interface import GaugeQuotientedInterfaceWork


STATUS = (
    "EXACT_NATIVE_MATERIAL_SERVICE_CAUSAL_QUOTIENT__"
    "KPHYS_MATERIAL_CROSSING_IS_SUBSET_BOUNDARY_FLUX__"
    "POSITIVE_SERVICE_OO_ON_NN_ARE_POST_SERVICE_RESTRICTIONS__"
    "RAW_MATERIAL_RELINK_AND_NEW_ANCESTRY_ARE_PREOWNER_LOCATORS"
)

RAW_MATERIAL_CAUSE_LABELS = frozenset(
    {
        PhysicalCause.MATERIAL_RELINK.value,
        PhysicalCause.NEW_COHERENT_ANCESTRY.value,
    }
)
MATERIAL_INTERFACE_SERVICE_PROVENANCE = "old_new_material_service_provenance"
MATERIAL_FRESH_SERVICE_PROVENANCE = "new_new_material_service_provenance"


def _scale(*values: float) -> float:
    return max(1.0, *(abs(float(v)) for v in values))


def require_native_service_owner_labels(labels: Iterable[str]) -> tuple[str, ...]:
    """Reject a raw material-state name where a native PDE owner is required.

    A material address exit can stop a carrier chart, but it is not itself an
    energy/service source.  Before entering the causal master it must resolve to
    an independently witnessed native owner (source, strain/dissipation, actual
    nonlinear work, or another already-certified physical supplier).
    """
    out = tuple(str(x) for x in labels if str(x))
    bad = tuple(sorted(set(out).intersection(RAW_MATERIAL_CAUSE_LABELS)))
    if bad:
        raise TypeError(
            "raw material relink/new-ancestry is a carrier/material-state locator, not an independent Navier-Stokes generation owner; "
            "resolve it first to conservative K_phys boundary flux or to the native physical source/service owner: "
            + ", ".join(bad)
        )
    return out


def require_native_service_cause_hits(hits: Sequence[CauseHit]) -> tuple[CauseHit, ...]:
    """CauseHit form of :func:`require_native_service_owner_labels`."""
    out = tuple(hits)
    require_native_service_owner_labels(h.cause.value for h in out)
    return out


@dataclass(frozen=True)
class MaterialBoundaryFluxCertificate:
    """Exact subset-divergence law for gauge-quotiented smooth K_phys work."""

    selected_roles: tuple[int, ...]
    complement_roles: tuple[int, ...]
    selected_signed_relink_work: float
    signed_boundary_flux_into_selected: float
    positive_boundary_inflow: float
    positive_boundary_outflow: float
    subset_divergence_residual: float
    pair_antisymmetry_residual: float
    row_binding_residual: float
    total_relink_work_residual: float
    native_relink_strain_split_residual: float
    same_physical_event: bool = True
    recursive_generation_created: bool = False
    physical_source_created: bool = False

    def __post_init__(self) -> None:
        if not self.selected_roles or not self.complement_roles:
            raise ValueError("a proper nonempty material-role subset is required")
        vals = (
            self.positive_boundary_inflow,
            self.positive_boundary_outflow,
            self.subset_divergence_residual,
            self.pair_antisymmetry_residual,
            self.row_binding_residual,
            self.total_relink_work_residual,
            self.native_relink_strain_split_residual,
        )
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in vals):
            raise ValueError("finite nonnegative material-boundary flux diagnostics required")
        if not math.isfinite(self.selected_signed_relink_work) or not math.isfinite(self.signed_boundary_flux_into_selected):
            raise ValueError("finite signed material-boundary flux required")
        if not self.same_physical_event or self.recursive_generation_created or self.physical_source_created:
            raise ValueError("conservative material boundary flux cannot create source or recursion depth")


def material_boundary_flux_from_smooth_relink(
    interface_work: GaugeQuotientedInterfaceWork,
    selected_roles: Sequence[int],
) -> MaterialBoundaryFluxCertificate:
    """Turn smooth material-role relink into its exact boundary-flux identity.

    For the already gauge-quotiented skew operator K_phys, the role-pair law is

        T_ab = -2 Re <eta_a u, K_phys eta_b u>,  T_ab=-T_ba,
        R_a = sum_b T_ab.

    Therefore for any proper material subset O,

        sum_{a in O} R_a = sum_{a in O, b notin O} T_ab.

    All O--O terms cancel pairwise.  The material crossing is thus a signed
    boundary flux of one same-time conservative interaction, not a new source.
    """
    if not isinstance(interface_work, GaugeQuotientedInterfaceWork):
        raise TypeError("typed gauge-quotiented smooth-interface work required")
    native = np.asarray(interface_work.signed_native_interface_atoms, dtype=float)
    relink = np.asarray(interface_work.signed_physical_relink_atoms, dtype=float)
    strain = np.asarray(interface_work.signed_existing_strain_atoms, dtype=float)
    T = np.asarray(interface_work.signed_physical_relink_pair_matrix, dtype=float)
    if (
        relink.ndim != 1
        or len(relink) < 2
        or native.shape != relink.shape
        or strain.shape != relink.shape
        or T.shape != (len(relink), len(relink))
    ):
        raise ValueError("matching nontrivial native/relink/strain role rows and bound K_phys pair matrix required")
    if np.any(~np.isfinite(native)) or np.any(~np.isfinite(relink)) or np.any(~np.isfinite(strain)) or np.any(~np.isfinite(T)):
        raise ValueError("finite native/relink/strain/K_phys pair law required")
    m = len(relink)
    selected = tuple(sorted(set(int(i) for i in selected_roles)))
    if any(i < 0 or i >= m for i in selected) or not selected or len(selected) == m:
        raise ValueError("selected material roles must form a proper nonempty subset")
    complement = tuple(i for i in range(m) if i not in set(selected))

    relink_scale = max(1.0, float(np.max(np.abs(T))), float(np.max(np.abs(relink))))
    native_scale = max(relink_scale, float(np.max(np.abs(native))), float(np.max(np.abs(strain))))
    antisym = float(np.max(np.abs(T + T.T)))
    row = float(np.max(np.abs(T.sum(axis=1) - relink)))
    total = abs(float(relink.sum()))
    native_split = float(np.max(np.abs(native - relink - strain)))
    if max(antisym, row, total) > 8e-11 * relink_scale:
        raise ValueError("supplied smooth relink law is not the exact conservative K_phys pair law")
    if native_split > 8e-11 * native_scale:
        raise ValueError("supplied smooth interface certificate lost native=relink+strain")

    cross = T[np.ix_(np.asarray(selected, dtype=int), np.asarray(complement, dtype=int))]
    boundary = float(cross.sum())
    selected_work = float(relink[list(selected)].sum())
    inflow = float(np.maximum(cross, 0.0).sum())
    outflow = float(np.maximum(-cross, 0.0).sum())
    subset_residual = abs(selected_work - boundary)
    if subset_residual > 8e-11 * relink_scale:
        raise AssertionError("K_phys material subset work failed exact boundary-flux cancellation")
    if abs(boundary - (inflow - outflow)) > 8e-11 * relink_scale:
        raise AssertionError("signed material boundary flux failed positive inflow/outflow reconstruction")

    return MaterialBoundaryFluxCertificate(
        selected_roles=selected,
        complement_roles=complement,
        selected_signed_relink_work=selected_work,
        signed_boundary_flux_into_selected=boundary,
        positive_boundary_inflow=inflow,
        positive_boundary_outflow=outflow,
        subset_divergence_residual=subset_residual,
        pair_antisymmetry_residual=antisym,
        row_binding_residual=row,
        total_relink_work_residual=total,
        native_relink_strain_split_residual=native_split,
    )


@dataclass(frozen=True)
class PositiveMaterialServiceQuotientCertificate:
    """Material ownership of an already-existing positive service law."""

    service_measure: str
    native_owner: str | None
    total_service: float
    old_old_service: float
    old_new_interface_service: float
    new_new_service: float
    partition_residual: float
    material_provenance: tuple[str, ...]
    native_owner_preserved: bool = True
    service_created_by_material_partition: bool = False
    physical_work_created_by_material_partition: bool = False
    recursive_generation_created_by_material_partition: bool = False

    def __post_init__(self) -> None:
        if not self.service_measure:
            raise ValueError("named positive service measure required")
        if self.native_owner is not None:
            if not self.native_owner:
                raise ValueError("native owner, when supplied, must be named")
            require_native_service_owner_labels((self.native_owner,))
        vals = (self.total_service, self.old_old_service, self.old_new_interface_service, self.new_new_service)
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in vals):
            raise ValueError("finite nonnegative positive-service partition required")
        if not math.isfinite(self.partition_residual):
            raise ValueError("finite service partition residual required")
        if abs(self.partition_residual) > 8e-12 * _scale(*vals):
            raise ValueError("material OO/ON/NN restrictions failed to conserve the supplied positive service law")
        expected: list[str] = []
        if self.old_new_interface_service > 0.0:
            expected.append(MATERIAL_INTERFACE_SERVICE_PROVENANCE)
        if self.new_new_service > 0.0:
            expected.append(MATERIAL_FRESH_SERVICE_PROVENANCE)
        if self.material_provenance != tuple(expected):
            raise ValueError("material service provenance does not match the exact positive restrictions")
        if (
            not self.native_owner_preserved
            or self.service_created_by_material_partition
            or self.physical_work_created_by_material_partition
            or self.recursive_generation_created_by_material_partition
        ):
            raise ValueError("post-service material restriction was promoted into a physical source/owner")


def positive_material_service_causal_quotient(
    *,
    service_measure: str,
    native_owner: str | None,
    edge_weights: Sequence[float],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
) -> PositiveMaterialServiceQuotientCertificate:
    """Restrict a positive service law by material endpoints without minting cause.

    The caller must already possess the positive physical law.  OO/ON/NN are then
    pointwise restrictions of that law.  A native causal owner is carried only if
    an upstream PDE theorem has independently certified one; this quotient never
    infers a supplier.  The restrictions may reveal interface/fresh provenance,
    capacity, or scale information, but they do not create another unit of service,
    work, or recursive event depth.
    """
    if native_owner is not None:
        require_native_service_owner_labels((native_owner,))
    p = partition_positive_edge_measure(edge_weights, old_here, old_neighbor)
    marks: list[str] = []
    if float(p["old_new_interface"]) > 0.0:
        marks.append(MATERIAL_INTERFACE_SERVICE_PROVENANCE)
    if float(p["new_new"]) > 0.0:
        marks.append(MATERIAL_FRESH_SERVICE_PROVENANCE)
    return PositiveMaterialServiceQuotientCertificate(
        service_measure=str(service_measure),
        native_owner=None if native_owner is None else str(native_owner),
        total_service=float(p["total"]),
        old_old_service=float(p["old_old"]),
        old_new_interface_service=float(p["old_new_interface"]),
        new_new_service=float(p["new_new"]),
        partition_residual=float(p["partition_residual"]),
        material_provenance=tuple(marks),
    )


def material_ownership_rereading_anti_theorem(
    *,
    service_measure: str,
    native_owner: str | None,
    edge_weights: Sequence[float],
    first_old_here: Sequence[bool],
    first_old_neighbor: Sequence[bool],
    second_old_here: Sequence[bool],
    second_old_neighbor: Sequence[bool],
) -> dict[str, object]:
    """Same positive physical law, two material readings, possibly different ON/NN.

    This is the service analogue of the selected-family anti-theorem.  The state
    and positive service weights are held fixed; only endpoint ownership is read
    differently.  Hence any change in ON/NN mass cannot by itself be a physical
    generation event.
    """
    a = positive_material_service_causal_quotient(
        service_measure=service_measure,
        native_owner=native_owner,
        edge_weights=edge_weights,
        old_here=first_old_here,
        old_neighbor=first_old_neighbor,
    )
    b = positive_material_service_causal_quotient(
        service_measure=service_measure,
        native_owner=native_owner,
        edge_weights=edge_weights,
        old_here=second_old_here,
        old_neighbor=second_old_neighbor,
    )
    tol = 8e-12 * _scale(a.total_service, b.total_service)
    if abs(a.total_service - b.total_service) > tol:
        raise AssertionError("material rereading changed the underlying positive service law")
    l1 = (
        abs(a.old_old_service - b.old_old_service)
        + abs(a.old_new_interface_service - b.old_new_interface_service)
        + abs(a.new_new_service - b.new_new_service)
    )
    return {
        "first": a,
        "second": b,
        "same_underlying_service": True,
        "same_native_owner": a.native_owner == b.native_owner,
        "ownership_partition_l1_change": l1,
        "material_generation_created": False,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "pde_subset_flux": (
            "after exact affine/Kelvin gauge quotient, K_phys is skew-adjoint; its role matrix T_ab is antisymmetric, so for every proper subset O of the certified smooth roles, "
            "sum_{a in O} R_a equals the signed O^c->O boundary flux and all O--O circulation cancels"
        ),
        "relink_semantics": "smooth K_phys material crossing is same-event conservative boundary flux: signed redistribution remains visible, but it creates no new net source/work charge, recursion depth, or scale progress",
        "positive_service_order": "the positive physical law must exist before coherent endpoints are restricted into OO/ON/NN; a native causal supplier is carried only when an upstream PDE theorem has independently certified one",
        "material_service_semantics": "ON and NN are provenance/restrictions of that already-existing positive law, not additional service, work, or causal generation; the quotient itself never infers a supplier",
        "anti_theorem": "holding the positive service weights fixed while changing only material ownership can change ON/NN strictly while total service and any independently supplied native owner stay fixed",
        "master_boundary": "raw material_relink/new_coherent_ancestry labels are unresolved carrier/material-state locators; they must resolve to conservative K_phys flux or an independently typed native physical owner before canonical recursion",
        "native_owners": "source/SGS, strain or critical dissipation, actual nonlinear work, and already-certified physical shell/service suppliers retain their own causal ownership; material provenance does not replace or duplicate them",
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class NativeMaterialServiceStress:
    samples: int
    maximum_subset_divergence_residual: float
    maximum_pair_antisymmetry_residual: float
    maximum_service_partition_residual: float
    rereadings_with_changed_material_partition: int
    raw_material_owner_rejections: int
    recursive_generation_creations: int


def stress(samples: int = 50_000, seed: int = 2026081307) -> NativeMaterialServiceStress:
    rng = np.random.default_rng(int(seed))
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    w_subset = w_skew = w_service = 0.0
    changed = rejected = recursive = 0
    for j in range(count):
        m = int(rng.integers(2, 18))
        A = rng.normal(size=(m, m))
        T = A - A.T
        relink = T.sum(axis=1)
        work = GaugeQuotientedInterfaceWork(
            signed_native_interface_atoms=tuple(float(x) for x in relink),
            signed_physical_relink_atoms=tuple(float(x) for x in relink),
            signed_existing_strain_atoms=tuple(0.0 for _ in range(m)),
            gauge_transport_operator_residual=0.0,
            skew_decomposition_residual=0.0,
            signed_physical_relink_pair_matrix=tuple(tuple(float(x) for x in row) for row in T),
        )
        mask = rng.random(m) < 0.5
        if bool(mask.all()):
            mask[0] = False
        if not bool(mask.any()):
            mask[0] = True
        flux = material_boundary_flux_from_smooth_relink(work, np.flatnonzero(mask).tolist())
        w_subset = max(w_subset, flux.subset_divergence_residual)
        w_skew = max(w_skew, flux.pair_antisymmetry_residual)
        recursive += int(flux.recursive_generation_created or flux.physical_source_created)

        n = int(rng.integers(4, 80))
        weights = rng.lognormal(mean=-0.5, sigma=1.2, size=n)
        a0 = rng.random(n) < 0.5
        a1 = rng.random(n) < 0.5
        b0 = ~a0 if j % 2 == 0 else (rng.random(n) < 0.5)
        b1 = ~a1 if j % 2 == 0 else (rng.random(n) < 0.5)
        anti = material_ownership_rereading_anti_theorem(
            service_measure="fixed_positive_law_with_pretyped_native_owner",
            native_owner=PhysicalCause.RESOLVED_SOURCE.value if j % 2 == 0 else PhysicalCause.HIGH_STRAIN_DISSIPATION.value,
            edge_weights=weights,
            first_old_here=a0,
            first_old_neighbor=a1,
            second_old_here=b0,
            second_old_neighbor=b1,
        )
        qa = anti["first"]
        qb = anti["second"]
        if not isinstance(qa, PositiveMaterialServiceQuotientCertificate) or not isinstance(qb, PositiveMaterialServiceQuotientCertificate):
            raise AssertionError("material-service anti-theorem lost typed quotient certificates")
        w_service = max(w_service, abs(qa.partition_residual), abs(qb.partition_residual))
        changed += int(float(anti["ownership_partition_l1_change"]) > 1e-12 * max(1.0, qa.total_service))
        recursive += int(qa.recursive_generation_created_by_material_partition or qb.recursive_generation_created_by_material_partition)

        if j % 7 == 0:
            for raw in RAW_MATERIAL_CAUSE_LABELS:
                try:
                    require_native_service_owner_labels((raw,))
                except TypeError:
                    rejected += 1
                else:
                    raise AssertionError("raw material locator entered native service owner boundary")

    if recursive:
        raise AssertionError("material quotient manufactured recursive generation")
    return NativeMaterialServiceStress(
        samples=count,
        maximum_subset_divergence_residual=w_subset,
        maximum_pair_antisymmetry_residual=w_skew,
        maximum_service_partition_residual=w_service,
        rereadings_with_changed_material_partition=changed,
        raw_material_owner_rejections=rejected,
        recursive_generation_creations=recursive,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081307)
    ap.add_argument("--outdir", type=Path, default=Path("results-native-material-service-causal-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples, args.seed)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "native_material_service_causal_quotient.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    md = f"""# Native material-service causal quotient

Status: **{STATUS}**.

After common affine/Kelvin gauge motion is removed, smooth `K_phys` material crossing is an exact subset boundary flux: internal role circulation cancels by antisymmetry. Separately, OO/ON/NN material ownership is read only after a positive service law already exists, so those restrictions preserve the native supplier and mint neither work nor recursion depth. Raw `material_relink` / `new_coherent_ancestry` names are therefore carrier/material-state locators until resolved to this conservative flux or to an independently witnessed native PDE owner.

Stress: `{out.samples}` states
- maximum K_phys subset-divergence residual: `{out.maximum_subset_divergence_residual:.3e}`
- maximum K_phys pair-antisymmetry residual: `{out.maximum_pair_antisymmetry_residual:.3e}`
- maximum positive-service partition residual: `{out.maximum_service_partition_residual:.3e}`
- material rereadings that changed OO/ON/NN: `{out.rereadings_with_changed_material_partition}`
- raw material-owner rejections: `{out.raw_material_owner_rejections}`
- recursive generation created by material quotient: `{out.recursive_generation_creations}`

No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
