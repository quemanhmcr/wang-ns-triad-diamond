from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Sequence

from src.heat_edge_material_ownership import partition_positive_edge_measure
from src.material_sidecar_stock_owner_decomposition import (
    MEMBERSHIP_PROVENANCE_CURRENCY,
    SELECTED_FAMILY_MOYAL_CURRENCY,
)
from src.physical_branch_compiler import PhysicalCause
from src.smooth_quadratic_carrier_interface import RELINK_OWNER
from src.smooth_relink_donor_quotient import (
    SMOOTH_RELINK_SAME_EVENT_RELAY,
    SmoothRelinkDonorCertificate,
)


STATUS = (
    "DRAFT_EXACT_POSITIVE_SERVICE_MATERIAL_OWNER_FACTORIZATION__"
    "MATERIAL_PARTITION_COMMUTES_WITH_NATIVE_OWNER_RESTRICTION__"
    "NO_STANDALONE_MATERIAL_OR_NEW_ANCESTRY_CAUSAL_ROOT__"
    "SMOOTH_KPHYS_RELINK_REMAINS_SAME_EVENT_PROVENANCE"
)


class MaterialEdgeClass(str, Enum):
    OO = "old_old"
    ON = "old_new_interface"
    NN = "new_new"


# These names may survive as provenance/legacy API labels, but this theorem does
# not admit them as the native owner of a positive service law.  Doing so would
# be circular: material classification is performed only after the positive law
# and its physical owner already exist.
FORBIDDEN_NATIVE_OWNER_ROOTS = frozenset(
    {
        PhysicalCause.MATERIAL_RELINK.value,
        PhysicalCause.NEW_COHERENT_ANCESTRY.value,
        RELINK_OWNER,
        MEMBERSHIP_PROVENANCE_CURRENCY,
        SELECTED_FAMILY_MOYAL_CURRENCY,
    }
)


@dataclass(frozen=True)
class NativePositiveServiceAtom:
    """One already-positive physical service atom with its pre-material owner."""

    weight: float
    native_owner: str
    old_here: bool
    old_neighbor: bool

    def __post_init__(self) -> None:
        w = float(self.weight)
        if not math.isfinite(w) or w < 0.0:
            raise ValueError("finite nonnegative positive-service atom weight required")
        if not self.native_owner:
            raise ValueError("every positive service atom requires a native physical owner")
        if self.native_owner in FORBIDDEN_NATIVE_OWNER_ROOTS:
            raise ValueError(
                "material/relink bookkeeping cannot be the native root of an already-positive service atom"
            )

    @property
    def material_class(self) -> MaterialEdgeClass:
        if self.old_here and self.old_neighbor:
            return MaterialEdgeClass.OO
        if self.old_here or self.old_neighbor:
            return MaterialEdgeClass.ON
        return MaterialEdgeClass.NN


def _sorted_mass_map(items: Iterable[tuple[str, float]]) -> tuple[tuple[str, float], ...]:
    acc: dict[str, float] = {}
    for key, value in items:
        acc[key] = acc.get(key, 0.0) + float(value)
    return tuple(sorted((key, value) for key, value in acc.items() if value > 0.0))


def _relative_residual(a: float, b: float) -> float:
    return abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), 1.0)


@dataclass(frozen=True)
class MaterialServiceOwnerFactorizationCertificate:
    """Exact finite-atom model of the Radon restriction/distributivity theorem."""

    total_service: float
    oo_service: float
    on_service: float
    nn_service: float
    owner_total: tuple[tuple[str, float], ...]
    owner_oo: tuple[tuple[str, float], ...]
    owner_on: tuple[tuple[str, float], ...]
    owner_nn: tuple[tuple[str, float], ...]
    material_native_owner_provenance: tuple[str, ...]
    ownership_partition_residual: float
    owner_disintegration_residual: float
    recursive_material_owner_created: bool = False
    recursive_new_ancestry_owner_created: bool = False
    new_recursive_vertex_created: bool = False
    owner_cloning_used: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        scalars = (
            self.total_service,
            self.oo_service,
            self.on_service,
            self.nn_service,
            self.ownership_partition_residual,
            self.owner_disintegration_residual,
        )
        if not all(math.isfinite(x) and x >= 0.0 for x in scalars):
            raise ValueError("finite nonnegative factorization diagnostics required")
        if self.ownership_partition_residual > 5.0e-13:
            raise ValueError("OO/ON/NN restrictions failed to partition the positive service law")
        if self.owner_disintegration_residual > 5.0e-13:
            raise ValueError("material restriction failed to commute with native-owner disintegration")
        for table in (self.owner_total, self.owner_oo, self.owner_on, self.owner_nn):
            if tuple(sorted(table)) != table:
                raise ValueError("owner mass tables must be sorted")
            if any(not owner or owner in FORBIDDEN_NATIVE_OWNER_ROOTS for owner, _ in table):
                raise ValueError("material bookkeeping appeared as a native owner root")
            if any(not math.isfinite(mass) or mass <= 0.0 for _, mass in table):
                raise ValueError("owner mass tables contain nonpositive/nonfinite mass")
        expected = tuple(sorted({owner for owner, _ in (*self.owner_on, *self.owner_nn)}))
        if self.material_native_owner_provenance != expected:
            raise ValueError("material-bearing service did not inherit exactly its pre-existing native owners")
        if (
            self.recursive_material_owner_created
            or self.recursive_new_ancestry_owner_created
            or self.new_recursive_vertex_created
            or self.owner_cloning_used
            or self.later_hahn_used
        ):
            raise ValueError("material restriction cannot mint event depth, clone charge, or Hahn again")


def factor_positive_service_by_material(
    atoms: Sequence[NativePositiveServiceAtom],
) -> MaterialServiceOwnerFactorizationCertificate:
    """Restrict one positive service law by OO/ON/NN after native ownership.

    Measure-level statement.  If ``sigma=sum_r sigma_r`` is an already-positive
    owner decomposition and ``chi_C`` is one of the nonnegative material
    indicators, then

        sigma_C = chi_C sigma = sum_r chi_C sigma_r.

    Thus ON/NN may mark interface/fresh provenance, but they cannot create a new
    causal root or event vertex.  The finite-atom implementation is the exact
    discrete algebra used for regression/adversarial checks; the Radon proof is
    pointwise distributivity of measurable restriction over a positive sum.
    """
    if not atoms:
        raise ValueError("at least one already-positive service atom required")
    if not all(isinstance(atom, NativePositiveServiceAtom) for atom in atoms):
        raise TypeError("typed NativePositiveServiceAtom inputs required")

    total = float(sum(atom.weight for atom in atoms))
    if total <= 0.0:
        raise ValueError("positive total physical service required")

    # Bind the candidate factorization to the already-certified native material
    # partition rather than reimplementing OO/ON/NN arithmetic locally.
    certified_partition = partition_positive_edge_measure(
        [atom.weight for atom in atoms],
        [atom.old_here for atom in atoms],
        [atom.old_neighbor for atom in atoms],
    )
    oo_mass = float(certified_partition["old_old"])
    on_mass = float(certified_partition["old_new_interface"])
    nn_mass = float(certified_partition["new_new"])
    partition_residual = _relative_residual(
        total, oo_mass + on_mass + nn_mass
    )

    owner_total = _sorted_mass_map((atom.native_owner, atom.weight) for atom in atoms if atom.weight > 0.0)
    owner_oo_items: list[tuple[str, float]] = []
    owner_on_items: list[tuple[str, float]] = []
    owner_nn_items: list[tuple[str, float]] = []
    for owner, _ in owner_total:
        owner_atoms = tuple(atom for atom in atoms if atom.native_owner == owner)
        part = partition_positive_edge_measure(
            [atom.weight for atom in owner_atoms],
            [atom.old_here for atom in owner_atoms],
            [atom.old_neighbor for atom in owner_atoms],
        )
        if part["old_old"] > 0.0:
            owner_oo_items.append((owner, float(part["old_old"])))
        if part["old_new_interface"] > 0.0:
            owner_on_items.append((owner, float(part["old_new_interface"])))
        if part["new_new"] > 0.0:
            owner_nn_items.append((owner, float(part["new_new"])))
    owner_oo = tuple(sorted(owner_oo_items))
    owner_on = tuple(sorted(owner_on_items))
    owner_nn = tuple(sorted(owner_nn_items))

    by_material_then_owner: dict[str, float] = {}
    for table in (owner_oo, owner_on, owner_nn):
        for owner, mass in table:
            by_material_then_owner[owner] = by_material_then_owner.get(owner, 0.0) + mass
    owner_disintegration_residual = max(
        (
            _relative_residual(mass, by_material_then_owner.get(owner, 0.0))
            for owner, mass in owner_total
        ),
        default=0.0,
    )

    material_owner_provenance = tuple(sorted({owner for owner, _ in (*owner_on, *owner_nn)}))
    return MaterialServiceOwnerFactorizationCertificate(
        total_service=total,
        oo_service=oo_mass,
        on_service=on_mass,
        nn_service=nn_mass,
        owner_total=owner_total,
        owner_oo=owner_oo,
        owner_on=owner_on,
        owner_nn=owner_nn,
        material_native_owner_provenance=material_owner_provenance,
        ownership_partition_residual=partition_residual,
        owner_disintegration_residual=owner_disintegration_residual,
    )


@dataclass(frozen=True)
class MaterialRecurrenceProjection:
    """Projection of material observations back to their already-existing roots."""

    native_owner_provenance: tuple[str, ...]
    same_event_relays: tuple[str, ...]
    sidecar_currencies: tuple[str, ...]
    selected_family_boundary_energy: float
    service_factorization: MaterialServiceOwnerFactorizationCertificate | None
    recursive_material_owner_created: bool = False
    recursive_new_ancestry_owner_created: bool = False
    new_recursive_vertex_created: bool = False

    def __post_init__(self) -> None:
        if tuple(sorted(set(self.native_owner_provenance))) != self.native_owner_provenance:
            raise ValueError("native owner provenance must be a sorted quotient")
        if any(owner in FORBIDDEN_NATIVE_OWNER_ROOTS for owner in self.native_owner_provenance):
            raise ValueError("material/new-ancestry bookkeeping survived as a primitive owner")
        if tuple(sorted(set(self.same_event_relays))) != self.same_event_relays:
            raise ValueError("same-event relays must be a sorted quotient")
        if tuple(sorted(set(self.sidecar_currencies))) != self.sidecar_currencies:
            raise ValueError("sidecar currencies must be a sorted quotient")
        boundary = float(self.selected_family_boundary_energy)
        if not math.isfinite(boundary) or boundary < 0.0:
            raise ValueError("finite nonnegative selected-family boundary energy required")
        if boundary > 0.0 and SELECTED_FAMILY_MOYAL_CURRENCY not in self.sidecar_currencies:
            raise ValueError("positive selected-family boundary energy lacks its Moyal sidecar currency")
        if (
            self.recursive_material_owner_created
            or self.recursive_new_ancestry_owner_created
            or self.new_recursive_vertex_created
        ):
            raise ValueError("material observation cannot mint a primitive causal root or event vertex")


def project_material_recurrence_to_native_owners(
    *,
    service_atoms: Sequence[NativePositiveServiceAtom] = (),
    membership_reread: bool = False,
    selected_family_boundary_energy: float = 0.0,
    smooth_relink: SmoothRelinkDonorCertificate | None = None,
    additional_native_owner_provenance: Sequence[str] = (),
    require_physical_role_change_owner: bool = False,
) -> MaterialRecurrenceProjection:
    """Remove material bookkeeping as a primitive recursive owner, fail-closed.

    This function never says that inherited provenance itself creates another
    event.  It only records which already-existing native owner remains attached
    after material restriction.
    """
    boundary = float(selected_family_boundary_energy)
    if not math.isfinite(boundary) or boundary < 0.0:
        raise ValueError("finite nonnegative selected-family boundary energy required")

    sidecars: set[str] = set()
    if membership_reread:
        sidecars.add(MEMBERSHIP_PROVENANCE_CURRENCY)
    if boundary > 0.0:
        sidecars.add(SELECTED_FAMILY_MOYAL_CURRENCY)

    relays: set[str] = set()
    if smooth_relink is not None:
        if not isinstance(smooth_relink, SmoothRelinkDonorCertificate):
            raise TypeError("typed smooth K_phys donor certificate required")
        if smooth_relink.recursive_generation_created or smooth_relink.new_causal_charge_created:
            raise ValueError("smooth relink certificate illegally created recursive generation")
        relays.add(SMOOTH_RELINK_SAME_EVENT_RELAY)

    factorization = factor_positive_service_by_material(service_atoms) if service_atoms else None
    native: set[str] = set(
        factorization.material_native_owner_provenance if factorization is not None else ()
    )
    for owner in additional_native_owner_provenance:
        if not owner:
            raise ValueError("empty additional native owner provenance")
        if owner in FORBIDDEN_NATIVE_OWNER_ROOTS:
            raise ValueError("material/new-ancestry bookkeeping cannot be supplied as an independent native owner")
        native.add(str(owner))

    if require_physical_role_change_owner and not native and smooth_relink is None:
        raise TypeError(
            "genuine role/probe change is fail-closed until native physical work/source/strain ownership or a bound smooth relink law is supplied"
        )

    return MaterialRecurrenceProjection(
        native_owner_provenance=tuple(sorted(native)),
        same_event_relays=tuple(sorted(relays)),
        sidecar_currencies=tuple(sorted(sidecars)),
        selected_family_boundary_energy=boundary,
        service_factorization=factorization,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "measure_identity": "for an already-positive owner law sigma=sum_r sigma_r and material class C, chi_C sigma=sum_r chi_C sigma_r; executable OO/ON/NN restriction is delegated to the certified heat_edge_material_ownership partition",
        "hahn_rule": "positivity is established before material classification; OO/ON/NN is positive restriction only and never a later Hahn split",
        "membership": "same-carrier material membership rereading is zero-charge provenance and creates no service/work/event vertex",
        "selected_family": "R_switch remains selected-family Moyal boundary currency with zero generation depth; it is not dW, stock, K_phys, or a causal stop",
        "freshness": "NN/fresh provenance is a restriction of existing positive service and inherits that service law's native owner provenance; NEW_COHERENT_ANCESTRY is not minted as a second root",
        "interface": "ON/interface provenance is a restriction of existing positive service and inherits the same native owner provenance; MATERIAL_RELINK is not minted as a second root",
        "smooth_role_relink": "certified gauge-quotiented K_phys relink remains finite same-event donor provenance; only independently registered native owners survive as event roots",
        "fail_closed": "a genuine role/probe change without independently registered physical work/source/strain ownership or a bound K_phys relay is rejected",
        "master_consequence": "material/new-ancestry labels do not form an independent letter in the recursive owner word once their native positive law is registered",
        "scope": "draft local factorization/reduction theorem only; source-HH-strain mixed recurrence, generic HH termination, initial-data and singular-time interfaces remain open; no global-regularity claim",
    }
