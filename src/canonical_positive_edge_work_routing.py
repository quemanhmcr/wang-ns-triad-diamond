from __future__ import annotations

import argparse
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.coherent_service_or_flat import FixedTransferLossGate, fixed_transfer_loss_gate
from src.continuum_helical_edge_measure_registration import (
    ContinuumEdgeMeasureLedger,
    ContinuumHelicalEdgeIdentity,
    ContinuumModalEdgeAtom,
    HelicalModeIdentity,
    _replay_physical_ledger,
    continuum_edge_measure_ledger,
    register_continuum_triad_fiber,
)
from src.helical import coupling_g, helical_basis
from src.joint_causal_stop_projection import JointStopProjection, joint_stop_master_projection
from src.physical_branch_compiler import (
    BlockWitness,
    MasterDisposition,
    PhysicalCause,
    PhysicalCurrency,
    TransferPartition,
    compile_transfer_measure,
)
from src.service_or_flat_rigidity import CURVATURE_DENOM
from src.single_edge_certificate import float_jstar
from src.triad_extremizer import symmetric_gamma, symmetric_rstar

STATUS = (
    "EXACT_CANONICAL_POSITIVE_EDGE_WORK_ROUTING__ACTUAL_DW_HAHN_RESTRICTION__"
    "GEOMETRY_BAD_SAME_LAW_STAGE_ZERO_TRANSFER_LOSS__GEOMETRY_GOOD_ONLY_YOUNG_ELIGIBLE__"
    "DETERMINISTIC_HARD_CELL_SIGNED_COMPRESSION__POSITIVE_PUSHFORWARD_NOT_REHAHN"
)

ETA0 = 1.0e-4
MAX_CERTIFIED_TAU = 0.1


def _finite_sum(values: Sequence[float], name: str) -> float:
    vals = tuple(float(v) for v in values)
    if not all(math.isfinite(v) for v in vals):
        raise ValueError(f"{name} terms must be finite")
    try:
        out = math.fsum(vals)
    except OverflowError as exc:
        raise ValueError(f"{name} left the finite physical range") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} left the finite physical range")
    return out


def _native_close(actual: float, expected: float, *, factor: float = 2.0e-10) -> bool:
    a = float(actual)
    b = float(expected)
    if not math.isfinite(a) or not math.isfinite(b):
        return False
    return abs(a - b) <= factor * max(abs(a), abs(b), 1.0e-300)


@dataclass(frozen=True, order=True)
class EdgeOccurrenceKey:
    """One occurrence in the registered quadrature of the physical edge Radon law."""

    fiber_index: int
    atom_index: int

    def __post_init__(self) -> None:
        if self.fiber_index < 0 or self.atom_index < 0:
            raise ValueError("nonnegative physical edge occurrence indices required")


@dataclass(frozen=True)
class CanonicalPositiveEdgeOccurrence:
    key: EdgeOccurrenceKey
    physical_edge_identity: ContinuumHelicalEdgeIdentity
    positive_work_mass: float
    capacity_mass: float
    signed_progress_mass: float
    signed_efficiency: float
    scale_progress: float
    geometry_good: bool

    def __post_init__(self) -> None:
        vals = (
            self.positive_work_mass,
            self.capacity_mass,
            self.signed_progress_mass,
            self.signed_efficiency,
            self.scale_progress,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite canonical edge occurrence data required")
        if self.positive_work_mass <= 0.0:
            raise ValueError("canonical positive occurrence must carry strictly positive dW+ mass")
        if self.capacity_mass <= 0.0:
            raise ValueError("positive physical work must have positive native capacity")
        if self.positive_work_mass > self.capacity_mass * (1.0 + 8.0e-10):
            raise AssertionError("positive physical edge work exceeded native capacity")
        expected_good = self.signed_efficiency > 1.0 - ETA0
        if self.geometry_good != expected_good:
            raise AssertionError("geometry-good flag was not read from the physical edge efficiency")


@dataclass(frozen=True, order=True)
class HardProductCell:
    """Deterministic hard event-role label for one unordered parent pair and child.

    These are observer labels only.  They may coarsen/refine the physical edge space,
    but they never create work.  Existing event-role registration supplies the Borel
    Fourier cells and pointwise orthogonal helicity projectors represented by the
    mode-role map used to construct this label.
    """

    parent_roles: tuple[str, str]
    child_role: str

    def __post_init__(self) -> None:
        if len(self.parent_roles) != 2:
            raise ValueError("hard product cell needs two unordered parent roles")
        parents = tuple(sorted(str(v) for v in self.parent_roles))
        child = str(self.child_role)
        if not all(parents) or not child:
            raise ValueError("nonempty hard role labels required")
        object.__setattr__(self, "parent_roles", parents)
        object.__setattr__(self, "child_role", child)


@dataclass(frozen=True)
class HardCellWork:
    cell: HardProductCell
    signed_work: float
    inherited_positive_work: float
    inherited_good_positive_work: float
    inherited_bad_positive_work: float
    fresh_cell_hahn_positive: float
    cancellation_gap: float
    fresh_cell_hahn_is_causal_law: bool = False

    def __post_init__(self) -> None:
        vals = (
            self.signed_work,
            self.inherited_positive_work,
            self.inherited_good_positive_work,
            self.inherited_bad_positive_work,
            self.fresh_cell_hahn_positive,
            self.cancellation_gap,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite hard-cell work required")
        if min(
            self.inherited_positive_work,
            self.inherited_good_positive_work,
            self.inherited_bad_positive_work,
            self.fresh_cell_hahn_positive,
            self.cancellation_gap,
        ) < -2.0e-12 * max(
            self.inherited_positive_work,
            self.fresh_cell_hahn_positive,
            abs(self.signed_work),
            1.0e-300,
        ):
            raise AssertionError("hard-cell positive mass/cancellation gap became negative")
        if self.fresh_cell_hahn_is_causal_law:
            raise ValueError("a later hard-cell Hahn split cannot be promoted to the master causal law")
        if not _native_close(
            self.inherited_good_positive_work + self.inherited_bad_positive_work,
            self.inherited_positive_work,
            factor=8.0e-10,
        ):
            raise AssertionError("hard cell lost the exact G/B restriction of inherited dW+")
        expected_hahn = max(self.signed_work, 0.0)
        if not _native_close(self.fresh_cell_hahn_positive, expected_hahn, factor=8.0e-10):
            raise AssertionError("hard-cell diagnostic Hahn mass is not the positive part of signed cell work")
        if self.fresh_cell_hahn_positive > self.inherited_positive_work + 8.0e-10 * max(
            self.fresh_cell_hahn_positive, self.inherited_positive_work, 1.0e-300
        ):
            raise AssertionError("positive part failed contraction under hard-cell aggregation")
        if not _native_close(
            self.cancellation_gap,
            self.inherited_positive_work - self.fresh_cell_hahn_positive,
            factor=8.0e-10,
        ):
            raise AssertionError("hard-cell cancellation gap lost pushforward provenance")


@dataclass(frozen=True)
class HardCellCompression:
    cells: tuple[HardCellWork, ...]
    canonical_positive_work: float
    inherited_positive_work: float
    inherited_good_positive_work: float
    inherited_bad_positive_work: float
    fresh_hahn_positive_work: float
    cancellation_gap: float
    deterministic_pushforward: bool = True
    fresh_hahn_is_causal_law: bool = False

    def __post_init__(self) -> None:
        vals = (
            self.canonical_positive_work,
            self.inherited_positive_work,
            self.inherited_good_positive_work,
            self.inherited_bad_positive_work,
            self.fresh_hahn_positive_work,
            self.cancellation_gap,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite hard-cell compression summary required")
        if not self.deterministic_pushforward or self.fresh_hahn_is_causal_law:
            raise ValueError("hard-cell compression must inherit dW+ and may not mint a new causal Hahn law")
        if not _native_close(self.canonical_positive_work, self.inherited_positive_work, factor=8.0e-10):
            raise AssertionError("deterministic hard-cell pushforward did not preserve canonical dW+ mass")
        if not _native_close(
            self.inherited_good_positive_work + self.inherited_bad_positive_work,
            self.inherited_positive_work,
            factor=8.0e-10,
        ):
            raise AssertionError("hard-cell compression lost exact geometry fate partition")
        if self.fresh_hahn_positive_work > self.inherited_positive_work + 8.0e-10 * max(
            self.fresh_hahn_positive_work, self.inherited_positive_work, 1.0e-300
        ):
            raise AssertionError("re-Hahn after hard aggregation exceeded inherited canonical cause")
        if not _native_close(
            self.cancellation_gap,
            self.inherited_positive_work - self.fresh_hahn_positive_work,
            factor=8.0e-10,
        ):
            raise AssertionError("compression cancellation gap changed")


@dataclass(frozen=True)
class BadPositiveWorkRoute:
    support: tuple[CanonicalPositiveEdgeOccurrence, ...]
    physical_work_mass: float
    capacity_mass: float
    signed_progress_mass: float
    deficit: float
    tau: float
    fixed_transfer_gate: FixedTransferLossGate
    transfer_partition: TransferPartition
    joint_projection: JointStopProjection

    def __post_init__(self) -> None:
        if not self.support or self.physical_work_mass <= 0.0 or self.capacity_mass <= 0.0:
            raise ValueError("nonempty bad positive-work restriction required")
        if any(edge.geometry_good for edge in self.support):
            raise AssertionError("geometry-good edge entered the bad causal sublaw")
        if not (0.0 < self.tau <= MAX_CERTIFIED_TAU):
            raise ValueError("bad route uses only the certified 0<tau<=0.1 range")
        expected_delta = self.tau * self.tau / CURVATURE_DENOM
        if not _native_close(self.fixed_transfer_gate.threshold, expected_delta, factor=5.0e-14):
            raise AssertionError("fixed-transfer threshold changed from the certified physical block interface")
        if not _native_close(self.fixed_transfer_gate.avg_transfer_deficit, self.deficit, factor=8.0e-10):
            raise AssertionError("fixed-transfer gate was not evaluated on this same bad physical restriction")
        if not self.fixed_transfer_gate.triggered or self.fixed_transfer_gate.cause != "physical_transfer_cost":
            raise AssertionError("bad physical restriction did not trigger the physical block transfer channel")
        if self.deficit < ETA0 - 2.0e-10:
            raise AssertionError("bad physical restriction failed its native eta0 transfer deficit")
        if not self.deficit > self.fixed_transfer_gate.threshold:
            raise AssertionError("bad physical restriction did not cross the certified fixed-transfer gate")
        if not _native_close(self.transfer_partition.total_mass, self.physical_work_mass, factor=8.0e-10):
            raise AssertionError("compiler was not bound to the same bad dW+ mass")
        if self.transfer_partition.xi_mass != 0.0 or not _native_close(
            self.transfer_partition.retained_mass, self.physical_work_mass, factor=8.0e-10
        ):
            raise AssertionError("bad stage-zero route altered the canonical causal sublaw before transfer loss")
        transfer_mass = self.transfer_partition.currency_mass.get(PhysicalCurrency.MULTIPLICATIVE_TRANSFER.value, 0.0)
        if not _native_close(transfer_mass, self.physical_work_mass, factor=8.0e-10):
            raise AssertionError("bad dW+ mass did not bind entirely to fixed transfer loss")
        if self.transfer_partition.first_time is not None:
            raise AssertionError("stage-zero fixed transfer loss invented a physical first-hit time")
        if self.transfer_partition.tied_causes != (PhysicalCause.TRANSFER_WORK_LOSS.value,):
            raise AssertionError("bad dW+ restriction lost TRANSFER_WORK_LOSS provenance")
        if self.joint_projection.first_time is not None:
            raise AssertionError("joint stage-zero transfer loss invented a physical first-hit time")
        if self.joint_projection.master_disposition != MasterDisposition.TRANSFER_COST.value:
            raise AssertionError("bad dW+ restriction did not terminate as TRANSFER_COST")
        if self.joint_projection.terminal_certificate_used != "stage_zero_fixed_transfer_loss":
            raise AssertionError("bad dW+ restriction lost its stage-zero terminal certificate")


@dataclass(frozen=True)
class YoungEligiblePositiveWork:
    support: tuple[CanonicalPositiveEdgeOccurrence, ...]
    physical_work_mass: float
    hard_cells: tuple[HardCellWork, ...]
    marking_good: bool = False
    registered_generated_continuation: bool = False
    young_certified: bool = False

    def __post_init__(self) -> None:
        if any(not edge.geometry_good for edge in self.support):
            raise AssertionError("geometry-bad edge entered Young eligibility")
        if self.marking_good or self.registered_generated_continuation or self.young_certified:
            raise ValueError("geometry-good work is only Young-eligible at this theorem layer")
        if not _native_close(
            self.physical_work_mass,
            _finite_sum([edge.positive_work_mass for edge in self.support], "Young-eligible positive work"),
            factor=8.0e-10,
        ):
            raise AssertionError("Young eligibility changed physical dW+ mass")


@dataclass(frozen=True)
class CanonicalPositiveEdgeWorkRouting:
    total_positive_work: float
    good_positive_work: float
    bad_positive_work: float
    mass_reconstruction_residual: float
    good_support: tuple[CanonicalPositiveEdgeOccurrence, ...]
    bad_support: tuple[CanonicalPositiveEdgeOccurrence, ...]
    bad_route: BadPositiveWorkRoute | None
    hard_cell_compression: HardCellCompression
    young_eligible: YoungEligiblePositiveWork
    capacity_used_as_causal_law: bool = False
    later_hahn_used_as_causal_law: bool = False

    def __post_init__(self) -> None:
        if self.capacity_used_as_causal_law or self.later_hahn_used_as_causal_law:
            raise ValueError("canonical routing may not replace dW+ by capacity or a later Hahn split")
        if not _native_close(self.good_positive_work + self.bad_positive_work, self.total_positive_work, factor=8.0e-10):
            raise AssertionError("canonical positive edge-work fate partition lost mass")
        if abs(self.mass_reconstruction_residual) > 8.0e-10 * max(self.total_positive_work, 1.0e-300):
            raise AssertionError("canonical positive edge-work reconstruction residual is too large")
        if (self.bad_positive_work > 0.0) != (self.bad_route is not None):
            raise AssertionError("bad route presence does not match actual bad causal mass")
        if not _native_close(self.hard_cell_compression.canonical_positive_work, self.total_positive_work, factor=8.0e-10):
            raise AssertionError("hard-cell compression is not a pushforward of this canonical dW+ law")
        if not _native_close(self.young_eligible.physical_work_mass, self.good_positive_work, factor=8.0e-10):
            raise AssertionError("Young-eligible handle is not the geometry-good dW+ restriction")


def _flatten_atoms(ledger: ContinuumEdgeMeasureLedger) -> tuple[tuple[EdgeOccurrenceKey, ContinuumModalEdgeAtom], ...]:
    out: list[tuple[EdgeOccurrenceKey, ContinuumModalEdgeAtom]] = []
    for fi, fiber in enumerate(ledger.physical_fibers):
        for ai, atom in enumerate(fiber.modal_atoms):
            out.append((EdgeOccurrenceKey(fi, ai), atom))
    if len(out) != ledger.modal_edges:
        raise AssertionError("physical edge occurrence count changed during replay")
    return tuple(out)


def _positive_occurrence(key: EdgeOccurrenceKey, atom: ContinuumModalEdgeAtom) -> CanonicalPositiveEdgeOccurrence | None:
    work = float(atom.signed_work_mass)
    if not work > 0.0:
        return None
    return CanonicalPositiveEdgeOccurrence(
        key=key,
        physical_edge_identity=atom.physical_edge_identity,
        positive_work_mass=work,
        capacity_mass=float(atom.capacity_mass),
        signed_progress_mass=float(atom.signed_progress_mass),
        signed_efficiency=float(atom.signed_efficiency),
        scale_progress=float(atom.scale_progress),
        geometry_good=bool(atom.signed_efficiency > 1.0 - ETA0),
    )


def exact_mode_role_map(ledger: ContinuumEdgeMeasureLedger) -> dict[HelicalModeIdentity, str]:
    """One deterministic hard role per exact wavevector/helicity mode identity."""
    replayed = _replay_physical_ledger(ledger)
    modes: set[HelicalModeIdentity] = set()
    for _, atom in _flatten_atoms(replayed):
        modes.update(atom.physical_edge_identity.parents)
        modes.add(atom.physical_edge_identity.child)
    return {
        mode: f"k={mode.wavevector};s={mode.helicity:+d}"
        for mode in sorted(modes)
    }


def single_hard_role_map(ledger: ContinuumEdgeMeasureLedger) -> dict[HelicalModeIdentity, str]:
    """Maximal deterministic coarsening, useful only to expose cancellation under aggregation."""
    return {mode: "all-hard-modes" for mode in exact_mode_role_map(ledger)}


def _validate_mode_roles(
    ledger: ContinuumEdgeMeasureLedger,
    mode_roles: Mapping[HelicalModeIdentity, str],
) -> dict[HelicalModeIdentity, str]:
    required: set[HelicalModeIdentity] = set()
    for _, atom in _flatten_atoms(ledger):
        required.update(atom.physical_edge_identity.parents)
        required.add(atom.physical_edge_identity.child)
    supplied = set(mode_roles)
    if supplied != required:
        missing = len(required - supplied)
        extra = len(supplied - required)
        raise ValueError(f"hard mode-role map must bind every and only physical mode identity (missing={missing}, extra={extra})")
    out = {mode: str(mode_roles[mode]) for mode in required}
    if any(not label for label in out.values()):
        raise ValueError("nonempty deterministic hard role labels required")
    return out


def compress_signed_edge_work_to_hard_cells(
    ledger: ContinuumEdgeMeasureLedger,
    mode_roles: Mapping[HelicalModeIdentity, str],
) -> HardCellCompression:
    """Compress the same signed edge law and push forward its already-fixed Hahn law.

    The deterministic hard role map is applied to physical Fourier/helicity mode
    identities.  For each product cell C we compute both

      T_C = (pi_# dW)(C)
      P_C = (pi_# dW+)(C).

    ``max(T_C,0)`` is retained only as a cancellation diagnostic.  It is never
    promoted to a second causal law.  This is precisely the signed hard-cell datum
    required before a downstream Young/Christ comparison.
    """
    replayed = _replay_physical_ledger(ledger)
    roles = _validate_mode_roles(replayed, mode_roles)
    buckets: dict[HardProductCell, list[tuple[float, float, bool]]] = {}
    for _, atom in _flatten_atoms(replayed):
        identity = atom.physical_edge_identity
        cell = HardProductCell(
            parent_roles=(roles[identity.parents[0]], roles[identity.parents[1]]),
            child_role=roles[identity.child],
        )
        work = float(atom.signed_work_mass)
        pos = max(work, 0.0)
        good = bool(pos > 0.0 and atom.signed_efficiency > 1.0 - ETA0)
        buckets.setdefault(cell, []).append((work, pos, good))

    cells: list[HardCellWork] = []
    for cell in sorted(buckets):
        rows = buckets[cell]
        signed = _finite_sum([row[0] for row in rows], "signed hard-cell work")
        inherited = _finite_sum([row[1] for row in rows], "inherited hard-cell dW+ mass")
        good = _finite_sum([row[1] for row in rows if row[2]], "hard-cell geometry-good dW+ mass")
        bad = _finite_sum([row[1] for row in rows if row[1] > 0.0 and not row[2]], "hard-cell geometry-bad dW+ mass")
        fresh = max(signed, 0.0)
        cells.append(
            HardCellWork(
                cell=cell,
                signed_work=signed,
                inherited_positive_work=inherited,
                inherited_good_positive_work=good,
                inherited_bad_positive_work=bad,
                fresh_cell_hahn_positive=fresh,
                cancellation_gap=inherited - fresh,
            )
        )

    canonical = float(replayed.positive_edge_work)
    inherited_total = _finite_sum([c.inherited_positive_work for c in cells], "hard-cell inherited dW+ total")
    good_total = _finite_sum([c.inherited_good_positive_work for c in cells], "hard-cell good dW+ total")
    bad_total = _finite_sum([c.inherited_bad_positive_work for c in cells], "hard-cell bad dW+ total")
    fresh_total = _finite_sum([c.fresh_cell_hahn_positive for c in cells], "hard-cell fresh Hahn diagnostic total")
    return HardCellCompression(
        cells=tuple(cells),
        canonical_positive_work=canonical,
        inherited_positive_work=inherited_total,
        inherited_good_positive_work=good_total,
        inherited_bad_positive_work=bad_total,
        fresh_hahn_positive_work=fresh_total,
        cancellation_gap=inherited_total - fresh_total,
    )


def route_canonical_positive_edge_work(
    ledger: ContinuumEdgeMeasureLedger,
    *,
    tau: float,
    mode_roles: Mapping[HelicalModeIdentity, str],
) -> CanonicalPositiveEdgeWorkRouting:
    """Route the canonical Hahn-positive physical edge law without changing measure.

    Capacity is consulted only after the actual dW+ restriction B has been selected,
    to certify why that same physical sublaw has fixed transfer loss.  Geometry-good
    dW+ is returned only as Young-eligible.  The hard representation carries both
    inherited dW+ and signed cell work; only the latter may feed Young/Christ.
    """
    replayed = _replay_physical_ledger(ledger)
    t = float(tau)
    if not math.isfinite(t) or not (0.0 < t <= MAX_CERTIFIED_TAU):
        raise ValueError("canonical routing is certified only for 0<tau<=0.1")

    positive: list[CanonicalPositiveEdgeOccurrence] = []
    for key, atom in _flatten_atoms(replayed):
        occurrence = _positive_occurrence(key, atom)
        if occurrence is not None:
            positive.append(occurrence)
    good = tuple(edge for edge in positive if edge.geometry_good)
    bad = tuple(edge for edge in positive if not edge.geometry_good)
    total_work = _finite_sum([e.positive_work_mass for e in positive], "canonical dW+ mass")
    good_work = _finite_sum([e.positive_work_mass for e in good], "geometry-good dW+ mass")
    bad_work = _finite_sum([e.positive_work_mass for e in bad], "geometry-bad dW+ mass")
    if not _native_close(total_work, replayed.positive_edge_work, factor=8.0e-10):
        raise AssertionError("positive edge occurrences failed to reconstruct canonical ledger dW+")
    reconstruction = (good_work + bad_work) - total_work

    bad_route: BadPositiveWorkRoute | None = None
    if bad_work > 0.0:
        bad_capacity = _finite_sum([e.capacity_mass for e in bad], "bad physical capacity")
        bad_progress = _finite_sum([e.signed_progress_mass for e in bad], "bad signed progress")
        if bad_capacity <= 0.0:
            raise AssertionError("nonzero bad dW+ restriction has zero native capacity")
        deficit = 1.0 - bad_progress / (float_jstar() * bad_capacity)
        transfer_gate = fixed_transfer_loss_gate(tau=t, avg_transfer_deficit=deficit)
        expected = t * t / CURVATURE_DENOM
        if not _native_close(transfer_gate.threshold, expected, factor=5.0e-14):
            raise AssertionError("physical fixed-transfer block interface threshold changed")
        if not transfer_gate.triggered or transfer_gate.cause != "physical_transfer_cost":
            raise AssertionError("same bad dW+ restriction failed to bind to the physical fixed-transfer block gate")
        transfer = compile_transfer_measure(
            total_mass=bad_work,
            xi_mass=0.0,
            witness=BlockWitness(fixed_transfer_loss=transfer_gate.triggered, kelvin_flat_certified=False),
        )
        projection = joint_stop_master_projection(fixed_transfer_loss=transfer_gate.triggered)
        bad_route = BadPositiveWorkRoute(
            support=bad,
            physical_work_mass=bad_work,
            capacity_mass=bad_capacity,
            signed_progress_mass=bad_progress,
            deficit=deficit,
            tau=t,
            fixed_transfer_gate=transfer_gate,
            transfer_partition=transfer,
            joint_projection=projection,
        )

    compression = compress_signed_edge_work_to_hard_cells(replayed, mode_roles)
    young_cells = tuple(cell for cell in compression.cells if cell.inherited_good_positive_work > 0.0)
    young = YoungEligiblePositiveWork(
        support=good,
        physical_work_mass=good_work,
        hard_cells=young_cells,
    )
    return CanonicalPositiveEdgeWorkRouting(
        total_positive_work=total_work,
        good_positive_work=good_work,
        bad_positive_work=bad_work,
        mass_reconstruction_residual=reconstruction,
        good_support=good,
        bad_support=bad,
        bad_route=bad_route,
        hard_cell_compression=compression,
        young_eligible=young,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "canonical_causal_law": "Hahn-positive part dW+ of the already reconstructed signed unordered Fourier/helicity NS edge measure",
        "routing_partition": "G={T_e>0,(J_e/J*)c_e>1-1e-4}; B={T_e>0,(J_e/J*)c_e<=1-1e-4}",
        "eta0": ETA0,
        "bad_deficit": "epsilon_B=1-F(B)/(J*A(B))>=1e-4 on the same measurable bad dW+ restriction",
        "fixed_transfer_threshold": f"delta_tau=tau^2/{CURVATURE_DENOM}, 0<tau<=0.1, hence delta_tau<1e-4",
        "bad_fate": "same dW+ restriction -> shared physical FixedTransferLossGate -> fixed_transfer_loss -> TRANSFER_WORK_LOSS -> TRANSFER_COST; first_time=None",
        "good_fate": "Young-eligible only; marking_good=False and generated continuation remain unproved",
        "hard_cell_handoff": "deterministic hard Fourier/helicity mode map carries pi_#dW+ as cause and pi_#dW as signed Young datum",
        "young_input": "T_C=(pi_#dW)(C), not gross (pi_#dW+)(C); physical cancellation retained before Young saturation",
        "later_hahn": "(pi_#dW)^+ is a cancellation diagnostic and never a second causal law; (pi_#dW)^+ <= pi_#dW+",
        "coherent_povm_scope": "general positive coherent localization remains an exact signed-work representation; no causal Hahn identification is made without a positive mass-preserving kernel from dW+",
        "capacity_is_causal_law": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class CanonicalPositiveEdgeRoutingStress:
    samples: int
    minimum_bad_deficit_margin: float
    minimum_fixed_transfer_margin: float
    worst_mass_reconstruction_relative: float
    worst_hard_pushforward_relative: float
    maximum_coarsened_cancellation_fraction: float
    nonforward_bad_cases: int
    geometry_good_marking_promotions: int
    first_time_sentinel_failures: int


def _pure_helical_fiber(
    *,
    x: np.ndarray,
    y: np.ndarray,
    sx: int,
    sy: int,
    sz: int,
    phase_sign: float = 1.0,
    quotient_measure_mass: float = 1.0,
) -> object:
    z = np.asarray(x, float) + np.asarray(y, float)
    g = coupling_g(np.asarray(x, float), np.asarray(y, float), -z, sx, sy, sz)
    signed_frequency = sx * float(np.linalg.norm(x)) - sy * float(np.linalg.norm(y))
    target_sign = 1.0 if signed_frequency >= 0.0 else -1.0
    az = float(phase_sign) * target_sign * np.exp(-1j * np.angle(g))
    return register_continuum_triad_fiber(
        x=np.asarray(x, float),
        y=np.asarray(y, float),
        z=z,
        ux=helical_basis(np.asarray(x, float), sx),
        uy=helical_basis(np.asarray(y, float), sy),
        uz=az * helical_basis(z, sz),
        quotient_measure_mass=quotient_measure_mass,
    )


def _near_extremal_positive_fiber(quotient_measure_mass: float) -> object:
    rstar = symmetric_rstar()
    gamma = symmetric_gamma(rstar)
    nx = math.exp(-gamma)
    ny = math.exp(-gamma)
    xx = 0.5 * (1.0 + nx * nx - ny * ny)
    yy = math.sqrt(nx * nx - xx * xx)
    return _pure_helical_fiber(
        x=np.array([xx, yy, 0.0]),
        y=np.array([1.0 - xx, -yy, 0.0]),
        sx=1,
        sy=-1,
        sz=1,
        quotient_measure_mass=quotient_measure_mass,
    )


def _nonforward_positive_fiber(quotient_measure_mass: float, *, phase_sign: float = 1.0) -> object:
    return _pure_helical_fiber(
        x=np.array([1.0, 0.0, 0.0]),
        y=np.array([-0.8, 0.6, 0.0]),
        sx=1,
        sy=-1,
        sz=1,
        phase_sign=phase_sign,
        quotient_measure_mass=quotient_measure_mass,
    )


def stress(samples: int = 50_000, seed: int = 20260812) -> CanonicalPositiveEdgeRoutingStress:
    if samples <= 0:
        raise ValueError("positive stress sample count required")
    rng = np.random.default_rng(seed)
    min_bad = math.inf
    min_transfer = math.inf
    worst_mass = 0.0
    worst_push = 0.0
    max_cancel = 0.0
    nonforward = 0
    marking = 0
    first_time_fail = 0

    for _ in range(samples):
        qg = float(np.exp(rng.uniform(-2.0, 2.0)))
        qb = float(np.exp(rng.uniform(-2.0, 2.0)))
        qn = float(np.exp(rng.uniform(-3.0, 1.0)))
        good = _near_extremal_positive_fiber(qg)
        bad = _nonforward_positive_fiber(qb)
        negative = _nonforward_positive_fiber(qn, phase_sign=-1.0)
        ledger = continuum_edge_measure_ledger((good, bad, negative))
        tau = float(np.exp(rng.uniform(math.log(1.0e-5), math.log(MAX_CERTIFIED_TAU))))
        exact_roles = exact_mode_role_map(ledger)
        out = route_canonical_positive_edge_work(ledger, tau=tau, mode_roles=exact_roles)
        if out.bad_route is None:
            raise AssertionError("physical nonforward positive edge disappeared from bad route")
        nonforward += 1
        min_bad = min(min_bad, out.bad_route.deficit - ETA0)
        min_transfer = min(min_transfer, out.bad_route.deficit - out.bad_route.fixed_transfer_gate.threshold)
        scale = max(out.total_positive_work, 1.0e-300)
        worst_mass = max(worst_mass, abs(out.mass_reconstruction_residual) / scale)
        worst_push = max(
            worst_push,
            abs(out.hard_cell_compression.inherited_positive_work - out.total_positive_work) / scale,
        )
        if out.young_eligible.marking_good or out.young_eligible.young_certified:
            marking += 1
        if out.bad_route.transfer_partition.first_time is not None:
            first_time_fail += 1

        coarse = compress_signed_edge_work_to_hard_cells(ledger, single_hard_role_map(ledger))
        if coarse.inherited_positive_work > 0.0:
            max_cancel = max(max_cancel, coarse.cancellation_gap / coarse.inherited_positive_work)
        if not _native_close(coarse.inherited_positive_work, out.total_positive_work, factor=8.0e-10):
            raise AssertionError("analyst coarsening changed inherited canonical dW+ mass")

    return CanonicalPositiveEdgeRoutingStress(
        samples=samples,
        minimum_bad_deficit_margin=min_bad,
        minimum_fixed_transfer_margin=min_transfer,
        worst_mass_reconstruction_relative=worst_mass,
        worst_hard_pushforward_relative=worst_push,
        maximum_coarsened_cancellation_fraction=max_cancel,
        nonforward_bad_cases=nonforward,
        geometry_good_marking_promotions=marking,
        first_time_sentinel_failures=first_time_fail,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--outdir", type=Path, default=Path("results-canonical-positive-edge-work-routing"))
    args = parser.parse_args()
    out = stress(samples=args.samples, seed=args.seed)
    cert = theorem_certificate()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (args.outdir / "stress.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    summary = f"""# Canonical positive edge-work routing\n\nStatus: **{STATUS}**.\n\n- samples: `{out.samples}`\n- minimum bad-deficit margin above eta0: `{out.minimum_bad_deficit_margin:.12g}`\n- minimum bad-deficit margin above fixed-transfer threshold: `{out.minimum_fixed_transfer_margin:.12g}`\n- worst exact G+B mass reconstruction relative residual: `{out.worst_mass_reconstruction_relative:.12g}`\n- worst deterministic hard-cell dW+ pushforward relative residual: `{out.worst_hard_pushforward_relative:.12g}`\n- maximum coarsened signed-cell cancellation fraction: `{out.maximum_coarsened_cancellation_fraction:.12g}`\n- retained positive nonforward bad cases: `{out.nonforward_bad_cases}`\n- geometry-good -> marking-good promotions: `{out.geometry_good_marking_promotions}`\n- fixed-transfer first-time sentinel failures: `{out.first_time_sentinel_failures}`\n\nThe causal law is the already-fixed Hahn-positive part of signed physical edge work. Capacity is not causal. Later hard cells inherit that positive law and separately retain signed cell work for Young/Christ. No global Navier--Stokes regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
