from __future__ import annotations

import argparse
import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.canonical_positive_edge_work_routing import ETA0, HardProductCell
from src.continuum_helical_edge_measure_registration import HelicalModeIdentity
from src.cyclic_helical_triad_donor_kernel import (
    ClosedHelicalTriadRegistration,
    CyclicTriadMeasureKernel,
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    register_closed_helical_triad,
    signed_good_integer_triad,
)
from src.helical import stable_norm3

STATUS = (
    "EXACT_CYCLIC_HARD_CELL_SINGLE_CHARGE_QUOTIENT__RESTRICTED_DW_MINUS_ROWS__"
    "CANONICAL_DW_PLUS_GOOD_BAD_COLUMNS__OVERLAPPING_RECIPIENT_SUPPORT_NO_DOUBLE_CHARGE__"
    "COARSE_SELF_LOOPS_ZERO_RECURSION_DEPTH"
)

FATE_GOOD = "geometry_good_young_eligible"
FATE_BAD = "geometry_bad_transfer_work_loss"


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


def _native_residual(actual: float, expected: float, native_scale: float) -> float:
    a = float(actual)
    b = float(expected)
    scale = float(native_scale)
    if not all(math.isfinite(v) for v in (a, b, scale)) or scale < 0.0:
        raise ValueError("finite nonnegative native work scale required")
    gap = abs(a - b)
    if scale == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / scale


def required_hard_modes(triad: ClosedHelicalTriadRegistration) -> tuple[HelicalModeIdentity, ...]:
    modes: set[HelicalModeIdentity] = set()
    for slot in triad.slots:
        modes.update(slot.edge_identity.parents)
        modes.add(slot.edge_identity.child)
    return tuple(sorted(modes))


def fine_hard_role_map(triad: ClosedHelicalTriadRegistration) -> dict[HelicalModeIdentity, str]:
    """Injective deterministic hard-role map on the six physical root/parent modes."""
    return {mode: f"hard-mode-{j}" for j, mode in enumerate(required_hard_modes(triad))}


def single_hard_role_map(triad: ClosedHelicalTriadRegistration) -> dict[HelicalModeIdentity, str]:
    """Maximal coarsening used to expose representation-level self-loops."""
    return {mode: "all-hard-modes" for mode in required_hard_modes(triad)}


def _validated_roles(
    triad: ClosedHelicalTriadRegistration,
    mode_roles: Mapping[HelicalModeIdentity, str],
) -> dict[HelicalModeIdentity, str]:
    required = set(required_hard_modes(triad))
    supplied = set(mode_roles)
    if supplied != required:
        raise ValueError(
            "hard role map must bind every and only physical mode used by the three cyclic roots "
            f"(missing={len(required-supplied)}, extra={len(supplied-required)})"
        )
    out = {mode: str(mode_roles[mode]) for mode in required}
    if any(not label for label in out.values()):
        raise ValueError("nonempty deterministic hard role labels required")
    return out


def _hard_cell_for_slot(slot, roles: Mapping[HelicalModeIdentity, str]) -> HardProductCell:
    identity = slot.edge_identity
    return HardProductCell(
        parent_roles=(roles[identity.parents[0]], roles[identity.parents[1]]),
        child_role=roles[identity.child],
    )


def _recipient_fate(slot) -> str:
    if slot.positive_work <= 0.0:
        raise ValueError("recipient fate is defined only on canonical positive child work")
    return FATE_GOOD if slot.signed_efficiency > 1.0 - ETA0 else FATE_BAD


@dataclass(frozen=True)
class HardCellDonorRecipientAtom:
    donor_closed_mode_index: int
    recipient_closed_mode_index: int
    donor_cell: HardProductCell
    recipient_cell: HardProductCell
    recipient_fate: str
    physical_work_mass: float
    coarse_self_loop: bool
    creates_new_event: bool = False
    creates_scale_progress: bool = False
    replaces_recipient_causal_law: bool = False

    def __post_init__(self) -> None:
        if self.donor_closed_mode_index == self.recipient_closed_mode_index:
            raise ValueError("physical cyclic donor and recipient roots must remain distinct")
        if self.recipient_fate not in (FATE_GOOD, FATE_BAD):
            raise ValueError("recipient fate must be inherited from canonical positive-edge routing")
        _finite_nonnegative(self.physical_work_mass, "hard-cell donor/recipient work mass")
        if self.physical_work_mass <= 0.0:
            raise ValueError("hard-cell donor/recipient atom must carry positive physical work")
        if self.coarse_self_loop != (self.donor_cell == self.recipient_cell):
            raise AssertionError("coarse self-loop flag changed under hard-cell pushforward")
        if self.creates_new_event or self.creates_scale_progress or self.replaces_recipient_causal_law:
            raise ValueError("hard-cell donor provenance may not create cause, event depth, or scale progress")


@dataclass(frozen=True)
class HardCellDonorCharge:
    cell: HardProductCell
    canonical_negative_work_mass: float
    outgoing_donor_mass: float
    marginal_native_residual: float

    def __post_init__(self) -> None:
        _finite_nonnegative(self.canonical_negative_work_mass, "canonical hard-cell negative work")
        _finite_nonnegative(self.outgoing_donor_mass, "outgoing hard-cell donor mass")
        _finite_nonnegative(self.marginal_native_residual, "donor marginal native residual")
        if self.canonical_negative_work_mass <= 0.0:
            raise ValueError("donor charge must carry positive canonical dW- mass")
        if self.marginal_native_residual > 5.0e-10:
            raise AssertionError("hard-cell donor row left canonical pi_#dW- on its native work scale")


@dataclass(frozen=True)
class HardCellRecipientCharge:
    cell: HardProductCell
    fate: str
    canonical_positive_work_mass: float
    incoming_donor_mass: float
    marginal_native_residual: float
    incoming_donor_cell_count: int
    charged_once_downstream: bool = True
    creates_new_causal_law: bool = False

    def __post_init__(self) -> None:
        if self.fate not in (FATE_GOOD, FATE_BAD):
            raise ValueError("recipient charge fate must be geometry-good or geometry-bad")
        _finite_nonnegative(self.canonical_positive_work_mass, "canonical recipient dW+ mass")
        _finite_nonnegative(self.incoming_donor_mass, "incoming donor provenance mass")
        _finite_nonnegative(self.marginal_native_residual, "recipient marginal native residual")
        if self.canonical_positive_work_mass <= 0.0 or self.incoming_donor_cell_count <= 0:
            raise ValueError("recipient charge must carry positive canonical dW+ with donor provenance")
        if self.marginal_native_residual > 5.0e-10:
            raise AssertionError("hard-cell donor columns did not reconstruct canonical pi_#dW+ fate charge")
        if not self.charged_once_downstream or self.creates_new_causal_law:
            raise ValueError("donor disintegration may not multiply or replace the canonical recipient cause")

    @property
    def downstream_semantics(self) -> str:
        if self.fate == FATE_BAD:
            return "existing_TRANSFER_WORK_LOSS_first_time_None"
        return "existing_geometry_good_Young_eligible_route"


@dataclass(frozen=True)
class RestrictedRecipientSubcharge:
    cell: HardProductCell
    fate: str
    submeasure_mass: float
    full_canonical_positive_work_mass: float
    native_work_mass_scale: float
    incoming_selected_donor_cell_count: int

    def __post_init__(self) -> None:
        if self.fate not in (FATE_GOOD, FATE_BAD):
            raise ValueError("restricted recipient fate must be inherited from canonical routing")
        _finite_nonnegative(self.submeasure_mass, "restricted recipient submeasure mass")
        _finite_nonnegative(self.full_canonical_positive_work_mass, "full canonical recipient dW+ mass")
        _finite_nonnegative(self.native_work_mass_scale, "restricted recipient native work-mass scale")
        if self.submeasure_mass <= 0.0 or self.native_work_mass_scale <= 0.0 or self.incoming_selected_donor_cell_count <= 0:
            raise ValueError("restricted recipient subcharge must carry positive donor provenance and native work scale")
        # The exact theorem is submeasure domination.  Floating certification must
        # not divide by an accidentally tiny realized recipient Hahn mass: near
        # phase cancellation that sign/mass is ill-conditioned.  Use only the
        # predecessor theorem's immutable native work-error envelope.
        if self.submeasure_mass - self.full_canonical_positive_work_mass > 5.0e-10 * self.native_work_mass_scale:
            raise AssertionError("restricted donor pushforward exceeded its canonical recipient charge on the native work scale")


@dataclass(frozen=True)
class RestrictedHardCellDonorPushforward:
    selected_donor_cells: tuple[HardProductCell, ...]
    selected_negative_work_mass: float
    recipient_subcharges: tuple[RestrictedRecipientSubcharge, ...]
    recipient_total_mass: float
    native_work_mass_scale: float
    mass_conservation_native_residual: float
    dominated_by_full_recipient_charge: bool
    creates_new_event: bool = False
    creates_new_owner: bool = False
    later_hahn_used: bool = False
    between_time_matching_used: bool = False

    def __post_init__(self) -> None:
        if not self.selected_donor_cells or len(set(self.selected_donor_cells)) != len(self.selected_donor_cells):
            raise ValueError("nonempty unique hard-cell donor restriction required")
        _finite_nonnegative(self.selected_negative_work_mass, "selected hard-cell negative work")
        _finite_nonnegative(self.recipient_total_mass, "restricted recipient mass")
        if self.selected_negative_work_mass <= 0.0 or self.native_work_mass_scale <= 0.0:
            raise ValueError("resolved donor restriction requires positive physical and native work mass")
        if self.mass_conservation_native_residual > 5.0e-10:
            raise AssertionError("restricted hard-cell dW- did not push forward mass-preservingly")
        if not self.dominated_by_full_recipient_charge:
            raise AssertionError("restricted donor pushforward exceeded canonical recipient dW+ charge")
        if self.creates_new_event or self.creates_new_owner or self.later_hahn_used or self.between_time_matching_used:
            raise ValueError("restricted same-time donor provenance changed causal or temporal semantics")


@dataclass(frozen=True)
class CyclicHardCellSingleChargeQuotient:
    atoms: tuple[HardCellDonorRecipientAtom, ...]
    donor_charges: tuple[HardCellDonorCharge, ...]
    recipient_charges: tuple[HardCellRecipientCharge, ...]
    total_negative_work_mass: float
    total_positive_work_mass: float
    good_recipient_mass: float
    bad_recipient_mass: float
    self_loop_mass: float
    self_loop_atom_count: int
    overlapping_recipient_charge_count: int
    native_work_mass_scale: float
    total_balance_native_residual: float
    worst_donor_marginal_native_residual: float
    worst_recipient_marginal_native_residual: float
    recipient_fate_partition_native_residual: float
    donor_partition_recombines_to_canonical_recipient: bool = True
    canonical_recipient_charged_once: bool = True
    self_loops_zero_recursion_depth: bool = True
    self_loops_create_scale_progress: bool = False
    donor_work_is_new_owner: bool = False
    negative_failure_payment_inferred: bool = False
    between_time_inventory_matching_inferred: bool = False

    def __post_init__(self) -> None:
        vals = (
            self.total_negative_work_mass,
            self.total_positive_work_mass,
            self.good_recipient_mass,
            self.bad_recipient_mass,
            self.self_loop_mass,
            self.native_work_mass_scale,
            self.total_balance_native_residual,
            self.worst_donor_marginal_native_residual,
            self.worst_recipient_marginal_native_residual,
            self.recipient_fate_partition_native_residual,
        )
        if not all(math.isfinite(float(v)) and float(v) >= 0.0 for v in vals):
            raise ValueError("finite nonnegative hard-cell single-charge data required")
        if self.native_work_mass_scale <= 0.0:
            raise ValueError("hard-cell quotient requires resolved positive native work mass")
        if max(
            self.total_balance_native_residual,
            self.worst_donor_marginal_native_residual,
            self.worst_recipient_marginal_native_residual,
            self.recipient_fate_partition_native_residual,
        ) > 5.0e-10:
            raise AssertionError("hard-cell single-charge quotient left the native physical work scale")
        if self.self_loop_atom_count != sum(atom.coarse_self_loop for atom in self.atoms):
            raise AssertionError("coarse self-loop count changed")
        expected_self = math.fsum(atom.physical_work_mass for atom in self.atoms if atom.coarse_self_loop)
        if _native_residual(expected_self, self.self_loop_mass, self.native_work_mass_scale) > 5.0e-10:
            raise AssertionError("coarse self-loop mass changed")
        keys = tuple((charge.cell, charge.fate) for charge in self.recipient_charges)
        if len(set(keys)) != len(keys):
            raise AssertionError("canonical recipient charge was duplicated after donor disintegration")
        if not (
            self.donor_partition_recombines_to_canonical_recipient
            and self.canonical_recipient_charged_once
            and self.self_loops_zero_recursion_depth
        ):
            raise ValueError("single-charge or same-event self-loop semantics were weakened")
        if (
            self.self_loops_create_scale_progress
            or self.donor_work_is_new_owner
            or self.negative_failure_payment_inferred
            or self.between_time_inventory_matching_inferred
        ):
            raise ValueError("hard-cell donor quotient invented progress, ownership, failure payment, or time matching")


def hard_cell_single_charge_quotient(
    triad: ClosedHelicalTriadRegistration,
    *,
    quotient_measure_mass: float,
    mode_roles: Mapping[HelicalModeIdentity, str],
) -> CyclicHardCellSingleChargeQuotient:
    """Push the certified cyclic donor kernel through a deterministic hard-cell map.

    Mathematically, for donor hard cell C and recipient hard cell/fate (D,F),

        K(C,D,F) = M_triangle({donor in C},{recipient in D and fate F}).

    Its row marginal is the inherited hard-cell Hahn-negative law.  Its column/fate
    marginal is the already-canonical Hahn-positive recipient law.  This is a
    positive disintegration of existing cause, not another causal measure.
    """
    roles = _validated_roles(triad, mode_roles)
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=quotient_measure_mass)
    if not kernel.numerically_resolved_transport or kernel.native_work_mass_scale <= 0.0:
        raise ValueError("hard-cell quotient refuses numerically unresolved near-zero cyclic work")

    cell_by_index = {
        slot.closed_mode_index: _hard_cell_for_slot(slot, roles)
        for slot in triad.slots
    }
    fate_by_index = {
        slot.closed_mode_index: _recipient_fate(slot)
        for slot in triad.slots
        if slot.positive_work > 0.0
    }

    atoms = tuple(
        HardCellDonorRecipientAtom(
            donor_closed_mode_index=atom.donor_closed_mode_index,
            recipient_closed_mode_index=atom.recipient_closed_mode_index,
            donor_cell=cell_by_index[atom.donor_closed_mode_index],
            recipient_cell=cell_by_index[atom.recipient_closed_mode_index],
            recipient_fate=fate_by_index[atom.recipient_closed_mode_index],
            physical_work_mass=atom.physical_work_mass,
            coarse_self_loop=(
                cell_by_index[atom.donor_closed_mode_index]
                == cell_by_index[atom.recipient_closed_mode_index]
            ),
        )
        for atom in kernel.atoms
    )
    if not atoms:
        raise AssertionError("resolved cyclic kernel lost all hard-cell donor atoms")

    native = kernel.native_work_mass_scale
    donor_cells = tuple(sorted({cell_by_index[i] for i, mass in enumerate(kernel.donor_edge_negative_masses) if mass > 0.0}))
    donor_charges: list[HardCellDonorCharge] = []
    for cell in donor_cells:
        expected = math.fsum(
            kernel.donor_edge_negative_masses[i]
            for i in range(3)
            if kernel.donor_edge_negative_masses[i] > 0.0 and cell_by_index[i] == cell
        )
        actual = math.fsum(atom.physical_work_mass for atom in atoms if atom.donor_cell == cell)
        donor_charges.append(
            HardCellDonorCharge(
                cell=cell,
                canonical_negative_work_mass=expected,
                outgoing_donor_mass=actual,
                marginal_native_residual=_native_residual(actual, expected, native),
            )
        )

    recipient_keys = tuple(sorted({
        (cell_by_index[i], fate_by_index[i])
        for i, mass in enumerate(kernel.recipient_edge_positive_masses)
        if mass > 0.0
    }, key=lambda item: (item[0], item[1])))
    recipient_charges: list[HardCellRecipientCharge] = []
    overlap_count = 0
    for cell, fate in recipient_keys:
        expected = math.fsum(
            kernel.recipient_edge_positive_masses[i]
            for i in range(3)
            if kernel.recipient_edge_positive_masses[i] > 0.0
            and cell_by_index[i] == cell
            and fate_by_index[i] == fate
        )
        selected_atoms = tuple(atom for atom in atoms if atom.recipient_cell == cell and atom.recipient_fate == fate)
        actual = math.fsum(atom.physical_work_mass for atom in selected_atoms)
        donor_cell_count = len({atom.donor_cell for atom in selected_atoms})
        if donor_cell_count > 1:
            overlap_count += 1
        recipient_charges.append(
            HardCellRecipientCharge(
                cell=cell,
                fate=fate,
                canonical_positive_work_mass=expected,
                incoming_donor_mass=actual,
                marginal_native_residual=_native_residual(actual, expected, native),
                incoming_donor_cell_count=donor_cell_count,
            )
        )

    total_n = math.fsum(charge.canonical_negative_work_mass for charge in donor_charges)
    total_p = math.fsum(charge.canonical_positive_work_mass for charge in recipient_charges)
    good = math.fsum(charge.canonical_positive_work_mass for charge in recipient_charges if charge.fate == FATE_GOOD)
    bad = math.fsum(charge.canonical_positive_work_mass for charge in recipient_charges if charge.fate == FATE_BAD)
    self_mass = math.fsum(atom.physical_work_mass for atom in atoms if atom.coarse_self_loop)
    return CyclicHardCellSingleChargeQuotient(
        atoms=atoms,
        donor_charges=tuple(donor_charges),
        recipient_charges=tuple(recipient_charges),
        total_negative_work_mass=total_n,
        total_positive_work_mass=total_p,
        good_recipient_mass=good,
        bad_recipient_mass=bad,
        self_loop_mass=self_mass,
        self_loop_atom_count=sum(atom.coarse_self_loop for atom in atoms),
        overlapping_recipient_charge_count=overlap_count,
        native_work_mass_scale=native,
        total_balance_native_residual=_native_residual(total_n, total_p, native),
        worst_donor_marginal_native_residual=max(charge.marginal_native_residual for charge in donor_charges),
        worst_recipient_marginal_native_residual=max(charge.marginal_native_residual for charge in recipient_charges),
        recipient_fate_partition_native_residual=_native_residual(good + bad, total_p, native),
    )


def aggregate_hard_cell_single_charge_quotients(
    quotients: Sequence[CyclicHardCellSingleChargeQuotient],
) -> CyclicHardCellSingleChargeQuotient:
    """Linearly assemble disjoint closed-triad measure pieces at hard resolution.

    Different closed-triad or quadrature pieces may land on the same hard recipient
    cell.  Their positive recipient masses are summed into one ``(cell,fate)``
    charge.  The input pieces remain physical measure contributions; aggregation
    creates no additional multiplicity or causal law.
    """
    pieces = tuple(quotients)
    if not pieces:
        raise ValueError("at least one hard-cell cyclic quotient is required")
    native = math.fsum(piece.native_work_mass_scale for piece in pieces)
    if not math.isfinite(native) or native <= 0.0:
        raise ValueError("aggregate hard-cell quotient requires positive native work mass")
    atoms = tuple(atom for piece in pieces for atom in piece.atoms)
    if not atoms:
        raise AssertionError("aggregate hard-cell quotient lost physical donor atoms")

    donor_expected: dict[HardProductCell, float] = {}
    for piece in pieces:
        for charge in piece.donor_charges:
            donor_expected[charge.cell] = donor_expected.get(charge.cell, 0.0) + charge.canonical_negative_work_mass
    donor_charges: list[HardCellDonorCharge] = []
    for cell in sorted(donor_expected):
        expected = donor_expected[cell]
        actual = math.fsum(atom.physical_work_mass for atom in atoms if atom.donor_cell == cell)
        donor_charges.append(
            HardCellDonorCharge(
                cell=cell,
                canonical_negative_work_mass=expected,
                outgoing_donor_mass=actual,
                marginal_native_residual=_native_residual(actual, expected, native),
            )
        )

    recipient_expected: dict[tuple[HardProductCell, str], float] = {}
    for piece in pieces:
        for charge in piece.recipient_charges:
            key = (charge.cell, charge.fate)
            recipient_expected[key] = recipient_expected.get(key, 0.0) + charge.canonical_positive_work_mass
    recipient_charges: list[HardCellRecipientCharge] = []
    overlap_count = 0
    for cell, fate in sorted(recipient_expected, key=lambda item: (item[0], item[1])):
        expected = recipient_expected[(cell, fate)]
        selected_atoms = tuple(
            atom for atom in atoms if atom.recipient_cell == cell and atom.recipient_fate == fate
        )
        actual = math.fsum(atom.physical_work_mass for atom in selected_atoms)
        donor_cell_count = len({atom.donor_cell for atom in selected_atoms})
        if donor_cell_count > 1:
            overlap_count += 1
        recipient_charges.append(
            HardCellRecipientCharge(
                cell=cell,
                fate=fate,
                canonical_positive_work_mass=expected,
                incoming_donor_mass=actual,
                marginal_native_residual=_native_residual(actual, expected, native),
                incoming_donor_cell_count=donor_cell_count,
            )
        )

    total_n = math.fsum(charge.canonical_negative_work_mass for charge in donor_charges)
    total_p = math.fsum(charge.canonical_positive_work_mass for charge in recipient_charges)
    good = math.fsum(charge.canonical_positive_work_mass for charge in recipient_charges if charge.fate == FATE_GOOD)
    bad = math.fsum(charge.canonical_positive_work_mass for charge in recipient_charges if charge.fate == FATE_BAD)
    self_mass = math.fsum(atom.physical_work_mass for atom in atoms if atom.coarse_self_loop)
    return CyclicHardCellSingleChargeQuotient(
        atoms=atoms,
        donor_charges=tuple(donor_charges),
        recipient_charges=tuple(recipient_charges),
        total_negative_work_mass=total_n,
        total_positive_work_mass=total_p,
        good_recipient_mass=good,
        bad_recipient_mass=bad,
        self_loop_mass=self_mass,
        self_loop_atom_count=sum(atom.coarse_self_loop for atom in atoms),
        overlapping_recipient_charge_count=overlap_count,
        native_work_mass_scale=native,
        total_balance_native_residual=_native_residual(total_n, total_p, native),
        worst_donor_marginal_native_residual=max(charge.marginal_native_residual for charge in donor_charges),
        worst_recipient_marginal_native_residual=max(charge.marginal_native_residual for charge in recipient_charges),
        recipient_fate_partition_native_residual=_native_residual(good + bad, total_p, native),
    )


def pushforward_restricted_hard_cell_donor_work(
    quotient: CyclicHardCellSingleChargeQuotient,
    *,
    donor_cells: Sequence[HardProductCell],
) -> RestrictedHardCellDonorPushforward:
    selected = tuple(sorted(set(donor_cells)))
    if not selected:
        raise ValueError("nonempty donor hard-cell restriction required")
    available = {charge.cell: charge for charge in quotient.donor_charges}
    if any(cell not in available for cell in selected):
        raise ValueError("selected hard cell carries no canonical negative donor charge")
    selected_mass = math.fsum(available[cell].canonical_negative_work_mass for cell in selected)
    native = quotient.native_work_mass_scale
    subcharges: list[RestrictedRecipientSubcharge] = []
    dominated = True
    for full in quotient.recipient_charges:
        selected_atoms = tuple(
            atom for atom in quotient.atoms
            if atom.donor_cell in selected
            and atom.recipient_cell == full.cell
            and atom.recipient_fate == full.fate
        )
        mass = math.fsum(atom.physical_work_mass for atom in selected_atoms)
        if mass <= 0.0:
            continue
        dominated = dominated and mass <= full.canonical_positive_work_mass + 5.0e-10 * native
        subcharges.append(
            RestrictedRecipientSubcharge(
                cell=full.cell,
                fate=full.fate,
                submeasure_mass=mass,
                full_canonical_positive_work_mass=full.canonical_positive_work_mass,
                native_work_mass_scale=native,
                incoming_selected_donor_cell_count=len({atom.donor_cell for atom in selected_atoms}),
            )
        )
    total = math.fsum(charge.submeasure_mass for charge in subcharges)
    return RestrictedHardCellDonorPushforward(
        selected_donor_cells=selected,
        selected_negative_work_mass=selected_mass,
        recipient_subcharges=tuple(subcharges),
        recipient_total_mass=total,
        native_work_mass_scale=native,
        mass_conservation_native_residual=_native_residual(total, selected_mass, native),
        dominated_by_full_recipient_charge=dominated,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "hard_cell_kernel": "K(C,D,F)=M_triangle(donor root in C, recipient root in D with inherited canonical fate F)",
        "row_marginal": "sum_{D,F} K(C,D,F) = (pi_# dW-)(C) on the closed-triad donor restriction",
        "column_fate_marginal": "sum_C K(C,D,G)= (pi_# dW_G+)(D) and sum_C K(C,D,B)= (pi_# dW_B+)(D)",
        "single_charge": "incoming donor atoms disintegrate one already-canonical recipient charge; overlapping donor support and multiple closed-triad measure pieces aggregate to one (cell,fate) charge rather than creating another dW+ cause",
        "bad_recipient": "geometry-bad recipient remains on the existing TRANSFER_WORK_LOSS stage-zero route with first_time=None",
        "good_recipient": "geometry-good recipient remains on the existing Young-eligible route; donor provenance does not certify Young/Christ",
        "coarse_self_loop": "pi(donor)=pi(recipient) is allowed; it is same-time physical redistribution with zero event depth and no scale progress",
        "failure_anti_shortcut": "the theorem does not infer that n_C pays failed good work g_C; domination inequalities are not causal maps",
        "time_anti_shortcut": "same-time cyclic provenance is not a between-time matching of earlier deposits to later withdrawals; no FIFO/LIFO or modal inventory matching is introduced",
        "capacity_is_causal_law": False,
        "later_hahn_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class CyclicHardCellSingleChargeStress:
    samples: int
    resolved_cases: int
    numerically_unresolved_cases: int
    one_donor_cases: int
    two_donor_cases: int
    cases_with_overlapping_recipient_charge: int
    cases_with_coarse_self_loop: int
    worst_balance_native_residual: float
    worst_donor_marginal_native_residual: float
    worst_recipient_marginal_native_residual: float
    worst_fate_partition_native_residual: float
    worst_restricted_pushforward_native_residual: float
    signed_good_good_recipient_mass: float
    signed_good_bad_recipient_mass: float
    signed_good_coarse_self_loop_fraction: float
    generic_two_donor_single_recipient_charge_verified: bool


def _random_closed_triad(rng: np.random.Generator) -> ClosedHelicalTriadRegistration:
    while True:
        k0 = rng.normal(size=3)
        k1 = rng.normal(size=3)
        k2 = -(k0 + k1)
        norms = tuple(stable_norm3(k) for k in (k0, k1, k2))
        if min(norms) > 0.08 and min(
            stable_norm3(k0-k1), stable_norm3(k1-k2), stable_norm3(k2-k0)
        ) > 0.05:
            break
    helicities = tuple(int(x) for x in rng.choice((-1, 1), size=3))
    amps = tuple(complex(x) for x in (rng.normal(size=3) + 1j * rng.normal(size=3)))
    return register_closed_helical_triad(
        wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=amps
    )


def _random_role_map(
    triad: ClosedHelicalTriadRegistration, rng: np.random.Generator
) -> dict[HelicalModeIdentity, str]:
    modes = required_hard_modes(triad)
    count = int(rng.integers(1, 5))
    return {mode: f"coarse-{int(rng.integers(0, count))}" for mode in modes}


def stress(samples: int = 75_000, seed: int = 2026081204) -> CyclicHardCellSingleChargeStress:
    if samples <= 0:
        raise ValueError("positive hard-cell single-charge stress sample count required")
    rng = np.random.default_rng(seed)
    resolved = unresolved = one = two = overlap = self_loop = 0
    wb = wd = wr = wf = wrestrict = 0.0
    for _ in range(int(samples)):
        triad = _random_closed_triad(rng)
        qmass = math.exp(float(rng.uniform(-8.0, 8.0)))
        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=qmass)
        if not kernel.numerically_resolved_transport:
            unresolved += 1
            continue
        quotient = hard_cell_single_charge_quotient(
            triad,
            quotient_measure_mass=qmass,
            mode_roles=_random_role_map(triad, rng),
        )
        resolved += 1
        if triad.donor_kernel.donor_count == 1:
            one += 1
        elif triad.donor_kernel.donor_count == 2:
            two += 1
        else:
            raise AssertionError("resolved cyclic donor count left {1,2}")
        overlap += int(quotient.overlapping_recipient_charge_count > 0)
        self_loop += int(quotient.self_loop_atom_count > 0)
        wb = max(wb, quotient.total_balance_native_residual)
        wd = max(wd, quotient.worst_donor_marginal_native_residual)
        wr = max(wr, quotient.worst_recipient_marginal_native_residual)
        wf = max(wf, quotient.recipient_fate_partition_native_residual)
        donor_cells = tuple(charge.cell for charge in quotient.donor_charges)
        chosen = donor_cells[:1] if len(donor_cells) == 1 or rng.random() < 0.5 else donor_cells
        restricted = pushforward_restricted_hard_cell_donor_work(quotient, donor_cells=chosen)
        wrestrict = max(wrestrict, restricted.mass_conservation_native_residual)

    signed_good, _ = signed_good_integer_triad()
    signed_q = hard_cell_single_charge_quotient(
        signed_good,
        quotient_measure_mass=1.0,
        mode_roles=fine_hard_role_map(signed_good),
    )
    coarse_q = hard_cell_single_charge_quotient(
        signed_good,
        quotient_measure_mass=1.0,
        mode_roles=single_hard_role_map(signed_good),
    )
    if signed_q.good_recipient_mass <= 0.0 or signed_q.bad_recipient_mass <= 0.0:
        raise AssertionError("signed-good triad did not split donor work into existing good and bad recipient fates")
    if coarse_q.self_loop_mass <= 0.0 or coarse_q.self_loop_atom_count != len(coarse_q.atoms):
        raise AssertionError("maximal hard coarsening did not expose the physical coarse self-loop anti-theorem")

    two_donor = generic_two_donor_counterexample()
    two_q = hard_cell_single_charge_quotient(
        two_donor,
        quotient_measure_mass=1.0,
        mode_roles=fine_hard_role_map(two_donor),
    )
    two_ok = (
        two_donor.donor_kernel.donor_count == 2
        and two_donor.donor_kernel.recipient_count == 1
        and len(two_q.recipient_charges) == 1
        and two_q.recipient_charges[0].incoming_donor_cell_count == 2
        and two_q.overlapping_recipient_charge_count == 1
    )
    if not two_ok:
        raise AssertionError("generic two-donor triad did not recombine to one canonical recipient charge")

    return CyclicHardCellSingleChargeStress(
        samples=int(samples),
        resolved_cases=resolved,
        numerically_unresolved_cases=unresolved,
        one_donor_cases=one,
        two_donor_cases=two,
        cases_with_overlapping_recipient_charge=overlap,
        cases_with_coarse_self_loop=self_loop,
        worst_balance_native_residual=wb,
        worst_donor_marginal_native_residual=wd,
        worst_recipient_marginal_native_residual=wr,
        worst_fate_partition_native_residual=wf,
        worst_restricted_pushforward_native_residual=wrestrict,
        signed_good_good_recipient_mass=signed_q.good_recipient_mass,
        signed_good_bad_recipient_mass=signed_q.bad_recipient_mass,
        signed_good_coarse_self_loop_fraction=coarse_q.self_loop_mass / coarse_q.total_positive_work_mass,
        generic_two_donor_single_recipient_charge_verified=two_ok,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=75_000)
    parser.add_argument("--seed", type=int, default=2026081204)
    parser.add_argument("--outdir", type=Path, default=Path("results-cyclic-hard-cell-single-charge"))
    args = parser.parse_args()
    out = stress(samples=args.samples, seed=args.seed)
    cert = theorem_certificate()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "certificate.json").write_text(json.dumps({"theorem": cert, "stress": asdict(out)}, indent=2, sort_keys=True) + "\n")
    summary = f"""# Cyclic hard-cell single-charge quotient

Status: **{STATUS}**.

The certified cyclic donor kernel is pushed through a deterministic hard Fourier/helicity map without changing cause.  For donor cell `C`, recipient cell `D`, and inherited recipient fate `F`, the positive table `K(C,D,F)` has row marginal `pi_#dW-` and column/fate marginal the already-canonical `pi_#dW+`.  Incoming atoms from several donor cells therefore disintegrate one recipient charge; they are not several causal charges.

Stress: `{out.samples}` physical closed triads
- resolved / numerically unresolved: `{out.resolved_cases}` / `{out.numerically_unresolved_cases}`
- one-donor / two-donor resolved cases: `{out.one_donor_cases}` / `{out.two_donor_cases}`
- overlapping-recipient-charge cases: `{out.cases_with_overlapping_recipient_charge}`
- coarse-self-loop cases: `{out.cases_with_coarse_self_loop}`
- worst balance native residual: `{out.worst_balance_native_residual:.3e}`
- worst donor marginal native residual: `{out.worst_donor_marginal_native_residual:.3e}`
- worst recipient marginal native residual: `{out.worst_recipient_marginal_native_residual:.3e}`
- worst recipient fate-partition native residual: `{out.worst_fate_partition_native_residual:.3e}`
- worst restricted-donor pushforward native residual: `{out.worst_restricted_pushforward_native_residual:.3e}`
- signed-good good/bad recipient masses: `{out.signed_good_good_recipient_mass:.12g}` / `{out.signed_good_bad_recipient_mass:.12g}`
- maximal-coarsening self-loop fraction: `{out.signed_good_coarse_self_loop_fraction:.12g}`
- generic two-donor -> one recipient single-charge anti-theorem: `{out.generic_two_donor_single_recipient_charge_verified}`

A geometry-bad recipient remains on the existing `TRANSFER_WORK_LOSS` stage-zero recursion route, but its modal energy remains in Navier--Stokes and may participate later.  A geometry-good recipient remains only Young-eligible.  A coarse self-loop is same-time physical redistribution with zero recursion depth and no scale progress.  The theorem does not say negative work pays failed good work and does not introduce any between-time deposit/withdrawal matching.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
