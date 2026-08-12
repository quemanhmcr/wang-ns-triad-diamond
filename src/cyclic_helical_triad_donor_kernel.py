from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

from src.canonical_positive_edge_work_routing import ETA0
from src.continuum_helical_edge_measure_registration import (
    ContinuumHelicalEdgeIdentity,
    HelicalModeIdentity,
    continuum_helical_edge_identity,
    unitary_fourier_convolution_factor,
)
from src.helical import coupling_g, stable_norm3
from src.helical_physical_edge_registration import (
    HelicalPhysicalEdgeRegistration,
    register_helical_physical_edge,
)

STATUS = (
    "EXACT_CYCLIC_HELICAL_TRIAD_DONOR_KERNEL__FULL_S3_CLOSED_TRIAD_QUOTIENT__"
    "THREE_ROOT_SIGNED_NS_ENERGY_CONSERVATION__CANONICAL_DW_MINUS_TO_DW_PLUS_MARGINALS__"
    "SIGNED_GOOD_UNIQUE_DONOR_AND_NONFORWARD_SIDE_RECIPIENT"
)

PARENT_RATIO_LO = Fraction(3, 5)
PARENT_RATIO_HI = Fraction(5, 8)
SIDE_TO_CHILD_LO = Fraction(3, 10)
SIDE_TO_CHILD_HI = Fraction(1, 3)
CHILD_TO_DONOR_LO = Fraction(3, 4)
CHILD_TO_DONOR_HI = Fraction(10, 13)
SIDE_TO_DONOR_LO = Fraction(3, 13)
SIDE_TO_DONOR_HI = Fraction(1, 4)
SAME_HELICITY_MULTIPLIER_UPPER = Fraction(9, 1600)
CLOSED_TRIAD_PERMUTATION_FACTOR = Fraction(1, 6)
EDGE_PARENT_SWAP_FACTOR = Fraction(1, 2)
SUM_RELATIVE_JACOBIAN = Fraction(1, 8)
CLOSED_TRIAD_ROOT_CHART_DENSITY = Fraction(1, 48)
EDGE_ROOT_CHART_DENSITY = Fraction(1, 16)


def _v3(value: Sequence[float], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must be one finite three-vector")
    if stable_norm3(out) == 0.0:
        raise ValueError(f"{name} must be nonzero")
    return out


def _complex(value: complex, name: str) -> complex:
    z = complex(value)
    if not math.isfinite(z.real) or not math.isfinite(z.imag):
        raise ValueError(f"finite {name} required")
    return z


def _relative_gap(actual: float, expected: float, scale: float) -> float:
    s = abs(float(scale))
    gap = abs(float(actual) - float(expected))
    if s == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / s


def _tuple_neg(v: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(-float(x) for x in v)


def cyclic_sum_relative_reroot(
    child: Sequence[float], relative_parent: Sequence[float]
) -> tuple[np.ndarray, np.ndarray]:
    """One cyclic re-root of the same closed triad in ``(z,r)`` coordinates.

    For the edge ``x+y=z`` with ``r=x-y``, re-root the closed triad at ``-x``
    with parents ``y,-z``.  Then

      z' = -(z+r)/2,
      r' =  (3z-r)/2.

    The corresponding 6D linear map has determinant +1.  This is a physical
    relabelling of the same closed triad, not a new event or scale step.
    """
    z = _v3(child, "child")
    r = np.asarray(relative_parent, dtype=float)
    if r.shape != (3,) or np.any(~np.isfinite(r)):
        raise ValueError("relative parent coordinate must be finite three-vector")
    zp = -0.5 * (z + r)
    rp = 0.5 * (3.0 * z - r)
    return zp, rp


def closed_triad_radon_certificate() -> dict[str, object]:
    """Exact finite-group multiplicity connecting closed triads to edge roots.

    Ordered closed triples ``k0+k1+k2=0`` are quotient by the full ``S3`` action,
    giving factor ``1/6`` a.e.  Marking one of the three root modes gives factor
    ``3/6=1/2``, exactly the existing unordered-parent edge quotient.  Passing
    from ordered parents to ``(z,r)`` contributes the already-certified ``1/8``
    inverse Jacobian, hence ``1/48`` per unmarked closed-triad root chart and
    ``3*(1/48)=1/16`` for the rooted edge law.
    """
    if 3 * CLOSED_TRIAD_PERMUTATION_FACTOR != EDGE_PARENT_SWAP_FACTOR:
        raise AssertionError("three closed-triad root marks do not recover the parent-swap quotient")
    if CLOSED_TRIAD_PERMUTATION_FACTOR * SUM_RELATIVE_JACOBIAN != CLOSED_TRIAD_ROOT_CHART_DENSITY:
        raise AssertionError("closed-triad root-chart density changed")
    if 3 * CLOSED_TRIAD_ROOT_CHART_DENSITY != EDGE_ROOT_CHART_DENSITY:
        raise AssertionError("closed-triad root marks do not recover the canonical 1/16 edge density")
    # One scalar-coordinate cyclic map is [[-1/2,-1/2],[3/2,-1/2]], det=1.
    det = Fraction(1, 4) + Fraction(3, 4)
    if det != 1:
        raise AssertionError("cyclic sum-relative re-root lost unit Jacobian")
    return {
        "ordered_closed_triad_S3_quotient": "1/6",
        "root_marks": 3,
        "rooted_parent_swap_factor": "1/2",
        "sum_relative_inverse_jacobian": "1/8",
        "closed_triad_root_chart_density": "1/48",
        "canonical_edge_density_after_three_roots": "1/16",
        "cyclic_reroot_abs_jacobian": "1",
        "fixed_loci": "S3 fixed sets require coincident wavevectors and are codimension-3 Lebesgue-null in the regular continuum",
        "reality_negation_quotiented": False,
        "radon": True,
    }


@dataclass(frozen=True)
class CyclicTriadSlotWork:
    closed_mode_index: int
    closed_mode: HelicalModeIdentity
    edge_identity: ContinuumHelicalEdgeIdentity
    parent_closed_indices: tuple[int, int]
    edge_registration: HelicalPhysicalEdgeRegistration
    expected_signed_work: float

    def __post_init__(self) -> None:
        if self.closed_mode_index not in (0, 1, 2):
            raise ValueError("closed triad root index must be 0,1,2")
        if len(self.parent_closed_indices) != 2 or self.closed_mode_index in self.parent_closed_indices:
            raise ValueError("cyclic root must use the other two closed modes as parents")
        if set((self.closed_mode_index, *self.parent_closed_indices)) != {0, 1, 2}:
            raise ValueError("cyclic root indices must resolve the full closed triad")
        actual = float(self.edge_registration.signed_child_energy_work)
        scale = max(
            abs(actual),
            abs(float(self.expected_signed_work)),
            float(self.edge_registration.native_modal_capacity),
            1.0e-300,
        )
        if _relative_gap(actual, self.expected_signed_work, scale) > 3.0e-10:
            raise AssertionError("cyclic root work disagrees with the common Waleffe triad factor")
        expected_child = _tuple_neg(self.closed_mode.wavevector)
        if self.edge_identity.child.wavevector != expected_child or self.edge_identity.child.helicity != self.closed_mode.helicity:
            raise AssertionError("cyclic root edge child is not the reality partner of its closed mode")

    @property
    def signed_work(self) -> float:
        return float(self.edge_registration.signed_child_energy_work)

    @property
    def positive_work(self) -> float:
        return max(self.signed_work, 0.0)

    @property
    def negative_work(self) -> float:
        return max(-self.signed_work, 0.0)

    @property
    def signed_efficiency(self) -> float:
        return float(self.edge_registration.normalized_multiplier * self.edge_registration.phase_alignment)


@dataclass(frozen=True)
class DonorRecipientFlow:
    donor_closed_mode_index: int
    recipient_closed_mode_index: int
    donor_child_mode: HelicalModeIdentity
    recipient_child_mode: HelicalModeIdentity
    physical_work: float

    def __post_init__(self) -> None:
        if self.donor_closed_mode_index == self.recipient_closed_mode_index:
            raise ValueError("donor and recipient must be different modes of the closed triad")
        if self.donor_closed_mode_index not in (0, 1, 2) or self.recipient_closed_mode_index not in (0, 1, 2):
            raise ValueError("donor/recipient closed-mode indices must be 0,1,2")
        if not math.isfinite(self.physical_work) or self.physical_work <= 0.0:
            raise ValueError("donor-recipient atom must carry positive physical work")


@dataclass(frozen=True)
class CyclicTriadDonorKernel:
    flows: tuple[DonorRecipientFlow, ...]
    total_positive_work: float
    total_negative_work: float
    donor_marginal_residual: float
    recipient_marginal_residual: float
    donor_count: int
    recipient_count: int
    transport_unique: bool
    canonical_positive_law_replaced: bool = False
    capacity_used: bool = False
    creates_new_event: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        vals = (
            self.total_positive_work,
            self.total_negative_work,
            self.donor_marginal_residual,
            self.recipient_marginal_residual,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite cyclic donor-kernel data required")
        if min(self.total_positive_work, self.total_negative_work) < 0.0:
            raise ValueError("nonnegative positive/negative work totals required")
        if self.canonical_positive_law_replaced or self.capacity_used or self.creates_new_event or self.later_hahn_used:
            raise ValueError("cyclic donor kernel may only add same-triad provenance to the canonical work law")
        if not self.transport_unique:
            raise ValueError("three-slot donor transport must be the unique work-marginal coupling")
        scale = max(self.total_positive_work, self.total_negative_work, 1.0e-300)
        if abs(self.total_positive_work - self.total_negative_work) > 3.0e-10 * scale:
            raise AssertionError("closed triad did not conserve nonlinear energy before Hahn splitting")
        if self.donor_marginal_residual > 4.0e-10 or self.recipient_marginal_residual > 4.0e-10:
            raise AssertionError("donor kernel failed the canonical dW-/dW+ marginals")
        if self.total_positive_work == 0.0:
            if self.flows or self.donor_count or self.recipient_count:
                raise AssertionError("zero-work triad cannot carry donor-recipient atoms")
        elif not (1 <= self.donor_count <= 2 and 1 <= self.recipient_count <= 2):
            raise AssertionError("nonzero three-slot conservative triad must have one/two donors and recipients")


@dataclass(frozen=True)
class DonorRecipientMeasureAtom:
    donor_closed_mode_index: int
    recipient_closed_mode_index: int
    donor_child_mode: HelicalModeIdentity
    recipient_child_mode: HelicalModeIdentity
    physical_work_mass: float

    def __post_init__(self) -> None:
        if self.donor_closed_mode_index == self.recipient_closed_mode_index:
            raise ValueError("measure donor and recipient must be distinct cyclic roots")
        if not math.isfinite(self.physical_work_mass) or self.physical_work_mass <= 0.0:
            raise ValueError("positive finite donor-recipient physical work mass required")


@dataclass(frozen=True)
class CyclicTriadMeasureKernel:
    quotient_measure_mass: float
    donor_edge_negative_masses: tuple[float, float, float]
    recipient_edge_positive_masses: tuple[float, float, float]
    atoms: tuple[DonorRecipientMeasureAtom, ...]
    total_mass: float
    donor_marginal_residual: float
    recipient_marginal_residual: float
    recipient_is_submeasure_of_canonical_dW_plus: bool = True
    donor_is_restriction_of_canonical_dW_minus: bool = True
    canonical_dW_plus_replaced: bool = False
    creates_new_event: bool = False

    def __post_init__(self) -> None:
        if not math.isfinite(self.quotient_measure_mass) or self.quotient_measure_mass < 0.0:
            raise ValueError("finite nonnegative closed-triad quotient measure mass required")
        vals = (*self.donor_edge_negative_masses, *self.recipient_edge_positive_masses, self.total_mass)
        if not all(math.isfinite(float(v)) and float(v) >= 0.0 for v in vals):
            raise ValueError("finite nonnegative cyclic work masses required")
        if (
            not self.recipient_is_submeasure_of_canonical_dW_plus
            or not self.donor_is_restriction_of_canonical_dW_minus
            or self.canonical_dW_plus_replaced
            or self.creates_new_event
        ):
            raise ValueError("measure kernel may only transport restricted dW- into inherited canonical dW+ provenance")
        dtotal = math.fsum(self.donor_edge_negative_masses)
        ptotal = math.fsum(self.recipient_edge_positive_masses)
        scale = max(dtotal, ptotal, self.total_mass, 1.0e-300)
        if max(abs(dtotal-self.total_mass), abs(ptotal-self.total_mass)) > 4.0e-10*scale:
            raise AssertionError("cyclic measure kernel changed total physical work mass")
        if self.donor_marginal_residual > 4.0e-10 or self.recipient_marginal_residual > 4.0e-10:
            raise AssertionError("cyclic measure kernel lost its dW-/dW+ marginals")


@dataclass(frozen=True)
class RestrictedNegativeRecipientPushforward:
    selected_donor_closed_mode_indices: tuple[int, ...]
    selected_negative_mass: float
    recipient_masses: tuple[float, float, float]
    recipient_total_mass: float
    recipient_dominated_by_full_canonical_positive_mass: tuple[bool, bool, bool]
    creates_new_event: bool = False
    capacity_used: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        if not self.selected_donor_closed_mode_indices:
            raise ValueError("nonempty negative donor restriction required")
        if len(set(self.selected_donor_closed_mode_indices)) != len(self.selected_donor_closed_mode_indices):
            raise ValueError("donor restriction contains duplicate cyclic roots")
        if any(i not in (0,1,2) for i in self.selected_donor_closed_mode_indices):
            raise ValueError("donor restriction indices must lie in {0,1,2}")
        if not math.isfinite(self.selected_negative_mass) or self.selected_negative_mass <= 0.0:
            raise ValueError("positive selected canonical negative-work mass required")
        if any(not math.isfinite(v) or v < 0.0 for v in self.recipient_masses):
            raise ValueError("finite nonnegative recipient pushforward masses required")
        if abs(math.fsum(self.recipient_masses)-self.recipient_total_mass) > 4.0e-10*max(self.recipient_total_mass,1.0e-300):
            raise AssertionError("restricted recipient pushforward total changed")
        if abs(self.recipient_total_mass-self.selected_negative_mass) > 4.0e-10*max(self.selected_negative_mass,1.0e-300):
            raise AssertionError("restricted dW- mass was not preserved by same-triad positive pushforward")
        if not all(self.recipient_dominated_by_full_canonical_positive_mass):
            raise AssertionError("restricted negative-work recipient law exceeded canonical dW+ on a cyclic root")
        if self.creates_new_event or self.capacity_used or self.later_hahn_used:
            raise ValueError("restricted negative-work pushforward changed physical event/cause semantics")


@dataclass(frozen=True)
class ClosedHelicalTriadRegistration:
    modes: tuple[HelicalModeIdentity, HelicalModeIdentity, HelicalModeIdentity]
    amplitudes: tuple[complex, complex, complex]
    common_waleffe_coupling: complex
    common_phase_work_factor: float
    slots: tuple[CyclicTriadSlotWork, CyclicTriadSlotWork, CyclicTriadSlotWork]
    signed_energy_conservation_residual: float
    cyclic_coupling_residual: float
    donor_kernel: CyclicTriadDonorKernel
    parent_permutation_quotiented: bool = True
    reality_negation_quotiented: bool = False

    def __post_init__(self) -> None:
        if len(self.modes) != 3 or len(self.amplitudes) != 3 or len(self.slots) != 3:
            raise ValueError("one closed helical triad requires exactly three modes/amplitudes/root slots")
        if tuple(sorted(self.modes)) != self.modes:
            raise ValueError("closed triad modes must use the canonical full-S3 storage order")
        if not self.parent_permutation_quotiented or self.reality_negation_quotiented:
            raise ValueError("closed triad quotient removes S3 ordering only, not global reality negation")
        vec_sum = np.sum(np.asarray([m.wavevector for m in self.modes], dtype=float), axis=0)
        scale = max(stable_norm3(np.asarray(m.wavevector, float)) for m in self.modes)
        if stable_norm3(vec_sum) > 3.0e-12 * scale:
            raise AssertionError("stored closed helical triad no longer sums to zero")
        works = tuple(slot.signed_work for slot in self.slots)
        variation = math.fsum(abs(w) for w in works)
        expected_res = 0.0 if variation == 0.0 else abs(math.fsum(works)) / variation
        if abs(self.signed_energy_conservation_residual - expected_res) > 3.0e-14 * max(1.0, expected_res):
            raise AssertionError("stored cyclic energy-conservation residual changed")
        if self.signed_energy_conservation_residual > 4.0e-10:
            raise AssertionError("one helical closed triad lost exact nonlinear energy conservation")
        if self.cyclic_coupling_residual > 4.0e-10:
            raise AssertionError("cyclic Waleffe triple-product coupling changed with the root")

    def slot_for_closed_mode_index(self, index: int) -> CyclicTriadSlotWork:
        for slot in self.slots:
            if slot.closed_mode_index == int(index):
                return slot
        raise KeyError(index)

    def slot_for_edge_child_wavevector(self, wavevector: Sequence[float]) -> CyclicTriadSlotWork:
        target = tuple(float(v) for v in _v3(wavevector, "edge child lookup"))
        for slot in self.slots:
            if slot.edge_identity.child.wavevector == target:
                return slot
        raise KeyError(target)


def _canonical_closed_data(
    wavevectors: Sequence[Sequence[float]],
    helicities: Sequence[int],
    amplitudes: Sequence[complex],
) -> tuple[tuple[HelicalModeIdentity, ...], tuple[complex, ...]]:
    if len(wavevectors) != 3 or len(helicities) != 3 or len(amplitudes) != 3:
        raise ValueError("closed helical triad requires three wavevectors, helicities and amplitudes")
    rows: list[tuple[HelicalModeIdentity, complex]] = []
    vectors = []
    for i, (k, s, a) in enumerate(zip(wavevectors, helicities, amplitudes)):
        kv = _v3(k, f"closed wavevector {i}")
        vectors.append(kv)
        rows.append((HelicalModeIdentity(tuple(float(v) for v in kv), int(s)), _complex(a, f"amplitude {i}")))
    scale = max(stable_norm3(k) for k in vectors)
    if stable_norm3(vectors[0] + vectors[1] + vectors[2]) > 3.0e-12 * scale:
        raise ValueError("closed helical wavevectors must sum to zero")
    if len({row[0].wavevector for row in rows}) != 3:
        raise ValueError("regular continuum closed triad requires three distinct wavevectors")
    rows.sort(key=lambda row: row[0])
    return tuple(row[0] for row in rows), tuple(row[1] for row in rows)


def _build_donor_kernel(slots: Sequence[CyclicTriadSlotWork]) -> CyclicTriadDonorKernel:
    positive = {slot.closed_mode_index: slot.positive_work for slot in slots if slot.positive_work > 0.0}
    negative = {slot.closed_mode_index: slot.negative_work for slot in slots if slot.negative_work > 0.0}
    total_p = math.fsum(positive.values())
    total_n = math.fsum(negative.values())
    scale = max(total_p, total_n, math.fsum(abs(slot.signed_work) for slot in slots), 1.0e-300)
    if abs(total_p - total_n) > 4.0e-10 * scale:
        raise AssertionError("cannot construct donor kernel before exact closed-triad energy conservation")
    if total_p == 0.0:
        return CyclicTriadDonorKernel((), 0.0, 0.0, 0.0, 0.0, 0, 0, True)
    if min(len(positive), len(negative)) != 1:
        raise AssertionError("a three-slot zero-sum work vector must have a singleton donor or recipient side")
    slots_by_index = {slot.closed_mode_index: slot for slot in slots}
    flows: list[DonorRecipientFlow] = []
    for donor, nwork in negative.items():
        for recipient, pwork in positive.items():
            work = nwork * pwork / total_p
            dslot = slots_by_index[donor]
            rslot = slots_by_index[recipient]
            flows.append(
                DonorRecipientFlow(
                    donor_closed_mode_index=donor,
                    recipient_closed_mode_index=recipient,
                    donor_child_mode=dslot.edge_identity.child,
                    recipient_child_mode=rslot.edge_identity.child,
                    physical_work=work,
                )
            )
    row_res = 0.0
    for donor, expected in negative.items():
        actual = math.fsum(flow.physical_work for flow in flows if flow.donor_closed_mode_index == donor)
        row_res = max(row_res, _relative_gap(actual, expected, max(expected, total_p)))
    col_res = 0.0
    for recipient, expected in positive.items():
        actual = math.fsum(flow.physical_work for flow in flows if flow.recipient_closed_mode_index == recipient)
        col_res = max(col_res, _relative_gap(actual, expected, max(expected, total_p)))
    return CyclicTriadDonorKernel(
        flows=tuple(flows),
        total_positive_work=total_p,
        total_negative_work=total_n,
        donor_marginal_residual=row_res,
        recipient_marginal_residual=col_res,
        donor_count=len(negative),
        recipient_count=len(positive),
        transport_unique=True,
    )


def register_closed_helical_triad(
    *,
    wavevectors: Sequence[Sequence[float]],
    helicities: Sequence[int],
    amplitudes: Sequence[complex],
) -> ClosedHelicalTriadRegistration:
    """Register all three cyclic child-work roots of one real-field closed triad.

    Stored closed modes satisfy ``k0+k1+k2=0``.  Root ``i`` is the physical edge
    with child ``-k_i`` and parents ``k_j,k_k``.  Reality gives child amplitude
    ``conj(a_i)`` and keeps the same helicity sign because the repository basis
    obeys ``h_s(-k)=conj(h_s(k))``.
    """
    modes, amps = _canonical_closed_data(wavevectors, helicities, amplitudes)
    k = tuple(np.asarray(mode.wavevector, dtype=float) for mode in modes)
    s = tuple(mode.helicity for mode in modes)
    g0 = coupling_g(k[1], k[2], k[0], s[1], s[2], s[0])
    cyclic_g = (
        g0,
        coupling_g(k[2], k[0], k[1], s[2], s[0], s[1]),
        coupling_g(k[0], k[1], k[2], s[0], s[1], s[2]),
    )
    gscale = max(abs(x) for x in cyclic_g)
    cyclic_res = 0.0 if gscale == 0.0 else max(abs(x - g0) for x in cyclic_g) / gscale
    if cyclic_res > 3.0e-10:
        raise AssertionError("Waleffe coupling is not cyclically invariant on the closed triad")
    common = 4.0 * float(np.real(amps[0] * amps[1] * amps[2] * np.conjugate(g0)))
    slots: list[CyclicTriadSlotWork] = []
    for i in range(3):
        j = (i + 1) % 3
        ell = (i + 2) % 3
        child = -k[i]
        reg = register_helical_physical_edge(
            x=k[j],
            y=k[ell],
            z=child,
            sx=s[j],
            sy=s[ell],
            sz=s[i],
            ax=amps[j],
            ay=amps[ell],
            az=np.conjugate(amps[i]),
        )
        edge_identity = continuum_helical_edge_identity(k[j], k[ell], child, s[j], s[ell], s[i])
        coeff = s[j] * stable_norm3(k[j]) - s[ell] * stable_norm3(k[ell])
        slots.append(
            CyclicTriadSlotWork(
                closed_mode_index=i,
                closed_mode=modes[i],
                edge_identity=edge_identity,
                parent_closed_indices=(j, ell),
                edge_registration=reg,
                expected_signed_work=coeff * common,
            )
        )
    variation = math.fsum(abs(slot.signed_work) for slot in slots)
    conservation = 0.0 if variation == 0.0 else abs(math.fsum(slot.signed_work for slot in slots)) / variation
    kernel = _build_donor_kernel(slots)
    return ClosedHelicalTriadRegistration(
        modes=(modes[0], modes[1], modes[2]),
        amplitudes=(amps[0], amps[1], amps[2]),
        common_waleffe_coupling=complex(g0),
        common_phase_work_factor=common,
        slots=(slots[0], slots[1], slots[2]),
        signed_energy_conservation_residual=conservation,
        cyclic_coupling_residual=cyclic_res,
        donor_kernel=kernel,
    )


def cyclic_triad_measure_kernel(
    triad: ClosedHelicalTriadRegistration, *, quotient_measure_mass: float
) -> CyclicTriadMeasureKernel:
    """Lift one density-level donor kernel to the physical Fourier work measure.

    ``quotient_measure_mass`` is mass of the full-S3 closed-triad base law.  The
    unitary Fourier factor is the same ``C_F`` already used by canonical edge
    registration.  Root marking is handled by the exact ``3*(1/48)=1/16``
    quotient theorem; no extra multiplicity is inserted here.
    """
    q = float(quotient_measure_mass)
    if not math.isfinite(q) or q < 0.0:
        raise ValueError("finite nonnegative closed-triad quotient measure mass required")
    factor = unitary_fourier_convolution_factor() * q
    donor = tuple(factor*slot.negative_work for slot in triad.slots)
    recipient = tuple(factor*slot.positive_work for slot in triad.slots)
    atoms = tuple(
        DonorRecipientMeasureAtom(
            donor_closed_mode_index=flow.donor_closed_mode_index,
            recipient_closed_mode_index=flow.recipient_closed_mode_index,
            donor_child_mode=flow.donor_child_mode,
            recipient_child_mode=flow.recipient_child_mode,
            physical_work_mass=factor*flow.physical_work,
        )
        for flow in triad.donor_kernel.flows
        if factor*flow.physical_work > 0.0
    )
    total = math.fsum(atom.physical_work_mass for atom in atoms)
    row_res = 0.0
    for i, expected in enumerate(donor):
        if expected <= 0.0:
            continue
        actual = math.fsum(atom.physical_work_mass for atom in atoms if atom.donor_closed_mode_index==i)
        row_res=max(row_res,_relative_gap(actual,expected,max(total,expected)))
    col_res = 0.0
    for i, expected in enumerate(recipient):
        if expected <= 0.0:
            continue
        actual = math.fsum(atom.physical_work_mass for atom in atoms if atom.recipient_closed_mode_index==i)
        col_res=max(col_res,_relative_gap(actual,expected,max(total,expected)))
    return CyclicTriadMeasureKernel(
        quotient_measure_mass=q,
        donor_edge_negative_masses=(donor[0],donor[1],donor[2]),
        recipient_edge_positive_masses=(recipient[0],recipient[1],recipient[2]),
        atoms=atoms,
        total_mass=total,
        donor_marginal_residual=row_res,
        recipient_marginal_residual=col_res,
    )


def pushforward_restricted_negative_work(
    kernel: CyclicTriadMeasureKernel, *, donor_closed_mode_indices: Sequence[int]
) -> RestrictedNegativeRecipientPushforward:
    """Push one canonical ``dW-`` restriction to same-triad positive recipients.

    The output is a positive submeasure of the already-existing canonical ``dW+``
    recipient law.  It is not a new Hahn decomposition and does not change event
    time.  This is the interface needed to read a hard-cell ``n_C`` obstruction
    as physical donor provenance rather than a cancellation budget.
    """
    selected=tuple(sorted(set(int(i) for i in donor_closed_mode_indices)))
    if not selected or any(i not in (0,1,2) for i in selected):
        raise ValueError("nonempty donor restriction inside {0,1,2} required")
    selected_mass=math.fsum(kernel.donor_edge_negative_masses[i] for i in selected)
    if selected_mass<=0.0:
        raise ValueError("selected cyclic roots carry no canonical negative work")
    recipient=[]
    dominated=[]
    for j in range(3):
        mass=math.fsum(
            atom.physical_work_mass
            for atom in kernel.atoms
            if atom.donor_closed_mode_index in selected and atom.recipient_closed_mode_index==j
        )
        recipient.append(mass)
        full=kernel.recipient_edge_positive_masses[j]
        dominated.append(mass <= full + 4.0e-10*max(full,mass,1.0e-300))
    return RestrictedNegativeRecipientPushforward(
        selected_donor_closed_mode_indices=selected,
        selected_negative_mass=selected_mass,
        recipient_masses=(recipient[0],recipient[1],recipient[2]),
        recipient_total_mass=math.fsum(recipient),
        recipient_dominated_by_full_canonical_positive_mass=(dominated[0],dominated[1],dominated[2]),
    )


@dataclass(frozen=True)
class SignedGoodSideRecipientCertificate:
    recipient_closed_mode_index: int
    energy_donor_closed_mode_index: int
    side_recipient_closed_mode_index: int
    recipient_work: float
    donor_negative_work: float
    side_positive_work: float
    side_to_recipient_ratio: float
    recipient_to_donor_ratio: float
    side_to_donor_ratio: float
    recipient_signed_efficiency: float
    parent_ratio_min: float
    parent_ratio_max: float
    side_forward_ratio: float
    side_geometric_multiplier: float
    parent_helicities_opposite: bool
    unique_energy_donor: bool
    side_is_positive_nonforward: bool
    interaction_parents_remain_two: bool = True
    side_terminal_transfer_loss_is_existing_router_consequence: bool = True
    creates_new_event: bool = False

    def __post_init__(self) -> None:
        if min(self.recipient_work, self.donor_negative_work, self.side_positive_work) <= 0.0:
            raise ValueError("signed-good donor/recipient certificate requires nonzero physical works")
        if not self.parent_helicities_opposite or not self.unique_energy_donor or not self.side_is_positive_nonforward:
            raise AssertionError("signed-good cyclic triad lost its rigid donor/side-recipient structure")
        if not self.interaction_parents_remain_two or not self.side_terminal_transfer_loss_is_existing_router_consequence:
            raise ValueError("energy-donor provenance may not erase the quadratic two-parent interaction")
        if self.creates_new_event:
            raise ValueError("cyclic donor/side-recipient registration is same-time provenance, not a new event")
        if not (float(SIDE_TO_CHILD_LO) < self.side_to_recipient_ratio < float(SIDE_TO_CHILD_HI)):
            raise AssertionError("signed-good side/recipient work ratio left the exact 3/10..1/3 interval")
        if not (float(CHILD_TO_DONOR_LO) < self.recipient_to_donor_ratio < float(CHILD_TO_DONOR_HI)):
            raise AssertionError("signed-good recipient/donor work ratio left the exact 3/4..10/13 interval")
        if not (float(SIDE_TO_DONOR_LO) < self.side_to_donor_ratio < float(SIDE_TO_DONOR_HI)):
            raise AssertionError("signed-good side/donor work ratio left the exact 3/13..1/4 interval")
        if not (float(PARENT_RATIO_LO) < self.parent_ratio_min <= self.parent_ratio_max < float(PARENT_RATIO_HI)):
            raise AssertionError("signed-good parent/child scale window changed")
        if self.side_forward_ratio >= 1.0 or self.side_geometric_multiplier != 0.0:
            raise AssertionError("side recipient is not the positive nonforward cyclic edge")


def signed_good_side_recipient_certificate(
    triad: ClosedHelicalTriadRegistration,
    *,
    recipient_closed_mode_index: int,
) -> SignedGoodSideRecipientCertificate:
    """Rigid cyclic energy geometry of a signed-good positive forward recipient.

    This theorem distinguishes *interaction parents* from *energy donors*.  Young
    and the quadratic source still use both parents.  On the signed-good forward
    core, however, the opposite-parent-helicity theorem makes exactly one of them
    an energy donor; the other is a simultaneous positive side recipient.
    """
    slot = triad.slot_for_closed_mode_index(recipient_closed_mode_index)
    row = slot.edge_registration
    efficiency = slot.signed_efficiency
    if not (slot.signed_work > 0.0 and efficiency > 1.0 - ETA0 and row.scale_progress > 0.0):
        raise ValueError("side-recipient theorem requires one actual signed-good positive forward root")
    child = float(row.child_frequency)
    ratios = (float(row.parent_x_frequency) / child, float(row.parent_y_frequency) / child)
    if not (float(PARENT_RATIO_LO) < min(ratios) <= max(ratios) < float(PARENT_RATIO_HI)):
        raise AssertionError("existing signed-good theorem no longer supplies the 3/5..5/8 parent window")
    if row.parent_x_helicity == row.parent_y_helicity:
        # For same parent helicities the exact sign factor ratio to the
        # sign-exhausted envelope is d(1+s)/[s(1+d)].  Signed-good gives
        # u=log(y/x)<=1/200 and the clean parent window gives s<5/4, so
        # d/s=tanh(u/2)<1/400 and the ratio is <(9/4)/400=9/1600.
        raise AssertionError(
            "signed-good edge cannot have same-helicity parents: J/J* would be <=9/1600"
        )
    if row.child_helicity not in (row.parent_x_helicity, row.parent_y_helicity):
        raise AssertionError("binary helicity signs lost the unique child-matching parent")
    parent_indices = slot.parent_closed_indices
    if row.parent_x_helicity == row.child_helicity:
        donor_index, side_index = parent_indices[0], parent_indices[1]
    else:
        donor_index, side_index = parent_indices[1], parent_indices[0]
    donor_slot = triad.slot_for_closed_mode_index(donor_index)
    side_slot = triad.slot_for_closed_mode_index(side_index)
    recipient_work = slot.signed_work
    donor_negative = donor_slot.negative_work
    side_positive = side_slot.positive_work
    if donor_negative <= 0.0 or side_positive <= 0.0:
        raise AssertionError("signed-good forward recipient did not produce one donor and one side recipient")
    scale = max(donor_negative, recipient_work + side_positive, 1.0e-300)
    if abs(donor_negative - recipient_work - side_positive) > 4.0e-10 * scale:
        raise AssertionError("signed-good donor loss did not split into child plus side-recipient work")
    side_row = side_slot.edge_registration
    if not (side_slot.signed_work > 0.0 and side_row.scale_progress == 0.0 and side_row.geometric_multiplier_J == 0.0):
        raise AssertionError("cyclic side recipient is not positive nonforward physical work")
    side_ratio = side_positive / recipient_work
    child_donor = recipient_work / donor_negative
    side_donor = side_positive / donor_negative
    return SignedGoodSideRecipientCertificate(
        recipient_closed_mode_index=recipient_closed_mode_index,
        energy_donor_closed_mode_index=donor_index,
        side_recipient_closed_mode_index=side_index,
        recipient_work=recipient_work,
        donor_negative_work=donor_negative,
        side_positive_work=side_positive,
        side_to_recipient_ratio=side_ratio,
        recipient_to_donor_ratio=child_donor,
        side_to_donor_ratio=side_donor,
        recipient_signed_efficiency=efficiency,
        parent_ratio_min=min(ratios),
        parent_ratio_max=max(ratios),
        side_forward_ratio=float(side_row.forward_ratio),
        side_geometric_multiplier=float(side_row.geometric_multiplier_J),
        parent_helicities_opposite=True,
        unique_energy_donor=True,
        side_is_positive_nonforward=True,
    )


def translate_closed_amplitudes(
    wavevectors: Sequence[Sequence[float]], amplitudes: Sequence[complex], translation: Sequence[float]
) -> tuple[complex, complex, complex]:
    if len(wavevectors) != 3 or len(amplitudes) != 3:
        raise ValueError("three closed modes/amplitudes required")
    x0 = np.asarray(translation, dtype=float)
    if x0.shape != (3,) or np.any(~np.isfinite(x0)):
        raise ValueError("finite translation vector required")
    return tuple(
        _complex(a, "translation amplitude") * np.exp(1j * float(np.dot(_v3(k, "translation wavevector"), x0)))
        for k, a in zip(wavevectors, amplitudes)
    )


def global_reality_negation(
    wavevectors: Sequence[Sequence[float]], amplitudes: Sequence[complex]
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], tuple[complex, complex, complex]]:
    if len(wavevectors) != 3 or len(amplitudes) != 3:
        raise ValueError("three closed modes/amplitudes required")
    return (
        tuple(-_v3(k, "reality wavevector") for k in wavevectors),
        tuple(np.conjugate(_complex(a, "reality amplitude")) for a in amplitudes),
    )


def generic_two_donor_counterexample() -> ClosedHelicalTriadRegistration:
    """Exact anti-theorem: generic positive child work need not have one donor."""
    k = (
        np.array([-1.0, -2.0, 0.0]),
        np.array([0.0, 2.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
    )
    s = (1, 1, 1)
    g = coupling_g(k[1], k[2], k[0], 1, 1, 1)
    # Choose product phase so the common factor is negative.  With canonical
    # magnitudes sqrt(5),2,1 this flips the difference cycle into two donors.
    a0 = -g / abs(g) if abs(g) > 0.0 else -1.0
    triad = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=(a0, 1.0, 1.0))
    if triad.donor_kernel.donor_count != 2 or triad.donor_kernel.recipient_count != 1:
        raise AssertionError("generic two-donor counterexample did not realize two donors")
    return triad


def signed_good_integer_triad() -> tuple[ClosedHelicalTriadRegistration, SignedGoodSideRecipientCertificate]:
    """The cutoff-7 integer near-extremal triad used by the evolved NS audit."""
    child = np.array([7.0, 6.0, 5.0])
    x = np.array([5.0, 0.0, 4.0])
    y = np.array([2.0, 6.0, 1.0])
    sx, sy, sz = 1, -1, 1
    g = coupling_g(x, y, -child, sx, sy, sz)
    ax = g / abs(g)
    # Closed modes are (-child,x,y).  The coefficient at -child is conj(az)=1.
    triad = register_closed_helical_triad(
        wavevectors=(-child, x, y),
        helicities=(sz, sx, sy),
        amplitudes=(1.0, ax, 1.0),
    )
    recipient_slot = triad.slot_for_edge_child_wavevector(child)
    cert = signed_good_side_recipient_certificate(
        triad, recipient_closed_mode_index=recipient_slot.closed_mode_index
    )
    return triad, cert


def theorem_certificate() -> dict[str, object]:
    quotient = closed_triad_radon_certificate()
    return {
        "status": STATUS,
        "closed_triad_quotient": quotient,
        "cyclic_work": "for k0+k1+k2=0, root i has child -ki and work T_i=(s_j|k_j|-s_l|k_l|) R with one common physical Waleffe phase factor R; therefore sum_i T_i=0 before Hahn",
        "donor_kernel": "Q=sum[T_i]_+=sum[-T_i]_+; M(i->j)=[-T_i]_+[T_j]_+/Q.  On three slots the zero-sum sign pattern has a singleton donor or recipient side, so this is the unique positive transport with dW- donor and dW+ recipient marginals",
        "causal_semantics": "the recipient marginal is the already-canonical dW+ law; the kernel adds same-time same-triad donor provenance and neither replaces cause nor creates an event",
        "negative_restriction": "every measurable/replayed restriction of cyclic dW- pushes through the same positive kernel to a mass-preserving submeasure of canonical recipient dW+; this rereads cancellation as donor provenance, not as a new owner or budget",
        "generic_anti_theorem": "a generic triad may have two energy donors and one positive recipient; unique donor is not a generic HH law",
        "signed_good_helicity": f"same-helicity parents would have J/J*<={SAME_HELICITY_MULTIPLIER_UPPER}; signed-good efficiency >1-1e-4 therefore forces opposite parent helicities",
        "signed_good_energy_roles": "for a positive forward signed-good child, exactly one interaction parent shares the child helicity and is the unique energy donor; the other interaction parent is a simultaneous positive side recipient",
        "side_work": "with both parent/child ratios in (3/5,5/8), 3/10<W_side+/W_child+<1/3, 3/4<W_child+/W_donor-<10/13, and 3/13<W_side+/W_donor-<1/4",
        "side_fate": "the side cyclic edge has the original high child among its parents, hence forward ratio<1, J=0, and its positive work is the already-certified positive-nonforward TRANSFER_WORK_LOSS sublaw",
        "interaction_ontology": "two quadratic interaction parents remain two even when only one is an energy donor; donor provenance does not linearize HH",
        "coherent_scope": "this cyclic kernel is a physical closed-Fourier-triad disintegration only; it is not a general coherent POVM kernel and fresh coherent Hahn remains noncausal",
        "reality_scope": "global k->-k reality negation is covariant but is not quotiented from the canonical edge law",
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class CyclicDonorKernelStress:
    samples: int
    one_donor_cases: int
    two_donor_cases: int
    zero_work_cases: int
    worst_energy_conservation_relative: float
    worst_cyclic_coupling_relative: float
    worst_donor_marginal_relative: float
    worst_recipient_marginal_relative: float
    worst_measure_donor_marginal_relative: float
    worst_measure_recipient_marginal_relative: float
    worst_restricted_negative_mass_residual: float
    worst_permutation_work_residual: float
    worst_translation_work_residual: float
    worst_wavevector_scaling_residual: float
    worst_amplitude_cubic_scaling_residual: float
    worst_reality_work_multiset_residual: float
    signed_good_side_ratio: float
    signed_good_child_donor_ratio: float
    signed_good_side_donor_ratio: float
    generic_two_donor_counterexample_passed: bool


def _work_multiset(triad: ClosedHelicalTriadRegistration) -> tuple[float, float, float]:
    return tuple(sorted(slot.signed_work for slot in triad.slots))


def _multiset_relative(a: Sequence[float], b: Sequence[float]) -> float:
    aa, bb = tuple(sorted(float(x) for x in a)), tuple(sorted(float(x) for x in b))
    scale = max(math.fsum(abs(x) for x in aa), math.fsum(abs(x) for x in bb), 1.0e-300)
    return max(abs(x - y) for x, y in zip(aa, bb)) / scale


def stress(samples: int = 75_000, seed: int = 2026081203) -> CyclicDonorKernelStress:
    if samples <= 0:
        raise ValueError("positive stress sample count required")
    rng = np.random.default_rng(seed)
    one = two = zero = 0
    we = wc = wd = wr = wmd = wmr = wrestrict = wp = wt = ws = wa = wreal = 0.0
    for i in range(int(samples)):
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
        triad = register_closed_helical_triad(
            wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=amps
        )
        kernel = triad.donor_kernel
        we = max(we, triad.signed_energy_conservation_residual)
        wc = max(wc, triad.cyclic_coupling_residual)
        wd = max(wd, kernel.donor_marginal_residual)
        wr = max(wr, kernel.recipient_marginal_residual)
        qmass = math.exp(float(rng.uniform(-8.0, 8.0)))
        mkernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=qmass)
        wmd = max(wmd, mkernel.donor_marginal_residual)
        wmr = max(wmr, mkernel.recipient_marginal_residual)
        donors = tuple(i for i,mass in enumerate(mkernel.donor_edge_negative_masses) if mass>0.0)
        if donors:
            chosen = donors[:1] if len(donors)==1 or rng.random()<0.5 else donors
            restricted = pushforward_restricted_negative_work(mkernel, donor_closed_mode_indices=chosen)
            wrestrict = max(
                wrestrict,
                abs(restricted.recipient_total_mass-restricted.selected_negative_mass)
                / max(restricted.selected_negative_mass,1.0e-300),
            )
        if kernel.total_positive_work == 0.0:
            zero += 1
        elif kernel.donor_count == 1:
            one += 1
        elif kernel.donor_count == 2:
            two += 1
        else:
            raise AssertionError("three-slot donor count left {1,2}")

        if i % 12 == 0:
            perm = tuple(int(x) for x in rng.permutation(3))
            ptriad = register_closed_helical_triad(
                wavevectors=tuple((k0, k1, k2)[j] for j in perm),
                helicities=tuple(helicities[j] for j in perm),
                amplitudes=tuple(amps[j] for j in perm),
            )
            wp = max(wp, _multiset_relative(_work_multiset(triad), _work_multiset(ptriad)))

            x0 = rng.normal(size=3)
            tamps = translate_closed_amplitudes((k0, k1, k2), amps, x0)
            ttriad = register_closed_helical_triad(
                wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=tamps
            )
            wt = max(wt, _multiset_relative(_work_multiset(triad), _work_multiset(ttriad)))

            lam = math.exp(float(rng.uniform(-5.0, 5.0)))
            striad = register_closed_helical_triad(
                wavevectors=(lam*k0, lam*k1, lam*k2), helicities=helicities, amplitudes=amps
            )
            expected_scale = tuple(lam*x for x in _work_multiset(triad))
            ws = max(ws, _multiset_relative(expected_scale, _work_multiset(striad)))

            amp_scale = math.exp(float(rng.uniform(-4.0, 4.0)))
            atriad = register_closed_helical_triad(
                wavevectors=(k0, k1, k2),
                helicities=helicities,
                amplitudes=tuple(amp_scale*a for a in amps),
            )
            expected_amp = tuple(amp_scale**3*x for x in _work_multiset(triad))
            wa = max(wa, _multiset_relative(expected_amp, _work_multiset(atriad)))

            nk, na = global_reality_negation((k0, k1, k2), amps)
            rtriad = register_closed_helical_triad(wavevectors=nk, helicities=helicities, amplitudes=na)
            wreal = max(wreal, _multiset_relative(_work_multiset(triad), _work_multiset(rtriad)))

    anti = generic_two_donor_counterexample()
    anti_ok = anti.donor_kernel.donor_count == 2 and anti.donor_kernel.recipient_count == 1
    if not anti_ok:
        raise AssertionError("generic two-donor anti-theorem regressed")
    if one == 0 or two == 0:
        raise AssertionError("random stress failed to exercise both one- and two-donor triads")
    _good, side = signed_good_integer_triad()
    return CyclicDonorKernelStress(
        samples=int(samples),
        one_donor_cases=one,
        two_donor_cases=two,
        zero_work_cases=zero,
        worst_energy_conservation_relative=we,
        worst_cyclic_coupling_relative=wc,
        worst_donor_marginal_relative=wd,
        worst_recipient_marginal_relative=wr,
        worst_measure_donor_marginal_relative=wmd,
        worst_measure_recipient_marginal_relative=wmr,
        worst_restricted_negative_mass_residual=wrestrict,
        worst_permutation_work_residual=wp,
        worst_translation_work_residual=wt,
        worst_wavevector_scaling_residual=ws,
        worst_amplitude_cubic_scaling_residual=wa,
        worst_reality_work_multiset_residual=wreal,
        signed_good_side_ratio=side.side_to_recipient_ratio,
        signed_good_child_donor_ratio=side.recipient_to_donor_ratio,
        signed_good_side_donor_ratio=side.side_to_donor_ratio,
        generic_two_donor_counterexample_passed=anti_ok,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=75_000)
    parser.add_argument("--seed", type=int, default=2026081203)
    parser.add_argument("--outdir", type=Path, default=Path("results-cyclic-helical-triad-donor-kernel"))
    args = parser.parse_args()
    out = stress(args.samples, args.seed)
    cert = theorem_certificate()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "certificate.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2, sort_keys=True) + "\n"
    )
    summary = f"""# Cyclic helical-triad donor/recipient kernel\n\nStatus: **{STATUS}**.\n\nOne regular closed helical triad is quotient by the full `S3` ordering symmetry.  Marking one of its three physical roots recovers exactly the existing parent-swap edge quotient: `3*(1/48)=1/16`.  The three cyclic root works use the same Waleffe phase factor and satisfy `T0+T1+T2=0` before Hahn.\n\nThe positive same-triad transport `M(i->j)=[-T_i]_+[T_j]_+/Q`, `Q=sum[T]_+=sum[-T]_+`, has donor marginal `dW-` and recipient marginal the already-canonical `dW+`.  On three slots this transport is unique because every nonzero zero-sum sign pattern has a singleton donor or recipient side.  It adds same-time donor provenance; it does not replace causality or create a new event.\n\nStress: `{out.samples}` generic physical closed triads\n- one-donor cases: `{out.one_donor_cases}`\n- two-donor cases: `{out.two_donor_cases}`\n- zero-work cases: `{out.zero_work_cases}`\n- worst cyclic energy-conservation relative residual: `{out.worst_energy_conservation_relative:.3e}`\n- worst cyclic-coupling relative residual: `{out.worst_cyclic_coupling_relative:.3e}`\n- worst density donor/recipient marginal residuals: `{out.worst_donor_marginal_relative:.3e}` / `{out.worst_recipient_marginal_relative:.3e}`\n- worst measure donor/recipient marginal residuals: `{out.worst_measure_donor_marginal_relative:.3e}` / `{out.worst_measure_recipient_marginal_relative:.3e}`\n- worst restricted dW- to canonical dW+ submeasure mass residual: `{out.worst_restricted_negative_mass_residual:.3e}`\n- worst permutation / translation work residuals: `{out.worst_permutation_work_residual:.3e}` / `{out.worst_translation_work_residual:.3e}`\n- worst wavevector / amplitude scaling residuals: `{out.worst_wavevector_scaling_residual:.3e}` / `{out.worst_amplitude_cubic_scaling_residual:.3e}`\n- worst global-reality work-multiset residual: `{out.worst_reality_work_multiset_residual:.3e}`\n\nOn the signed-good integer triad, the interaction parents remain two but exactly one is the energy donor.  The other is a simultaneous positive nonforward side recipient.  The certified work ratios are\n\n`3/10 < W_side+/W_child+ < 1/3`,\n`3/4 < W_child+/W_donor- < 10/13`,\n`3/13 < W_side+/W_donor- < 1/4`.\n\nThe side work is real physical energy redistribution and already lies on the existing positive-nonforward transfer-loss route.  It is not dissipation, a reset budget, or a proof that the good child itself terminates.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
