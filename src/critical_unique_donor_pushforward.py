from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from src.continuum_helical_edge_measure_registration import HelicalModeIdentity
from src.critical_parabolic_donor_threshold import theorem_certificate as parabolic_threshold_certificate
from src.cyclic_helical_triad_donor_kernel import ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel
from src.helical import stable_norm3
from src.helical_energy_helicity_barycentric_rigidity import (
    certify_helical_energy_helicity_rigidity,
    critical_helicity_pair_balance,
)

STATUS = (
    "DRAFT_CRITICAL_EFFICIENT_UNIQUE_ENERGY_DONOR_PUSHFORWARD__"
    "CANONICAL_RECIPIENT_DW_PLUS_RESTRICTED_THROUGH_EXISTING_CYCLIC_KERNEL__"
    "MASS_PRESERVING_DONOR_MODE_SCALE_MEASURE__ALL_DONOR_RATIOS_LT_5_OVER_8__NO_ARGMAX"
)


def _radius(mode: HelicalModeIdentity) -> float:
    return stable_norm3(mode.wavevector)


def _close(a: float, b: float, scale: float) -> bool:
    return abs(float(a)-float(b)) <= 8.0e-10*max(abs(float(scale)),abs(float(a)),abs(float(b)),1.0e-300)


@dataclass(frozen=True)
class CriticalEfficientRecipientDonorPushforward:
    recipient_closed_mode_index: int
    donor_closed_mode_index: int
    recipient_mode: HelicalModeIdentity
    donor_mode: HelicalModeIdentity
    canonical_recipient_positive_mass: float
    pushed_donor_mass: float
    child_frequency: float
    donor_frequency: float
    donor_child_ratio: float
    critical_efficiency: float
    critical_efficiency_threshold: float
    mass_native_residual: float
    donor_is_unique_median_energy_donor: bool = True
    cyclic_kernel_reused: bool = True
    later_hahn_used: bool = False
    argmax_scale_selector_used: bool = False
    creates_new_event: bool = False

    def __post_init__(self) -> None:
        vals=(
            self.canonical_recipient_positive_mass,self.pushed_donor_mass,
            self.child_frequency,self.donor_frequency,self.donor_child_ratio,
            self.critical_efficiency,self.critical_efficiency_threshold,self.mass_native_residual,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite critical donor pushforward data required")
        if min(self.canonical_recipient_positive_mass,self.pushed_donor_mass,self.child_frequency,self.donor_frequency)>0.0:
            pass
        else:
            raise ValueError("positive recipient/donor mass and frequencies required")
        if not self.critical_efficiency > self.critical_efficiency_threshold:
            raise ValueError("recipient is not on the super-threshold critical-efficient law")
        if not self.donor_child_ratio < 5.0/8.0:
            raise AssertionError("critical-efficient donor pushforward escaped D<5/8")
        if self.mass_native_residual>8.0e-10:
            raise AssertionError("critical recipient->donor pushforward changed canonical work mass")
        if not self.donor_is_unique_median_energy_donor or not self.cyclic_kernel_reused:
            raise ValueError("critical donor provenance must come from the existing unique cyclic donor law")
        if self.later_hahn_used or self.argmax_scale_selector_used or self.creates_new_event:
            raise ValueError("critical donor pushforward may not re-Hahn, argmax-select a scale, or mint an event")


def critical_efficient_recipient_donor_pushforward(
    triad: ClosedHelicalTriadRegistration,
    kernel: CyclicTriadMeasureKernel,
    *,
    recipient_closed_mode_index: int,
) -> CriticalEfficientRecipientDonorPushforward:
    i=int(recipient_closed_mode_index)
    if i not in (0,1,2):
        raise ValueError("recipient index must be one closed-triad root")
    if not kernel.numerically_resolved_transport:
        raise ValueError("critical donor pushforward refuses unresolved cyclic work signs")
    rigidity=certify_helical_energy_helicity_rigidity(triad)
    critical=critical_helicity_pair_balance(triad)
    if rigidity.transfer_orientation != "mean_preserving_spread":
        raise ValueError("critical-efficient recipient donor theorem requires a physical spread")
    if i not in rigidity.strict_uv_frontier_positive_slots:
        raise ValueError("recipient must be a strict positive UV-frontier root")
    slot=triad.slot_for_closed_mode_index(i)
    M=float(slot.edge_registration.child_frequency)
    A=float(slot.edge_registration.native_modal_capacity)
    if M<=0.0 or A<=0.0:
        raise ValueError("positive native child scale/capacity required")
    efficiency=float(critical.absolute_critical_rate)/(M*A)
    threshold=parabolic_threshold_certificate().critical_efficiency_threshold
    if not efficiency>threshold:
        raise ValueError("strict UV recipient is below the critical parabolic threshold")

    donor_index=rigidity.ordered_slot_indices[1]
    donor_slot=triad.slot_for_closed_mode_index(donor_index)
    if donor_slot.signed_work>=0.0:
        raise AssertionError("median curl mode is not the unique energy donor on a spread")
    recipient_mass=float(kernel.recipient_edge_positive_masses[i])
    if recipient_mass<=0.0:
        raise ValueError("selected critical-efficient recipient carries no canonical dW+ mass")
    atoms=tuple(a for a in kernel.atoms if a.recipient_closed_mode_index==i)
    if not atoms:
        raise AssertionError("canonical recipient dW+ lost cyclic donor provenance")
    if any(a.donor_closed_mode_index!=donor_index for a in atoms):
        raise AssertionError("critical-efficient recipient has more than its unique median energy donor")
    pushed=math.fsum(a.physical_work_mass for a in atoms)
    native=max(kernel.native_work_mass_scale,recipient_mass,pushed,1.0e-300)
    if not _close(pushed,recipient_mass,native):
        raise AssertionError("recipient restriction failed to push back mass-preservingly to its unique donor")
    donor_mode=atoms[0].donor_child_mode
    recipient_mode=atoms[0].recipient_child_mode
    if any(a.donor_child_mode!=donor_mode or a.recipient_child_mode!=recipient_mode for a in atoms):
        raise AssertionError("one closed recipient restriction changed physical donor/recipient identities")
    Nd=_radius(donor_mode)
    Nc=_radius(recipient_mode)
    ratio=Nd/Nc
    if not ratio<5.0/8.0:
        raise AssertionError("actual critical-efficient unique donor violates the exact parabolic theorem")
    return CriticalEfficientRecipientDonorPushforward(
        recipient_closed_mode_index=i,
        donor_closed_mode_index=donor_index,
        recipient_mode=recipient_mode,
        donor_mode=donor_mode,
        canonical_recipient_positive_mass=recipient_mass,
        pushed_donor_mass=pushed,
        child_frequency=Nc,
        donor_frequency=Nd,
        donor_child_ratio=ratio,
        critical_efficiency=efficiency,
        critical_efficiency_threshold=threshold,
        mass_native_residual=abs(pushed-recipient_mass)/native,
    )


@dataclass(frozen=True)
class CriticalDonorScaleMeasure:
    atoms: tuple[CriticalEfficientRecipientDonorPushforward,...]
    total_recipient_mass: float
    total_donor_pushforward_mass: float
    donor_scale_weighted_mass: float
    maximum_donor_child_ratio: float
    mass_native_residual: float
    argmax_scale_selector_used: bool = False

    def __post_init__(self)->None:
        if not self.atoms:
            raise ValueError("nonempty critical donor-scale measure required")
        if self.mass_native_residual>8.0e-10:
            raise AssertionError("critical donor-scale pushforward changed total canonical mass")
        if not self.maximum_donor_child_ratio<5.0/8.0:
            raise AssertionError("critical donor-scale law escaped the uniform 5/8 ratio")
        if self.argmax_scale_selector_used:
            raise ValueError("critical donor-scale law must retain the full pushforward distribution")


def critical_donor_scale_measure(
    atoms: Sequence[CriticalEfficientRecipientDonorPushforward],
)->CriticalDonorScaleMeasure:
    rows=tuple(atoms)
    if not rows:
        raise ValueError("nonempty critical donor pushforwards required")
    rec=math.fsum(a.canonical_recipient_positive_mass for a in rows)
    don=math.fsum(a.pushed_donor_mass for a in rows)
    weighted=math.fsum(a.donor_frequency*a.pushed_donor_mass for a in rows)
    native=max(rec,don,1.0e-300)
    return CriticalDonorScaleMeasure(
        atoms=rows,
        total_recipient_mass=rec,
        total_donor_pushforward_mass=don,
        donor_scale_weighted_mass=weighted,
        maximum_donor_child_ratio=max(a.donor_child_ratio for a in rows),
        mass_native_residual=abs(rec-don)/native,
    )


def theorem_certificate()->dict[str,object]:
    return {
        "status":STATUS,
        "domain":"strict-UV positive heterochiral spreads whose native critical-production efficiency exceeds 1/(8 sqrt2)",
        "donor":"energy+helicity barycentric rigidity makes the median curl mode the unique physical energy donor",
        "measure":"restrict the already-canonical recipient dW+ and use the existing cyclic donor kernel; its preimage is the unique donor and the pushforward preserves work mass exactly",
        "scale":"every pushed donor atom has N_d/N_c<5/8 by the critical parabolic theorem",
        "distribution":"retain the complete positive donor-mode/scale pushforward; no argmax, largest cell, or output-shell selector is introduced",
        "causal":"no later Hahn, no new event, no capacity-as-causality, and no replacement of canonical dW+",
        "scope":"this solves same-event donor-scale disintegration; between-time continuation into an earlier generated donor event remains a separate state-continuity problem",
    }
