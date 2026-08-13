from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.critical_shell_service_reentry import (
    critical_shell_bounded_service_lower,
    critical_shell_integrated_service_lower,
)
from src.cyclic_helical_triad_donor_kernel import (
    ClosedHelicalTriadRegistration,
    CyclicTriadMeasureKernel,
    cyclic_triad_measure_kernel,
    register_closed_helical_triad,
)
from src.hard_tail_true_upward_supply import (
    HardTailUpwardSupplySplit,
    UpwardSupplyAtom,
    deep_upward_resolved_contact_fixture,
    hard_tail_upward_supply_split,
)
from src.high_tail_natural_window_reentry import natural_window_geometry, temporal_concentration_statistics
from src.high_tail_ultraviolet_locality import ultraviolet_hh_work_constant


STATUS = (
    "EXACT_RESOLVED_CONTACT_NATIVE_BINDING__SIGNED_VH_HH_REPARTITION_BEFORE_HAHN__"
    "CANONICAL_DWPLUS_RESTRICTS_WITHOUT_CLONING__MIXED_TO_K_OR_STRAIN_COVER__"
    "HH_COMPLEMENT_EDGE_TOTAL_VARIATION_BEFORE_AGGREGATE_HAHN__"
    "R_5_OVER_4_DIRECT_NATURAL_WINDOW__COMMON_N_DW_UNIT"
)

K_RELAY = "conservative_same_event_K_relink"
S_STRAIN = "existing_symmetric_strain_deformation"
HH_WINDOW = "resolved_contact_HH_complement_natural_window"
CONTACT_PARENT_UPPER_RATIO = 5.0 / 4.0


def _finite_positive(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"positive finite {name} required")
    return out


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


@dataclass(frozen=True)
class ResolvedContactSmoothBinding:
    boundary: float
    recipient_shell_scale: float
    recipient_shell_index: int
    donor_closed_mode_index: int
    recipient_closed_mode_index: int
    recipient_radius: float
    resolved_parent_radius: float
    uv_parent_radius: float
    donor_radius: float
    donor_is_resolved_parent: bool
    resolved_parent_cutoff_value: float
    canonical_positive_mass: float
    mixed_vh_submeasure_mass: float
    hh_complement_submeasure_mass: float
    common_unit_scale: float
    maximum_parent_to_shell_ratio: float
    canonical_cause_preserved: bool = True
    same_time_donor_provenance_preserved: bool = True
    later_hahn_used: bool = False
    low_low_mass_created: bool = False
    owner_mass_cloned: bool = False
    recipient_shell_reweighting_used: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        M = _finite_positive(self.recipient_shell_scale, "recipient shell scale")
        rr = _finite_positive(self.recipient_radius, "recipient radius")
        rl = _finite_positive(self.resolved_parent_radius, "resolved parent radius")
        rh = _finite_positive(self.uv_parent_radius, "UV parent radius")
        rd = _finite_positive(self.donor_radius, "donor radius")
        mu = _finite_positive(self.canonical_positive_mass, "canonical positive mass")
        q = float(self.resolved_parent_cutoff_value)
        if not math.isfinite(q) or not (0.0 <= q <= 1.0):
            raise ValueError("positive smooth cutoff value must lie in [0,1]")
        if self.recipient_shell_index < 1:
            raise ValueError("resolved-contact recipient shell must lie above the tail boundary")
        if abs(M - (2.0 ** self.recipient_shell_index) * N) > 5.0e-12 * M:
            raise AssertionError("recipient shell scale changed from boundary-anchored dyadic geometry")
        if not (0.5 * M < rr <= M):
            raise AssertionError("recipient mode left its hard output shell")
        tol = 5.0e-12 * M
        if rl > 0.25 * M + tol:
            raise AssertionError("resolved-contact parent left B_(M/4)")
        if abs(rl - 0.25 * M) <= tol and q > 5.0e-12:
            raise ValueError("a smooth multiplier supported in B_(M/4) must vanish at boundary contact")
        if not rh > 0.25 * M - tol:
            raise AssertionError("resolved-contact edge must have exactly one parent in B_(M/4)")
        if rh > rr + rl + tol or rh > CONTACT_PARENT_UPPER_RATIO * M + tol:
            raise AssertionError("resolved-contact UV parent lost exact <=5M/4 triangle bound")
        if rd > N + 5.0e-12 * N:
            raise AssertionError("upward donor left the low side of the physical tail boundary")
        if abs(self.maximum_parent_to_shell_ratio - max(rl, rh) / M) > 5.0e-12:
            raise AssertionError("stored parent/shell ratio changed from physical geometry")
        mixed = _finite_nonnegative(self.mixed_vh_submeasure_mass, "mixed V-h submeasure mass")
        hh = _finite_nonnegative(self.hh_complement_submeasure_mass, "HH-complement submeasure mass")
        scale = max(1.0, mu)
        if abs(mixed - q * mu) > 2.0e-12 * scale:
            raise AssertionError("mixed submeasure is not q times the canonical positive cause")
        if abs(hh - (1.0 - q) * mu) > 2.0e-12 * scale:
            raise AssertionError("HH complement is not (1-q) times the canonical positive cause")
        if abs(mixed + hh - mu) > 2.0e-12 * scale:
            raise AssertionError("smooth repartition cloned or lost canonical dW+")
        if self.common_unit_scale != N:
            raise ValueError("resolved-contact binding must retain the parent-tail N dW unit")
        if (
            not self.canonical_cause_preserved
            or not self.same_time_donor_provenance_preserved
            or self.later_hahn_used
            or self.low_low_mass_created
            or self.owner_mass_cloned
            or self.recipient_shell_reweighting_used
        ):
            raise ValueError("resolved-contact binding changed certified causal semantics")

    @property
    def mixed_common_unit_mass(self) -> float:
        return self.common_unit_scale * self.mixed_vh_submeasure_mass

    @property
    def hh_common_unit_mass(self) -> float:
        return self.common_unit_scale * self.hh_complement_submeasure_mass


def resolved_contact_smooth_binding(
    atom: UpwardSupplyAtom,
    *,
    resolved_parent_cutoff_value: float,
) -> ResolvedContactSmoothBinding:
    """Restrict one canonical upward dW+ atom through the actual u=V+h split.

    On a resolved-contact recipient edge, the shell lower bound |z|>M/2 makes
    two low parents impossible.  Hence exactly one interaction parent has
    radius <=M/4 and the other is >M/4.  For any real nonnegative smooth
    low-pass S supported in B_(M/4), the UV parent has S=0.  Writing
    q=S(k_low), bilinearity gives the signed edge identity

        dW = q dW|_{V-h mixed} + (1-q) dW|_{h-h}.

    The incoming atom is already a positive donor-restricted submeasure of the
    canonical dW+.  Multiplication by q and 1-q is therefore a positive
    restriction/pushforward of that same cause, not another Hahn split.
    """
    if not atom.resolved_scale_parent_contact or atom.pure_uv_hh_by_support:
        raise ValueError("resolved-contact canonical upward atom required")
    q = float(resolved_parent_cutoff_value)
    if not math.isfinite(q) or not (0.0 <= q <= 1.0):
        raise ValueError("resolved-parent cutoff value must lie in [0,1]")
    M = atom.recipient_shell_scale
    radii = tuple(float(r) for r in atom.interaction_parent_radii)
    if len(radii) != 2:
        raise ValueError("exactly two interaction parents required")
    low_i = 0 if radii[0] <= radii[1] else 1
    low = radii[low_i]
    high = radii[1 - low_i]
    tol = 5.0e-12 * M
    if low > 0.25 * M + tol:
        raise AssertionError("resolved-contact atom lost its low parent")
    if not high > 0.25 * M - tol:
        raise AssertionError("two resolved parents cannot feed a shell with |z|>M/2")
    donor_is_low = abs(atom.donor_radius - low) <= tol
    mu = atom.physical_work_mass
    return ResolvedContactSmoothBinding(
        boundary=atom.boundary,
        recipient_shell_scale=M,
        recipient_shell_index=atom.recipient_shell_index,
        donor_closed_mode_index=atom.donor_closed_mode_index,
        recipient_closed_mode_index=atom.recipient_closed_mode_index,
        recipient_radius=atom.recipient_radius,
        resolved_parent_radius=low,
        uv_parent_radius=high,
        donor_radius=atom.donor_radius,
        donor_is_resolved_parent=donor_is_low,
        resolved_parent_cutoff_value=q,
        canonical_positive_mass=mu,
        mixed_vh_submeasure_mass=q * mu,
        hh_complement_submeasure_mass=(1.0 - q) * mu,
        common_unit_scale=atom.common_unit_scale,
        maximum_parent_to_shell_ratio=max(radii) / M,
    )


@dataclass(frozen=True)
class SingleChargedRecipientMixedCause:
    recipient_closed_mode_index: int
    donor_closed_mode_indices: tuple[int, ...]
    canonical_mixed_submeasure_mass: float
    common_unit_scale: float
    recipient_submeasure_single_charged: bool = True
    donor_provenance_is_sidecar: bool = True

    def __post_init__(self) -> None:
        if self.recipient_closed_mode_index not in (0, 1, 2):
            raise ValueError("recipient cyclic root must lie in {0,1,2}")
        if not self.donor_closed_mode_indices or len(set(self.donor_closed_mode_indices)) != len(self.donor_closed_mode_indices):
            raise ValueError("nonempty unique donor sidecar indices required")
        if any(i not in (0, 1, 2) or i == self.recipient_closed_mode_index for i in self.donor_closed_mode_indices):
            raise ValueError("donor sidecars must be distinct cyclic roots from the recipient")
        _finite_nonnegative(self.canonical_mixed_submeasure_mass, "coalesced canonical mixed mass")
        _finite_positive(self.common_unit_scale, "coalesced common N scale")
        if not self.recipient_submeasure_single_charged or not self.donor_provenance_is_sidecar:
            raise ValueError("coalesced recipient cause must be single-charged with donor provenance only as sidecar")


def coalesce_recipient_mixed_cause(
    bindings: Sequence[ResolvedContactSmoothBinding],
) -> SingleChargedRecipientMixedCause:
    """Coalesce all donor provenance for one recipient before K/S ownership.

    The cyclic donor kernel may split one canonical positive recipient edge among
    two same-time donors.  Those donor labels survive as provenance, but the
    recipient physical mixed work has only one K/S decomposition.  Summing the
    donor-restricted mixed masses first prevents repeated owner charging.
    """
    rows = tuple(bindings)
    if not rows:
        raise ValueError("nonempty recipient mixed binding family required")
    recipient = rows[0].recipient_closed_mode_index
    N = rows[0].common_unit_scale
    q = rows[0].resolved_parent_cutoff_value
    shell = rows[0].recipient_shell_scale
    if any(row.recipient_closed_mode_index != recipient for row in rows):
        raise ValueError("K/S coalescence cannot mix distinct recipient physical edges")
    if any(abs(row.common_unit_scale - N) > 5e-12 * N for row in rows):
        raise ValueError("recipient donor sidecars changed the common N dW unit")
    if any(abs(row.recipient_shell_scale - shell) > 5e-12 * shell for row in rows):
        raise ValueError("recipient donor sidecars changed the physical recipient shell")
    if any(abs(row.resolved_parent_cutoff_value - q) > 5e-12 for row in rows):
        raise ValueError("same recipient edge acquired inconsistent smooth cutoff values")
    donors = tuple(sorted(row.donor_closed_mode_index for row in rows))
    if len(set(donors)) != len(donors):
        raise ValueError("duplicate donor sidecar would double-charge canonical recipient cause")
    mass = math.fsum(row.mixed_vh_submeasure_mass for row in rows)
    return SingleChargedRecipientMixedCause(
        recipient_closed_mode_index=recipient,
        donor_closed_mode_indices=donors,
        canonical_mixed_submeasure_mass=mass,
        common_unit_scale=N,
    )


@dataclass(frozen=True)
class SignedResolvedKSAtom:
    signed_mixed_work: float
    signed_skew_work: float
    signed_strain_work: float
    same_physical_atom: bool = True
    same_resolved_operator: bool = True
    observer_gauge_quotiented_or_fixed_event: bool = True

    def __post_init__(self) -> None:
        I, K, S = (float(self.signed_mixed_work), float(self.signed_skew_work), float(self.signed_strain_work))
        if not all(math.isfinite(v) for v in (I, K, S)):
            raise ValueError("finite signed resolved K/S work required")
        scale = max(1.0, abs(I), abs(K), abs(S))
        if abs(I - K - S) > 8.0e-12 * scale:
            raise ValueError("signed I=K+S identity must hold before positive ownership")
        if not (self.same_physical_atom and self.same_resolved_operator and self.observer_gauge_quotiented_or_fixed_event):
            raise ValueError("K/S cover refuses representation substitution or unquotiented observer motion")


@dataclass(frozen=True)
class PositiveKSOwnerCover:
    canonical_mixed_submeasure_mass: float
    positive_mixed_work: float
    positive_skew_work: float
    positive_strain_work: float
    common_unit_scale: float
    signed_identity_residual: float
    positive_cover_margin: float
    half_owner_threshold: float
    joint_owner_witnesses: tuple[str, ...]
    canonical_cause_unsplit: bool = True
    recipient_submeasure_single_charged: bool = True
    donor_provenance_is_sidecar: bool = True
    later_hahn_on_canonical_cause: bool = False
    owner_fraction_matching_used: bool = False
    own_shell_reweighting_used: bool = False

    def __post_init__(self) -> None:
        mu = _finite_nonnegative(self.canonical_mixed_submeasure_mass, "canonical mixed submeasure mass")
        Ip = _finite_nonnegative(self.positive_mixed_work, "positive mixed work")
        Kp = _finite_nonnegative(self.positive_skew_work, "positive skew work")
        Sp = _finite_nonnegative(self.positive_strain_work, "positive strain work")
        N = _finite_positive(self.common_unit_scale, "common N scale")
        del N
        tol = 2.0e-11 * max(1.0, mu, Ip, Kp, Sp)
        if mu > Ip + tol:
            raise AssertionError("canonical mixed cause exceeds same-atom positive mixed work")
        if Kp + Sp + tol < Ip or self.positive_cover_margin < -tol:
            raise AssertionError("positive mixed work escaped K+/S+ cover")
        if self.signed_identity_residual > 2.0e-11:
            raise AssertionError("K/S owner cover left the signed physical identity")
        if mu > 0.0 and not self.joint_owner_witnesses:
            raise AssertionError("positive mixed canonical cause has no K/S owner witness")
        if (
            not self.canonical_cause_unsplit
            or not self.recipient_submeasure_single_charged
            or not self.donor_provenance_is_sidecar
            or self.later_hahn_on_canonical_cause
            or self.owner_fraction_matching_used
            or self.own_shell_reweighting_used
        ):
            raise ValueError("K/S cover requires one single-charged recipient submeasure; donor provenance cannot clone owner mass")


def cover_canonical_mixed_submeasure_by_ks(
    cause: SingleChargedRecipientMixedCause,
    signed_atom: SignedResolvedKSAtom,
) -> PositiveKSOwnerCover:
    """Keep the canonical cause unsplit and cover it by existing K+/S+ laws.

    From the same-atom signed identity I=K+S one has I+<=K++S+.  Since the
    coalesced donor-restricted mixed cause mu is already dominated by I+, it follows that
    mu<=K++S+.  No proportional allocation is manufactured.  The half-owner
    statement is only a witness alternative: at least one existing component
    positive law has mass >=mu/2.
    """
    if not isinstance(cause, SingleChargedRecipientMixedCause):
        raise TypeError("K/S ownership requires a coalesced single-charged recipient cause")
    mu = _finite_nonnegative(cause.canonical_mixed_submeasure_mass, "canonical mixed submeasure mass")
    N = _finite_positive(cause.common_unit_scale, "common N scale")
    I = float(signed_atom.signed_mixed_work)
    K = float(signed_atom.signed_skew_work)
    S = float(signed_atom.signed_strain_work)
    Ip, Kp, Sp = max(I, 0.0), max(K, 0.0), max(S, 0.0)
    scale = max(1.0, abs(I), abs(K), abs(S), mu)
    if mu > Ip + 8.0e-12 * scale:
        raise ValueError("canonical mixed submeasure is not dominated by same-atom mixed dW+")
    cover = Kp + Sp - Ip
    if cover < -8.0e-12 * scale:
        raise AssertionError("signed K/S identity failed its positive cover")
    threshold = 0.5 * mu
    witnesses: list[str] = []
    tol = 8.0e-12 * scale
    if mu > 0.0:
        if Kp + tol >= threshold:
            witnesses.append(K_RELAY)
        if Sp + tol >= threshold:
            witnesses.append(S_STRAIN)
        if not witnesses:
            raise AssertionError("K+/S+ positive cover lost the half-owner alternative")
    return PositiveKSOwnerCover(
        canonical_mixed_submeasure_mass=mu,
        positive_mixed_work=Ip,
        positive_skew_work=Kp,
        positive_strain_work=Sp,
        common_unit_scale=N,
        signed_identity_residual=abs(I - K - S) / scale,
        positive_cover_margin=cover,
        half_owner_threshold=threshold,
        joint_owner_witnesses=tuple(witnesses),
        recipient_submeasure_single_charged=True,
        donor_provenance_is_sidecar=True,
    )


def resolved_contact_component_route(
    *,
    required_contact_common_work_lower: float,
    actual_contact_common_work: float,
    hh_complement_common_work: float,
    mixed_common_work: float,
    positive_skew_common_work: float,
    positive_strain_common_work: float,
) -> dict[str, object]:
    """Integrated owner alternative in the unchanged N dW unit.

    Contact = HH_complement + mixed exactly.  If contact>=L, then either
    HH_complement>=L/2 or mixed>=L/2.  On the latter branch the signed K/S
    cover gives K+ or S+ at least L/4.  K is same-event conservative provenance;
    S is the existing strain/deformation owner.
    """
    L = _finite_positive(required_contact_common_work_lower, "required contact lower")
    C = _finite_positive(actual_contact_common_work, "actual contact common work")
    H = _finite_nonnegative(hh_complement_common_work, "HH complement common work")
    I = _finite_nonnegative(mixed_common_work, "mixed common work")
    K = _finite_nonnegative(positive_skew_common_work, "positive skew common work")
    S = _finite_nonnegative(positive_strain_common_work, "positive strain common work")
    tol = 2.0e-11 * max(1.0, L, C, H, I, K, S)
    if C + tol < L:
        raise ValueError("actual resolved-contact work is below the supplied owner lower")
    if abs(H + I - C) > tol:
        raise ValueError("resolved-contact common work did not split exactly into HH plus mixed")
    if K + S + tol < I:
        raise ValueError("integrated mixed common work escaped the existing K+/S+ cover")
    owners: list[str] = []
    if H + tol >= 0.5 * L:
        owners.append(HH_WINDOW)
    if K + tol >= 0.25 * L:
        owners.append(K_RELAY)
    if S + tol >= 0.25 * L:
        owners.append(S_STRAIN)
    if not owners:
        raise AssertionError("resolved-contact owner cover lost all physical continuations")
    return {
        "required_contact_common_work_lower": L,
        "actual_contact_common_work": C,
        "hh_complement_common_work": H,
        "mixed_common_work": I,
        "positive_skew_common_work": K,
        "positive_strain_common_work": S,
        "hh_owner_threshold": 0.5 * L,
        "ks_owner_threshold": 0.25 * L,
        "joint_physical_continuations": tuple(owners),
        "canonical_cause_replaced": False,
        "later_hahn_used": False,
        "recipient_shell_reweighting_used": False,
    }




def hard_tail_resolved_contact_route(
    *,
    physical_tail_dissipation: float,
    viscosity: float,
    actual_contact_common_work: float,
    hh_complement_common_work: float,
    mixed_common_work: float,
    positive_skew_common_work: float,
    positive_strain_common_work: float,
) -> dict[str, object]:
    """Attach the binding to the certified true-upward contact-owner lower.

    Once the upstream true-upward support alternative selects resolved contact,
    its common-unit mass is at least nu D_tail/2.  Therefore the native binding
    gives exactly one or more of

      contact HH complement >= nu D_tail/4,
      K+ relink witness       >= nu D_tail/8,
      S+ strain witness       >= nu D_tail/8.
    """
    D = _finite_positive(physical_tail_dissipation, "physical tail dissipation")
    nu = _finite_positive(viscosity, "viscosity")
    clean_contact = 0.5 * nu * D
    out = resolved_contact_component_route(
        required_contact_common_work_lower=clean_contact,
        actual_contact_common_work=actual_contact_common_work,
        hh_complement_common_work=hh_complement_common_work,
        mixed_common_work=mixed_common_work,
        positive_skew_common_work=positive_skew_common_work,
        positive_strain_common_work=positive_strain_common_work,
    )
    return {
        **out,
        "clean_contact_owner_lower": clean_contact,
        "clean_contact_HH_lower": 0.25 * nu * D,
        "clean_K_or_S_lower": 0.125 * nu * D,
        "physical_tail_dissipation": D,
        "viscosity": nu,
    }


def canonical_contact_hh_shell_law(positive_shell_common_work: Mapping[int, float]) -> dict[str, object]:
    """Positive shell pushforward of the same donor-restricted canonical HH submeasure.

    No signed shell aggregation and no downstream Hahn operation occurs here.
    Each input mass is already a positive restriction of canonical edge dW+.
    """
    if not positive_shell_common_work:
        raise ValueError("nonempty canonical contact-HH shell law required")
    items = sorted((int(j), float(w)) for j, w in positive_shell_common_work.items())
    if any(j < 1 or not math.isfinite(w) or w < 0.0 for j, w in items):
        raise ValueError("high-tail shell levels j>=1 with finite nonnegative canonical masses required")
    total = math.fsum(w for _, w in items)
    if total <= 0.0:
        raise ValueError("positive canonical contact-HH common work required")
    positive = [(j, w) for j, w in items if w > 0.0]
    j_star, w_star = max(positive, key=lambda row: row[1])
    p_star = w_star / total
    return {
        "total_canonical_contact_HH_common_work": total,
        "selected_shell_level": j_star,
        "selected_shell_common_work": w_star,
        "p_max": p_star,
        "H_inf_output_scale": -math.log(p_star),
        "later_hahn_used": False,
    }


@dataclass(frozen=True)
class CanonicalHHTotalVariationYoungBridge:
    unitary_fourier_factor: float
    helical_l1_over_vector_l2_factor: float
    unordered_parent_orbit_factor: float
    native_edge_capacity_factor: float
    exact_edge_variation_prefactor_over_A3: float
    clean_pair_young_prefactor_over_A3: float
    edge_variation_to_clean_young_ratio: float
    canonical_positive_submeasure_dominated_by_edge_variation: bool = True
    edge_variation_dominated_by_capacity_measure: bool = True
    aggregate_hahn_used: bool = False

    def __post_init__(self) -> None:
        C_F = _finite_positive(self.unitary_fourier_factor, "unitary Fourier convolution factor")
        h = _finite_positive(self.helical_l1_over_vector_l2_factor, "helical l1/l2 factor")
        orbit = _finite_positive(self.unordered_parent_orbit_factor, "unordered-parent orbit factor")
        native = _finite_positive(self.native_edge_capacity_factor, "native edge-capacity factor")
        exact = _finite_positive(self.exact_edge_variation_prefactor_over_A3, "edge-variation Young prefactor")
        clean = _finite_positive(self.clean_pair_young_prefactor_over_A3, "clean Young prefactor")
        ratio = _finite_positive(self.edge_variation_to_clean_young_ratio, "edge/clean Young ratio")
        expected_exact = native * orbit * (h**3) * C_F
        scale = max(1.0, exact, expected_exact)
        if abs(exact - expected_exact) > 5.0e-14 * scale:
            raise AssertionError("helicity-resolved edge total variation lost its quotient/Young prefactor")
        if abs(clean - native) > 5.0e-14 * max(1.0, clean, native):
            raise AssertionError("clean pair Young prefactor changed from 4 A_3")
        if abs(ratio - exact / clean) > 5.0e-14 * max(1.0, ratio):
            raise AssertionError("stored edge/clean Young ratio changed")
        if not ratio < 1.0:
            raise AssertionError("clean Young constant no longer dominates helicity edge total variation")
        if (
            not self.canonical_positive_submeasure_dominated_by_edge_variation
            or not self.edge_variation_dominated_by_capacity_measure
            or self.aggregate_hahn_used
        ):
            raise ValueError("canonical HH capacity bridge must pass through signed edge total variation before any aggregate Hahn")


def canonical_hh_edge_total_variation_young_bridge() -> CanonicalHHTotalVariationYoungBridge:
    """Bind canonical edge ``dW+`` to the clean sharp-Young window constant.

    This is the missing representation step between an atomwise canonical
    positive law and the older aggregate-HH natural-window estimate.

    For the continuum helicity-resolved edge measure, native capacity is

        A_e = 4 |z| |a_x a_y a_z|.

    Passing from ordered parents to the unordered physical orbit contributes
    ``1/2``.  At each frequency the two helical coefficients obey

        sum_s |a_s(k)| <= sqrt(2) |u_hat(k)|.

    Therefore Young on the *total variation* of the edge law has coefficient

        4 * (1/2) * (sqrt(2))^3 * C_F * A_3
        = 4 sqrt(2) C_F A_3.

    The clean physical pair inequality already uses ``4 A_3``.  Hence the
    exact helicity-edge variation / clean-Young ratio is

        sqrt(2) C_F < 1,

    with ``C_F=(2pi)^(-3/2)``.  Consequently every donor-restricted canonical
    positive submeasure satisfies

        dmu_HH <= d|W_HH| <= dA_HH <= clean Young envelope,

    before any signed shell aggregation or downstream Hahn operation.
    """
    C_F = unitary_fourier_convolution_factor()
    helical = math.sqrt(2.0)
    orbit = 0.5
    native = 4.0
    exact = native * orbit * (helical**3) * C_F
    clean = native
    return CanonicalHHTotalVariationYoungBridge(
        unitary_fourier_factor=C_F,
        helical_l1_over_vector_l2_factor=helical,
        unordered_parent_orbit_factor=orbit,
        native_edge_capacity_factor=native,
        exact_edge_variation_prefactor_over_A3=exact,
        clean_pair_young_prefactor_over_A3=clean,
        edge_variation_to_clean_young_ratio=exact / clean,
    )


def canonical_contact_hh_natural_window_capacity_upper(
    window_peak_child_mass: float,
    parent_frequency: float,
    global_energy: float,
    scaled_lifetime: float,
    locality_radius: float = CONTACT_PARENT_UPPER_RATIO,
) -> float:
    """Absolute-work capacity for the canonical positive contact-HH submeasure.

    The donor-restricted canonical HH law is first dominated by the total
    variation of the same signed helicity-resolved edge measure.  The certified
    quotient/helicity bridge proves that its exact Young prefactor
    ``4 sqrt(2) C_F A_3`` is strictly below the clean ``4 A_3`` constant used
    here.  Thus no aggregate signed shell work or aggregate Hahn law is inserted.
    The chosen resolved multiplier satisfies 0<=S<=1, so the complement
    multiplier also satisfies 0<=1-S<=1 and ||h||_2<=||u||_2.  Relative to the
    older |S|<=1 capacity this removes the factor 2 on each of two HH parents.

      N W_window <= 3 sqrt(pi) R c N E_global sqrt(mu_window).
    """
    mu = _finite_nonnegative(window_peak_child_mass, "window peak child mass")
    N = _finite_positive(parent_frequency, "parent frequency")
    E = _finite_positive(global_energy, "global energy")
    c = _finite_positive(scaled_lifetime, "scaled lifetime")
    R = _finite_positive(locality_radius, "parent upper-comparability ratio")
    if R <= 1.0:
        raise ValueError("contact-HH upper-comparability ratio must exceed one")
    bridge = canonical_hh_edge_total_variation_young_bridge()
    if not bridge.edge_variation_to_clean_young_ratio < 1.0:
        raise AssertionError("clean natural-window capacity lost edge-total-variation domination")
    return R * ultraviolet_hh_work_constant() * c * N * E * math.sqrt(mu)


def contact_hh_direct_natural_window_reentry(
    *,
    positive_hh_shell_common_work: Mapping[int, float],
    required_total_hh_common_work_lower: float,
    parent_frequency: float,
    global_energy: float,
    scaled_lifetime: float,
    viscosity: float,
    maximum_window_common_work: float,
    window_length: float,
    window_peak_child_mass: float,
) -> dict[str, object]:
    """Send the actual contact-HH remainder directly to the natural-window law.

    Every resolved-contact HH remainder has both physical parent frequencies
    <=(5/4)M.  No ultraviolet-locality theorem is needed to discover that upper
    comparability.  The input shell masses are the same positive HH submeasure,
    already in the common N dW unit.
    """
    lower = _finite_positive(required_total_hh_common_work_lower, "required HH common-work lower")
    N = _finite_positive(parent_frequency, "parent frequency")
    E = _finite_positive(global_energy, "global energy")
    c = _finite_positive(scaled_lifetime, "scaled lifetime")
    nu = _finite_positive(viscosity, "viscosity")
    Ww = _finite_positive(maximum_window_common_work, "maximum window common work")
    Tw = _finite_positive(window_length, "window length")
    mu = _finite_positive(window_peak_child_mass, "window peak child mass")
    shell = canonical_contact_hh_shell_law(positive_hh_shell_common_work)
    H = float(shell["total_canonical_contact_HH_common_work"])
    if H + 6.0e-13 * max(H, lower) < lower:
        raise ValueError("actual contact-HH submeasure is below its supplied owner lower")
    j = int(shell["selected_shell_level"])
    W = float(shell["selected_shell_common_work"])
    p_s = float(shell["p_max"])
    geom = natural_window_geometry(N, j, c)
    Tnatural = float(geom["selected_natural_window"])
    if abs(Tw - Tnatural) > 8.0e-13 * max(Tw, Tnatural):
        raise ValueError("contact-HH temporal window is not the selected shell natural window")
    temporal = temporal_concentration_statistics(W, Ww, Tw)
    p_t = float(temporal["p_time"])
    capacity = canonical_contact_hh_natural_window_capacity_upper(
        mu, N, E, c, CONTACT_PARENT_UPPER_RATIO
    )
    if Ww > capacity + 8.0e-13 * max(Ww, capacity):
        raise ValueError("contact-HH natural-window work exceeds the physical comparable-parent capacity")
    weighted = math.sqrt(mu) / (p_s * p_t)
    clean = lower / (
        ultraviolet_hh_work_constant() * CONTACT_PARENT_UPPER_RATIO * c * N * E
    )
    margin = weighted - clean
    if margin < -1.0e-11 * max(weighted, clean):
        raise AssertionError("direct contact-HH scale-time concentration bound failed")
    return {
        "selected_shell_level": j,
        "selected_shell_frequency": float(geom["selected_shell_frequency"]),
        "forward_scale_ratio": float(geom["forward_scale_ratio"]),
        "natural_time_ratio": float(geom["natural_time_ratio"]),
        "selected_natural_window": Tnatural,
        "parent_upper_comparability_ratio": CONTACT_PARENT_UPPER_RATIO,
        "p_scale": p_s,
        "H_inf_output_scale": float(shell["H_inf_output_scale"]),
        "p_time": p_t,
        "H_inf_time": float(temporal["H_inf_time"]),
        "total_contact_HH_common_work": H,
        "selected_shell_common_work": W,
        "maximum_window_common_work": Ww,
        "window_peak_child_mass": mu,
        "natural_window_common_work_capacity": capacity,
        "edge_variation_to_clean_young_ratio": canonical_hh_edge_total_variation_young_bridge().edge_variation_to_clean_young_ratio,
        "weighted_sqrt_child_mass": weighted,
        "clean_weighted_sqrt_child_mass_lower": clean,
        "scale_time_tradeoff_margin": margin,
        "next_owner": "generic_critical_shell_first_stop",
        "full_survivor_own_scale_service_lower": critical_shell_bounded_service_lower(mu, c, nu),
        "full_survivor_integrated_service_lower": critical_shell_integrated_service_lower(mu, c, nu),
        "full_survivor_service_is_conditional": True,
        "output_scale_locality_theorem_used": False,
        "later_hahn_used": False,
        "packet_persistence_used": False,
        "time_partition_used": False,
        "status": STATUS,
    }


def coarse_hahn_cancellation_counterexample() -> dict[str, float | bool]:
    first, second = 1.0, -0.9
    coarse = first + second
    return {
        "canonical_positive_first_atom": first,
        "second_signed_atom": second,
        "coarse_signed_work": coarse,
        "coarse_positive_hahn_mass": max(coarse, 0.0),
        "atomic_cause_exceeds_coarse_hahn": first > max(coarse, 0.0),
    }


def boundary_contact_counterexample() -> dict[str, float | bool]:
    """Actual closed triad with resolved contact but zero mixed mass at support edge."""
    _, _, split = deep_upward_resolved_contact_fixture()
    boundary_atoms = [
        atom
        for atom in split.atoms
        if atom.resolved_scale_parent_contact
        and abs(min(atom.interaction_parent_radii) - 0.25 * atom.recipient_shell_scale)
        <= 5.0e-12 * atom.recipient_shell_scale
    ]
    if not boundary_atoms:
        raise AssertionError("certified boundary-contact fixture lost its M/4 atom")
    atom = max(boundary_atoms, key=lambda a: a.physical_work_mass)
    out = resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=0.0)
    return {
        "recipient_shell_index": atom.recipient_shell_index,
        "resolved_parent_to_shell_ratio": out.resolved_parent_radius / out.recipient_shell_scale,
        "mixed_fraction": out.mixed_vh_submeasure_mass / out.canonical_positive_mass,
        "hh_fraction": out.hh_complement_submeasure_mass / out.canonical_positive_mass,
        "contact_is_not_interface_owner": out.mixed_vh_submeasure_mass == 0.0,
    }


def _find_interior_contact_fixture() -> tuple[
    ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel, HardTailUpwardSupplySplit
]:
    """Finite physical search for one M=4N interior-contact upward edge."""
    N = math.sqrt(2.0)
    wavevectors = (
        np.asarray((-3.0, -2.0, 0.0)),
        np.asarray((1.0, 0.0, 0.0)),
        np.asarray((2.0, 2.0, 0.0)),
    )
    phases = (
        1.0 + 0.0j,
        -1.0 + 0.0j,
        0.0 + 1.0j,
        0.0 - 1.0j,
        complex(2.0**-0.5, 2.0**-0.5),
        complex(2.0**-0.5, -2.0**-0.5),
    )
    for helicities in ((a, b, c) for a in (-1, 1) for b in (-1, 1) for c in (-1, 1)):
        for phase in phases:
            for position in range(3):
                amps = [1.0 + 0.0j] * 3
                amps[position] = phase
                triad = register_closed_helical_triad(
                    wavevectors=wavevectors, helicities=helicities, amplitudes=tuple(amps)
                )
                kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
                if not kernel.numerically_resolved_transport:
                    continue
                try:
                    split = hard_tail_upward_supply_split(triad, kernel, boundary=N)
                except ValueError:
                    continue
                for atom in split.atoms:
                    M = atom.recipient_shell_scale
                    low = min(atom.interaction_parent_radii)
                    if atom.resolved_scale_parent_contact and atom.recipient_shell_index == 2 and 0.125 * M < low < 0.25 * M:
                        return triad, kernel, split
    raise AssertionError("finite physical search found no interior resolved-contact upward fixture")


def interior_contact_fixture() -> tuple[
    ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel, HardTailUpwardSupplySplit
]:
    return _find_interior_contact_fixture()


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "signed_first": "on each resolved-contact edge, dW=dW_mixed+dW_HH with dW_mixed=q dW and dW_HH=(1-q)dW before any new Hahn operation",
        "canonical_positive_push": "the donor-restricted canonical dW+ atom is restricted by q and 1-q; its same-time donor provenance is retained on both labels and total mass is unchanged",
        "contact_geometry": "for every resolved-contact upward shell j>=1, |z|>M/2 and one parent<=M/4 force the other parent into (M/4,5M/4]; low-low work cannot feed the recipient shell",
        "hard_tail_corollary": "if the certified contact owner carries at least nu D_tail/2, then contact HH carries nu D_tail/4 or existing K+ relink / S+ strain carries nu D_tail/8",
        "cutoff_uniformity": "the theorem uses only 0<=q<=1 for the resolved parent and S=0 on the >M/4 parent; no plateau radius or transition profile is part of the causal law",
        "mixed_owner": "after donor provenance is coalesced to one single-charged recipient submeasure, the same resolved operator is split signed-first as I=K+S; canonical mixed dW+ is covered once by existing K+ and strain S+ without fraction matching",
        "donor_sidecar": "multiple same-time donors retain provenance but cannot cause repeated K/S owner charging on the same recipient physical work",
        "K_semantics": "K is conservative same-event relink/circulation provenance and creates zero recursive depth",
        "S_semantics": "S is the existing symmetric strain/deformation owner, not a new interface currency",
        "HH_semantics": "the h-h complement remains a distinct canonical positive submeasure; it is dominated by signed edge total variation, then by native dA, whose exact helicity/unordered-parent Young prefactor is below the clean 4 A_3 envelope; both parents are <=5M/4, giving the direct R=5/4 natural-window capacity without aggregate Hahn",
        "edge_total_variation_bridge": asdict(canonical_hh_edge_total_variation_young_bridge()),
        "positive_cutoff_gain": "0<=S<=1 is required for positive cause restriction and simultaneously gives ||h||_2<=||u||_2, removing the older factor four from two |S|<=1 complement bounds",
        "anti_theorem": "resolved-frequency contact alone does not imply interface/mixed ownership; an M/4 boundary atom has q=0 and can be 100% h-h",
        "causal_unit": "all measures remain in the parent-tail common unit N dW; recipient shell M is geometry only",
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class ResolvedContactNativeBindingStress:
    status: str
    samples: int
    maximum_partition_residual: float
    maximum_parent_upper_ratio: float
    minimum_positive_ks_cover_margin: float
    minimum_half_owner_margin: float
    boundary_contact_hh_fraction: float
    coarse_hahn_counterexample_gap: float
    edge_variation_to_clean_young_ratio: float


def stress(samples: int = 50_000, seed: int = 2026081302) -> ResolvedContactNativeBindingStress:
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = np.random.default_rng(seed)
    max_partition = 0.0
    max_ratio = 0.0
    min_cover = math.inf
    min_half = math.inf
    for _ in range(count):
        mu = float(np.exp(rng.uniform(-5.0, 5.0)))
        q = float(rng.uniform(0.0, 1.0))
        mixed = q * mu
        hh = (1.0 - q) * mu
        max_partition = max(max_partition, abs(mixed + hh - mu) / max(mu, 1.0e-300))
        low_ratio = float(rng.uniform(1.0e-6, 0.25))
        recipient_ratio = float(rng.uniform(0.500001, 1.0))
        high_ratio = float(rng.uniform(max(0.250001, recipient_ratio - low_ratio), min(1.25, recipient_ratio + low_ratio)))
        max_ratio = max(max_ratio, high_ratio)

        I = float(np.exp(rng.uniform(-5.0, 5.0)))
        K = float(rng.uniform(-2.0 * I, 2.0 * I))
        S = I - K
        cause_mass = float(rng.uniform(0.0, 1.0)) * I
        cause = SingleChargedRecipientMixedCause(
            recipient_closed_mode_index=0,
            donor_closed_mode_indices=(1,),
            canonical_mixed_submeasure_mass=cause_mass,
            common_unit_scale=1.0,
        )
        cover = cover_canonical_mixed_submeasure_by_ks(cause, SignedResolvedKSAtom(I, K, S))
        min_cover = min(min_cover, cover.positive_cover_margin)
        half_margin = max(cover.positive_skew_work, cover.positive_strain_work) - 0.5 * cause_mass
        min_half = min(min_half, half_margin)

    boundary = boundary_contact_counterexample()
    coarse = coarse_hahn_cancellation_counterexample()
    return ResolvedContactNativeBindingStress(
        status=STATUS,
        samples=count,
        maximum_partition_residual=max_partition,
        maximum_parent_upper_ratio=max_ratio,
        minimum_positive_ks_cover_margin=min_cover,
        minimum_half_owner_margin=min_half,
        boundary_contact_hh_fraction=float(boundary["hh_fraction"]),
        coarse_hahn_counterexample_gap=float(coarse["canonical_positive_first_atom"] - coarse["coarse_positive_hahn_mass"]),
        edge_variation_to_clean_young_ratio=canonical_hh_edge_total_variation_young_bridge().edge_variation_to_clean_young_ratio,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081302)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-contact-native-binding"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples, args.seed)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "resolved_contact_native_binding.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    md = f"""# Resolved-contact native binding\n\nStatus: **{STATUS}**.\n\nThe theorem is uniform over every nonnegative smooth resolved cutoff supported in `B_(M/4)`.  It does not select a causal plateau.  On a resolved-contact upward edge there is exactly one parent at `<=M/4`; the other lies in `(M/4,5M/4]`.  If the low-parent cutoff value is `q`, signed work first satisfies `dW=q dW+(1-q)dW`, with the two terms being actual mixed `V-h` and actual `h-h` work.  Only then is the already-canonical donor-restricted `dW+` restricted by the same nonnegative weights.\n\nThe mixed positive submeasure remains one canonical cause.  The same physical resolved operator is split signed-first as `I=K+S`, giving `mu_mixed<=I+<=K++S+`; no proportional owner matching is made.  `K` is conservative same-event relink provenance and `S` is existing strain/deformation.  The separate HH complement has both parent frequencies `<=5M/4`.  Because `0<=S<=1` also gives `||h||_2<=||u||_2`, its canonical positive submeasure first passes through the signed helicity-edge total variation; the exact variation/clean-Young ratio is `sqrt(2) C_F<1`, so the clean `R=5/4` sharp-Young capacity applies without an aggregate shell Hahn or an output-scale locality theorem.\n\nStress: `{out.samples}` admissible geometry/K-S states\n- maximum positive repartition residual: `{out.maximum_partition_residual:.3e}`\n- maximum parent/shell upper ratio: `{out.maximum_parent_upper_ratio:.12g}`\n- minimum K/S positive-cover margin: `{out.minimum_positive_ks_cover_margin:.3e}`\n- minimum half-owner margin: `{out.minimum_half_owner_margin:.3e}`\n- boundary-contact HH fraction: `{out.boundary_contact_hh_fraction:.12g}`\n- coarse-Hahn anti-theorem gap: `{out.coarse_hahn_counterexample_gap:.12g}`\n- edge-total-variation / clean-Young ratio: `{out.edge_variation_to_clean_young_ratio:.12g}`\n\nNo later Hahn split of the canonical cause, temporal deposit matching, recipient-shell reweighting, new event clock, or Navier--Stokes global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
