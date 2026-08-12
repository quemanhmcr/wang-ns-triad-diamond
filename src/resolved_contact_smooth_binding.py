from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.cyclic_helical_triad_donor_kernel import (
    CyclicTriadMeasureKernel,
    ClosedHelicalTriadRegistration,
    cyclic_triad_measure_kernel,
    register_closed_helical_triad,
)
from src.hard_tail_true_upward_supply import (
    HardTailUpwardSupplySplit,
    UpwardSupplyAtom,
    deep_upward_resolved_contact_fixture,
    hard_tail_upward_supply_split,
)

STATUS = (
    "EXACT_RESOLVED_CONTACT_SMOOTH_BINDING__POSITIVE_CORE_CUTOFF_BEFORE_CAUSAL_PUSH__"
    "DEEP_M_GE_8N_IS_PURE_MIXED_VH__M_EQ_4N_RETAINS_TRUE_TRANSITION_HH__"
    "DONOR_RESTRICTED_DWPLUS_BINDS_BEFORE_COARSE_HAHN__NONCLONING_KS_OWNER_SPLIT__"
    "COMMON_N_DW_UNIT"
)

CORE_FRACTION = 1.0 / 8.0
SUPPORT_FRACTION = 1.0 / 4.0


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


def canonical_positive_resolved_cutoff(radius: float, shell_scale: float) -> float:
    """One canonical C-infinity radial low-pass value.

    The represented multiplier is real, 0<=S<=1, equal to one on B_(M/8),
    and zero on and outside |k|=M/4.  Only these pointwise facts are used by the
    theorem.  The transition formula is the standard flat C-infinity step.
    """
    r = _finite_nonnegative(radius, "frequency radius")
    M = _finite_positive(shell_scale, "recipient shell scale")
    core = CORE_FRACTION * M
    outer = SUPPORT_FRACTION * M
    if r <= core:
        return 1.0
    if r >= outer:
        return 0.0
    t = (r - core) / (outer - core)
    a = math.exp(-1.0 / t)
    b = math.exp(-1.0 / (1.0 - t))
    q = b / (a + b)
    if not (0.0 < q < 1.0):
        raise AssertionError("canonical transition cutoff left (0,1)")
    return q


@dataclass(frozen=True)
class DeepContactSmoothRepartition:
    boundary: float
    recipient_shell_scale: float
    recipient_shell_index: int
    donor_radius: float
    other_parent_radius: float
    recipient_radius: float
    donor_cutoff_value: float
    other_parent_cutoff_value: float
    canonical_positive_mass: float
    low_low_bound_mass: float
    mixed_vh_bound_mass: float
    transition_hh_bound_mass: float
    common_unit_scale: float
    strict_deep_fully_mixed: bool
    borderline_transition_shell: bool
    canonical_cause_preserved: bool = True
    later_hahn_used: bool = False
    coarse_hahn_push_used: bool = False
    cloned_owner_mass: bool = False
    recipient_shell_reweighting_used: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        M = _finite_positive(self.recipient_shell_scale, "recipient shell scale")
        rd = _finite_positive(self.donor_radius, "donor radius")
        ro = _finite_positive(self.other_parent_radius, "other parent radius")
        rr = _finite_positive(self.recipient_radius, "recipient radius")
        m = _finite_positive(self.canonical_positive_mass, "canonical positive submeasure mass")
        qd = float(self.donor_cutoff_value)
        qo = float(self.other_parent_cutoff_value)
        if not all(math.isfinite(q) and 0.0 <= q <= 1.0 for q in (qd, qo)):
            raise ValueError("positive resolved cutoff must satisfy 0<=S<=1")
        if self.recipient_shell_index < 2 or abs(M - (2.0 ** self.recipient_shell_index) * N) > 5e-12 * M:
            raise ValueError("deep-contact binding requires a boundary-anchored shell M>=4N")
        if not (rd <= N and 0.5 * M < rr <= M):
            raise ValueError("deep-contact atom left its low-to-high donor/recipient geometry")
        if rd > 0.25 * M + 5e-12 * M:
            raise AssertionError("deep donor left the resolved-support radius")
        if not ro > 0.25 * M - 5e-12 * M:
            raise AssertionError("the non-donor parent of a deep upward edge must lie above M/4")
        if qo != 0.0:
            raise AssertionError("strict resolved cutoff must vanish on the other deep parent")
        for name, value in (
            ("low-low bound mass", self.low_low_bound_mass),
            ("mixed V-h bound mass", self.mixed_vh_bound_mass),
            ("transition HH bound mass", self.transition_hh_bound_mass),
        ):
            _finite_nonnegative(value, name)
        tol = 2e-12 * max(1.0, m)
        if abs(self.low_low_bound_mass) > tol:
            raise AssertionError("high recipient shell acquired impossible low-low canonical mass")
        if abs(self.mixed_vh_bound_mass - qd * m) > tol:
            raise AssertionError("mixed positive submeasure is not the cutoff restriction q_d dW+")
        if abs(self.transition_hh_bound_mass - (1.0 - qd) * m) > tol:
            raise AssertionError("transition HH submeasure is not (1-q_d)dW+")
        if abs(self.mixed_vh_bound_mass + self.transition_hh_bound_mass - m) > tol:
            raise AssertionError("smooth repartition cloned or lost canonical positive cause")
        if ro > rr + rd + 5.0e-12 * M:
            raise AssertionError("interaction-parent geometry left triad closure")
        if self.transition_hh_bound_mass > tol:
            if not (0.125 * M < rd <= 0.25 * M + 5.0e-12 * M):
                raise AssertionError("nonzero transition-HH mass must come from the genuine cutoff annulus M/8<|kd|<=M/4")
            if ro > 1.25 * M + 5.0e-12 * M:
                raise AssertionError("M=4N transition-HH other parent lost automatic <=5M/4 comparability")
        if self.strict_deep_fully_mixed != (self.recipient_shell_index >= 3):
            raise AssertionError("strict-deep label changed from M>=8N")
        if self.borderline_transition_shell != (self.recipient_shell_index == 2):
            raise AssertionError("borderline transition label changed from M=4N")
        if self.strict_deep_fully_mixed:
            if rd > 0.125 * M + 5e-12 * M:
                raise AssertionError("M>=8N donor should lie in the resolved core B_(M/8)")
            if qd != 1.0 or abs(self.mixed_vh_bound_mass - m) > tol or self.transition_hh_bound_mass > tol:
                raise AssertionError("strict-deep upward atom failed to become pure mixed V-h work")
            if not ro > 0.25 * M:
                raise AssertionError("strict-deep other parent failed the exact >M/4 UV separation")
        if self.common_unit_scale != N:
            raise ValueError("causal unit must remain the parent-tail N dW unit")
        if (
            not self.canonical_cause_preserved
            or self.later_hahn_used
            or self.coarse_hahn_push_used
            or self.cloned_owner_mass
            or self.recipient_shell_reweighting_used
        ):
            raise ValueError("smooth binding may only restrict/push the canonical cause in the unchanged N dW unit")

    @property
    def mixed_common_unit_mass(self) -> float:
        return self.common_unit_scale * self.mixed_vh_bound_mass

    @property
    def transition_hh_common_unit_mass(self) -> float:
        return self.common_unit_scale * self.transition_hh_bound_mass


def deep_contact_smooth_repartition(atom: UpwardSupplyAtom) -> DeepContactSmoothRepartition:
    """Push one donor-restricted canonical upward atom through the actual V+h split.

    On a deep upward edge the donor is the unique parent at or below M/4.  The
    other parent is strictly above M/4 by triangle closure and therefore is pure
    h.  For the canonical positive cutoff, the edge work has the signed identity

        dW = q_d dW + (1-q_d) dW,

    where the first term is the actual mixed V-h interaction and the second is
    actual h-h work.  Since the input atom is already a positive submeasure of
    canonical dW+, multiplying by q_d and 1-q_d is restriction/pushforward, not
    a new Hahn decomposition.
    """
    if not atom.deep_upward_shell or not atom.resolved_scale_parent_contact:
        raise ValueError("deep resolved-contact upward atom required")
    M = atom.recipient_shell_scale
    parent_radii = tuple(float(r) for r in atom.interaction_parent_radii)
    low = min(parent_radii)
    high = max(parent_radii)
    tol = 5e-12 * M
    if abs(low - atom.donor_radius) > tol:
        raise AssertionError("deep cyclic energy donor is not the unique resolved-scale interaction parent")
    if not high > 0.25 * M:
        raise AssertionError("triangle closure should force the other deep parent above M/4")
    qd = canonical_positive_resolved_cutoff(atom.donor_radius, M)
    qo = canonical_positive_resolved_cutoff(high, M)
    m = atom.physical_work_mass
    return DeepContactSmoothRepartition(
        boundary=atom.boundary,
        recipient_shell_scale=M,
        recipient_shell_index=atom.recipient_shell_index,
        donor_radius=atom.donor_radius,
        other_parent_radius=high,
        recipient_radius=atom.recipient_radius,
        donor_cutoff_value=qd,
        other_parent_cutoff_value=qo,
        canonical_positive_mass=m,
        low_low_bound_mass=0.0,
        mixed_vh_bound_mass=qd * m,
        transition_hh_bound_mass=(1.0 - qd) * m,
        common_unit_scale=atom.common_unit_scale,
        strict_deep_fully_mixed=(atom.recipient_shell_index >= 3),
        borderline_transition_shell=(atom.recipient_shell_index == 2),
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
        vals = tuple(float(v) for v in (self.signed_mixed_work, self.signed_skew_work, self.signed_strain_work))
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("finite signed resolved K/S work required")
        scale = max(1.0, *(abs(v) for v in vals))
        if abs(vals[0] - vals[1] - vals[2]) > 8e-12 * scale:
            raise ValueError("signed identity must be verified before any positive K/S binding")
        if not (self.same_physical_atom and self.same_resolved_operator and self.observer_gauge_quotiented_or_fixed_event):
            raise ValueError("K/S binding refuses representation substitution or unquotiented observer motion")


@dataclass(frozen=True)
class PositiveKSSubmeasureBinding:
    canonical_mixed_submeasure_mass: float
    available_positive_mixed_work: float
    available_positive_skew_work: float
    available_positive_strain_work: float
    skew_bound_mass: float
    strain_bound_mass: float
    common_unit_scale: float
    signed_identity_residual: float
    canonical_mass_residual: float
    maximum_domination_excess: float
    canonical_cause_replaced: bool = False
    later_hahn_on_canonical_cause: bool = False
    owner_mass_cloned: bool = False
    own_shell_reweighting_used: bool = False

    def __post_init__(self) -> None:
        mu = _finite_nonnegative(self.canonical_mixed_submeasure_mass, "canonical mixed submeasure mass")
        I = _finite_nonnegative(self.available_positive_mixed_work, "positive mixed work")
        K = _finite_nonnegative(self.available_positive_skew_work, "positive skew work")
        S = _finite_nonnegative(self.available_positive_strain_work, "positive strain work")
        mk = _finite_nonnegative(self.skew_bound_mass, "skew-bound canonical mass")
        ms = _finite_nonnegative(self.strain_bound_mass, "strain-bound canonical mass")
        scale = max(1.0, mu, I, K, S)
        tol = 2e-11 * scale
        if mu > I + tol:
            raise AssertionError("donor-restricted canonical cause exceeds the same-atom positive mixed work")
        if abs(mk + ms - mu) > tol or self.canonical_mass_residual > 2e-11:
            raise AssertionError("K/S binding cloned or lost canonical mixed cause")
        if mk > K + tol or ms > S + tol or self.maximum_domination_excess > tol:
            raise AssertionError("bound canonical submeasure exceeds an existing component positive law")
        _finite_positive(self.common_unit_scale, "common N work scale")
        if self.signed_identity_residual > 2e-11:
            raise AssertionError("K/S binding left the signed physical identity")
        if self.canonical_cause_replaced or self.later_hahn_on_canonical_cause or self.owner_mass_cloned or self.own_shell_reweighting_used:
            raise ValueError("K/S provenance may not replace, re-Hahn, clone, or M-reweight canonical dW+")


def bind_canonical_mixed_submeasure_to_ks(
    canonical_mixed_submeasure_mass: float,
    signed_atom: SignedResolvedKSAtom,
    *,
    common_unit_scale: float,
) -> PositiveKSSubmeasureBinding:
    """Non-cloning positive domination of canonical dW+ by existing K+/S+ laws.

    This function does *not* Hahn-split the canonical cause.  The input ``mu`` is
    already a positive restriction of canonical dW+.  K and S have their own
    signed physical laws on the same refined work atom.  From I=K+S,
    I+<=K++S+.  The symmetric common-unit allocation below restricts exactly mu
    into submeasures dominated by those two existing positive component laws.
    """
    mu = _finite_nonnegative(canonical_mixed_submeasure_mass, "canonical mixed submeasure mass")
    N = _finite_positive(common_unit_scale, "common N work scale")
    I = float(signed_atom.signed_mixed_work)
    K = float(signed_atom.signed_skew_work)
    S = float(signed_atom.signed_strain_work)
    ip = max(I, 0.0)
    kp = max(K, 0.0)
    sp = max(S, 0.0)
    scale = max(1.0, abs(I), abs(K), abs(S), mu)
    if mu > ip + 8e-12 * scale:
        raise ValueError("canonical mixed submeasure is not dominated by same-atom mixed dW+")
    capacity = kp + sp
    if mu == 0.0:
        mk = ms = 0.0
    else:
        if capacity <= 0.0:
            raise AssertionError("positive mixed cause has no K+/S+ cover")
        theta = mu / capacity
        if theta > 1.0 + 8e-12:
            raise AssertionError("signed K/S positive cover is smaller than canonical mixed submeasure")
        mk = theta * kp
        ms = theta * sp
    identity_residual = abs(I - K - S) / scale
    mass_residual = abs(mk + ms - mu) / scale
    domination_excess = max(0.0, mk - kp, ms - sp)
    return PositiveKSSubmeasureBinding(
        canonical_mixed_submeasure_mass=mu,
        available_positive_mixed_work=ip,
        available_positive_skew_work=kp,
        available_positive_strain_work=sp,
        skew_bound_mass=mk,
        strain_bound_mass=ms,
        common_unit_scale=N,
        signed_identity_residual=identity_residual,
        canonical_mass_residual=mass_residual,
        maximum_domination_excess=domination_excess,
    )


def coarse_hahn_cancellation_counterexample() -> dict[str, float | bool]:
    """Atomic positive cause can exceed the Hahn mass after downstream aggregation."""
    first = 1.0
    second = -0.9
    coarse = first + second
    return {
        "canonical_positive_first_atom": first,
        "second_signed_atom": second,
        "coarse_signed_work": coarse,
        "coarse_positive_hahn_mass": max(coarse, 0.0),
        "atomic_cause_exceeds_coarse_hahn": first > max(coarse, 0.0),
    }


def interior_transition_resolved_contact_fixture() -> tuple[
    ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel, HardTailUpwardSupplySplit
]:
    """Physical M=4N fixture with the radial donor strictly inside the smooth transition.

    Fix N=sqrt(2), hence M=4sqrt(2), and the closed geometry

        (-3,-2,0) + (1,0,0) + (2,2,0) = 0.

    The radius-one donor satisfies M/8<1<M/4, while the recipient radius
    sqrt(13) lies strictly in (M/2,M].  A finite helicity/phase search only
    orients the actual Waleffe work so that this low root donates to the deep
    recipient; it does not change cutoff geometry.
    """
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
        complex(-2.0**-0.5, 2.0**-0.5),
        complex(-2.0**-0.5, -2.0**-0.5),
    )
    for s0 in (-1, 1):
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                for phase in phases:
                    for position in range(3):
                        amps = [1.0 + 0.0j, 1.0 + 0.0j, 1.0 + 0.0j]
                        amps[position] = phase
                        triad = register_closed_helical_triad(
                            wavevectors=wavevectors,
                            helicities=(s0, s1, s2),
                            amplitudes=tuple(amps),
                        )
                        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
                        if not kernel.numerically_resolved_transport:
                            continue
                        try:
                            split = hard_tail_upward_supply_split(triad, kernel, boundary=N)
                        except ValueError:
                            continue
                        targets = [
                            atom
                            for atom in split.atoms
                            if atom.recipient_shell_index == 2
                            and abs(atom.donor_radius - 1.0) <= 5.0e-12
                        ]
                        if targets:
                            for atom in targets:
                                binding = deep_contact_smooth_repartition(atom)
                                if not (0.0 < binding.donor_cutoff_value < 1.0):
                                    raise AssertionError("interior-transition fixture failed 0<q_d<1")
                                if not (binding.mixed_vh_bound_mass > 0.0 and binding.transition_hh_bound_mass > 0.0):
                                    raise AssertionError("interior-transition canonical cause did not split into two positive physical pieces")
                            return triad, kernel, split
    raise AssertionError("finite physical helicity/phase search found no interior-transition M=4N upward atom")


def strict_deep_resolved_mixed_fixture() -> tuple[
    ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel, HardTailUpwardSupplySplit
]:
    """Deterministic finite search for a true upward atom into M=8N.

    Geometry is fixed before phases/helicities are searched:

        (-5,-4,0) + (1,0,0) + (4,4,0) = 0.

    At N=1 the radius-one root is eligible as a radial donor while the two high
    radii sqrt(32), sqrt(41) lie in the recipient shell M=8.  We search only the
    finite physical helicity/phase orientation needed to make the radius-one
    donor feed a positive high recipient; support geometry is unchanged.
    """
    wavevectors = (
        np.asarray((-5.0, -4.0, 0.0)),
        np.asarray((1.0, 0.0, 0.0)),
        np.asarray((4.0, 4.0, 0.0)),
    )
    phases = (
        1.0 + 0.0j,
        -1.0 + 0.0j,
        0.0 + 1.0j,
        0.0 - 1.0j,
        complex(2.0**-0.5, 2.0**-0.5),
        complex(2.0**-0.5, -2.0**-0.5),
        complex(-2.0**-0.5, 2.0**-0.5),
        complex(-2.0**-0.5, -2.0**-0.5),
    )
    for s0 in (-1, 1):
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                for phase in phases:
                    for position in range(3):
                        amps = [1.0 + 0.0j, 1.0 + 0.0j, 1.0 + 0.0j]
                        amps[position] = phase
                        triad = register_closed_helical_triad(
                            wavevectors=wavevectors,
                            helicities=(s0, s1, s2),
                            amplitudes=tuple(amps),
                        )
                        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
                        if not kernel.numerically_resolved_transport:
                            continue
                        try:
                            split = hard_tail_upward_supply_split(triad, kernel, boundary=1.0)
                        except ValueError:
                            continue
                        strict = [
                            atom
                            for atom in split.atoms
                            if atom.deep_upward_shell
                            and atom.recipient_shell_index >= 3
                            and abs(atom.donor_radius - 1.0) <= 5.0e-12
                        ]
                        if strict:
                            for atom in strict:
                                binding = deep_contact_smooth_repartition(atom)
                                if not binding.strict_deep_fully_mixed:
                                    raise AssertionError("strict-deep fixture lost pure mixed V-h binding")
                            return triad, kernel, split
    raise AssertionError("finite physical helicity/phase search found no M=8N low-to-high strict-deep atom")


def boundary_contact_counterexample() -> dict[str, float | bool]:
    triad, kernel, split = deep_upward_resolved_contact_fixture()
    del triad, kernel
    deep = [a for a in split.atoms if a.deep_upward_shell]
    if not deep:
        raise AssertionError("certified deep fixture lost its deep atom")
    atom = min(deep, key=lambda a: abs(a.donor_radius - 0.25 * a.recipient_shell_scale))
    binding = deep_contact_smooth_repartition(atom)
    return {
        "donor_to_shell_ratio": atom.donor_radius / atom.recipient_shell_scale,
        "cutoff_value": binding.donor_cutoff_value,
        "mixed_fraction": binding.mixed_vh_bound_mass / binding.canonical_positive_mass,
        "hh_fraction": binding.transition_hh_bound_mass / binding.canonical_positive_mass,
        "support_contact_is_not_interface": binding.mixed_vh_bound_mass == 0.0 and binding.transition_hh_bound_mass > 0.0,
    }


def interior_transition_fixture_observation() -> dict[str, float | bool]:
    triad, kernel, split = interior_transition_resolved_contact_fixture()
    del triad, kernel
    targets = [
        atom
        for atom in split.atoms
        if atom.recipient_shell_index == 2 and abs(atom.donor_radius - 1.0) <= 5.0e-12
    ]
    if not targets:
        raise AssertionError("interior-transition fixture lost its target atom")
    atom = max(targets, key=lambda a: a.physical_work_mass)
    binding = deep_contact_smooth_repartition(atom)
    return {
        "donor_to_shell_ratio": atom.donor_radius / atom.recipient_shell_scale,
        "cutoff_value": binding.donor_cutoff_value,
        "mixed_fraction": binding.mixed_vh_bound_mass / binding.canonical_positive_mass,
        "hh_fraction": binding.transition_hh_bound_mass / binding.canonical_positive_mass,
        "both_positive": binding.mixed_vh_bound_mass > 0.0 and binding.transition_hh_bound_mass > 0.0,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "signed_first": "for each deep upward recipient edge, dW_mix=q_d dW and dW_HH=(1-q_d)dW, so dW=dW_mix+dW_HH before any new Hahn operation",
        "positive_push": "donor-restricted canonical dW+ is multiplied only by q_d and 1-q_d; the two positive submeasures sum to the original cause",
        "strict_deep": "M>=8N forces donor<=M/8, q_d=1, other parent>M/4, hence the entire canonical upward atom is actual mixed V-h work",
        "borderline": "M=4N is a genuine transition shell; contact at M/4 has q_d=0 and can remain 100% HH; whenever transition-HH mass is nonzero its donor lies in (M/8,M/4] and the other parent lies in (M/4,5M/4], so this HH piece is already comparable",
        "ks_binding": "only after same-atom signed I=K+S is verified does a non-cloning positive domination bind canonical mixed dW+ into existing K+/S+ owner laws",
        "coarse_hahn_forbidden": "canonical atomic dW+ is bound before downstream aggregation; a later coarse Hahn mass may be smaller by cancellation and cannot receive the cause",
        "common_unit": "all bound masses retain N dW; recipient M is geometry only",
        "cutoff": "choose a canonical smooth real radial cutoff 0<=S<=1, S=1 on B_(M/8), S=0 on and outside M/4",
        "later_hahn_used": False,
        "cause_cloning_used": False,
        "recipient_shell_reweighting_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class ResolvedContactSmoothBindingStress:
    status: str
    samples: int
    strict_deep_samples: int
    borderline_samples: int
    maximum_partition_residual: float
    maximum_strict_deep_hh_fraction: float
    minimum_strict_deep_mixed_fraction: float
    maximum_ks_mass_residual: float
    maximum_ks_domination_excess: float
    coarse_hahn_counterexample_gap: float
    boundary_contact_mixed_fraction: float
    boundary_contact_hh_fraction: float
    interior_transition_mixed_fraction: float
    interior_transition_hh_fraction: float


def stress(samples: int = 50_000, seed: int = 2026081301) -> ResolvedContactSmoothBindingStress:
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = np.random.default_rng(int(seed))
    strict = borderline = 0
    max_partition = 0.0
    max_strict_hh = 0.0
    min_strict_mixed = 1.0
    max_ks_mass = 0.0
    max_domination = 0.0
    for _ in range(count):
        N = float(2.0 ** rng.uniform(-4.0, 4.0))
        j = int(rng.integers(2, 7))
        M = (2.0 ** j) * N
        rd = float(rng.uniform(0.02, 1.0)) * N
        q = canonical_positive_resolved_cutoff(rd, M)
        m = float(10.0 ** rng.uniform(-6.0, 3.0))
        mixed = q * m
        hh = (1.0 - q) * m
        max_partition = max(max_partition, abs(mixed + hh - m) / max(1.0, m))
        if j >= 3:
            strict += 1
            max_strict_hh = max(max_strict_hh, hh / m)
            min_strict_mixed = min(min_strict_mixed, mixed / m)
            if q != 1.0:
                raise AssertionError("M>=8N random donor escaped canonical cutoff core")
        else:
            borderline += 1

        I = float(10.0 ** rng.uniform(-5.0, 3.0))
        K = float(rng.uniform(-2.5, 3.5)) * I
        S = I - K
        mu = float(rng.uniform(0.0, 1.0)) * I
        bind = bind_canonical_mixed_submeasure_to_ks(
            mu,
            SignedResolvedKSAtom(I, K, S),
            common_unit_scale=N,
        )
        max_ks_mass = max(max_ks_mass, bind.canonical_mass_residual)
        max_domination = max(max_domination, bind.maximum_domination_excess)

    coarse = coarse_hahn_cancellation_counterexample()
    boundary = boundary_contact_counterexample()
    interior = interior_transition_fixture_observation()
    return ResolvedContactSmoothBindingStress(
        status=STATUS,
        samples=count,
        strict_deep_samples=strict,
        borderline_samples=borderline,
        maximum_partition_residual=max_partition,
        maximum_strict_deep_hh_fraction=max_strict_hh,
        minimum_strict_deep_mixed_fraction=min_strict_mixed if strict else 1.0,
        maximum_ks_mass_residual=max_ks_mass,
        maximum_ks_domination_excess=max_domination,
        coarse_hahn_counterexample_gap=float(coarse["canonical_positive_first_atom"] - coarse["coarse_positive_hahn_mass"]),
        boundary_contact_mixed_fraction=float(boundary["mixed_fraction"]),
        boundary_contact_hh_fraction=float(boundary["hh_fraction"]),
        interior_transition_mixed_fraction=float(interior["mixed_fraction"]),
        interior_transition_hh_fraction=float(interior["hh_fraction"]),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081301)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-contact-smooth-binding"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples, args.seed)
    cert = theorem_certificate()
    (args.outdir / "resolved_contact_smooth_binding.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2, sort_keys=True) + "\n"
    )
    boundary = boundary_contact_counterexample()
    interior = interior_transition_fixture_observation()
    coarse = coarse_hahn_cancellation_counterexample()
    md = f"""# Resolved-contact smooth binding\n\nStatus: **{STATUS}**.\n\nThe actual positive cutoff is chosen with `0<=S<=1`, `S=1` on `B_(M/8)`, and `S=0` on and outside `M/4`. For a deep upward atom the cyclic donor is the unique parent at or below `M/4`; the other parent is strictly above `M/4`. Thus before any new Hahn operation,\n\n`dW = S(k_d) dW_mixed + (1-S(k_d)) dW_HH`,\n\nand the donor-restricted canonical `dW+` is pushed by these same two nonnegative weights. No mass is cloned and the common causal unit remains `N dW`.\n\nFor `M>=8N`, the donor lies in `B_(M/8)`, so `S(k_d)=1`: the whole canonical upward atom is actual mixed `V-h` work. The shell `M=4N` is genuinely different. The certified boundary-contact fixture has donor/shell ratio `{boundary['donor_to_shell_ratio']:.12g}`, cutoff value `{boundary['cutoff_value']:.12g}`, mixed fraction `{boundary['mixed_fraction']:.12g}`, and HH fraction `{boundary['hh_fraction']:.12g}`. Hence support contact alone is not interface ownership.\n\nOn an actual mixed atom, first verify the signed identity `I=K+S` from the same resolved operator. Only then bind the already-canonical positive submeasure into positive K/S owner submeasures by domination; their masses sum exactly to the incoming canonical cause. A downstream coarse Hahn is forbidden as the transport map: the fixed cancellation counterexample leaves coarse positive mass `{coarse['coarse_positive_hahn_mass']:.12g}` from an atomic canonical cause `{coarse['canonical_positive_first_atom']:.12g}`.\n\nStress: `{out.samples}` states\n- strict-deep samples: `{out.strict_deep_samples}`\n- borderline samples: `{out.borderline_samples}`\n- maximum smooth partition residual: `{out.maximum_partition_residual:.3e}`\n- maximum strict-deep HH fraction: `{out.maximum_strict_deep_hh_fraction:.3e}`\n- minimum strict-deep mixed fraction: `{out.minimum_strict_deep_mixed_fraction:.12g}`\n- maximum K/S canonical-mass residual: `{out.maximum_ks_mass_residual:.3e}`\n- maximum K/S domination excess: `{out.maximum_ks_domination_excess:.3e}`\n\nNo output-scale causal reweighting, new clock, new event depth, or Navier--Stokes global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(json.dumps(asdict(out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
