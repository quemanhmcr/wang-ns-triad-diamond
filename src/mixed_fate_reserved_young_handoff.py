from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.canonical_positive_edge_work_routing import (
    ETA0,
    CanonicalPositiveEdgeWorkRouting,
    HardCellWork,
)
from src.continuum_helical_edge_measure_registration import (
    ContinuumEdgeMeasureLedger,
    _replay_physical_ledger,
)
from src.single_edge_certificate import RSTAR_HI, RSTAR_LO

STATUS = (
    "EXACT_MIXED_FATE_RESERVED_YOUNG_HANDOFF__"
    "CANONICAL_DW_PLUS_INHERITED__FULL_SIGNED_HARD_CELL_ONLY__"
    "TERMINAL_BAD_ASSISTANCE_RESERVED_WITHOUT_NEW_TRILINEAR_WORK"
)


def _close(a: float, b: float, *, factor: float = 2.0e-10) -> bool:
    return abs(float(a) - float(b)) <= factor * max(abs(float(a)), abs(float(b)), 1.0e-300)


def inherited_negative_work(cell: HardCellWork) -> float:
    """Pushforward mass of the already-fixed canonical Hahn-negative edge law.

    Since ``pi_# dW = pi_# dW+ - pi_# dW-`` on every hard product cell,

      n_C = (pi_# dW-)(C) = P_C - T_C.

    This is not a new Hahn split after aggregation.  It is the negative mass of
    the original edge-space Hahn decomposition, recovered from its two verified
    pushforwards.
    """
    raw = float(cell.inherited_positive_work) - float(cell.signed_work)
    tol = 3.0e-10 * max(
        abs(float(cell.inherited_positive_work)),
        abs(float(cell.signed_work)),
        1.0e-300,
    )
    if raw < -tol:
        raise AssertionError("hard cell violates pi_#dW = pi_#dW+ - pi_#dW-")
    return 0.0 if raw < 0.0 else raw


@dataclass(frozen=True)
class ReservedYoungHandoff:
    signed_full_cell_work: float
    inherited_good_positive_work: float
    inherited_bad_positive_work: float
    inherited_negative_work: float
    certified_full_cell_young_upper: float
    reservation_value: float
    full_signed_young_deficit: float
    terminal_bad_assistance_fraction: float
    reserved_young_deficit: float
    normalized_symbol_freezing_error: float
    christ_modulus: float
    full_signed_christ_gate: bool
    reserved_christ_gate: bool
    young_handoff_certified: bool
    full_signed_trilinear_work_used: bool = True
    reservation_value_used_as_trilinear_work: bool = False
    inherited_good_causal_mass_changed: bool = False
    later_hahn_used_as_causal_law: bool = False

    def __post_init__(self) -> None:
        numeric = (
            self.signed_full_cell_work,
            self.inherited_good_positive_work,
            self.inherited_bad_positive_work,
            self.inherited_negative_work,
            self.certified_full_cell_young_upper,
            self.reservation_value,
            self.full_signed_young_deficit,
            self.terminal_bad_assistance_fraction,
            self.reserved_young_deficit,
            self.normalized_symbol_freezing_error,
            self.christ_modulus,
        )
        if not all(math.isfinite(float(v)) for v in numeric):
            raise ValueError("finite mixed-fate Young handoff data required")
        if self.certified_full_cell_young_upper <= 0.0:
            raise ValueError("positive full-cell Young upper required")
        if min(
            self.inherited_good_positive_work,
            self.inherited_bad_positive_work,
            self.inherited_negative_work,
            self.terminal_bad_assistance_fraction,
            self.normalized_symbol_freezing_error,
            self.christ_modulus,
        ) < 0.0:
            raise ValueError("nonnegative inherited masses/errors and positive modulus required")
        if not (0.0 < self.christ_modulus < 1.0):
            raise ValueError("external Christ near-extremal modulus must lie in (0,1)")
        if (
            not self.full_signed_trilinear_work_used
            or self.reservation_value_used_as_trilinear_work
            or self.inherited_good_causal_mass_changed
            or self.later_hahn_used_as_causal_law
        ):
            raise ValueError("reserved handoff changed signed-work or canonical-causality semantics")
        expected_reservation = self.signed_full_cell_work - self.inherited_bad_positive_work
        if not _close(self.reservation_value, expected_reservation, factor=4.0e-10):
            raise AssertionError("reservation value is not T_C-b_C")
        expected_signed = (
            self.inherited_good_positive_work
            + self.inherited_bad_positive_work
            - self.inherited_negative_work
        )
        if not _close(self.signed_full_cell_work, expected_signed, factor=5.0e-10):
            raise AssertionError("mixed-fate Hahn pushforwards do not reconstruct full signed cell work")
        Y = self.certified_full_cell_young_upper
        if abs(self.signed_full_cell_work) > Y + 5.0e-10 * max(Y, abs(self.signed_full_cell_work)):
            raise AssertionError("certified Young upper does not bound the full signed hard-cell work")
        full_deficit = 1.0 - self.signed_full_cell_work / Y
        bad_fraction = self.inherited_bad_positive_work / Y
        reserved_deficit = 1.0 - self.reservation_value / Y
        if not _close(self.full_signed_young_deficit, full_deficit, factor=5.0e-10):
            raise AssertionError("full signed Young deficit changed")
        if not _close(self.terminal_bad_assistance_fraction, bad_fraction, factor=5.0e-10):
            raise AssertionError("bad assistance fraction changed")
        if not _close(self.reserved_young_deficit, reserved_deficit, factor=5.0e-10):
            raise AssertionError("reserved Young deficit changed")
        if not _close(
            self.reserved_young_deficit,
            self.full_signed_young_deficit + self.terminal_bad_assistance_fraction,
            factor=8.0e-10,
        ):
            raise AssertionError("reserved deficit failed exact full-deficit + bad-assistance identity")
        expected_full_gate = (
            self.signed_full_cell_work > 0.0
            and self.full_signed_young_deficit + self.normalized_symbol_freezing_error <= self.christ_modulus
        )
        expected_reserved_gate = (
            self.inherited_good_positive_work > 0.0
            and self.reservation_value > 0.0
            and self.reserved_young_deficit + self.normalized_symbol_freezing_error <= self.christ_modulus
        )
        if self.full_signed_christ_gate != expected_full_gate:
            raise AssertionError("full signed Christ gate changed")
        if self.reserved_christ_gate != expected_reserved_gate:
            raise AssertionError("reserved Christ gate changed")
        if self.young_handoff_certified != self.reserved_christ_gate:
            raise AssertionError("only the bad-assistance-free reserved gate may certify mixed-fate Young handoff")
        if self.reserved_christ_gate and not self.full_signed_christ_gate:
            raise AssertionError("reserved gate cannot pass unless the actual full signed cell also passes Christ")


def certify_reserved_young_handoff(
    cell: HardCellWork,
    *,
    certified_full_cell_young_upper: float,
    normalized_symbol_freezing_error: float,
    christ_modulus: float,
) -> ReservedYoungHandoff:
    """Certify that terminal bad-positive work is not needed for Young/Christ.

    The scalar ``reservation_value=T_C-b_C`` is *never* passed to Young and is
    never called physical trilinear work.  It is only a sufficient-condition
    certificate.  Young/Christ still sees the actual full signed hard-cell form
    ``T_C``.  The canonical good ``dW+`` mass is inherited unchanged.
    """
    Y = float(certified_full_cell_young_upper)
    xi = float(normalized_symbol_freezing_error)
    modulus = float(christ_modulus)
    if not all(math.isfinite(v) for v in (Y, xi, modulus)):
        raise ValueError("finite Young upper, symbol error and Christ modulus required")
    if Y <= 0.0 or xi < 0.0 or not (0.0 < modulus < 1.0):
        raise ValueError("positive Young upper, nonnegative symbol error, and Christ modulus in (0,1) required")
    if not cell.inherited_good_positive_work > 0.0:
        raise ValueError("Young handoff requires inherited geometry-good canonical dW+ mass")
    n = inherited_negative_work(cell)
    T = float(cell.signed_work)
    b = float(cell.inherited_bad_positive_work)
    reservation = T - b
    full_deficit = 1.0 - T / Y
    bad_fraction = b / Y
    reserved_deficit = 1.0 - reservation / Y
    full_gate = T > 0.0 and full_deficit + xi <= modulus
    reserved_gate = reservation > 0.0 and reserved_deficit + xi <= modulus
    return ReservedYoungHandoff(
        signed_full_cell_work=T,
        inherited_good_positive_work=float(cell.inherited_good_positive_work),
        inherited_bad_positive_work=b,
        inherited_negative_work=n,
        certified_full_cell_young_upper=Y,
        reservation_value=reservation,
        full_signed_young_deficit=full_deficit,
        terminal_bad_assistance_fraction=bad_fraction,
        reserved_young_deficit=reserved_deficit,
        normalized_symbol_freezing_error=xi,
        christ_modulus=modulus,
        full_signed_christ_gate=full_gate,
        reserved_christ_gate=reserved_gate,
        young_handoff_certified=reserved_gate,
    )


CERTIFIED_GOOD_WORK_CAPACITY_DENSITY_LOWER = Fraction(19, 100)
ACTUAL_NS_GOOD_RATIO_SQUARED = Fraction(41, 110)


def good_work_capacity_density_lower(eta: float = ETA0) -> float:
    """Clean rigorous pointwise lower ``dW/dA >= 19/100`` on the eta0 core.

    The analytic quantity is ``J*(1-eta0)/(gamma*+1/80)``.  The dedicated Arb
    certificate below proves that it is strictly larger than the conservative
    rational 19/100 using the already-certified bracket for ``r*``.  Runtime
    work comparisons use only this clean lower, never an optimizer decimal.
    """
    e = float(eta)
    if not math.isfinite(e) or not _close(e, ETA0, factor=2.0e-14):
        raise ValueError("clean physical RN lower is certified specifically on the canonical eta0=1e-4 core")
    return float(CERTIFIED_GOOD_WORK_CAPACITY_DENSITY_LOWER)


def low_deficit_contamination_ratio_upper(block_transfer_deficit: float, eta: float = ETA0) -> float:
    """Universal actual-work bound ``(B+N)/G`` on one low-deficit block.

    Capacity is only the reference used in the proof:

      A(G^c)/A <= epsilon/eta0,
      |dW| <= dA on G^c,
      dW/dA >= 19/100 on G.

    Thus, for ``epsilon<eta0``, with ``q=epsilon/eta0``,

      (bad-positive work + negative work) / good-positive work
      <= q / ((19/100)(1-q)).

    This is not the forbidden inference "capacity majority => work majority";
    the conversion uses an explicit pointwise physical RN lower on the same good
    edge law and the native total-variation bound outside it.
    """
    eps = float(block_transfer_deficit)
    e = float(eta)
    if not math.isfinite(eps) or eps < 0.0 or not math.isfinite(e) or not _close(e, ETA0, factor=2.0e-14):
        raise ValueError("finite nonnegative deficit on the canonical eta0 core required")
    if eps >= e:
        return math.inf
    q = eps / e
    if q == 0.0:
        return 0.0
    return q / (good_work_capacity_density_lower(e) * (1.0 - q))


@dataclass(frozen=True)
class LowDeficitWorkContamination:
    block_transfer_deficit: float
    good_work: float
    bad_positive_work: float
    negative_work: float
    contamination_work: float
    actual_contamination_to_good_ratio: float
    universal_contamination_to_good_ratio_upper: float
    canonical_good_mass_changed: bool = False
    capacity_used_as_causal_law: bool = False

    def __post_init__(self) -> None:
        numeric = (
            self.block_transfer_deficit,
            self.good_work,
            self.bad_positive_work,
            self.negative_work,
            self.contamination_work,
            self.actual_contamination_to_good_ratio,
            self.universal_contamination_to_good_ratio_upper,
        )
        if not all(math.isfinite(float(v)) for v in numeric):
            raise ValueError("finite low-deficit actual-work contamination data required")
        if min(numeric) < 0.0 or self.good_work <= 0.0:
            raise ValueError("nonnegative work/deficit data and positive good work required")
        if self.canonical_good_mass_changed or self.capacity_used_as_causal_law:
            raise ValueError("low-deficit comparison may not change canonical dW+ or promote capacity")
        if not _close(self.contamination_work, self.bad_positive_work + self.negative_work, factor=5.0e-10):
            raise AssertionError("contamination work is not bad-positive plus canonical negative work")
        if not _close(
            self.actual_contamination_to_good_ratio,
            self.contamination_work / self.good_work,
            factor=5.0e-10,
        ):
            raise AssertionError("actual contamination/good ratio changed")
        if self.actual_contamination_to_good_ratio > self.universal_contamination_to_good_ratio_upper + 2.0e-8 * max(
            self.universal_contamination_to_good_ratio_upper, 1.0e-300
        ):
            raise AssertionError("realized actual-work contamination exceeded the native low-deficit theorem")


def certify_low_deficit_work_contamination(
    ledger: ContinuumEdgeMeasureLedger,
    routing: CanonicalPositiveEdgeWorkRouting,
) -> LowDeficitWorkContamination:
    """Bind the low-native-deficit estimate to the same canonical physical law."""
    replayed = _replay_physical_ledger(ledger)
    if not _close(routing.total_positive_work, replayed.positive_edge_work, factor=8.0e-10):
        raise AssertionError("routing and low-deficit certificate are not bound to the same dW+ law")
    if not _close(routing.good_positive_work, routing.hard_cell_compression.inherited_good_positive_work, factor=8.0e-10):
        raise AssertionError("routing good dW+ and hard pushforward disagree")
    if not _close(routing.bad_positive_work, routing.hard_cell_compression.inherited_bad_positive_work, factor=8.0e-10):
        raise AssertionError("routing bad dW+ and hard pushforward disagree")

    negative = math.fsum(inherited_negative_work(cell) for cell in routing.hard_cell_compression.cells)
    if not _close(negative, replayed.negative_edge_work, factor=8.0e-10):
        raise AssertionError("hard cells do not push forward the canonical edge Hahn-negative mass")
    good = float(routing.good_positive_work)
    bad = float(routing.bad_positive_work)
    if good <= 0.0:
        raise ValueError("low-deficit contamination theorem requires nonzero good canonical work")
    contamination = bad + negative
    eps = float(replayed.block_transfer_deficit)
    ratio_upper = low_deficit_contamination_ratio_upper(eps, ETA0)
    if not math.isfinite(ratio_upper):
        raise ValueError("quantitative contamination theorem requires native block deficit < eta0")
    return LowDeficitWorkContamination(
        block_transfer_deficit=eps,
        good_work=good,
        bad_positive_work=bad,
        negative_work=negative,
        contamination_work=contamination,
        actual_contamination_to_good_ratio=contamination / good,
        universal_contamination_to_good_ratio_upper=ratio_upper,
    )


def reserved_failure_good_work_upper(
    handoff: ReservedYoungHandoff,
    *,
    full_signed_christ_margin_floor: float,
) -> float:
    """Bound good work in a cell that still fails reservation despite margin.

    Put ``mu=margin_floor``.  If the ordinary full-signed Young/Christ premise has
    spare margin at least ``mu`` but the reserved gate fails, then the exact
    reservation identity forces ``b_C/Y_C > mu``.  Since ``T_C<=Y_C`` and
    ``T_C=g_C+b_C-n_C``, one gets

      g_C < n_C + (1/mu - 1) b_C.

    The right side consists only of already-existing physical work laws.
    """
    mu = float(full_signed_christ_margin_floor)
    if not math.isfinite(mu) or not (0.0 < mu < 1.0):
        raise ValueError("Christ safety margin floor must lie in (0,1)")
    full_margin = handoff.christ_modulus - (
        handoff.full_signed_young_deficit + handoff.normalized_symbol_freezing_error
    )
    if full_margin + 5.0e-12 < mu:
        raise ValueError("handoff does not have the claimed full-signed Christ safety margin")
    if handoff.young_handoff_certified:
        return 0.0
    # With christ_modulus<1, a failed reserved gate cannot fail solely through
    # reservation_value<=0 while its deficit inequality passes: that value would
    # make delta_res>=1.  Hence failure plus the full margin forces b/Y>mu.
    if not (
        handoff.reserved_young_deficit + handoff.normalized_symbol_freezing_error
        > handoff.christ_modulus
    ):
        raise AssertionError("failed reserved gate did not fail its physical deficit inequality")
    if not handoff.terminal_bad_assistance_fraction > mu - 5.0e-12:
        raise AssertionError("reserved failure did not consume the claimed Christ safety margin")
    upper = handoff.inherited_negative_work + (1.0 / mu - 1.0) * handoff.inherited_bad_positive_work
    if handoff.inherited_good_positive_work > upper + 2.0e-9 * max(
        handoff.inherited_good_positive_work, upper, 1.0e-300
    ):
        raise AssertionError("reserved-failure good work escaped bad/backscatter domination")
    return upper


def reserved_failure_good_fraction_upper(
    block_transfer_deficit: float,
    *,
    full_signed_christ_margin_floor: float,
) -> float:
    """Universal fraction of good work that can fail reservation at margin ``mu``.

    Summing the cell inequality gives

      G_fail <= c_mu (B+N),  c_mu=max(1,1/mu-1),

    then the low-native-deficit actual-work theorem converts this to a fraction of
    the same canonical good ``dW+`` law.  ``mu`` is not an analyst fate threshold;
    it is the spare modulus in the existing full-signed Young/Christ premise.
    """
    mu = float(full_signed_christ_margin_floor)
    if not math.isfinite(mu) or not (0.0 < mu < 1.0):
        raise ValueError("Christ safety margin floor must lie in (0,1)")
    ratio = low_deficit_contamination_ratio_upper(block_transfer_deficit, ETA0)
    if not math.isfinite(ratio):
        return 1.0
    coefficient = max(1.0, 1.0 / mu - 1.0)
    return min(1.0, coefficient * ratio)


def arb_clean_constant_certificate() -> dict[str, str]:
    """Rigorous clean constants for the RN lower and the integer actual-NS triad."""
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover - exercised in Actions
        raise RuntimeError("python-flint is required for the mixed-fate rigorous certificate") from exc

    ctx.prec = 180

    def aq(q: Fraction):
        return arb(f"{q.numerator}/{q.denominator}")

    rlo = aq(RSTAR_LO)
    rhi = aq(RSTAR_HI)
    rstar = rlo.union(rhi)
    gamma = -rstar.log()
    sqrt2 = arb(2).sqrt()
    jstar = ((4 * rstar * rstar - 1).sqrt() * gamma) / (4 * sqrt2 * rstar)
    eta = aq(Fraction(1, 10_000))
    gap = aq(Fraction(1, 80))
    density_lower = jstar * (1 - eta) / (gamma + gap)
    clean = aq(CERTIFIED_GOOD_WORK_CAPACITY_DENSITY_LOWER)
    if not density_lower > clean:
        raise AssertionError(f"good-work/capacity clean lower 19/100 failed: {density_lower}")

    r0 = aq(ACTUAL_NS_GOOD_RATIO_SQUARED).sqrt()
    j0 = ((4 * r0 * r0 - 1).sqrt() * (-r0.log())) / (4 * sqrt2 * r0)
    normalized = j0 / jstar
    good_threshold = 1 - eta
    if not normalized > good_threshold:
        raise AssertionError(f"integer actual-NS triad missed eta0 good geometry: {normalized}")

    return {
        "rstar": str(rstar),
        "jstar": str(jstar),
        "gamma_star": str(gamma),
        "exact_good_density_lower_interval": str(density_lower),
        "certified_clean_good_density_lower": "19/100",
        "actual_ns_parent_child_ratio": str(r0),
        "actual_ns_normalized_geometric_multiplier": str(normalized),
        "canonical_good_threshold": "9999/10000",
    }


def phase_fate_role_refinement_counterexample(eta: float = ETA0) -> dict[str, float | bool]:
    """Same Fourier/helicity geometry, two positive phases, opposite G/B fate."""
    e = float(eta)
    if not (0.0 < e < 0.25):
        raise ValueError("small positive fate threshold required")
    r0 = math.sqrt(float(ACTUAL_NS_GOOD_RATIO_SQUARED))
    rstar = 0.5 * (float(RSTAR_LO) + float(RSTAR_HI))
    j0 = math.sqrt(4.0 * r0 * r0 - 1.0) * math.log(1.0 / r0) / (4.0 * math.sqrt(2.0) * r0)
    js = math.sqrt(4.0 * rstar * rstar - 1.0) * math.log(1.0 / rstar) / (4.0 * math.sqrt(2.0) * rstar)
    multiplier = j0 / js
    good_phase = 1.0
    bad_target_efficiency = 1.0 - 1.5 * e
    bad_phase = bad_target_efficiency / multiplier
    if not (0.0 < bad_phase < 1.0):
        raise AssertionError("phase counterexample left the positive-work sector")
    good_eff = multiplier * good_phase
    bad_eff = multiplier * bad_phase
    return {
        "eta": e,
        "same_normalized_geometric_multiplier": multiplier,
        "good_phase": good_phase,
        "bad_phase": bad_phase,
        "good_signed_efficiency": good_eff,
        "bad_signed_efficiency": bad_eff,
        "good_is_good": good_eff > 1.0 - e,
        "bad_is_bad": 0.0 < bad_eff <= 1.0 - e,
        "mode_role_geometry_changed": False,
    }


def positive_subtraction_not_trilinear_counterexample() -> dict[str, float | bool]:
    """Hahn-positive subtraction is nonlinear, so T-b may never be called work."""
    t1 = 1.0
    t2 = -1.0
    b1 = max(t1, 0.0)
    b2 = max(t2, 0.0)
    tsum = t1 + t2
    bsum_state = max(tsum, 0.0)
    return {
        "t1": t1,
        "t2": t2,
        "b_of_t1": b1,
        "b_of_t2": b2,
        "b_of_sum": bsum_state,
        "sum_of_b": b1 + b2,
        "positive_part_is_additive": math.isclose(bsum_state, b1 + b2),
        "reservation_is_trilinear_work": False,
    }


def coherent_fresh_hahn_kernel_counterexample() -> dict[str, float | bool]:
    """A non-diagonal localization can create cell cross-work from zero edge work.

    This is a finite exact trilinear countermodel to the shortcut
    ``fresh coherent Hahn = positive kernel applied to canonical edge Hahn``.
    It does *not* rule out a separately proved physical kernel whose output is
    different from fresh coherent Hahn.
    """

    def P(x: tuple[float, float]) -> tuple[float, float]:
        m = 0.5 * (x[0] + x[1])
        return (m, m)

    def tri(a: tuple[float, float], b: tuple[float, float], c: tuple[float, float]) -> float:
        return a[0] * b[0] * c[0] + a[1] * b[1] * c[1]

    u = (1.0, 0.0)
    v = (0.0, 1.0)
    w = (1.0, 1.0)
    edge_atoms = (u[0] * v[0] * w[0], u[1] * v[1] * w[1])
    localized = tri(P(u), P(v), P(w))
    canonical_positive = sum(max(x, 0.0) for x in edge_atoms)
    return {
        "canonical_edge_atom_0": edge_atoms[0],
        "canonical_edge_atom_1": edge_atoms[1],
        "canonical_positive_mass": canonical_positive,
        "non_diagonal_localized_signed_atom": localized,
        "fresh_localized_hahn_positive": max(localized, 0.0),
        "fresh_hahn_can_be_positive_kernel_pushforward_of_zero_canonical_mass": False,
    }


def theorem_certificate(*, rigorous: bool = False) -> dict[str, object]:
    phase = phase_fate_role_refinement_counterexample()
    nonlinear = positive_subtraction_not_trilinear_counterexample()
    coherent = coherent_fresh_hahn_kernel_counterexample()
    out: dict[str, object] = {
        "status": STATUS,
        "cell_hahn_identity": "T_C=g_C+b_C-n_C with g=pi_#dW_G+, b=pi_#dW_B+, n=pi_#dW- from the single edge-space Hahn decomposition",
        "reservation_identity": "delta_res=1-(T_C-b_C)/Y_C = (1-T_C/Y_C)+b_C/Y_C = delta_full+bad_assistance",
        "young_semantics": "Young/Christ receives only the actual full signed trilinear T_C; T_C-b_C is a sufficient-condition certificate only and is never a work law or causal law",
        "safe_handoff": "delta_res+xi<=delta_Christ implies delta_full+xi<=delta_Christ and certifies the inherited good dW+ branch without terminal bad-positive assistance",
        "fate_pure_recovery": "b_C=0 gives delta_res=delta_full, exactly recovering the previous fate-pure handoff",
        "low_deficit_actual_work_bound": "(B+N)/G <= q/[(19/100)(1-q)], q=epsilon/eta0; capacity is reference only and 19/100 is an Arb-certified pointwise physical RN lower",
        "reserved_failure_mass_bound": "for full-signed Christ safety margin mu in (0,1), failed reservation gives g<n+(1/mu-1)b; therefore G_fail/G <= max(1,1/mu-1) q/[(19/100)(1-q)]",
        "mode_role_no_go": phase,
        "nonlinear_subtraction_no_go": nonlinear,
        "coherent_fresh_hahn_no_go": coherent,
        "general_coherent_kernel_scope": "the countermodel rules out identifying fresh non-diagonal coherent Hahn atoms with a kernel pushforward; it does not rule out a different separately proved physical positive kernel",
        "claims_global_regularity": False,
    }
    if rigorous:
        out["arb_clean_constants"] = arb_clean_constant_certificate()
    return out


@dataclass(frozen=True)
class ReservedYoungStress:
    samples: int
    safe_handoff_cases: int
    mixed_fate_cases: int
    reserved_failure_domination_cases: int
    worst_reservation_identity_relative: float
    worst_hahn_identity_relative: float
    worst_positive_scaling_deficit_residual: float
    minimum_reserved_implies_full_margin: float
    phase_role_counterexample_passed: bool
    nonlinear_subtraction_counterexample_passed: bool
    coherent_fresh_hahn_counterexample_passed: bool


def _make_cell(g: float, b: float, n: float, label: str = "stress") -> HardCellWork:
    from src.canonical_positive_edge_work_routing import HardProductCell

    T = float(g) + float(b) - float(n)
    P = float(g) + float(b)
    fresh = max(T, 0.0)
    return HardCellWork(
        cell=HardProductCell(parent_roles=(f"{label}-p1", f"{label}-p2"), child_role=f"{label}-c"),
        signed_work=T,
        inherited_positive_work=P,
        inherited_good_positive_work=float(g),
        inherited_bad_positive_work=float(b),
        fresh_cell_hahn_positive=fresh,
        cancellation_gap=P - fresh,
    )


def stress(samples: int = 50_000, seed: int = 2026081202) -> ReservedYoungStress:
    rng = np.random.default_rng(seed)
    worst_res = 0.0
    worst_hahn = 0.0
    worst_scale = 0.0
    min_implication = math.inf
    safe = 0
    mixed = 0
    dominated = 0
    for i in range(int(samples)):
        scale0 = math.exp(float(rng.uniform(-12.0, 12.0)))
        if i % 2 == 0:
            Y = scale0
            full_deficit = float(rng.uniform(0.0, 0.025))
            bad_fraction = float(rng.uniform(0.0, 0.025))
            T = (1.0 - full_deficit) * Y
            b = bad_fraction * Y
            n = float(rng.uniform(0.0, 0.5)) * Y
            g = T - b + n
            if g <= 0.0:
                g = 0.05 * Y
                n = g + b - T
                if n < 0.0:
                    n = 0.0
                    T = g + b
                    Y = max(Y, 1.01 * T)
        else:
            g = scale0 * math.exp(float(rng.uniform(-2.0, 2.0)))
            b = scale0 * float(rng.uniform(0.0, 1.5))
            n = scale0 * float(rng.uniform(0.0, 2.5))
            T = g + b - n
            Y = abs(T) + scale0 * float(rng.uniform(0.02, 2.0))
        if b > 0.0:
            mixed += 1
        cell = _make_cell(g, b, n, label=f"s{i}")
        xi = float(rng.uniform(0.0, 0.01))
        modulus = float(rng.uniform(0.015, 0.08))
        cert = certify_reserved_young_handoff(
            cell,
            certified_full_cell_young_upper=Y,
            normalized_symbol_freezing_error=xi,
            christ_modulus=modulus,
        )
        safe += int(cert.young_handoff_certified)
        full_margin = modulus - (cert.full_signed_young_deficit + xi)
        if not cert.young_handoff_certified and full_margin > 2.0e-6:
            mu = 0.5 * full_margin
            upper = reserved_failure_good_work_upper(cert, full_signed_christ_margin_floor=mu)
            if cert.inherited_good_positive_work > upper + 3.0e-9 * max(cert.inherited_good_positive_work, upper, 1.0e-300):
                raise AssertionError("reserved-failure margin domination regressed")
            dominated += 1
        rscale = max(Y, abs(cert.reservation_value), abs(cert.signed_full_cell_work), 1.0e-300)
        worst_res = max(
            worst_res,
            abs(
                cert.reserved_young_deficit
                - cert.full_signed_young_deficit
                - cert.terminal_bad_assistance_fraction
            ) / max(1.0, abs(cert.reserved_young_deficit)),
        )
        worst_hahn = max(
            worst_hahn,
            abs(
                cert.signed_full_cell_work
                - cert.inherited_good_positive_work
                - cert.inherited_bad_positive_work
                + cert.inherited_negative_work
            ) / rscale,
        )
        lam = math.exp(float(rng.uniform(-8.0, 8.0)))
        scaled = _make_cell(lam * g, lam * b, lam * n, label=f"q{i}")
        scert = certify_reserved_young_handoff(
            scaled,
            certified_full_cell_young_upper=lam * Y,
            normalized_symbol_freezing_error=xi,
            christ_modulus=modulus,
        )
        worst_scale = max(
            worst_scale,
            abs(scert.full_signed_young_deficit - cert.full_signed_young_deficit),
            abs(scert.terminal_bad_assistance_fraction - cert.terminal_bad_assistance_fraction),
            abs(scert.reserved_young_deficit - cert.reserved_young_deficit),
        )
        if cert.reserved_christ_gate:
            margin = modulus - (cert.full_signed_young_deficit + xi)
            min_implication = min(min_implication, margin)
            if margin < -2.0e-12:
                raise AssertionError("reserved mixed-fate gate passed while full signed Christ gate failed")

    phase = phase_fate_role_refinement_counterexample()
    nonlinear = positive_subtraction_not_trilinear_counterexample()
    coherent = coherent_fresh_hahn_kernel_counterexample()
    phase_ok = bool(phase["good_is_good"] and phase["bad_is_bad"] and not phase["mode_role_geometry_changed"])
    nonlinear_ok = bool(not nonlinear["positive_part_is_additive"] and not nonlinear["reservation_is_trilinear_work"])
    coherent_ok = bool(
        coherent["canonical_positive_mass"] == 0.0
        and coherent["fresh_localized_hahn_positive"] > 0.0
        and not coherent["fresh_hahn_can_be_positive_kernel_pushforward_of_zero_canonical_mass"]
    )
    if not (phase_ok and nonlinear_ok and coherent_ok):
        raise AssertionError("one mixed-fate anti-shortcut counterexample regressed")
    if safe == 0 or mixed == 0:
        raise AssertionError("mixed-fate stress failed to exercise safe and mixed branches")
    if dominated == 0:
        raise AssertionError("mixed-fate stress failed to exercise Christ-margin failure domination")
    return ReservedYoungStress(
        samples=int(samples),
        reserved_failure_domination_cases=dominated,
        safe_handoff_cases=safe,
        mixed_fate_cases=mixed,
        worst_reservation_identity_relative=worst_res,
        worst_hahn_identity_relative=worst_hahn,
        worst_positive_scaling_deficit_residual=worst_scale,
        minimum_reserved_implies_full_margin=0.0 if math.isinf(min_implication) else min_implication,
        phase_role_counterexample_passed=phase_ok,
        nonlinear_subtraction_counterexample_passed=nonlinear_ok,
        coherent_fresh_hahn_counterexample_passed=coherent_ok,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=50_000)
    parser.add_argument("--seed", type=int, default=2026081202)
    parser.add_argument("--outdir", type=Path, default=Path("results-mixed-fate-reserved-young-handoff"))
    args = parser.parse_args()
    out = stress(args.samples, args.seed)
    cert = theorem_certificate(rigorous=True)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "certificate.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2, sort_keys=True) + "\n"
    )
    summary = f"""# Mixed-fate reserved Young/Christ handoff\n\nStatus: **{STATUS}**.\n\nThe canonical edge law is unchanged.  On one hard product cell, write\n\n`T_C = g_C + b_C - n_C`,\n\nwhere `g_C` is inherited geometry-good `dW+`, `b_C` is inherited already-terminal geometry-bad `dW+`, and `n_C` is the pushforward of the original edge Hahn-negative law.  For any certified full-cell Young upper `Y_C`,\n\n`delta_res = 1-(T_C-b_C)/Y_C = (1-T_C/Y_C) + b_C/Y_C`.\n\n`T_C-b_C` is **not** a new trilinear work and is never sent to Young.  It is only a sufficient-condition certificate: if `delta_res + xi <= delta_Christ`, then the actual full signed `T_C` satisfies the ordinary complex Young/Christ gate without needing terminal bad-positive assistance.\n\nStress: `{out.samples}` algebra/scaling states\n- safe mixed/fate-pure handoffs exercised: `{out.safe_handoff_cases}`\n- mixed-fate states exercised: `{out.mixed_fate_cases}`\n- worst reservation identity relative residual: `{out.worst_reservation_identity_relative:.3e}`\n- worst Hahn pushforward identity relative residual: `{out.worst_hahn_identity_relative:.3e}`\n- worst positive scaling deficit residual: `{out.worst_positive_scaling_deficit_residual:.3e}`\n- phase-only hard-role fate counterexample: `{out.phase_role_counterexample_passed}`\n- nonlinear Hahn-subtraction counterexample: `{out.nonlinear_subtraction_counterexample_passed}`\n- fresh coherent-Hahn kernel shortcut counterexample: `{out.coherent_fresh_hahn_counterexample_passed}`\n\nNo global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
