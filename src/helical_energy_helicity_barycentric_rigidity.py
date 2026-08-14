from __future__ import annotations

import math
from dataclasses import dataclass

from src.cyclic_helical_triad_donor_kernel import ClosedHelicalTriadRegistration
from src.helical import stable_norm3

STATUS = (
    "DRAFT_HELICAL_ENERGY_HELICITY_BARYCENTRIC_RIGIDITY__"
    "T_EQUALS_R_LAMBDA_CROSS_ONE__ENERGY_AND_HELICITY_CONSERVED__"
    "MEDIAN_CURL_MODE_SINGLETON__CONVEX_SPREAD_OR_CONTRACTION__"
    "QUADRATIC_MOMENT_IS_ENSTROPHY_PRODUCTION"
)


def curl_eigenvalue(mode) -> float:
    return float(mode.helicity) * stable_norm3(mode.wavevector)


def _native_ratio(value: float, scale: float) -> float:
    if scale <= 0.0:
        return 0.0 if value == 0.0 else math.inf
    return abs(float(value)) / float(scale)


@dataclass(frozen=True)
class HelicalEnergyHelicityRigidity:
    lambdas_by_slot: tuple[float, float, float]
    works_by_slot: tuple[float, float, float]
    energy_conservation_native_residual: float
    helicity_conservation_native_residual: float
    ordered_slot_indices: tuple[int, int, int]
    ordered_lambdas: tuple[float, float, float]
    ordered_works: tuple[float, float, float]
    lambda_gaps: tuple[float, float]
    transfer_orientation: str
    quadratic_moment_production: float
    quadratic_identity_native_residual: float
    median_is_unique_singleton_side: bool
    strict_uv_frontier_positive_slots: tuple[int, ...]
    strict_uv_frontier_slots_are_spreads: bool

def certify_helical_energy_helicity_rigidity(
    triad: ClosedHelicalTriadRegistration,
) -> HelicalEnergyHelicityRigidity:
    lambdas = tuple(curl_eigenvalue(slot.closed_mode) for slot in triad.slots)
    works = tuple(float(slot.signed_work) for slot in triad.slots)
    capacities = tuple(float(slot.edge_registration.native_modal_capacity) for slot in triad.slots)

    work_scale = math.fsum(capacities)
    helicity_scale = math.fsum(abs(lam) * cap for lam, cap in zip(lambdas, capacities))
    quadratic_scale = math.fsum(lam * lam * cap for lam, cap in zip(lambdas, capacities))

    energy_res = _native_ratio(math.fsum(works), work_scale)
    helicity_work = math.fsum(lam * work for lam, work in zip(lambdas, works))
    helicity_res = _native_ratio(helicity_work, helicity_scale)
    if energy_res > 5.0e-10:
        raise AssertionError("closed helical triad lost nonlinear energy conservation")
    if helicity_res > 8.0e-10:
        raise AssertionError("closed helical triad lost nonlinear helicity conservation")

    order = tuple(sorted(range(3), key=lambda i: (lambdas[i], i)))
    lo, mid, hi = order
    ordered_lambdas = (lambdas[lo], lambdas[mid], lambdas[hi])
    ordered_works = (works[lo], works[mid], works[hi])
    gap_lo = ordered_lambdas[1] - ordered_lambdas[0]
    gap_hi = ordered_lambdas[2] - ordered_lambdas[1]
    lam_scale = max(abs(v) for v in ordered_lambdas)
    gap_tol = 8.0e-12 * max(lam_scale, 1.0)
    work_tol = 8.0e-10 * max(work_scale, 1.0e-300)

    if gap_lo <= gap_tol or gap_hi <= gap_tol:
        orientation = "curl_eigenvalue_degenerate"
        singleton = False
    else:
        tlo, tm, thi = ordered_works
        if abs(tm) <= work_tol:
            orientation = "numerically_unresolved_or_zero"
            singleton = False
        elif tm < 0.0:
            orientation = "mean_preserving_spread"
            singleton = tlo > -work_tol and thi > -work_tol
        else:
            orientation = "mean_preserving_contraction"
            singleton = tlo < work_tol and thi < work_tol

        expected_lo = -(gap_hi / (gap_lo + gap_hi)) * tm
        expected_hi = -(gap_lo / (gap_lo + gap_hi)) * tm
        bary_scale = max(work_scale, abs(tm), abs(expected_lo), abs(expected_hi), 1.0e-300)
        if abs(tlo - expected_lo) > 1.2e-9 * bary_scale:
            raise AssertionError("low curl-eigenvalue work left the energy-helicity barycentric law")
        if abs(thi - expected_hi) > 1.2e-9 * bary_scale:
            raise AssertionError("high curl-eigenvalue work left the energy-helicity barycentric law")
        if not singleton and orientation in {"mean_preserving_spread", "mean_preserving_contraction"}:
            raise AssertionError("median curl-eigenvalue mode is not the unique singleton work side")

    quadratic = math.fsum(lam * lam * work for lam, work in zip(lambdas, works))
    expected_quadratic = -ordered_works[1] * gap_lo * gap_hi
    quadratic_res = _native_ratio(quadratic - expected_quadratic, quadratic_scale)
    if quadratic_res > 1.2e-9:
        raise AssertionError("quadratic curl-spectral moment left the exact enstrophy identity")

    frontier: list[int] = []
    radii = tuple(abs(lam) for lam in lambdas)
    radial_tol = 8.0e-12 * max(max(radii), 1.0)
    for slot in triad.slots:
        i = slot.closed_mode_index
        j, ell = slot.parent_closed_indices
        if works[i] > work_tol and radii[i] > max(radii[j], radii[ell]) + radial_tol:
            frontier.append(i)
    frontier_spread = not frontier or orientation == "mean_preserving_spread"
    if not frontier_spread:
        raise AssertionError("strict UV-frontier positive child was not a barycentric spread")

    return HelicalEnergyHelicityRigidity(
        lambdas_by_slot=lambdas,
        works_by_slot=works,
        energy_conservation_native_residual=energy_res,
        helicity_conservation_native_residual=helicity_res,
        ordered_slot_indices=order,
        ordered_lambdas=ordered_lambdas,
        ordered_works=ordered_works,
        lambda_gaps=(gap_lo, gap_hi),
        transfer_orientation=orientation,
        quadratic_moment_production=quadratic,
        quadratic_identity_native_residual=quadratic_res,
        median_is_unique_singleton_side=singleton,
        strict_uv_frontier_positive_slots=tuple(frontier),
        strict_uv_frontier_slots_are_spreads=frontier_spread,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "curl_coordinate": "lambda=s|k| is the physical curl eigenvalue of each helical Fourier mode",
        "normal_form": "T_triangle=R_triangle(lambda cross 1) on the already-certified closed-triad rooted works",
        "two_invariants": "sum_i T_i=0 and sum_i lambda_i T_i=0: nonlinear energy and helicity conservation",
        "dimension": "unless all curl eigenvalues coincide, the three-work vector lies on the unique one-dimensional intersection of the energy/helicity conservation planes",
        "median_law": "for distinct ordered lambda_-<lambda_0<lambda_+, the median slot has the opposite work sign to both extremes and is the unique singleton donor/recipient side",
        "barycenter": "a median donor splits energy to the two extremes with the unique weights making lambda_0 their barycenter; reverse orientation is contraction",
        "convex_order": "every convex curl-spectral moment increases on a spread and decreases on the reverse contraction",
        "enstrophy": "sum lambda_i^2 T_i=-T_median(lambda_median-lambda_-)(lambda_+-lambda_median), the exact triad contribution to nonlinear enstrophy/vortex-stretching production",
        "uv_frontier": "a positive child strictly beyond both interaction-parent radii must be an extreme recipient, hence a spread; a contraction recipient always has a same-event donor at radius at least as large",
        "signed_good": "the existing signed-good unique-donor/side-recipient work ratios are barycentric corollaries after opposite parent helicities are certified",
        "critical_pair": "for K_+=sum_{s=+}|k|E and K_-=sum_{s=-}|k|E, every closed triad has dK_+^NL=dK_-^NL; homochiral triads give zero, while a heterochiral triad source is the |k|-weighted work of its singleton-helicity mode",
        "critical_side": "positive H1/2 production is a heterochiral spread and carries a compulsory positive nonforward side recipient on the same triad",
        "good_edge_critical": "on canonical good positive work, signed-good geometry gives 18/49 < (d/dt ||u||_{Hdot1/2}^2)_triad / (M W_child+) < 20/49; this is a state-balance source, not a causal reweighting",
        "causal_scope": "this does not replace canonical dW+ or create a new probability/currency; it identifies invariant structure and native critical-state production of the same signed physical work",
        "claims_mixed_recurrence_closed": False,
        "claims_global_regularity": False,
    }

@dataclass(frozen=True)
class HelicalCriticalPairBalance:
    positive_helicity_critical_rate: float
    negative_helicity_critical_rate: float
    absolute_critical_rate: float
    pair_source_native_residual: float
    homochiral: bool
    singleton_helicity_closed_mode_index: int | None
    singleton_helicity_sign: int | None
    singleton_helicity_weighted_work: float
    singleton_identity_native_residual: float
    critical_growth_has_nonforward_side: bool
    side_recipient_closed_mode_index: int | None
    side_weighted_positive_work: float
    side_tax_native_margin: float

    def __post_init__(self) -> None:
        vals = (
            self.positive_helicity_critical_rate,
            self.negative_helicity_critical_rate,
            self.absolute_critical_rate,
            self.pair_source_native_residual,
            self.singleton_helicity_weighted_work,
            self.singleton_identity_native_residual,
            self.side_weighted_positive_work,
            self.side_tax_native_margin,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite critical helicity-pair balance required")
        if self.pair_source_native_residual > 1.2e-9 or self.singleton_identity_native_residual > 1.2e-9:
            raise AssertionError("critical helicity-pair law left its native weighted-work scale")
        if self.homochiral:
            if self.singleton_helicity_closed_mode_index is not None or self.singleton_helicity_sign is not None:
                raise ValueError("homochiral triad may not invent a singleton helicity mode")
        else:
            if self.singleton_helicity_closed_mode_index not in (0,1,2) or self.singleton_helicity_sign not in (-1,1):
                raise ValueError("heterochiral triad requires its unique singleton-helicity mode")
        if self.absolute_critical_rate > 0.0 and not self.homochiral:
            if not self.critical_growth_has_nonforward_side or self.side_recipient_closed_mode_index not in (0,1,2):
                raise AssertionError("positive critical growth lost its compulsory same-triad nonforward side")
            if self.side_tax_native_margin < -1.2e-9:
                raise AssertionError("critical growth exceeded the weighted nonforward side tax")


def critical_helicity_pair_balance(
    triad: ClosedHelicalTriadRegistration,
) -> HelicalCriticalPairBalance:
    """Exact nonlinear balance of the two positive helicity magnitudes.

    K_+=sum_{s=+}|k|E and K_-=sum_{s=-}|k|E.  Nonlinear helicity
    conservation forces dK_+^NL=dK_-^NL on every closed triad.  For a
    heterochiral triad one helicity sign occurs on exactly one mode, so the
    common source is just that singleton mode's |k|-weighted rooted work.
    """
    slots = tuple(triad.slots)
    lambdas = tuple(curl_eigenvalue(slot.closed_mode) for slot in slots)
    works = tuple(float(slot.signed_work) for slot in slots)
    capacities = tuple(float(slot.edge_registration.native_modal_capacity) for slot in slots)
    weighted_scale = math.fsum(abs(lam) * cap for lam, cap in zip(lambdas, capacities))
    weighted_scale = max(weighted_scale, 1.0e-300)

    kp = math.fsum(max(lam, 0.0) * work for lam, work in zip(lambdas, works))
    km = math.fsum(max(-lam, 0.0) * work for lam, work in zip(lambdas, works))
    kab = kp + km
    pair_res = abs(kp-km) / weighted_scale
    if pair_res > 1.2e-9:
        raise AssertionError("nonlinear helicity conservation did not pair the two critical helicity sectors")

    signs = tuple(slot.closed_mode.helicity for slot in slots)
    homochiral = len(set(signs)) == 1
    singleton_index: int | None = None
    singleton_sign: int | None = None
    singleton_weighted = 0.0
    singleton_res = 0.0
    side_index: int | None = None
    side_weighted = 0.0
    side_margin = 0.0
    has_side = False

    if homochiral:
        if abs(kab) > 1.2e-9 * weighted_scale:
            raise AssertionError("homochiral triad changed the absolute-helicity/H1/2 critical moment")
    else:
        counts = {s: signs.count(s) for s in (-1,1)}
        singleton_sign = -1 if counts[-1] == 1 else 1
        singleton_index = signs.index(singleton_sign)
        singleton_weighted = abs(lambdas[singleton_index]) * works[singleton_index]
        singleton_res = abs(0.5*kab-singleton_weighted) / weighted_scale
        if singleton_res > 1.2e-9:
            raise AssertionError("critical source is not the singleton-helicity rooted work")

        if kab > 1.2e-12 * weighted_scale:
            positive_indices = [i for i,w in enumerate(works) if w > 8.0e-10*max(math.fsum(capacities),1.0e-300)]
            if len(positive_indices) != 2:
                raise AssertionError("critical-growing spread did not have exactly two positive recipients")
            side_index = min(positive_indices, key=lambda i: abs(lambdas[i]))
            side_slot = slots[side_index]
            side_weighted = abs(lambdas[side_index]) * works[side_index]
            if not (side_slot.edge_registration.scale_progress == 0.0 and side_slot.edge_registration.geometric_multiplier_J == 0.0):
                raise AssertionError("smaller-radial spread recipient is not the exact nonforward side edge")
            side_margin = (side_weighted - 0.5*kab) / weighted_scale
            has_side = True
            if side_margin < -1.2e-9:
                raise AssertionError("nonforward side does not dominate half the positive critical production")

    return HelicalCriticalPairBalance(
        positive_helicity_critical_rate=kp,
        negative_helicity_critical_rate=km,
        absolute_critical_rate=kab,
        pair_source_native_residual=pair_res,
        homochiral=homochiral,
        singleton_helicity_closed_mode_index=singleton_index,
        singleton_helicity_sign=singleton_sign,
        singleton_helicity_weighted_work=singleton_weighted,
        singleton_identity_native_residual=singleton_res,
        critical_growth_has_nonforward_side=has_side,
        side_recipient_closed_mode_index=side_index,
        side_weighted_positive_work=side_weighted,
        side_tax_native_margin=side_margin,
    )

CRITICAL_GOOD_EDGE_RATIO_LO = 18.0 / 49.0
CRITICAL_GOOD_EDGE_RATIO_HI = 20.0 / 49.0


@dataclass(frozen=True)
class GoodEdgeCriticalSourceBridge:
    recipient_closed_mode_index: int
    child_frequency: float
    recipient_positive_work: float
    critical_source_rate: float
    critical_source_to_child_scale_work_ratio: float
    ratio_lower: float = CRITICAL_GOOD_EDGE_RATIO_LO
    ratio_upper: float = CRITICAL_GOOD_EDGE_RATIO_HI
    canonical_causality_reweighted: bool = False
    critical_source_declared_causal_probability: bool = False

    def __post_init__(self) -> None:
        if self.recipient_closed_mode_index not in (0,1,2):
            raise ValueError("good-edge critical bridge requires one closed-triad recipient root")
        if min(self.child_frequency, self.recipient_positive_work, self.critical_source_rate) <= 0.0:
            raise ValueError("good-edge critical bridge requires positive physical scale/work/source")
        if not (self.ratio_lower < self.critical_source_to_child_scale_work_ratio < self.ratio_upper):
            raise AssertionError("good-edge critical source left the exact 18/49..20/49 corridor")
        if self.canonical_causality_reweighted or self.critical_source_declared_causal_probability:
            raise ValueError("critical H1/2 source is a state-balance law, not a replacement causal measure")


def good_edge_critical_source_bridge(
    triad: ClosedHelicalTriadRegistration,
    *,
    recipient_closed_mode_index: int,
) -> GoodEdgeCriticalSourceBridge:
    """Canonical good positive edge -> native critical H1/2 production.

    A good edge has positive work and signed efficiency >1-eta0.  The existing
    cyclic theorem then makes it a signed-good forward main recipient with
    parent/child ratios D,S in (3/5,5/8), one median energy donor and one
    positive nonforward side recipient.  The new energy-helicity law identifies
    the same triad as a heterochiral spread, so its nonlinear critical source is

        dK/dt = 2 S M W_side.

    Dividing by M W_child and using W_side/W_child=(1-D)/(D+S) gives the exact
    open interval 18/49 .. 20/49.
    """
    from src.cyclic_helical_triad_donor_kernel import signed_good_side_recipient_certificate

    cert = signed_good_side_recipient_certificate(
        triad, recipient_closed_mode_index=recipient_closed_mode_index
    )
    rigidity = certify_helical_energy_helicity_rigidity(triad)
    critical = critical_helicity_pair_balance(triad)
    if rigidity.transfer_orientation != "mean_preserving_spread":
        raise AssertionError("canonical good positive edge is not a barycentric spread")
    if recipient_closed_mode_index not in rigidity.strict_uv_frontier_positive_slots:
        raise AssertionError("canonical good positive edge is not the unique forward spread recipient")
    if critical.absolute_critical_rate <= 0.0:
        raise AssertionError("canonical good positive edge failed to create positive H1/2 critical source")
    if cert.side_recipient_closed_mode_index != critical.side_recipient_closed_mode_index:
        raise AssertionError("signed-good side recipient disagrees with critical side-tax law")
    slot = triad.slot_for_closed_mode_index(recipient_closed_mode_index)
    M = float(slot.edge_registration.child_frequency)
    W = float(slot.signed_work)
    ratio = critical.absolute_critical_rate / (M * W)
    return GoodEdgeCriticalSourceBridge(
        recipient_closed_mode_index=recipient_closed_mode_index,
        child_frequency=M,
        recipient_positive_work=W,
        critical_source_rate=critical.absolute_critical_rate,
        critical_source_to_child_scale_work_ratio=ratio,
    )
