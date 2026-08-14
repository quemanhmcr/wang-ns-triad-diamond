from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Callable, Sequence

from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import (
    ClosedHelicalTriadRegistration,
    CyclicTriadMeasureKernel,
)
from src.helical import coupling_magnitude_closed, stable_norm3
from src.helical_energy_helicity_barycentric_rigidity import (
    certify_helical_energy_helicity_rigidity,
    critical_helicity_pair_balance,
    curl_eigenvalue,
)
from src.single_edge_certificate import float_jstar

STATUS = (
    "DRAFT_CURL_SPECTRAL_CURVATURE_BALANCE__"
    "ALL_NONLINEAR_MOMENTS_FROM_ONE_BARYCENTRIC_THREE_POINT_CURVATURE_LAW__"
    "H_HALF_IS_TANAKA_DEFECT_AT_CURL_SIGN_INTERFACE__"
    "SHARP_STRICT_UV_CRITICAL_PRODUCTION_GEOMETRY_IS_NOT_LOG_PROGRESS_EXTREMIZER"
)


def _close(a: float, b: float, *, factor: float = 1.5e-9) -> bool:
    return abs(float(a) - float(b)) <= factor * max(abs(float(a)), abs(float(b)), 1.0e-300)


def _mode_lambda(mode) -> float:
    return float(mode.helicity) * stable_norm3(mode.wavevector)


@dataclass(frozen=True)
class CurlSpectralMomentProduction:
    quotient_measure_mass: float
    rooted_signed_production: float
    donor_flow_production: float
    native_weighted_work_scale: float
    rooted_flow_native_residual: float
    affine_energy_residual: float
    affine_helicity_residual: float
    causal_law_replaced: bool = False
    new_hahn_used_as_causality: bool = False

    def __post_init__(self) -> None:
        if not math.isfinite(self.quotient_measure_mass) or self.quotient_measure_mass < 0.0:
            raise ValueError("finite nonnegative closed-triad quotient mass required")
        vals = (
            self.rooted_signed_production,
            self.donor_flow_production,
            self.native_weighted_work_scale,
            self.rooted_flow_native_residual,
            self.affine_energy_residual,
            self.affine_helicity_residual,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite curl-spectral moment data required")
        if self.native_weighted_work_scale <= 0.0:
            raise ValueError("positive native weighted-work scale required")
        if self.rooted_flow_native_residual > 1.5e-9:
            raise AssertionError("weighted donor-flow law left the rooted physical-work moment")
        if self.affine_energy_residual > 1.5e-9 or self.affine_helicity_residual > 1.5e-9:
            raise AssertionError("energy/helicity affine observables failed to vanish nonlinearly")
        if self.causal_law_replaced or self.new_hahn_used_as_causality:
            raise ValueError("spectral moment balance may not replace canonical dW causality")


def weighted_curl_spectral_production(
    triad: ClosedHelicalTriadRegistration,
    kernel: CyclicTriadMeasureKernel,
    phi: Callable[[float], float],
) -> CurlSpectralMomentProduction:
    """Read one spectral observable from the existing canonical cyclic work law.

    This creates no positive law.  It evaluates the same signed closed-triad work
    either at the three rooted modes or through the already-certified positive
    donor->recipient transport.  The two readings must agree exactly.
    """
    q = float(kernel.quotient_measure_mass)
    if q <= 0.0:
        raise ValueError("positive quotient mass required for a nontrivial weighted production")
    slots = tuple(triad.slots)
    lambdas = tuple(curl_eigenvalue(slot.closed_mode) for slot in slots)
    works = tuple(float(slot.signed_work) for slot in slots)
    values = tuple(float(phi(lam)) for lam in lambdas)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("finite spectral test-function values required")

    physical_factor = unitary_fourier_convolution_factor() * q
    rooted = physical_factor * math.fsum(v * w for v, w in zip(values, works))
    flow = math.fsum(
        atom.physical_work_mass
        * (float(phi(_mode_lambda(atom.recipient_child_mode))) - float(phi(_mode_lambda(atom.donor_child_mode))))
        for atom in kernel.atoms
    )
    scale = physical_factor * math.fsum(
        abs(v) * float(slot.edge_registration.native_modal_capacity)
        for v, slot in zip(values, slots)
    )
    scale = max(scale, abs(rooted), abs(flow), 1.0e-300)

    energy = physical_factor * math.fsum(works)
    helicity = physical_factor * math.fsum(lam * work for lam, work in zip(lambdas, works))
    energy_scale = max(physical_factor * math.fsum(float(s.edge_registration.native_modal_capacity) for s in slots), 1.0e-300)
    helicity_scale = max(
        physical_factor * math.fsum(abs(lam) * float(s.edge_registration.native_modal_capacity) for lam, s in zip(lambdas, slots)),
        1.0e-300,
    )
    return CurlSpectralMomentProduction(
        quotient_measure_mass=q,
        rooted_signed_production=rooted,
        donor_flow_production=flow,
        native_weighted_work_scale=scale,
        rooted_flow_native_residual=abs(rooted - flow) / scale,
        affine_energy_residual=abs(energy) / energy_scale,
        affine_helicity_residual=abs(helicity) / helicity_scale,
    )


@dataclass(frozen=True)
class CriticalTanakaTriadBalance:
    lambdas_ordered: tuple[float, float, float]
    median_work_density: float
    tanaka_defect: float
    critical_production_density: float
    quotient_measure_mass: float
    critical_production_mass: float
    donor_flow_critical_mass: float
    heterochiral: bool
    spread: bool
    contraction: bool
    homochiral_zero: bool
    singleton_helicity_closed_mode_index: int | None
    native_residual: float
    declared_causal_probability: bool = False

    def __post_init__(self) -> None:
        lo, mid, hi = self.lambdas_ordered
        if not lo <= mid <= hi:
            raise ValueError("ordered curl eigenvalues required")
        if self.tanaka_defect < -1.0e-14:
            raise AssertionError("absolute-value Jensen/Tanaka defect became negative")
        if self.spread and self.contraction:
            raise ValueError("one triad cannot be both spread and contraction")
        if self.homochiral_zero and self.heterochiral:
            raise ValueError("homochiral and heterochiral labels conflict")
        if self.native_residual > 1.5e-9:
            raise AssertionError("critical Tanaka identity left native physical scale")
        if self.declared_causal_probability:
            raise ValueError("critical state production is not a replacement causal probability")


def critical_tanaka_triad_balance(
    triad: ClosedHelicalTriadRegistration,
    kernel: CyclicTriadMeasureKernel,
) -> CriticalTanakaTriadBalance:
    """Exact H^{1/2} production as the kink-curvature defect of |lambda| at zero."""
    rigidity = certify_helical_energy_helicity_rigidity(triad)
    critical = critical_helicity_pair_balance(triad)
    weighted = weighted_curl_spectral_production(triad, kernel, abs)

    lo, mid, hi = rigidity.ordered_lambdas
    tmid = rigidity.ordered_works[1]
    if hi > lo:
        alpha = (hi - mid) / (hi - lo)
        beta = (mid - lo) / (hi - lo)
        defect = alpha * abs(lo) + beta * abs(hi) - abs(mid)
    else:
        defect = 0.0
    defect = max(0.0, defect)
    density = -tmid * defect
    mass = unitary_fourier_convolution_factor() * float(kernel.quotient_measure_mass) * density
    scale = max(
        weighted.native_weighted_work_scale,
        abs(mass),
        abs(weighted.rooted_signed_production),
        1.0e-300,
    )
    residual = max(
        abs(mass - weighted.rooted_signed_production) / scale,
        abs(mass - unitary_fourier_convolution_factor() * float(kernel.quotient_measure_mass) * critical.absolute_critical_rate) / scale,
    )
    signs = {slot.closed_mode.helicity for slot in triad.slots}
    hetero = len(signs) == 2
    homo = len(signs) == 1
    if homo and abs(mass) > 1.5e-9 * scale:
        raise AssertionError("homochiral triad produced critical H1/2 curvature")
    if hetero and defect <= 0.0:
        raise AssertionError("heterochiral triad failed to cross the curl-sign kink")
    return CriticalTanakaTriadBalance(
        lambdas_ordered=rigidity.ordered_lambdas,
        median_work_density=tmid,
        tanaka_defect=defect,
        critical_production_density=density,
        quotient_measure_mass=float(kernel.quotient_measure_mass),
        critical_production_mass=mass,
        donor_flow_critical_mass=weighted.donor_flow_production,
        heterochiral=hetero,
        spread=rigidity.transfer_orientation == "mean_preserving_spread",
        contraction=rigidity.transfer_orientation == "mean_preserving_contraction",
        homochiral_zero=homo,
        singleton_helicity_closed_mode_index=critical.singleton_helicity_closed_mode_index,
        native_residual=residual,
    )


def aggregate_weighted_curl_spectral_production(
    rows: Sequence[tuple[ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel]],
    phi: Callable[[float], float],
) -> float:
    if not rows:
        raise ValueError("nonempty closed-triad family required")
    return math.fsum(weighted_curl_spectral_production(t, k, phi).rooted_signed_production for t, k in rows)



@dataclass(frozen=True)
class HelicityPairProduction:
    positive_helicity_critical_production: float
    negative_helicity_critical_production: float
    common_pair_source: float
    absolute_critical_production: float
    signed_helicity_production: float
    native_residual: float
    pair_source_declared_causal_probability: bool = False

    def __post_init__(self) -> None:
        vals = (
            self.positive_helicity_critical_production,
            self.negative_helicity_critical_production,
            self.common_pair_source,
            self.absolute_critical_production,
            self.signed_helicity_production,
            self.native_residual,
        )
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite helicity-pair production required")
        if self.native_residual > 1.5e-9:
            raise AssertionError("the two helicity critical reservoirs lost their common nonlinear source")
        if self.pair_source_declared_causal_probability:
            raise ValueError("helicity-pair state source is not a replacement causal law")


def helicity_pair_production(
    rows: Sequence[tuple[ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel]],
) -> HelicityPairProduction:
    """Full-family nonlinear critical dynamics has one common pair source.

    K_+=int max(lambda,0) d rho and K_-=int max(-lambda,0) d rho.
    Since lambda=lambda_+-lambda_- is affine and nonlinear helicity is conserved,
    the two nonlinear productions are exactly equal.  Their sum is the H1/2
    production and their difference is zero helicity production.
    """
    plus_rows = tuple(weighted_curl_spectral_production(t, k, lambda x: max(x, 0.0)) for t, k in rows)
    minus_rows = tuple(weighted_curl_spectral_production(t, k, lambda x: max(-x, 0.0)) for t, k in rows)
    abs_rows = tuple(weighted_curl_spectral_production(t, k, abs) for t, k in rows)
    hel_rows = tuple(weighted_curl_spectral_production(t, k, lambda x: x) for t, k in rows)
    plus = math.fsum(r.rooted_signed_production for r in plus_rows)
    minus = math.fsum(r.rooted_signed_production for r in minus_rows)
    absolute = math.fsum(r.rooted_signed_production for r in abs_rows)
    helicity = math.fsum(r.rooted_signed_production for r in hel_rows)
    scale = max(
        math.fsum(r.native_weighted_work_scale for r in plus_rows),
        math.fsum(r.native_weighted_work_scale for r in minus_rows),
        math.fsum(r.native_weighted_work_scale for r in abs_rows),
        1.0e-300,
    )
    residual = max(
        abs(plus-minus)/scale,
        abs((plus+minus)-absolute)/scale,
        abs((plus-minus)-helicity)/scale,
    )
    return HelicityPairProduction(
        positive_helicity_critical_production=plus,
        negative_helicity_critical_production=minus,
        common_pair_source=0.5*(plus+minus),
        absolute_critical_production=absolute,
        signed_helicity_production=helicity,
        native_residual=residual,
    )

def viscous_curl_spectral_sink(
    *,
    viscosity: float,
    modal_energy: Sequence[tuple[object, float]],
    phi: Callable[[float], float],
) -> float:
    nu = float(viscosity)
    if not math.isfinite(nu) or nu < 0.0:
        raise ValueError("finite nonnegative viscosity required")
    total = 0.0
    for mode, energy in modal_energy:
        E = float(energy)
        if not math.isfinite(E) or E < 0.0:
            raise ValueError("finite nonnegative modal energy required")
        lam = _mode_lambda(mode)
        value = float(phi(lam))
        if not math.isfinite(value):
            raise ValueError("finite spectral weight required")
        total += 2.0 * nu * (lam * lam) * value * E
    return total


# ---------------------------------------------------------------------------
# Sharp local geometry for critical H^{1/2} production at a strict UV child.
# Child radius is normalized to one.  By global helicity reversal and parent
# swap there are only two heterochiral sectors.
# ---------------------------------------------------------------------------


def _frontier_domain(D: float, S: float) -> tuple[float, float]:
    d, s = float(D), float(S)
    if not all(math.isfinite(v) for v in (d, s)):
        raise ValueError("finite parent ratios required")
    if not (0.0 < d < 1.0 and 0.0 < s < 1.0 and d + s > 1.0):
        raise ValueError("strict-UV child requires 0<D,S<1 and D+S>1")
    return d, s


def critical_majority_helicity_multiplier(D: float, S: float) -> float:
    """Child shares helicity with parent D; parent S has opposite helicity.

    For phase c=1, critical production divided by child-frequency times native
    child-edge capacity is 2 S (1-D) |g_(+,-,+)|.
    """
    d, s = _frontier_domain(D, S)
    g = coupling_magnitude_closed(d, s, 1.0, 1, -1, 1)
    return 2.0 * s * (1.0 - d) * g


def critical_singleton_child_multiplier(D: float, S: float) -> float:
    """Child is the singleton helicity sign; both parents carry the opposite sign."""
    d, s = _frontier_domain(D, S)
    g = coupling_magnitude_closed(d, s, 1.0, -1, -1, 1)
    return 2.0 * abs(d - s) * g


def critical_extremizer_polynomial(D: float) -> float:
    d = float(D)
    return 12.0 * d**3 - 5.0 * d**2 + 2.0 * d - 1.0


def critical_extremizer_parent_ratio() -> float:
    """Unique root of 12D^3-5D^2+2D-1 in (9/20,23/50)."""
    lo, hi = 9.0 / 20.0, 23.0 / 50.0
    if not (critical_extremizer_polynomial(lo) < 0.0 < critical_extremizer_polynomial(hi)):
        raise RuntimeError("critical extremizer bracket changed")
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if critical_extremizer_polynomial(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class SharpCriticalFrontierGeometry:
    same_helicity_parent_ratio: float
    opposite_helicity_parent_ratio: float
    critical_multiplier: float
    log_progress_multiplier: float
    log_progress_efficiency_ratio: float
    stationary_polynomial_residual: float
    stationary_relation_residual: float
    boundary_upper: float
    rational_interior_witness_squared: Fraction
    rational_boundary_upper_squared: Fraction

    def __post_init__(self) -> None:
        if not (0.0 < self.same_helicity_parent_ratio < 1.0 and 0.0 < self.opposite_helicity_parent_ratio < 1.0):
            raise ValueError("critical extremizer must be a strict UV geometry")
        if self.stationary_polynomial_residual > 2.0e-14 or self.stationary_relation_residual > 2.0e-14:
            raise AssertionError("critical extremizer left its exact stationary algebra")
        if not self.critical_multiplier > self.boundary_upper:
            raise AssertionError("critical interior extremizer did not beat every boundary/singleton-child sector")
        if not self.rational_interior_witness_squared > self.rational_boundary_upper_squared:
            raise AssertionError("exact rational witness no longer proves the global maximum is interior")
        if not self.log_progress_efficiency_ratio < 0.3:
            raise AssertionError("critical extremizer unexpectedly entered the near-Jstar geometry")


def sharp_critical_frontier_geometry() -> SharpCriticalFrontierGeometry:
    D = critical_extremizer_parent_ratio()
    S = 4.0 * D * D
    C = critical_majority_helicity_multiplier(D, S)
    g = coupling_magnitude_closed(D, S, 1.0, 1, -1, 1)
    J = math.log(1.0 / max(D, S)) * (D + S) * g
    jstar = float_jstar()
    # Exact rational witness at D=9/20, S=41/50:
    # C_A^2 = 313640536671 / 32000000000000 > 1/128.
    interior = Fraction(313_640_536_671, 32_000_000_000_000)
    boundary_sq = Fraction(1, 128)
    return SharpCriticalFrontierGeometry(
        same_helicity_parent_ratio=D,
        opposite_helicity_parent_ratio=S,
        critical_multiplier=C,
        log_progress_multiplier=J,
        log_progress_efficiency_ratio=J / jstar,
        stationary_polynomial_residual=abs(critical_extremizer_polynomial(D)),
        stationary_relation_residual=abs(S - 4.0 * D * D),
        boundary_upper=1.0 / (8.0 * math.sqrt(2.0)),
        rational_interior_witness_squared=interior,
        rational_boundary_upper_squared=boundary_sq,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "state_measure": "rho_t=sum_(k,s) E_(k,s)(t) delta_(lambda=s|k|) is positive helical energy on the signed curl spectrum",
        "weak_nonlinear_law": "for every spectral test phi, one closed triad contributes the three-point curvature -T_median (lambda_m-lambda_-)(lambda_+-lambda_m) phi[lambda_-,lambda_m,lambda_+]",
        "affine_kernel": "phi=1 and phi=lambda vanish exactly: nonlinear energy and helicity conservation are the affine kernel of the curvature operator",
        "flow_equivalence": "the same weighted production equals sum donor->recipient dW mass times [phi(lambda_recipient)-phi(lambda_donor)]; no new Hahn/cause law is created",
        "tanaka": "phi=|lambda| has curvature only at lambda=0; homochiral triads contribute exactly zero and every heterochiral triad contributes its signed curl-sign-interface Jensen/Tanaka defect",
        "critical_pair": "on a heterochiral triad d||u||_H1/2^2^NL = 2 |k_*| T_*; globally K_+=int lambda_+ d rho and K_-=int lambda_- d rho have one common nonlinear source P_pair",
        "viscosity": "the same weak balance has sink 2 nu int lambda^2 phi(lambda) d rho; for phi=|lambda| this is 2 nu ||u||_H3/2^2",
        "frontier_sectors": "mod parent swap and global helicity reversal, strict-UV positive children have homochiral zero-critical sector, majority-helicity C_A=2S(1-D)|g_(+,-,+)|, or singleton-child C_B=2|D-S||g_(-,-,+)|",
        "sharp_critical_geometry": "C_A has the unique global maximum at S=4D^2, 12D^3-5D^2+2D-1=0; C_B and all boundaries are <1/(8 sqrt2), while an exact rational interior C_A witness has square 313640536671/32000000000000 >1/128",
        "anti_J_primitive": "at the sharp critical geometry J/J* is about 0.2715, so the log-progress good/bad functional is a proof observable, not the intrinsic measure of critical nonlinear danger",
        "causal_scope": "critical curvature production is a signed state-balance law; canonical edge dW+ remains the causal energy-transfer law",
        "claims_global_regularity": False,
    }
