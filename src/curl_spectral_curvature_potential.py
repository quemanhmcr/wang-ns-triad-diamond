from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel
from src.helical_energy_helicity_barycentric_rigidity import curl_eigenvalue

STATUS = (
    "DRAFT_EXACT_CURL_SPECTRAL_CURVATURE_POTENTIAL__"
    "CLOSED_TRIAD_SOURCE_EQUALS_SECOND_DERIVATIVE_OF_SIGNED_TENT__"
    "ENERGY_HELICITY_AFFINE_KERNEL__TANAKA_CRITICAL_VALUE__ENSTROPHY_AREA__"
    "RADIAL_TAIL_SLOPE__CRITICAL_TANAKA_SCALE_COCYCLE"
)


def _close(a: float, b: float, *, scale: float, factor: float = 1.5e-9) -> bool:
    return abs(float(a) - float(b)) <= factor * max(abs(float(scale)), abs(float(a)), abs(float(b)), 1.0e-300)


def _theta(a: float, m: float, b: float, x: float) -> float:
    if not a < m < b:
        raise ValueError("strict ordered curl eigenvalues required")
    if x < a or x > b:
        return 0.0
    if x <= m:
        return (b - m) * (x - a) / (b - a)
    return (m - a) * (b - x) / (b - a)


def _theta_slope(a: float, m: float, b: float, x: float) -> float:
    if not a < m < b:
        raise ValueError("strict ordered curl eigenvalues required")
    if x < a or x > b:
        return 0.0
    if x < m:
        return (b - m) / (b - a)
    if x > m:
        return -(m - a) / (b - a)
    # The slope jumps at m.  Tail flux is evaluated away from atoms; choosing
    # the symmetric representative prevents a fake one-sided convention.
    return 0.5 * ((b - m) - (m - a)) / (b - a)


@dataclass(frozen=True)
class TriadCurlSpectralPotential:
    ordered_lambdas: tuple[float, float, float]
    ordered_source_masses: tuple[float, float, float]
    signed_tent_mass: float
    quotient_measure_mass: float
    physical_fourier_factor: float
    native_work_mass_scale: float
    source_reconstruction_native_residual: float
    uses_later_hahn: bool = False
    creates_causal_law: bool = False
    creates_owner: bool = False

    def __post_init__(self) -> None:
        a, m, b = self.ordered_lambdas
        if not a < m < b:
            raise ValueError("curl-potential theorem requires three distinct ordered curl eigenvalues")
        if not all(math.isfinite(v) for v in (*self.ordered_source_masses, self.signed_tent_mass)):
            raise ValueError("finite signed curl-potential masses required")
        if self.quotient_measure_mass < 0.0 or self.physical_fourier_factor < 0.0 or self.native_work_mass_scale <= 0.0:
            raise ValueError("valid physical continuum normalization required")
        if self.source_reconstruction_native_residual > 1.5e-9:
            raise AssertionError("tent second derivative left the canonical rooted physical-work source")
        if self.uses_later_hahn or self.creates_causal_law or self.creates_owner:
            raise ValueError("curl potential is a state-balance representation, not new causality")

    @property
    def a(self) -> float:
        return self.ordered_lambdas[0]

    @property
    def m(self) -> float:
        return self.ordered_lambdas[1]

    @property
    def b(self) -> float:
        return self.ordered_lambdas[2]

    def potential(self, lam: float) -> float:
        return self.signed_tent_mass * _theta(self.a, self.m, self.b, float(lam))

    def slope(self, lam: float) -> float:
        return self.signed_tent_mass * _theta_slope(self.a, self.m, self.b, float(lam))

    def weighted_nonlinear_production(self, phi: Callable[[float], float]) -> float:
        vals = tuple(float(phi(lam)) for lam in self.ordered_lambdas)
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("finite spectral observable required")
        return math.fsum(w * v for w, v in zip(self.ordered_source_masses, vals))

    @property
    def potential_integral(self) -> float:
        height = (self.b - self.m) * (self.m - self.a) / (self.b - self.a)
        return self.signed_tent_mass * 0.5 * (self.b - self.a) * height

    @property
    def critical_production(self) -> float:
        # |lambda|'' = 2 delta_0 in distributions.
        return 2.0 * self.potential(0.0)

    @property
    def enstrophy_production(self) -> float:
        # (lambda^2)'' = 2.
        return 2.0 * self.potential_integral

    def radial_tail_net_production(self, radius: float) -> float:
        """Nonlinear d/dt of energy in {|lambda|>R}.

        With j=-partial_lambda kappa, the continuity equation is
        partial_t rho + partial_lambda j = 0.  Hence the exterior radial stock
        changes by j(R)-j(-R)=kappa'(-R)-kappa'(R).
        """
        r = float(radius)
        if not math.isfinite(r) or r < 0.0:
            raise ValueError("finite nonnegative radial boundary required")
        return self.slope(-r) - self.slope(r)

    @property
    def integrated_radial_first_moment_production(self) -> float:
        """Exact layer-cake integral of net radial tail production over R>=0."""
        breaks = sorted({0.0, abs(self.a), abs(self.m), abs(self.b)})
        total = 0.0
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            if hi <= lo:
                continue
            mid = 0.5 * (lo + hi)
            total += (hi - lo) * self.radial_tail_net_production(mid)
        return total


@dataclass(frozen=True)
class CurlSpectralCurvatureFamily:
    triads: tuple[TriadCurlSpectralPotential, ...]

    def __post_init__(self) -> None:
        if not self.triads:
            raise ValueError("nonempty curl-spectral potential family required")

    def potential(self, lam: float) -> float:
        return math.fsum(row.potential(lam) for row in self.triads)

    def slope(self, lam: float) -> float:
        return math.fsum(row.slope(lam) for row in self.triads)

    def weighted_nonlinear_production(self, phi: Callable[[float], float]) -> float:
        return math.fsum(row.weighted_nonlinear_production(phi) for row in self.triads)

    @property
    def critical_production(self) -> float:
        return 2.0 * self.potential(0.0)

    @property
    def enstrophy_production(self) -> float:
        return 2.0 * math.fsum(row.potential_integral for row in self.triads)

    def radial_tail_net_production(self, radius: float) -> float:
        r = float(radius)
        if r < 0.0 or not math.isfinite(r):
            raise ValueError("finite nonnegative radial boundary required")
        return self.slope(-r) - self.slope(r)


@dataclass(frozen=True)
class CriticalTanakaScaleCocycle:
    donor_frequency: float
    recipient_frontier_frequency: float
    physical_donor_work_mass: float
    critical_production: float
    normalized_critical_production: float
    radial_gap_fraction: float
    log_scale_displacement: float
    log_parabolic_lifetime_expansion: float
    critical_to_gap_margin: float
    gap_to_log_margin: float
    uses_log_progress_J: bool = False
    uses_capacity_as_causality: bool = False
    creates_budget: bool = False

    def __post_init__(self) -> None:
        if min(self.donor_frequency, self.recipient_frontier_frequency, self.physical_donor_work_mass, self.critical_production) <= 0.0:
            raise ValueError("positive spread donor/frontier/work/critical production required")
        if not self.donor_frequency < self.recipient_frontier_frequency:
            raise ValueError("critical-growing spread must have a genuine radial frontier beyond its median donor")
        if self.normalized_critical_production <= 0.0:
            raise ValueError("positive normalized critical production required")
        if self.critical_to_gap_margin < -1.5e-9 or self.gap_to_log_margin < -1.5e-12:
            raise AssertionError("critical Tanaka-scale cocycle inequality failed")
        expected = 2.0 * self.log_scale_displacement
        if not _close(self.log_parabolic_lifetime_expansion, expected, scale=max(1.0, abs(expected)), factor=5.0e-13):
            raise AssertionError("parabolic lifetime expansion is not twice log scale displacement")
        if self.uses_log_progress_J or self.uses_capacity_as_causality or self.creates_budget:
            raise ValueError("critical cocycle is intrinsic state/scale geometry, not proof currency")


def triad_curl_spectral_potential(
    triad: ClosedHelicalTriadRegistration,
    kernel: CyclicTriadMeasureKernel,
) -> TriadCurlSpectralPotential:
    """Lift one already-registered closed triad to its intrinsic curl potential.

    The canonical continuum normalization is inherited from the cyclic measure
    theorem: dW=C_F T dLambda_closed with C_F=(2pi)^(-3/2).  No new positive
    part is taken here; the potential is signed before any causal Hahn routing.
    """
    lambdas = tuple(float(curl_eigenvalue(slot.closed_mode)) for slot in triad.slots)
    order = tuple(sorted(range(3), key=lambda i: (lambdas[i], i)))
    a, m, b = (lambdas[i] for i in order)
    scale_lam = max(abs(a), abs(m), abs(b), 1.0)
    if not (m - a > 8e-12 * scale_lam and b - m > 8e-12 * scale_lam):
        raise ValueError("curl-potential theorem fails closed on curl-eigenvalue degeneracy")

    works = tuple(float(slot.signed_work) for slot in triad.slots)
    C_F = unitary_fourier_convolution_factor()
    factor = C_F * float(kernel.quotient_measure_mass)
    source = tuple(factor * works[i] for i in order)
    signed_tent_mass = -source[1]
    native = max(float(kernel.native_work_mass_scale), factor * math.fsum(float(slot.edge_registration.native_modal_capacity) for slot in triad.slots), 1.0e-300)

    wa = signed_tent_mass * (b - m) / (b - a)
    wm = -signed_tent_mass
    wb = signed_tent_mass * (m - a) / (b - a)
    residual = max(abs(source[0] - wa), abs(source[1] - wm), abs(source[2] - wb)) / native

    # Bind back to the existing cyclic recipient/donor marginals.  This catches
    # a wrong quotient/root-chart factor without using realized cancellation as
    # the numerical scale.
    for slot in triad.slots:
        i = slot.closed_mode_index
        expected_pos = factor * max(float(slot.signed_work), 0.0)
        expected_neg = factor * max(-float(slot.signed_work), 0.0)
        if not _close(kernel.recipient_edge_positive_masses[i], expected_pos, scale=native):
            raise AssertionError("curl potential lost canonical recipient dW+ normalization")
        if not _close(kernel.donor_edge_negative_masses[i], expected_neg, scale=native):
            raise AssertionError("curl potential lost canonical donor dW- normalization")

    return TriadCurlSpectralPotential(
        ordered_lambdas=(a, m, b),
        ordered_source_masses=source,
        signed_tent_mass=signed_tent_mass,
        quotient_measure_mass=float(kernel.quotient_measure_mass),
        physical_fourier_factor=factor,
        native_work_mass_scale=native,
        source_reconstruction_native_residual=residual,
    )


def critical_tanaka_scale_cocycle(
    potential: TriadCurlSpectralPotential,
) -> CriticalTanakaScaleCocycle:
    """Positive critical pair creation forces native donor-scale displacement.

    For a spread, Q=-T_m is the actual total donor work.  The two extreme
    recipients are the barycentric image of the median donor.  Therefore

      Pcrit/Q = E_recipient |lambda| - |lambda_m|
              <= N_c-N_d,

    where N_c=max(|lambda_-|,|lambda_+|) and N_d=|lambda_m|.  Dividing by N_c
    gives eta_crit <= 1-N_d/N_c <= log(N_c/N_d).  With parabolic lifetimes
    T_N=c N^-2, the last logarithm is one half of log(T_d/T_c).
    """
    if potential.signed_tent_mass <= 0.0:
        raise ValueError("critical Tanaka-scale cocycle requires a physical barycentric spread")
    pcrit = float(potential.critical_production)
    if pcrit <= 0.0:
        raise ValueError("homochiral/noncritical spread carries no positive curl-sign Tanaka production")
    Nd = abs(float(potential.m))
    Nc = max(abs(float(potential.a)), abs(float(potential.b)))
    if Nd <= 0.0 or not Nd < Nc:
        raise ValueError("positive critical spread requires a nonzero median donor below the extreme radial frontier")
    Q = float(potential.signed_tent_mass)
    eta = pcrit / (Nc * Q)
    gap = 1.0 - Nd / Nc
    logscale = math.log(Nc / Nd)
    return CriticalTanakaScaleCocycle(
        donor_frequency=Nd,
        recipient_frontier_frequency=Nc,
        physical_donor_work_mass=Q,
        critical_production=pcrit,
        normalized_critical_production=eta,
        radial_gap_fraction=gap,
        log_scale_displacement=logscale,
        log_parabolic_lifetime_expansion=2.0 * logscale,
        critical_to_gap_margin=gap - eta,
        gap_to_log_margin=logscale - gap,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "potential": "for ordered curl eigenvalues a<m<b, one closed triad has signed energy source d rho_NL = partial_lambda^2(q Theta_{a,m,b}), where Theta is the barycentric tent and q=-T_m with the canonical C_F dLambda_closed normalization",
        "full_weak_form": "summing closed triads gives partial_t rho = partial_lambda^2 kappa - 2 nu lambda^2 rho in the finite/Galerkin law and locally on continuum curl intervals where the registered Radon moments exist",
        "intrinsic_reconstruction": "energy and helicity make the nonlinear source have zero total mass and zero first lambda moment; the decaying potential is therefore uniquely reconstructed from the physical source by kappa(lambda)=1/2 int |lambda-mu| dS_NL(mu)",
        "observables": "phi=1 and lambda are the affine kernel; phi=|lambda| gives 2 kappa(0), phi=lambda^2 gives 2 int kappa, and radial exterior energy has nonlinear rate kappa'(-R)-kappa'(R)",
        "pair_source": "K_+ and K_- each receive the same nonlinear source kappa(0); critical H1/2 production, opposite-helicity pair creation and integrated net radial first-moment flux are the same state-balance quantity",
        "cocycle": "for every positive critical-growing spread, eta_crit=Pcrit/(N_c Q)<=1-N_d/N_c<=log(N_c/N_d)=1/2 log(T_d/T_c); no J, capacity, phase classifier or event budget is used",
        "constitutive_frontier": "the balance law alone is not regularity: the remaining primitive question is how the actual Waleffe quadratic constitutive law constrains the signed curvature current kappa relative to viscous killing",
        "causal_scope": "kappa is a signed state-current potential reconstructed before Hahn; canonical dW+ remains the physical causal work law",
        "claims_global_regularity": False,
    }
