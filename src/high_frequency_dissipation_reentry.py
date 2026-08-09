from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.critical_shell_service_reentry import critical_shell_bounded_service_lower

STATUS = (
    "EXACT_SQUARE_LP_HIGH_SERVICE_TO_PHYSICAL_TAIL_ENERGY__"
    "INHERITED_CRITICAL_SHELL_OR_PHYSICAL_REGENERATION__NO_RESOLVED_DV_RELABEL"
)


def dyadic_high_ratios(count: int) -> np.ndarray:
    """M_j/N=2^j for hard annuli N<|xi|<=2N, 2N<|xi|<=4N, ... ."""
    n = int(count)
    if n <= 0:
        raise ValueError("positive dyadic shell count required")
    return 2.0 ** np.arange(1, n + 1, dtype=float)


def integrated_hard_annular_currency(shell_mass_integrals: Sequence[float]) -> float:
    """D_>^hard=sum_j 2^j int mu_j d tau, mu_j=M_j||P_j u||_2^2.

    Here tau=N^2 t and M_j=2^j N.  This is an auxiliary **hard-annular**
    comparison currency.  The coherent increment theorem uses smooth LP bands;
    the two are not identified.
    """
    a = np.asarray(tuple(float(x) for x in shell_mass_integrals), dtype=float)
    if a.ndim != 1 or len(a) == 0 or np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError("finite nonnegative shell-mass integrals required")
    return float(np.dot(dyadic_high_ratios(len(a)), a))


def scaled_tail_energy_from_shell_masses(shell_critical_masses: Sequence[float]) -> float:
    """N||P_>N u||_2^2=sum_j 2^{-j} mu_j at one time.

    Since sum_{j>=1}2^{-j}=1, the scaled tail energy is a convex dyadic
    average of critical shell masses.  Hence some shell has mu_j at least this
    tail energy.  This is the key summability that high enstrophy itself lacks.
    """
    a = np.asarray(tuple(float(x) for x in shell_critical_masses), dtype=float)
    if a.ndim != 1 or len(a) == 0 or np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError("finite nonnegative shell critical masses required")
    weights = 1.0 / dyadic_high_ratios(len(a))
    return float(np.dot(weights, a))


def high_tail_scaled_gradient_bounds(shell_mass_integrals: Sequence[float]) -> tuple[float, float]:
    """Hard-annulus spectral bounds for N int ||grad P_>N u||_2^2 dt.

    On M_j/2<|xi|<=M_j,

      (1/4) 2^j mu_j <= N^{-1}||grad P_j u||_2^2 <= 2^j mu_j

    in scaled time.  After integration,

      D_>^hard/4 <= N int||grad P_>N u||_2^2 dt <= D_>^hard.
    """
    D = integrated_hard_annular_currency(shell_mass_integrals)
    return D / 4.0, D


def lp_tail_comparison_constant(
    lower_support_fraction: float = 0.5,
    square_bessel_upper: float = 1.0,
) -> float:
    """Exact LP-to-hard-tail constant from support plus the L2 square/Bessel law.

    Suppose the same smooth high LP multipliers used upstream obey

      supp phi_j subset {|xi|>=a M_j},   sum_(j>=1)|phi_j(xi)|^2 <= B,

    with `a>=1/2`, `M_j=2^jN`.  Then every high multiplier lies in `|xi|>=N` and

      sum_j (M_j^2/N)||phi_j(D)u||_2^2
      <= (B/a^2) N^-1 ||grad P_>N u||_2^2.

    Thus the physical hard-tail dissipation dominates the LP high currency by
    `c_LP=a^2/B`.  A square-normalized smooth dyadic partition has `B=1`; with
    the standard lower annular support `a=1/2`, the canonical constant is 1/4.
    """
    a = float(lower_support_fraction)
    B = float(square_bessel_upper)
    if a < 0.5 or B <= 0 or not math.isfinite(a + B):
        raise ValueError("require finite lower support fraction >=1/2 and positive square-Bessel upper")
    return a * a / B


def canonical_square_lp_tail_comparison_constant() -> float:
    """Canonical smooth square-LP choice: lower annular support 1/2, Bessel upper 1."""
    return lp_tail_comparison_constant(0.5, 1.0)


def physical_tail_dissipation_lower_from_lp(
    lp_high_currency: float,
    lp_to_physical_tail_lower: float,
) -> float:
    """Certified supplier from a chosen smooth LP high currency to physical tail dissipation.

    The coherent-increment `d_high` comes from a standard smooth LP partition.
    Its comparison with the orthogonal hard-tail gradient is partition-dependent:

      N int||grad P_>N u||_2^2 dt >= c_LP D_high.

    This function keeps `c_LP` explicit.  The canonical smooth square-normalized
    LP analysis frame registered by coherent increment service has lower support
    fraction 1/2 and L2 square-Bessel constant 1, hence `c_LP=1/4`.  Other LP
    frames may enter only with their own certified comparison constant.
    """
    D = float(lp_high_currency)
    c = float(lp_to_physical_tail_lower)
    if D <= 0 or c <= 0 or not math.isfinite(D + c):
        raise ValueError("positive finite LP currency and comparison lower required")
    return c * D


def direct_high_enstrophy_shell_counterexample(level: int, high_currency: float) -> dict[str, float]:
    """One-shell family proving D_high alone has no scale-independent mass floor.

    Put all high currency at M_j=2^j N and choose mu_j=D/2^j.  Then
    2^j mu_j=D while mu_j ->0 as j->infinity.
    """
    j = int(level)
    D = float(high_currency)
    if j < 1 or D <= 0 or not math.isfinite(D):
        raise ValueError("level>=1 and positive finite high currency required")
    mu = D / (2.0**j)
    return {
        "level": float(j),
        "frequency_ratio": 2.0**j,
        "critical_shell_mass": mu,
        "high_currency": (2.0**j) * mu,
    }


def high_tail_energy_owner_threshold(physical_tail_dissipation: float, viscosity: float) -> float:
    """Clean half-pigeonhole threshold `nu D_tail`.

    The native physical currency is

      D_tail = N int ||grad P_>N u||_2^2 dt.

    Exact hard-tail energy gives

      N E_>(s) + N W_>^+ >= 2 nu D_tail.

    Hence inherited scaled tail energy or positive scaled nonlinear tail work is
    at least `nu D_tail`.  Exact ties are retained jointly.
    """
    D = float(physical_tail_dissipation)
    nu = float(viscosity)
    if D <= 0 or nu <= 0 or not math.isfinite(D + nu):
        raise ValueError("positive finite physical tail dissipation and viscosity required")
    return nu * D


def classify_high_tail_energy_owners(
    physical_tail_dissipation_lower: float,
    viscosity: float,
    inherited_scaled_tail_energy: float,
    positive_scaled_tail_work: float,
) -> dict[str, object]:
    """Classify the two native owners forced by the exact hard-tail energy law.

    The first input is any certified lower bound for the **physical** currency
    `D_tail=N int||grad P_>N u||_2^2 dt`.  Inputs two and three are
    `N||P_>N u(s)||_2^2` and `N W_>^+`, with

      W_>^+=int 2[Re<P_>N u,-P_>N P div(u tensor u)>]_+ dt.

    No smooth LP scalar is identified with the hard tail: an LP supplier must
    first pass through `physical_tail_dissipation_lower_from_lp`.
    """
    D = float(physical_tail_dissipation_lower)
    nu = float(viscosity)
    E = float(inherited_scaled_tail_energy)
    W = float(positive_scaled_tail_work)
    if min(D, nu) <= 0 or min(E, W) < 0 or not all(math.isfinite(x) for x in (D, nu, E, W)):
        raise ValueError("valid positive physical-tail data required")
    required_total = 2.0 * nu * D
    tol = 4e-13 * max(1.0, required_total, E + W)
    if E + W + tol < required_total:
        raise ValueError("data violate the physical hard-tail energy lower required by the theorem")
    threshold = high_tail_energy_owner_threshold(D, nu)
    owners: list[str] = []
    if E + tol >= threshold:
        owners.append("inherited_tail_energy")
    if W + tol >= threshold:
        owners.append("positive_nonlinear_regeneration")
    if not owners:
        raise AssertionError("two-owner hard-tail pigeonhole failed")
    return {
        "required_total_scaled_energy_or_work": required_total,
        "owner_threshold": threshold,
        "joint_owners": tuple(owners),
        "resolved_DV_relabel": False,
    }


def inherited_tail_shell_witness(shell_critical_masses: Sequence[float]) -> dict[str, float]:
    """Return a hard high shell whose critical mass dominates scaled tail energy.

    For a finite truncation the dyadic weights sum to <1, so max(mu_j) is at
    least the represented tail energy.  Infinite tails follow by monotone
    convergence.  The shell frequency is M_j=2^j N and therefore lies strictly
    above the source block scale.
    """
    a = np.asarray(tuple(float(x) for x in shell_critical_masses), dtype=float)
    if a.ndim != 1 or len(a) == 0 or np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError("finite nonnegative shell critical masses required")
    tail = scaled_tail_energy_from_shell_masses(a)
    j0 = int(np.argmax(a))
    mu = float(a[j0])
    tol = 3e-13 * max(1.0, mu, tail)
    if mu + tol < tail:
        raise AssertionError("dyadic tail energy failed to expose a critical shell")
    j = j0 + 1
    return {
        "dyadic_level": float(j),
        "shell_to_block_frequency_ratio": 2.0**j,
        "critical_shell_mass": mu,
        "scaled_tail_energy": tail,
        "critical_shell_margin": mu - tail,
    }


def inherited_branch_clean_shell_mass_lower(physical_tail_dissipation_lower: float, viscosity: float) -> float:
    """If inherited tail owns the gate, some high shell has `mu>=nu D_tail`."""
    return high_tail_energy_owner_threshold(physical_tail_dissipation_lower, viscosity)


def lp_inherited_branch_clean_shell_mass_lower(
    lp_high_currency: float,
    lp_to_physical_tail_lower: float,
    viscosity: float,
) -> float:
    """Smooth-LP supplier: `mu>=nu c_LP D_high` on the inherited owner branch."""
    D_tail = physical_tail_dissipation_lower_from_lp(
        lp_high_currency, lp_to_physical_tail_lower
    )
    return inherited_branch_clean_shell_mass_lower(D_tail, viscosity)


def inherited_branch_full_survivor_service_lower(
    physical_tail_dissipation_lower: float,
    viscosity: float,
    scaled_lifetime: float,
) -> float:
    """Existing generic-shell service lower, conditional on a full natural survivor.

    This is deliberately not an automatic outcome of the high-tail theorem: the
    generic critical-shell theorem still enforces its own observed-history guard
    and may instead hit strain/interface/HH or t=0.
    """
    mu0 = inherited_branch_clean_shell_mass_lower(physical_tail_dissipation_lower, viscosity)
    return critical_shell_bounded_service_lower(mu0, scaled_lifetime, viscosity)


def positive_shell_work_disintegration(
    positive_scaled_tail_work: float,
    positive_scaled_shell_works: Sequence[float],
) -> dict[str, float]:
    """Tail positive work -> actual hard-shell positive-work law.

    `positive_scaled_shell_works[j-1]=N W_j^+`.  Orthogonality gives

      W_>^+ <= sum_j W_j^+.

    Since M_j/N=2^j>=2,

      sum_j M_j W_j^+ = sum_j 2^j (N W_j^+) >= 2 N W_>^+.

    No single-shell mass conclusion is drawn from this possibly infinite work
    law.  It remains actual positive physical work.
    """
    W = float(positive_scaled_tail_work)
    a = np.asarray(tuple(float(x) for x in positive_scaled_shell_works), dtype=float)
    if W < 0 or not math.isfinite(W) or a.ndim != 1 or len(a) == 0 or np.any(~np.isfinite(a)) or np.any(a < 0):
        raise ValueError("finite nonnegative tail/shell positive works required")
    tol = 4e-13 * max(1.0, W, float(a.sum()))
    if float(a.sum()) + tol < W:
        raise ValueError("hard-shell positive works do not dominate positive tail work")
    ratios = dyadic_high_ratios(len(a))
    ownscale = float(np.dot(ratios, a))
    if ownscale + tol < 2.0 * W:
        raise AssertionError("own-scale hard-shell positive-work lower failed")
    return {
        "scaled_tail_positive_work": W,
        "scaled_shell_positive_work_sum": float(a.sum()),
        "own_scale_positive_shell_work": ownscale,
        "clean_own_scale_lower": 2.0 * W,
    }


def classify_regeneration_work_owners(
    own_scale_positive_shell_work: float,
    positive_hh_work: float,
    positive_resolved_interface_work: float,
) -> dict[str, object]:
    """Low-low-free shell work -> HH or resolved cross/interface physical work.

    On shell M_j with V=S_{M_j/4}u, the hard shell lies in |xi|>M_j/2 while
    B(V,V) lies in |xi|<=M_j/2, so low-low work vanishes.  Expanding u=V+h gives
    the exact signed shell work as HH plus the mixed resolved-cross/interface
    work.  Therefore positive parts obey

      W_shell^+ <= W_HH^+ + W_interface^+.

    The inputs here are the aggregate own-scale positive measures after summing
    shell-time atoms.  Large interface work is not renamed HH work, and large HH
    work is not yet declared generated-energy productivity without its own energy
    gate.
    """
    S = float(own_scale_positive_shell_work)
    H = float(positive_hh_work)
    I = float(positive_resolved_interface_work)
    if S < 0 or H < 0 or I < 0 or not all(math.isfinite(x) for x in (S, H, I)):
        raise ValueError("finite nonnegative work measures required")
    tol = 4e-13 * max(1.0, S, H + I)
    if H + I + tol < S:
        raise ValueError("HH/interface positive works do not cover shell positive work")
    threshold = S / 2.0
    owners: list[str] = []
    if H + tol >= threshold:
        owners.append("positive_HH_regeneration")
    if I + tol >= threshold:
        owners.append("positive_resolved_cross_interface")
    if S > 0 and not owners:
        raise AssertionError("low-low-free work owner pigeonhole failed")
    return {
        "owner_threshold": threshold,
        "joint_owners": tuple(owners),
        "HH_is_productivity_generated_branch": False,
        "interface_is_free": False,
    }


def high_tail_clean_reentry_thresholds(physical_tail_dissipation_lower: float, viscosity: float) -> dict[str, float | str]:
    """Clean constants from a physical tail-dissipation lower bound."""
    D = float(physical_tail_dissipation_lower)
    nu = float(viscosity)
    owner = high_tail_energy_owner_threshold(D, nu)
    return {
        "physical_tail_dissipation_lower": D,
        "energy_or_tail_work_owner": owner,
        "inherited_critical_shell_mass": owner,
        "positive_tail_work_scaled": owner,
        "own_scale_shell_work_if_regeneration": 2.0 * owner,
        "HH_or_interface_work_if_regeneration": owner,
        "master_semantics": "RECURSE_CRITICAL",
    }


def lp_high_clean_reentry_thresholds(
    lp_high_currency: float,
    lp_to_physical_tail_lower: float,
    viscosity: float,
) -> dict[str, float | str]:
    """Clean source-facing thresholds with the LP/hard-tail comparison explicit."""
    D_tail = physical_tail_dissipation_lower_from_lp(
        lp_high_currency, lp_to_physical_tail_lower
    )
    out = high_tail_clean_reentry_thresholds(D_tail, viscosity)
    return {
        "lp_high_currency": float(lp_high_currency),
        "lp_to_physical_tail_lower": float(lp_to_physical_tail_lower),
        **out,
    }


def theorem_certificate() -> dict[str, object]:
    counter = direct_high_enstrophy_shell_counterexample(40, 1.0)
    if counter["critical_shell_mass"] >= 1e-10:
        raise AssertionError("counterexample did not expose absence of a direct high-enstrophy mass floor")
    c_lp = canonical_square_lp_tail_comparison_constant()
    if abs(c_lp - 0.25) > 1e-15:
        raise AssertionError("canonical square-LP tail comparison lost its one-quarter constant")
    return {
        "status": STATUS,
        "native_currency": "D_tail=N int||grad P_>N u||_2^2 dt is the physical orthogonal-tail dissipation currency; it is neither smooth-LP d_high nor resolved low-pass D_V",
        "lp_supplier": "if high LP multipliers satisfy supp phi_j subset {|xi|>=a M_j}, a>=1/2, and sum|phi_j|^2<=B, then D_tail>=(a^2/B)D_high; the canonical square-normalized smooth dyadic choice a=1/2,B=1 gives c_LP=1/4 exactly",
        "spectral_bridge": "for hard annular D_>^hard=sum 2^j int mu_j d tau, D_>^hard/4 <= D_tail <= D_>^hard",
        "energy_gate": "N||P_>N u(s)||_2^2 + N W_>^+ >= 2 nu D_tail, so inherited tail energy or actual positive nonlinear tail work is >=nu D_tail",
        "inherited_branch": "N||P_>N u||^2=sum 2^-j mu_j with sum 2^-j=1; therefore some M_j>=2N has mu_j>=nu D_tail (or >=nu c_LP D_high for an LP supplier) and enters the generic critical-shell theorem",
        "regeneration_branch": "positive tail work disintegrates into hard-shell positive work; own-scale shell work is at least twice N W_>^+; low-low is support-excluded at every shell, so positive shell work is covered by HH plus resolved-cross/interface positive work",
        "no_false_productivity": "the HH regeneration owner is actual positive work but is not automatically the generated-energy branch W_HH>=8E1/15; its own physical energy gate is still required before KL productivity",
        "anti_relabel": "D_high alone has no critical-shell floor: at level j choose mu_j=2^-j D_high; smooth LP d_high is not identified with the hard tail without c_LP; no high-frequency owner is renamed resolved low-pass D_V",
        "master_rule": "both inherited shell and physical regeneration/interface are recursive scale-sensitive owners; no additive reset is created",
    }


@dataclass(frozen=True)
class HighFrequencyReentryStress:
    samples: int
    minimum_gradient_lower_margin: float
    minimum_gradient_upper_margin: float
    minimum_energy_owner_margin: float
    minimum_inherited_shell_margin: float
    minimum_shell_work_disintegration_margin: float
    minimum_regeneration_owner_margin: float
    maximum_joint_energy_owner_count: int
    maximum_joint_regeneration_owner_count: int


def stress(samples: int = 50_000, seed: int = 20260809) -> HighFrequencyReentryStress:
    rng = np.random.default_rng(seed)
    mgl = mgu = me = mis = msw = mro = float("inf")
    max_e = max_r = 0
    for _ in range(samples):
        n = int(rng.integers(1, 18))
        U = np.exp(rng.uniform(-12.0, 3.0, size=n))
        ratios = dyadic_high_ratios(n)
        D = integrated_hard_annular_currency(U)

        # Exact shell spectral radii r_j=|xi|/M_j in (1/2,1] give the true
        # scaled gradient integral G=sum 2^j r_j^2 int mu_j d tau.
        radial = rng.uniform(0.5, 1.0, size=n)
        G = float(np.dot(ratios * radial * radial, U))
        lo, hi = high_tail_scaled_gradient_bounds(U)
        mgl = min(mgl, G - lo)
        mgu = min(mgu, hi - G)
        if G + 3e-12 * max(1.0, D) < lo or G > hi + 3e-12 * max(1.0, D):
            raise AssertionError("hard high-tail spectral gradient bounds failed")

        nu = float(rng.uniform(0.03, 3.0))
        # Exact energy law only makes the available E0+W larger than 2nu G.
        total = 2.0 * nu * G * float(rng.uniform(1.0, 2.5))
        q = float(rng.uniform(0.0, 1.0))
        E = q * total
        W = (1.0 - q) * total
        gate = classify_high_tail_energy_owners(G, nu, E, W)
        threshold = float(gate["owner_threshold"])
        hard_supplier = physical_tail_dissipation_lower_from_lp(D, 0.25)
        if hard_supplier > G + 3e-12 * max(1.0, D):
            raise AssertionError("hard-annulus LP supplier exceeded physical tail dissipation")
        me = min(me, max(E, W) - threshold)
        max_e = max(max_e, len(tuple(gate["joint_owners"])))
        if max(E, W) + 3e-12 * max(1.0, threshold) < threshold:
            raise AssertionError("high-tail energy owner threshold failed")

        if E > 0:
            # Build shell masses with exact represented tail energy E.
            p = rng.dirichlet(np.ones(n))
            masses = E * p * ratios
            tail = scaled_tail_energy_from_shell_masses(masses)
            if abs(tail - E) > 3e-12 * max(1.0, E):
                raise AssertionError("synthetic inherited tail energy identity failed")
            wit = inherited_tail_shell_witness(masses)
            mis = min(mis, float(wit["critical_shell_margin"]))

        if W > 0:
            p = rng.dirichlet(np.ones(n))
            # Shell positive work dominates tail positive work; allow slack.
            shell_scaled = W * float(rng.uniform(1.0, 2.0)) * p
            law = positive_shell_work_disintegration(W, shell_scaled)
            own = float(law["own_scale_positive_shell_work"])
            msw = min(msw, own - 2.0 * W)

            r = float(rng.uniform(0.0, 1.0))
            slack = own * float(rng.uniform(1.0, 1.8))
            H = r * slack
            I = (1.0 - r) * slack
            split = classify_regeneration_work_owners(own, H, I)
            mro = min(mro, max(H, I) - own / 2.0)
            max_r = max(max_r, len(tuple(split["joint_owners"])))

        level = int(rng.integers(10, 60))
        cex = direct_high_enstrophy_shell_counterexample(level, D)
        if abs(float(cex["high_currency"]) - D) > 2e-12 * max(1.0, D):
            raise AssertionError("direct high-enstrophy counterexample lost its currency")

    return HighFrequencyReentryStress(
        samples,
        mgl,
        mgu,
        me,
        mis,
        msw,
        mro,
        max_e,
        max_r,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-frequency-dissipation-reentry"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "high_frequency_dissipation_reentry.json").write_text(
        json.dumps(data, indent=2), encoding="utf-8"
    )
    md = f"""# High-frequency service dissipation -> inherited critical shell or physical regeneration

Status: **{cert['status']}**.

The physical theorem is formulated first in the PDE's own currency

`D_tail=N int ||grad P_>N u||_2^2 dt`.

This orthogonal hard-tail dissipation is **not** resolved low-pass `D_V`, and it is not silently identified with the smooth LP `d_high` used by coherent increment service.  For high LP analysis multipliers supported above `a M_j` with square-Bessel upper `B`, Plancherel gives

`D_tail >= c_LP D_high`,  `c_LP=a^2/B`.

The canonical smooth square-normalized analysis--synthesis frame registered upstream has `a=1/2`, `B=1`, hence `c_LP=1/4`.  Independently, for the auxiliary hard annuli `M_j=2^jN`, `M_j/2<|xi|<=M_j`, the same one-quarter spectral lower holds.  Writing `mu_j=M_j||P_j u||_2^2`, the hard annular currency is `D_>^hard=int sum_j2^j mu_j d tau`.

Indeed `D_high` alone cannot force a critical shell: placing all currency at level `j` with `mu_j=2^-j D_high` keeps `2^j mu_j=D_high` while `mu_j->0`.

The hard-annulus comparison gives

`D_>^hard/4 <= D_tail <= D_>^hard`.

The missing physics is then viscosity plus the exact hard-tail energy law.  With `w=P_>N u` and actual positive nonlinear tail work

`W_>^+=int 2[Re<w,-P_>N P div(u tensor u)>]_+ dt`,

the Navier--Stokes energy identity gives

`N||w(s)||_2^2 + N W_>^+ >= 2 nu D_tail`.

Therefore at least one native owner carries `nu D_tail`.  With an LP supplier this is at least `nu c_LP D_high`; for the hard-annulus lower it is `nu D_>^hard/4` (exact ties remain joint):

1. **Inherited tail energy.**  Since

   `N||w(s)||_2^2=sum_j 2^-j mu_j(s)`

   and `sum_j2^-j=1`, some actual high shell `M_j>=2N` obeys

   `M_j||P_j u(s)||_2^2 >= nu D_tail`,

   hence at least `nu c_LP D_high` for an LP supplier.

   This is a genuine critical-shell seed and enters the existing generic shell first-stop/service theorem.  That theorem's observed-history guard is unchanged: the present theorem does not manufacture a full natural survivor.

2. **Positive nonlinear regeneration.**  Orthogonality disintegrates tail work into hard shell positive works `W_j^+` with `sum W_j^+>=W_>^+`.  Because every `M_j/N>=2`,

   `sum_j M_j W_j^+ >= 2N W_>^+`.

   At shell `M_j`, choose `V=S_(M_j/4)u`.  The hard shell is strictly above `M_j/2`, whereas `B(V,V)` is supported at or below `M_j/2`; low--low work is absent.  Expanding `u=V+h`, the signed shell work is exactly HH work plus resolved mixed/cross-interface work.  Hence positive parts satisfy

   `W_shell^+ <= W_HH^+ + W_interface^+`.

   On the clean regeneration branch the aggregate own-scale shell work is at least `2 nu D_tail`, so HH or resolved-interface positive work carries at least `nu D_tail` (at least `nu c_LP D_high` for the LP supplier).

The last HH statement is **not** the generated-energy gate.  Actual positive HH work still must pass its own energy comparison before the physical KL productivity theorem may be invoked.  Likewise interface work remains interface/strain provenance rather than being declared free.

Thus the anonymous high-frequency enstrophy exit has been replaced by native physical owners:

`smooth-LP D_high -> physical D_tail -> inherited critical shell OR actual positive HH/interface regeneration`.

No `D_high -> D_V` relabel, no packet persistence, no additive reset.

Stress: `{out.samples}` high-tail/shell/work states
- minimum hard-tail gradient lower margin: `{out.minimum_gradient_lower_margin:.3e}`
- minimum hard-tail gradient upper margin: `{out.minimum_gradient_upper_margin:.3e}`
- minimum two-owner energy margin: `{out.minimum_energy_owner_margin:.3e}`
- minimum inherited-shell margin: `{out.minimum_inherited_shell_margin:.3e}`
- minimum shell-work disintegration margin: `{out.minimum_shell_work_disintegration_margin:.3e}`
- minimum HH/interface owner margin: `{out.minimum_regeneration_owner_margin:.3e}`
- maximum sampled joint energy-owner count: `{out.maximum_joint_energy_owner_count}`
- maximum sampled joint regeneration-owner count: `{out.maximum_joint_regeneration_owner_count}`

This theorem closes the unit mismatch in the coherent-service `D_high` branch at the physical-energy level.  Supplier-specific continuation of a positive regeneration owner and the low-frequency pressure-reservoir lineage remain separate master-facing questions.  No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
