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
    "EXACT_HIGH_FREQUENCY_SERVICE_DISSIPATION_TO_INHERITED_CRITICAL_SHELL_OR_PHYSICAL_REGENERATION__"
    "NO_RESOLVED_DV_RELABEL"
)


def dyadic_high_ratios(count: int) -> np.ndarray:
    """M_j/N=2^j for hard annuli N<|xi|<=2N, 2N<|xi|<=4N, ... ."""
    n = int(count)
    if n <= 0:
        raise ValueError("positive dyadic shell count required")
    return 2.0 ** np.arange(1, n + 1, dtype=float)


def integrated_high_lp_currency(shell_mass_integrals: Sequence[float]) -> float:
    """D_>^LP=sum_j 2^j int mu_j d tau, mu_j=M_j||P_j u||_2^2.

    Here tau=N^2 t and M_j=2^j N.  This is the native high-frequency
    normalized-enstrophy currency used by the coherent increment service route.
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

      D_>^LP/4 <= N int||grad P_>N u||_2^2 dt <= D_>^LP.
    """
    D = integrated_high_lp_currency(shell_mass_integrals)
    return D / 4.0, D


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


def high_tail_energy_owner_threshold(high_currency: float, viscosity: float) -> float:
    """Clean half-pigeonhole threshold nu D_high/4.

    Hard-tail energy gives

      N E_>(s) + N W_>^+ >= 2 nu N int||grad w||^2 dt
                              >= (nu/2) D_>^LP.

    Hence inherited scaled tail energy or positive scaled nonlinear tail work is
    at least nu D_>^LP/4.  Exact ties are retained jointly.
    """
    D = float(high_currency)
    nu = float(viscosity)
    if D <= 0 or nu <= 0 or not math.isfinite(D + nu):
        raise ValueError("positive finite high currency and viscosity required")
    return nu * D / 4.0


def classify_high_tail_energy_owners(
    high_currency: float,
    viscosity: float,
    inherited_scaled_tail_energy: float,
    positive_scaled_tail_work: float,
) -> dict[str, object]:
    """Classify the two native owners forced by the exact hard-tail energy law.

    Inputs are `N||P_>N u(s)||_2^2` and `N W_>^+`, with

      W_>^+=int 2[Re<P_>N u,-P_>N P div(u tensor u)>]_+ dt.

    The function checks only the exact consequence needed by the theorem.  It
    does not infer work from a coefficient impulse and it does not call D_high
    resolved D_V.
    """
    D = float(high_currency)
    nu = float(viscosity)
    E = float(inherited_scaled_tail_energy)
    W = float(positive_scaled_tail_work)
    if min(D, nu) <= 0 or min(E, W) < 0 or not all(math.isfinite(x) for x in (D, nu, E, W)):
        raise ValueError("valid positive high-tail data required")
    required_total = nu * D / 2.0
    tol = 4e-13 * max(1.0, required_total, E + W)
    if E + W + tol < required_total:
        raise ValueError("data violate the hard-tail energy lower required by the theorem")
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


def inherited_branch_clean_shell_mass_lower(high_currency: float, viscosity: float) -> float:
    """If inherited tail owns the energy gate, some high shell has mu>=nu D/4."""
    return high_tail_energy_owner_threshold(high_currency, viscosity)


def inherited_branch_full_survivor_service_lower(
    high_currency: float,
    viscosity: float,
    scaled_lifetime: float,
) -> float:
    """Existing generic-shell service lower, conditional on a full natural survivor.

    This is deliberately not an automatic outcome of the high-tail theorem: the
    generic critical-shell theorem still enforces its own observed-history guard
    and may instead hit strain/interface/HH or t=0.
    """
    mu0 = inherited_branch_clean_shell_mass_lower(high_currency, viscosity)
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


def high_tail_clean_reentry_thresholds(high_currency: float, viscosity: float) -> dict[str, float | str]:
    """Clean constants for the complete high-tail dichotomy."""
    D = float(high_currency)
    nu = float(viscosity)
    owner = high_tail_energy_owner_threshold(D, nu)
    return {
        "high_currency": D,
        "energy_or_tail_work_owner": owner,
        "inherited_critical_shell_mass": owner,
        "positive_tail_work_scaled": owner,
        "own_scale_shell_work_if_regeneration": 2.0 * owner,
        "HH_or_interface_work_if_regeneration": owner,
        "master_semantics": "RECURSE_CRITICAL",
    }


def theorem_certificate() -> dict[str, object]:
    counter = direct_high_enstrophy_shell_counterexample(40, 1.0)
    if counter["critical_shell_mass"] >= 1e-10:
        raise AssertionError("counterexample did not expose absence of a direct high-enstrophy mass floor")
    return {
        "status": STATUS,
        "native_currency": "D_>^LP=int sum_(j>=1) 2^j mu_j d tau on hard annuli M_j=2^j N; it is high-frequency normalized enstrophy, not resolved D_V",
        "spectral_bridge": "D_>^LP/4 <= N int ||grad P_>N u||_2^2 dt <= D_>^LP",
        "energy_gate": "N||P_>N u(s)||_2^2 + N W_>^+ >= (nu/2)D_>^LP, so inherited tail energy or actual positive nonlinear tail work is >=nu D_>^LP/4",
        "inherited_branch": "N||P_>N u||^2=sum 2^-j mu_j with sum 2^-j=1; therefore some M_j>=2N has mu_j>=nu D_>^LP/4 and enters the existing generic critical-shell theorem",
        "regeneration_branch": "positive tail work disintegrates into hard-shell positive work; own-scale shell work is at least twice N W_>^+; low-low is support-excluded at every shell, so positive shell work is covered by HH plus resolved-cross/interface positive work",
        "no_false_productivity": "the HH regeneration owner is actual positive work but is not automatically the generated-energy branch W_HH>=8E1/15; its own physical energy gate is still required before KL productivity",
        "anti_relabel": "D_high alone has no critical-shell floor: at level j choose mu_j=2^-j D_high; also no high-frequency owner is renamed resolved low-pass D_V",
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
        D = integrated_high_lp_currency(U)

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
        gate = classify_high_tail_energy_owners(D, nu, E, W)
        threshold = float(gate["owner_threshold"])
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

The high-frequency branch of coherent increment service is kept in its native units.  On hard dyadic annuli

`M_j=2^j N`, `M_j/2<|xi|<=M_j`, `j>=1`,

write `mu_j=M_j||P_j u||_2^2` and

`D_>^LP=int sum_j 2^j mu_j d tau`, `tau=N^2 t`.

This is **not** resolved low-pass `D_V`.  Indeed `D_high` alone cannot force a critical shell: placing all currency at level `j` with `mu_j=2^-j D_high` keeps `2^j mu_j=D_high` while `mu_j->0`.

The missing physics is viscosity plus the exact hard-tail energy law.  For `w=P_>N u`, hard-annulus support gives

`D_>^LP/4 <= N int ||grad w||_2^2 dt <= D_>^LP`.

With actual positive nonlinear tail work

`W_>^+=int 2[Re<w,-P_>N P div(u tensor u)>]_+ dt`,

the Navier--Stokes energy identity gives

`N||w(s)||_2^2 + N W_>^+ >= 2 nu N int||grad w||_2^2 dt >= (nu/2)D_>^LP`.

Therefore at least one native owner carries `nu D_>^LP/4` (exact ties remain joint):

1. **Inherited tail energy.**  Since

   `N||w(s)||_2^2=sum_j 2^-j mu_j(s)`

   and `sum_j2^-j=1`, some actual high shell `M_j>=2N` obeys

   `M_j||P_j u(s)||_2^2 >= nu D_>^LP/4`.

   This is a genuine critical-shell seed and enters the existing generic shell first-stop/service theorem.  That theorem's observed-history guard is unchanged: the present theorem does not manufacture a full natural survivor.

2. **Positive nonlinear regeneration.**  Orthogonality disintegrates tail work into hard shell positive works `W_j^+` with `sum W_j^+>=W_>^+`.  Because every `M_j/N>=2`,

   `sum_j M_j W_j^+ >= 2N W_>^+`.

   At shell `M_j`, choose `V=S_(M_j/4)u`.  The hard shell is strictly above `M_j/2`, whereas `B(V,V)` is supported at or below `M_j/2`; low--low work is absent.  Expanding `u=V+h`, the signed shell work is exactly HH work plus resolved mixed/cross-interface work.  Hence positive parts satisfy

   `W_shell^+ <= W_HH^+ + W_interface^+`.

   On the clean regeneration branch the aggregate own-scale shell work is at least `nu D_>^LP/2`, so HH or resolved-interface positive work carries at least `nu D_>^LP/4`.

The last HH statement is **not** the generated-energy gate.  Actual positive HH work still must pass its own energy comparison before the physical KL productivity theorem may be invoked.  Likewise interface work remains interface/strain provenance rather than being declared free.

Thus the anonymous high-frequency enstrophy exit has been replaced by native physical owners:

`D_high -> inherited critical shell OR actual positive HH/interface regeneration`.

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
