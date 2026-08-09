from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower

TRANSPORTER_RADIUS = 1.0 / 4.0
LOW_STRAIN_THRESHOLD = 1.0 / 30.0


def dyadic_resolved_shell_upper(child_frequency: float, shell_index: int) -> float:
    """Upper radius M_j=(N/4)2^{-j} of the j-th resolved dyadic annulus."""
    N = float(child_frequency)
    j = int(shell_index)
    if N <= 0 or not math.isfinite(N) or j < 0:
        raise ValueError("positive child frequency and nonnegative shell index required")
    return TRANSPORTER_RADIUS * N * (2.0 ** (-j))


def infinite_resolved_shell_upper_sum(child_frequency: float) -> float:
    """sum_{j>=0} M_j=N/2."""
    N = float(child_frequency)
    if N <= 0 or not math.isfinite(N):
        raise ValueError("positive finite child frequency required")
    return 0.5 * N


def low_mass_dissipation_upper(
    child_frequency: float,
    scaled_lifetime: float,
    critical_mass_threshold: float,
) -> float:
    """Upper normalized D_V mass carried by low-critical-mass dyadic atoms.

    Let A_j={M_j/2<|xi|<=M_j}, E_j(t)=||P_j u(t)||_2^2 and
    mu_j(t)=M_j E_j(t).  For the standard strict low-pass multiplier
    V=S_(N/4)u with |s_N(xi)|<=1,

      N ||grad P_j V||_2^2 <= N M_j mu_j(t).

    Therefore on {mu_j<mu_*}, over T=cN^-2,

      D_bad <= N T mu_* sum_j M_j = c mu_*/2.
    """
    N = float(child_frequency)
    c = float(scaled_lifetime)
    mu = float(critical_mass_threshold)
    if N <= 0 or c <= 0 or mu < 0 or not all(math.isfinite(x) for x in (N, c, mu)):
        raise ValueError("positive finite scale/lifetime and nonnegative mass threshold required")
    T = c / (N * N)
    return N * T * mu * infinite_resolved_shell_upper_sum(N)


def high_strain_ancestor_mass_threshold(scaled_lifetime: float) -> float:
    """Critical shell mass threshold D_*/c used for the half-dissipation law."""
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    return clean_high_strain_dissipation_lower(c) / c


def retained_dissipation_lower(
    total_normalized_dissipation: float,
    scaled_lifetime: float,
    critical_mass_threshold: float,
) -> float:
    D = float(total_normalized_dissipation)
    bad = low_mass_dissipation_upper(1.0, scaled_lifetime, critical_mass_threshold)
    if D < 0 or not math.isfinite(D):
        raise ValueError("finite nonnegative normalized dissipation required")
    return max(0.0, D - bad)


def retained_fraction_lower(
    total_normalized_dissipation: float,
    scaled_lifetime: float,
    critical_mass_threshold: float,
) -> float:
    D = float(total_normalized_dissipation)
    if D <= 0:
        return 0.0
    return retained_dissipation_lower(D, scaled_lifetime, critical_mass_threshold) / D


def ancestor_lifetime_ratio(child_frequency: float, ancestor_frequency: float) -> float:
    """Parabolic natural-lifetime ratio T_M/T_N=(N/M)^2."""
    N = float(child_frequency)
    M = float(ancestor_frequency)
    if min(N, M) <= 0 or not math.isfinite(N + M):
        raise ValueError("positive finite frequencies required")
    return (N / M) ** 2


def shell_dissipation_law(
    child_frequency: float,
    scaled_lifetime: float,
    shell_energies_u: np.ndarray,
    resolved_energy_fractions: np.ndarray,
    shell_gradient_fractions: np.ndarray,
    time_weights: Sequence[float] | None = None,
    critical_mass_threshold: float | None = None,
) -> dict[str, float]:
    """Finite regression of the continuum dissipation-weighted ancestor law.

    Rows are dyadic shells, columns physical-time cells.  `resolved_energy_fractions`
    is ||P_j V||^2/||P_j u||^2 in [0,1].  `shell_gradient_fractions` is the
    normalized mean |xi|^2/M_j^2 in [1/4,1] for the resolved shell energy.
    """
    N = float(child_frequency)
    c = float(scaled_lifetime)
    E = np.asarray(shell_energies_u, float)
    a = np.asarray(resolved_energy_fractions, float)
    rho = np.asarray(shell_gradient_fractions, float)
    if E.ndim != 2 or a.shape != E.shape or rho.shape != E.shape or E.size == 0:
        raise ValueError("matching nonempty shell/time arrays required")
    if np.any(~np.isfinite(E)) or np.any(~np.isfinite(a)) or np.any(~np.isfinite(rho)):
        raise ValueError("finite shell data required")
    if np.any(E < 0) or np.any(a < 0) or np.any(a > 1) or np.any(rho < 0) or np.any(rho > 1):
        raise ValueError("nonnegative energy, contraction a in [0,1], rho in [0,1] required")
    if N <= 0 or c <= 0:
        raise ValueError("positive scale and lifetime required")
    nt = E.shape[1]
    if time_weights is None:
        wt = np.full(nt, c / (N * N * nt), float)
    else:
        wt = np.asarray(time_weights, float)
        if wt.shape != (nt,) or np.any(wt < 0) or not np.all(np.isfinite(wt)):
            raise ValueError("valid physical time weights required")
        target = c / (N * N)
        if abs(float(wt.sum()) - target) > 2e-12 * max(1.0, target):
            raise ValueError("time weights must sum to c N^-2")

    M = np.array([dyadic_resolved_shell_upper(N, j) for j in range(E.shape[0])], float)[:, None]
    mu = M * E
    density = N * rho * (M * M) * a * E
    D = float(np.sum(density * wt[None, :]))
    threshold = high_strain_ancestor_mass_threshold(c) if critical_mass_threshold is None else float(critical_mass_threshold)
    if threshold < 0 or not math.isfinite(threshold):
        raise ValueError("finite nonnegative threshold required")
    good = mu >= threshold
    Dgood = float(np.sum(np.where(good, density, 0.0) * wt[None, :]))
    Dbad = D - Dgood
    analytic_bad = low_mass_dissipation_upper(N, c, threshold)
    return {
        "total_normalized_dissipation": D,
        "good_ancestor_dissipation": Dgood,
        "bad_low_mass_dissipation": Dbad,
        "analytic_bad_upper": analytic_bad,
        "bad_upper_margin": analytic_bad - Dbad,
        "retained_fraction": Dgood / D if D > 0 else 0.0,
        "critical_mass_threshold": threshold,
        "largest_ancestor_scale_ratio": TRANSPORTER_RADIUS,
        "minimum_ancestor_lifetime_ratio": 1.0 / (TRANSPORTER_RADIUS**2),
    }


def theorem_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    Dstar = clean_high_strain_dissipation_lower(c)
    mustar = high_strain_ancestor_mass_threshold(c)
    bad = low_mass_dissipation_upper(1.0, c, mustar)
    if bad > Dstar / 2.0 + 1e-14 * max(1.0, Dstar):
        raise AssertionError("clean half-dissipation ancestor threshold failed")
    return {
        "status": "EXACT_HIGH_STRAIN_DISSIPATION_WEIGHTED_RESOLVED_ANCESTOR__HALF_DV_ON_CRITICAL_LOW_SHELLS__MATERIAL_RENEWAL_REMAINS",
        "dyadic_geometry": "A_j={M_j/2<|xi|<=M_j}, M_j=(N/4)2^-j, sum_j M_j=N/2",
        "shell_bound": "with |S_(N/4)(xi)|<=1, N||grad P_j V||_2^2 <= N M_j [M_j||P_j u||_2^2]",
        "low_mass_bound": "on mu_j(t)=M_j||P_j u||_2^2<mu_*, total normalized dissipation is <=c mu_*/2",
        "high_strain": f"K>=1/30 gives D_V>=D_*={Dstar:.12g}; choosing mu_*=D_*/c={mustar:.12g} leaves at least half of the actual D_V measure on critical resolved-shell ancestors",
        "causal_law": "the retained law is the actual positive normalized resolved dissipation N|xi|^2|Vhat|^2 dt dxi, not a packet probability or shell argmax",
        "scale": "every retained ancestor shell has M<=N/4, hence its parabolic natural lifetime is at least 16 times the child lifetime",
        "currency": "D_V is still O(1/N) in physical viscous cost and is not promoted to an additive finite reset; the theorem adds a physical low-frequency ancestor mark to its recursive route",
        "scope": "a resolved-shell critical-mass ancestor is not yet a selected coherent transfer parent; material/coherent localization and universal slab renewal of this dissipation-seeded route remain",
    }


@dataclass(frozen=True)
class ResolvedAncestorStress:
    samples: int
    minimum_bad_upper_margin: float
    minimum_half_law_margin: float
    minimum_critical_mass_threshold: float
    minimum_lifetime_ratio: float
    maximum_retained_fraction: float


def stress(samples: int = 50_000, seed: int = 20260809) -> ResolvedAncestorStress:
    rng = np.random.default_rng(seed)
    mbad = mhalf = mmu = mlife = float("inf")
    maxq = 0.0
    for _ in range(samples):
        c = float(math.exp(rng.uniform(-2.0, 2.0)))
        N = float(math.exp(rng.uniform(0.0, 8.0)))
        ns = int(rng.integers(2, 18))
        nt = int(rng.integers(2, 20))
        mu_star = high_strain_ancestor_mass_threshold(c)
        Dstar = clean_high_strain_dissipation_lower(c)

        # Construct arbitrary physical shell states, then scale their energy so
        # total D_V lies on/above the high-strain lower branch.
        E = rng.lognormal(mean=0.0, sigma=2.0, size=(ns, nt))
        a = rng.random((ns, nt))
        rho = 0.25 + 0.75 * rng.random((ns, nt))
        base = shell_dissipation_law(N, c, E, a, rho, critical_mass_threshold=mu_star)
        if base["total_normalized_dissipation"] <= 1e-300:
            continue
        scale = float(rng.uniform(1.0, 5.0)) * Dstar / float(base["total_normalized_dissipation"])
        E *= scale
        out = shell_dissipation_law(N, c, E, a, rho, critical_mass_threshold=mu_star)
        D = float(out["total_normalized_dissipation"])
        if D + 5e-12 * max(1.0, Dstar) < Dstar:
            raise AssertionError("constructed high-strain dissipation state fell below D_*")
        mbad = min(mbad, float(out["bad_upper_margin"]))
        half_margin = float(out["good_ancestor_dissipation"]) - 0.5 * D
        mhalf = min(mhalf, half_margin)
        mmu = min(mmu, mu_star)
        mlife = min(mlife, float(out["minimum_ancestor_lifetime_ratio"]))
        maxq = max(maxq, float(out["retained_fraction"]))
        tol = 8e-12 * max(1.0, D, Dstar)
        if float(out["bad_upper_margin"]) < -tol:
            raise AssertionError("low-critical-mass shells exceeded geometric dissipation upper")
        if half_margin < -tol:
            raise AssertionError("high strain failed to leave half D_V on critical resolved ancestors")

    return ResolvedAncestorStress(samples, mbad, mhalf, mmu, mlife, maxq)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-strain-resolved-ancestor"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "high_strain_resolved_ancestor.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# High strain is a dissipation-weighted resolved-ancestor event\n\nStatus: **{cert['status']}**.\n\nThe normalized dissipation `D_V=N int ||grad V||_2^2 dt` should not be treated as an abstract reset currency.  The strict transporter itself reveals where that dissipation lives.  Decompose its Fourier ball into deterministic dyadic annuli\n\n`A_j={{M_j/2<|xi|<=M_j}},  M_j=(N/4)2^(-j)`.\n\nThey satisfy `sum_j M_j=N/2`.  Put `E_j(t)=||P_j u(t)||_2^2` and the actual critical shell mass `mu_j(t)=M_j E_j(t)`.  Since the standard low-pass multiplier is an L2 contraction and `|xi|<=M_j` on `A_j`,\n\n`N ||grad P_j V||_2^2 <= N M_j mu_j(t)`.\n\nOn a natural lifetime `T=cN^-2`, the part of the **actual resolved dissipation measure** lying where `mu_j(t)<mu_*` is therefore at most\n\n`N T mu_* sum_j M_j = c mu_*/2`.\n\nThe high-strain collision gives\n\n`K>=1/30 => D_V>=D_* = 32 pi^2/(75 c)`.\n\nChoose `mu_*=D_*/c=32 pi^2/(75 c^2)`.  Then the low-mass part is at most `D_*/2`, so for every `D_V>=D_*`,\n\n`D_V({{mu_j>=mu_*}}) >= D_V-D_*/2 >= D_V/2`.\n\nThus at least **half of the actual normalized dissipation law** already carries a simultaneous low-frequency ancestor with fixed critical mass.  No shell is selected by argmax, and there is no logarithmic loss from the infinitely many low shells.  The law is the physical density `N|xi|^2|Vhat|^2 dt dxi` disintegrated by dyadic frequency.\n\nEvery such ancestor has `M_j<=N/4`; its parabolic natural lifetime is therefore at least `16` child lifetimes.  This is a much stronger scale separation than the signed-good HH parent ratio, but it is a different physical object: a resolved reservoir/ancestor state, not yet a transfer-generated hard parent.\n\nStress: `{out.samples}` random shell/time dissipation laws\n- minimum low-mass dissipation upper margin: `{out.minimum_bad_upper_margin:.3e}`\n- minimum half-law margin: `{out.minimum_half_law_margin:.3e}`\n- minimum sampled clean critical-mass threshold: `{out.minimum_critical_mass_threshold:.3e}`\n- minimum ancestor/child lifetime ratio: `{out.minimum_lifetime_ratio:.3f}`\n- maximum retained dissipation fraction: `{out.maximum_retained_fraction:.6f}`\n\nThis does **not** promote `D_V` to a globally finite reset: its physical viscous cost is still `nu D_V/N` and remains summable on geometric high-frequency chains.  What changes is its recursive meaning.  A high-strain stop is now accompanied, on at least half of its own physical dissipation law, by a critical resolved-shell ancestor at a genuinely lower scale.  The remaining bridge is to attach these dissipation-seeded shell ancestors to the existing material/coherent reservoir/reuse or renewed-slab machinery without inventing a packet selector.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
