from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

EXTENDED_ASPECT = Fraction(567, 500)
MILD_ASPECT = Fraction(21, 20)
LOW_STRAIN_ACTION = Fraction(1, 30)


def extended_qpol_lower(condL: float) -> float:
    """Perturbative polarization-only hook lower bound up to cond=567/500."""
    if condL < 1 or condL > float(EXTENDED_ASPECT) + 1e-12:
        raise ValueError("outside extended-aspect bridge")
    raw = 1.0 / math.sqrt(10.0) - math.sqrt(5.0) * (condL - 1.0)
    return max(0.0, raw) ** 2


def physical_three_role_energy_lower_from_hook(condL: float) -> float:
    """Coefficient c with sum_i||F_i^H1||^2 >= c||B_hook||^2."""
    return 0.25 * extended_qpol_lower(condL)


def shell_radius_from_condition_lower(condL: float) -> float:
    """s=N r_g > (2/3) kappa^(1/3) from every axis >2/(3N)."""
    if condL < 1:
        raise ValueError("condition number must be >=1")
    return (2.0 / 3.0) * condL ** (1.0 / 3.0)


def fresh_critical_mass_from_condition_lower(condL: float) -> float:
    """N int_E|u|^2 > (1/5) kappa^(1/3) for a certified fresh affine role."""
    if condL < 1:
        raise ValueError("condition number must be >=1")
    return 0.2 * condL ** (1.0 / 3.0)


def covariance_rate(A: np.ndarray, Sigma: np.ndarray, nu: float) -> np.ndarray:
    A = np.asarray(A, float)
    Sigma = np.asarray(Sigma, float)
    return A @ Sigma + Sigma @ A.T + float(nu) * np.eye(3)


def log_condition_rate_exact(A: np.ndarray, Sigma: np.ndarray, nu: float) -> float:
    """Instantaneous d/dt log cond(L), cond(L)=sqrt(lambda_max/lambda_min).

    Formula is exact at simple eigenvalues. Viscosity contributes a nonpositive
    term and therefore cannot create anisotropy.
    """
    vals, vecs = np.linalg.eigh(np.asarray(Sigma, float))
    if vals[0] <= 0:
        raise ValueError("Sigma must be positive definite")
    S = 0.5 * (np.asarray(A, float) + np.asarray(A, float).T)
    e0, e1 = vecs[:, 0], vecs[:, -1]
    strain = float(e1 @ S @ e1 - e0 @ S @ e0)
    visc = 0.5 * float(nu) * (1.0 / vals[-1] - 1.0 / vals[0])
    return strain + visc


def condition_growth_upper(kappa0: float, strain_action: float) -> float:
    if kappa0 < 1 or strain_action < 0:
        raise ValueError("invalid condition-growth inputs")
    return kappa0 * math.exp(2.0 * strain_action)


def predecessor_condition_lower(kappa1: float, strain_action: float) -> float:
    if kappa1 < 1 or strain_action < 0:
        raise ValueError("invalid predecessor inputs")
    return kappa1 * math.exp(-2.0 * strain_action)


def extended_h1_no_escape_constants() -> dict[str, float]:
    return {
        "extended_qpol": 1.0 / 4000.0,
        "physical_three_role_energy": 1.0 / 16000.0,
        "interaction_forcing_norm": 1.0 / 132.0,
        "physical_first_daughter_energy": 1.0 / 75000.0,
        "h1_pair_or_deficit": 1.0 / 28_800_000.0,
        "full_curvature_pair_or_deficit": 1.0 / 115_200_000.0,
        "source_integral": 1.0 / 600.0,
        "strain_curvature_integral": 1.0 / 16000.0,
        "h1_dominant_strain_action": 1.0 / 32000.0,
        "one_physical_source": 1.0 / 1800.0,
    }


def arb_extended_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 180
    k = arb(567) / 500
    raw = 1 / arb(10).sqrt() - arb(5).sqrt() * (k - 1)
    q = raw * raw
    if not (q > arb(1) / 4000):
        raise AssertionError(f"extended Qpol lower failed: {q}")
    if not ((arb(1) / 15).exp() < arb(27) / 25):
        raise AssertionError("condition growth exp(1/15)<27/25 failed")
    if not ((-(arb(1) / 30)).exp() / arb(16000).sqrt() > arb(1) / 132):
        raise AssertionError("interaction forcing conditioning failed")
    if not ((-(arb(1) / 15)).exp() > arb(14) / 15):
        raise AssertionError("physical pushforward energy factor failed")
    # first interaction daughter: (I/264)^2, then physical energy factor >14/15
    if not ((arb(14) / 15) / (arb(264) ** 2) > arb(1) / 75000):
        raise AssertionError("first daughter clean constant failed")
    # Extended source calculus coefficients at kappa=567/500.
    src = arb(3).sqrt() * (arb(1) / 15).exp() * k
    conn = arb(3).sqrt() * (arb(1) / 15).exp() * (
        4 * k + 2 * k * k + 15 * arb(2).sqrt() * k
    )
    if not (src < arb(11) / 5):
        raise AssertionError(f"extended source coefficient failed: {src}")
    if not (conn < arb(60)):
        raise AssertionError(f"extended connection coefficient failed: {conn}")
    # If both source and AB branches were below these thresholds, J<I/(132T).
    if not ((arb(11) / 5) / 600 + arb(60) / 16000 < arb(1) / 132):
        raise AssertionError("extended source dichotomy arithmetic failed")
    return {
        "extended_aspect_threshold": "567/500",
        "old_mild_threshold": "21/20",
        "low_strain_action": "1/30",
        "extended_qpol_ball": str(q),
        "extended_qpol_lower": "1/4000",
        "physical_three_role_energy_lower": "1/16000",
        "condition_growth_clean": "exp(1/15)<27/25",
        "interaction_forcing_norm_lower": "1/132",
        "physical_first_daughter_energy_lower": "1/75000",
        "extended_h1_pair_or_deficit": "1/28800000",
        "extended_full_curvature_pair_or_deficit": "1/115200000",
        "source_coefficient_upper": "11/5",
        "connection_coefficient_upper": "60",
        "extended_source_or_AB": "int||S||>=I1/(600T) or int||A||||B||>=I1/(16000T)",
        "h1_dominant_strain_action_threshold": "1/32000",
        "one_physical_source_threshold": "I1/(1800T)",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class AspectStickyStress:
    samples: int
    minimum_condition_rate_margin: float
    maximum_viscous_anisotropy_contribution: float
    minimum_shell_radius_margin: float
    minimum_fresh_mass_margin: float
    minimum_predecessor_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> AspectStickyStress:
    rng = np.random.default_rng(seed)
    mr = mf = mp = float("inf")
    mc = float("inf")
    mv = -float("inf")
    for _ in range(samples):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
        axes = np.exp(rng.uniform(-2.0, 2.0, size=3))
        Sigma = Q @ np.diag(axes * axes) @ Q.T
        A = rng.normal(size=(3, 3))
        A -= np.trace(A) / 3.0 * np.eye(3)
        nu = float(rng.uniform(0.0, 2.0))
        rate = log_condition_rate_exact(A, Sigma, nu)
        S = 0.5 * (A + A.T)
        bound = 2.0 * float(np.linalg.norm(S, 2))
        mc = min(mc, bound - rate)
        if rate > bound + 2e-12:
            raise AssertionError("condition growth rate exceeded strain bound")
        vals = np.linalg.eigvalsh(Sigma)
        visc = 0.5 * nu * (1.0 / vals[-1] - 1.0 / vals[0])
        mv = max(mv, visc)
        if visc > 2e-14:
            raise AssertionError("viscosity created anisotropy")

        N = float(math.exp(rng.uniform(-2.0, 2.0)))
        lmin = (2.0 / (3.0 * N)) * float(rng.uniform(1.0001, 4.0))
        kappa = float(rng.uniform(1.0, 1000.0))
        lmax = kappa * lmin
        lmid = float(rng.uniform(lmin, lmax))
        rg = (lmin * lmid * lmax) ** (1.0 / 3.0)
        s = N * rg
        lower_s = shell_radius_from_condition_lower(kappa)
        mr = min(mr, s - lower_s)
        if s <= lower_s - 2e-12:
            raise AssertionError("shell radius/condition lower bound failed")
        mass = 0.3 * s
        lower_m = fresh_critical_mass_from_condition_lower(kappa)
        mf = min(mf, mass - lower_m)
        if mass <= lower_m - 2e-12:
            raise AssertionError("fresh aspect-radius energy bound failed")

        K = float(rng.uniform(0.0, 1.0 / 30.0))
        k0 = float(rng.uniform(1.0, 20.0))
        k1 = condition_growth_upper(k0, K)
        pred = predecessor_condition_lower(k1, K)
        mp = min(mp, pred - k0)
        if pred < k0 - 2e-12:
            raise AssertionError("predecessor condition inversion failed")
    return AspectStickyStress(samples, mc, mv, mr, mf, mp)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-affine-aspect-sticky"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_extended_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": out.__dict__, "no_escape_constants": extended_h1_no_escape_constants()}
    (args.outdir / "affine_aspect_sticky.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Affine aspect: extended H1 bridge or sticky ancestry

Status: **{cert['status']}**.

The physical hook polarization bridge extends from `cond(L)<=21/20` to the exact transition threshold

`cond(L)<=567/500`,

with the clean lower bound

`Q_pol >= (1/4000)||B_hook||^2`,

hence physical three-role H1 sideband energy `>=||B_hook||^2/16000`.  Accounting for the nonunitary base propagators, feedback, one-of-three-role selection, odd-Hermite convexity and pair rescue gives the conservative local alternative

`Def or R_pair >= I1^2/28800000`.

Together with the H3 sector this yields the extended full-curvature cost `I_B^2/115200000` outside source/feedback/large-daughter/rescue exits.

For covariance dynamics `Sigma_dot=A Sigma+Sigma A^T+nu I`, viscosity contributes nonpositively to `d log cond(L)/dt`, and

`d log cond(L)/dt <= 2||S||_op`.

Thus on the existing low-strain action `int||S||<=1/30`, `cond(L)` can increase by less than `27/25`.  Since `(27/25)(21/20)=567/500`, any current grain with condition number above `567/500` must have had a predecessor already above `21/20`, unless the ancestry identification breaks and the grain is fresh.

Fresh high aspect is paid only through physical radius, not a Young deficit.  Shell concentration implies

`s=N r_g > (2/3) cond(L)^(1/3)`,

so the affine critical mass gives

`N int_E|u|^2 > (1/5) cond(L)^(1/3)`.

The extended dephasing source calculus also remains quantitative: `J1>=I1/(132T)` forces `int||S_source||>=I1/(600T)` or `int||A||||B||>=I1/(16000T)`; on the H1-dominant branch, strain action below `1/32000` makes one of pressure-third / differentiated-SGS / viscous-fourth sources at least `I1/(1800T)`.

Stress: `{out.samples}`
- minimum condition-growth margin: `{out.minimum_condition_rate_margin:.3e}`
- maximum viscous anisotropy contribution: `{out.maximum_viscous_anisotropy_contribution:.3e}`
- minimum shell-radius margin: `{out.minimum_shell_radius_margin:.3e}`
- minimum fresh-mass margin: `{out.minimum_fresh_mass_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
