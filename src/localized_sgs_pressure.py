from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ResolvedWindowStep:
    """Integrated localized resolved-energy bookkeeping on one time interval.

    Exact PDE identity (with the definitions in the accompanying note):
        combined_work + dissipation = energy_in - energy_out + leakage.
    """
    energy_in: float
    energy_out: float
    combined_work: float
    dissipation: float
    leakage: float

    @property
    def residual(self) -> float:
        return self.combined_work + self.dissipation - (
            self.energy_in - self.energy_out + self.leakage
        )


def infer_combined_work(energy_in: float, energy_out: float, dissipation: float, leakage: float) -> float:
    return energy_in - energy_out + leakage - dissipation


def pressure_cancellation_trichotomy(raw_sgs: float, combined_work: float) -> dict[str, float | str]:
    """Exact algebraic dichotomy for positive raw SGS transfer.

    Write combined_work = raw_sgs + pressure_divergence_work.  If raw_sgs>0,
    either combined work retains at least half, or the pressure boundary term
    opposes it with magnitude at least half.
    """
    if raw_sgs < 0.0:
        raise ValueError("the trichotomy is stated for nonnegative raw SGS transfer")
    pressure_work = combined_work - raw_sgs
    if combined_work >= 0.5 * raw_sgs:
        branch = "combined_work"
        retained = combined_work
        pressure_cancel = max(0.0, -pressure_work)
    else:
        branch = "pressure_boundary"
        retained = max(0.0, combined_work)
        pressure_cancel = -pressure_work
    return {
        "branch": branch,
        "raw_sgs": float(raw_sgs),
        "combined_work": float(combined_work),
        "pressure_divergence_work": float(pressure_work),
        "combined_retained": float(retained),
        "pressure_cancellation": float(pressure_cancel),
        "half_flux": float(0.5 * raw_sgs),
    }


def critical_boundary_charge_lower_bound(raw_sgs: float, grad_chi_sup: float) -> float:
    """CKN-type annular charge forced by the pressure-cancellation branch.

    If |int p u.grad chi| >= raw_sgs/2, Holder and Young give
      raw_sgs/(2 ||grad chi||_inf)
        <= ||p||_{3/2} ||u||_3
        <= (2/3) int |p|^{3/2} + (1/3) int |u|^3
        <= int (|p|^{3/2}+|u|^3).
    """
    if raw_sgs < 0.0 or grad_chi_sup <= 0.0:
        raise ValueError("require nonnegative raw flux and positive cutoff gradient")
    return raw_sgs / (2.0 * grad_chi_sup)


def weighted_chain_identity(steps: list[ResolvedWindowStep], weights: list[float]) -> dict[str, float]:
    """Summation-by-parts identity for a chain with matching endpoint energies.

    Requires E_out[j]=E_in[j+1].  Then
      sum a_j(W_j+D_j-L_j)
       = a_0 E_0 - a_{n-1}E_n + sum_{j=1}^{n-1}(a_j-a_{j-1})E_j.
    For nonincreasing nonnegative weights and E_j>=0, this is <=a_0 E_0.
    """
    if len(steps) == 0 or len(steps) != len(weights):
        raise ValueError("matching nonempty step/weight lists required")
    if any(a < 0.0 for a in weights):
        raise ValueError("weights must be nonnegative")
    for j in range(len(steps) - 1):
        scale = max(1.0, abs(steps[j].energy_out), abs(steps[j + 1].energy_in))
        if abs(steps[j].energy_out - steps[j + 1].energy_in) > 1e-11 * scale:
            raise ValueError("chain endpoint energies do not match")
    lhs = sum(a * (s.combined_work + s.dissipation - s.leakage) for a, s in zip(weights, steps))
    rhs = weights[0] * steps[0].energy_in - weights[-1] * steps[-1].energy_out
    for j in range(1, len(steps)):
        rhs += (weights[j] - weights[j - 1]) * steps[j].energy_in
    return {
        "lhs": float(lhs),
        "rhs": float(rhs),
        "residual": float(lhs - rhs),
        "initial_budget": float(weights[0] * steps[0].energy_in),
        "nonincreasing": bool(all(weights[j] <= weights[j - 1] + 1e-15 for j in range(1, len(weights)))),
    }


def positive_work_depletion_bound(steps: list[ResolvedWindowStep], weights: list[float]) -> dict[str, float]:
    """Budget positive combined work by initial energy, leakage and backscatter.

    For nonincreasing weights and exact steps:
      sum a W_+ + sum a D
        <= a_0 E_0 + sum a L_+ + sum a W_-.
    """
    cert = weighted_chain_identity(steps, weights)
    if not cert["nonincreasing"]:
        raise ValueError("depletion bound requires nonincreasing weights")
    positive_work = sum(a * max(s.combined_work, 0.0) for a, s in zip(weights, steps))
    backscatter = sum(a * max(-s.combined_work, 0.0) for a, s in zip(weights, steps))
    diss = sum(a * s.dissipation for a, s in zip(weights, steps))
    positive_leak = sum(a * max(s.leakage, 0.0) for a, s in zip(weights, steps))
    rhs = weights[0] * steps[0].energy_in + positive_leak + backscatter
    return {
        "positive_work": float(positive_work),
        "dissipation": float(diss),
        "backscatter": float(backscatter),
        "positive_leakage": float(positive_leak),
        "rhs_budget": float(rhs),
        "margin": float(rhs - positive_work - diss),
    }


def random_chain_stress(samples: int = 50000, seed: int = 20260807) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    worst_identity = 0.0
    worst_depletion = float("inf")
    pressure_branch_checks = 0
    for _ in range(samples):
        n = int(rng.integers(1, 15))
        energies = rng.lognormal(mean=0.0, sigma=1.0, size=n + 1)
        weights = np.sort(rng.uniform(0.05, 1.0, size=n))[::-1]
        steps: list[ResolvedWindowStep] = []
        for j in range(n):
            diss = float(rng.uniform(0.0, 0.4))
            leak = float(rng.normal(0.0, 0.25))
            work = infer_combined_work(float(energies[j]), float(energies[j + 1]), diss, leak)
            steps.append(ResolvedWindowStep(float(energies[j]), float(energies[j + 1]), work, diss, leak))
        cert = weighted_chain_identity(steps, weights.tolist())
        worst_identity = max(worst_identity, abs(cert["residual"]))
        dep = positive_work_depletion_bound(steps, weights.tolist())
        worst_depletion = min(worst_depletion, dep["margin"])

        raw = float(rng.lognormal(mean=-1.0, sigma=1.0))
        combined = float(rng.normal(0.4 * raw, raw))
        tri = pressure_cancellation_trichotomy(raw, combined)
        if tri["branch"] == "combined_work":
            if float(tri["combined_work"]) + 1e-14 < 0.5 * raw:
                raise AssertionError("combined-work branch failed")
        else:
            pressure_branch_checks += 1
            if float(tri["pressure_cancellation"]) + 1e-14 < 0.5 * raw:
                raise AssertionError("pressure branch failed")
    return {
        "samples": samples,
        "worst_chain_identity_error": float(worst_identity),
        "worst_depletion_margin": float(worst_depletion),
        "pressure_branch_checks": int(pressure_branch_checks),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50000)
    ap.add_argument("--outdir", type=Path, default=Path("results-localized-sgs-pressure"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = random_chain_stress(args.samples)
    (args.outdir / "localized_sgs_pressure.json").write_text(json.dumps(out, indent=2))
    md = f'''# Localized SGS / pressure-work ledger

The identities are analytic; random traces only validate implementation.

- random finite chains: `{out['samples']}`
- worst weighted summation-by-parts residual: `{out['worst_chain_identity_error']:.3e}`
- worst positive-work depletion margin: `{out['worst_depletion_margin']:.3e}`
- pressure-cancellation branches checked: `{out['pressure_branch_checks']}`

For a positive raw SGS flux `S`, either combined work retains at least `S/2`,
or pressure boundary work has magnitude at least `S/2`. In the latter branch,

`int_A (|u|^3+|p|^(3/2)) >= S/(2 ||grad chi||_inf)`.
'''
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
