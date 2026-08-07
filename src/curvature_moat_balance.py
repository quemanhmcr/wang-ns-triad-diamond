from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def localization_error(commutator_coeff: float, curvature_coeff: float, kappa: float, moat_width: float) -> float:
    """a/M + b*kappa*M: filter commutator plus affine-Taylor curvature error."""
    if min(commutator_coeff, curvature_coeff, kappa) < 0 or moat_width <= 0:
        raise ValueError("invalid localization parameters")
    return commutator_coeff / moat_width + curvature_coeff * kappa * moat_width


def optimal_moat_width(commutator_coeff: float, curvature_coeff: float, kappa: float) -> float:
    if commutator_coeff <= 0 or curvature_coeff <= 0 or kappa <= 0:
        raise ValueError("positive coefficients and curvature are required")
    return math.sqrt(commutator_coeff / (curvature_coeff * kappa))


def optimal_localization_error(commutator_coeff: float, curvature_coeff: float, kappa: float) -> float:
    if min(commutator_coeff, curvature_coeff, kappa) < 0:
        raise ValueError("invalid localization parameters")
    if kappa == 0 or commutator_coeff == 0 or curvature_coeff == 0:
        return 0.0
    return 2.0 * math.sqrt(commutator_coeff * curvature_coeff * kappa)


def affine_window_curvature_coeff(lifetime_c: float, deformation_factor: float, grad_chi0: float) -> float:
    """Dimensionless b in b*kappa*M for T=c N^-2 and R=M/N.

    The affine-window bound is (H/2) R Q |grad chi0|.  Integrating for
    T=c N^-2 and writing kappa=N^-3 H gives b=(c/2)Q|grad chi0|.
    """
    if min(lifetime_c, deformation_factor, grad_chi0) < 0:
        raise ValueError("parameters must be nonnegative")
    return 0.5 * lifetime_c * deformation_factor * grad_chi0


def scalar_hessian_bernstein_constant(lambda_cutoff: float = 1.0) -> float:
    """Unitary-Fourier scalar constant for a ball |xi|<=lambda*N.

    Using |xi_i xi_j|<=|xi|^2 and Cauchy-Schwarz,
      N^-3 ||partial_ij f||_inf <= C_B (N||f||_2^2)^(1/2),
    C_B=(2pi)^(-3/2) sqrt(4pi/7) lambda^(7/2).
    Vector/operator Hessian norms can be absorbed by a dimension factor outside
    this scalar certificate.
    """
    if lambda_cutoff <= 0:
        raise ValueError("lambda_cutoff must be positive")
    return (2.0 * math.pi) ** (-1.5) * math.sqrt(4.0 * math.pi / 7.0) * lambda_cutoff ** 3.5


def critical_mass_lower_from_curvature(kappa: float, bernstein_constant: float) -> float:
    """If kappa<=C_B sqrt(mu), then mu>=(kappa/C_B)^2."""
    if kappa < 0 or bernstein_constant <= 0:
        raise ValueError("invalid curvature/mass parameters")
    return (kappa / bernstein_constant) ** 2


def critical_mass_lower_from_unavoidable_error(
    error_threshold: float,
    commutator_coeff: float,
    curvature_coeff: float,
    bernstein_constant: float,
) -> float:
    """Mass forced when even the optimally balanced moat costs >=eta.

    2 sqrt(a b kappa)>=eta => kappa>=eta^2/(4ab), and
    kappa<=C_B sqrt(mu) then gives
      mu>=eta^4/(16 a^2 b^2 C_B^2).
    """
    if error_threshold < 0 or min(commutator_coeff, curvature_coeff, bernstein_constant) <= 0:
        raise ValueError("invalid parameters")
    return error_threshold ** 4 / (
        16.0 * commutator_coeff ** 2 * curvature_coeff ** 2 * bernstein_constant ** 2
    )


def old_quadratic_schedule_partial(n_terms: int) -> tuple[float, float]:
    """Countermodel: M_j=(j+3)^2, kappa_j=(j+3)^-3.

    Sum 1/M converges while sum M*kappa is harmonic and diverges.
    """
    comm = 0.0
    curv = 0.0
    for j in range(n_terms):
        n = j + 3.0
        M = n * n
        kappa = n ** -3
        comm += 1.0 / M
        curv += M * kappa
    return comm, curv


def balanced_schedule_partial(n_terms: int) -> tuple[float, float, float]:
    """For the same kappa_j=(j+3)^-3 use M_j=kappa_j^-1/2.

    Both commutator and curvature terms become (j+3)^(-3/2), summable.
    Returns finite partial sums and a rigorous elementary tail upper bound from
    the integral test: tail after n_terms <=2/sqrt(n_terms+2).
    """
    comm = 0.0
    curv = 0.0
    for j in range(n_terms):
        n = j + 3.0
        kappa = n ** -3
        M = kappa ** -0.5
        comm += 1.0 / M
        curv += M * kappa
    tail_upper = 2.0 / math.sqrt(n_terms + 2.0)
    return comm, curv, tail_upper


def stress(samples: int = 100_000, seed: int = 20260807) -> dict[str, float]:
    import numpy as np

    rng = np.random.default_rng(seed)
    worst_opt_residual = 0.0
    worst_mass_residual = 0.0
    for _ in range(samples):
        a = float(10 ** rng.uniform(-3, 2))
        b = float(10 ** rng.uniform(-3, 2))
        kappa = float(10 ** rng.uniform(-8, 1))
        M = optimal_moat_width(a, b, kappa)
        direct = localization_error(a, b, kappa, M)
        formula = optimal_localization_error(a, b, kappa)
        worst_opt_residual = max(worst_opt_residual, abs(direct - formula) / max(1.0, formula))

        C = float(10 ** rng.uniform(-2, 2))
        eta = float(10 ** rng.uniform(-5, 0))
        mu = critical_mass_lower_from_unavoidable_error(eta, a, b, C)
        kappa_threshold = eta * eta / (4 * a * b)
        direct_mu = critical_mass_lower_from_curvature(kappa_threshold, C)
        worst_mass_residual = max(worst_mass_residual, abs(mu - direct_mu) / max(1.0, direct_mu))

    c_old, h_old = old_quadratic_schedule_partial(100_000)
    c_bal, h_bal, tail = balanced_schedule_partial(100_000)
    return {
        "samples": samples,
        "worst_optimality_residual": worst_opt_residual,
        "worst_mass_identity_residual": worst_mass_residual,
        "old_schedule_commutator_partial_100k": c_old,
        "old_schedule_curvature_partial_100k": h_old,
        "balanced_commutator_partial_100k": c_bal,
        "balanced_curvature_partial_100k": h_bal,
        "balanced_tail_upper_after_100k": tail,
        "scalar_bernstein_constant_lambda1": scalar_hessian_bernstein_constant(1.0),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=100_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-curvature-moat"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = stress(args.samples)
    (args.outdir / "curvature_moat_balance.json").write_text(json.dumps(result, indent=2))
    md = f"""# Curvature-balanced spatial moat\n\n- exact localization law: `E(M)=a/M+b*kappa*M`\n- optimizer: `M*=sqrt(a/(b kappa))`\n- optimum: `E*=2 sqrt(a b kappa)`\n- algebra stress checks: `{result['samples']}`\n- worst optimality residual: `{result['worst_optimality_residual']:.3e}`\n- scalar unitary-Fourier Hessian Bernstein constant at lambda=1: `{result['scalar_bernstein_constant_lambda1']:.12f}`\n- old quadratic schedule, commutator partial sum at 100k: `{result['old_schedule_commutator_partial_100k']:.6f}`\n- old quadratic schedule, curvature partial sum at 100k: `{result['old_schedule_curvature_partial_100k']:.6f}` (harmonic divergence)\n- balanced schedule, each partial sum at 100k: `{result['balanced_commutator_partial_100k']:.6f}`\n- balanced schedule tail upper after 100k: `{result['balanced_tail_upper_after_100k']:.6f}`\n\nThis module records a correction to the previous localization heuristic: an\nexpanding moat cannot be chosen from commutator considerations alone.  The moat\nwidth must balance filter nonlocality against the curvature of the transported\nvelocity field.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
