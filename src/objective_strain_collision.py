from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

PRESSURE_HESSIAN_KERNEL_CONSTANT = 150.0


def derivative_bernstein_constant(order: int, lambda_cutoff: float = 1.0) -> float:
    """Scalar unitary-Fourier L2->Linf derivative constant in R^3.

    For supp fhat in |xi|<=lambda*N and critical mass mu=N||f||_2^2,
      N^{-(order+1)} ||partial^order f||_inf <= C_order sqrt(mu),
    using |monomial|<=|xi|^order and Cauchy-Schwarz.
    """
    if order < 0 or lambda_cutoff <= 0:
        raise ValueError("invalid derivative Bernstein parameters")
    return (
        (2.0 * math.pi) ** (-1.5)
        * math.sqrt(4.0 * math.pi / (2 * order + 3))
        * lambda_cutoff ** (order + 1.5)
    )


def quadratic_objective_source_mass_lower(source_level: float, gradient_constant: float) -> float:
    """Mass forced by Q=-S^2-Omega^2+[S,Omega].

    ||Q|| <=4||A||^2.  If N^-2||A||<=C_A sqrt(mu), then
      N^-4||Q||<=4 C_A^2 mu.
    """
    if source_level < 0 or gradient_constant <= 0:
        raise ValueError("invalid source parameters")
    return source_level / (4.0 * gradient_constant * gradient_constant)


def viscous_objective_source_mass_lower(source_level: float, nu: float, third_derivative_constant: float) -> float:
    """Mass forced by nu Delta S in a band-limited packet model."""
    if source_level < 0 or nu < 0 or third_derivative_constant <= 0:
        raise ValueError("invalid viscous source parameters")
    if source_level == 0:
        return 0.0
    if nu == 0:
        return math.inf
    return (source_level / (nu * third_derivative_constant)) ** 2


def pressure_hessian_kernel_component_bound() -> float:
    """Clean elementary far-field bound for Hessian of the pressure kernel.

    K_ij=(4pi)^-1(3 z_i z_j r^-5-delta_ij r^-3).  A direct product-rule
    estimate bounds each scalar second derivative by 204/(4pi) r^-5.  Taking
    the 3x3 Hessian Frobenius norm and sum_ij |v_i v_j|<=3|v|^2 gives
    9*204/(4pi)=146.1... <150.
    """
    raw = 9.0 * 204.0 / (4.0 * math.pi)
    if raw >= PRESSURE_HESSIAN_KERNEL_CONSTANT:
        raise AssertionError("clean pressure-Hessian constant no longer dominates")
    return PRESSURE_HESSIAN_KERNEL_CONSTANT


def far_pressure_hessian_no_fresh_coefficient(
    packing_constant: float,
    first_shell: int = 3,
    kernel_constant: float = PRESSURE_HESSIAN_KERNEL_CONSTANT,
) -> float:
    """C_far in N^-4 |Hess p_far| <= C_far mu_max.

    A packet of critical mass mu at distance 2^n/N contributes O(mu 2^-5n).
    Three-dimensional packing gives at most C_geom 2^{3n} packets, leaving
    sum 2^-2n=(4/3)4^-n0.
    """
    if packing_constant <= 0 or first_shell < 0 or kernel_constant <= 0:
        raise ValueError("invalid pressure Hessian packing parameters")
    return kernel_constant * packing_constant * (4.0 / 3.0) * (4.0 ** (-first_shell))


def pressure_hessian_fresh_mass_threshold(
    source_level: float,
    near_coefficient: float,
    far_coefficient: float,
) -> float:
    """Fresh/local critical mass forced by a pressure-Hessian source level.

    If absence of mass >mu implies H_near<=C_near mu and H_far<=C_far mu,
    then source>=rho contradicts that absence once both channels are below
    half the budget.
    """
    if source_level <= 0 or near_coefficient < 0 or far_coefficient < 0:
        raise ValueError("invalid pressure Hessian collision parameters")
    candidates = []
    if near_coefficient > 0:
        candidates.append(source_level / (2.0 * near_coefficient))
    if far_coefficient > 0:
        candidates.append(source_level / (2.0 * far_coefficient))
    if not candidates:
        raise ValueError("at least one pressure channel is required")
    return min(candidates)


def objective_source_channel_level(strain_number: float, lifetime_c: float) -> float:
    """Per-channel normalized source forced by 5% coherence failure.

    Let d=sigma N^2 and T=c N^-2.  Coherence failure gives
      int ||S_circ|| dt >= d/20.
    Since S_circ=Q-Hess p+nu Delta S, one of three channels has average
    N^-4 magnitude at least sigma/(60 c).
    """
    if strain_number < 0 or lifetime_c <= 0:
        raise ValueError("invalid strain/lifetime parameters")
    return strain_number / (60.0 * lifetime_c)


def source_collision_thresholds(
    strain_number: float,
    lifetime_c: float,
    gradient_constant: float,
    nu: float,
    third_derivative_constant: float,
    pressure_near_coefficient: float,
    pressure_far_coefficient: float,
) -> dict[str, float]:
    rho = objective_source_channel_level(strain_number, lifetime_c)
    return {
        "channel_source_level": rho,
        "quadratic_mass_lower": quadratic_objective_source_mass_lower(rho, gradient_constant),
        "viscous_mass_lower": viscous_objective_source_mass_lower(rho, nu, third_derivative_constant),
        "pressure_mass_lower": pressure_hessian_fresh_mass_threshold(rho, pressure_near_coefficient, pressure_far_coefficient),
    }


def stress(samples: int = 50_000, seed: int = 20260807) -> dict[str, float]:
    import numpy as np

    rng = np.random.default_rng(seed)
    worst_quad_residual = 0.0
    worst_visc_residual = 0.0
    worst_pressure_margin = float("inf")
    kernel = pressure_hessian_kernel_component_bound()
    C1 = derivative_bernstein_constant(1, 1.0)
    C3 = derivative_bernstein_constant(3, 1.0)
    cfar = far_pressure_hessian_no_fresh_coefficient(3.0, 3)

    for _ in range(samples):
        rho = float(10 ** rng.uniform(-6, -0.2))
        CA = float(10 ** rng.uniform(-2, 1))
        muq = quadratic_objective_source_mass_lower(rho, CA)
        directq = 4 * CA * CA * muq
        worst_quad_residual = max(worst_quad_residual, abs(directq - rho) / max(1.0, rho))

        nu = float(10 ** rng.uniform(-2, 1))
        C = float(10 ** rng.uniform(-2, 1))
        muv = viscous_objective_source_mass_lower(rho, nu, C)
        directv = nu * C * math.sqrt(muv)
        worst_visc_residual = max(worst_visc_residual, abs(directv - rho) / max(1.0, rho))

        cnear = float(10 ** rng.uniform(-1, 2))
        muf = 0.9 * pressure_hessian_fresh_mass_threshold(rho, cnear, cfar)
        upper = cnear * muf + cfar * muf
        # The split-half certificate implies each term <rho/2, hence total<rho.
        worst_pressure_margin = min(worst_pressure_margin, rho - upper)
        if upper >= rho + 1e-14:
            raise AssertionError("pressure Hessian fresh threshold failed")

    return {
        "samples": samples,
        "pressure_hessian_kernel_constant": kernel,
        "gradient_bernstein_constant_lambda1": C1,
        "third_derivative_bernstein_constant_lambda1": C3,
        "far_pressure_hessian_coefficient_example": cfar,
        "worst_quadratic_identity_residual": worst_quad_residual,
        "worst_viscous_identity_residual": worst_visc_residual,
        "minimum_pressure_collision_margin": worst_pressure_margin,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-objective-strain"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = stress(args.samples)
    (args.outdir / "objective_strain_collision.json").write_text(json.dumps(result, indent=2))
    md = f"""# Objective-strain source collision\n\n- exact source split: `S_circ = Q(A) - Hess p + nu Delta S`\n- clean far pressure-Hessian kernel constant: `< {result['pressure_hessian_kernel_constant']:.0f}`\n- far shell exponent after 3D packing: `5-3=2`\n- scalar gradient Bernstein constant (lambda=1): `{result['gradient_bernstein_constant_lambda1']:.12f}`\n- scalar third-derivative Bernstein constant (lambda=1): `{result['third_derivative_bernstein_constant_lambda1']:.12f}`\n- packet collision stress cases: `{result['samples']}`\n- worst quadratic threshold residual: `{result['worst_quadratic_identity_residual']:.3e}`\n- worst viscous threshold residual: `{result['worst_viscous_identity_residual']:.3e}`\n- minimum pressure-Hessian collision margin: `{result['minimum_pressure_collision_margin']:.3e}`\n\nA 5% strain-coherence failure over `T=cN^-2`, with initial non-conformal strain\n`d=sigma N^2`, forces at least one normalized source channel at level\n`sigma/(60c)`.  In the stated band-limited packet model the quadratic and\nviscous channels force critical mass directly; the pressure-Hessian far field\nhas a stronger `2^-2n` summable packing gain, while its near field is passed to\nthe local critical-mass coefficient.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
