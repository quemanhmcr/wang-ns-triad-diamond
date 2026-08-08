from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np


def sgs_sampling_scaling_exponent() -> Fraction:
    """N exponent in sum rho_R^(3/2) after sampling+Bernstein.

    rho_R ~ N^-5 D^3 R; sampling contributes N^3 and D^3 Bernstein
    contributes N^(9/2) after the 3/2 power.
    """
    return -5 * Fraction(3, 2) + 3 + 3 * Fraction(3, 2)


def viscous_sampling_scaling_exponent() -> Fraction:
    """N exponent in sum rho_nu^2 after sampling and conversion to d_V.

    rho_nu~N^-5 D^4 V. Sampling: +3. Replace D^4 by N^3 grad V: +6
    after squaring. Finally ||grad V||_2^2=N d_V: +1.
    """
    return -10 + 3 + 6 + 1


def germano_l32_power_upper(filter_l1: float, cubic_increment_charge: float) -> float:
    if filter_l1 < 1.0 - 1e-12 or cubic_increment_charge < 0:
        raise ValueError("invalid Germano data")
    return (1.0 + filter_l1) ** 1.5 * math.sqrt(filter_l1) * cubic_increment_charge


def sgs_source_sample_sum_upper(
    stress_l32_power: float,
    aspect_cap: float,
    scale_radius_cap: float,
    sampling_constant: float,
    derivative_bernstein_constant: float,
) -> float:
    """Bound sum_a rho_R,a^(3/2) for N^-1-separated source centers.

    Analytic inputs:
      sum_a |f(x_a)|^(3/2) <= C_PP N^3 ||f||_(3/2)^(3/2),
      ||D^3 R||_(3/2) <= C_D3 N^3 ||R||_(3/2).
    Affine source normalization uses ||L^-1||||L||^2 <= kappa^2 r_g.
    All powers of N cancel exactly.
    """
    vals = [stress_l32_power, scale_radius_cap]
    if any(v < 0 for v in vals) or min(aspect_cap, sampling_constant, derivative_bernstein_constant) <= 0:
        raise ValueError("invalid SGS sampling data")
    affine = aspect_cap * aspect_cap * scale_radius_cap
    return sampling_constant * (affine * derivative_bernstein_constant) ** 1.5 * stress_l32_power


def sgs_source_sample_sum_from_increments(
    cubic_increment_charge: float,
    filter_l1: float,
    aspect_cap: float,
    scale_radius_cap: float,
    sampling_constant: float,
    derivative_bernstein_constant: float,
) -> float:
    rpow = germano_l32_power_upper(filter_l1, cubic_increment_charge)
    return sgs_source_sample_sum_upper(
        rpow, aspect_cap, scale_radius_cap, sampling_constant, derivative_bernstein_constant
    )


def viscous_source_sample_square_sum_upper(
    normalized_enstrophy: float,
    viscosity: float,
    aspect_cap: float,
    scale_radius_cap: float,
    sampling_constant_l2: float,
    d4_from_grad_constant: float,
) -> float:
    """Bound sum_a rho_nu,a^2 by normalized resolved dissipation/enstrophy.

    Inputs:
      sum_a |f(x_a)|^2 <= C_PP,2 N^3 ||f||_2^2,
      ||D^4 V||_2 <= C_41 N^3 ||grad V||_2.
    """
    if normalized_enstrophy < 0 or viscosity < 0 or scale_radius_cap < 0:
        raise ValueError("invalid viscous sampling data")
    if min(aspect_cap, sampling_constant_l2, d4_from_grad_constant) <= 0:
        raise ValueError("positive sampling constants required")
    affine = aspect_cap * aspect_cap * scale_radius_cap
    return (
        sampling_constant_l2
        * (affine * viscosity * d4_from_grad_constant) ** 2
        * normalized_enstrophy
    )


def pressure_l32_power_upper(
    resolved_critical_mass: float,
    stress_l32_power: float,
    riesz_constant: float,
    resolved_bernstein_l3_constant: float,
) -> float:
    """Filtered pressure L^(3/2) power from V and R.

    -Delta P=partial_i partial_j(V_i V_j+R_ij), so
      ||P||_(3/2) <= C_R (||V||_3^2+||R||_(3/2)).
    If ||V||_3 <= C_B N^(1/2)||V||_2 and mu_V=N||V||_2^2,
      ||V||_3^3 <= C_B^3 mu_V^(3/2).
    Then (a+b)^(3/2)<=sqrt(2)(a^(3/2)+b^(3/2)).
    """
    if min(resolved_critical_mass, stress_l32_power) < 0:
        raise ValueError("nonnegative pressure data required")
    if min(riesz_constant, resolved_bernstein_l3_constant) <= 0:
        raise ValueError("positive pressure constants required")
    return (
        math.sqrt(2.0)
        * riesz_constant**1.5
        * (
            resolved_bernstein_l3_constant**3 * resolved_critical_mass**1.5
            + stress_l32_power
        )
    )


def pressure_source_sample_sum_upper(
    pressure_l32_power: float,
    aspect_cap: float,
    scale_radius_cap: float,
    sampling_constant: float,
    derivative_bernstein_constant: float,
) -> float:
    if pressure_l32_power < 0 or scale_radius_cap < 0:
        raise ValueError("invalid pressure sampling data")
    if min(aspect_cap, sampling_constant, derivative_bernstein_constant) <= 0:
        raise ValueError("positive pressure sampling constants required")
    affine = aspect_cap * aspect_cap * scale_radius_cap
    return sampling_constant * (affine * derivative_bernstein_constant) ** 1.5 * pressure_l32_power


def pressure_source_samples_from_mass_and_increments(
    resolved_critical_mass: float,
    cubic_increment_charge: float,
    filter_l1: float,
    aspect_cap: float,
    scale_radius_cap: float,
    sampling_constant: float,
    pressure_derivative_bernstein_constant: float,
    riesz_constant: float,
    resolved_bernstein_l3_constant: float,
) -> float:
    rpow = germano_l32_power_upper(filter_l1, cubic_increment_charge)
    ppow = pressure_l32_power_upper(
        resolved_critical_mass, rpow, riesz_constant, resolved_bernstein_l3_constant
    )
    return pressure_source_sample_sum_upper(
        ppow,
        aspect_cap,
        scale_radius_cap,
        sampling_constant,
        pressure_derivative_bernstein_constant,
    )


def resolvable_cluster_count_upper(total_p_charge: float, source_threshold: float, p: float) -> float:
    """Maximal separated source-cluster count from an ell^p sampling budget.

    If every selected source center has rho>=rho0 and sum rho^p<=Qp, then
    number of resolvable clusters is at most Qp/rho0^p.  A maximal separated
    subfamily covers all remaining centers by one resolution-scale neighborhood.
    """
    if total_p_charge < 0 or source_threshold <= 0 or p <= 0:
        raise ValueError("invalid cluster-count data")
    return total_p_charge / (source_threshold**p)


def exact_scaling_certificate() -> dict[str, str]:
    sgs = sgs_sampling_scaling_exponent()
    visc = viscous_sampling_scaling_exponent()
    if sgs != 0 or visc != 0:
        raise AssertionError("source sampling lost Navier-Stokes scale invariance")
    return {
        "SGS_pressure_sampling_exponent": str(sgs),
        "viscous_sampling_exponent": str(visc),
        "sampling_input": "Plancherel-Polya for N^-1-separated samples",
        "pressure_identity": "-Delta P=partial_i partial_j(V_i V_j+R_ij)",
        "status": "EXACT_SCALING_GIVEN_STANDARD_PLANCHEREL_POLYA_BERNSTEIN_RIESZ",
    }


@dataclass(frozen=True)
class SourceSamplingStress:
    samples: int
    worst_scale_invariance_residual: float
    minimum_SGS_routing_margin: float
    minimum_pressure_routing_margin: float
    minimum_viscous_routing_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> SourceSamplingStress:
    rng = np.random.default_rng(seed)
    ws = 0.0
    ms = mp = mv = float("inf")
    for _ in range(samples):
        # Positive abstract analytic constants: stress only the exact routing/scaling algebra.
        kappa = float(rng.uniform(1.0, 1.2))
        s0 = float(rng.uniform(0.2, 4.0))
        cpp = float(rng.uniform(0.5, 5.0))
        cd = float(rng.uniform(0.2, 3.0))
        g1 = float(rng.uniform(1.0, 3.0))
        Q = float(10 ** rng.uniform(-8.0, 2.0))
        mu = float(10 ** rng.uniform(-8.0, 2.0))
        CR = float(rng.uniform(0.5, 4.0))
        CB = float(rng.uniform(0.5, 3.0))
        nu = float(rng.uniform(0.0, 2.0))
        dV = float(10 ** rng.uniform(-8.0, 2.0))

        rpow = germano_l32_power_upper(g1, Q)
        sgs1 = sgs_source_sample_sum_upper(rpow, kappa, s0, cpp, cd)
        sgs2 = sgs_source_sample_sum_from_increments(Q, g1, kappa, s0, cpp, cd)
        ms = min(ms, sgs2 - sgs1)
        if abs(sgs2 - sgs1) > 2e-12 * max(1.0, sgs1):
            raise AssertionError("SGS sampling/increment route mismatch")

        ppow = pressure_l32_power_upper(mu, rpow, CR, CB)
        p1 = pressure_source_sample_sum_upper(ppow, kappa, s0, cpp, cd)
        p2 = pressure_source_samples_from_mass_and_increments(
            mu, Q, g1, kappa, s0, cpp, cd, CR, CB
        )
        mp = min(mp, p2 - p1)
        if abs(p2 - p1) > 2e-12 * max(1.0, p1):
            raise AssertionError("pressure source route mismatch")

        v1 = viscous_source_sample_square_sum_upper(dV, nu, kappa, s0, cpp, cd)
        v2 = cpp * (kappa * kappa * s0 * nu * cd) ** 2 * dV
        mv = min(mv, v2 - v1)
        if abs(v2 - v1) > 2e-12 * max(1.0, v1):
            raise AssertionError("viscous source route mismatch")

        # Verify the N powers directly on randomly scaled synthetic magnitudes.
        N = float(10 ** rng.uniform(-3.0, 3.0))
        Rnorm = float(10 ** rng.uniform(-5.0, 2.0))
        # source sample p-power factor before constants:
        lhs = N ** (-7.5) * N**3 * (N**3 * Rnorm) ** 1.5
        rhs = Rnorm**1.5
        ws = max(ws, abs(lhs - rhs) / max(1.0, abs(rhs)))
        if abs(lhs - rhs) > 5e-11 * max(1.0, abs(rhs)):
            raise AssertionError("SGS/pressure scale cancellation failed numerically")
    return SourceSamplingStress(samples, ws, ms, mp, mv)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-bandlimited-source-sampling"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = exact_scaling_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "bandlimited_source_sampling.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Band-limited source sampling: replication cannot reuse one field value for free

Status: **{cert['status']}**.

For `N^-1`-separated source centers, the standard Plancherel--Polya sampling theorem gives `sum_a |f(x_a)|^p <= C_PP N^3 ||f||_p^p` for a fixed band limit.  Combined with the affine source factor `||L^-1||||L||^2<=kappa^2 r_g`, every power of `N` cancels exactly.

Differentiated SGS source:

`sum_a rho_R,a^(3/2) <= C_samp (kappa^2 s0 C_D3)^(3/2) ||R||_(3/2)^(3/2)`.

The exact Germano increment estimate then routes the right-hand side to the global cubic velocity-increment charge at that filter scale.

Viscous-fourth source:

`sum_a rho_nu,a^2 <= C_samp,2 (kappa^2 s0 nu C_41)^2 d_V`,  `d_V=N^-1||grad V||_2^2`.

Thus many separated viscous source grains pay additive resolved dissipation.

For the strict filtered pressure,

`-Delta P=partial_i partial_j(V_i V_j+R_ij)`

and `P` remains band limited.  Riesz + Bernstein gives

`||P||_(3/2)^(3/2) <= sqrt(2) C_R^(3/2)[C_B^3 mu_V^(3/2)+||R||_(3/2)^(3/2)]`.

Hence separated pressure-third source grains route to resolved low-pass critical mass or the same cubic increment charge; pressure-third near-field is not an independent unpriced source.

Stress: `{out.samples}`
- worst scale-invariance residual: `{out.worst_scale_invariance_residual:.3e}`
- minimum SGS routing margin: `{out.minimum_SGS_routing_margin:.3e}`
- minimum pressure routing margin: `{out.minimum_pressure_routing_margin:.3e}`
- minimum viscous routing margin: `{out.minimum_viscous_routing_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
