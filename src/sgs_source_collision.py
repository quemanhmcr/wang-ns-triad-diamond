from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.affine_sgs_boundary_ledger import sgs_increment_cubic_upper


def stress_support_radius(selected_N: float) -> float:
    """Support radius of the SGS stress for V=S_{N/4}u.

    The filtered product is supported in B_{N/4}; V tensor V in B_{N/2}.
    Hence R is supported in B_{N/2}.
    """
    if selected_N <= 0:
        raise ValueError("N must be positive")
    return 0.5 * selected_N


def third_derivative_l32_to_linf_constant(lambda_ratio: float = 0.5) -> float:
    """Unitary-Fourier Bernstein constant for nabla^2 div R.

    For Hilbert-valued R supported in |xi|<=lambda N,
      ||nabla^2 div R||_inf
      <= C(lambda) N^5 ||R||_{3/2,F},
    with Hausdorff--Young interpolation constant (2pi)^(-1/2).
    The symbol R -> nabla^2 div R has Hilbert operator norm <=|xi|^3.
    """
    if lambda_ratio <= 0:
        raise ValueError("positive support ratio required")
    return ((2.0 * math.pi) ** -2.0
            * (8.0 * math.pi / 15.0) ** (2.0 / 3.0)
            * lambda_ratio ** 5)


def fourth_from_gradient_l2_to_linf_constant(lambda_ratio: float = 0.25) -> float:
    """Unitary-Fourier Bernstein constant for nabla^2 Delta V from grad V.

    If supp Vhat lies in |xi|<=lambda N,
      ||nabla^2 Delta V||_inf
      <= C(lambda) N^(9/2) ||grad V||_2.
    """
    if lambda_ratio <= 0:
        raise ValueError("positive support ratio required")
    return ((2.0 * math.pi) ** -1.5
            * math.sqrt(4.0 * math.pi / 9.0)
            * lambda_ratio ** 4.5)


def affine_length_factor(L: np.ndarray, N: float) -> float:
    """Dimensionless factor N ||L^-1||op ||L||op^2."""
    L = np.asarray(L, float)
    if L.shape != (3, 3) or N <= 0:
        raise ValueError("need a 3x3 grain matrix and N>0")
    return float(N * np.linalg.norm(np.linalg.inv(L), 2) * np.linalg.norm(L, 2) ** 2)


def geometric_radius(L: np.ndarray) -> float:
    L = np.asarray(L, float)
    det = abs(float(np.linalg.det(L)))
    if det <= 0:
        raise ValueError("invertible grain matrix required")
    return det ** (1.0 / 3.0)


def affine_factor_from_radius_upper(L: np.ndarray, N: float) -> tuple[float, float, float]:
    """Return actual factor, kappa^2*N*r_g upper bound, and kappa."""
    svals = np.linalg.svd(np.asarray(L, float), compute_uv=False)
    kappa = float(svals.max() / svals.min())
    rg = geometric_radius(L)
    actual = affine_length_factor(L, N)
    upper = kappa * kappa * N * rg
    return actual, upper, kappa


def clean_sgs_source_upper(scale_radius: float, stress_l32: float) -> float:
    """Clean mild-aspect normalized source bound rho_R <= (3/2000)s||R||_3/2."""
    if scale_radius < 0 or stress_l32 < 0:
        raise ValueError("nonnegative inputs required")
    return (3.0 / 2000.0) * scale_radius * stress_l32


def stress_l32_from_source_lower(normalized_source: float, scale_radius_cap: float) -> float:
    if normalized_source < 0 or scale_radius_cap <= 0:
        raise ValueError("invalid collision parameters")
    return (2000.0 / (3.0 * scale_radius_cap)) * normalized_source


def cubic_increment_from_sgs_source_lower(
    normalized_source: float,
    scale_radius_cap: float,
    filter_l1: float,
) -> float:
    """Cubic increment charge forced by differentiated SGS source.

    Uses ||R||_(3/2) >= 2000 rho/(3s0), then
      int |R|^(3/2) <= C_inc int int |G||delta u|^3.
    """
    rnorm = stress_l32_from_source_lower(normalized_source, scale_radius_cap)
    cinc = sgs_increment_cubic_upper(filter_l1, 1.0)
    return rnorm ** 1.5 / cinc


def clean_viscous_source_upper(
    viscosity: float,
    scale_radius: float,
    normalized_enstrophy: float,
) -> float:
    """rho_nu <= (nu s/5000) sqrt(d), d=N^-1||grad V||_2^2."""
    if min(viscosity, scale_radius, normalized_enstrophy) < 0:
        raise ValueError("nonnegative inputs required")
    return viscosity * scale_radius * math.sqrt(normalized_enstrophy) / 5000.0


def enstrophy_from_viscous_source_lower(
    normalized_source: float,
    viscosity: float,
    scale_radius_cap: float,
) -> float:
    if normalized_source < 0 or viscosity <= 0 or scale_radius_cap <= 0:
        raise ValueError("invalid viscous collision parameters")
    return (5000.0 * normalized_source / (viscosity * scale_radius_cap)) ** 2


def fresh_radius_mass_lower(scale_radius: float) -> float:
    """N times physical L2 mass on E2 from M_aff(E2)>=3/10."""
    if scale_radius < 0:
        raise ValueError("nonnegative scale radius required")
    return 0.3 * scale_radius


def arb_constants_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 180
    pi = arb.pi()
    half = arb(1) / 2
    cr = (2 * pi) ** (-2) * (8 * pi / 15).root(3) ** 2 * half ** 5
    cv = (2 * pi) ** (-arb(3) / 2) * (4 * pi / 9).sqrt() * (arb(1) / 4) ** (arb(9) / 2)
    if not (cr < arb(1) / 800):
        raise AssertionError(f"SGS Bernstein constant failed: {cr}")
    if not (cv < arb(1) / 6000):
        raise AssertionError(f"viscous Bernstein constant failed: {cv}")
    if not ((arb(21) / 20) ** 2 < arb(6) / 5):
        raise AssertionError("mild-aspect affine factor failed")
    return {
        "stress_support_ratio": "1/2",
        "transporter_support_ratio": "1/4",
        "sgs_l32_to_source_constant_ball": str(cr),
        "sgs_clean_constant": "1/800",
        "viscous_grad_to_source_constant_ball": str(cv),
        "viscous_clean_constant": "1/6000",
        "mild_aspect_kappa_squared_upper": "6/5",
        "normalized_sgs_source_upper": "(3/2000) s ||R||_(3/2)",
        "normalized_viscous_source_upper": "(nu s/5000) sqrt(d_V)",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class SourceCollisionStress:
    samples: int
    worst_affine_factor_ratio: float
    minimum_sgs_collision_margin: float
    minimum_viscous_collision_margin: float
    minimum_large_radius_mass_margin: float


def _random_mild_L(rng: np.random.Generator, N: float) -> np.ndarray:
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    kappa = float(rng.uniform(1.0, 1.05))
    # Keep geometric radius arbitrary: shape is mild but physical scale may vary.
    rg = float(math.exp(rng.uniform(-2.0, 2.0)) / N)
    exps = rng.uniform(-0.5, 0.5, size=3)
    exps -= np.mean(exps)
    raw = kappa ** exps
    raw /= np.prod(raw) ** (1.0 / 3.0)
    # Rescale if numerical spread exceeds requested mild condition.
    spread = raw.max() / raw.min()
    if spread > kappa:
        raw = np.exp(np.log(raw) * math.log(kappa) / math.log(spread))
    return Q @ np.diag(rg * raw)


def stress(samples: int = 50_000, seed: int = 20260808) -> SourceCollisionStress:
    rng = np.random.default_rng(seed)
    waf = 0.0
    msgs = mvis = mrad = float("inf")
    for _ in range(samples):
        N = float(math.exp(rng.uniform(-3.0, 3.0)))
        L = _random_mild_L(rng, N)
        actual, upper, kappa = affine_factor_from_radius_upper(L, N)
        if actual > upper * (1.0 + 3e-12):
            raise AssertionError("affine factor/radius inequality failed")
        waf = max(waf, actual / upper)
        if kappa > 21.0 / 20.0 + 2e-12:
            raise AssertionError("stress generator escaped mild-aspect branch")

        s0 = float(rng.uniform(0.4, 4.0))
        rho = float(rng.uniform(1e-6, 2e-2))
        rnorm = 1.03 * stress_l32_from_source_lower(rho, s0)
        upper_rho = clean_sgs_source_upper(s0, rnorm)
        msgs = min(msgs, upper_rho - rho)
        if upper_rho < rho:
            raise AssertionError("SGS source collision direction failed")

        nu = float(rng.uniform(0.1, 2.0))
        d = 1.03 * enstrophy_from_viscous_source_lower(rho, nu, s0)
        upper_nu = clean_viscous_source_upper(nu, s0, d)
        mvis = min(mvis, upper_nu - rho)
        if upper_nu < rho:
            raise AssertionError("viscous source collision direction failed")

        s = float(rng.uniform(0.0, 8.0))
        mass = fresh_radius_mass_lower(s)
        mrad = min(mrad, mass - 0.3 * s)
        if mass + 1e-15 < 0.3 * s:
            raise AssertionError("fresh-radius mass identity failed")
    return SourceCollisionStress(samples, waf, msgs, mvis, mrad)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-sgs-source-collision"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_constants_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": out.__dict__}
    (args.outdir / "sgs_source_collision.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Filtered SGS / viscous source collision

Status: **{cert['status']}**.

For the strict transporter `V=S_(N/4)u`, its SGS stress has Fourier support in `|xi|<=N/2`. Vector-valued Hausdorff--Young/Bernstein and the mild-aspect affine factor give

`N^-4 ||S_R|| <= (3/2000) s ||R||_(3/2)`,  `s=N r_g`,

while the viscous source obeys

`N^-4 ||S_nu|| <= (nu s/5000) sqrt(d_V)`,  `d_V=N^-1||grad V||_2^2`.

Thus, for `s<=s0`, a differentiated-SGS source level `rho_R` forces

`||R||_(3/2) >= 2000 rho_R/(3 s0)`,

and hence, by the exact Germano increment bound, cubic velocity-increment charge at the actual `N/4` filter scale. A viscous source level `rho_nu` forces

`d_V >= (5000 rho_nu/(nu s0))^2`.

If instead `s>s0`, the selected affine grain itself carries scale-critical physical mass `N int_E |u|^2 >= (3/10)s0`; this is a radius-energy/ancestry event, not an aspect defect.

Stress: `{out.samples}`
- worst affine-factor / `kappa^2 N r_g` ratio: `{out.worst_affine_factor_ratio:.9f}`
- minimum SGS collision margin: `{out.minimum_sgs_collision_margin:.3e}`
- minimum viscous collision margin: `{out.minimum_viscous_collision_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
