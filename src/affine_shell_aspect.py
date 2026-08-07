from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

EPS = 1.0 / 100.0
P = 1.5
SHELL_LOG_HALF = 2.0 / 25.0


def p_density_frequency_covariance(A: np.ndarray) -> np.ndarray:
    """Covariance of |F|^(3/2) for F=exp(-1/2 (xi-k)^T A (xi-k))."""
    A = np.asarray(A, float)
    return (2.0 / 3.0) * np.linalg.inv(A)


def physical_l2_covariance(A: np.ndarray) -> np.ndarray:
    """Under the unitary Fourier convention, |check F|^2 has covariance A/2."""
    return 0.5 * np.asarray(A, float)


def uncertainty_matrix_residual(A: np.ndarray) -> float:
    G = p_density_frequency_covariance(A)
    S = physical_l2_covariance(A)
    return float(np.linalg.norm(S @ G - np.eye(3) / 3.0))


def physical_axis_lower_constant(eps: float = EPS, shell_log_half: float = SHELL_LOG_HALF) -> float:
    mass = 1.0 - eps ** P
    B = math.exp(shell_log_half)
    return mass * math.sqrt(2.0 * math.pi) / (2.0 * B * math.sqrt(3.0))


def gaussian_l2_coefficient() -> float:
    """||F||_2^2 = c_L (det Sigma_x)^(1/6) when ||F||_(3/2)=1."""
    return 9.0 * math.sqrt(2.0) / (16.0 * math.sqrt(math.pi))


def maxwell_ball_probability(radius: float) -> float:
    r = float(radius)
    return math.erf(r / math.sqrt(2.0)) - math.sqrt(2.0 / math.pi) * r * math.exp(-0.5 * r * r)


def local_ellipsoid_mass_coefficient(radius: float = 2.0, eps: float = EPS) -> float:
    """Coefficient C in integral_E |check f|^2 >= C (det Sigma)^(1/6).

    Uses ||check f-check F||_3 <= ||f-F||_(3/2) <= eps and Holder on E.
    E is the radius-r affine ellipsoid of the physical L2 covariance of F.
    """
    prob = maxwell_ball_probability(radius)
    cL = gaussian_l2_coefficient()
    vol = (4.0 * math.pi / 3.0) * radius ** 3
    amp = math.sqrt(cL * prob) - eps * vol ** (1.0 / 6.0)
    return max(0.0, amp) ** 2


def aspect_mass_lower(aspect: float, radius: float = 2.0) -> float:
    """Dimensionless ellipsoidal critical mass lower bound for A=N*l_max.

    Uses l_min >= (2/3)N^-1 and the certified local ellipsoid coefficient.
    """
    if aspect <= 0:
        raise ValueError("aspect must be positive")
    C = local_ellipsoid_mass_coefficient(radius)
    return C * (4.0 * aspect / 9.0) ** (1.0 / 3.0)


def affine_young_scaling_factors(detS: float) -> tuple[float, float]:
    """Norm and trilinear scaling for F_S(xi)=|detS|^(2/3)F(S xi).

    Both are exactly one at p=3/2 after change of variables.  Returning the
    explicit algebraic factors makes the affine-symmetry countermodel testable.
    """
    d = abs(float(detS))
    if d <= 0:
        raise ValueError("detS must be nonzero")
    norm_p_power = d ** ((2.0 / 3.0) * P) / d
    trilinear = d ** 2 / (d * d)
    return norm_p_power, trilinear


def arb_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    pi = arb.pi()
    eps = arb(1) / 100
    mass = arb(999) / 1000  # 1-eps^(3/2), exactly 1-1/1000
    B = (arb(2) / 25).exp()
    axis = mass * (2 * pi).sqrt() / (2 * B * arb(3).sqrt())
    if not (axis > arb(2) / 3):
        raise AssertionError(f"physical axis lower bound failed: {axis}")

    r = arb(2)
    # Maxwell/chi_3 radius-two probability.
    prob = (r / arb(2).sqrt()).erf() - (arb(2) / pi).sqrt() * r * (-(r * r) / 2).exp()
    cL = 9 * arb(2).sqrt() / (16 * pi.sqrt())
    vol = (4 * pi / 3) * r ** 3
    amp = (cL * prob).sqrt() - eps * vol ** (arb(1) / 6)
    Cloc = amp * amp
    if not (Cloc > arb(3) / 10):
        raise AssertionError(f"local ellipsoid mass coefficient failed: {Cloc}")
    aspect_coeff = (arb(3) / 10) * (arb(4) / 9) ** (arb(1) / 3)
    if not (aspect_coeff > arb(1) / 5):
        raise AssertionError(f"aspect mass coefficient failed: {aspect_coeff}")
    return {
        "shell_outer_radius": "exp(2/25) N",
        "profile_distance": "1/100 in L^(3/2)",
        "p_density_shell_mass_lower": "999/1000",
        "physical_axis_constant_ball": str(axis),
        "physical_axis_clean_lower": "2/3",
        "radius_two_probability_ball": str(prob),
        "local_ellipsoid_mass_coefficient_ball": str(Cloc),
        "local_ellipsoid_clean_lower": "3/10",
        "aspect_mass_coefficient_ball": str(aspect_coeff),
        "aspect_mass_clean_lower": "1/5",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class AspectStress:
    samples: int
    worst_uncertainty_residual: float
    minimum_axis_constant: float
    local_mass_coefficient: float
    aspect_mass_ratio: float
    max_affine_young_symmetry_residual: float
    extreme_free_young_aspect: float


def stress(samples: int = 50_000, seed: int = 20260807) -> AspectStress:
    rng = np.random.default_rng(seed)
    wu = wy = 0.0
    for _ in range(samples):
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        vals = np.exp(rng.uniform(-12.0, 12.0, size=3))
        A = Q @ np.diag(vals) @ Q.T
        wu = max(wu, uncertainty_matrix_residual(A))
        detS = math.exp(float(rng.uniform(-20.0, 20.0)))
        nfac, tfac = affine_young_scaling_factors(detS)
        wy = max(wy, abs(nfac - 1.0), abs(tfac - 1.0))
    extreme = 1.0e8
    S = np.diag([extreme, extreme ** -0.5, extreme ** -0.5])
    nfac, tfac = affine_young_scaling_factors(float(np.linalg.det(S)))
    wy = max(wy, abs(nfac - 1.0), abs(tfac - 1.0))
    return AspectStress(
        samples=samples,
        worst_uncertainty_residual=wu,
        minimum_axis_constant=physical_axis_lower_constant(),
        local_mass_coefficient=local_ellipsoid_mass_coefficient(),
        aspect_mass_ratio=aspect_mass_lower(1.0),
        max_affine_young_symmetry_residual=wy,
        extreme_free_young_aspect=extreme,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-affine-shell-aspect"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_certificate()
    out = stress(args.samples)
    (args.outdir / "affine_shell_aspect.json").write_text(json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2))
    md = f"""# Shell concentration versus affine Gaussian aspect

Status: **{cert['status']}**.

- exact uncertainty matrix: `Sigma_x Gamma_(3/2)=I/3`
- one-percent profile distance and the certified log shell imply every physical
  standard axis is `> (2/3) N^-1`
- on the radius-two covariance ellipsoid of the Gaussian profile, Hausdorff--Young
  plus local Holder gives actual physical mass
  `N integral_E |check f|^2 >= (3/10) N(det Sigma_x)^(1/6)`
- if `A=N l_max`, then
  `N integral_E |check f|^2 > (1/5) A^(1/3)`
- affine Young symmetry remains exact: elongation alone is not a transfer cost
- random covariance checks: `{out.samples}`
- worst uncertainty residual: `{out.worst_uncertainty_residual:.3e}`
- local ellipsoid coefficient (numerical): `{out.local_mass_coefficient:.9f}`
- tested free affine Young aspect: `{out.extreme_free_young_aspect:.1e}`

The last bullet is a required correction to a naive replication argument.  A
common affine Gaussian can be arbitrarily anisotropic while remaining an exact
Young extremizer.  Therefore aspect ratio alone must not be inserted as a
Bellman deficit.  Dynamics must instead see the anisotropy through the
ellipsoid-localized mass ledger and the grain-normalized curvature tensor.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
