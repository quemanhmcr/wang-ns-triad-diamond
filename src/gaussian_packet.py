from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import ndtri
from scipy.stats import qmc

from .helical import coupling_magnitude_closed
from .triad_extremizer import symmetric_jstar, symmetric_rstar


def packet_center() -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[int, int, int], float]:
    """Return the symmetric single-edge optimum normalized by |z|=1."""
    x = symmetric_rstar()
    h = math.sqrt(x * x - 0.25)
    p = np.array([0.5, h, 0.0])
    q = np.array([0.5, -h, 0.0])
    z = p + q
    signs = (1, -1, -1)
    return p, q, z, signs, symmetric_jstar(x)


def edge_efficiency_batch(p: np.ndarray, q: np.ndarray, signs: tuple[int, int, int]) -> np.ndarray:
    """Vectorized magnitude-only log-scale efficiency for z=p+q."""
    sx, sy, sz = signs
    z = p + q
    nx = np.linalg.norm(p, axis=1)
    ny = np.linalg.norm(q, axis=1)
    nz = np.linalg.norm(z, axis=1)
    cross = np.linalg.norm(np.cross(p, q), axis=1)
    area = 0.5 * cross
    safe = (nx > 1e-12) & (ny > 1e-12) & (nz > 1e-12)
    g = np.zeros_like(nx)
    g[safe] = area[safe] * np.abs(sx * nx[safe] + sy * ny[safe] + sz * nz[safe]) / (
        2.0 * math.sqrt(2.0) * nx[safe] * ny[safe] * nz[safe]
    )
    forward = np.zeros_like(nx)
    forward[safe] = nz[safe] / np.maximum(nx[safe], ny[safe])
    progress = np.maximum(0.0, np.log(np.maximum(forward, 1e-300)))
    raw = np.zeros_like(nx)
    raw[safe] = np.abs(sx * nx[safe] - sy * ny[safe]) * g[safe] / nz[safe]
    return progress * raw


def scalar_frequency_ratio(delta_k: np.ndarray, sigma: float) -> float:
    delta_k = np.asarray(delta_k, float)
    return float(math.exp(-float(delta_k @ delta_k) / (12.0 * sigma * sigma)))


def scalar_position_ratio(xs: np.ndarray, sigma: float) -> float:
    xs = np.asarray(xs, float)
    if xs.shape != (3, 3):
        raise ValueError("xs must have shape (3,3)")
    pair = sum(float(np.dot(xs[i] - xs[j], xs[i] - xs[j])) for i in range(3) for j in range(i + 1, 3))
    return float(math.exp(-(sigma * sigma / 3.0) * pair))


def scalar_width_ratio(sigmas: np.ndarray, dimension: int = 3) -> float:
    sigmas = np.asarray(sigmas, float)
    if sigmas.shape != (3,) or np.min(sigmas) <= 0:
        raise ValueError("sigmas must be three positive numbers")
    a = sigmas * sigmas / float(np.sum(sigmas * sigmas))
    return float((27.0 * np.prod(a)) ** (dimension / 4.0))


def sample_overlap_measure(
    k1: np.ndarray,
    k2: np.ndarray,
    k3: np.ndarray,
    sigma: float,
    power: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample the normalized Gaussian measure proportional to f1(p)f2(q)f3(p+q)."""
    sampler = qmc.Sobol(d=6, scramble=True, seed=seed)
    u = np.clip(sampler.random_base2(power), 1e-12, 1.0 - 1e-12)
    normals = ndtri(u).reshape(-1, 2, 3)

    mean_p = (2.0 * k1 - k2 + k3) / 3.0
    mean_q = (-k1 + 2.0 * k2 + k3) / 3.0
    covariance_2 = (2.0 * sigma * sigma / 3.0) * np.array([[2.0, -1.0], [-1.0, 2.0]])
    chol = np.linalg.cholesky(covariance_2)
    correlated = np.einsum("ab,nbc->nac", chol, normals)
    p = mean_p + correlated[:, 0, :]
    q = mean_q + correlated[:, 1, :]
    return p, q


def monte_carlo_cap_ratio(sigma: float, power: int = 17, seed: int = 0) -> dict[str, float]:
    k1, k2, k3, signs, j_star = packet_center()
    p, q = sample_overlap_measure(k1, k2, k3, sigma, power, seed)
    vals = edge_efficiency_batch(p, q, signs)
    ratio = vals / j_star
    return {
        "sigma": float(sigma),
        "samples": int(len(vals)),
        "mean_ratio": float(np.mean(ratio)),
        "std_ratio": float(np.std(ratio)),
        "q01": float(np.quantile(ratio, 0.01)),
        "q50": float(np.quantile(ratio, 0.50)),
        "q99": float(np.quantile(ratio, 0.99)),
        "center_efficiency": float(j_star),
    }



def sample_young_overlap_measure(
    k1: np.ndarray,
    k2: np.ndarray,
    k3: np.ndarray,
    sigma: float,
    power: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Overlap measure for symmetric L^(3/2) Gaussian Young extremizers."""
    sampler = qmc.Sobol(d=6, scramble=True, seed=seed)
    u = np.clip(sampler.random_base2(power), 1e-12, 1.0 - 1e-12)
    normals = ndtri(u).reshape(-1, 2, 3)
    mean_p = (2.0 * k1 - k2 + k3) / 3.0
    mean_q = (-k1 + 2.0 * k2 + k3) / 3.0
    covariance_2 = (sigma * sigma / 3.0) * np.array([[2.0, -1.0], [-1.0, 2.0]])
    chol = np.linalg.cholesky(covariance_2)
    correlated = np.einsum("ab,nbc->nac", chol, normals)
    return mean_p + correlated[:, 0, :], mean_q + correlated[:, 1, :]


def young_gaussian_ratio(sigma: float, power: int = 17, seed: int = 0) -> dict[str, float]:
    k1, k2, k3, signs, j_star = packet_center()
    p, q = sample_young_overlap_measure(k1, k2, k3, sigma, power, seed)
    vals = edge_efficiency_batch(p, q, signs) / j_star
    return {
        "sigma": float(sigma),
        "samples": int(len(vals)),
        "weighted_young_ratio": float(np.mean(vals)),
        "deficit_over_sigma": float((1.0 - np.mean(vals)) / sigma),
        "std_ratio": float(np.std(vals)),
    }


def sharp_young_constant(dimension: int = 3) -> float:
    """Sharp symmetric p=q=r=3/2 trilinear Young constant in R^d."""
    return float((math.sqrt(3.0) / 2.0) ** dimension)

def optimize_coherence_score() -> dict[str, object]:
    """Optimize the exact scalar coherence factors; maximum should be aligned/equal-width."""
    # Variables: normalized frequency mismatch, two independent spatial offsets,
    # and log-width perturbations (third log width fixed by zero-sum gauge).
    def objective(v: np.ndarray) -> float:
        dk = np.array([v[0], v[1], v[2]])
        xs = np.array([[0.0, 0.0, 0.0], [v[3], v[4], 0.0], [v[5], v[6], 0.0]])
        logs = np.array([v[7], v[8], -v[7] - v[8]])
        widths = np.exp(logs)
        score = scalar_frequency_ratio(dk, 1.0) * scalar_position_ratio(xs, 1.0) * scalar_width_ratio(widths)
        return -score

    bounds = [(-2.0, 2.0)] * 3 + [(-2.0, 2.0)] * 4 + [(-1.0, 1.0)] * 2
    res = differential_evolution(objective, bounds, seed=13, popsize=12, maxiter=180, tol=1e-10, polish=True)
    return {"score": float(-res.fun), "params": [float(x) for x in res.x], "success": bool(res.success)}


def symmetric_curvature() -> dict[str, float]:
    _, _, _, _, j_star = packet_center()
    x_star = packet_center()[0]
    x = float(np.linalg.norm(x_star))

    def j(t: float) -> float:
        return math.sqrt(max(0.0, 4.0 * t * t - 1.0)) * math.log(1.0 / t) / (4.0 * math.sqrt(2.0) * t)

    h = 1e-5
    second = (j(x + h) - 2.0 * j(x) + j(x - h)) / (h * h)
    return {"x_star": x, "j_star": j_star, "second_derivative": float(second), "normalized_curvature": float(-second / j_star)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--power", type=int, default=17)
    parser.add_argument("--outdir", type=Path, default=Path("results-packet"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    sigmas = [0.005, 0.01, 0.02, 0.04, 0.08, 0.12]
    cap = [monte_carlo_cap_ratio(s, power=args.power, seed=100 + i) for i, s in enumerate(sigmas)]
    young = [young_gaussian_ratio(s, power=args.power, seed=300 + i) for i, s in enumerate(sigmas)]

    mismatch = []
    for a in [0.0, 0.25, 0.5, 1.0, 2.0]:
        mismatch.append({"delta_over_sigma": a, "ratio": scalar_frequency_ratio(np.array([a, 0.0, 0.0]), 1.0)})

    separation = []
    for a in [0.0, 0.25, 0.5, 1.0, 2.0]:
        xs = np.array([[0.0, 0.0, 0.0], [a, 0.0, 0.0], [0.0, 0.0, 0.0]])
        separation.append({"sigma_times_separation": a, "ratio": scalar_position_ratio(xs, 1.0)})

    widths = []
    for triplet in ([1, 1, 1], [0.8, 1.0, 1.2], [0.6, 1.0, 1.4], [0.5, 1.0, 1.5]):
        widths.append({"sigmas": triplet, "ratio": scalar_width_ratio(np.array(triplet, float))})

    result = {
        "model": "L2-normalized isotropic Fourier Gaussian packets",
        "dimension": 3,
        "cap_coefficient_stability": cap,
        "weighted_young_gaussians": young,
        "sharp_young_constant": sharp_young_constant(),
        "frequency_mismatch_exact": mismatch,
        "spatial_separation_exact": separation,
        "width_balance_exact": widths,
        "coherence_optimization": optimize_coherence_score(),
        "single_edge_curvature": symmetric_curvature(),
    }
    (args.outdir / "packet_inverse.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# Gaussian packet inverse experiment",
        "",
        "The exact scalar envelope factor separates into frequency resonance, spatial overlap, and width balance.",
        "",
        f"Coherence optimizer score: `{result['coherence_optimization']['score']:.12f}` (ideal is 1).",
        "",
        "## Helical coefficient stability over narrow caps",
        "",
        "| sigma / child frequency | mean J/J* | deficit | deficit / sigma |",
        "|---:|---:|---:|---:|",
    ]
    for row in cap:
        deficit = 1.0 - row["mean_ratio"]
        lines.append(f"| {row['sigma']:.4f} | {row['mean_ratio']:.9f} | {deficit:.9f} | {deficit/row['sigma']:.6f} |")
    lines += ["", "## Scale-critical weighted Young Gaussian packets", "", f"Sharp scalar constant in R^3: `{result['sharp_young_constant']:.12f}`", "", "| sigma / child frequency | weighted ratio | deficit / sigma |", "|---:|---:|---:|"]
    for row in young:
        lines.append(f"| {row['sigma']:.4f} | {row['weighted_young_ratio']:.9f} | {row['deficit_over_sigma']:.6f} |")
    curv = result["single_edge_curvature"]
    lines += [
        "",
        "## Symmetric single-edge stability",
        "",
        f"- x*: `{curv['x_star']:.12f}`",
        f"- J*: `{curv['j_star']:.12f}`",
        f"- J''(x*): `{curv['second_derivative']:.9f}`",
        f"- normalized local curvature -J''/J*: `{curv['normalized_curvature']:.9f}`",
        "",
        "Numerical cap averages are experimental. The scalar Gaussian overlap formulas are exact for the stated model.",
    ]
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
