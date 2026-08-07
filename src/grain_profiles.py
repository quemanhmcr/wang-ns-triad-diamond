from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

from .gaussian_packet import packet_center, scalar_frequency_ratio, scalar_position_ratio
from .helical import edge_metrics

P = 1.5
ALPHA = 2.0 / 3.0  # 1/P


def probability_vector(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or np.any(x < 0) or not np.isfinite(x).all():
        raise ValueError("expected a finite nonnegative vector")
    s = float(np.sum(x))
    if s <= 0:
        raise ValueError("vector must have positive sum")
    return x / s


def component_score(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Sharp component Bellman functional for p=3/2 mass probabilities."""
    x, y, z = map(probability_vector, (x, y, z))
    if not (len(x) == len(y) == len(z)):
        raise ValueError("vectors must have equal length")
    return float(np.sum((x * y * z) ** ALPHA))


def holder_component_bound(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Hölder upper bound S <= (sum x^2 sum y^2 sum z^2)^(1/3) <= 1."""
    x, y, z = map(probability_vector, (x, y, z))
    return float((np.sum(x*x) * np.sum(y*y) * np.sum(z*z)) ** (1.0 / 3.0))


@dataclass(frozen=True)
class DominantCertificate:
    score: float
    index: int
    x_mass: float
    y_mass: float
    z_mass: float
    x_lower: float
    yz_lower_each: float


def dominant_component_certificate(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> DominantCertificate:
    """A rigorous certificate extracted from the proof of the component lemma.

    If S is the component score and C* maximizes x, then
      x_C* >= S^3,
      y_C*, z_C* >= max(0, S^3 - sqrt(6(1-S))).

    The shared-component bound follows from the exact AM--GM stability
      (x^2+y^2+z^2)/3-(xyz)^(2/3)
        >= ((x-y)^2+(y-z)^2+(z-x)^2)/6.
    """
    x, y, z = map(probability_vector, (x, y, z))
    s = component_score(x, y, z)
    i = int(np.argmax(x))
    x_lower = s ** 3
    yz_lower = max(0.0, x_lower - math.sqrt(max(0.0, 6.0 * (1.0 - s))))
    return DominantCertificate(
        score=s,
        index=i,
        x_mass=float(x[i]),
        y_mass=float(y[i]),
        z_mass=float(z[i]),
        x_lower=float(x_lower),
        yz_lower_each=float(yz_lower),
    )


def softmax_pair(a: float) -> np.ndarray:
    if a >= 0:
        e = math.exp(-a)
        return np.array([1.0 / (1.0 + e), e / (1.0 + e)])
    e = math.exp(a)
    return np.array([e / (1.0 + e), 1.0 / (1.0 + e)])


def kernel_transfer(kernel: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> float:
    """Coarse grain transfer model with p-mass vectors x,y,z."""
    kernel = np.asarray(kernel, float)
    x, y, z = map(probability_vector, (x, y, z))
    if kernel.shape != (len(x), len(y), len(z)):
        raise ValueError("kernel shape mismatch")
    weights = (x[:, None, None] * y[None, :, None] * z[None, None, :]) ** ALPHA
    return float(np.sum(kernel * weights))


def generic_rotation(theta: float) -> np.ndarray:
    """Rodrigues rotation about a fixed generic axis that moves every center."""
    axis = np.array([1.0, 2.0, 3.0], dtype=float)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    eye = np.eye(3)
    return eye + math.sin(theta) * K + (1.0 - math.cos(theta)) * (K @ K)


def two_branch_kernel(sigma: float, angle: float, spatial_scaled: float) -> np.ndarray:
    """Exact scalar Gaussian cross-overlap times center helical efficiency.

    Branch 1 is a rigid rotation of branch 0, preserving its internal triad geometry.
    spatial_scaled means sigma * physical branch separation.
    """
    p0, q0, z0, signs, j_star = packet_center()
    rot = generic_rotation(angle)
    ps = [p0, rot @ p0]
    qs = [q0, rot @ q0]
    zs = [z0, rot @ z0]
    sep = spatial_scaled / sigma
    positions = [np.zeros(3), np.array([sep, 0.0, 0.0])]
    out = np.zeros((2, 2, 2), dtype=float)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                mismatch = ps[i] + qs[j] - zs[k]
                freq = scalar_frequency_ratio(mismatch, sigma)
                pos = scalar_position_ratio(np.array([positions[i], positions[j], positions[k]]), sigma)
                eff = edge_metrics(ps[i], qs[j], ps[i] + qs[j], *signs).efficiency / j_star
                # The h-packet is centered at zs[k], so mismatch is handled by the exact Gaussian factor.
                out[i, j, k] = freq * pos * max(0.0, min(1.0, eff))
    return out


def optimize_two_branch_kernel(kernel: np.ndarray, seed: int = 0) -> dict[str, object]:
    def objective(v: np.ndarray) -> float:
        x, y, z = (softmax_pair(float(v[i])) for i in range(3))
        return -kernel_transfer(kernel, x, y, z)

    res = differential_evolution(
        objective,
        [(-12.0, 12.0)] * 3,
        seed=seed,
        popsize=14,
        maxiter=220,
        tol=1e-11,
        polish=True,
    )
    x, y, z = (softmax_pair(float(res.x[i])) for i in range(3))
    return {
        "score": float(-res.fun),
        "x": x.tolist(),
        "y": y.tolist(),
        "z": z.tolist(),
        "dominant_common_mass": float(min(np.max(x), np.max(y), np.max(z))),
        "success": bool(res.success),
    }


def random_component_experiment(n: int, samples: int, seed: int) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    best_nonconcentrated = {"score": -1.0}
    max_violation = 0.0
    for _ in range(samples):
        concentration = 10.0 ** rng.uniform(-1.2, 1.2)
        x = rng.dirichlet(np.full(n, concentration))
        y = rng.dirichlet(np.full(n, concentration))
        z = rng.dirichlet(np.full(n, concentration))
        s = component_score(x, y, z)
        b = holder_component_bound(x, y, z)
        max_violation = max(max_violation, s - b, b - 1.0)
        common = max(min(x[i], y[i], z[i]) for i in range(n))
        if common <= 0.8 and s > best_nonconcentrated["score"]:
            best_nonconcentrated = {
                "score": float(s),
                "common_mass": float(common),
                "x": x.tolist(),
                "y": y.tolist(),
                "z": z.tolist(),
            }
    return {
        "components": n,
        "samples": samples,
        "max_numerical_inequality_violation": float(max_violation),
        "best_with_no_common_mass_above_0.8": best_nonconcentrated,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("results-grain-profiles"))
    parser.add_argument("--samples", type=int, default=200000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    random_results = [random_component_experiment(n, args.samples, 100+n) for n in (2, 3, 4, 8)]

    equal_two = component_score(np.array([0.5, 0.5]), np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    dominant_examples = []
    for a in [0.8, 0.9, 0.95, 0.99]:
        x = np.array([a, 1-a])
        cert = dominant_component_certificate(x, x, x)
        dominant_examples.append({
            "mass": a,
            "score": cert.score,
            "certificate": cert.__dict__,
        })

    branch_rows = []
    for sigma in [0.01, 0.02, 0.04]:
        for angle_over_sigma in [2.0, 4.0, 8.0, 12.0]:
            angle = angle_over_sigma * sigma
            for spatial_scaled in [0.0, 1.0, 2.0]:
                kernel = two_branch_kernel(sigma, angle, spatial_scaled)
                matched = [float(kernel[0,0,0]), float(kernel[1,1,1])]
                cross = [float(kernel[i,j,k]) for i in range(2) for j in range(2) for k in range(2) if not (i==j==k)]
                opt = optimize_two_branch_kernel(kernel, seed=int(1e5*sigma + 10*angle_over_sigma + spatial_scaled))
                branch_rows.append({
                    "sigma": sigma,
                    "angle_over_sigma": angle_over_sigma,
                    "sigma_times_spatial_separation": spatial_scaled,
                    "matched_min": min(matched),
                    "cross_max": max(cross),
                    "cross_sum": sum(cross),
                    "optimized": opt,
                })

    result = {
        "theorem": {
            "functional": "S=sum_C (X_C Y_C Z_C)^(2/3)",
            "sharp_upper_bound": 1.0,
            "two_equal_disconnected_components": equal_two,
            "dominant_certificate_formula": {
                "X_Cstar": "at least S^3",
                "Y_Cstar_and_Z_Cstar": "at least max(0,S^3-sqrt(6(1-S)))",
            },
        },
        "dominant_examples": dominant_examples,
        "random_checks": random_results,
        "two_branch_gaussian_model": branch_rows,
    }
    (args.outdir / "grain_profiles.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# Transfer-preserving grain profile experiment",
        "",
        "## Exact component Bellman lemma",
        "",
        f"Two equal disconnected components give normalized transfer `{equal_two:.12f}` (one half of the single-component optimum).",
        "",
        "| dominant p-mass per side | component score | certified X mass | certified Y,Z mass |",
        "|---:|---:|---:|---:|",
    ]
    for row in dominant_examples:
        c = row["certificate"]
        lines.append(f"| {row['mass']:.3f} | {row['score']:.9f} | {c['x_lower']:.9f} | {c['yz_lower_each']:.9f} |")
    lines += ["", "## Gaussian two-branch interaction", "", "| sigma | angle/sigma | sigma*space-sep | max cross edge | optimized score | common dominant mass |", "|---:|---:|---:|---:|---:|---:|"]
    for row in branch_rows:
        if row["angle_over_sigma"] in (4.0, 8.0, 12.0) and row["sigma_times_spatial_separation"] in (0.0, 2.0):
            lines.append(
                f"| {row['sigma']:.3f} | {row['angle_over_sigma']:.1f} | {row['sigma_times_spatial_separation']:.1f} | "
                f"{row['cross_max']:.3e} | {row['optimized']['score']:.9f} | {row['optimized']['dominant_common_mass']:.9f} |"
            )
    lines += [
        "",
        "The component lemma is exact. Gaussian cross-edge values and optimizations belong to the stated coarse-grain model; they are not a Navier--Stokes theorem.",
    ]
    (args.outdir / "summary.md").write_text("\n".join(lines)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()
