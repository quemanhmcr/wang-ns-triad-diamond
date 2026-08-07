from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .triad_extremizer import symmetric_rstar


@dataclass(frozen=True)
class Triad:
    left: str
    right: str
    child: str


def vertices_of(triads: Sequence[Triad]) -> list[str]:
    return sorted({v for t in triads for v in (t.left, t.right, t.child)})


def arcs_of(triads: Sequence[Triad]) -> list[tuple[str, str, int]]:
    arcs: list[tuple[str, str, int]] = []
    for i, t in enumerate(triads):
        arcs.append((t.left, t.child, i))
        arcs.append((t.right, t.child, i))
    return arcs


def incidence_matrix(triads: Sequence[Triad]) -> tuple[np.ndarray, list[str], list[tuple[str, str, int]]]:
    verts = vertices_of(triads)
    index = {v: i for i, v in enumerate(verts)}
    arcs = arcs_of(triads)
    B = np.zeros((len(arcs), len(verts)), dtype=float)
    for r, (tail, head, _) in enumerate(arcs):
        B[r, index[tail]] = -1.0
        B[r, index[head]] = 1.0
    return B, verts, arcs


def cycle_rank(triads: Sequence[Triad]) -> int:
    B, _, _ = incidence_matrix(triads)
    return int(B.shape[0] - np.linalg.matrix_rank(B, tol=1e-10))


def weighted_hodge_energy(
    triads: Sequence[Triad],
    gamma: float,
    weights: Sequence[float] | None = None,
) -> dict[str, object]:
    """Minimize sum_a w_a (ell_head-ell_tail-gamma)^2 by weighted least squares."""
    B, verts, arcs = incidence_matrix(triads)
    q = B.shape[0]
    w = np.ones(q) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != (q,) or np.min(w) <= 0:
        raise ValueError("positive weight per parent-child arc required")
    target = np.full(q, float(gamma))
    A = np.sqrt(w)[:, None] * B
    b = np.sqrt(w) * target
    ell, *_ = np.linalg.lstsq(A, b, rcond=None)
    residual = B @ ell - target
    energy = float(np.dot(w * residual, residual))

    # The optimal dual cycle vector is z=W r. Since B^T W r=0, it lies in cycle space.
    z = w * residual
    dual_stationarity = float(np.linalg.norm(B.T @ z))
    denom = float(np.dot(z / w, z))
    dual = 0.0 if denom <= 1e-30 else float(np.dot(z, target) ** 2 / denom)
    return {
        "vertices": verts,
        "arcs": arcs,
        "levels": {v: float(ell[i] / gamma) if gamma != 0 else float(ell[i]) for i, v in enumerate(verts)},
        "residual": residual.tolist(),
        "energy": energy,
        "dual_energy": dual,
        "dual_stationarity": dual_stationarity,
        "cycle_rank": int(q - np.linalg.matrix_rank(B, tol=1e-10)),
        "flat": bool(energy <= 1e-18),
    }


def triad_residual_identity(log_left: float, log_right: float, log_child: float, gamma: float) -> dict[str, float]:
    u = log_left - log_right
    v = log_child - 0.5 * (log_left + log_right) - gamma
    r_left = log_child - log_left - gamma
    r_right = log_child - log_right - gamma
    return {
        "u": u,
        "v": v,
        "r_left": r_left,
        "r_right": r_right,
        "lhs": r_left * r_left + r_right * r_right,
        "rhs": 0.5 * u * u + 2.0 * v * v,
    }


def nonflat_reuse_motif() -> list[Triad]:
    return [
        Triad("a", "b", "m"),
        Triad("m", "c", "d"),
        Triad("b", "c", "n"),
    ]


def flat_butterfly_motif() -> list[Triad]:
    return [
        Triad("a", "b", "m"),
        Triad("a", "c", "n"),
        Triad("m", "n", "d"),
    ]


def optimal_geometry(r_star: float) -> dict[str, float]:
    c = 1.0 / (2.0 * r_star * r_star) - 1.0
    theta = math.acos(c)
    R = 1.0 / r_star
    return {"c": c, "theta": theta, "R": R}


def butterfly_vectors(cos_theta: float) -> dict[str, np.ndarray]:
    if not (-1.0 < cos_theta < 1.0):
        raise ValueError("cos_theta must be in (-1,1)")
    s = math.sqrt(1.0 - cos_theta * cos_theta)
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([cos_theta, s, 0.0])
    c = np.array([cos_theta, -s, 0.0])
    R = math.sqrt(2.0 + 2.0 * cos_theta)
    m = (a + b) / R
    n = (a + c) / R
    d = (m + n) / R
    return {"a": a, "b": b, "c": c, "m": m, "n": n, "d": d}


def butterfly_rigidity_certificate(cos_theta: float) -> dict[str, float]:
    V = butterfly_vectors(cos_theta)
    a, b, c, m, n, d = (V[k] for k in ("a", "b", "c", "m", "n", "d"))
    R = math.sqrt(2.0 + 2.0 * cos_theta)
    theta = math.acos(cos_theta)
    return {
        "a_dot_b": float(a @ b),
        "a_dot_c": float(a @ c),
        "b_dot_c": float(b @ c),
        "required_b_dot_c": float(2.0 * cos_theta * cos_theta - 1.0),
        "m_dot_n": float(m @ n),
        "d_dot_a": float(d @ a),
        "d_error": float(np.linalg.norm(d - a)),
        "m_angle_from_d": float(math.acos(np.clip(m @ d, -1.0, 1.0))),
        "required_companion_angle": theta,
        "internal_angle_gap": float(theta - math.acos(np.clip(m @ d, -1.0, 1.0))),
        "R": R,
    }


def perturb_butterfly(cos_theta: float, eps: float, samples: int = 20000, seed: int = 0) -> dict[str, float | int | None]:
    """Numerical near-rigidity probe; not used as a theorem."""
    rng = np.random.default_rng(seed)
    V = butterfly_vectors(cos_theta)
    a, b0, c0 = V["a"], V["b"], V["c"]
    accepted = 0
    max_d_error = 0.0
    min_next_gap = math.inf
    theta_req = math.acos(cos_theta)
    for _ in range(samples):
        b = b0 + rng.normal(scale=eps, size=3)
        c = c0 + rng.normal(scale=eps, size=3)
        b /= np.linalg.norm(b)
        c /= np.linalg.norm(c)
        if abs(float(a @ b) - cos_theta) > eps or abs(float(a @ c) - cos_theta) > eps:
            continue
        m = (a + b) / np.linalg.norm(a + b)
        n = (a + c) / np.linalg.norm(a + c)
        if abs(float(m @ n) - cos_theta) > eps:
            continue
        d = (m + n) / np.linalg.norm(m + n)
        accepted += 1
        max_d_error = max(max_d_error, float(np.linalg.norm(d - a)))
        gaps = [abs(math.acos(np.clip(float(d @ x), -1.0, 1.0)) - theta_req) for x in (m, n)]
        min_next_gap = min(min_next_gap, min(gaps))
    return {
        "eps": eps,
        "samples": samples,
        "accepted": accepted,
        "max_d_error": max_d_error,
        "min_internal_companion_gap": None if accepted == 0 else min_next_gap,
        "exact_gap": 0.5 * theta_req,
    }



def planar_midpoint_children(angles: Sequence[float], theta: float, tolerance: float = 1e-12) -> list[float]:
    """All lifted-angle midpoints of pairs separated by theta.

    Angles live on the real line, not modulo 2*pi. This is the exact planar
    flat-extremal model inside a chosen angular chart.
    """
    a = sorted(float(x) for x in angles)
    out: list[float] = []
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if abs((a[j] - a[i]) - theta) <= tolerance:
                out.append(0.5 * (a[i] + a[j]))
    return sorted(out)


def planar_erosion_step(
    internal_angles: Sequence[float],
    theta: float,
    fresh_angles: Sequence[float] = (),
    pair_tolerance: float = 0.0,
    midpoint_error: float = 0.0,
) -> dict[str, float]:
    """Exact/near-exact angular erosion ledger.

    If parent pairs have separation at least theta-pair_tolerance and children
    are within midpoint_error of their angular midpoint, then
      diameter(next) <= diameter(internal union fresh)
                        -(theta-pair_tolerance)+2 midpoint_error.
    The returned bound is theorem-level; `children` are generated exactly only
    for the zero-tolerance diagnostic model.
    """
    internal = sorted(float(x) for x in internal_angles)
    if not internal:
        raise ValueError("nonempty internal angle set required")
    augmented = sorted(internal + [float(x) for x in fresh_angles])
    d_internal = internal[-1] - internal[0]
    d_augmented = augmented[-1] - augmented[0]
    fresh_expansion = d_augmented - d_internal
    bound = max(0.0, d_augmented - (theta - pair_tolerance) + 2.0 * midpoint_error)
    children = planar_midpoint_children(augmented, theta) if pair_tolerance == 0.0 and midpoint_error == 0.0 else []
    actual = 0.0 if len(children) <= 1 else children[-1] - children[0]
    return {
        "internal_diameter": d_internal,
        "augmented_diameter": d_augmented,
        "fresh_expansion": fresh_expansion,
        "next_diameter_bound": bound,
        "actual_exact_next_diameter": actual,
        "child_count": float(len(children)),
    }


def planar_fresh_span_lower_bound(
    initial_diameter: float,
    final_diameter: float,
    depth: int,
    theta: float,
    pair_tolerances: Sequence[float] | None = None,
    midpoint_errors: Sequence[float] | None = None,
) -> float:
    """Minimum total boundary expansion needed by a depth-L near-flat cascade."""
    if depth < 0:
        raise ValueError("depth must be nonnegative")
    pt = [0.0] * depth if pair_tolerances is None else list(pair_tolerances)
    me = [0.0] * depth if midpoint_errors is None else list(midpoint_errors)
    if len(pt) != depth or len(me) != depth:
        raise ValueError("one tolerance and midpoint error per level required")
    erosion = sum(theta - pt[j] - 2.0 * me[j] for j in range(depth))
    return max(0.0, final_diameter - initial_diameter + erosion)

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("results-cycle-hodge"))
    parser.add_argument("--samples", type=int, default=30000)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    r_star = symmetric_rstar()
    geom = optimal_geometry(r_star)
    gamma = math.log(1.0 / r_star)
    nonflat = weighted_hodge_energy(nonflat_reuse_motif(), gamma)
    butterfly = weighted_hodge_energy(flat_butterfly_motif(), gamma)
    cert = butterfly_rigidity_certificate(geom["c"])
    probes = [perturb_butterfly(geom["c"], eps, samples=args.samples, seed=100 + i) for i, eps in enumerate((1e-3, 3e-3, 1e-2, 3e-2))]

    result = {
        "r_star": r_star,
        "gamma": gamma,
        "geometry": geom,
        "nonflat_reuse": nonflat,
        "flat_butterfly": butterfly,
        "butterfly_certificate": cert,
        "near_butterfly_probes": probes,
    }
    (args.outdir / "cycle_hodge.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# Cycle Hodge and flat-cell rigidity",
        "",
        f"Nonflat reuse motif: cycle rank `{nonflat['cycle_rank']}`, Hodge energy `{nonflat['energy']:.12f}`.",
        f"Flat butterfly: cycle rank `{butterfly['cycle_rank']}`, Hodge energy `{butterfly['energy']:.3e}`.",
        "",
        "## Exact butterfly geometry",
        "",
        f"- cos(theta*): `{geom['c']:.12f}`",
        f"- theta*: `{geom['theta']:.12f}` rad",
        f"- scale ratio R*: `{geom['R']:.12f}`",
        f"- d returns to a with error: `{cert['d_error']:.3e}`",
        f"- internal next-companion angular gap: `{cert['internal_angle_gap']:.12f}` rad",
        "",
        "## Near-cell numerical probe",
        "",
        "| eps | accepted | max |d-a| | min internal companion gap |",
        "|---:|---:|---:|---:|",
    ]
    for row in probes:
        gap = row["min_internal_companion_gap"]
        lines.append(f"| {row['eps']:.4g} | {row['accepted']} | {row['max_d_error']:.6g} | {('n/a' if gap is None else f'{gap:.6g}')} |")
    lines += [
        "",
        "The Hodge identities and exact butterfly certificate are theorem-level finite-dimensional statements.",
        "The perturbation table is only a numerical stability probe.",
    ]
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
