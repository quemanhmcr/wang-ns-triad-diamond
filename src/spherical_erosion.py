from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import minimize

R_STAR = 1.6369181300426772
C_STAR = R_STAR / 2.0
THETA_STAR = 2.0 * math.acos(C_STAR)
KAPPA_STAR = -math.log(C_STAR)
TETRAHEDRAL_DIAMETER = math.acos(-1.0 / 3.0)


def normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    n = float(np.linalg.norm(x))
    if n <= 1e-14:
        raise ValueError("cannot normalize zero vector")
    return x / n


def spherical_angle(x: np.ndarray, y: np.ndarray) -> float:
    x = normalize(x)
    y = normalize(y)
    return float(math.acos(float(np.clip(x @ y, -1.0, 1.0))))


def spherical_midpoint(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    p = normalize(p)
    q = normalize(q)
    if float(p @ q) <= -1.0 + 1e-13:
        raise ValueError("antipodal points have no unique short-arc midpoint")
    return normalize(p + q)


def hemisphere_barrier(x: np.ndarray, pole: np.ndarray) -> float:
    x = normalize(x)
    pole = normalize(pole)
    support = float(pole @ x)
    if support <= 0.0:
        raise ValueError("point is outside the open hemisphere")
    return -math.log(support)


def midpoint_barrier_margin(p: np.ndarray, q: np.ndarray, pole: np.ndarray) -> float:
    """Return RHS-LHS in the exact spherical barrier inequality.

    For theta=d(p,q), kappa(theta)=-log cos(theta/2),
      Phi(mid(p,q)) <= (Phi(p)+Phi(q))/2 - kappa(theta).
    The return value is nonnegative up to floating-point error.
    """
    theta = spherical_angle(p, q)
    midpoint = spherical_midpoint(p, q)
    kappa = -math.log(math.cos(theta / 2.0))
    lhs = hemisphere_barrier(midpoint, pole)
    rhs = 0.5 * (hemisphere_barrier(p, pole) + hemisphere_barrier(q, pole)) - kappa
    return float(rhs - lhs)


def make_companion(parent: np.ndarray, tangent: np.ndarray, theta: float = THETA_STAR) -> np.ndarray:
    parent = normalize(parent)
    tangent = np.asarray(tangent, dtype=float)
    tangent = tangent - float(tangent @ parent) * parent
    tangent = normalize(tangent)
    return normalize(math.cos(theta) * parent + math.sin(theta) * tangent)


def lineage_barrier_certificate(
    x0: np.ndarray,
    companions: Iterable[np.ndarray],
    pole: np.ndarray,
    theta: float = THETA_STAR,
) -> dict[str, object]:
    """Build the exact midpoint lineage and verify the additive companion ledger."""
    pole = normalize(pole)
    states = [normalize(x0)]
    phis = [hemisphere_barrier(states[0], pole)]
    companion_phis: list[float] = []
    step_margins: list[float] = []
    for q in companions:
        q = normalize(q)
        actual_theta = spherical_angle(states[-1], q)
        if abs(actual_theta - theta) > 1e-8:
            raise ValueError("companion does not have the requested angular separation")
        companion_phis.append(hemisphere_barrier(q, pole))
        step_margins.append(midpoint_barrier_margin(states[-1], q, pole))
        nxt = spherical_midpoint(states[-1], q)
        states.append(nxt)
        phis.append(hemisphere_barrier(nxt, pole))

    L = len(companion_phis)
    kappa = -math.log(math.cos(theta / 2.0))
    exact_rhs = 2.0 * kappa * L - phis[0] + sum(phis[1:-1]) + (2.0 * phis[-1] if L else 0.0)
    simple_rhs = 2.0 * kappa * L - phis[0]
    total = float(sum(companion_phis))
    return {
        "states": [x.tolist() for x in states],
        "state_barriers": phis,
        "companion_barriers": companion_phis,
        "step_margins": step_margins,
        "total_companion_barrier": total,
        "exact_rhs": float(exact_rhs),
        "simple_rhs": float(simple_rhs),
        "exact_margin": float(total - exact_rhs),
        "simple_margin": float(total - simple_rhs),
    }


def next_cap_support(current_support: float, theta: float = THETA_STAR) -> float:
    """Exact lower support after one no-fresh midpoint generation."""
    if not 0.0 < current_support <= 1.0:
        raise ValueError("support must lie in (0,1]")
    return current_support / math.cos(theta / 2.0)


def no_fresh_lifespan_bound(cap_radius: float, theta: float = THETA_STAR) -> int:
    """Maximum number of nonempty exact generations allowed by cap support."""
    if not 0.0 <= cap_radius < math.pi / 2.0:
        raise ValueError("cap radius must be in [0,pi/2)")
    support = math.cos(cap_radius)
    kappa = -math.log(math.cos(theta / 2.0))
    return int(math.floor(math.log(1.0 / support) / kappa + 1e-12))


def convex_hull_distance(points: np.ndarray) -> dict[str, object]:
    """Numerically compute dist(0,conv(points)) and a simplex witness."""
    points = np.asarray(points, dtype=float)
    points = np.array([normalize(x) for x in points])
    n = len(points)
    x0 = np.full(n, 1.0 / n)

    def objective(lam: np.ndarray) -> float:
        v = lam @ points
        return 0.5 * float(v @ v)

    cons = {"type": "eq", "fun": lambda lam: float(np.sum(lam) - 1.0)}
    res = minimize(objective, x0, method="SLSQP", bounds=[(0.0, 1.0)] * n, constraints=cons,
                   options={"ftol": 1e-14, "maxiter": 3000})
    lam = np.clip(res.x, 0.0, 1.0)
    lam /= np.sum(lam)
    witness = lam @ points
    return {
        "distance": float(np.linalg.norm(witness)),
        "weights": lam.tolist(),
        "witness": witness.tolist(),
        "success": bool(res.success),
    }


def angular_diameter(points: np.ndarray) -> float:
    points = np.asarray(points, dtype=float)
    if len(points) <= 1:
        return 0.0
    return max(spherical_angle(points[i], points[j]) for i in range(len(points)) for j in range(i + 1, len(points)))


def regular_tetrahedron() -> np.ndarray:
    pts = np.array([
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ])
    return np.array([normalize(x) for x in pts])


def equal_marginal_barycenter_test(pairs: list[tuple[np.ndarray, np.ndarray]], weights: np.ndarray) -> dict[str, object]:
    """Use the symmetric coupling (p,q),(q,p), so both parent marginals agree."""
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.sum(weights)
    parent_b = np.zeros(3)
    child_b = np.zeros(3)
    cos_halves = []
    for w, (p, q) in zip(weights, pairs):
        p = normalize(p)
        q = normalize(q)
        theta = spherical_angle(p, q)
        cos_half = math.cos(theta / 2.0)
        cos_halves.append(cos_half)
        parent_b += 0.5 * w * (p + q)
        child_b += w * spherical_midpoint(p, q)
    if max(cos_halves) - min(cos_halves) > 1e-10:
        raise ValueError("all pairs must have one common angle")
    c = float(np.mean(cos_halves))
    predicted = parent_b / c
    return {
        "parent_barycenter": parent_b.tolist(),
        "child_barycenter": child_b.tolist(),
        "predicted_child": predicted.tolist(),
        "error": float(np.linalg.norm(child_b - predicted)),
        "amplification": float(1.0 / c),
    }



def collision_entropy(weights: np.ndarray) -> float:
    weights = np.asarray(weights, dtype=float)
    if np.min(weights) < 0 or float(np.sum(weights)) <= 0:
        raise ValueError("weights must be nonnegative and nontrivial")
    weights = weights / np.sum(weights)
    return -math.log(float(np.sum(weights * weights)))


def barycenter_collision_certificate(points: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    """Exact atomic inequality H_2 >= log(2/(1+|b|))."""
    points = np.asarray(points, dtype=float)
    points = np.array([normalize(x) for x in points])
    weights = np.asarray(weights, dtype=float)
    weights = weights / np.sum(weights)
    barycenter = weights @ points
    bary_norm = float(np.linalg.norm(barycenter))
    q = float(np.sum(weights * weights))
    upper_q = 0.5 * (1.0 + bary_norm)
    return {
        "barycenter_norm": bary_norm,
        "collision_probability": q,
        "collision_upper_bound": upper_q,
        "entropy": -math.log(q),
        "entropy_lower_bound": math.log(2.0 / (1.0 + bary_norm)),
        "margin": upper_q - q,
    }


def balanced_chain_entropy_bound(depth: int, cos_half: float = C_STAR) -> dict[str, float]:
    """Entropy forced by a depth-L equal-marginal flat chain.

    If the chain survives from level j to depth L, |b_j| <= cos_half^(L-j).
    Hence H_2(mu_j) >= log(2/(1+cos_half^(L-j))).
    """
    if depth < 0 or not 0.0 < cos_half < 1.0:
        raise ValueError("invalid depth or cosine")
    terms = [math.log(2.0 / (1.0 + cos_half ** r)) for r in range(1, depth + 1)]
    total = float(sum(terms))
    deficit = float(depth * math.log(2.0) - total)
    return {
        "depth": depth,
        "total_entropy_lower_bound": total,
        "linear_main_term": float(depth * math.log(2.0)),
        "bounded_deficit": deficit,
        "per_level_terms": terms,
    }

def random_unit(rng: np.random.Generator) -> np.ndarray:
    return normalize(rng.normal(size=3))


def random_barrier_probe(samples: int, seed: int = 0) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    minimum_margin = float("inf")
    accepted = 0
    for _ in range(samples):
        p = random_unit(rng)
        tangent = random_unit(rng)
        q = make_companion(p, tangent)
        midpoint = spherical_midpoint(p, q)
        pole = normalize(midpoint + 0.35 * random_unit(rng))
        if min(float(pole @ p), float(pole @ q), float(pole @ midpoint)) <= 1e-9:
            continue
        accepted += 1
        minimum_margin = min(minimum_margin, midpoint_barrier_margin(p, q, pole))
    return {"samples": samples, "accepted": accepted, "minimum_margin": float(minimum_margin)}


def build_lineage(length: int, seed: int = 1) -> tuple[np.ndarray, list[np.ndarray], np.ndarray]:
    rng = np.random.default_rng(seed)
    pole = np.array([0.0, 0.0, 1.0])
    state = normalize(np.array([0.5, 0.0, 1.0]))
    x0 = state.copy()
    companions: list[np.ndarray] = []
    for _ in range(length):
        found = False
        # Bias the tangent so the companion stays in the northern hemisphere.
        for _attempt in range(10000):
            tangent = random_unit(rng) + 1.5 * pole
            q = make_companion(state, tangent)
            nxt = spherical_midpoint(state, q)
            if float(pole @ q) > 0.03 and float(pole @ nxt) > 0.03:
                companions.append(q)
                state = nxt
                found = True
                break
        if not found:
            break
    return x0, companions, pole


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=50000)
    parser.add_argument("--outdir", type=Path, default=Path("results-spherical-erosion"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    probe = random_barrier_probe(args.samples, seed=71)
    x0, companions, pole = build_lineage(12, seed=19)
    lineage = lineage_barrier_certificate(x0, companions, pole)

    cap_rows = []
    for degrees in [30, 45, 60, 75, 85]:
        radius = math.radians(degrees)
        cap_rows.append({
            "radius_degrees": degrees,
            "support": math.cos(radius),
            "lifespan_bound": no_fresh_lifespan_bound(radius),
        })

    tetra = regular_tetrahedron()
    tetra_hull = convex_hull_distance(tetra)
    tetra_diameter = angular_diameter(tetra)

    rng = np.random.default_rng(101)
    pairs = []
    for _ in range(40):
        p = random_unit(rng)
        q = make_companion(p, random_unit(rng))
        pairs.append((p, q))
    bary_weights = rng.random(len(pairs))
    bary = equal_marginal_barycenter_test(pairs, bary_weights)
    atomic_points = np.array([p for pair in pairs for p in pair])
    atomic_weights = np.repeat(bary_weights / np.sum(bary_weights) / 2.0, 2)
    collision = barycenter_collision_certificate(atomic_points, atomic_weights)
    entropy_rows = [balanced_chain_entropy_bound(L) for L in [4, 8, 12, 24, 48]]
    entropy_constant = sum(math.log(1.0 + C_STAR ** r) for r in range(1, 20000))

    result = {
        "constants": {
            "r_star": 1.0 / R_STAR,
            "R_star": R_STAR,
            "theta_star": THETA_STAR,
            "theta_star_degrees": math.degrees(THETA_STAR),
            "cos_half": C_STAR,
            "kappa_star": KAPPA_STAR,
            "tetrahedral_diameter": TETRAHEDRAL_DIAMETER,
            "tetrahedral_diameter_degrees": math.degrees(TETRAHEDRAL_DIAMETER),
        },
        "barrier_probe": probe,
        "lineage": lineage,
        "cap_lifespans": cap_rows,
        "tetrahedron": {
            "convex_hull": tetra_hull,
            "diameter": tetra_diameter,
            "diameter_error": tetra_diameter - TETRAHEDRAL_DIAMETER,
        },
        "barycenter_amplification": bary,
        "barycenter_collision": collision,
        "balanced_chain_entropy": entropy_rows,
        "balanced_entropy_deficit_constant": entropy_constant,
    }
    (args.outdir / "spherical_erosion.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# Spherical flat-network erosion",
        "",
        f"- theta*: `{THETA_STAR:.12f}` rad (`{math.degrees(THETA_STAR):.6f}` deg)",
        f"- kappa*=-log cos(theta*/2): `{KAPPA_STAR:.12f}`",
        f"- tetrahedral non-hemisphere threshold: `{math.degrees(TETRAHEDRAL_DIAMETER):.6f}` deg",
        "",
        "## Exact barrier probe",
        "",
        f"Accepted samples: `{probe['accepted']}` / `{probe['samples']}`",
        f"Minimum RHS-LHS margin: `{probe['minimum_margin']:.3e}`",
        "",
        "## Exact lineage ledger",
        "",
        f"Lineage length: `{len(companions)}`",
        f"Total companion barrier: `{lineage['total_companion_barrier']:.9f}`",
        f"Exact telescoping RHS: `{lineage['exact_rhs']:.9f}`",
        f"Exact margin: `{lineage['exact_margin']:.3e}`",
        "",
        "## No-fresh hemispherical lifespan",
        "",
        "| cap radius | support | maximum exact generations |",
        "|---:|---:|---:|",
    ]
    for row in cap_rows:
        lines.append(f"| {row['radius_degrees']} deg | {row['support']:.6f} | {row['lifespan_bound']} |")
    lines += [
        "",
        "## Balanced exception",
        "",
        f"Regular tetrahedron hull distance: `{tetra_hull['distance']:.3e}`",
        f"Regular tetrahedron diameter: `{math.degrees(tetra_diameter):.9f}` deg",
        f"Equal-marginal barycenter amplification error: `{bary['error']:.3e}`",
        f"Atomic collision inequality margin: `{collision['margin']:.3e}`",
        f"Asymptotic entropy deficit constant: `{entropy_constant:.9f}`",
        "",
        "## Balanced-chain entropy",
        "",
        "| depth | entropy lower bound | L log 2 - bound |",
        "|---:|---:|---:|",
    ]
    for row in entropy_rows:
        lines.append(f"| {row['depth']} | {row['total_entropy_lower_bound']:.9f} | {row['bounded_deficit']:.9f} |")
    lines += [
        "",
        "The midpoint barrier, lifespan bound, tetrahedral diameter certificate, and barycenter identity are exact finite-dimensional statements. The random probe only tests the implementation.",
    ]
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
