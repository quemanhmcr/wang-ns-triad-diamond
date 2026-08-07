from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

from .helical import edge_metrics

Array = np.ndarray


def _prob(v: Array, axis: int | None = None) -> Array:
    v = np.asarray(v, dtype=float)
    if np.any(v < 0):
        raise ValueError("probability weights must be nonnegative")
    s = np.sum(v, axis=axis, keepdims=True)
    if np.any(s <= 0):
        raise ValueError("probability weights must have positive sum")
    return v / s


def component_score(x: Array, y: Array, z: Array) -> float:
    """Bellman score sum_C (X_C Y_C Z_C)^(2/3)."""
    x, y, z = map(lambda a: _prob(np.asarray(a, float)), (x, y, z))
    if not (x.shape == y.shape == z.shape):
        raise ValueError("component arrays must have equal shape")
    return float(np.sum((x * y * z) ** (2.0 / 3.0)))


@dataclass(frozen=True)
class StepCertificate:
    parent_score: float
    exact_ratio: float
    entropy_bound: float
    reuse_average: float
    combined_bound: float
    entropy_cost: float
    reuse_cost: float
    total_cost: float
    cross_error: float
    error_correction: float
    observed_ratio_bound: float


def refinement_certificate(
    x: Array,
    y: Array,
    z: Array,
    ax: Array,
    ay: Array,
    az: Array,
    rho: Array | None = None,
    cross_error: float = 0.0,
) -> StepCertificate:
    """Certify one transfer-adapted refinement step.

    x,y,z are parent p=3/2 mass distributions. Rows of ax,ay,az are
    conditional child distributions inside each parent component. rho[v]
    is a local reuse-efficiency factor in [0,1].
    """
    x, y, z = map(lambda a: _prob(np.asarray(a, float)), (x, y, z))
    ax, ay, az = map(lambda a: _prob(np.asarray(a, float), axis=1), (ax, ay, az))
    if not (x.shape == y.shape == z.shape):
        raise ValueError("parent distributions must have equal shape")
    if not (ax.shape == ay.shape == az.shape and ax.shape[0] == x.size):
        raise ValueError("conditional arrays must have shape (parents, children)")
    if rho is None:
        rho = np.ones_like(x)
    rho = np.asarray(rho, float)
    if rho.shape != x.shape or np.any((rho < 0) | (rho > 1)):
        raise ValueError("rho must lie in [0,1] and match parent shape")
    if cross_error < 0:
        raise ValueError("cross_error must be nonnegative")

    parent_terms = (x * y * z) ** (2.0 / 3.0)
    s0 = float(np.sum(parent_terms))
    if s0 <= 0:
        raise ValueError("parent Bellman score is zero")
    lam = parent_terms / s0

    cx = np.sum(ax * ax, axis=1)
    cy = np.sum(ay * ay, axis=1)
    cz = np.sum(az * az, axis=1)
    local_bellman = np.sum((ax * ay * az) ** (2.0 / 3.0), axis=1)
    geometric_collision = (cx * cy * cz) ** (1.0 / 3.0)

    exact_ratio = float(np.sum(lam * rho * local_bellman))
    bx = float(np.sum(lam * cx))
    by = float(np.sum(lam * cy))
    bz = float(np.sum(lam * cz))
    entropy_bound = float((bx * by * bz) ** (1.0 / 3.0))

    denom = float(np.sum(lam * geometric_collision))
    if denom <= 0:
        reuse_average = 0.0
    else:
        reuse_average = float(np.sum(lam * rho * geometric_collision) / denom)
    combined_bound = reuse_average * entropy_bound

    # exact_ratio <= sum lam*rho*g <= reuse_average*(sum lam*g)
    # and sum lam*g <= (sum lam*cx sum lam*cy sum lam*cz)^(1/3).
    tol = 5e-12
    if exact_ratio > combined_bound + tol:
        raise AssertionError((exact_ratio, combined_bound))

    entropy_cost = -math.log(max(entropy_bound, 1e-300))
    reuse_cost = -math.log(max(reuse_average, 1e-300))
    total_cost = entropy_cost + reuse_cost

    observed_ratio_bound = min(1.0, combined_bound + cross_error)
    error_correction = math.log1p(cross_error / max(combined_bound, 1e-300))
    return StepCertificate(
        parent_score=s0,
        exact_ratio=exact_ratio,
        entropy_bound=entropy_bound,
        reuse_average=reuse_average,
        combined_bound=combined_bound,
        entropy_cost=entropy_cost,
        reuse_cost=reuse_cost,
        total_cost=total_cost,
        cross_error=float(cross_error),
        error_correction=float(error_correction),
        observed_ratio_bound=float(observed_ratio_bound),
    )


def refine_masses(parent: Array, conditional: Array) -> Array:
    parent = _prob(np.asarray(parent, float))
    conditional = _prob(np.asarray(conditional, float), axis=1)
    if conditional.shape[0] != parent.size:
        raise ValueError("row count must match parent count")
    return (parent[:, None] * conditional).reshape(-1)


def holonomy_convex_cost(gamma: float, a_linear: float, b_quadratic: float) -> float:
    """Exact minimum of A|s| + B t^2 subject to s+t=gamma.

    This is the abstract cost obtained from a linear cusp defect plus a
    quadratic scale-shift defect and the log-scale holonomy identity.
    """
    if gamma <= 0 or a_linear <= 0 or b_quadratic <= 0:
        raise ValueError("parameters must be positive")
    if 2.0 * b_quadratic * gamma <= a_linear:
        return b_quadratic * gamma * gamma
    return a_linear * gamma - a_linear * a_linear / (4.0 * b_quadratic)


def single_edge_optimum() -> tuple[float, float, float]:
    fun = lambda r: -math.sqrt(max(0.0, 4.0 * r * r - 1.0)) * math.log(1.0 / r) / (4.0 * math.sqrt(2.0) * r)
    res = minimize_scalar(fun, bounds=(0.500000001, 0.999999), method="bounded", options={"xatol": 1e-14})
    r = float(res.x)
    return r, -float(res.fun), math.log(1.0 / r)


def efficiency_xy(x: float, y: float) -> float:
    """Exact helical edge efficiency for child |z|=1 and signs (+,-,-)."""
    if x <= 0 or y <= 0 or x + y <= 1 or abs(x - y) >= 1:
        return 0.0
    cos_theta = (1.0 - x * x - y * y) / (2.0 * x * y)
    if abs(cos_theta) > 1.0:
        return 0.0
    sin_theta = math.sqrt(max(0.0, 1.0 - cos_theta * cos_theta))
    p = np.array([x, 0.0, 0.0])
    q = np.array([y * cos_theta, y * sin_theta, 0.0])
    z = p + q
    return edge_metrics(p, q, z, 1, -1, -1).efficiency


def search_local_stability(
    box: tuple[float, float] = (0.54, 0.72),
    seed: int = 17,
) -> dict[str, float]:
    """Numerically search conservative A,B in d >= A|u|+Bv^2.

    This is experimental, not an interval-arithmetic certificate.
    """
    rstar, jstar, gamma = single_edge_optimum()

    # Optimize A and B with a small safety objective: maximize A+0.08B
    # subject to sampled/differential-evolution minimum margin >= 0.
    candidates: list[tuple[float, float, float]] = []
    for a in np.linspace(0.04, 0.30, 14):
        for b in np.linspace(0.5, 8.0, 16):
            def margin(v: Array) -> float:
                x, y = float(v[0]), float(v[1])
                j = efficiency_xy(x, y)
                deficit = 1.0 - j / jstar
                u = math.log(x / y)
                w = -0.5 * (math.log(x) + math.log(y)) - gamma
                return deficit - a * abs(u) - b * w * w

            res = differential_evolution(
                margin,
                bounds=[box, box],
                seed=seed,
                maxiter=120,
                popsize=10,
                polish=True,
                tol=1e-9,
            )
            m = float(res.fun)
            if m >= -2e-6:
                candidates.append((a + 0.08 * b, a, b))
    if not candidates:
        return {"A": 0.0, "B": 0.0, "minimum_margin": -1.0, "rstar": rstar, "jstar": jstar, "gamma": gamma}
    _, a, b = max(candidates)

    def final_margin(v: Array) -> float:
        x, y = map(float, v)
        deficit = 1.0 - efficiency_xy(x, y) / jstar
        u = math.log(x / y)
        w = -0.5 * (math.log(x) + math.log(y)) - gamma
        return deficit - a * abs(u) - b * w * w

    final = differential_evolution(
        final_margin,
        bounds=[box, box],
        seed=seed + 1,
        maxiter=500,
        popsize=20,
        polish=True,
        tol=1e-11,
    )
    cost = holonomy_convex_cost(gamma, a, b)
    return {
        "A": float(a),
        "B": float(b),
        "minimum_margin": float(final.fun),
        "argmin_x": float(final.x[0]),
        "argmin_y": float(final.x[1]),
        "rstar": rstar,
        "jstar": jstar,
        "gamma": gamma,
        "holonomy_cost": float(cost),
        "reuse_factor_candidate": float(math.exp(-cost)),
    }


def random_step_checks(samples: int, seed: int = 0) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    worst = -1e9
    worst_log = -1e9
    for _ in range(samples):
        parents = int(rng.integers(1, 7))
        children = int(rng.integers(2, 7))
        x = rng.dirichlet(np.ones(parents))
        y = rng.dirichlet(np.ones(parents))
        z = rng.dirichlet(np.ones(parents))
        ax = rng.dirichlet(np.ones(children), size=parents)
        ay = rng.dirichlet(np.ones(children), size=parents)
        az = rng.dirichlet(np.ones(children), size=parents)
        rho = rng.uniform(0.72, 1.0, size=parents)
        cert = refinement_certificate(x, y, z, ax, ay, az, rho)
        worst = max(worst, cert.exact_ratio - cert.combined_bound)
        if cert.exact_ratio > 0:
            worst_log = max(worst_log, cert.total_cost + math.log(cert.exact_ratio))
    return {"samples": samples, "max_bound_violation": float(worst), "max_log_violation": float(worst_log)}


def equal_branch_cascade(depth: int, branches: int, reuse_factor: float = 1.0) -> dict[str, float]:
    x = y = z = np.array([1.0])
    product_exact = 1.0
    sum_entropy = 0.0
    sum_reuse = 0.0
    for _ in range(depth):
        m = x.size
        conditional = np.full((m, branches), 1.0 / branches)
        rho = np.full(m, reuse_factor)
        cert = refinement_certificate(x, y, z, conditional, conditional, conditional, rho)
        product_exact *= cert.exact_ratio
        sum_entropy += cert.entropy_cost
        sum_reuse += cert.reuse_cost
        x = refine_masses(x, conditional)
        y = refine_masses(y, conditional)
        z = refine_masses(z, conditional)
    return {
        "depth": depth,
        "branches": branches,
        "reuse_factor": reuse_factor,
        "product_exact": product_exact,
        "predicted_product": math.exp(-(sum_entropy + sum_reuse)),
        "entropy_cost": sum_entropy,
        "reuse_cost": sum_reuse,
        "final_components": int(x.size),
        "final_component_score": component_score(x, y, z),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=50000)
    parser.add_argument("--outdir", type=Path, default=Path("results-multiscale"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    stability = search_local_stability()
    checks = random_step_checks(args.samples, seed=20260807)
    cascades = [
        equal_branch_cascade(depth=12, branches=2, reuse_factor=1.0),
        equal_branch_cascade(depth=12, branches=1, reuse_factor=0.9116174203759786),
        equal_branch_cascade(depth=12, branches=2, reuse_factor=0.9116174203759786),
    ]
    result = {"local_stability_experimental": stability, "random_checks": checks, "cascades": cascades}
    (args.outdir / "multiscale_bellman.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    lines = [
        "# Multiscale transfer Bellman experiment",
        "",
        "## Exact refinement inequality",
        "",
        f"Random checks: `{checks['samples']}`",
        f"Maximum numerical bound violation: `{checks['max_bound_violation']:.3e}`",
        f"Maximum logarithmic violation: `{checks['max_log_violation']:.3e}`",
        "",
        "## Experimental local triad stability",
        "",
        f"A: `{stability['A']:.9f}`",
        f"B: `{stability['B']:.9f}`",
        f"minimum searched margin: `{stability['minimum_margin']:.3e}`",
        f"gamma*: `{stability['gamma']:.12f}`",
        f"convex holonomy cost: `{stability.get('holonomy_cost', 0.0):.12f}`",
        f"candidate reuse factor: `{stability.get('reuse_factor_candidate', 1.0):.12f}`",
        "",
        "The refinement/Bellman inequality and convex minimization formula are exact.",
        "The A,B local stability constants are numerical candidates, not interval-certified.",
        "",
        "## Model cascades",
        "",
    ]
    for c in cascades:
        lines += [
            f"- depth={c['depth']}, branches={c['branches']}, rho={c['reuse_factor']:.9f}: "
            f"product={c['product_exact']:.6e}, total cost={c['entropy_cost'] + c['reuse_cost']:.6f}"
        ]
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
