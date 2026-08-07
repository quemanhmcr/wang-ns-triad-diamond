from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar

from .helical import all_signs, diamond_metrics, edge_metrics


def vectors_from_params(x: np.ndarray):
    log_rb, log_rc, theta, phi, psi = x
    rb, rc = math.exp(log_rb), math.exp(log_rc)
    a = np.array([1.0, 0.0, 0.0])
    b = rb * np.array([math.cos(theta), math.sin(theta), 0.0])
    c = rc * np.array([
        math.sin(phi) * math.cos(psi),
        math.sin(phi) * math.sin(psi),
        math.cos(phi),
    ])
    return a, b, c


def single_symmetric_objective(x: float) -> float:
    if not (0.500001 < x < 0.999999):
        return 1e6
    val = math.sqrt(4*x*x - 1)/(4*math.sqrt(2)*x) * math.log(1/x)
    return -val


def single_reference() -> Dict[str, float]:
    r = minimize_scalar(single_symmetric_objective, bounds=(0.500001, 0.999999), method="bounded")
    return {"x": float(r.x), "J": float(-r.fun)}


def objective(params: np.ndarray, signs, lambda_phase: float, j_ref: float) -> float:
    a, b, c = vectors_from_params(params)
    try:
        m = diamond_metrics(a, b, c, signs)
    except (ValueError, FloatingPointError):
        return 1e4
    edges = list(m["edges"].values())
    ratios = np.array([e.efficiency / j_ref for e in edges])
    # Require genuine forward progress on every edge.
    if any(e.forward_ratio <= 1.0005 for e in edges):
        return 100.0 + sum(max(0.0, 1.0005-e.forward_ratio) for e in edges)
    soft_min = -math.log(np.sum(np.exp(-20.0 * ratios))) / 20.0
    phase_penalty = lambda_phase * (m["phase_frustration"] / math.pi) ** 2
    # Mild regularizer against extreme scale ratios.
    scale_penalty = 0.002 * (params[0] ** 2 + params[1] ** 2)
    return -soft_min + phase_penalty + scale_penalty


def serialise_candidate(params, signs, lambda_phase, j_ref):
    a, b, c = vectors_from_params(np.asarray(params))
    m = diamond_metrics(a, b, c, signs)
    return {
        "params": [float(v) for v in params],
        "signs": [int(v) for v in signs],
        "lambda_phase": float(lambda_phase),
        "min_edge_ratio": float(m["min_efficiency"] / j_ref),
        "geom_mean_edge_ratio": float(m["geom_mean_efficiency"] / j_ref),
        "mean_edge_ratio": float(m["mean_efficiency"] / j_ref),
        "phase_frustration": float(m["phase_frustration"]),
        "geom_holonomy": float(m["geom_holonomy"]),
        "target_holonomy": float(m["target_holonomy"]),
        "vectors": {k: [float(t) for t in v] for k, v in m["vectors"].items()},
        "edges": {k: asdict(v) for k, v in m["edges"].items()},
    }


def run(quick: bool, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    ref = single_reference()
    j_ref = ref["J"]
    bounds = [
        (math.log(0.35), math.log(2.5)),
        (math.log(0.35), math.log(2.5)),
        (0.08, math.pi-0.08),
        (0.08, math.pi-0.08),
        (0.0, 2*math.pi),
    ]
    lambdas = [0.0, 1.0, 10.0, 100.0]
    maxiter = 80 if quick else 280
    popsize = 10 if quick else 20
    workers = 1
    candidates: List[Dict] = []
    sign_patterns = list(all_signs())

    for li, lam in enumerate(lambdas):
        for si, signs in enumerate(sign_patterns):
            seed = 1729 + 1000*li + si
            result = differential_evolution(
                objective,
                bounds=bounds,
                args=(signs, lam, j_ref),
                seed=seed,
                maxiter=maxiter,
                popsize=popsize,
                tol=2e-7 if not quick else 2e-5,
                polish=True,
                updating="immediate",
                workers=workers,
            )
            cand = serialise_candidate(result.x, signs, lam, j_ref)
            cand["optimizer_fun"] = float(result.fun)
            cand["success"] = bool(result.success)
            candidates.append(cand)
            print(json.dumps({
                "lambda": lam,
                "sign_index": si,
                "min_ratio": cand["min_edge_ratio"],
                "phase": cand["phase_frustration"],
            }), flush=True)

    candidates.sort(key=lambda c: (-c["min_edge_ratio"], c["phase_frustration"]))
    best_by_lambda = {}
    for lam in lambdas:
        subset = [c for c in candidates if c["lambda_phase"] == lam]
        # rank by the actual penalised score used conceptually
        subset.sort(key=lambda c: (-(c["min_edge_ratio"] - lam*(c["phase_frustration"]/math.pi)**2), c["phase_frustration"]))
        best_by_lambda[str(lam)] = subset[:10]

    feasible = {}
    for threshold in [0.5, 0.25, 0.1, 0.05, 0.02, 0.01]:
        subset = [c for c in candidates if c["phase_frustration"] <= threshold]
        feasible[str(threshold)] = max(subset, key=lambda c: c["min_edge_ratio"], default=None)

    payload = {
        "single_symmetric_reference": ref,
        "best_by_lambda": best_by_lambda,
        "best_under_phase_threshold": feasible,
        "top_unconstrained": candidates[:30],
        "quick": quick,
    }
    (outdir / "results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    lines = [
        "# Triad diamond optimisation results",
        "",
        f"Single-triad symmetric reference: `J={j_ref:.12f}` at `x={ref['x']:.12f}`.",
        "",
        "## Best candidates under phase-frustration thresholds",
        "",
        "| phase threshold | min edge / reference | actual phase | signs |",
        "|---:|---:|---:|:---|",
    ]
    for threshold, c in feasible.items():
        if c is None:
            lines.append(f"| {threshold} | none | none | none |")
        else:
            lines.append(f"| {threshold} | {c['min_edge_ratio']:.8f} | {c['phase_frustration']:.8g} | `{c['signs']}` |")
    lines += ["", "## Interpretation guardrail", "", "This is a finite-dimensional numerical probe, not a theorem about the PDE. A numerical gap is evidence for a candidate rigidity lemma only after interval-arithmetic or analytic certification."]
    (outdir / "summary.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--quick", action="store_true")
    p.add_argument("--outdir", default="results")
    args = p.parse_args()
    run(args.quick, Path(args.outdir))
