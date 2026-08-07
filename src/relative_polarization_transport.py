from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

J2 = np.array([[0.0, 1.0], [-1.0, 0.0]])
LAMBDA_CHILD = np.array([1.0, -1.0])


def tracefree_symmetric(delta: float, beta: float) -> np.ndarray:
    return np.array([[delta, beta], [beta, -delta]], dtype=float)


def sl2_symmetric_step(D: np.ndarray, dt: float) -> np.ndarray:
    """Exact exp(-dt D) for a real symmetric trace-free 2x2 matrix."""
    D = np.asarray(D, float)
    if np.linalg.norm(D - D.T) > 1e-12 or abs(float(np.trace(D))) > 1e-12:
        raise ValueError("D must be symmetric trace free")
    rho2 = 0.5 * float(np.sum(D * D))
    if rho2 < 1e-28:
        return np.eye(2) - dt * D
    rho = math.sqrt(rho2)
    x = rho * dt
    return math.cosh(x) * np.eye(2) - (math.sinh(x) / rho) * D


def parent_wedge(U: np.ndarray, V: np.ndarray) -> complex:
    return complex(np.asarray(U, complex).T @ J2 @ np.asarray(V, complex))


def child_factor(Z: np.ndarray) -> complex:
    return complex(LAMBDA_CHILD.T @ np.asarray(Z, complex))


def polarization_numerator(U: np.ndarray, V: np.ndarray, Z: np.ndarray) -> complex:
    return parent_wedge(U, V) * child_factor(Z)


def wedge_rhs(U: np.ndarray, V: np.ndarray, D1: np.ndarray, D2: np.ndarray) -> complex:
    """Exact d/dt(U^T J V) for Udot=-D1 U, Vdot=-D2 V."""
    return complex(np.asarray(U, complex).T @ J2 @ (np.asarray(D1, float) - np.asarray(D2, float)) @ np.asarray(V, complex))


def polarization_rhs(U: np.ndarray, V: np.ndarray, Z: np.ndarray, D1: np.ndarray, D2: np.ndarray, D3: np.ndarray) -> complex:
    W = parent_wedge(U, V)
    L = child_factor(Z)
    rel = wedge_rhs(U, V, D1, D2)
    child = complex(LAMBDA_CHILD.T @ np.asarray(D3, float) @ np.asarray(Z, complex))
    return rel * L - W * child


def pointwise_capacity_bound(U: np.ndarray, V: np.ndarray, Z: np.ndarray, D1: np.ndarray, D2: np.ndarray, D3: np.ndarray) -> float:
    """2 sqrt(||D1-D2||_F^2+||D3||_F^2) ||U||||V||||Z||."""
    q = float(np.sum((np.asarray(D1) - np.asarray(D2)) ** 2) + np.sum(np.asarray(D3) ** 2))
    amp = float(np.linalg.norm(U) * np.linalg.norm(V) * np.linalg.norm(Z))
    return 2.0 * math.sqrt(max(0.0, q)) * amp


def evolve_piecewise(U: np.ndarray, V: np.ndarray, Z: np.ndarray, rows: list[tuple[float, np.ndarray, np.ndarray, np.ndarray]]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    U = np.asarray(U, complex).copy(); V = np.asarray(V, complex).copy(); Z = np.asarray(Z, complex).copy()
    for dt, D1, D2, D3 in rows:
        U = sl2_symmetric_step(D1, dt) @ U
        V = sl2_symmetric_step(D2, dt) @ V
        Z = sl2_symmetric_step(D3, dt) @ Z
    return U, V, Z


def common_hyperbolic_countermodel(a: float = 8.0) -> dict[str, float]:
    """Large common SL(2) motion with exactly invariant parent wedge."""
    D = tracefree_symmetric(1.0, 0.0)
    M = sl2_symmetric_step(D, a)
    U = np.array([1.0 + 0.2j, -0.4 + 0.7j])
    V = np.array([0.3 - 0.1j, 0.8 + 0.5j])
    w0 = parent_wedge(U, V)
    w1 = parent_wedge(M @ U, M @ V)
    return {
        "a": a,
        "propagator_distance": float(np.linalg.norm(M - np.eye(2))),
        "condition_number": float(np.linalg.cond(M)),
        "wedge_relative_residual": float(abs(w1 - w0) / max(1.0, abs(w0))),
    }


@dataclass(frozen=True)
class TransportStress:
    samples: int
    worst_wedge_rhs_residual: float
    worst_polarization_rhs_residual: float
    worst_pointwise_bound_margin: float
    worst_common_timeordered_wedge_residual: float


def stress(samples: int = 50_000, seed: int = 20260807) -> TransportStress:
    rng = np.random.default_rng(seed)
    ww = wp = wc = 0.0
    margin = float("inf")
    eps = 2e-7
    for _ in range(samples):
        U = rng.normal(size=2) + 1j * rng.normal(size=2)
        V = rng.normal(size=2) + 1j * rng.normal(size=2)
        Z = rng.normal(size=2) + 1j * rng.normal(size=2)
        D1 = tracefree_symmetric(*rng.normal(size=2))
        D2 = tracefree_symmetric(*rng.normal(size=2))
        D3 = tracefree_symmetric(*rng.normal(size=2))

        # Differential identities by centered finite difference of exact steps.
        Up = sl2_symmetric_step(D1, eps) @ U
        Vp = sl2_symmetric_step(D2, eps) @ V
        Zp = sl2_symmetric_step(D3, eps) @ Z
        Um = sl2_symmetric_step(D1, -eps) @ U
        Vm = sl2_symmetric_step(D2, -eps) @ V
        Zm = sl2_symmetric_step(D3, -eps) @ Z
        numw = (parent_wedge(Up, Vp) - parent_wedge(Um, Vm)) / (2 * eps)
        anaw = wedge_rhs(U, V, D1, D2)
        ww = max(ww, abs(numw - anaw) / max(1.0, abs(anaw)))
        nump = (polarization_numerator(Up, Vp, Zp) - polarization_numerator(Um, Vm, Zm)) / (2 * eps)
        anap = polarization_rhs(U, V, Z, D1, D2, D3)
        wp = max(wp, abs(nump - anap) / max(1.0, abs(anap)))
        bnd = pointwise_capacity_bound(U, V, Z, D1, D2, D3)
        margin = min(margin, bnd - abs(anap))
        if abs(anap) > bnd + 2e-10:
            raise AssertionError("pointwise capacity bound violated")

        # Arbitrary noncommuting common time-ordered SL(2) history is neutral.
        rows = []
        for _j in range(4):
            D = tracefree_symmetric(*rng.normal(scale=1.4, size=2))
            rows.append((float(rng.uniform(0.01, 0.15)), D, D, np.zeros((2, 2))))
        U1, V1, _ = evolve_piecewise(U, V, Z, rows)
        wc = max(wc, abs(parent_wedge(U1, V1) - parent_wedge(U, V)) / max(1.0, abs(parent_wedge(U, V))))
    return TransportStress(samples, ww, wp, margin, wc)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-relative-polarization"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cm = common_hyperbolic_countermodel()
    payload = {"stress": out.__dict__, "common_hyperbolic_countermodel": cm}
    (args.outdir / "relative_polarization_transport.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Exact relative-polarization transport

The identities are exact; random checks validate the implementation.

- exact parent wedge law: `Wdot = U^T J (D1-D2) V`
- exact polarization numerator law includes only relative parent generator and child generator
- pointwise capacity bound: `|Pdot| <= 2 sqrt(||D1-D2||_F^2+||D3||_F^2)||U||||V||||Z||`
- random differential/time-ordered checks: `{out.samples}`
- worst wedge RHS residual: `{out.worst_wedge_rhs_residual:.3e}`
- worst polarization RHS residual: `{out.worst_polarization_rhs_residual:.3e}`
- minimum pointwise bound margin: `{out.worst_pointwise_bound_margin:.3e}`
- worst arbitrary common time-ordered wedge residual: `{out.worst_common_timeordered_wedge_residual:.3e}`
- hyperbolic common-gauge countermodel: `||M-I||={cm['propagator_distance']:.3e}`, `cond(M)={cm['condition_number']:.3e}`, wedge residual `{cm['wedge_relative_residual']:.3e}`

Thus the full time-ordered parent observable does not require a Magnus expansion:
common `SL(2)` motion cancels pointwise in the symplectic wedge.  Euclidean
propagator distance from the identity is not a physical polarization defect.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
