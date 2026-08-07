from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from .affine_grain_dynamics import extremal_parent_directions, sym, tracefree_2x2
from .single_edge_certificate import float_rstar

INTRINSIC_COERCIVITY = 43.0 / 100.0


def gram_rhs(K: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Exact intrinsic Gram evolution for K_dot=-B K.

    K is 3x2 with parent carriers as columns.  G=K^T K obeys
      G_dot=-K^T(B+B^T)K.
    All Euclidean triangle side lengths are functions of G alone.
    """
    K = np.asarray(K, dtype=float)
    B = np.asarray(B, dtype=float)
    if K.shape != (3, 2) or B.shape != (3, 3):
        raise ValueError("expected K 3x2 and B 3x3")
    return -(K.T @ (B + B.T) @ K)


def shape_from_gram(G: np.ndarray, rstar: float | None = None) -> tuple[float, float, float]:
    """Signed (u,v,H) from the two-parent Gram matrix."""
    G = np.asarray(G, dtype=float)
    if G.shape != (2, 2):
        raise ValueError("expected a 2x2 Gram matrix")
    if rstar is None:
        rstar = float_rstar()
    aa, bb = float(G[0, 0]), float(G[1, 1])
    cc = float(G[0, 0] + 2 * G[0, 1] + G[1, 1])
    if min(aa, bb, cc) <= 0:
        raise ValueError("degenerate triangle")
    la, lb, lc = 0.5 * math.log(aa), 0.5 * math.log(bb), 0.5 * math.log(cc)
    u = lb - la
    v = lc - 0.5 * (la + lb) + math.log(rstar)
    H = 0.5 * u * u + 2.0 * v * v
    return u, v, H


def shape_rates_from_gram(G: np.ndarray, Gdot: np.ndarray) -> tuple[float, float, float]:
    G = np.asarray(G, dtype=float)
    Gdot = np.asarray(Gdot, dtype=float)
    aa, bb = float(G[0, 0]), float(G[1, 1])
    cc = float(G[0, 0] + 2 * G[0, 1] + G[1, 1])
    aad, bbd = float(Gdot[0, 0]), float(Gdot[1, 1])
    ccd = float(Gdot[0, 0] + 2 * Gdot[0, 1] + Gdot[1, 1])
    ra, rb, rc = 0.5 * aad / aa, 0.5 * bbd / bb, 0.5 * ccd / cc
    udot = rb - ra
    vdot = rc - 0.5 * (ra + rb)
    speed2 = 0.5 * udot * udot + 2.0 * vdot * vdot
    return udot, vdot, speed2


def oriented_extremal_parents(E: np.ndarray, rstar: float | None = None) -> np.ndarray:
    """3x2 extremal parent matrix in an arbitrary orthonormal plane basis E."""
    E = np.asarray(E, dtype=float)
    if E.shape != (3, 2) or not np.allclose(E.T @ E, np.eye(2), atol=1e-12):
        raise ValueError("E must be a 3x2 orthonormal frame")
    if rstar is None:
        rstar = float_rstar()
    na, nb, _ = extremal_parent_directions(rstar)
    return E @ np.column_stack((rstar * na, rstar * nb))


def restricted_shape_driver(B: np.ndarray, E: np.ndarray) -> np.ndarray:
    """Trace-free symmetric restriction of the common carrier driver to plane E."""
    B = np.asarray(B, dtype=float)
    E = np.asarray(E, dtype=float)
    H = E.T @ sym(B) @ E
    return tracefree_2x2(H)


def intrinsic_extremal_shape_rates(B: np.ndarray, E: np.ndarray, rstar: float | None = None) -> tuple[float, float, float]:
    K = oriented_extremal_parents(E, rstar)
    G = K.T @ K
    Gdot = gram_rhs(K, B)
    return shape_rates_from_gram(G, Gdot)


def intrinsic_coercivity_ratio(B: np.ndarray, E: np.ndarray, rstar: float | None = None) -> float:
    D = restricted_shape_driver(B, E)
    denom = float(np.sum(D * D))
    if denom == 0.0:
        return math.inf
    return intrinsic_extremal_shape_rates(B, E, rstar)[2] / denom


def plane_normal_rhs(n: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Exact normal evolution for the carrier plane under K_dot=-B K.

    If n.K=0, then n_dot=(I-nn^T) B^T?  Here carrier equation is Kdot=-B K;
    differentiating n^T K=0 gives n_dot-B^T n parallel n, so after normalization
      n_dot=(I-nn^T) B^T n.
    """
    n = np.asarray(n, dtype=float)
    B = np.asarray(B, dtype=float)
    n = n / np.linalg.norm(n)
    return (np.eye(3) - np.outer(n, n)) @ B.T @ n


def gaussian_effective_carrier_driver(A: np.ndarray, P: np.ndarray, nu: float) -> np.ndarray:
    """B in kappa_dot=-B kappa for a common Gaussian precision P."""
    A = np.asarray(A, dtype=float)
    P = np.asarray(P, dtype=float)
    if nu < 0:
        raise ValueError("nu must be nonnegative")
    return A.T + 2.0 * nu * np.linalg.inv(P)


def random_orthonormal_plane(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q[:, :2], Q[:, 2]


def stress(samples: int = 50_000, seed: int = 20260807) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    worst_ratio = float("inf")
    worst_gram = 0.0
    worst_shape_rate = 0.0
    worst_normal_constraint = 0.0
    worst_isotropic_viscous_intrinsic = 0.0
    r = float_rstar()

    for _ in range(samples):
        E, n = random_orthonormal_plane(rng)
        A = rng.normal(size=(3, 3))
        A -= np.trace(A) / 3 * np.eye(3)
        P0 = rng.normal(size=(3, 3)); P = P0 @ P0.T + 0.5 * np.eye(3)
        nu = float(rng.uniform(0.0, 1.5))
        B = gaussian_effective_carrier_driver(A, P, nu)

        K = oriented_extremal_parents(E, r)
        G = K.T @ K
        Gdot = gram_rhs(K, B)
        Kdot = -B @ K
        direct = Kdot.T @ K + K.T @ Kdot
        worst_gram = max(worst_gram, np.linalg.norm(Gdot - direct) / max(1.0, np.linalg.norm(direct)))

        ud, vd, sp2 = shape_rates_from_gram(G, Gdot)
        u2, v2, sp22 = intrinsic_extremal_shape_rates(B, E, r)
        worst_shape_rate = max(worst_shape_rate, abs(ud-u2), abs(vd-v2), abs(sp2-sp22))
        ratio = intrinsic_coercivity_ratio(B, E, r)
        if math.isfinite(ratio):
            worst_ratio = min(worst_ratio, ratio)
            if ratio + 2e-11 < INTRINSIC_COERCIVITY:
                raise AssertionError(("full-3D intrinsic coercivity failed", ratio))

        ndot = plane_normal_rhs(n, B)
        # Differentiate n^T K=0 exactly.
        constraint = ndot @ K + n @ Kdot
        worst_normal_constraint = max(worst_normal_constraint, float(np.linalg.norm(constraint)))

        # Isotropic viscosity is a scalar B contribution and vanishes from the
        # trace-free restricted driver in every plane orientation.
        p = float(rng.uniform(0.2, 4.0))
        Bnu = A.T + 2.0 * nu / p * np.eye(3)
        Dnu = restricted_shape_driver(Bnu, E)
        D0 = restricted_shape_driver(A.T, E)
        worst_isotropic_viscous_intrinsic = max(worst_isotropic_viscous_intrinsic, float(np.linalg.norm(Dnu-D0)))

    return {
        "samples": samples,
        "worst_full3d_intrinsic_coercivity": worst_ratio,
        "worst_gram_identity_residual": worst_gram,
        "worst_shape_rate_residual": worst_shape_rate,
        "worst_plane_normal_constraint_residual": worst_normal_constraint,
        "worst_isotropic_viscosity_intrinsic_residual": worst_isotropic_viscous_intrinsic,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=50_000)
    ap.add_argument('--outdir', type=Path, default=Path('results-intrinsic-plane'))
    args = ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    result=stress(args.samples)
    (args.outdir/'intrinsic_triad_plane.json').write_text(json.dumps(result,indent=2))
    md=f"""# Intrinsic 3D triad-plane dynamics\n\n- exact Gram law: `Gdot=-K^T(B+B^T)K`\n- scalar multiplier geometry depends only on the intrinsic Gram matrix\n- full-3D instantaneous extremal coercivity: `>=43/100` times the squared\n  trace-free symmetric restriction of the carrier driver to the evolving plane\n- random 3D plane/strain/Gaussian checks: `{result['samples']}`\n- worst full-3D coercivity seen: `{result['worst_full3d_intrinsic_coercivity']:.9f}`\n- worst Gram residual: `{result['worst_gram_identity_residual']:.3e}`\n- worst shape-rate residual: `{result['worst_shape_rate_residual']:.3e}`\n- worst plane-normal constraint residual: `{result['worst_plane_normal_constraint_residual']:.3e}`\n- worst isotropic-viscosity intrinsic residual: `{result['worst_isotropic_viscosity_intrinsic_residual']:.3e}`\n\nExtrinsic tilt of the common triad plane is therefore a gauge for scalar\nside-length/multiplier geometry.  The remaining genuinely three-dimensional\nissue is coherent transport of helical polarization and of the spatial packet\nframe, not the existence of a fixed Euclidean plane.\n"""
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
