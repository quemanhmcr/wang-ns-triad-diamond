from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

from src.helical_spin_transport import forward_normal_coupling
from src.full_strain_observability import extremal_parent_directions, restriction_tracefree_norm2, tracefree_2x2, transverse_frame

RSTAR_LO = Fraction(61090410158, 100_000_000_000)
RSTAR_HI = Fraction(61090410160, 100_000_000_000)
J2 = np.array([[0.0, 1.0], [-1.0, 0.0]])


def isosceles_geometry(r: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not (0.5 < r < 10.0):
        raise ValueError("need a nondegenerate isosceles parent ratio r>1/2")
    c = 1.0 / (2.0 * r)
    ss = math.sqrt(1.0 - c * c)
    x = r * np.array([c, ss, 0.0])
    y = r * np.array([c, -ss, 0.0])
    z = x + y
    return x, y, z


def child_energy_helicity_tensor(r: float) -> np.ndarray:
    """Tensor A*g on an isosceles forward triad, indices (+,-) for each role."""
    x, y, z = isosceles_geometry(r)
    signs = (1, -1)
    T = np.zeros((2, 2, 2), dtype=complex)
    for i, sx in enumerate(signs):
        for j, sy in enumerate(signs):
            for k, sz in enumerate(signs):
                A = sx * r - sy * r
                T[i, j, k] = A * forward_normal_coupling(x, y, z, sx, sy, sz)
    return T


def tensor_factor_coefficient(r: float) -> float:
    """Positive C(r) in T_{s1 s2 s3}=-i C s3 epsilon_{s1s2}."""
    x, y, z = isosceles_geometry(r)
    area = 0.5 * np.linalg.norm(np.cross(x, y))
    return float(area / (math.sqrt(2.0) * r))


def factorized_tensor(r: float) -> np.ndarray:
    C = tensor_factor_coefficient(r)
    signs = (1, -1)
    T = np.zeros((2, 2, 2), dtype=complex)
    for i, _sx in enumerate(signs):
        for j, _sy in enumerate(signs):
            eps = J2[i, j]
            for k, sz in enumerate(signs):
                T[i, j, k] = -1j * C * sz * eps
    return T


def parent_wedge(u: np.ndarray, v: np.ndarray) -> complex:
    return complex(np.asarray(u, complex).T @ J2 @ np.asarray(v, complex))


def relative_parent_wedge(u: np.ndarray, v: np.ndarray, M1: np.ndarray, M2: np.ndarray) -> complex:
    return parent_wedge(np.asarray(M1, complex) @ u, np.asarray(M2, complex) @ v)


def relative_matrix_formula(u: np.ndarray, v: np.ndarray, M1: np.ndarray, M2: np.ndarray) -> complex:
    M1 = np.asarray(M1, complex); M2 = np.asarray(M2, complex)
    if abs(np.linalg.det(M1) - 1.0) > 1e-9:
        raise ValueError("M1 must have determinant one")
    R = np.linalg.inv(M1) @ M2
    return complex(np.asarray(u, complex).T @ J2 @ R @ np.asarray(v, complex))


@lru_cache(maxsize=32)
def _transfer_observability_frames(rstar: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Immutable transverse frames for one fixed physical parent ratio."""
    k1, k2 = extremal_parent_directions(float(rstar))
    Epi = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    E1 = transverse_frame(k1)
    E2 = transverse_frame(k2)
    E3 = transverse_frame(np.array([1.0, 0.0, 0.0]))
    for E in (Epi, E1, E2, E3):
        E.setflags(write=False)
    return Epi, E1, E2, E3


def _transfer_relevant_strain_observability_validated(
    S: np.ndarray, rstar: float = 0.610904101586766
) -> tuple[float, float]:
    """Evaluate Q_rel after symmetry/trace provenance was already validated."""
    Epi, E1, E2, E3 = _transfer_observability_frames(float(rstar))
    Dp = tracefree_2x2(Epi.T @ S @ Epi)
    D1 = tracefree_2x2(E1.T @ S @ E1)
    D2 = tracefree_2x2(E2.T @ S @ E2)
    D3 = tracefree_2x2(E3.T @ S @ E3)
    Q = float(np.sum(Dp * Dp) + np.sum((D1 - D2) ** 2) + np.sum(D3 * D3))
    return Q, float(np.sum(S * S))


def _transfer_relevant_strain_observability_validated_batch(
    S: np.ndarray, rstar: float = 0.610904101586766
) -> tuple[np.ndarray, np.ndarray]:
    """Bitwise-equivalent batched Q_rel on already-validated trace-free strains."""
    Epi, E1, E2, E3 = _transfer_observability_frames(float(rstar))
    Ds = []
    for E in (Epi, E1, E2, E3):
        M = np.matmul(np.matmul(E.T, S), E)
        tr = 0.5 * (M[:, 0, 0] + M[:, 1, 1])
        D = M.copy()
        D[:, 0, 0] -= tr
        D[:, 1, 1] -= tr
        Ds.append(D)
    Dp, D1, D2, D3 = Ds
    Q = (
        np.sum(Dp * Dp, axis=(1, 2))
        + np.sum((D1 - D2) ** 2, axis=(1, 2))
        + np.sum(D3 * D3, axis=(1, 2))
    )
    N = np.sum(S * S, axis=(1, 2))
    return Q, N


def transfer_relevant_strain_observability(S: np.ndarray, rstar: float = 0.610904101586766) -> tuple[float, float]:
    """Q_rel=||D_Pi||^2+||D1-D2||^2+||D_child||^2 and ||S||^2."""
    S = np.asarray(S, dtype=float)
    if np.linalg.norm(S - S.T) > 1e-10 or abs(float(np.trace(S))) > 1e-10:
        raise ValueError("S must be symmetric trace free")
    return _transfer_relevant_strain_observability_validated(S, rstar)


def arb_transfer_relevant_observability_certificate() -> dict[str, str]:
    """Certify Q_rel >= 1/2 ||S||^2 on the r* bracket.

    Put C=cos(phi)^2=1/(4r^2) and
      S=[[a,b,x],[b,d,y],[x,y,-a-d]].
    Direct expansion gives
      Q_rel-1/2||S||^2
       = (3/2)d^2 +(1+8C-8C^2)b^2 +(7-8C)x^2 + y^2.
    """
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    def aq(q: Fraction):
        return arb(q.numerator) / q.denominator
    r = aq(RSTAR_LO).union(aq(RSTAR_HI))
    C = 1 / (4 * r * r)
    bcoef = 1 + 8 * C - 8 * C * C
    xcoef = 7 - 8 * C
    if not (bcoef > arb(0) and xcoef > arb(0)):
        raise AssertionError(f"transfer-relevant strain coefficients failed: b={bcoef}, x={xcoef}")
    return {
        "rstar_ball": str(r),
        "cos2_half_angle_ball": str(C),
        "relative_observability_lower": "1/2",
        "remainder_identity": "(3/2)d^2+(1+8C-8C^2)b^2+(7-8C)x^2+y^2",
        "b_coefficient_ball": str(bcoef),
        "x_coefficient_ball": str(xcoef),
        "equality_mode": "S=diag(a,0,-a) in child-aligned coordinates",
        "status": "CERTIFIED",
    }


def random_sl2(rng: np.random.Generator) -> np.ndarray:
    # exponentiate a random real trace-free 2x2 matrix by eigen decomposition
    A = rng.normal(size=(2, 2))
    A -= 0.5 * np.trace(A) * np.eye(2)
    vals, vecs = np.linalg.eig(A)
    M = vecs @ np.diag(np.exp(vals)) @ np.linalg.inv(vecs)
    M = np.real_if_close(M, tol=1000).astype(float)
    # numerical determinant normalization
    M /= math.sqrt(abs(float(np.linalg.det(M))))
    if np.linalg.det(M) < 0:
        M[:, 0] *= -1
    return M


@dataclass(frozen=True)
class SymplecticStress:
    samples: int
    worst_tensor_factorization_residual: float
    worst_common_sl2_invariance_residual: float
    worst_relative_matrix_residual: float
    worst_transfer_observability_ratio: float


def stress(samples: int = 50_000, seed: int = 20260807) -> SymplecticStress:
    rng = np.random.default_rng(seed)
    wt = ws = wr = 0.0
    wo = float("inf")
    r = 0.610904101586766
    for _ in range(samples):
        rr = float(rng.uniform(0.52, 1.2))
        T = child_energy_helicity_tensor(rr)
        F = factorized_tensor(rr)
        wt = max(wt, float(np.linalg.norm(T - F)) / max(1.0, float(np.linalg.norm(F))))

        u = rng.normal(size=2) + 1j * rng.normal(size=2)
        v = rng.normal(size=2) + 1j * rng.normal(size=2)
        M = random_sl2(rng)
        before = parent_wedge(u, v)
        after = relative_parent_wedge(u, v, M, M)
        ws = max(ws, abs(after - before) / max(1.0, abs(before)))

        M2 = random_sl2(rng)
        direct = relative_parent_wedge(u, v, M, M2)
        rel = relative_matrix_formula(u, v, M, M2)
        wr = max(wr, abs(direct - rel) / max(1.0, abs(direct)))

        X = rng.normal(size=(3, 3))
        S = 0.5 * (X + X.T)
        S -= np.trace(S) / 3.0 * np.eye(3)
        Q, N = transfer_relevant_strain_observability(S, r)
        if N > 1e-16:
            wo = min(wo, Q / N)
            if Q + 2e-12 < 0.5 * N:
                raise AssertionError("transfer-relevant strain observability failed")
    return SymplecticStress(samples, wt, ws, wr, wo)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-helicity-symplectic"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_transfer_relevant_observability_certificate()
    out = stress(args.samples)
    (args.outdir / "extremal_helicity_symplectic.json").write_text(json.dumps({"certificate": cert, "stress": out.__dict__}, indent=2))
    md = f"""# Extremal helicity tensor and symplectic parent gauge

Status: **{cert['status']}** for the transfer-relevant strain observable.

- exact isosceles tensor: `T_(s1,s2,s3)=-i C(r) s3 epsilon_(s1,s2)`
- exact parent symmetry: `(M u)^T J (M v)=u^T J v` for every `M in SL(2)`
- exact relative formula: `(M1 u)^T J (M2 v)=u^T J (M1^-1 M2)v`
- transfer-relevant strain tomography:
  `||D_Pi||^2+||D1-D2||^2+||D_child||^2 >= 1/2 ||S||^2`
- exact positive remainder: `{cert['remainder_identity']}`
- random checks: `{out.samples}`
- worst tensor-factorization residual: `{out.worst_tensor_factorization_residual:.3e}`
- worst common-SL2 invariance residual: `{out.worst_common_sl2_invariance_residual:.3e}`
- worst relative-matrix residual: `{out.worst_relative_matrix_residual:.3e}`
- worst observed transfer-relevant ratio: `{out.worst_transfer_observability_ratio:.9f}`

This corrects a tempting but false interpretation: absolute helicity conversion
is not by itself a transfer cost.  At equal-parent geometry a common determinant-
one deformation of both parent helicity spinors is an exact symmetry of the
nonlinear parent wedge.  The physical variables are relative parent polarization,
child polarization, and scalar triad shape.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
