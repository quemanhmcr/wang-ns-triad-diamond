from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


def symmetrize_rank3(B: np.ndarray) -> np.ndarray:
    B = np.asarray(B, dtype=float)
    if B.shape != (3, 3, 3):
        raise ValueError("B must be 3x3x3")
    out = np.zeros_like(B)
    for perm in itertools.permutations(range(3)):
        out += np.transpose(B, perm)
    return out / 6.0


def whitened_velocity_hessian(L: np.ndarray, H: np.ndarray) -> np.ndarray:
    """B=L^{-1} H[L.,L.] in grain coordinates.

    H[i,j,k]=partial_j partial_k U_i, symmetric in j,k.  L is any invertible
    grain frame with physical covariance Sigma=L L^T.
    """
    L = np.asarray(L, dtype=float)
    H = np.asarray(H, dtype=float)
    if L.shape != (3, 3) or H.shape != (3, 3, 3):
        raise ValueError("wrong shapes")
    Linv = np.linalg.inv(L)
    return np.einsum("ai,ijk,jb,kc->abc", Linv, H, L, L)


def normalized_carrier(L: np.ndarray, k: np.ndarray) -> np.ndarray:
    return np.asarray(L, dtype=float).T @ np.asarray(k, dtype=float)


def quadratic_phase_tensor(B: np.ndarray, q: np.ndarray) -> np.ndarray:
    """C_bc=q_a B_abc, the quadratic wavefront/chirp tangent."""
    return np.einsum("a,abc->bc", np.asarray(q, float), np.asarray(B, float))


def cubic_trace(T: np.ndarray) -> np.ndarray:
    return np.einsum("iik->k", np.asarray(T, float))


def full_quadratic_advection_residual_sq(B: np.ndarray, q: np.ndarray) -> float:
    """Exact ||R_2 . grad psi||_2^2 / ||psi||_2^2 for an unchirped Gaussian.

    psi(Lz)=const exp(-|z|^2/4) exp(i q.z), z~N(0,I) under |psi|^2, and
    R_2=(1/2)H[Lz,Lz].  The phase-quadratic and real cubic pieces are pointwise
    orthogonal, so Wick's formula is exact.
    """
    B = np.asarray(B, float)
    C = quadratic_phase_tensor(B, q)
    T = symmetrize_rank3(B)
    trC = float(np.trace(C))
    t = cubic_trace(T)
    phase = 0.25 * (2.0 * float(np.sum(C * C)) + trC * trC)
    cubic = (6.0 * float(np.sum(T * T)) + 9.0 * float(np.dot(t, t))) / 16.0
    return phase + cubic


def osculating_transverse_residual_sq(B: np.ndarray) -> float:
    """Exact residual after projecting to Gaussian tangent parameters.

    Quadratic phase q.B is absorbed by chirp/covariance.  The cubic polynomial
    T_abc z_a z_b z_c decomposes into its tangent linear piece 3 tr(T).z and the
    third Hermite chaos.  The latter has squared Gaussian norm 6||T||_F^2.
    """
    T = symmetrize_rank3(B)
    return 3.0 / 8.0 * float(np.sum(T * T))


def osculating_transverse_bound(B: np.ndarray) -> float:
    return math.sqrt(6.0) / 4.0 * float(np.linalg.norm(B))


def gaussian_laplacian_multiplier(G: np.ndarray, k: np.ndarray, y: np.ndarray) -> complex:
    """Exact multiplier (Delta psi)/psi for a complex/chirped Gaussian.

    psi(y)=exp(-1/2 y^T G y + i k.y), with symmetric complex G.  The result is
    a polynomial of degree at most two, hence tangent to the Gaussian manifold.
    """
    G=np.asarray(G,complex); k=np.asarray(k,float); y=np.asarray(y,float)
    return complex(y @ (G @ G) @ y - 2j*(k @ G @ y) - np.dot(k,k) - np.trace(G))


def _third_hermite_value_with_trace(T: np.ndarray, trace_T: np.ndarray, z: np.ndarray) -> float:
    """Evaluate the same H3 polynomial with the tensor trace already computed."""
    p = float(np.einsum("abc,a,b,c", T, z, z, z))
    return p - 3.0 * float(np.dot(trace_T, z))


def third_hermite_value(T: np.ndarray, z: np.ndarray) -> float:
    T = np.asarray(T, float)
    z = np.asarray(z, float)
    return _third_hermite_value_with_trace(T, cubic_trace(T), z)


def transform_hessian(S: np.ndarray, H: np.ndarray) -> np.ndarray:
    """Hessian under x'=Sx, U'(x')=S U(S^{-1}x')."""
    S = np.asarray(S, float)
    Sinv = np.linalg.inv(S)
    return np.einsum("ia,abc,bj,ck->ijk", S, np.asarray(H, float), Sinv, Sinv)


def random_invertible(rng: np.random.Generator, log_condition: float = 8.0) -> np.ndarray:
    q1, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q2, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    logs = rng.uniform(-0.5 * log_condition, 0.5 * log_condition, size=3)
    return q1 @ np.diag(np.exp(logs)) @ q2.T


@dataclass(frozen=True)
class ForcingStress:
    samples: int
    worst_affine_B_residual: float
    worst_affine_q_residual: float
    worst_orthogonal_gauge_residual: float
    worst_hermite_moment_relative_error: float
    worst_transverse_bound_ratio: float
    extreme_frame_condition: float
    extreme_affine_residual: float


def stress(samples: int = 50_000, seed: int = 20260807) -> ForcingStress:
    rng = np.random.default_rng(seed)
    waB = waq = wog = wh = wt = 0.0
    hermite_checks = min(500, samples)
    for n in range(samples):
        L = random_invertible(rng, log_condition=6.0)
        H = rng.normal(size=(3, 3, 3))
        H = 0.5 * (H + H.swapaxes(1, 2))
        k = rng.normal(size=3)
        B = whitened_velocity_hessian(L, H)
        q = normalized_carrier(L, k)

        S = random_invertible(rng, log_condition=5.0)
        Lp = S @ L
        Hp = transform_hessian(S, H)
        kp = np.linalg.solve(S.T, k)
        Bp = whitened_velocity_hessian(Lp, Hp)
        qp = normalized_carrier(Lp, kp)
        waB = max(waB, float(np.linalg.norm(Bp - B)) / max(1.0, float(np.linalg.norm(B))))
        waq = max(waq, float(np.linalg.norm(qp - q)) / max(1.0, float(np.linalg.norm(q))))

        O, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        Bg = whitened_velocity_hessian(L @ O, H)
        qg = normalized_carrier(L @ O, k)
        # Coordinate change z=O z' gives B'_{abc}=O_{pa} B_{pqr} O_{qb} O_{rc}.
        Bexpected = np.einsum("pa,pqr,qb,rc->abc", O, B, O, O)
        qexpected = O.T @ q
        wog = max(wog, float(np.linalg.norm(Bg - Bexpected)), float(np.linalg.norm(qg - qexpected)))

        exact = math.sqrt(osculating_transverse_residual_sq(B))
        bound = osculating_transverse_bound(B)
        if exact > bound + 2e-10:
            raise AssertionError("osculating transverse bound failed")
        if bound > 1e-14:
            wt = max(wt, exact / bound)

        if n < hermite_checks:
            T = symmetrize_rank3(B)
            trace_T = cubic_trace(T)
            # Same 3000-point Monte Carlo check; only the T-invariant trace is hoisted.
            zz = rng.normal(size=(3000, 3))
            vals = np.array([_third_hermite_value_with_trace(T, trace_T, z) for z in zz])
            empirical = float(np.mean(vals * vals))
            analytic = 6.0 * float(np.sum(T * T))
            if analytic > 1e-12:
                wh = max(wh, abs(empirical - analytic) / analytic)

    # Explicit huge-condition common affine coordinate transform.  The physical
    # tensor and grain transform together, so the intrinsic forcing is unchanged.
    A = 1.0e5
    S = np.diag([A, 1.0, 1.0 / A])
    L = np.diag([1.2, 0.8, 1.1])
    H = rng.normal(size=(3, 3, 3)); H = 0.5 * (H + H.swapaxes(1, 2))
    k = rng.normal(size=3)
    B = whitened_velocity_hessian(L, H); q = normalized_carrier(L, k)
    Bp = whitened_velocity_hessian(S @ L, transform_hessian(S, H))
    qp = normalized_carrier(S @ L, np.linalg.solve(S.T, k))
    extreme = max(float(np.linalg.norm(Bp - B)), float(np.linalg.norm(qp - q)))
    return ForcingStress(
        samples=samples,
        worst_affine_B_residual=waB,
        worst_affine_q_residual=waq,
        worst_orthogonal_gauge_residual=wog,
        worst_hermite_moment_relative_error=wh,
        worst_transverse_bound_ratio=wt,
        extreme_frame_condition=float(np.linalg.cond(S @ L)),
        extreme_affine_residual=extreme,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-affine-gaussian-forcing"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    payload = {
        "theorem": {
            "whitened_hessian": "B=L^-1 H[L.,L.]",
            "normalized_carrier": "q=L^T k",
            "full_residual_sq": "1/4(2||q.B||_F^2+tr(q.B)^2)+1/16(6||Sym B||_F^2+9||tr Sym B||^2)",
            "osculating_transverse_sq": "3/8 ||Sym B||_F^2",
            "osculating_bound": "sqrt(6)/4 ||B||_F",
            "interpretation": "quadratic phase is chirp tangent; cubic trace is center tangent; first transverse mode is third Hermite chaos",
        },
        "stress": asdict(out),
    }
    (args.outdir / "affine_gaussian_forcing.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Affine-covariant Gaussian forcing and Hermite projection

- intrinsic curvature tensor: `B=L^-1 H[L.,L.]`
- intrinsic carrier: `q=L^T k`
- exact unprojected quadratic-advection residual is the Wick formula recorded in the JSON
- after allowing Gaussian center/carrier/covariance/chirp to osculate the flow,
  the first transverse forcing is exactly third Hermite chaos:
  `||F_perp||/||psi|| = sqrt(6)/4 ||Sym B||_F <= sqrt(6)/4 ||B||_F`
- random affine checks: `{out.samples}`
- worst affine B residual: `{out.worst_affine_B_residual:.3e}`
- worst affine q residual: `{out.worst_affine_q_residual:.3e}`
- worst orthogonal grain-gauge residual: `{out.worst_orthogonal_gauge_residual:.3e}`
- worst sampled transverse/bound ratio: `{out.worst_transverse_bound_ratio:.9f}`
- extreme transformed grain condition number: `{out.extreme_frame_condition:.3e}`
- extreme affine-invariance residual: `{out.extreme_affine_residual:.3e}`

A large Euclidean aspect ratio is therefore not itself a forcing cost.  The
physical curvature variable is the Hessian expressed in the grain's own affine
metric.  Quadratic wavefront curvature belongs to the Gaussian tangent manifold;
calling it residual forcing would double-count a packet degree of freedom.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
