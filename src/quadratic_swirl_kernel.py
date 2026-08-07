from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.affine_gaussian_forcing import symmetrize_rank3


def levi_civita() -> np.ndarray:
    e = np.zeros((3, 3, 3))
    e[0, 1, 2] = e[1, 2, 0] = e[2, 0, 1] = 1.0
    e[1, 0, 2] = e[2, 1, 0] = e[0, 2, 1] = -1.0
    return e


EPS3 = levi_civita()


def swirl_tensor(M: np.ndarray) -> np.ndarray:
    """B_abc=eps_abd M_dc + eps_acd M_db, symmetric in b,c."""
    M = np.asarray(M, float)
    return np.einsum("abd,dc->abc", EPS3, M) + np.einsum("acd,db->abc", EPS3, M)


def swirl_velocity(M: np.ndarray, z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, float)
    return np.cross(z, np.asarray(M, float) @ z)


def quadratic_velocity(B: np.ndarray, z: np.ndarray) -> np.ndarray:
    return 0.5 * np.einsum("abc,b,c->a", np.asarray(B, float), z, z)


def reconstruct_M_from_kernel(B: np.ndarray) -> np.ndarray:
    """For a divergence-free Sym(B)=0 kernel, M_dc=(1/3)eps_abd B_abc."""
    return (1.0 / 3.0) * np.einsum("abd,abc->dc", EPS3, np.asarray(B, float))


def divergence_gradient(B: np.ndarray) -> np.ndarray:
    return np.einsum("aac->c", np.asarray(B, float))


def gaussian_envelope_advection(M: np.ndarray, z: np.ndarray) -> float:
    """V.grad exp(-|z|^2/4) divided by the envelope itself."""
    V = swirl_velocity(M, z)
    return -0.5 * float(np.dot(V, z))


def carrier_chirp_matrix(M: np.ndarray, q: np.ndarray) -> np.ndarray:
    """Symmetric C with q.(z x Mz)=z^T C z."""
    B = swirl_tensor(M)
    return 0.5 * np.einsum("a,abc->bc", np.asarray(q, float), B)


@dataclass(frozen=True)
class SwirlStress:
    samples: int
    worst_full_symmetry_residual: float
    worst_divergence_residual: float
    worst_reconstruction_residual: float
    worst_velocity_residual: float
    worst_radial_advection_residual: float
    worst_chirp_residual: float


def stress(samples: int = 50_000, seed: int = 20260807) -> SwirlStress:
    rng = np.random.default_rng(seed)
    ws = wd = wr = wv = wg = wc = 0.0
    for _ in range(samples):
        X = rng.normal(size=(3, 3))
        M = 0.5 * (X + X.T)
        M -= np.trace(M) / 3.0 * np.eye(3)
        B = swirl_tensor(M)
        ws = max(ws, float(np.linalg.norm(symmetrize_rank3(B))))
        wd = max(wd, float(np.linalg.norm(divergence_gradient(B))))
        Mr = reconstruct_M_from_kernel(B)
        wr = max(wr, float(np.linalg.norm(Mr - M)))
        z = rng.normal(size=3); q = rng.normal(size=3)
        V1 = quadratic_velocity(B, z); V2 = swirl_velocity(M, z)
        wv = max(wv, float(np.linalg.norm(V1 - V2)))
        wg = max(wg, abs(gaussian_envelope_advection(M, z)))
        C = carrier_chirp_matrix(M, q)
        direct = float(np.dot(q, V2))
        quad = float(z @ C @ z)
        wc = max(wc, abs(direct - quad))
    return SwirlStress(samples, ws, wd, wr, wv, wg, wc)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-quadratic-swirl-kernel"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    payload = {
        "theorem": {
            "kernel": "Sym_abc B=0, B_abc=B_acb, B_aac=0",
            "representation": "B_abc=eps_abd M_dc+eps_acd M_db with M symmetric tracefree",
            "velocity": "V(z)=z cross (M z)",
            "reconstruction": "M_dc=(1/3)eps_abd B_abc",
            "envelope": "V.grad exp(-|z|^2/4)=0",
            "carrier": "q.V is a quadratic chirp tangent",
        },
        "stress": asdict(out),
    }
    (args.outdir / "quadratic_swirl_kernel.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Quadratic incompressible swirl kernel

For a quadratic normalized velocity `V_a=(1/2)B_abc z_b z_c`, the divergence-free
third-Hermite kernel has the exact representation

`B_abc=eps_abd M_dc+eps_acd M_db`, `M=M^T`, `tr M=0`,

hence `V(z)=z cross (M z)`.  It is tangent to Gaussian level spheres:
`V.z=0`; a carrier sees only the quadratic chirp `q.V`.

- random checks: `{out.samples}`
- worst full-symmetry residual: `{out.worst_full_symmetry_residual:.3e}`
- worst divergence residual: `{out.worst_divergence_residual:.3e}`
- worst M reconstruction residual: `{out.worst_reconstruction_residual:.3e}`
- worst velocity representation residual: `{out.worst_velocity_residual:.3e}`
- worst radial Gaussian-advection residual: `{out.worst_radial_advection_residual:.3e}`
- worst chirp reconstruction residual: `{out.worst_chirp_residual:.3e}`

This five-dimensional mode is a real dynamical symmetry of the scalar Gaussian
envelope, not missing coercivity.  It must be routed to vector/polarization
transport rather than charged as scalar non-affine forcing.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
