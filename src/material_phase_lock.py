from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


def kelvin_rhs(k: np.ndarray, A: np.ndarray, residual: np.ndarray | None = None) -> np.ndarray:
    out = -np.asarray(A, float).T @ np.asarray(k, float)
    if residual is not None:
        out = out + np.asarray(residual, float)
    return out


def hessian_phase_rhs(
    K: np.ndarray,
    k: np.ndarray,
    A: np.ndarray,
    H: np.ndarray,
    residual: np.ndarray | None = None,
) -> np.ndarray:
    """Material Hessian equation for D_t phi=0 plus an optional source.

    H[i,j,k]=partial_j partial_k U_i and H[kcarrier]_{jk}=k_i H[i,j,k].
    """
    A = np.asarray(A, float); K = np.asarray(K, float); k = np.asarray(k, float)
    Hk = np.einsum("i,ijk->jk", k, np.asarray(H, float))
    out = -A.T @ K - K @ A - Hk
    if residual is not None:
        out = out + np.asarray(residual, float)
    return 0.5 * (out + out.T)


def phase_lock_source(rho1: float, rho2: float, rho3: float) -> float:
    return float(rho1 + rho2 - rho3)


def gradient_lock_source(r1: np.ndarray, r2: np.ndarray, r3: np.ndarray) -> np.ndarray:
    return np.asarray(r1, float) + np.asarray(r2, float) - np.asarray(r3, float)


def hessian_lock_source(G1: np.ndarray, G2: np.ndarray, G3: np.ndarray) -> np.ndarray:
    return np.asarray(G1, float) + np.asarray(G2, float) - np.asarray(G3, float)


def rk4_step(k: np.ndarray, K: np.ndarray, t: float, dt: float, A_fun, H_fun):
    def rhs(kk, KK, tt):
        A=A_fun(tt); H=H_fun(tt)
        return kelvin_rhs(kk,A), hessian_phase_rhs(KK,kk,A,H)
    k1,K1=rhs(k,K,t)
    k2,K2=rhs(k+.5*dt*k1,K+.5*dt*K1,t+.5*dt)
    k3,K3=rhs(k+.5*dt*k2,K+.5*dt*K2,t+.5*dt)
    k4,K4=rhs(k+dt*k3,K+dt*K3,t+dt)
    return k+dt*(k1+2*k2+2*k3+k4)/6, K+dt*(K1+2*K2+2*K3+K4)/6


@dataclass(frozen=True)
class PhaseLockStress:
    samples: int
    worst_gradient_lock_residual: float
    worst_hessian_lock_residual: float
    worst_common_hessian_source_residual: float
    minimum_differential_source_margin: float


def stress(samples: int = 50_000, seed: int = 20260807) -> PhaseLockStress:
    rng=np.random.default_rng(seed)
    wg=wh=wH=0.0
    margin=float("inf")
    # Algebraic random checks.
    for _ in range(samples):
        A=rng.normal(size=(3,3))
        H=rng.normal(size=(3,3,3)); H=.5*(H+H.swapaxes(1,2))
        k1=rng.normal(size=3); k2=rng.normal(size=3); k3=k1+k2
        K1=rng.normal(size=(3,3)); K1=.5*(K1+K1.T)
        K2=rng.normal(size=(3,3)); K2=.5*(K2+K2.T); K3=K1+K2
        g=kelvin_rhs(k1,A)+kelvin_rhs(k2,A)-kelvin_rhs(k3,A)
        h=hessian_phase_rhs(K1,k1,A,H)+hessian_phase_rhs(K2,k2,A,H)-hessian_phase_rhs(K3,k3,A,H)
        wg=max(wg,float(np.linalg.norm(g))); wh=max(wh,float(np.linalg.norm(h)))
        common=np.einsum("i,ijk->jk",k1+k2-k3,H)
        wH=max(wH,float(np.linalg.norm(common)))

        r1=rng.normal(size=3);r2=rng.normal(size=3);r3=rng.normal(size=3)
        G1=rng.normal(size=(3,3));G1=.5*(G1+G1.T)
        G2=rng.normal(size=(3,3));G2=.5*(G2+G2.T)
        G3=rng.normal(size=(3,3));G3=.5*(G3+G3.T)
        lhs=np.linalg.norm(gradient_lock_source(r1,r2,r3))
        rhs=np.linalg.norm(r1)+np.linalg.norm(r2)+np.linalg.norm(r3)
        margin=min(margin,rhs-lhs)
        lhs2=np.linalg.norm(hessian_lock_source(G1,G2,G3))
        rhs2=np.linalg.norm(G1)+np.linalg.norm(G2)+np.linalg.norm(G3)
        margin=min(margin,rhs2-lhs2)

    # Time-ordered nonlinear-gradient/Hessian history: the three ODEs are solved
    # separately, yet exact initial resonance remains at numerical RK4 error.
    coeffA=rng.normal(size=(3,3,4)); coeffH=rng.normal(size=(3,3,3,4)); coeffH=.5*(coeffH+coeffH.swapaxes(1,2))
    def A_fun(t): return sum(coeffA[:,:,j]*t**j for j in range(4))*.08
    def H_fun(t): return sum(coeffH[:,:,:,j]*t**j for j in range(4))*.05
    k1=rng.normal(size=3); k2=rng.normal(size=3); k3=k1+k2
    K1=rng.normal(size=(3,3)); K1=.5*(K1+K1.T)
    K2=rng.normal(size=(3,3)); K2=.5*(K2+K2.T); K3=K1+K2
    dt=1/2000
    for n in range(2000):
        t=n*dt
        k1,K1=rk4_step(k1,K1,t,dt,A_fun,H_fun)
        k2,K2=rk4_step(k2,K2,t,dt,A_fun,H_fun)
        k3,K3=rk4_step(k3,K3,t,dt,A_fun,H_fun)
    wg=max(wg,float(np.linalg.norm(k1+k2-k3)))
    wh=max(wh,float(np.linalg.norm(K1+K2-K3)))
    return PhaseLockStress(samples,wg,wh,wH,margin)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--outdir",type=Path,default=Path("results-material-phase-lock"))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    payload={
        "theorem":{
            "material_phase":"D_t(phi1+phi2-phi3)=rho1+rho2-rho3",
            "gradient":"D_t(k1+k2-k3)=-A^T(k1+k2-k3)+(r1+r2-r3)",
            "hessian":"D_t Klock=-A^T Klock-Klock A-H_U[klock]+(G1+G2-G3)",
            "consequence":"common nonlinear advection, including common velocity Hessian chirp, is an exact triad phase gauge",
        },
        "stress":asdict(out),
    }
    (args.outdir/"material_phase_lock.json").write_text(json.dumps(payload,indent=2))
    md=f"""# Material triad phase-lock gauge

For three phases transported by the same resolved velocity, the signed triad
phase `Phi=phi1+phi2-phi3` is itself a materially transported scalar.  Exact
carrier resonance and quadratic chirp lock therefore persist under arbitrary
common non-affine advection.  Only differential packet/resolved-flow sources can
dephase the transfer.

- random algebraic checks: `{out.samples}`
- worst gradient-lock residual: `{out.worst_gradient_lock_residual:.3e}`
- worst Hessian/chirp-lock residual: `{out.worst_hessian_lock_residual:.3e}`
- worst common Hessian-source residual at exact resonance: `{out.worst_common_hessian_source_residual:.3e}`
- minimum differential-source triangle margin: `{out.minimum_differential_source_margin:.3e}`

This explains why the quadratic `q.B` term in the affine Gaussian forcing module
belongs to a common phase gauge rather than to the transfer-facing residual.
"""
    (args.outdir/"summary.md").write_text(md); print(md)

if __name__=="__main__": main()
