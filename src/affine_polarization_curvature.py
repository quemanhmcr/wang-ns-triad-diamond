from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.extremal_helicity_symplectic import transfer_relevant_strain_observability
from src.quadratic_swirl_kernel import swirl_tensor


def physical_strain_gradient(H: np.ndarray, L: np.ndarray) -> np.ndarray:
    """C_ijc=sym_ij(H_ijk L_kc): strain variation per normalized grain coordinate."""
    H=np.asarray(H,float);L=np.asarray(L,float)
    G=np.einsum("ijk,kc->ijc",H,L)
    return 0.5*(G+G.swapaxes(0,1))


def divergence_gradient(H: np.ndarray) -> np.ndarray:
    return np.einsum("iik->k",np.asarray(H,float))


def rms_transfer_relevant_curvature(H: np.ndarray, L: np.ndarray, rstar: float=0.610904101586766) -> tuple[float,float]:
    """E_z Q_rel(S(z)) and E_z ||S(z)||^2 for z~N(0,I)."""
    C=physical_strain_gradient(H,L)
    qsum=0.0;nsum=0.0
    for c in range(3):
        Q,N=transfer_relevant_strain_observability(C[:,:,c],rstar)
        qsum+=Q;nsum+=N
    return qsum,nsum


def hessian_from_grain_tensor(B: np.ndarray, L: np.ndarray) -> np.ndarray:
    """Invert B=L^-1 H[L,L]."""
    B=np.asarray(B,float);L=np.asarray(L,float);Li=np.linalg.inv(L)
    return np.einsum("ia,abc,bj,ck->ijk",L,B,Li,Li)


@dataclass(frozen=True)
class CurvatureStress:
    samples:int
    worst_rms_observability_ratio:float
    worst_gaussian_expectation_residual:float
    worst_swirl_scalar_kernel:float
    minimum_swirl_polarization_signal:float


def stress(samples:int=50_000,seed:int=20260807)->CurvatureStress:
    rng=np.random.default_rng(seed)
    worst=float("inf");we=0.0;wsk=0.0;minsig=float("inf")
    monte=min(500,samples)
    for n in range(samples):
        # Physical divergence-free Hessian: symmetric in j,k and H_iik=0.
        H=rng.normal(size=(3,3,3));H=.5*(H+H.swapaxes(1,2))
        # Enforce differentiated incompressibility by solving H_22k from the trace.
        for k in range(3): H[2,2,k]=H[2,k,2]=-(H[0,0,k]+H[1,1,k])
        Q,_=np.linalg.qr(rng.normal(size=(3,3)));axes=np.exp(rng.uniform(-2,2,size=3));L=Q@np.diag(axes)@Q.T
        qv,nv=rms_transfer_relevant_curvature(H,L)
        if nv>1e-14:
            ratio=qv/nv;worst=min(worst,ratio)
            if ratio<.5-2e-10: raise AssertionError("integrated transfer-relevant tomography failed")
        if n<monte:
            C=physical_strain_gradient(H,L)
            zz=rng.normal(size=(4000,3))
            vals=[]; norms=[]
            for z in zz:
                S=np.einsum("ijc,c->ij",C,z)
                qx,nx=transfer_relevant_strain_observability(S)
                vals.append(qx);norms.append(nx)
            we=max(we,abs(float(np.mean(vals))-qv)/max(1.,qv),abs(float(np.mean(norms))-nv)/max(1.,nv))

        # Explicit scalar-forcing kernel transformed through an affine grain.
        X=rng.normal(size=(3,3));M=.5*(X+X.T);M-=np.trace(M)/3*np.eye(3)
        B=swirl_tensor(M); wsk=max(wsk,float(np.linalg.norm(sum(np.transpose(B,p) for p in __import__('itertools').permutations(range(3)))/6)))
        Hs=hessian_from_grain_tensor(B,L)
        qs,ns=rms_transfer_relevant_curvature(Hs,L)
        if np.linalg.norm(B)>1e-12:
            minsig=min(minsig,qs/(np.linalg.norm(B)**2))
            # Signal may deteriorate with aspect but cannot be negative.
            if qs<-1e-12: raise AssertionError("swirl polarization signal negative")
    return CurvatureStress(samples,worst,we,wsk,minsig)


def main()->None:
    ap=argparse.ArgumentParser();ap.add_argument("--samples",type=int,default=50_000);ap.add_argument("--outdir",type=Path,default=Path("results-affine-polarization-curvature"));args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    payload={"theorem":{
        "strain_gradient":"C_ijc=sym_ij(H_ijk L_kc)",
        "gaussian_rms":"E||S(z)||_F^2=||C||_F^2",
        "transfer_rms":"E Q_rel(S(z))=sum_c Q_rel(C_c)>=1/2||C||_F^2",
        "interpretation":"spatial strain/helical-generator curvature is the vector channel complementary to scalar third-Hermite forcing",
    },"stress":asdict(out)}
    (args.outdir/"affine_polarization_curvature.json").write_text(json.dumps(payload,indent=2))
    md=f"""# Affine-grain polarization curvature

For `A(x)=A0+H[Lz]`, let `C_ijc=sym_ij(H_ijk L_kc)`.  Differentiated
incompressibility makes each matrix `C_c` trace free.  The certified extremal
relative-polarization tomography theorem applies to every `C_c`, and Gaussian
orthogonality gives

`E_z Q_rel(S(z)) = sum_c Q_rel(C_c) >= (1/2)||C||_F^2`.

- random affine Hessian checks: `{out.samples}`
- worst observed RMS observability ratio: `{out.worst_rms_observability_ratio:.9f}`
- worst Monte-Carlo expectation residual: `{out.worst_gaussian_expectation_residual:.3e}`
- worst scalar third-Hermite symmetry residual for tested swirl kernels: `{out.worst_swirl_scalar_kernel:.3e}`
- minimum sampled swirl polarization signal / `||B||^2`: `{out.minimum_swirl_polarization_signal:.3e}`

Thus quadratic curvature has two distinct packet channels: full-symmetric
curvature creates third-Hermite envelope forcing, while physical symmetric
gradient variation creates transfer-distinguishable shape/polarization forcing.
The five-dimensional swirl kernel belongs to the second channel rather than the
first.  No claim is made here of an aspect-independent lower bound comparing the
second channel directly to the affine-normalized `||B||` norm.
"""
    (args.outdir/"summary.md").write_text(md);print(md)

if __name__=="__main__":main()
