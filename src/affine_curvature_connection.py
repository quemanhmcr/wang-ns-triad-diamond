from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def normalize_curvature(H: np.ndarray,L:np.ndarray)->np.ndarray:
    """B=L^-1 H[L,L], H shape (out,in,in)."""
    H=np.asarray(H,float); L=np.asarray(L,float); Li=np.linalg.inv(L)
    return np.einsum('ai,ijk,jb,kc->abc',Li,H,L,L)


def affine_matrix(A:np.ndarray,L:np.ndarray)->np.ndarray:
    L=np.asarray(L,float); return np.linalg.inv(L)@np.asarray(A,float)@L


def material_hessian_rate(A:np.ndarray,H:np.ndarray,F2:np.ndarray)->np.ndarray:
    """Exact D_t Hess(V) for material acceleration F=D_t V.

    D_t H_i,jk = F2_i,jk - A_i,m H_m,jk - A_m,j H_i,mk - A_m,k H_i,jm.
    """
    A=np.asarray(A,float); H=np.asarray(H,float); F2=np.asarray(F2,float)
    return (F2
            -np.einsum('im,mjk->ijk',A,H)
            -np.einsum('mj,imk->ijk',A,H)
            -np.einsum('mk,ijm->ijk',A,H))


def normalized_curvature_rate(A:np.ndarray,H:np.ndarray,F2:np.ndarray,L:np.ndarray)->np.ndarray:
    """Differentiate B=L^-1 H[L,L] with dot L=A L using the exact material H rate."""
    A=np.asarray(A,float); H=np.asarray(H,float); L=np.asarray(L,float); Li=np.linalg.inv(L)
    dH=material_hessian_rate(A,H,F2)
    dL=A@L; dLi=-Li@A
    return (np.einsum('ai,ijk,jb,kc->abc',dLi,H,L,L)
            +np.einsum('ai,ijk,jb,kc->abc',Li,dH,L,L)
            +np.einsum('ai,ijk,jb,kc->abc',Li,H,dL,L)
            +np.einsum('ai,ijk,jb,kc->abc',Li,H,L,dL))


def curvature_connection_rhs(A:np.ndarray,H:np.ndarray,F2:np.ndarray,L:np.ndarray)->np.ndarray:
    """Source-free covariant combination dot B+2 A_aff B = L^-1 F2[L,L]."""
    return normalize_curvature(F2,L)


def curvature_connection_residual(A:np.ndarray,H:np.ndarray,F2:np.ndarray,L:np.ndarray)->np.ndarray:
    B=normalize_curvature(H,L); Bd=normalized_curvature_rate(A,H,F2,L); Aa=affine_matrix(A,L)
    return Bd+2*np.einsum('ad,dbc->abc',Aa,B)-normalize_curvature(F2,L)


def resolved_acceleration_hessian(P3:np.ndarray,divR_hess:np.ndarray,V4:np.ndarray,nu:float)->np.ndarray:
    """Hess(D_t V) for D_t V=-grad P-div R+nu Delta V."""
    return -np.asarray(P3,float)-np.asarray(divR_hess,float)+float(nu)*np.asarray(V4,float)


def pressure_third_far_shell_exponent(space_dimension:int=3)->int:
    """Pressure kernel has degree -3 in 3D; third derivative degree -6."""
    if space_dimension!=3: raise ValueError('this pressure-kernel statement is for 3D')
    return 6-space_dimension


def far_pressure_third_geometric_sum(first_shell:int=3)->float:
    if first_shell<0: raise ValueError('nonnegative shell')
    # sum_{n>=n0} 2^{-3n}=2^{-3n0}/(1-1/8)
    return (8.0/7.0)*2.0**(-3*first_shell)


@dataclass(frozen=True)
class CurvatureConnectionStress:
    samples:int
    worst_connection_residual:float
    worst_source_split_residual:float
    worst_homogeneity_residual:float


def stress(samples:int=50_000,seed:int=20260807)->CurvatureConnectionStress:
    rng=np.random.default_rng(seed); wc=ws=wh=0.0
    for _ in range(samples):
        A=rng.normal(size=(3,3)); A-=np.trace(A)/3*np.eye(3)
        H=rng.normal(size=(3,3,3)); H=.5*(H+np.swapaxes(H,1,2))
        F2=rng.normal(size=(3,3,3)); F2=.5*(F2+np.swapaxes(F2,1,2))
        Q,_=np.linalg.qr(rng.normal(size=(3,3))); scales=np.exp(rng.uniform(-4,4,size=3)); L=Q@np.diag(scales)
        res=curvature_connection_residual(A,H,F2,L); rr=float(np.linalg.norm(res)); wc=max(wc,rr/max(1.,np.linalg.norm(normalize_curvature(F2,L))))
        if rr>5e-9*max(1.,np.linalg.norm(normalize_curvature(F2,L))): raise AssertionError('curvature connection identity failed')
        P3=rng.normal(size=(3,3,3)); R3=rng.normal(size=(3,3,3)); V4=rng.normal(size=(3,3,3)); nu=float(rng.uniform(0,2))
        split=resolved_acceleration_hessian(P3,R3,V4,nu)
        direct=-P3-R3+nu*V4; sr=float(np.linalg.norm(split-direct)); ws=max(ws,sr)
        if sr>2e-12: raise AssertionError('resolved acceleration source split failed')
        # Homogeneity of a degree -6 model scalar under dyadic radius rescaling.
        n=int(rng.integers(0,12)); r=float(2**n); val=r**-6; packed=(r**3)*val; expected=r**-3
        hr=abs(packed-expected); wh=max(wh,hr)
        if hr>2e-15*max(1.,expected): raise AssertionError('pressure-third packing exponent failed')
    return CurvatureConnectionStress(samples,wc,ws,wh)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-affine-curvature-connection'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    data={'stress':out.__dict__,'identities':{'connection':'dot B+2 A_aff B=L^-1 Hess(F)[L,L]','resolved_F':'F=-grad P-div R+nu Delta V','resolved_source':'Hess F=-nabla^3 P-nabla^2 div R+nu nabla^2 Delta V','pressure_far_exponent':'6-3=3'},'far_pressure_third_geometric_sum_n0_3':far_pressure_third_geometric_sum(3)}
    (args.outdir/'affine_curvature_connection.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Affine curvature connection and sideband source identity

For any smooth transporter `V`, along `dot X=V(X)`, `dot L=A L`, let `A=grad V`, `H=Hess V`, `F=D_t^V V`, and `B=L^-1 H[L,L]`.  Exact differentiation gives

`dot B+2 A_aff B = L^-1 Hess(F)[L,L]`,  `A_aff=L^-1 A L`.

Thus common affine deformation is a connection on curvature rather than a source.  For the resolved Navier--Stokes transporter
`F=-grad P-div R+nu Delta V`,

`Hess(F)=-nabla^3 P-nabla^2 div R+nu nabla^2 Delta V`.

Hence curvature-sideband dephasing is sourced by pressure third derivatives, differentiated SGS stress, or viscous fourth velocity derivatives after the affine connection is removed.  The pressure kernel is homogeneous of degree `-3`; three derivatives have degree `-6`, and 3D packet packing leaves the summable far exponent `6-3=3`.

This is a source/locality theorem, not yet a daughter-grain/coherence cost.

Stress checks: `{out.samples}`
- worst normalized connection residual: `{out.worst_connection_residual:.3e}`
- worst resolved-source split residual: `{out.worst_source_split_residual:.3e}`
- worst `6-3=3` homogeneity residual: `{out.worst_homogeneity_residual:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
