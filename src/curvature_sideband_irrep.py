from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.hermite_helicity_ledger import sym3


def divergence_trace(B:np.ndarray)->np.ndarray:
    return np.einsum('aac->c',np.asarray(B,float))


def symmetric_trace(T:np.ndarray)->np.ndarray:
    T=np.asarray(T,float)
    return np.einsum('aac->c',T)


def envelope_reconstruction(B:np.ndarray)->np.ndarray:
    """Orthogonal divergence-free reconstruction from T=Sym B.

    For B[a,b,c]=B[a,c,b] and sum_a B[a,a,c]=0, let T=Sym B and t=T[a,a,c].
    Then
      B_E=T-1/2(delta_ab t_c+delta_ac t_b)+delta_bc t_a
    has the same full symmetrization and is divergence free.
    """
    T=sym3(np.asarray(B,float)); t=symmetric_trace(T); I=np.eye(3)
    return (T
            -0.5*np.einsum('ab,c->abc',I,t)
            -0.5*np.einsum('ac,b->abc',I,t)
            +np.einsum('bc,a->abc',I,t))


def hook_component(B:np.ndarray)->np.ndarray:
    return np.asarray(B,float)-envelope_reconstruction(B)


def levi_civita()->np.ndarray:
    e=np.zeros((3,3,3))
    e[0,1,2]=e[1,2,0]=e[2,0,1]=1
    e[0,2,1]=e[2,1,0]=e[1,0,2]=-1
    return e


def hook_M(K:np.ndarray)->np.ndarray:
    """Inverse of K_abc=eps_abd M_dc+eps_acd M_db on the divergence-free hook sector."""
    e=levi_civita(); K=np.asarray(K,float)
    return (1.0/3.0)*np.einsum('abe,abc->ec',e,K)


def hook_from_M(M:np.ndarray)->np.ndarray:
    e=levi_civita(); M=np.asarray(M,float)
    return np.einsum('abd,dc->abc',e,M)+np.einsum('acd,db->abc',e,M)


def hook_strain_sideband(K:np.ndarray)->np.ndarray:
    """Normalized grain-coordinate H1 strain coefficient C_ad,c=sym_ad K_ad,c."""
    K=np.asarray(K,float)
    return 0.5*(K+np.swapaxes(K,0,1))


def irrep_norm_ledger(B:np.ndarray)->dict[str,float]:
    B=np.asarray(B,float); T=sym3(B); t=symmetric_trace(T); E=envelope_reconstruction(B); H=B-E; M=hook_M(H); C=hook_strain_sideband(H)
    return {
        'B2':float(np.sum(B*B)),
        'T2':float(np.sum(T*T)),
        'trace2':float(np.dot(t,t)),
        'envelope2':float(np.sum(E*E)),
        'hook2':float(np.sum(H*H)),
        'M2':float(np.sum(M*M)),
        'hook_H1_2':float(np.sum(C*C)),
        'orthogonality':float(np.einsum('abc,abc',E,H)),
    }


def clean_combined_observability(B:np.ndarray)->float:
    """Return ||Sym B||^2+||C_hook||^2-(1/6)||B||^2."""
    B=np.asarray(B,float); T=sym3(B); H=hook_component(B); C=hook_strain_sideband(H)
    return float(np.sum(T*T)+np.sum(C*C)-np.sum(B*B)/6.0)


_DIVFREE_COORDS=[(a,b,c) for a in range(3) for b in range(3) for c in range(b,3)]

def _full_curvature_from_coords(x:np.ndarray)->np.ndarray:
    B=np.zeros((3,3,3))
    for val,(a,b,c) in zip(x,_DIVFREE_COORDS):
        B[a,b,c]=B[a,c,b]=val
    return B

def _build_divfree_nullspace()->np.ndarray:
    n=len(_DIVFREE_COORDS); C=np.zeros((3,n)); eye=np.eye(n)
    for j in range(n): C[:,j]=divergence_trace(_full_curvature_from_coords(eye[j]))
    _,_,vh=np.linalg.svd(C)
    return vh[3:].T

_DIVFREE_NULLSPACE=_build_divfree_nullspace()

def random_divfree_curvature(rng:np.random.Generator)->np.ndarray:
    """Generate jk-symmetric B with divergence trace zero by cached nullspace coordinates."""
    x=_DIVFREE_NULLSPACE@rng.normal(size=_DIVFREE_NULLSPACE.shape[1])
    return _full_curvature_from_coords(x)


@dataclass(frozen=True)
class CurvatureIrrepStress:
    samples:int
    worst_envelope_sym_residual:float
    worst_hook_sym_residual:float
    worst_hook_divergence_residual:float
    worst_orthogonality_residual:float
    worst_envelope_norm_identity_residual:float
    worst_hook_M_reconstruction_residual:float
    worst_hook_norm_identity_residual:float
    worst_H1_norm_identity_residual:float
    minimum_combined_observability_margin:float
    maximum_trace_ratio:float


def stress(samples:int=50_000,seed:int=20260807)->CurvatureIrrepStress:
    rng=np.random.default_rng(seed)
    wes=whs=whd=wo=wen=wm=whn=wc=0.0; minobs=float('inf'); maxtr=0.0
    for _ in range(samples):
        B=random_divfree_curvature(rng); T=sym3(B); t=symmetric_trace(T); E=envelope_reconstruction(B); H=B-E
        wes=max(wes,float(np.linalg.norm(sym3(E)-T)))
        whs=max(whs,float(np.linalg.norm(sym3(H))))
        whd=max(whd,float(np.linalg.norm(divergence_trace(H))))
        wo=max(wo,abs(float(np.einsum('abc,abc',E,H))))
        en=abs(float(np.sum(E*E)-(np.sum(T*T)+3*np.dot(t,t)))); wen=max(wen,en)
        M=hook_M(H); rec=hook_from_M(M); wm=max(wm,float(np.linalg.norm(rec-H)))
        whn=max(whn,abs(float(np.sum(H*H)-6*np.sum(M*M))))
        C=hook_strain_sideband(H); wc=max(wc,abs(float(np.sum(C*C)-0.25*np.sum(H*H))))
        obs=clean_combined_observability(B); minobs=min(minobs,obs)
        if np.sum(T*T)>1e-14: maxtr=max(maxtr,float(np.dot(t,t)/np.sum(T*T)))
        scale=max(1.,float(np.linalg.norm(B)))
        if max(wes,whs,whd,wm)>2e-11*scale: raise AssertionError('irrep tensor identity failed')
        if abs(float(np.einsum('abc,abc',E,H)))>2e-11*scale*scale: raise AssertionError('irrep orthogonality failed')
        if en>2e-10*scale*scale or whn>2e-10*scale*scale or wc>2e-10*scale*scale: raise AssertionError('irrep norm identity failed')
        if obs<-2e-11*scale*scale: raise AssertionError('combined H3/H1 observability failed')
        if np.sum(T*T)>1e-14 and np.dot(t,t)>(5/3)*np.sum(T*T)+2e-11*scale*scale: raise AssertionError('sharp symmetric trace bound failed')
    return CurvatureIrrepStress(samples,wes,whs,whd,wo,wen,wm,whn,wc,minobs,maxtr)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-curvature-sideband-irrep'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    data={'stress':out.__dict__,'identities':{
        'envelope':'B_E=T-1/2(delta_ab t_c+delta_ac t_b)+delta_bc t_a',
        'orthogonal_split':'B=B_E+B_H, <B_E,B_H>=0',
        'envelope_norm':'||B_E||^2=||T||^2+3||t||^2',
        'sharp_trace':'||t||^2<=(5/3)||T||^2',
        'hook':'B_H=eps_abd M_dc+eps_acd M_db, M symmetric tracefree',
        'hook_norm':'||B_H||^2=6||M||^2',
        'hook_H1':'||C_hook||^2=(1/4)||B_H||^2',
        'combined':'||Sym B||^2+||C_hook||^2 >= (1/6)||B||^2',
    }}
    (args.outdir/'curvature_sideband_irrep.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Divergence-free curvature irreducible sideband split

For `B[a,b,c]=B[a,c,b]` with differentiated incompressibility `B[a,a,c]=0`, let `T=Sym B`, `t_c=T[a,a,c]`.  The exact divergence-free tensor carrying the envelope symmetrization is

`B_E = T - 1/2(delta_ab t_c+delta_ac t_b) + delta_bc t_a`.

Then `B_H=B-B_E` is orthogonal to `B_E`, has `Sym B_H=0` and remains divergence free.  Hence it is exactly the five-dimensional quadratic-swirl sector

`B_H[a,b,c]=eps[a,b,d] M[d,c]+eps[a,c,d] M[d,b]`, `M=M^T`, `tr M=0`.

Exact norm identities:
- `||B_E||^2=||T||^2+3||t||^2`;
- sharp symmetric trace bound `||t||^2<=(5/3)||T||^2`, hence `||B_E||^2<=6||T||^2`;
- `||B_H||^2=6||M||^2`;
- normalized hook H1 strain coefficient `C_hook=1/2(B_H+swap_output_input(B_H))` obeys `||C_hook||^2=(1/4)||B_H||^2`.

Therefore the two intrinsic sideband channels obey the aspect-free curvature observability

`||Sym B||^2 + ||C_hook||^2 >= (1/6)||B||^2`.

This is the representation-theoretic statement `15=(7+3)_envelope + 5_swirl`: normalized non-affine curvature cannot disappear simultaneously from the H3 envelope and H1/swirl sectors.  The `C_hook` norm is an intrinsic grain-coordinate sideband norm; comparison with physical Euclidean helicity curvature still uses the existing polarization/ancestry bridge and is not claimed aspect-uniformly here.

Stress: `{out.samples}`
- worst envelope-sym residual: `{out.worst_envelope_sym_residual:.3e}`
- worst hook-sym residual: `{out.worst_hook_sym_residual:.3e}`
- worst hook-divergence residual: `{out.worst_hook_divergence_residual:.3e}`
- worst orthogonality residual: `{out.worst_orthogonality_residual:.3e}`
- worst envelope-norm identity residual: `{out.worst_envelope_norm_identity_residual:.3e}`
- worst hook-M reconstruction residual: `{out.worst_hook_M_reconstruction_residual:.3e}`
- worst hook-norm identity residual: `{out.worst_hook_norm_identity_residual:.3e}`
- worst H1-norm identity residual: `{out.worst_H1_norm_identity_residual:.3e}`
- minimum combined-observability margin: `{out.minimum_combined_observability_margin:.3e}`
- maximum sampled trace ratio: `{out.maximum_trace_ratio:.9f}` (sharp bound `5/3`)
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
