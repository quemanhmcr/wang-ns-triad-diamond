from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.transfer_profile_extraction import trilinear_replacement_loss


def base_child_parent_remainder_loss(epsilon: float) -> float:
    """Work-level bound for T(f,g,H)-T(F,G,H).

    Assume ||f||=||g||=1, ||f-F||,||g-G||,||h-H||<=eps and hence
    ||F||,||G||,||H||<=1+eps.  For a trilinear form of norm one,
      T(f,g,H)-T(F,G,H)
       =T(f-F,g,H)+T(F,g-G,H),
    giving eps(1+eps)+(1+eps)^2 eps.
    """
    if epsilon < 0: raise ValueError("epsilon must be nonnegative")
    e=epsilon
    return 2*e+3*e*e+e**3


def child_representation_loss(epsilon: float) -> float:
    """T(f,g,h)-T(f,g,H) <= eps for unit actual parents."""
    if epsilon < 0: raise ValueError("epsilon must be nonnegative")
    return epsilon


def split_replacement_identity(epsilon: float) -> float:
    return base_child_parent_remainder_loss(epsilon)+child_representation_loss(epsilon)


def rank_one_trilinear(a: np.ndarray,b: np.ndarray,c: np.ndarray,x: np.ndarray,y: np.ndarray,z: np.ndarray) -> complex:
    return np.vdot(a,x)*np.vdot(b,y)*np.vdot(c,z)


def unit(v: np.ndarray)->np.ndarray:
    n=np.linalg.norm(v)
    if n<1e-14: raise ValueError("zero")
    return v/n


@dataclass(frozen=True)
class CrossForcingStress:
    samples: int
    worst_base_source_ratio: float
    worst_child_mismatch_ratio: float
    worst_split_identity_residual: float


def stress(samples:int=50_000,seed:int=20260807)->CrossForcingStress:
    rng=np.random.default_rng(seed); wb=wc=wi=0.0
    dim=5
    for _ in range(samples):
        a=unit(rng.normal(size=dim)+1j*rng.normal(size=dim)); b=unit(rng.normal(size=dim)+1j*rng.normal(size=dim)); c=unit(rng.normal(size=dim)+1j*rng.normal(size=dim))
        f=unit(rng.normal(size=dim)+1j*rng.normal(size=dim)); g=unit(rng.normal(size=dim)+1j*rng.normal(size=dim)); h=unit(rng.normal(size=dim)+1j*rng.normal(size=dim))
        eps=float(rng.uniform(1e-7,.05))
        rf=unit(rng.normal(size=dim)+1j*rng.normal(size=dim))*float(rng.uniform(0,eps))
        rg=unit(rng.normal(size=dim)+1j*rng.normal(size=dim))*float(rng.uniform(0,eps))
        rh=unit(rng.normal(size=dim)+1j*rng.normal(size=dim))*float(rng.uniform(0,eps))
        F=f-rf; G=g-rg; H=h-rh
        main=rank_one_trilinear(a,b,c,f,g,H)-rank_one_trilinear(a,b,c,F,G,H)
        bound=base_child_parent_remainder_loss(eps)
        if abs(main)>bound+2e-12: raise AssertionError("base-child cross forcing bound failed")
        wb=max(wb,abs(main)/bound if bound else 0.0)
        child=rank_one_trilinear(a,b,c,f,g,h)-rank_one_trilinear(a,b,c,f,g,H)
        cb=child_representation_loss(eps)
        if abs(child)>cb+2e-12: raise AssertionError("child representation bound failed")
        wc=max(wc,abs(child)/cb if cb else 0.0)
        ident=split_replacement_identity(eps)-trilinear_replacement_loss(eps)
        wi=max(wi,abs(ident))
        if abs(ident)>2e-14: raise AssertionError("replacement split identity failed")
    return CrossForcingStress(samples,wb,wc,wi)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-packet-cross-forcing'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    e=.01
    data={'stress':out.__dict__,'at_one_percent':{'base_child_parent_remainder_loss':base_child_parent_remainder_loss(e),'child_representation_loss':child_representation_loss(e),'full_replacement_loss':trilinear_replacement_loss(e)}}
    (args.outdir/'packet_cross_forcing.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Base-spinor cross forcing is the existing trilinear cross-error

For normalized actual parents `f,g` and Gaussian roles `F,G,H` with distance at most `eps`, the degree-zero/base-child source error obeys

`|T(f,g,H)-T(F,G,H)| <= (2 eps+3 eps^2+eps^3) ||T||`.

The child representation mismatch obeys

`|T(f,g,h)-T(f,g,H)| <= eps ||T||`.

Their sum is exactly the one-shot replacement polynomial
`3 eps+3 eps^2+eps^3` already used in the profile ledger.
At `eps=1%` this splits as `{base_child_parent_remainder_loss(e):.6f} + {e:.6f} = {trilinear_replacement_loss(e):.6f}`.

Thus the **work-level degree-zero forcing** produced by parent remainders/cross components is already an omitted trilinear cross interaction and belongs to the existing `eta_j` / `Xi` ledger.  This does not claim an `L^2` bound for the entire nonlinear residual; orthogonal Hermite sidebands are handled separately.

Stress checks: `{out.samples}` rank-one norm-one complex trilinear forms
- worst base-source/bound ratio: `{out.worst_base_source_ratio:.9f}`
- worst child-mismatch/bound ratio: `{out.worst_child_mismatch_ratio:.9f}`
- worst replacement-split identity residual: `{out.worst_split_identity_residual:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
