from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def split_pair_rescue(weights:np.ndarray,edges:np.ndarray,ancestry:np.ndarray)->dict[str,float]:
    weights=np.asarray(weights,float); edges=np.asarray(edges,int); ancestry=np.asarray(ancestry,int)
    if len(weights)!=len(edges) or np.any(weights<0): raise ValueError('invalid weighted graph')
    total=float(np.sum(weights))
    if total<=0: raise ValueError('positive pair-rescue mass required')
    same=0.0
    for a,(u,v) in zip(weights,edges):
        if u==v: raise ValueError('pair edge needs distinct odd atoms')
        if ancestry[u]==ancestry[v]: same+=float(a)
    return {'total':total,'same':same,'cross':total-same,'same_fraction':same/total,'cross_fraction':(total-same)/total}


def same_ancestry_endpoint_law(weights:np.ndarray,edges:np.ndarray,ancestry:np.ndarray)->tuple[np.ndarray,np.ndarray,float]:
    """Endpoint law on same-ancestry rescue edges.

    w_v=d_v/(2 W_same).  The induced ancestry weights equal the fraction of
    same-edge rescue mass in each ancestry class exactly.
    """
    weights=np.asarray(weights,float); edges=np.asarray(edges,int); ancestry=np.asarray(ancestry,int)
    n=len(ancestry); d=np.zeros(n); W=0.0
    for a,(u,v) in zip(weights,edges):
        if ancestry[u]==ancestry[v]:
            d[u]+=a; d[v]+=a; W+=float(a)
    if W<=0: raise ValueError('no same-ancestry rescue mass')
    w=d/(2*W)
    labels=np.unique(ancestry)
    WA=np.array([np.sum(w[ancestry==A]) for A in labels])
    return w,WA,W


def collision_chain(w:np.ndarray,WA:np.ndarray)->dict[str,float]:
    w=np.asarray(w,float); WA=np.asarray(WA,float)
    if abs(np.sum(w)-1)>1e-10 or abs(np.sum(WA)-1)>1e-10: raise ValueError('probability laws required')
    Qat=float(np.sum(w*w)); Qanc=float(np.sum(WA*WA))
    if Qanc+1e-13<Qat: raise AssertionError('ancestry collision must dominate atomic collision')
    return {'Q_atomic':Qat,'Q_ancestry':Qanc,'H_atomic':-math.log(Qat),'H_ancestry':-math.log(Qanc),'distinct_same_ancestry_pair_mass':Qanc-Qat}


def pair_rescue_ancestry_route(weights:np.ndarray,edges:np.ndarray,ancestry:np.ndarray,cross_threshold:float=0.5,h:float=0.7,alpha:float=0.5)->dict[str,float|str]:
    """Exact finite-graph routing of odd-sideband pair rescue.

    Branches:
      cross Xi; dominant reused atom; ancestry Bellman entropy; same-ancestry pair/cycle mass.
    """
    if not (0<cross_threshold<1 and h>0 and 0<alpha<1): raise ValueError('invalid thresholds')
    split=split_pair_rescue(weights,edges,ancestry)
    if split['cross_fraction']>=cross_threshold:
        return {'branch':'cross_Xi',**split,'cross_mass_lower':cross_threshold*split['total']}
    w,WA,Wsame=same_ancestry_endpoint_law(weights,edges,ancestry)
    c=collision_chain(w,WA); H=c['H_atomic']; Ha=c['H_ancestry']
    if H<h:
        vmax=float(np.max(w))
        if vmax+1e-13<math.exp(-h): raise AssertionError('dominant endpoint theorem failed')
        branch='dominant_reused_daughter'
        extra={'dominant_endpoint_weight':vmax,'dominant_lower':math.exp(-h)}
    elif Ha>=alpha*h:
        branch='ancestry_Bellman_entropy'
        extra={'ancestry_entropy_lower':alpha*h}
    else:
        lower=math.exp(-alpha*h)-math.exp(-h)
        if c['distinct_same_ancestry_pair_mass']+1e-12<lower: raise AssertionError('hidden pair-cycle mass theorem failed')
        branch='same_ancestry_pair_cycle'
        extra={'pair_cycle_mass_lower':lower}
    return {'branch':branch,**split,**c,**extra}


def ancestry_weights_equal_edge_class_mass(weights:np.ndarray,edges:np.ndarray,ancestry:np.ndarray)->float:
    w,WA,W=same_ancestry_endpoint_law(weights,edges,ancestry)
    labels=np.unique(ancestry); direct=[]
    for A in labels:
        mass=0.0
        for a,(u,v) in zip(weights,edges):
            if ancestry[u]==A and ancestry[v]==A: mass+=float(a)
        direct.append(mass/W)
    return float(np.linalg.norm(WA-np.asarray(direct)))


@dataclass(frozen=True)
class PairAncestryStress:
    samples:int
    worst_class_mass_identity_residual:float
    worst_chain_identity_residual:float
    minimum_branch_margin:float
    branch_counts:dict[str,int]


def random_graph(rng:np.random.Generator,n:int)->tuple[np.ndarray,np.ndarray,np.ndarray]:
    ancestry=rng.integers(0,max(2,n//4),size=n)
    m=int(rng.integers(n,4*n+1)); es=set()
    while len(es)<m:
        u,v=rng.integers(0,n,size=2)
        if u!=v: es.add(tuple(sorted((int(u),int(v)))))
        if len(es)>=n*(n-1)//2: break
    edges=np.asarray(list(es),int); weights=np.exp(rng.uniform(-4,1,size=len(edges)))
    return weights,edges,ancestry


def stress(samples:int=50_000,seed:int=20260807)->PairAncestryStress:
    rng=np.random.default_rng(seed); wi=wc=0.0; mm=float('inf'); counts={}
    for _ in range(samples):
        n=int(rng.integers(4,18)); weights,edges,anc=random_graph(rng,n)
        split=split_pair_rescue(weights,edges,anc)
        if split['same']<=1e-14:
            # force one same edge for identity testing if possible
            found=False
            for u in range(n):
                vs=np.where((anc==anc[u]) & (np.arange(n)!=u))[0]
                if len(vs):
                    edges=np.vstack([edges,[u,int(vs[0])]]); weights=np.r_[weights,1.0]; found=True; break
            if not found: continue
        wi=max(wi,ancestry_weights_equal_edge_class_mass(weights,edges,anc))
        w,WA,_=same_ancestry_endpoint_law(weights,edges,anc); c=collision_chain(w,WA)
        # direct chain identity Qanc-Qat equals distinct same-class endpoint pair mass
        direct=0.0
        for A in np.unique(anc):
            ix=np.where(anc==A)[0]
            for i in ix:
                for j in ix:
                    if i!=j: direct+=w[i]*w[j]
        cr=abs(direct-c['distinct_same_ancestry_pair_mass']); wc=max(wc,cr)
        if wi>2e-12 or cr>3e-12: raise AssertionError('pair ancestry exact identity failed')
        out=pair_rescue_ancestry_route(weights,edges,anc); counts[out['branch']]=counts.get(out['branch'],0)+1
        if out['branch']=='cross_Xi': margin=out['cross_fraction']-.5
        elif out['branch']=='dominant_reused_daughter': margin=out['dominant_endpoint_weight']-out['dominant_lower']
        elif out['branch']=='ancestry_Bellman_entropy': margin=out['H_ancestry']-out['ancestry_entropy_lower']
        else: margin=out['distinct_same_ancestry_pair_mass']-out['pair_cycle_mass_lower']
        mm=min(mm,float(margin))
        if margin<-2e-12: raise AssertionError('pair ancestry branch margin failed')
    return PairAncestryStress(samples,wi,wc,mm,counts)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-sideband-pair-ancestry'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    data={'stress':out.__dict__,'theorems':{
        'parity_graph':'nonzero odd-sideband rescue edges contain exactly two odd vertices and one even/base role',
        'same_endpoint_law':'w_v=d_v/(2W_same)',
        'ancestry_edge_mass':'W_A equals same-edge pair-rescue mass fraction in ancestry A',
        'chain':'Q_anc-Q_at=sum_A sum_{u!=v in A} w_u w_v',
        'route':'cross Xi / dominant reused daughter / ancestry Bellman entropy / same-ancestry pair-cycle mass',
    }}
    (args.outdir/'sideband_pair_ancestry.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Odd-sideband pair rescue to spacetime ancestry

Parity makes every nonzero rescue interaction containing odd Hermite modes use exactly **two odd daughter endpoints** and one even/base role.  Therefore pair rescue projects canonically to an ordinary weighted graph on odd daughter atoms.

Split rescue mass into `W_cross` (odd endpoints in different ancestry components) and `W_same`.  A cross fraction at least `1/2` is already an existing `Xi`/cross-component branch.

On same-ancestry edges define

`w_v=d_v/(2 W_same)`.

For each ancestry class `A`, `W_A=sum_(v in A) w_v` equals **exactly** the fraction of same-edge rescue mass lying in `A`.  Hence

`Q_anc-Q_at = sum_A sum_(u!=v in A) w_u w_v`.

With thresholds `h>0`, `alpha in (0,1)`, the same-edge graph has the exact trichotomy:

1. `H_at<h` -> some odd daughter carries endpoint weight `> exp(-h)`: dominant reused daughter;
2. `H_at>=h` and `H_anc>=alpha h` -> ancestry/component Bellman collision entropy;
3. `H_at>=h` and `H_anc<alpha h` -> distinct same-ancestry pair mass at least `exp(-alpha h)-exp(-h)`, hence repeated attachments which become cycle-rank gain after contraction by the existing atomic-to-ancestry theorem.

Together with the cross branch, pair-sideband rescue is routed into the existing currencies `Xi`, dominant reuse, Bellman entropy, or ancestry cycles.  No spatial disjointness of Hermite modes is assumed.

Stress: `{out.samples}`
- branch counts: `{out.branch_counts}`
- worst ancestry-class mass identity residual: `{out.worst_class_mass_identity_residual:.3e}`
- worst collision-chain identity residual: `{out.worst_chain_identity_residual:.3e}`
- minimum branch margin: `{out.minimum_branch_margin:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
