from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict,dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.causal_binary_ancestry import root_count_energy_upper


def entropy(weights:Sequence[float])->float:
    w=np.asarray(weights,float)
    if np.any(w<0) or w.sum()<=0: raise ValueError('positive probability mass required')
    w=w/w.sum(); nz=w[w>0]
    return float(-np.dot(nz,np.log(nz)))


def parent_pushforward(child_weights:Sequence[float], parent_pairs:np.ndarray, n_parents:int|None=None)->np.ndarray:
    w=np.asarray(child_weights,float)
    if np.any(w<0) or w.sum()<=0: raise ValueError('positive child law required')
    w=w/w.sum(); pairs=np.asarray(parent_pairs,int)
    if pairs.shape!=(len(w),2) or np.any(pairs<0): raise ValueError('parent pairs must have shape (children,2)')
    if n_parents is None: n_parents=int(pairs.max())+1
    if n_parents<=int(pairs.max()): raise ValueError('n_parents too small')
    out=np.zeros(n_parents,float)
    for c,wc in enumerate(w):
        out[pairs[c,0]] += .5*wc
        out[pairs[c,1]] += .5*wc
    if abs(out.sum()-1)>2e-13: raise AssertionError('parent pushforward lost mass')
    return out


def layer_reuse_information(child_weights:Sequence[float], parent_pairs:np.ndarray, n_parents:int|None=None)->dict[str,float|np.ndarray]:
    w=np.asarray(child_weights,float); w=w/w.sum()
    p=parent_pushforward(w,parent_pairs,n_parents)
    hc=entropy(w); hp=entropy(p)
    reuse=hc+math.log(2.0)-hp
    if reuse<-2e-13: raise AssertionError('conditional reuse information became negative')
    return {'child_entropy':hc,'parent_entropy':hp,'reuse_information':max(0.0,reuse),'parent_weights':p}


def layered_reuse_information(parent_maps:Sequence[np.ndarray],terminal_weights:Sequence[float]|None=None)->dict[str,object]:
    """Maps are ordered root->terminal: map[j] has shape (n_{j+1},2) into level-j parents."""
    maps=[np.asarray(x,int) for x in parent_maps]
    if not maps: raise ValueError('at least one layer required')
    nL=maps[-1].shape[0]
    if terminal_weights is None:
        if nL!=1: raise ValueError('default terminal law requires one terminal packet')
        w=np.ones(1)
    else:
        w=np.asarray(terminal_weights,float); w=w/w.sum()
        if len(w)!=nL: raise ValueError('terminal weight size mismatch')
    infos=[None]*len(maps); laws=[None]*(len(maps)+1); laws[-1]=w.copy()
    for j in range(len(maps)-1,-1,-1):
        npar=int(maps[j].max())+1
        row=layer_reuse_information(w,maps[j],npar)
        infos[j]=float(row['reuse_information']); w=np.asarray(row['parent_weights']); laws[j]=w.copy()
    total=float(sum(infos)); exact=entropy(laws[-1])+len(maps)*math.log(2.0)-entropy(laws[0])
    if abs(total-exact)>5e-13*max(1.0,abs(exact)): raise AssertionError('weighted reuse information telescope failed')
    return {'reuse_information_by_layer':infos,'total_reuse_information':total,'root_entropy':entropy(laws[0]),'terminal_entropy':entropy(laws[-1]),'root_weights':laws[0],'laws':laws}


def physical_reuse_information_lower(depth:int,base_frequency:float,global_energy:float,critical_mass:float=.2,frame_budget:float=1.0)->float:
    if depth<0 or min(base_frequency,global_energy,critical_mass,frame_budget)<=0: raise ValueError('bad physical data')
    # H(root)<=log n0 and n0<=P E Nbase (25/24)^L / eta.
    return depth*math.log(48.0/25.0)-math.log(frame_budget*global_energy*base_frequency/critical_mass)


def rich_layer_threshold()->float:
    return math.log(4.0/3.0)


def first_information_rich_depth(base_frequency:float,global_energy:float,critical_mass:float=.2,frame_budget:float=1.0)->int:
    B=frame_budget*global_energy*base_frequency/critical_mass
    if B<=0: raise ValueError('positive budget required')
    if B<1: return 0
    # L log(48/25)-log B > L log(4/3) iff L log(36/25)>log B.
    x=math.log(B)/math.log(36.0/25.0)
    L=int(math.floor(x))+1
    while L*math.log(36.0/25.0)<=math.log(B): L+=1
    return L


def clean_rich_layer_consequence(depth:int,base_frequency:float,global_energy:float,critical_mass:float=.2,frame_budget:float=1.0)->dict[str,float|int|str]:
    lo=physical_reuse_information_lower(depth,base_frequency,global_energy,critical_mass,frame_budget)
    avg=lo/depth if depth>0 else 0.0
    th=rich_layer_threshold()
    if lo>depth*th:
        return {'branch':'transfer_weighted_reuse_information_rich','total_reuse_information_lower':lo,'one_layer_reuse_information_lower':avg,'clean_layer_threshold':th}
    return {'branch':'causal_depth_not_yet_forced','total_reuse_information_lower':lo,'average_lower':avg,'first_forced_depth':first_information_rich_depth(base_frequency,global_energy,critical_mass,frame_budget)}


def random_surjective_map(rng:np.random.Generator,nchild:int,nparent:int)->np.ndarray:
    if not (1<=nparent<=2*nchild): raise ValueError('bad layer sizes')
    slots=np.empty(2*nchild,dtype=int)
    slots[:nparent]=np.arange(nparent)
    if nparent<2*nchild: slots[nparent:]=rng.integers(0,nparent,size=2*nchild-nparent)
    rng.shuffle(slots)
    return slots.reshape(nchild,2)


def random_layered_maps(rng:np.random.Generator,depth:int)->list[np.ndarray]:
    counts=[0]*(depth+1); counts[-1]=1
    for j in range(depth-1,-1,-1): counts[j]=int(rng.integers(1,2*counts[j+1]+1))
    return [random_surjective_map(rng,counts[j+1],counts[j]) for j in range(depth)]

@dataclass(frozen=True)
class WeightedCausalStress:
    samples:int
    worst_telescope_residual:float
    minimum_conditional_information:float
    minimum_depth_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260808)->WeightedCausalStress:
    rng=np.random.default_rng(seed); wt=0.0; mi=float('inf'); md=float('inf'); branches={}
    for _ in range(samples):
        L=int(rng.integers(1,16)); maps=random_layered_maps(rng,L); out=layered_reuse_information(maps); exact=L*math.log(2)-float(out['root_entropy']); wt=max(wt,abs(float(out['total_reuse_information'])-exact)); mi=min(mi,min(out['reuse_information_by_layer']))
        if min(out['reuse_information_by_layer'])<-2e-13: raise AssertionError('negative conditional reuse information')
        N=float(rng.lognormal(0,1)); E=float(rng.lognormal(0,1)); c=clean_rich_layer_consequence(L,N,E); b=str(c['branch']); branches[b]=branches.get(b,0)+1
        if b=='transfer_weighted_reuse_information_rich':
            md=min(md,float(c['total_reuse_information_lower'])-L*rich_layer_threshold())
    return WeightedCausalStress(samples,wt,mi,md,branches)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-weighted-causal-reuse'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); out=stress(args.samples)
    data={'status':'EXACT_TRANSFER_WEIGHTED_CAUSAL_REUSE_INFORMATION','theorem':{
        'layer_identity':'R_j=H(child)+log2-H(parent)=H(child,role|parent)>=0',
        'telescope':'sum R_j=L log2-H(root) for one terminal atom',
        'physical_lower':'sum R_j >= L log(48/25)-log(P E_global N_base/eta)',
        'rich_layer':'L log(36/25)>log(P E_global N_base/eta) forces some R_j>log(4/3)',
        'role_baseline':'the log2 two-parent-role entropy is baseline causality and is not itself charged as Bellman cost',
    },'stress':asdict(out)}
    (args.outdir/'weighted_causal_reuse.json').write_text(json.dumps(data,indent=2))
    md=f'''# Transfer-weighted causal reuse information

Status: **EXACT_TRANSFER_WEIGHTED_CAUSAL_REUSE_INFORMATION**.

Let `w_(j+1)` be the normalized physical causal-transfer law on child packets of one synchronized layer.  Give each of its two parent-role slots weight `w/2`, and push this joint `(child,role)` law through the physical parent-label map to obtain `w_j`.

Define

`R_j = H(child,role | parent) = H(w_(j+1)) + log 2 - H(w_j) >=0`.

The `log 2` is the exact cost-free two-parent-role baseline.  `R_j`, not `log 2`, measures transfer-weighted causal merging/reuse.

For a one-terminal ancestry the identities telescope exactly:

`sum_j R_j = L log 2 - H(w_0)`.

Since `H(w_0)<=log n_0` and the coherent root-energy / signed-good scale theorem gives `n_0<=(P E_global N_base/eta)(25/24)^L`,

`sum_j R_j >= L log(48/25)-log(P E_global N_base/eta)`.

Therefore the same clean depth condition

`L log(36/25)>log(P E_global N_base/eta)`

forces at least one physical-transfer-weighted layer with

`R_j > log(4/3)`.

This removes the main caveat of the raw cycle-count theorem: the forced reuse signal can be measured with the causal transfer law itself.  What remains is a local conversion of this conditional reuse information into the existing collision-pair / Hodge-resistance / component-entropy currencies; one must not charge the baseline two parent roles as branching loss.

Stress: `{out.samples}`
- worst information-telescope residual: `{out.worst_telescope_residual:.3e}`
- minimum sampled conditional reuse information: `{out.minimum_conditional_information:.3e}`
- minimum rich-depth margin: `{out.minimum_depth_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
