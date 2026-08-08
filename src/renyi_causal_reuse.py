from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict,dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.atomic_component_entropy import collision_chain
from src.weighted_causal_reuse import parent_pushforward, random_layered_maps


def collision(weights:Sequence[float])->float:
    w=np.asarray(weights,float)
    if np.any(w<0) or w.sum()<=0: raise ValueError('positive law required')
    w=w/w.sum(); return float(np.dot(w,w))


def layer_collision_reuse(child_weights:Sequence[float],parent_pairs:np.ndarray,n_parents:int|None=None)->dict[str,float|np.ndarray]:
    w=np.asarray(child_weights,float); w=w/w.sum(); p=parent_pushforward(w,parent_pairs,n_parents)
    qc=collision(w); qp=collision(p); hidden=qp-.5*qc
    if hidden<-3e-14: raise AssertionError('pushforward reduced collision below binary-slot baseline')
    hidden=max(0.0,hidden); theta=2.0*hidden/qc
    lhs=qp; rhs=.5*qc*(1+theta)
    if abs(lhs-rhs)>3e-14*max(1.0,lhs): raise AssertionError('collision reuse recurrence failed')
    return {'q_child':qc,'q_parent':qp,'hidden_parent_slot_pair_mass':hidden,'theta':theta,'parent_weights':p,'child_h2':-math.log(qc),'parent_h2':-math.log(qp)}


def layered_collision_reuse(parent_maps:Sequence[np.ndarray])->dict[str,object]:
    maps=[np.asarray(x,int) for x in parent_maps]
    if not maps or maps[-1].shape[0]!=1: raise ValueError('one-terminal layered maps required')
    w=np.ones(1); rows=[None]*len(maps); actions=[0.0]*len(maps); qlaws=[None]*(len(maps)+1); qlaws[-1]=1.0
    for j in range(len(maps)-1,-1,-1):
        row=layer_collision_reuse(w,maps[j],int(maps[j].max())+1); rows[j]=row; actions[j]=math.log1p(float(row['theta'])); w=np.asarray(row['parent_weights']); qlaws[j]=float(row['q_parent'])
    total=float(sum(actions)); exact=len(maps)*math.log(2.0)+math.log(qlaws[0])
    if abs(total-exact)>5e-13*max(1.0,abs(exact)): raise AssertionError('Renyi causal reuse action telescope failed')
    return {'rows':rows,'action_by_layer':actions,'total_action':total,'root_collision':qlaws[0],'root_weights':w,'q_by_layer':qlaws}


def physical_action_lower(depth:int,base_frequency:float,global_energy:float,critical_mass:float=.2,frame_budget:float=1.0)->float:
    if depth<0 or min(base_frequency,global_energy,critical_mass,frame_budget)<=0: raise ValueError('bad physical data')
    return depth*math.log(48.0/25.0)-math.log(frame_budget*global_energy*base_frequency/critical_mass)


def rich_theta_threshold()->float:
    return 1.0/3.0


def rich_layer_route(child_weights:Sequence[float],parent_pairs:np.ndarray,child_ancestry_labels:Sequence[object]|None=None)->dict[str,float|str]:
    row=layer_collision_reuse(child_weights,parent_pairs)
    theta=float(row['theta']); qc=float(row['q_child']); hidden=float(row['hidden_parent_slot_pair_mass']); hchild=float(row['child_h2'])
    if theta<=1/3: return {'branch':'not_renyi_reuse_rich','theta':theta}
    if hchild<math.log(2.0):
        # qc>1/2 and hidden=(theta/2)qc >1/12.
        if hidden<=1/12-2e-14: raise AssertionError('clean hidden parent-slot pair mass failed')
        return {'branch':'weighted_parent_slot_reuse_pair','theta':theta,'hidden_pair_mass':hidden,'clean_hidden_pair_lower':1/12,'child_h2':hchild}
    if child_ancestry_labels is None:
        return {'branch':'child_atomic_collision_entropy','theta':theta,'child_h2':hchild,'clean_child_entropy_lower':math.log(2.0)}
    w=np.asarray(child_weights,float); w=w/w.sum(); chain=collision_chain(w,child_ancestry_labels)
    if chain['h_ancestry']>=.5*math.log(2.0)-1e-13:
        return {'branch':'child_component_Bellman_entropy','theta':theta,'child_h2':hchild,'ancestry_h2':chain['h_ancestry'],'clean_ancestry_entropy_lower':.5*math.log(2.0)}
    req=1/math.sqrt(2.0)-.5
    if chain['hidden_pair_mass']<req-2e-13: raise AssertionError('child same-ancestry pair lower failed')
    return {'branch':'child_same_ancestry_pair_cycle','theta':theta,'child_h2':hchild,'hidden_pair_mass':chain['hidden_pair_mass'],'clean_hidden_pair_lower':.2,'sharp_clean_formula':req}


def depth_forces_renyi_rich(depth:int,base_frequency:float,global_energy:float,critical_mass:float=.2,frame_budget:float=1.0)->bool:
    B=frame_budget*global_energy*base_frequency/critical_mass
    return depth*math.log(36.0/25.0)>math.log(B)


def random_ancestry_labels(rng:np.random.Generator,n:int)->list[int]:
    k=int(rng.integers(1,max(2,n+1))); return rng.integers(0,k,size=n).tolist()

@dataclass(frozen=True)
class RenyiCausalStress:
    samples:int
    worst_recurrence_residual:float
    worst_action_telescope_residual:float
    minimum_rich_action_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260808)->RenyiCausalStress:
    rng=np.random.default_rng(seed); wr=wa=0.0; mm=float('inf'); branches={}
    for _ in range(samples):
        L=int(rng.integers(1,15)); maps=random_layered_maps(rng,L); out=layered_collision_reuse(maps); exact=L*math.log(2)+math.log(float(out['root_collision'])); wa=max(wa,abs(float(out['total_action'])-exact))
        for row in out['rows']:
            wr=max(wr,abs(float(row['q_parent'])-.5*float(row['q_child'])*(1+float(row['theta']))))
        # stress local routing by finding max-theta layer and adding random ancestry labels
        j=int(np.argmax([float(r['theta']) for r in out['rows']])); row=out['rows'][j]; childn=maps[j].shape[0]
        # reconstruct child law from q_by recursion by rerunning from terminal up to j
        w=np.ones(1)
        laws=[None]*(L+1); laws[L]=w
        for k in range(L-1,-1,-1): w=parent_pushforward(w,maps[k],int(maps[k].max())+1); laws[k]=w
        route=rich_layer_route(laws[j+1],maps[j],random_ancestry_labels(rng,childn)); b=str(route['branch']); branches[b]=branches.get(b,0)+1
        N=float(rng.lognormal(0,1)); E=float(rng.lognormal(0,1));
        if depth_forces_renyi_rich(L,N,E):
            lo=physical_action_lower(L,N,E); mm=min(mm,lo-L*math.log(4/3))
    return RenyiCausalStress(samples,wr,wa,mm,branches)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-renyi-causal-reuse'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); out=stress(args.samples)
    data={'status':'EXACT_RENYI_CAUSAL_REUSE_TO_EXISTING_CURRENCIES','theorem':{
        'layer_recurrence':'Q_parent=(1/2)Q_child+R_hidden=(1/2)Q_child(1+theta)',
        'action_telescope':'sum log(1+theta_j)=L log2+log Q_root',
        'physical_lower':'action >= L log(48/25)-log(P E_global N_base/eta)',
        'rich_layer':'depth threshold forces theta_j>1/3',
        'rich_route':'H2_child>=log2 -> component entropy >=(1/2)log2 or same-ancestry pair >1/5; else hidden parent-slot pair >1/12',
        'role_baseline':'the factor 1/2 is the exact cost-free two-parent-role baseline',
    },'stress':asdict(out)}
    (args.outdir/'renyi_causal_reuse.json').write_text(json.dumps(data,indent=2))
    md=f'''# Renyi causal reuse: weighted binary ancestry reaches existing Bellman/cycle currencies

Status: **EXACT_RENYI_CAUSAL_REUSE_TO_EXISTING_CURRENCIES**.

For a causal child law `w`, duplicating each event into its two parent-role slots gives collision probability `Q_slot=Q_child/2`.  Pushforward through the physical parent label map has

`Q_parent = Q_child/2 + R_hidden`,

where `R_hidden` is exactly the weighted hidden pair mass of distinct causal slots sharing a parent.  Define

`theta = 2 R_hidden / Q_child`.

Then

`Q_parent=(Q_child/2)(1+theta)`.

Across a one-terminal depth-`L` ancestry,

`sum_j log(1+theta_j)=L log2+log Q_root`.

Using `Q_root>=1/n_0` and the coherent root-energy / signed-good scale bound gives

`sum_j log(1+theta_j) >= L log(48/25)-log(P E_global N_base/eta)`.

Therefore the clean depth condition `L log(36/25)>log(P E_global N_base/eta)` forces some layer with

`theta_j>1/3`.

At that layer there is a uniform existing-currency route:

- if `H2_child<log2`, then `Q_child>1/2` and `R_hidden>1/12`: a physical transfer-weighted parent-slot reuse pair;
- if `H2_child>=log2`, apply the existing atomic-to-ancestry collision chain.  Either ancestry/component collision entropy is at least `(1/2)log2`, or child same-ancestry hidden pair mass exceeds `1/sqrt(2)-1/2>1/5`.

Thus the baseline two parent roles remain free, but any sufficiently deep causal ancestry must eventually pay a **transfer-weighted** Bellman/component entropy or pair/cycle currency.  The remaining PDE issue is constructing the synchronized causal layers and registering these parent labels with the same coherent/nested ancestry used by the master.

Stress: `{out.samples}`
- worst layer recurrence residual: `{out.worst_recurrence_residual:.3e}`
- worst action telescope residual: `{out.worst_action_telescope_residual:.3e}`
- minimum rich-depth action margin: `{out.minimum_rich_action_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
