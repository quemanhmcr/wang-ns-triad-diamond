from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict,dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

from src.single_edge_certificate import RSTAR_LO,RSTAR_HI

PARENT_LO=Fraction(3,5)
PARENT_HI=Fraction(5,8)
ROOT_MASS=Fraction(1,5)
REUSE_FRACTION=Fraction(1,4)


def reuse_savings(counts:Sequence[int])->list[int]:
    n=list(map(int,counts))
    if len(n)<2 or n[-1]!=1 or min(n)<=0: raise ValueError('positive layered counts ending at one terminal required')
    r=[]
    for a,b in zip(n[:-1],n[1:]):
        x=2*b-a
        if x<0: raise ValueError('causal layer has more distinct parents than two slots per child')
        r.append(x)
    return r


def reuse_fractions(counts:Sequence[int])->list[float]:
    n=list(map(int,counts)); r=reuse_savings(n)
    return [rr/(2.0*n[j+1]) for j,rr in enumerate(r)]


def incidence_cycle_rank(counts:Sequence[int])->int:
    """Cycle rank of a connected layered 3-uniform incidence graph."""
    n=list(map(int,counts)); reuse_savings(n)
    m=sum(n[1:]); packets=sum(n)
    return 2*m-packets+1


def causal_cycle_identity(counts:Sequence[int])->tuple[int,int]:
    r=reuse_savings(counts)
    return incidence_cycle_rank(counts),sum(r)


def binary_product_ratio(counts:Sequence[int])->float:
    n=list(map(int,counts)); rho=reuse_fractions(n)
    prod=float(np.prod([1.0-x for x in rho]))
    exact=n[0]/(2.0**(len(n)-1))
    if abs(prod-exact)>2e-13*max(1.0,exact): raise AssertionError('binary reuse product identity failed')
    return prod


def binary_reuse_action(counts:Sequence[int])->float:
    rho=reuse_fractions(counts)
    if any(x>=1 for x in rho): return math.inf
    return float(sum(-math.log1p(-x) for x in rho))


def root_frequency_upper(base_frequency:float,depth:int)->float:
    """Using child/base <=(5/3)^L and root/terminal <=(5/8)^L."""
    if base_frequency<=0 or depth<0: raise ValueError('bad scale data')
    return base_frequency*(25.0/24.0)**depth


def root_count_energy_upper(base_frequency:float,depth:int,global_energy:float,critical_mass:float=float(ROOT_MASS),frame_budget:float=1.0)->float:
    if min(base_frequency,global_energy,critical_mass,frame_budget)<=0 or depth<0: raise ValueError('positive data required')
    return frame_budget*global_energy*root_frequency_upper(base_frequency,depth)/critical_mass


def reuse_action_lower(base_frequency:float,depth:int,global_energy:float,critical_mass:float=float(ROOT_MASS),frame_budget:float=1.0)->float:
    """Lower action from binary formal leaves versus Moyal/Riesz root energy budget."""
    B=frame_budget*global_energy*base_frequency/critical_mass
    if B<=0 or depth<0: raise ValueError('bad data')
    return depth*math.log(48.0/25.0)-math.log(B)


def quarter_reuse_depth(base_frequency:float,global_energy:float,critical_mass:float=float(ROOT_MASS),frame_budget:float=1.0)->int:
    """First integer L with L log(36/25)>log(P E N0/eta)."""
    B=frame_budget*global_energy*base_frequency/critical_mass
    if B<=0: raise ValueError('positive budget required')
    if B<1: return 0
    x=math.log(B)/math.log(36.0/25.0)
    L=int(math.floor(x))+1
    while L*math.log(36.0/25.0)<=math.log(B): L+=1
    while L>0 and (L-1)*math.log(36.0/25.0)>math.log(B): L-=1
    return L


def classify_causal_ancestry(counts:Sequence[int],base_frequency:float,global_energy:float,critical_mass:float=float(ROOT_MASS),frame_budget:float=1.0)->dict[str,float|int|str]:
    n=list(map(int,counts)); L=len(n)-1; rho=reuse_fractions(n); action=binary_reuse_action(n); beta=incidence_cycle_rank(n)
    root_cap=root_count_energy_upper(base_frequency,L,global_energy,critical_mass,frame_budget)
    if n[0]>root_cap+1e-12*max(1.0,root_cap):
        return {'branch':'root_energy_exceeded','root_count':n[0],'root_count_upper':root_cap,'cycle_rank':beta,'reuse_action':action}
    if max(rho,default=0.0)>=.25-1e-15:
        j=int(np.argmax(rho)); return {'branch':'reuse_rich_layer','layer':j,'reuse_fraction':rho[j],'reuse_savings':reuse_savings(n)[j],'cycle_rank':beta,'reuse_action':action}
    stop=quarter_reuse_depth(base_frequency,global_energy,critical_mass,frame_budget)
    return {'branch':'shallow_binary_tree','depth':L,'quarter_reuse_forced_by_depth':stop,'root_count':n[0],'root_count_upper':root_cap,'cycle_rank':beta,'reuse_action':action}


def arb_signed_good_certificate()->dict[str,str]:
    try:
        from flint import arb,ctx
    except ImportError as exc: raise RuntimeError('python-flint required') from exc
    ctx.prec=180
    def aq(q): return arb(q.numerator)/q.denominator
    r=aq(RSTAR_LO).union(aq(RSTAR_HI)); lo=r*(-arb(1)/80).exp(); hi=r*(arb(1)/80).exp()
    if not (lo>arb(3)/5): raise AssertionError(f'parent lower ratio failed {lo}')
    if not (hi<arb(5)/8): raise AssertionError(f'parent upper ratio failed {hi}')
    if Fraction(2,1)/PARENT_HI != Fraction(16,5): raise AssertionError('binary/root shrink ratio algebra failed')
    if Fraction(16,5)/Fraction(5,3) != Fraction(48,25): raise AssertionError('distinguished-scale reuse slope failed')
    if Fraction(48,25)/Fraction(4,3) != Fraction(36,25): raise AssertionError('quarter-reuse stopping ratio failed')
    return {
        'signed_good_parent_child_ratio':'3/5 < N_parent/N_child < 5/8',
        'forward_scale_progress':'8/5 < N_child/N_parent < 5/3',
        'root_scale_relative_to_distinguished_base':'N_root,max < (25/24)^L N_base',
        'cycle_identity':'beta=sum_j(2 n_(j+1)-n_j)',
        'product_identity':'n_0/2^L=product_j(1-rho_j)',
        'reuse_action':'sum -log(1-rho_j) >= L log(48/25)-log(P E N0/eta)',
        'quarter_reuse_depth':'L log(36/25)>log(P E N0/eta) forces some rho_j>=1/4',
        'clean_root_mass':'eta=1/5 for separated selected affine grains',
        'status':'CERTIFIED_GIVEN_SYNCHRONIZED_CAUSAL_LAYERING_AND_COHERENT_ROOT_BUDGET',
    }


def random_counts(rng:np.random.Generator,L:int)->list[int]:
    n=[0]*(L+1); n[L]=1
    for j in range(L-1,-1,-1):
        n[j]=int(rng.integers(1,2*n[j+1]+1))
    return n

@dataclass(frozen=True)
class CausalAncestryStress:
    samples:int
    worst_cycle_identity_residual:int
    worst_product_residual:float
    minimum_action_margin:float
    minimum_quarter_stop_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260808)->CausalAncestryStress:
    rng=np.random.default_rng(seed); wc=0; wp=0.0; ma=mq=float('inf'); branches={}
    for _ in range(samples):
        L=int(rng.integers(1,18)); n=random_counts(rng,L); beta,sr=causal_cycle_identity(n); wc=max(wc,abs(beta-sr)); prod=binary_product_ratio(n); wp=max(wp,abs(prod-n[0]/2**L)); action=binary_reuse_action(n); exact=L*math.log(2)-math.log(n[0]); ma=min(ma,action-exact)
        if abs(action-exact)>3e-13*max(1.0,exact): raise AssertionError('reuse action telescope failed')
        N0=float(rng.lognormal(0,1)); E=float(rng.lognormal(0,1)); out=classify_causal_ancestry(n,N0,E); b=str(out['branch']); branches[b]=branches.get(b,0)+1
        stop=quarter_reuse_depth(N0,E)
        if L>=stop and n[0]<=root_count_energy_upper(N0,L,E)+1e-12 and max(reuse_fractions(n),default=0)<.25-1e-14:
            # theorem says this cannot happen once the strict stopping inequality holds
            if L*math.log(36/25)>math.log(5*N0*E)+1e-12:
                raise AssertionError('quarter reuse stopping theorem failed')
        mq=min(mq, L*math.log(36/25)-math.log(5*N0*E) if L>=stop else 0.0)
    return CausalAncestryStress(samples,wc,wp,ma,mq,branches)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-causal-binary-ancestry'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=arb_signed_good_certificate(); out=stress(args.samples)
    (args.outdir/'causal_binary_ancestry.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f'''# Causal binary ancestry: branching outruns scale dilation unless reuse is sticky

Status: **{cert['status']}**.

On the signed-good physical core, Arb certifies

`3/5 < N_parent/N_child < 5/8`, hence `8/5 < N_child/N_parent < 5/3`.

Consider a synchronized layered backward ancestry of depth `L`, with one terminal packet and every packet at level `j+1` generated by one triad using two parent incidences at level `j`.  Let `n_j` be the number of distinct packet ancestors.  Define

`r_j=2 n_(j+1)-n_j >=0`,  `rho_j=r_j/(2 n_(j+1))`.

Then the connected 3-uniform incidence graph has the exact cycle identity

`beta = sum_j r_j`,

and binary branching has the exact multiplicative reuse law

`n_0/2^L = product_j (1-rho_j)`.

Thus the cumulative reuse action is

`A_reuse=sum_j -log(1-rho_j)=L log 2-log n_0`.

If distinct root coherent cells each carry critical mass `N E>=eta` and their positive Moyal/frame energy has budget `P E_global`, signed-good scale synchronization gives

`n_0 <= (P E_global N_base/eta)(25/24)^L`.

Therefore

`A_reuse >= L log(48/25)-log(P E_global N_base/eta)`.

For the clean Moyal/fresh-grain constants `P=1`, `eta=1/5`, if every layer had `rho_j<1/4`, then `A_reuse<L log(4/3)`.  Hence once

`L log(36/25)>log(5 E_global N_base)`,

some layer must obey

`rho_j>=1/4`.

So a causally complete fresh ancestry cannot remain an almost-binary tree indefinitely: binary parent creation outruns the allowed scale dilation and finite coherent energy.  A sufficiently deep ancestry must contain a reuse-rich layer in which at least one quarter of parent slots are saved by merging/reuse; those savings are exact cycle-rank units of the layered incidence graph.

This theorem is deliberately conditional on a synchronized causal layering.  Establishing such a layering for the full Duhamel/packet PDE extraction, with summable asynchronous interfaces, remains a continuum step.

Stress: `{out.samples}`
- worst cycle-identity residual: `{out.worst_cycle_identity_residual}`
- worst binary-product residual: `{out.worst_product_residual:.3e}`
- minimum reuse-action telescope margin: `{out.minimum_action_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
