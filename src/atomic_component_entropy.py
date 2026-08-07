from __future__ import annotations

import argparse, json, math
from collections import defaultdict, deque
from pathlib import Path
from typing import Iterable
import numpy as np


def normalize(weights):
    w=np.asarray(weights,float)
    if np.any(w<0) or w.sum()<=0: raise ValueError('bad weights')
    return w/w.sum()


def collision_chain(weights, labels):
    w=normalize(weights)
    labels=list(labels)
    if len(w)!=len(labels): raise ValueError('length mismatch')
    masses=defaultdict(float)
    intern=defaultdict(float)
    for wi,a in zip(w,labels):
        masses[a]+=float(wi)
        intern[a]+=float(wi*wi)
    q_atom=float(np.dot(w,w))
    q_anc=float(sum(x*x for x in masses.values()))
    hidden=q_anc-q_atom
    # exact decomposition hidden = sum_A W_A^2 (1-q_A)
    decomp=0.0
    q_cond={}
    for a,W in masses.items():
        qa=intern[a]/(W*W)
        q_cond[a]=qa
        decomp += W*W*(1.0-qa)
    h_atom=-math.log(q_atom)
    h_anc=-math.log(q_anc)
    hidden_entropy=h_atom-h_anc
    conditional_mean=q_atom/q_anc
    return dict(q_atom=q_atom,q_ancestry=q_anc,hidden_pair_mass=hidden,
                hidden_pair_decomposition=decomp,h_atom=h_atom,h_ancestry=h_anc,
                hidden_entropy=hidden_entropy,conditional_collision_mean=conditional_mean,
                masses=dict(masses),q_cond=q_cond)


def entropy_or_pair_bound(h: float, alpha: float=0.5):
    if h<0 or not (0<alpha<1): raise ValueError
    return math.exp(-alpha*h)-math.exp(-h)


def verify_entropy_pair_dichotomy(weights, labels, alpha=0.5):
    d=collision_chain(weights,labels)
    h=d['h_atom']
    if d['h_ancestry'] >= alpha*h - 1e-14:
        return True, 'component-entropy', 0.0
    req=entropy_or_pair_bound(h,alpha)
    return d['hidden_pair_mass']+1e-14 >= req, 'pair-cycle', req


def effective_support(weights):
    w=normalize(weights)
    return 1.0/float(np.dot(w,w))


def core_mass_bound(weights, lam=4.0):
    """Size-biased Markov core: mass on atoms w_i <= lam*Q is >= 1-1/lam."""
    w=normalize(weights); q=float(np.dot(w,w)); mask=w<=lam*q+1e-15
    return float(w[mask].sum()), int(mask.sum()), q


def graph_cycle_rank(vertices: Iterable[str], edges: Iterable[tuple[str,str]]) -> int:
    V=set(vertices); E=list(edges)
    adj={v:[] for v in V}
    for a,b in E:
        V.add(a); V.add(b)
        adj.setdefault(a,[]).append(b); adj.setdefault(b,[]).append(a)
    seen=set(); comps=0
    for v in V:
        if v in seen: continue
        comps+=1; dq=[v]; seen.add(v)
        while dq:
            x=dq.pop()
            for y in adj.get(x,[]):
                if y not in seen: seen.add(y); dq.append(y)
    return len(E)-len(V)+comps


def ancestry_contracted_rank(triads, ancestry_labels, fresh_token='FRESH'):
    """Exact cycle-rank gain after contracting each old ancestry component.

    triads: iterable of triples of packet names. ancestry_labels maps packet->old component label;
    fresh packets should have unique labels or fresh_token+packet name.
    Each triad gets its own incidence node. Parallel edges are retained.
    """
    triads=list(triads)
    labels={}
    for p,a in ancestry_labels.items():
        labels[p]=(fresh_token+':'+p) if a==fresh_token else ('A:'+str(a))
    vertices=set(labels.values())
    edges=[]
    for j,t in enumerate(triads):
        tn='T:'+str(j); vertices.add(tn)
        for p in t:
            if p not in labels: raise ValueError(f'missing label {p}')
            edges.append((tn,labels[p]))
    return graph_cycle_rank(vertices,edges)


def raw_incidence_rank(triads):
    triads=list(triads); vertices=set(); edges=[]
    for j,t in enumerate(triads):
        tn='T:'+str(j); vertices.add(tn)
        for p in t:
            pn='P:'+p; vertices.add(pn); edges.append((tn,pn))
    return graph_cycle_rank(vertices,edges)


def ancestry_cycle_gain(triads, ancestry_labels, fresh_token='FRESH'):
    raw=raw_incidence_rank(triads)
    contracted=ancestry_contracted_rank(triads,ancestry_labels,fresh_token)
    groups=defaultdict(set)
    for p,a in ancestry_labels.items():
        if a!=fresh_token: groups[a].add(p)
    attachment_gain=sum(max(0,len(s)-1) for s in groups.values())
    return dict(raw_rank=raw,contracted_rank=contracted,
                rank_gain=contracted-raw,attachment_lower_bound=attachment_gain)



def pair_biased_multiplicity_certificate(weights, labels, lam=2.0):
    """Exact entropy-to-reused-multiplicity certificate.

    Under alpha_A proportional to W_A^2, E_alpha q_A = exp(-(H_at-H_anc)).
    Hence alpha{q_A <= lam exp(-d)} >= 1-1/lam; and since q_A >= 1/k_A,
    those ancestry classes contain k_A >= exp(d)/lam distinct atoms.
    """
    if lam <= 1: raise ValueError('lam must exceed 1')
    d=collision_chain(weights,labels)
    masses=d['masses']; q=d['q_cond']; qanc=d['q_ancestry']; gap=d['hidden_entropy']
    alpha={a:(W*W/qanc) for a,W in masses.items()}
    threshold=lam*math.exp(-gap)
    good=[a for a in alpha if q[a] <= threshold + 1e-15]
    good_mass=sum(alpha[a] for a in good)
    # actual atom multiplicities
    counts=defaultdict(int)
    for a in labels: counts[a]+=1
    min_actual=min((counts[a] for a in good), default=math.inf)
    return dict(hidden_entropy=gap,pair_biased_good_mass=good_mass,
                theorem_good_mass=1.0-1.0/lam,
                multiplicity_lower_bound=math.exp(gap)/lam,
                minimum_actual_multiplicity=min_actual,good_labels=good)

def random_trial(rng, n=24, k=6):
    w=rng.dirichlet(np.full(n,0.45))
    labels=rng.integers(0,k,size=n).tolist()
    d=collision_chain(w,labels)
    ok,mode,req=verify_entropy_pair_dichotomy(w,labels,0.5)
    core=core_mass_bound(w,4.0)
    return dict(ok=ok, mode=mode, required=req, **{x:d[x] for x in ['q_atom','q_ancestry','hidden_pair_mass','h_atom','h_ancestry']}, core_mass=core[0])


def example_graphs():
    # Three reused atoms from ancestry A are attached by a current triad chain.
    triads=[('a','x','m'),('b','m','n'),('c','n','d')]
    labels={p:'FRESH' for t in triads for p in t}
    for p in ('a','b','c'): labels[p]='oldA'
    ex1=ancestry_cycle_gain(triads,labels)
    # Two old ancestry groups, two attachments each.
    triads2=[('a','x','m'),('b','m','n'),('c','y','r'),('d','r','s'),('n','s','z')]
    labels2={p:'FRESH' for t in triads2 for p in t}
    for p in ('a','b'): labels2[p]='oldA'
    for p in ('c','d'): labels2[p]='oldB'
    ex2=ancestry_cycle_gain(triads2,labels2)
    return ex1,ex2


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50000); ap.add_argument('--outdir',type=Path,default=Path('results-atomic-component')); args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20260807)
    rows=[random_trial(rng,n=int(rng.integers(8,50)),k=int(rng.integers(2,10))) for _ in range(args.samples)]
    violations=sum(not r['ok'] for r in rows)
    min_core=min(r['core_mass'] for r in rows)
    modes=defaultdict(int)
    for r in rows: modes[r['mode']]+=1
    ex1,ex2=example_graphs()
    # Equal atom, grouped examples exhibit exact hidden pair mass.
    examples=[]
    for groupsize in [1,2,4,8]:
        n=16; labels=[i//groupsize for i in range(n)]; w=np.ones(n)/n; d=collision_chain(w,labels)
        examples.append(dict(groupsize=groupsize, groups=len(set(labels)), q_atom=d['q_atom'],q_ancestry=d['q_ancestry'],hidden=d['hidden_pair_mass'],h_atom=d['h_atom'],h_ancestry=d['h_ancestry']))
    out=dict(samples=args.samples,violations=violations,min_core_mass_lambda4=min_core,modes=dict(modes),equal_examples=examples,graph_example_one=ex1,graph_example_two=ex2)
    (args.outdir/'atomic_component.json').write_text(json.dumps(out,indent=2))
    lines=['# Atomic-to-component entropy transfer','',f'Random dichotomy checks: `{args.samples}`; violations: `{violations}`.',f'Minimum size-biased core mass for lambda=4: `{min_core:.12f}` (theorem bound 0.75).','', '## Equal-weight grouped atoms','', '| atoms per ancestry | ancestry groups | Q_atom | Q_ancestry | hidden pair mass |', '|---:|---:|---:|---:|---:|']
    for x in examples: lines.append(f"| {x['groupsize']} | {x['groups']} | {x['q_atom']:.6f} | {x['q_ancestry']:.6f} | {x['hidden']:.6f} |")
    lines += ['', '## Ancestry contraction cycle certificates','',f"Example 1: `{ex1}`",f"Example 2: `{ex2}`",'', 'The collision chain rule and ancestry-contraction rank identities are exact. Random checks only test implementation.']
    (args.outdir/'summary.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))

if __name__=='__main__': main()
