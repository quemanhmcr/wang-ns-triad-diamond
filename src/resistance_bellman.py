from __future__ import annotations

import argparse, json, math
from dataclasses import dataclass
from pathlib import Path
import numpy as np

@dataclass(frozen=True)
class Edge:
    u: int
    v: int
    c: float


def normalize_edges(edges):
    s=sum(e.c for e in edges)
    if s<=0: raise ValueError('positive total conductance required')
    return [Edge(e.u,e.v,e.c/s) for e in edges]


def adjacency(n, edges):
    a=[[] for _ in range(n)]
    for k,e in enumerate(edges):
        r=1.0/e.c
        a[e.u].append((e.v,r,k)); a[e.v].append((e.u,r,k))
    return a


def tree_resistance_matrix(n, edges):
    a=adjacency(n,edges)
    R=np.zeros((n,n))
    for s in range(n):
        stack=[(s,-1,0.0)]
        while stack:
            u,par,d=stack.pop(); R[s,u]=d
            for v,r,_ in a[u]:
                if v!=par: stack.append((v,u,d+r))
    return R


def component_collision(n, edges, cutmask, p):
    parent=list(range(n))
    def find(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b: parent[b]=a
    cross=0.0
    for k,e in enumerate(edges):
        if (cutmask>>k)&1: cross += e.c
        else: union(e.u,e.v)
    masses={}
    for i,w in enumerate(p): masses[find(i)]=masses.get(find(i),0.0)+w
    q=sum(x*x for x in masses.values())
    return q,cross


def poisson_certificate(edges, p, R, lam):
    p=np.asarray(p,float); p=p/p.sum()
    pair=np.outer(p,p)
    F=float(np.sum(pair*(1.0-np.exp(-R/lam))))
    expQ=1.0-F
    expC=sum(e.c*(1.0-math.exp(-1.0/(e.c*lam))) for e in edges)
    rho=sum(min(e.c,1.0/lam) for e in edges)
    if F<=0:
        K=float('inf'); qbound=1.0; cbound=float('inf')
    else:
        K=2.0*(2.0-F)/F
        qbound=1.0-F/2.0
        cbound=K*expC
    edgeQ=sum(e.c*e.c for e in edges)
    edgeH=-math.log(edgeQ)
    edgeQbound=max(0.0,min(1.0,1.0-rho+1.0/lam))
    edgeHbound=float('inf') if edgeQbound==0 else -math.log(edgeQbound)
    return dict(F=F, expQ=expQ, expC=expC, rho=rho, K=K,
                qbound=qbound, cbound=cbound, edgeQ=edgeQ, edgeH=edgeH,
                edgeQbound=edgeQbound, edgeHbound=edgeHbound)


def weighted_median_resistance(R,p):
    p=np.asarray(p,float); p=p/p.sum(); vals=[]
    for i in range(len(p)):
        for j in range(len(p)):
            vals.append((float(R[i,j]),float(p[i]*p[j])))
    vals.sort()
    acc=0.0
    for r,w in vals:
        acc+=w
        if acc>=0.5: return max(r,1e-12)
    return vals[-1][0]


def exhaustive_witness(n,edges,p,cert):
    if len(edges)>18: return None
    best=None
    for mask in range(1<<len(edges)):
        q,c=component_collision(n,edges,mask,p)
        if q<=cert['qbound']+1e-12 and c<=cert['cbound']+1e-12:
            best=(q,c,mask); break
    return best


def random_tree(rng,n):
    edges=[]
    for v in range(1,n):
        u=int(rng.integers(0,v)); c=float(np.exp(rng.normal()))
        edges.append(Edge(u,v,c))
    return normalize_edges(edges)


def theorem_random_checks(samples=1000,seed=0):
    rng=np.random.default_rng(seed); violations=0; min_entropy_margin=1e9; min_witness_margin=1e9
    for _ in range(samples):
        n=int(rng.integers(3,10)); edges=random_tree(rng,n)
        p=rng.dirichlet(np.ones(n)); R=tree_resistance_matrix(n,edges)
        lam=weighted_median_resistance(R,p)
        cert=poisson_certificate(edges,p,R,lam)
        # exact edge entropy inequality Q_edge <= 1-rho+1/lambda
        min_entropy_margin=min(min_entropy_margin, cert['edgeQbound']-cert['edgeQ'])
        wit=exhaustive_witness(n,edges,p,cert)
        if wit is None:
            violations+=1
        else:
            q,c,_=wit
            min_witness_margin=min(min_witness_margin, cert['qbound']-q, cert['cbound']-c)
    return dict(samples=samples,violations=violations,min_edge_entropy_margin=min_entropy_margin,min_witness_margin=min_witness_margin)


def toy(kind,n=16):
    if kind=='star':
        edges=[Edge(0,i,1.0) for i in range(1,n)]
        p=np.array([0.0]+[1/(n-1)]*(n-1))
    elif kind=='path':
        edges=[Edge(i,i+1,1.0) for i in range(n-1)]
        p=np.ones(n)/n
    elif kind=='bottleneck':
        # two short dense-in-tree arms joined by a tiny conductance edge
        edges=[]
        mid=n//2
        for i in range(mid-1): edges.append(Edge(i,i+1,10.0))
        edges.append(Edge(mid-1,mid,0.05))
        for i in range(mid,n-1): edges.append(Edge(i,i+1,10.0))
        p=np.ones(n)/n
    else: raise ValueError(kind)
    edges=normalize_edges(edges); R=tree_resistance_matrix(n,edges); lam=weighted_median_resistance(R,p)
    cert=poisson_certificate(edges,p,R,lam)
    return dict(kind=kind,n=n,lambda_median=lam,**cert)


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=2000); ap.add_argument('--outdir',type=Path,default=Path('results-resistance-bellman'))
    a=ap.parse_args(); a.outdir.mkdir(parents=True,exist_ok=True)
    rnd=theorem_random_checks(a.samples,17)
    toys=[toy(k,18) for k in ['star','path','bottleneck']]
    F0=0.5*(1.0-math.exp(-1.0)); h0=-math.log(1.0-F0/2.0); K0=2*(2-F0)/F0
    out=dict(random_checks=rnd,toys=toys,median_constants=dict(F0=F0,h0=h0,K0=K0))
    (a.outdir/'resistance_bellman.json').write_text(json.dumps(out,indent=2))
    lines=['# Resistance-to-Bellman Poisson stopping','',
           f"Random exact checks: `{rnd['samples']}`; violations: `{rnd['violations']}`.",
           f"Minimum edge-entropy margin: `{rnd['min_edge_entropy_margin']:.3e}`.",
           f"Minimum exhaustive-witness margin: `{rnd['min_witness_margin']:.3e}`.",'',
           '## Median constants','',
           f"F0: `{F0:.12f}`",f"Bellman entropy floor h0: `{h0:.12f}`",f"simultaneous-cut multiplier K0: `{K0:.12f}`",'',
           '## Toy trees','',
           '| type | median R | F | rho | cut bound | component-Q bound | edge H | edge-H lower bound |',
           '|:--|--:|--:|--:|--:|--:|--:|--:|']
    for t in toys:
        lines.append(f"| {t['kind']} | {t['lambda_median']:.6g} | {t['F']:.6f} | {t['rho']:.6f} | {t['cbound']:.6f} | {t['qbound']:.6f} | {t['edgeH']:.6f} | {t['edgeHbound']:.6f} |")
    lines += ['', 'The Poisson collision identity, simultaneous stopping bound, and edge-conductance collision bound are exact finite-dimensional statements. Random/exhaustive checks only validate implementation.']
    (a.outdir/'summary.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))
if __name__=='__main__': main()
