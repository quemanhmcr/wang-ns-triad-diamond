from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import asdict
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

from .helical import coupling_magnitude_closed, diamond_metrics


def single_J(v, signs):
    x,y=v
    if not (0 < x <= y < 1 and x+y > 1):
        return 0.0
    sx,sy,sz=signs
    area_term=(x+y+1)*(-x+y+1)*(x-y+1)*(x+y-1)
    if area_term <= 0:
        return 0.0
    g=coupling_magnitude_closed(x,y,1.0,sx,sy,sz)
    return math.log(1/y)*abs(sx*x-sy*y)*g


def find_single_global(quick=False):
    best=None
    for signs in itertools.product((-1,1), repeat=3):
        def obj(v):
            x,y=v
            penalty=0.0
            if x>y: penalty += 20*(x-y)**2
            if x+y<1: penalty += 20*(1-x-y)**2
            return -single_J((min(x,y),max(x,y)), signs)+penalty
        r=differential_evolution(obj,[(0.05,0.9999),(0.05,0.9999)],seed=410+hash(signs)%1000,
                                 maxiter=70 if quick else 350,popsize=10 if quick else 24,
                                 tol=1e-10,polish=True)
        x,y=sorted(r.x)
        rec={'J':single_J((x,y),signs),'x':float(x),'y':float(y),'signs':list(signs)}
        if best is None or rec['J']>best['J']:
            best=rec
    return best


def vectors(params):
    log_rb,log_rc,theta,phi,psi=params
    rb,rc=math.exp(log_rb),math.exp(log_rc)
    a=np.array([1.0,0.0,0.0])
    b=rb*np.array([math.cos(theta),math.sin(theta),0.0])
    c=rc*np.array([math.sin(phi)*math.cos(psi),math.sin(phi)*math.sin(psi),math.cos(phi)])
    return a,b,c


def softmin(vals,tau=30.0):
    vals=np.asarray(vals,float)
    return -math.log(np.exp(-tau*vals).sum())/tau


def optimise_motif(jstar,edge_names,quick=False):
    bounds=[(math.log(.25),math.log(4.0)),(math.log(.25),math.log(4.0)),
            (.03,math.pi-.03),(.03,math.pi-.03),(0,2*math.pi)]
    best=None
    for idx,signs in enumerate(itertools.product((-1,1),repeat=6)):
        def obj(p):
            a,b,c=vectors(p)
            try: dm=diamond_metrics(a,b,c,signs)
            except Exception: return 100.0
            es=[dm['edges'][name] for name in edge_names]
            if any(e.forward_ratio <= 1.0002 for e in es):
                return 20+sum(max(0,1.0002-e.forward_ratio) for e in es)
            ratios=[e.efficiency/jstar for e in es]
            return -softmin(ratios)+0.001*(p[0]**2+p[1]**2)
        r=differential_evolution(obj,bounds,seed=8000+idx,maxiter=55 if quick else 260,
                                 popsize=8 if quick else 18,tol=1e-8,polish=True)
        a,b,c=vectors(r.x); dm=diamond_metrics(a,b,c,signs)
        es=[dm['edges'][name] for name in edge_names]
        rec={
            'signs':list(signs),'params':[float(x) for x in r.x],
            'min_ratio':float(min(e.efficiency for e in es)/jstar),
            'ratios':{name:float(dm['edges'][name].efficiency/jstar) for name in edge_names},
            'forward_ratios':{name:float(dm['edges'][name].forward_ratio) for name in edge_names},
            'lengths':{k:float(np.linalg.norm(v)) for k,v in dm['vectors'].items()},
            'phase_frustration':float(dm['phase_frustration']),
        }
        if best is None or rec['min_ratio']>best['min_ratio']:
            best=rec
        print(json.dumps({'motif':edge_names,'i':idx,'best':rec['min_ratio']}),flush=True)
    return best


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--quick',action='store_true'); ap.add_argument('--outdir',default='results-scale')
    args=ap.parse_args(); out=Path(args.outdir); out.mkdir(exist_ok=True,parents=True)
    single=find_single_global(args.quick); jstar=single['J']
    tri=optimise_motif(jstar,['ab_m','mc_d','bc_n'],args.quick)
    diamond=optimise_motif(jstar,['ab_m','mc_d','bc_n','an_d'],args.quick)
    payload={'single_global':single,'three_edge_reuse':tri,'four_edge_diamond':diamond,'quick':args.quick}
    (out/'scale_gap.json').write_text(json.dumps(payload,indent=2))
    gap3=1-tri['min_ratio']; gap4=1-diamond['min_ratio']
    md=f'''# Scale-holonomy experiment

- Global single-edge numerical maximum: `{jstar:.12f}` at parent ratios `x={single['x']:.9f}`, `y={single['y']:.9f}`, signs `{single['signs']}`.
- Best three-edge reuse motif: min ratio `{tri['min_ratio']:.9f}`; numerical deficit `{gap3:.9f}`.
- Best four-edge diamond: min ratio `{diamond['min_ratio']:.9f}`; numerical deficit `{gap4:.9f}`.

These are non-certified numerical values. They identify candidate constants for a later interval-arithmetic proof; they are not theorem statements.
'''
    (out/'summary.md').write_text(md)
    print(md)

if __name__=='__main__': main()
