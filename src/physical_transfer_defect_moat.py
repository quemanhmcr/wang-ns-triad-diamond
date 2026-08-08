from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Hashable, Sequence

import numpy as np

from src.nested_grains import TriadEdge, incidence_components

GOOD_THRESHOLD = 1e-4
CLEAN_CHANGE_OF_MEASURE = 53.0 / 50.0
CLEAN_GOOD_CORE_MASS = 0.5
PHYSICAL_DEFECT_MEAN_FACTOR = 106.0 / 25.0
SQUARE_SCHEDULE_SUM_CLEAN = 13.0 / 20.0  # sum_{j>=0}(j+2)^-2 < 0.65


class UnionFind:
    def __init__(self, items: Sequence[Hashable]):
        self.parent={x:x for x in items}
        self.rank={x:0 for x in items}
    def find(self,x):
        p=self.parent[x]
        if p!=x:
            self.parent[x]=self.find(p)
        return self.parent[x]
    def union(self,a,b):
        a=self.find(a); b=self.find(b)
        if a==b: return
        if self.rank[a]<self.rank[b]: a,b=b,a
        self.parent[b]=a
        if self.rank[a]==self.rank[b]: self.rank[a]+=1


def physical_good_core_defect_mean_upper(block_deficit: float) -> float:
    """Mean scale defect under normalized physical child-transfer law on G_1e-4.

    Capacity law: sum w D <=2 epsilon from Def>=D/2.
    If epsilon<1/20000, good-core capacity mass >=1/2.  After normalization
    the mean is <=4 epsilon.  Physical/capacity density ratio <=53/50, hence
    <=(106/25) epsilon.
    """
    if block_deficit<0:
        raise ValueError('nonnegative block deficit required')
    return PHYSICAL_DEFECT_MEAN_FACTOR*block_deficit


def moat_cross_fraction_upper(mean_defect: float, radius: float, bins: int) -> float:
    if mean_defect<0 or radius<=0 or bins<1:
        raise ValueError('invalid moat parameters')
    return 1.0/bins + 2.0*mean_defect/radius


def square_schedule_total_upper(mean_defect_upper: float, M0: float, R0: float) -> float:
    if mean_defect_upper<0 or M0<=0 or R0<=0:
        raise ValueError('positive schedule bases required')
    return SQUARE_SCHEDULE_SUM_CLEAN*(1.0/M0+2.0*mean_defect_upper/R0)


def logarithmic_xi_upper(mean_defect_upper: float, M0: float, R0: float, cost_cap: float) -> float:
    """On low-cost blocks xi_j=log(1+eta_j e^Cj)<=e^Ccap eta_j."""
    if cost_cap<0:
        raise ValueError('nonnegative cost cap required')
    return math.exp(cost_cap)*square_schedule_total_upper(mean_defect_upper,M0,R0)


def choose_defect_moat(defects: Sequence[float], weights: Sequence[float], radius: float, bins: int) -> dict[str,float|int]:
    d=np.asarray(defects,float); w=np.asarray(weights,float)
    if d.shape!=w.shape or d.ndim!=1 or np.any(d<0) or np.any(w<0) or radius<=0 or bins<1:
        raise ValueError('invalid defect law')
    total=float(w.sum())
    if total<=0:
        raise ValueError('positive total transfer required')
    lo=radius/2.0; width=(radius/2.0)/bins
    masses=[]
    for j in range(bins):
        a=lo+j*width; b=a+width
        mask=(d>=a)&(d<b if j<bins-1 else d<=b)
        masses.append(float(w[mask].sum()))
    j=int(np.argmin(masses)); lower=lo+j*width; upper=lower+width
    return {'bin':j,'lower':lower,'upper':upper,'moat_mass':masses[j],'total':total}


def cross_component_mass(
    vertices: Sequence[tuple[Hashable,Hashable,Hashable]],
    defects: Sequence[float],
    weights: Sequence[float],
    lower: float,
    upper: float,
) -> dict[str,object]:
    if not (len(vertices)==len(defects)==len(weights)):
        raise ValueError('edge arrays must match')
    allv=list({x for e in vertices for x in e})
    uf=UnionFind(allv)
    for e,d in zip(vertices,defects):
        if d<lower:
            a,b,c=e; uf.union(a,b); uf.union(b,c)
    cross=moat=tail=0.0
    internal=[]
    for e,d,w in zip(vertices,defects,weights):
        roots={uf.find(x) for x in e}
        if len(roots)>1:
            cross+=w
            if lower<=d<upper: moat+=w
            elif d>=upper: tail+=w
            else: raise AssertionError('short cross edge should have connected its vertices')
        else:
            internal.append(TriadEdge(tuple(str(x) for x in e),float(d),float(w),float(w)))
    comps=incidence_components(internal) if internal else []
    return {'cross_mass':cross,'moat_cross_mass':moat,'tail_cross_mass':tail,'internal_edges':len(internal),'incidence':comps}


def defect_moat_certificate(
    vertices: Sequence[tuple[Hashable,Hashable,Hashable]],
    defects: Sequence[float],
    weights: Sequence[float],
    radius: float,
    bins: int,
) -> dict[str,object]:
    d=np.asarray(defects,float); w=np.asarray(weights,float); total=float(w.sum())
    if total<=0: raise ValueError('positive transfer required')
    mean=float(np.dot(d,w)/total)
    moat=choose_defect_moat(d,w,radius,bins)
    row=cross_component_mass(vertices,d,w,float(moat['lower']),float(moat['upper']))
    bound=moat_cross_fraction_upper(mean,radius,bins)*total
    if float(row['cross_mass'])>bound+3e-12*max(1.0,total):
        raise AssertionError('defect moat cross bound failed')
    return {'mean_defect':mean,'bound':bound,**moat,**row}


def theorem_certificate() -> dict[str,object]:
    return {
        'status':'EXACT_TRANSFER_WEIGHTED_DEFECT_MOAT_GIVEN_PHYSICAL_GOOD_CORE',
        'physical_mean':'E_phys D <= (106/25) epsilon on epsilon<1/20000 good core',
        'one_depth':'eta_cross <= 1/M + 2 E_phys[D]/R',
        'schedule':'M_j=M0(j+2)^2, R_j=R0(j+2)^2',
        'sum':'sum eta_j <= (13/20)(1/M0+2 Dbar/R0)',
        'log_xi':'on C_j<=Ccap, sum xi_j <= exp(Ccap) sum eta_j',
        'percolation':'retained connected components route through exact fresh-or-cycle incidence identity',
    }


@dataclass(frozen=True)
class MoatStress:
    samples:int
    minimum_cross_bound_margin:float
    minimum_markov_tail_margin:float
    minimum_moat_pigeonhole_margin:float
    minimum_schedule_margin:float
    incidence_checks:int


def stress(samples:int=50_000,seed:int=20260808)->MoatStress:
    rng=np.random.default_rng(seed)
    mc=mt=mm=ms=float('inf'); ic=0
    # Exact clean schedule sum check on large partial sums.
    partial=sum(1.0/(j+2)**2 for j in range(200000))
    if partial>=SQUARE_SCHEDULE_SUM_CLEAN:
        raise AssertionError('13/20 square schedule bound failed')
    for q in range(samples):
        nvert=int(rng.integers(5,35)); nedge=int(rng.integers(2,80))
        verts=[]
        for _ in range(nedge):
            e=tuple(rng.choice(nvert,size=3,replace=False).tolist())
            verts.append(e)
        # Wide adversarial defect distribution and highly nonuniform positive transfer.
        defects=rng.lognormal(mean=-2.0,sigma=2.0,size=nedge)
        weights=rng.lognormal(mean=-1.0,sigma=1.8,size=nedge)
        R=float(10**rng.uniform(-2,2)); M=int(rng.integers(1,80))
        row=defect_moat_certificate(verts,defects,weights,R,M)
        total=float(row['total']); mean=float(row['mean_defect'])
        mc=min(mc,float(row['bound'])-float(row['cross_mass']))
        # Tail cross is bounded by all tail mass, then Markov at upper>=R/2.
        tail_all=float(weights[defects>=float(row['upper'])].sum())
        mark=2*mean*total/R
        mt=min(mt,mark-tail_all)
        if tail_all>mark+3e-12*max(1.0,total): raise AssertionError('Markov tail bound failed')
        mm=min(mm,total/M-float(row['moat_mass']))
        if float(row['moat_mass'])>total/M+3e-12*max(1.0,total): raise AssertionError('moat pigeonhole failed')
        if q<5000:
            for comp in row['incidence']:
                ic+=1
                if int(comp['fresh_units'])+int(comp['cycle_rank'])!=2*int(comp['triads']):
                    raise AssertionError('fresh-cycle incidence identity failed after moat')
        Dbar=float(10**rng.uniform(-5,0)); M0=float(rng.uniform(2,100)); R0=float(rng.uniform(.2,100))
        finite=sum((1/M0+2*Dbar/R0)/(j+2)**2 for j in range(500))
        clean=square_schedule_total_upper(Dbar,M0,R0)
        ms=min(ms,clean-finite)
        if finite>clean+2e-12: raise AssertionError('summable cross schedule failed')
    return MoatStress(samples,mc,mt,mm,ms,ic)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-physical-transfer-defect-moat'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); out=stress(args.samples)
    (args.outdir/'physical_transfer_defect_moat.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f'''# Physical transfer-weighted defect moat\n\nStatus: **{cert['status']}**.\n\nOn the `eta_0=10^-4` signed-good core, capacity and actual positive child-transfer laws have density ratio at most `53/50`.  If the block deficit is `epsilon<1/20000`, capacity good-core mass is at least `1/2`, while single-edge stability gives `sum w_e D_e<=2 epsilon`. Hence the normalized **physical** good-core law satisfies\n\n`E_phys[D] <= (106/25) epsilon`.\n\nAt one recursive coherent-cell depth choose a scalar defect radius `R`, split `[R/2,R]` into `M` bins, and delete the bin with least physical transfer.  Connect packet/coherent vertices using all edges below the lower moat boundary. Every cross-component edge is then either in the chosen moat or has `D>=R/2`. Therefore\n\n`eta_cross <= 1/M + 2 E_phys[D]/R`.\n\nChoose `M_j=M0(j+2)^2` and `R_j=R0(j+2)^2`. Since `sum_(j>=0)(j+2)^-2<13/20`,\n\n`sum_j eta_cross,j <= (13/20)[1/M0+2 Dbar/R0]`.\n\nOn the low-cost extraction branch `C_j<=Ccap`, `xi_j=log(1+eta_j exp(C_j))<=exp(Ccap)eta_j`, so this gives a summable and tunably small physical-transfer `Xi`.\n\nNo phase-space packet count or Gaussian synthesis tail is used. If the low-defect hypergraph remains connected instead of splitting, that is not an interface failure: each retained connected component obeys the exact incidence identity `(n-1)+beta=2m` and therefore routes to fresh-rich or cycle-rich ancestry.\n\nStress: `{out.samples}` random weighted triad graphs\n- minimum cross-bound margin: `{out.minimum_cross_bound_margin:.3e}`\n- minimum Markov-tail margin: `{out.minimum_markov_tail_margin:.3e}`\n- minimum moat-pigeonhole margin: `{out.minimum_moat_pigeonhole_margin:.3e}`\n- minimum square-schedule margin: `{out.minimum_schedule_margin:.3e}`\n- incidence checks: `{out.incidence_checks}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
