from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.single_edge_certificate import RSTAR_LO,RSTAR_HI
from src.spherical_erosion import C_STAR,KAPPA_STAR,normalize

LOCAL=Fraction(2,25)
TAU_CLEAN=Fraction(1,100)
KAPPA0_CLEAN=Fraction(17,100)


def half_cosine_from_uv(rstar:float,u:float,v:float)->float:
    """cos(theta/2) from normalized child length 1 and signed-good x,y."""
    if rstar<=0: raise ValueError('positive rstar required')
    c2=math.exp(2*v)/(4*rstar*rstar)-math.sinh(u/2)**2
    if c2<=0: raise ValueError('invalid triangle data')
    return math.sqrt(c2)


def half_cosine_error_clean(u:float,v:float)->float:
    return 1.2*abs(v)+0.2*u*u


def child_midpoint_error_clean(u:float)->float:
    return 0.5*abs(u)


def barycenter_step_error_upper(hodge_energy:float,parent_barycenter_mismatch:float)->float:
    """Clean error e in ||b_child-b_old/c_*||.

    H=E(2v^2+u^2/2) under normalized positive transfer weights.
    """
    if min(hodge_energy,parent_barycenter_mismatch)<0: raise ValueError('nonnegative inputs required')
    return 2.0*math.sqrt(hodge_energy)+0.5*hodge_energy+0.625*parent_barycenter_mismatch


def exact_potential_zeta_upper(step_error:float)->float:
    if not (0<=step_error<1): raise ValueError('step error must be in [0,1)')
    return -math.log1p(-step_error)


def clean_potential_zeta_upper(hodge_energy:float,parent_barycenter_mismatch:float)->float:
    e=barycenter_step_error_upper(hodge_energy,parent_barycenter_mismatch)
    if e>0.5+1e-15: raise ValueError('clean 2e bound requires e<=1/2')
    return 2.0*e


def tau_flat_zeta_upper(tau:float,parent_mismatch_ratio:float=1.0)->float:
    """If sqrt(H)<=tau/3 and Delta_b<=ratio*tau."""
    if tau<0 or parent_mismatch_ratio<0: raise ValueError('nonnegative inputs required')
    H=(tau/3.0)**2; d=parent_mismatch_ratio*tau
    return clean_potential_zeta_upper(H,d)


def tau_flat_kappa0_lower(tau:float,parent_mismatch_ratio:float=1.0)->float:
    return KAPPA_STAR-tau_flat_zeta_upper(tau,parent_mismatch_ratio)


def barycentric_potential_step_upper(parent_bary_norm:float,hodge_energy:float,parent_barycenter_mismatch:float)->float:
    """Upper P_child-P_parent on the concentrated branch |b_parent|>=c_*.

    Returns -kappa_*+zeta.  The exact zeta uses c_* e / |b|, while the clean
    bound below replaces it by e because |b|>=c_*.
    """
    if not (C_STAR<=parent_bary_norm<=1): raise ValueError('concentrated branch requires |b|>=c_*')
    e=barycenter_step_error_upper(hodge_energy,parent_barycenter_mismatch)
    if e>=parent_bary_norm/C_STAR: raise ValueError('error too large for positive child barycenter lower bound')
    z=-math.log1p(-C_STAR*e/parent_bary_norm)
    return -KAPPA_STAR+z


def _random_frame(rng:np.random.Generator)->tuple[np.ndarray,np.ndarray]:
    m=normalize(rng.normal(size=3)); d=rng.normal(size=3); d-=float(d@m)*m; d=normalize(d); return m,d


def actual_coupling_barycenters(rng:np.random.Generator,n:int,concentrated:bool=False):
    w=rng.dirichlet(np.ones(n)); p=[];q=[];child=[]; us=[];vs=[]
    r=float(rng.uniform(float(RSTAR_LO),float(RSTAR_HI)))
    for _ in range(n):
        u=float(rng.uniform(-float(LOCAL),float(LOCAL))); v=float(rng.uniform(-float(LOCAL),float(LOCAL))); c=half_cosine_from_uv(r,u,v); s=math.sqrt(max(0,1-c*c))
        if concentrated:
            # midpoints in a narrow cap around +e1
            m=normalize(np.array([1.0,float(rng.normal(scale=.08)),float(rng.normal(scale=.08))])); d=rng.normal(size=3); d-=float(d@m)*m; d=normalize(d)
        else: m,d=_random_frame(rng)
        pu=c*m+s*d; qu=c*m-s*d
        x=r*math.exp(-v-u/2); y=r*math.exp(-v+u/2)
        cu=normalize(x*pu+y*qu)
        p.append(pu);q.append(qu);child.append(cu);us.append(u);vs.append(v)
    p=np.array(p);q=np.array(q);child=np.array(child);us=np.array(us);vs=np.array(vs)
    b1=w@p; b2=w@q; bc=w@child; H=float(np.dot(w,2*vs*vs+.5*us*us)); d=float(np.linalg.norm(b2-b1))
    return b1,b2,bc,H,d,w,p,q,child,us,vs,r


def arb_geometry_certificate()->dict[str,str]:
    try:
        from flint import arb,ctx
    except ImportError as exc: raise RuntimeError('python-flint required') from exc
    ctx.prec=200
    def aq(q): return arb(q.numerator)/q.denominator
    r=aq(RSTAR_LO).union(aq(RSTAR_HI)); cs=1/(2*r); umax=arb(2)/25; vmax=arb(2)/25
    # local minimum half-cosine: v=-vmax, |u|=umax, worst r endpoint is enclosed automatically
    c2=(-2*vmax).exp()/(4*r*r)-(umax/2).sinh()**2
    if not (cs>arb(4)/5): raise AssertionError('c_* >4/5 failed')
    if not (c2>arb(9)/16): raise AssertionError(f'local c^2>9/16 failed: {c2}')
    denom=arb(31)/20 # c+c_*>3/4+4/5
    vcoef=(arb(7)/10)*2*(arb(4)/25).exp()/denom
    ucoef=(arb(1)/4)*(arb(1)/25).cosh()**2/denom
    if not (vcoef<arb(6)/5): raise AssertionError(f'half-angle v coefficient failed {vcoef}')
    if not (ucoef<arb(1)/5): raise AssertionError(f'half-angle u coefficient failed {ucoef}')
    # tau=1/100: zeta <=31 tau/12 + tau^2/9; kappa_*=-log c_*.
    tau=arb(1)/100; zeta=arb(31)/12*tau+tau*tau/9; kappa=(-cs.log())
    if not (kappa-zeta>arb(17)/100): raise AssertionError(f'physical flat erosion kappa0 failed: {kappa-zeta}')
    return {
        'cstar_lower':'4/5','local_half_cosine_lower':'3/4',
        'half_cosine_error':'|c_e-c_*| <= (6/5)|v|+(1/5)u^2',
        'child_midpoint_error':'||n_child-m_e|| <= |u|/2',
        'barycenter_step_error':'e <=2 sqrt(E_H)+(1/2)E_H+(5/8)Delta_b',
        'clean_zeta':'zeta <=4 sqrt(E_H)+E_H+(5/4)Delta_b',
        'tau':'1/100','tau_marginal_mismatch':'Delta_b<=tau',
        'zeta_tau':'<=31 tau/12+tau^2/9',
        'master_erosion_rate':'kappa0>17/100',
        'status':'CERTIFIED_PHYSICAL_BARYCENTRIC_PERTURBATION',
    }


@dataclass(frozen=True)
class PhysicalFlatEpisodeStress:
    samples:int
    worst_barycenter_error_ratio:float
    worst_potential_margin:float
    minimum_tau_kappa_margin:float


def stress(samples:int=50_000,seed:int=20260808)->PhysicalFlatEpisodeStress:
    rng=np.random.default_rng(seed); rb=0.0; pm=float('inf'); km=float('inf')
    for _ in range(samples):
        n=int(rng.integers(2,24)); b1,b2,bc,H,d,*_=actual_coupling_barycenters(rng,n,False); e=barycenter_step_error_upper(H,d); actual=float(np.linalg.norm(bc-b1/C_STAR)); rb=max(rb,actual/max(e,1e-300))
        if actual>e+3e-12: raise AssertionError('physical barycenter step bound failed')
    for _ in range(min(samples,8000)):
        n=int(rng.integers(2,20)); b1,b2,bc,H,d,*_=actual_coupling_barycenters(rng,n,True); bn=float(np.linalg.norm(b1)); e=barycenter_step_error_upper(H,d)
        if bn>=C_STAR and e<=.45 and np.linalg.norm(bc)>1e-14:
            P1=-math.log(bn); Pc=-math.log(float(np.linalg.norm(bc))); up=P1+barycentric_potential_step_upper(bn,H,d); pm=min(pm,up-Pc)
            if Pc>up+4e-12: raise AssertionError('barycentric potential perturbation failed')
    tau=.01; k=tau_flat_kappa0_lower(tau,1.0); km=k-.17
    if k<=.17: raise AssertionError('tau-flat erosion lost clean 0.17 rate')
    return PhysicalFlatEpisodeStress(samples,rb,pm,km)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-physical-flat-episode'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=arb_geometry_certificate(); out=stress(args.samples)
    (args.outdir/'physical_flat_episode.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f'''# Physical flat episode: explicit master perturbation and erosion rate

Status: **{cert['status']}**.

For a signed-good transfer edge, the exact half-angle formula is

`c_e^2 = c_*^2 exp(2v) - sinh(u/2)^2`.

Arb on the full local box certifies

`|c_e-c_*| <= (6/5)|v|+(1/5)u^2`,

while the exact unequal-parent midpoint decomposition gives

`||n_child-m_e|| <= |u|/2`.

For a normalized positive-transfer coupling with physical Hodge energy

`H=E[2v^2+u^2/2]`

and parent-barycenter mismatch `Delta_b=|b_2-b_1|`, these identities imply

`||b_child-b_1/c_*|| <= e`,

`e <= 2 sqrt(H)+(1/2)H+(5/8)Delta_b`.

On the concentrated master branch `|b_1|>=c_*`, the barycentric potential obeys

`P_child <= P_parent-kappa_* -log(1-e)`.

For `e<=1/2`, one may take the clean perturbation

`zeta <=4 sqrt(H)+H+(5/4)Delta_b`.

Combining with the service-or-flat theorem at `tau=1/100`, `sqrt(H)<=tau/3`, and the synchronized-marginal threshold `Delta_b<=tau`, gives

`zeta <=31 tau/12+tau^2/9`,

hence

`kappa_0=kappa_*-zeta >17/100`.

Thus a physical `1%`-Kelvin-flat, parent-synchronized block is already a quantitative master flat step with a uniform barycentric erosion rate; the abstract per-step `zeta_j` is no longer free.

Stress: `{out.samples}`
- worst actual barycenter error / clean bound: `{out.worst_barycenter_error_ratio:.9f}`
- minimum potential inequality margin: `{out.worst_potential_margin:.3e}`
- margin above clean `kappa_0=0.17`: `{out.minimum_tau_kappa_margin:.3e}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
