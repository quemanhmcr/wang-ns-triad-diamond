from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.affine_coherent_moyal import periodic_discrete_stft
from src.atomic_component_entropy import collision_chain


def cubic_to_square_threshold(cubic_charge: float, filter_l1: float, lp_constant: float, bernstein_constant: float) -> float:
    """Y=(Q/g1)^(2/3)/(C_LP C_B)^2.

    Since Q/g1 is an average of ||delta_r u||_3^3, some filter displacement r
    has ||delta_r u||_3^2 >= (Q/g1)^(2/3).  LP+Bernstein then force the total
    square service sum_j M_j ||delta_r u_j||_2^2 to be at least Y.
    """
    if min(cubic_charge, filter_l1, lp_constant, bernstein_constant) <= 0:
        raise ValueError('positive parameters required')
    return (cubic_charge / filter_l1) ** (2.0 / 3.0) / (lp_constant * bernstein_constant) ** 2


def low_square_service_lower(square_threshold: float, high_normalized_enstrophy: float) -> float:
    """Low-band actual increment service after high-band estimate <=2 d_high."""
    if min(square_threshold, high_normalized_enstrophy) < 0:
        raise ValueError('nonnegative inputs required')
    return max(0.0, square_threshold - 2.0 * high_normalized_enstrophy)


def periodic_increment_covariance_residual(f: np.ndarray, g: np.ndarray, shift: int) -> float:
    """Exact discrete analogue of V_g(delta_r f)(X,k)=e^{-ikr}V_g f(X-r,k)-V_g f(X,k)."""
    f=np.asarray(f,complex); g=np.asarray(g,complex); n=len(f)
    if f.ndim!=1 or g.shape!=f.shape: raise ValueError('same 1D shape required')
    A=periodic_discrete_stft(f,g)
    D=periodic_discrete_stft(np.roll(f,shift)-f,g)
    k=np.arange(n)
    phase=np.exp(-2j*np.pi*k*shift/n)
    pred=np.empty_like(A)
    for m in range(n): pred[m]=phase*A[(m-shift)%n]-A[m]
    return float(np.max(np.abs(D-pred)))


def discrete_cell_increment_energies(f: np.ndarray, g: np.ndarray, shift: int, labels: np.ndarray) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    """Return exact increment energy and the two neighboring coherent energies per cell."""
    f=np.asarray(f,complex); g=np.asarray(g,complex); labels=np.asarray(labels); n=len(f)
    A=periodic_discrete_stft(f,g); D=periodic_discrete_stft(np.roll(f,shift)-f,g)
    if labels.shape!=A.shape or np.any(labels<0): raise ValueError('bad labels')
    mmax=int(labels.max())+1
    inc=np.zeros(mmax); here=np.zeros(mmax); shifted=np.zeros(mmax)
    for a in range(mmax):
        mask=(labels==a)
        inc[a]=float(np.sum(np.abs(D[mask])**2)/n)
        here[a]=float(np.sum(np.abs(A[mask])**2)/n)
        # sum |A[m-shift,k]|^2 over (m,k) in the current cell
        nbr=np.empty_like(A.real)
        for m in range(n): nbr[m]=np.abs(A[(m-shift)%n])**2
        shifted[a]=float(np.sum(nbr[mask])/n)
    return inc,here,shifted


def coherent_service_route(
    square_threshold: float,
    high_normalized_enstrophy: float,
    old_pool_capacity: float,
    old_old_service: float,
    cross_interface_service: float,
    new_new_edge_services: Sequence[float],
    ancestry_labels: Sequence[object] | None = None,
    dominant_fraction: float = 0.25,
    ancestry_alpha: float = 0.5,
) -> dict[str,float|str]:
    """Route an actual low-band coherent increment service measure.

    If d_high>=Y/4, dissipation wins.  Otherwise S_low>=Y/2.  Once old-pool
    capacity is <=Y/8, old-old service is too small.  Then either cross-interface
    service is >=Y/8 (Xi) or new-new service is >=Y/4.  The latter either has a
    dominant service edge, forcing coherent critical mass >=Y/32 at theta=1/4,
    or pays log2 ancestry entropy / 1/4 hidden same-ancestry pair mass.
    """
    vals=[square_threshold,high_normalized_enstrophy,old_pool_capacity,old_old_service,cross_interface_service]
    if min(vals)<0: raise ValueError('nonnegative inputs required')
    w=np.asarray(new_new_edge_services,float)
    if np.any(w<0): raise ValueError('nonnegative new-edge services required')
    if not (0<dominant_fraction<1 and 0<ancestry_alpha<1): raise ValueError('bad routing thresholds')
    Y=float(square_threshold)
    if high_normalized_enstrophy >= Y/4.0:
        return {'branch':'high_frequency_dissipation','threshold':Y/4.0,'branch_value':high_normalized_enstrophy,'margin':high_normalized_enstrophy-Y/4.0}
    slow=low_square_service_lower(Y,high_normalized_enstrophy)
    total=old_old_service+cross_interface_service+float(w.sum())
    if total+1e-12*max(1.0,Y)<slow:
        raise ValueError('cell service weights do not realize the forced low-band service')
    if old_pool_capacity > Y/8.0:
        return {'branch':'old_pool_not_yet_eroded','threshold':Y/8.0,'branch_value':old_pool_capacity,'low_service_lower':slow}
    if old_old_service > old_pool_capacity+1e-12*max(1.0,Y):
        raise ValueError('old-old service exceeds certified old-pool capacity')
    if cross_interface_service >= Y/8.0:
        return {'branch':'selected_interface_Xi','threshold':Y/8.0,'branch_value':cross_interface_service,'margin':cross_interface_service-Y/8.0}
    fresh=float(w.sum())
    if fresh+2e-12*max(1.0,Y) < Y/4.0:
        raise AssertionError('new-new service lower bound failed')
    if fresh<=0: raise AssertionError('positive new service required')
    p=w/fresh
    imax=int(np.argmax(p)); pmax=float(p[imax])
    if pmax>=dominant_fraction:
        edge=float(w[imax])
        # s_e <= 2 M(E_C+E_{C-r}); the RHS critical coherent cluster mass is >=s_e/2.
        mass=edge/2.0
        clean=dominant_fraction*Y/8.0
        if mass+2e-13*max(1.0,Y)<clean: raise AssertionError('dominant coherent cluster mass lower failed')
        return {'branch':'dominant_new_coherent_cluster','new_service':fresh,'dominant_edge_service':edge,'coherent_critical_mass_lower':mass,'clean_mass_lower':clean,'margin':mass-clean}
    q=float(np.dot(p,p)); h=-math.log(q); h0=-math.log(dominant_fraction)
    if h+1e-13<h0: raise AssertionError('service-edge collision entropy failed')
    if ancestry_labels is None:
        return {'branch':'new_service_collision_entropy','H_atomic':h,'entropy_lower':h0,'new_service':fresh}
    if len(ancestry_labels)!=len(w): raise ValueError('ancestry label length mismatch')
    chain=collision_chain(p,ancestry_labels)
    if chain['h_ancestry']>=ancestry_alpha*h0-1e-13:
        return {'branch':'new_service_Bellman_entropy','H_atomic':h,'H_ancestry':chain['h_ancestry'],'ancestry_entropy_lower':ancestry_alpha*h0,'new_service':fresh}
    pair_lower=dominant_fraction**ancestry_alpha-dominant_fraction
    if chain['hidden_pair_mass']+2e-13<pair_lower: raise AssertionError('new-service ancestry pair bound failed')
    return {'branch':'new_service_same_ancestry_cycle','H_atomic':h,'H_ancestry':chain['h_ancestry'],'hidden_pair_mass':chain['hidden_pair_mass'],'hidden_pair_lower':pair_lower,'new_service':fresh}


def exact_certificate() -> dict[str,str]:
    return {
        'filter_pigeonhole':'some r has ||delta_r u||_3^2 >= (Q/g1)^(2/3)',
        'actual_low_service':'S_low(r)=sum_{j<=0} M_j ||delta_r u_j||_2^2 >= Y-2 d_high',
        'coherent_edge_measure':'s_jC=M_j int_C |V_g delta_r u_j|^2, sum s_jC=S_low',
        'translation_covariance':'V_g(delta_r u)(X,k)=e^{-ik.r}V_g u(X-r,k)-V_g u(X,k)',
        'local_energy_capacity':'s_jC <= 2 M_j[E_j(C)+E_j(C-r)]',
        'clean_route':'d_high>=Y/4 OR oldcap>Y/8 OR Xi>=Y/8 OR fresh coherent mass>=Y/32 OR log2 entropy OR 1/4 cycle mass',
        'status':'EXACT_MOYAL_CELL_ROUTING_GIVEN_STANDARD_LP_BERNSTEIN',
    }


@dataclass(frozen=True)
class CoherentIncrementStress:
    samples:int
    worst_translation_covariance_residual:float
    worst_local_capacity_ratio:float
    minimum_route_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260808)->CoherentIncrementStress:
    rng=np.random.default_rng(seed); wt=wc=0.0; mm=float('inf'); branches={}
    # expensive STFT identities on a representative subset
    for _ in range(min(samples,3000)):
        n=int(rng.integers(8,36)); f=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g); r=int(rng.integers(-n//3,n//3+1))
        res=periodic_increment_covariance_residual(f,g,r); wt=max(wt,res)
        if res>2e-10*max(1.0,np.linalg.norm(f)): raise AssertionError('increment STFT covariance failed')
        A=periodic_discrete_stft(f,g); labels=rng.integers(0,int(rng.integers(2,10)),size=A.shape)
        inc,here,nbr=discrete_cell_increment_energies(f,g,r,labels)
        ratio=float(np.max(inc/np.maximum(2*(here+nbr),1e-300))); wc=max(wc,ratio)
        if np.any(inc>2*(here+nbr)+2e-11*max(1.0,float(inc.max()))): raise AssertionError('coherent increment local capacity failed')
    for _ in range(samples):
        Y=float(rng.lognormal(-1,.8)); d=float(rng.uniform(0,Y/3))
        if d>=Y/4:
            route=coherent_service_route(Y,d,0,0,0,[])
        else:
            oldcap=float(rng.uniform(0,Y/6))
            if oldcap>Y/8:
                slow=low_square_service_lower(Y,d); route=coherent_service_route(Y,d,oldcap,min(oldcap,slow),0,[max(0.0,slow-min(oldcap,slow))])
            else:
                old=float(rng.uniform(0,oldcap)); cross=float(rng.uniform(0,Y/6))
                slow=low_square_service_lower(Y,d); need=max(0.0,slow-old-cross)
                if cross>=Y/8:
                    w=[need]
                else:
                    n=int(rng.integers(3,20)); probs=rng.dirichlet(np.ones(n)); w=(probs*max(need,Y/4)).tolist()
                labels=rng.integers(0,max(2,len(w)//4+1),size=len(w)).tolist()
                route=coherent_service_route(Y,d,oldcap,old,cross,w,labels)
        b=str(route['branch']); branches[b]=branches.get(b,0)+1
        if 'margin' in route: mm=min(mm,float(route['margin']))
        else: mm=min(mm,0.0)
    return CoherentIncrementStress(samples,wt,wc,mm,branches)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-coherent-increment-service'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=exact_certificate(); out=stress(args.samples)
    (args.outdir/'coherent_increment_service.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f'''# Coherent increment service: SGS charge becomes phase-space edges

Status: **{cert['status']}**.

Let `Q` be the cubic increment charge and `g1=||G||_1`.  Some filter displacement `r` satisfies

`||delta_r u||_3^2 >= (Q/g1)^(2/3)`.

With `Y=(Q/g1)^(2/3)/(C_LP C_B)^2`, standard LP/Bernstein gives

`S_low(r)=sum_(j<=0) M_j ||delta_r u_j||_2^2 >= Y-2 d_high`.

For any normalized coherent window and phase-space partition,

`s_(j,C)=M_j int_C |V_g delta_r u_j|^2 dmu`

is a positive **actual increment-service measure** and sums exactly to `S_low`. Translation covariance gives

`V_g(delta_r u)(X,k)=exp(-ik.r)V_g u(X-r,k)-V_g u(X,k)`,

hence

`s_(j,C) <= 2 M_j[E_j(C)+E_j(C-r)]`.

On `d_high<Y/4`, `S_low>=Y/2`.  Once whole-old-pool erosion gives `old_capacity<=Y/8`, either selected old/new interface service is at least `Y/8`, or new/new coherent service is at least `Y/4`.  With the clean quarter-dominance threshold, the latter yields

- a dominant new coherent cluster with critical Moyal mass at least `Y/32`; or
- ancestry Bellman entropy at least `log 2`; or
- same-ancestry hidden pair/cycle mass at least `1/4`.

Thus the old low-band reservoir selected by the Onsager collision is no longer merely an aggregate global mass: the **actual increment itself** generates nearby coherent-cell edges, and after old-pool erosion those edges must cross an interface or create new coherent ancestry.

Stress: `{out.samples}`
- worst STFT translation-covariance residual: `{out.worst_translation_covariance_residual:.3e}`
- worst local increment-capacity ratio: `{out.worst_local_capacity_ratio:.9f}`
- minimum routing margin: `{out.minimum_route_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
