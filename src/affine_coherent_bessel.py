from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

CLEAN_SEPARATION = 4.0
CLEAN_BESSEL = 25.0 / 4.0
RIESZ_SEPARATION = 5.0
RIESZ_OFFDIAGONAL = 3.0 / 50.0


def intrinsic_phase_point(X: np.ndarray, k: np.ndarray, L: np.ndarray) -> np.ndarray:
    X=np.asarray(X,float); k=np.asarray(k,float); L=np.asarray(L,float)
    if X.shape!=(3,) or k.shape!=(3,) or L.shape!=(3,3): raise ValueError('bad phase point')
    return np.concatenate([0.5*np.linalg.solve(L,X), L.T@k])


def coherent_overlap_magnitude(X:np.ndarray,k:np.ndarray,Y:np.ndarray,l:np.ndarray,L:np.ndarray)->float:
    """Exact L2 overlap magnitude of equal-covariance affine Gaussian packets.

    g_(X,k,L)(x)=C det(L)^(-1/2) exp(-|L^-1(x-X)|^2/4) exp(i k.x).
    """
    dx=np.linalg.solve(np.asarray(L,float),np.asarray(X,float)-np.asarray(Y,float))
    dq=np.asarray(L,float).T@(np.asarray(k,float)-np.asarray(l,float))
    return math.exp(-0.125*float(dx@dx)-0.5*float(dq@dq))


def overlap_from_intrinsic_points(zeta:np.ndarray,eta:np.ndarray)->float:
    d=np.asarray(zeta,float)-np.asarray(eta,float)
    return math.exp(-0.5*float(d@d))


def shell_packing_count_upper(n:int,dimension:int=6)->int:
    """Coarse delta-separated packing bound for shell n delta <= d < (n+1)delta."""
    if n<1 or dimension<1: raise ValueError('bad shell')
    return (2*n+3)**dimension


def gram_row_sum_bound(delta:float=CLEAN_SEPARATION,dimension:int=6,terms:int=50)->float:
    if delta<=0 or dimension<1 or terms<2: raise ValueError('bad data')
    return 1.0+sum(shell_packing_count_upper(n,dimension)*math.exp(-0.5*(n*delta)**2) for n in range(1,terms))


def gram_offdiagonal_row_bound(delta:float=RIESZ_SEPARATION,dimension:int=6,terms:int=50)->float:
    return gram_row_sum_bound(delta,dimension,terms)-1.0


def synthesis_coefficient_energy_upper(field_energy:float)->float:
    """For a 5-separated equal-covariance coherent synthesis: sum|c|^2 <=50/47 ||f||^2."""
    if field_energy<0: raise ValueError('nonnegative energy required')
    return (50.0/47.0)*field_energy


def arb_clean_bessel_certificate()->dict[str,str]:
    """Certify the infinite 6D packing row sum is <25/4 at delta=4.

    t_n=(2n+3)^6 exp(-8n^2). For n>=2 the ratio decreases; certify r_2<1/1000,
    then sum_{n>=2}t_n<=t_2/(1-r_2).
    """
    try:
        from flint import arb,ctx
    except ImportError as exc: raise RuntimeError('python-flint required') from exc
    ctx.prec=160
    t1=arb(5)**6*(-arb(8)).exp()
    t2=arb(7)**6*(-arb(32)).exp()
    r2=(arb(9)/7)**6*(-arb(40)).exp()
    if not(r2<arb(1)/1000): raise AssertionError(f'tail ratio failed {r2}')
    upper=arb(1)+t1+t2/(1-r2)
    if not(upper<arb(25)/4): raise AssertionError(f'Bessel row bound failed {upper}')
    # At delta=5, certify the off-diagonal Gram row is <3/50.
    u1=arb(5)**6*(-arb(25)/2).exp()
    u2=arb(7)**6*(-arb(50)).exp()
    rr=(arb(9)/7)**6*(-arb(125)/2).exp()  # t3/t2, tail ratio decreases thereafter
    if not(rr<arb(1)/1000000): raise AssertionError(f'Riesz tail ratio failed {rr}')
    off=u1+u2/(1-rr)
    if not(off<arb(3)/50): raise AssertionError(f'Riesz off-diagonal bound failed {off}')
    return {
        'intrinsic_phase_coordinate':'zeta=(L^-1 X/2,L^T k)',
        'overlap':'|<g_a,g_b>|=exp(-|zeta_a-zeta_b|^2/2)',
        'separation':'delta>=4',
        'packing_row_sum_ball':str(upper),
        'clean_Bessel_constant':'25/4',
        'Riesz_separation':'delta>=5',
        'Riesz_offdiagonal_ball':str(off),
        'Riesz_Gram_spectrum':'[47/50,53/50]',
        'synthesis_coefficient_budget':'sum |c_a|^2 <= (50/47)||sum c_a g_a||_2^2',
        'status':'CERTIFIED_EQUAL_COVARIANCE_COHERENT_BESSEL_RIESZ',
    }


def gram_matrix(points:np.ndarray)->np.ndarray:
    Z=np.asarray(points,float)
    D=Z[:,None,:]-Z[None,:,:]
    return np.exp(-0.5*np.sum(D*D,axis=2))


def min_pair_distance(points:np.ndarray)->float:
    Z=np.asarray(points,float)
    if len(Z)<2:return math.inf
    D=Z[:,None,:]-Z[None,:,:]
    d2=np.sum(D*D,axis=2); d2+=np.eye(len(Z))*1e100
    return math.sqrt(float(np.min(d2)))


def random_separated_points(rng:np.random.Generator,n:int,delta:float=4.0,dimension:int=6)->np.ndarray:
    pts=[]; attempts=0
    radius=max(6.0,2.5*n**(1/dimension)*delta)
    while len(pts)<n and attempts<200000:
        z=rng.uniform(-radius,radius,size=dimension); attempts+=1
        if all(np.linalg.norm(z-p)>=delta for p in pts): pts.append(z)
    if len(pts)<n: raise RuntimeError('could not generate separated set')
    return np.array(pts)


@dataclass(frozen=True)
class CoherentBesselStress:
    samples:int
    worst_affine_invariance_residual:float
    worst_overlap_coordinate_residual:float
    maximum_gram_operator_norm:float
    maximum_gram_row_sum:float
    minimum_separation_margin:float
    minimum_Riesz_eigenvalue:float
    maximum_Riesz_eigenvalue:float


def stress(samples:int=50_000,seed:int=20260808)->CoherentBesselStress:
    rng=np.random.default_rng(seed); wa=wo=0.; mg=mr=0.; ms=float('inf'); rlo=float('inf'); rhi=0.
    # 50k arbitrary affine/pair identity checks.
    for _ in range(samples):
        L=rng.normal(size=(3,3))
        while abs(np.linalg.det(L))<.1: L=rng.normal(size=(3,3))
        X=rng.normal(size=3);Y=rng.normal(size=3);k=rng.normal(size=3);l=rng.normal(size=3)
        z=intrinsic_phase_point(X,k,L); e=intrinsic_phase_point(Y,l,L)
        a=coherent_overlap_magnitude(X,k,Y,l,L); b=overlap_from_intrinsic_points(z,e)
        wo=max(wo,abs(a-b))
        S=rng.normal(size=(3,3))
        while abs(np.linalg.det(S))<.1: S=rng.normal(size=(3,3))
        Lp=S@L;Xp=S@X;Yp=S@Y;kp=np.linalg.solve(S.T,k);lp=np.linalg.solve(S.T,l)
        zp=intrinsic_phase_point(Xp,kp,Lp); ep=intrinsic_phase_point(Yp,lp,Lp)
        wa=max(wa,float(np.linalg.norm(zp-z)),float(np.linalg.norm(ep-e)))
    # Adversarial finite Gram probes.
    for _ in range(1500):
        n=int(rng.integers(2,28)); Z=random_separated_points(rng,n)
        sep=min_pair_distance(Z); ms=min(ms,sep-CLEAN_SEPARATION)
        G=gram_matrix(Z); rows=float(np.max(np.sum(np.abs(G),axis=1)))
        op=float(np.linalg.eigvalsh(G)[-1]); mr=max(mr,rows);mg=max(mg,op)
        if op>CLEAN_BESSEL+2e-10 or rows>CLEAN_BESSEL+2e-10: raise AssertionError('finite Gram exceeded clean Schur bound')
    # Riesz probes at separation 5.
    for _ in range(1500):
        n=int(rng.integers(2,28)); Z=random_separated_points(rng,n,RIESZ_SEPARATION)
        ev=np.linalg.eigvalsh(gram_matrix(Z)); rlo=min(rlo,float(ev[0]));rhi=max(rhi,float(ev[-1]))
        if ev[0]<47/50-2e-10 or ev[-1]>53/50+2e-10: raise AssertionError('finite Gram violated clean Riesz spectrum')
    if wa>2e-11 or wo>2e-13: raise AssertionError('coherent affine identity failed')
    return CoherentBesselStress(samples,wa,wo,mg,mr,ms,rlo,rhi)


def main()->None:
    ap=argparse.ArgumentParser();ap.add_argument('--samples',type=int,default=50_000);ap.add_argument('--outdir',type=Path,default=Path('results-affine-coherent-bessel'));args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_clean_bessel_certificate();out=stress(args.samples)
    (args.outdir/'affine_coherent_bessel.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f"""# Affine coherent-state Bessel budget

Status: **{cert['status']}**.

For equal physical covariance factor `L`, normalized affine Gaussian packets satisfy exactly

`|<g_(X,k,L),g_(Y,l,L)>| = exp[-|L^-1(X-Y)|^2/8-|L^T(k-l)|^2/2]`.

With intrinsic phase point `zeta=(L^-1 X/2,L^T k)`, this is simply `exp(-|zeta_a-zeta_b|^2/2)` and is invariant under a common physical affine change.

If the intrinsic phase points are `4`-separated in R^6, disjoint-ball packing gives at most `(2n+3)^6` centers in shell `[4n,4(n+1))`. Arb certifies the infinite absolute Gram row sum is `<25/4`. Schur therefore gives

`sum_a |<f,g_a>|^2 <= (25/4)||f||_2^2`.

At separation `5`, Arb further certifies the off-diagonal Gram row `<3/50`, hence the Gram spectrum lies in `[47/50,53/50]`. Therefore arbitrary synthesis coefficients in that separated equal-covariance family satisfy `sum|c_a|^2 <= (50/47)||sum c_a g_a||_2^2`. This directly supplies the Bessel/frame coefficient budget needed by old-pool erosion inside one covariance cell. Changes of covariance cell and transfer extraction before the separated coherent synthesis remain iterative-interface issues.

Stress: `{out.samples}` affine/pair checks plus 1500 finite Gram probes
- worst affine-coordinate invariance residual: `{out.worst_affine_invariance_residual:.3e}`
- worst overlap-coordinate residual: `{out.worst_overlap_coordinate_residual:.3e}`
- maximum sampled Gram operator norm: `{out.maximum_gram_operator_norm:.9f}`
- maximum sampled absolute row sum: `{out.maximum_gram_row_sum:.9f}`
- minimum separation margin: `{out.minimum_separation_margin:.3e}`
- minimum sampled 5-separated Gram eigenvalue: `{out.minimum_Riesz_eigenvalue:.9f}`
- maximum sampled 5-separated Gram eigenvalue: `{out.maximum_Riesz_eigenvalue:.9f}`
"""
    (args.outdir/'summary.md').write_text(md);print(md)
if __name__=='__main__':main()
