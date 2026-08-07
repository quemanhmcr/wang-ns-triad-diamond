from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def cubic_boundary_pointwise_bound(U: np.ndarray, W: np.ndarray, R: np.ndarray) -> tuple[float,float]:
    """Bound |e W|+|R U| by critical cubic densities.

    e=|U|^2/2.  Frobenius norm is used for R.
    Young gives
      .5 |U|^2 |W| + |R||U|
      <= (2/3)|U|^3 +(1/6)|W|^3 +(2/3)|R|^(3/2).
    """
    U=np.asarray(U,float); W=np.asarray(W,float); R=np.asarray(R,float)
    u=float(np.linalg.norm(U)); w=float(np.linalg.norm(W)); r=float(np.linalg.norm(R,'fro'))
    lhs=0.5*u*u*w+r*u
    rhs=(2.0/3.0)*u**3+(1.0/6.0)*w**3+(2.0/3.0)*r**1.5
    return lhs,rhs


def cubic_annular_charge_lower_bound(leakage_abs: float, grad_chi_sup: float) -> float:
    """Scale-critical charge forced by resolved transport+RU boundary leakage."""
    if leakage_abs < 0 or grad_chi_sup <= 0:
        raise ValueError("require nonnegative leakage and positive gradient")
    return leakage_abs / grad_chi_sup


def affine_gradient_upper(N: float, Cchi: float, M: float) -> float:
    if min(N,Cchi,M)<=0: raise ValueError("positive inputs")
    return 1.5*N*Cchi/M


def affine_cubic_charge_lower_bound(leakage_abs: float, N: float, Cchi: float, M: float) -> float:
    """Use ||grad chi||<=3 N Cchi/(2M)."""
    return 2.0*M*leakage_abs/(3.0*N*Cchi)


def affine_pressure_charge_lower_bound(raw_sgs: float, N: float, Cchi: float, M: float) -> float:
    """Pressure-cancellation branch: charge >= S/(2||grad chi||) >= S M/(3 N Cchi)."""
    if raw_sgs < 0: raise ValueError("nonnegative raw SGS required")
    return raw_sgs*M/(3.0*N*Cchi)


def viscous_boundary_lifetime_bound(c: float, nu: float, delta: float, Cchi: float, M: float, energy_sup: float) -> float:
    """Time-integrated resolved viscous boundary flux on T=c N^-2.

    If supp Uhat lies in |xi|<=exp(delta)N and
    ||grad chi||<=3NCchi/(2M), then
      nu int_0^T |int grad chi.grad e|
      <= (3/2) exp(delta) c nu Cchi M^-1 sup_t ||U||_2^2.
    """
    if min(c,nu,Cchi,M,energy_sup)<0 or M<=0: raise ValueError("invalid inputs")
    return 1.5*math.exp(delta)*c*nu*Cchi*energy_sup/M


def clean_viscous_boundary_lifetime_bound(c: float, nu: float, Cchi: float, M: float, energy_sup: float) -> float:
    """For delta<=1/20, exp(delta)<11/10: coefficient 33/20."""
    return (33.0/20.0)*c*nu*Cchi*energy_sup/M


def partition_flux_residual(etas_dt: np.ndarray, grad_etas: np.ndarray, energy_density: np.ndarray, flux: np.ndarray) -> float:
    """Pointwise sum of partition boundary terms.

    If sum eta_alpha=1 at all times, then sum dt eta=0 and sum grad eta=0,
    so sum_alpha [e dt eta_alpha + grad eta_alpha . F]=0 exactly.
    Inputs grad_etas shape (n,d), flux shape (d,), energy_density scalar/length-1.
    """
    e=float(np.asarray(energy_density).reshape(-1)[0])
    return e*float(np.sum(etas_dt))+float(np.sum(np.asarray(grad_etas,float),axis=0)@np.asarray(flux,float))


def arb_boundary_certificate() -> dict[str,str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec=160
    if not ((arb(1)/20).exp() < arb(11)/10):
        raise AssertionError("exp(1/20)<11/10 failed")
    if Fraction(3,2)*Fraction(11,10) != Fraction(33,20):
        raise AssertionError("exact rational viscous coefficient arithmetic failed")
    return {"smooth_filter_high_side":"delta<=1/20", "exp_delta_bound":"exp(delta)<11/10", "clean_viscous_coefficient":"33/20", "status":"CERTIFIED"}


@dataclass(frozen=True)
class BoundaryStress:
    samples: int
    worst_cubic_ratio: float
    minimum_cubic_margin: float
    worst_partition_residual: float
    minimum_clean_viscous_margin: float


def stress(samples:int=50_000,seed:int=20260807)->BoundaryStress:
    rng=np.random.default_rng(seed)
    wr=wp=0.0; mm=float('inf'); vm=float('inf')
    for _ in range(samples):
        U=rng.normal(size=3); W=rng.normal(size=3); R=rng.normal(size=(3,3))
        lhs,rhs=cubic_boundary_pointwise_bound(U,W,R)
        if lhs>rhs+2e-12: raise AssertionError("cubic boundary Young bound failed")
        mm=min(mm,rhs-lhs)
        if rhs>1e-14: wr=max(wr,lhs/rhs)
        n=int(rng.integers(2,9)); d=3
        dt=rng.normal(size=n); dt-=np.mean(dt)
        G=rng.normal(size=(n,d)); G-=np.mean(G,axis=0)
        e=np.array([abs(float(rng.normal()))]); F=rng.normal(size=d)
        pr=abs(partition_flux_residual(dt,G,e,F)); wp=max(wp,pr)
        if pr>2e-12: raise AssertionError("partition boundary cancellation failed")
        c=float(rng.uniform(.01,1.)); nu=float(rng.uniform(0.,2.)); C=float(rng.uniform(.1,3.)); M=float(rng.uniform(1.,40.)); E=float(rng.uniform(.01,10.)); delta=float(rng.uniform(0.,.05))
        exact=viscous_boundary_lifetime_bound(c,nu,delta,C,M,E)
        clean=clean_viscous_boundary_lifetime_bound(c,nu,C,M,E)
        vm=min(vm,clean-exact)
        if exact>clean+2e-12: raise AssertionError("clean viscous boundary bound failed")
    return BoundaryStress(samples,wr,mm,wp,vm)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-affine-sgs-boundary'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_boundary_certificate(); out=stress(args.samples)
    (args.outdir/'affine_sgs_boundary_ledger.json').write_text(json.dumps({'certificate':cert,'stress':out.__dict__},indent=2),encoding='utf-8')
    md=f"""# Affine SGS boundary ledger

Status: **{cert['status']}**.

For an affine window transported by the strict low-pass affine jet, write
`W=U-V_aff`.  The nonviscous combined-work leakage is
`int grad chi . (e W + R U)`.  Pointwise Young gives

`.5 |U|^2 |W| + |R||U| <= (2/3)|U|^3 +(1/6)|W|^3 +(2/3)|R|^(3/2)`.

Hence large differential-advection/SGS-transport leakage forces the scale-critical annular charge
`int_A(|U|^3+|W|^3+|R|^(3/2)) >= |L_cubic|/||grad chi||`.
For the affine shell window this is at least `2M |L_cubic|/(3 N Cchi)`.
Pressure cancellation separately forces `S M/(3 N Cchi)` in `|U|^3+|P|^(3/2)`.

The resolved viscous boundary term is not a new source: on a parabolic lifetime and the smooth spectral support,
`|L_nu| <= (33/20)c nu Cchi M^-1 sup ||U||_2^2` for `delta<=1/20`.
Thus it renormalizes the existing `1/M` localization coefficient.

A quadratic spatial partition `sum eta_alpha=1` has exact total boundary cancellation:
`sum_alpha(e partial_t eta_alpha + grad eta_alpha.F)=0`.  Overlap can matter after selecting a lineage, but there is no global packet-count loss from the partition itself.

Stress checks: `{out.samples}`
- worst cubic Young ratio: `{out.worst_cubic_ratio:.9f}`
- minimum cubic margin: `{out.minimum_cubic_margin:.3e}`
- worst partition cancellation residual: `{out.worst_partition_residual:.3e}`
- minimum clean viscous-bound margin: `{out.minimum_clean_viscous_margin:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
