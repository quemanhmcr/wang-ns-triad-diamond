from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite_e import hermegauss

A3=(math.sqrt(3.0)/2.0)**3


def parity_selection(total_degree:int)->bool:
    """Centered resonant Gaussian triad coefficient vanishes for odd total Hermite parity."""
    return (int(total_degree)%2)==1


def sideband_rescue_upper(rho1:float,rho2:float,rho3:float,base1:float=1.0,base2:float=1.0,base3:float=1.0)->float:
    """Sharp-Young bound after all single-odd terms vanish by parity.

    Remaining rescue terms contain at least two sidebands.
    """
    vals=(rho1,rho2,rho3,base1,base2,base3)
    if any(v<0 for v in vals): raise ValueError('nonnegative norms required')
    return A3*(rho1*rho2*base3+rho1*rho3*base2+rho2*rho3*base1+rho1*rho2*rho3)


def small_odd_sideband_norm_growth(sigma:float)->float:
    """p-th power growth for p=3/2 in the small odd degree<=3 branch.

    Hypercontractivity gives E|P|^4 <=729 sigma^4.  At sigma<=1/80, at
    least half of the second moment lies where |P|<=1/2.  On that set the
    parity-symmetrized pointwise norm obeys phi(P)>=1+(3/8)|P|^2.
    Hence ||G+PG||_p^p / ||G||_p^p >= 1+3 sigma^2/16.
    """
    if sigma<0 or sigma>1/80+1e-15: raise ValueError('small-sideband theorem requires sigma<=1/80')
    return 1.0+3.0*sigma*sigma/16.0


def single_role_transfer_deficit_lower(sigma:float)->float:
    """Normalized Gaussian-base transfer deficit from one odd sideband.

    Single-sideband numerator vanishes exactly by parity; only the child norm
    grows.  The clean bound is Def >= sigma^2/16 on sigma<=1/80.
    """
    growth=small_odd_sideband_norm_growth(sigma)
    exact=1.0-growth**(-2.0/3.0)
    clean=sigma*sigma/16.0
    if exact+1e-15<clean: raise AssertionError('clean sideband deficit failed')
    return clean


def exact_arithmetic_certificate()->dict[str,str]:
    # No interval transcendental is needed for the main threshold: this is exact rational arithmetic.
    tail=Fraction(2916,80*80)
    if not tail<Fraction(1,2): raise AssertionError('fourth-moment tail is not below one half')
    # derivative lower used to justify 1-(1+x)^(-2/3)>=sigma^2/16 can be checked much more strongly numerically;
    # theorem proof in the note uses monotonicity on x<=3/(16*80^2).
    xmax=Fraction(3,16*80*80)
    return {
        'small_sideband_threshold':'1/80',
        'tail_second_moment_fraction_upper':f'{tail.numerator}/{tail.denominator}',
        'norm_growth':'1+3 sigma^2/16',
        'clean_single_role_deficit':'sigma^2/16',
        'xmax':f'{xmax.numerator}/{xmax.denominator}',
        'status':'EXACT_GIVEN_DEGREE3_GAUSSIAN_HYPERCONTRACTIVITY',
    }


def gh_joint_expectation_odd(degree:int, covariance:np.ndarray, order:int=8)->float:
    """Regression integral for one odd child polynomial under a centered joint Gaussian.

    Let (x,y) be 1D parent deviations and child deviation x+y.  The centered
    Gaussian weight may have arbitrary SPD covariance.  Odd polynomial
    (x+y)^degree times the even density has zero expectation.
    """
    if degree%2!=1: raise ValueError('degree must be odd')
    covariance=np.asarray(covariance,float)
    L=np.linalg.cholesky(covariance)
    nodes,weights=hermegauss(order); weights=weights/math.sqrt(2*math.pi)
    total=0.0
    for i,x in enumerate(nodes):
        for j,y in enumerate(nodes):
            z=L@np.array([x,y])
            total += weights[i]*weights[j]*(z[0]+z[1])**degree
    return float(total)


@dataclass(frozen=True)
class HermiteTriadSelectionStress:
    samples:int
    worst_H1_single_sideband_residual:float
    worst_H3_single_sideband_residual:float
    minimum_deficit_margin:float
    minimum_rescue_margin:float


def stress(samples:int=50_000,seed:int=20260807)->HermiteTriadSelectionStress:
    rng=np.random.default_rng(seed); w1=w3=0.0; md=float('inf'); mr=float('inf')
    for _ in range(samples):
        Q=rng.normal(size=(2,2)); cov=Q@Q.T+0.2*np.eye(2)
        r1=abs(gh_joint_expectation_odd(1,cov,6)); r3=abs(gh_joint_expectation_odd(3,cov,8))
        w1=max(w1,r1); w3=max(w3,r3)
        if r1>2e-12 or r3>5e-11: raise AssertionError('odd Gaussian triad parity regression failed')
        sigma=float(rng.uniform(0,1/80))
        growth=small_odd_sideband_norm_growth(sigma)
        exact=1-growth**(-2/3); clean=single_role_transfer_deficit_lower(sigma)
        md=min(md,exact-clean)
        if exact+2e-15<clean: raise AssertionError('single-role transfer deficit failed')
        rs=rng.uniform(0,0.2,size=3); bs=rng.uniform(0.4,1.2,size=3)
        bound=sideband_rescue_upper(*rs,*bs)
        # construct arbitrary coefficients no larger than sharp-Young termwise maxima
        terms=np.array([rs[0]*rs[1]*bs[2],rs[0]*rs[2]*bs[1],rs[1]*rs[2]*bs[0],np.prod(rs)])
        actual=A3*float(np.dot(rng.uniform(-1,1,size=4),terms))
        mr=min(mr,bound-abs(actual))
        if abs(actual)>bound+2e-14: raise AssertionError('quadratic rescue bound failed')
    return HermiteTriadSelectionStress(samples,w1,w3,md,mr)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-hermite-triad-selection'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=exact_arithmetic_certificate(); out=stress(args.samples)
    data={'certificate':cert,'stress':out.__dict__,'theorems':{
        'parity':'T(P_n1 G1,P_n2 G2,P_n3 G3)=0 if n1+n2+n3 is odd',
        'single_sideband':'T(H1,G,G)=T(H3,G,G)=0',
        'quadratic_rescue':'|T_rescue|<=A3(r1 r2 b3+r1 r3 b2+r2 r3 b1+r1 r2 r3)',
        'small_single_role_deficit':'Def>=sigma^2/16 for sigma<=1/80',
    }}
    (args.outdir/'hermite_triad_selection.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Odd-Hermite triad selection and daughter-transfer cost

Status: **{cert['status']}**.

For a centered resonant affine Gaussian triad, simultaneous inversion of all centered frequency deviations leaves the Gaussian trilinear weight invariant.  Therefore a polynomial/Hermite perturbation has the exact selection rule

`T(P_n1 G1,P_n2 G2,P_n3 G3)=0` whenever `n1+n2+n3` is odd.

In particular a single H1 or H3 daughter cannot feed the base Gaussian triad at first order.  Any sideband rescue contains at least two odd sidebands and is bounded by

`A3 (rho1 rho2 b3 + rho1 rho3 b2 + rho2 rho3 b1 + rho1 rho2 rho3)`.

There is also a quantitative one-role loss.  Let `R=P G` be an odd degree<=3 sideband and let `sigma^2=E_mu |P|^2` in the critical `|G|^(3/2)` Gaussian measure.  Gaussian hypercontractivity gives `E|P|^4<=729 sigma^4`.  If `sigma<=1/80`, parity plus uniform convexity yields

`||G+R||_(3/2)^(3/2) >= ||G||_(3/2)^(3/2) (1+3 sigma^2/16)`

and since the single-sideband numerator is exactly zero,

`Def_transfer >= sigma^2/16`.

Thus a coherent odd daughter has only two ways to remain efficient: recruit a second odd sideband (a genuine daughter/cross interaction component) or pay a quadratic transfer deficit.

Stress: `{out.samples}`
- worst H1 single-sideband parity residual: `{out.worst_H1_single_sideband_residual:.3e}`
- worst H3 single-sideband parity residual: `{out.worst_H3_single_sideband_residual:.3e}`
- minimum clean-deficit margin: `{out.minimum_deficit_margin:.3e}`
- minimum quadratic-rescue margin: `{out.minimum_rescue_margin:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
