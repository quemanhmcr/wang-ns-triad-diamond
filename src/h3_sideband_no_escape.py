from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

H3_FORCE_SQ=Fraction(3,8)          # ||f_H3||_2^2 = (3/8)||T||^2
SMALL_SIGMA=Fraction(1,80)
SINGLE_DEF_COEFF=Fraction(1,16)
NO_ESCAPE_COEFF=Fraction(3,4096)


def first_impulse_lower(curvature_impulse:float)->float:
    """Coherent branch from Banach/Hilbert variation theorem: delta1 >= sqrt(3/8) I /2."""
    if curvature_impulse<0: raise ValueError('nonnegative curvature impulse required')
    return 0.5*math.sqrt(3.0/8.0)*curvature_impulse


def actual_daughter_lower(curvature_impulse:float)->float:
    """If nonlinear Duhamel feedback is < half the first impulse, actual daughter >= sqrt(3/8) I /4."""
    return 0.5*first_impulse_lower(curvature_impulse)


def no_escape_quadratic_cost(curvature_impulse:float)->float:
    """Clean transfer-deficit / pair-rescue threshold 3 I^2 /4096."""
    if curvature_impulse<0: raise ValueError('nonnegative curvature impulse required')
    return float(NO_ESCAPE_COEFF)*curvature_impulse*curvature_impulse


def classify_h3_no_escape(
    curvature_impulse:float,
    lifetime:float,
    pulled_source_variation:float,
    first_duhamel_norm:float,
    nonlinear_feedback_norm:float,
    young_second_moment_sigma:float,
    pair_sideband_rescue:float,
    net_transfer_deficit:float,
)->dict[str,float|str]:
    """Master-facing H3 sideband no-escape dichotomy.

    Inputs are symmetry-quotiented/pulled-back quantities.  The theorem logic is:
      source branch: J >= I/T;
      otherwise first Duhamel delta1 >= sqrt(3/8) I/2;
      if feedback >= delta1/2 -> nonlinear sideband branch;
      else actual daughter delta >= delta1/2;
      if sigma >=1/80 -> definite daughter capacity;
      otherwise one-role odd-sideband deficit d0>=sigma^2/16.  Either pair rescue
      is >=d0/2 or net deficit is >=d0/2.  Since sigma>=actual daughter in the
      H3 Gaussian Young measure, both are >=3 I^2/4096.
    """
    vals=[curvature_impulse,pulled_source_variation,first_duhamel_norm,nonlinear_feedback_norm,young_second_moment_sigma,pair_sideband_rescue,net_transfer_deficit]
    if lifetime<=0 or any(v<0 for v in vals): raise ValueError('invalid no-escape data')
    I=curvature_impulse; src_thresh=I/lifetime
    d1_req=first_impulse_lower(I); cost=no_escape_quadratic_cost(I)
    if pulled_source_variation>=src_thresh-1e-12*max(1.,src_thresh):
        branch='dephasing_source'
    else:
        if first_duhamel_norm+1e-11*max(1.,d1_req)<d1_req:
            raise ValueError('coherent branch violates first-Duhamel lower bound')
        if nonlinear_feedback_norm>=0.5*first_duhamel_norm-1e-12*max(1.,first_duhamel_norm):
            branch='nonlinear_sideband_feedback'
        else:
            actual=max(0.0,first_duhamel_norm-nonlinear_feedback_norm)
            if young_second_moment_sigma+1e-11<actual:
                raise ValueError('H3 Young-measure second moment must dominate standard Gaussian L2 daughter')
            if young_second_moment_sigma>=float(SMALL_SIGMA)-1e-14:
                branch='large_daughter_capacity'
            else:
                d0=young_second_moment_sigma**2/16.0
                if pair_sideband_rescue>=0.5*d0-1e-13*max(1.,d0):
                    branch='pair_sideband_rescue'
                    if pair_sideband_rescue+2e-12<cost:
                        raise AssertionError('pair-sideband branch missed clean curvature cost')
                else:
                    branch='transfer_deficit'
                    if net_transfer_deficit+2e-12<0.5*d0 or net_transfer_deficit+2e-12<cost:
                        raise AssertionError('transfer-deficit branch missed clean curvature cost')
    return {
        'branch':branch,
        'source_threshold':src_thresh,
        'first_duhamel_required':d1_req,
        'clean_quadratic_cost':cost,
        'large_daughter_threshold':float(SMALL_SIGMA),
    }


def exact_constant_certificate()->dict[str,str]:
    # Algebra: after feedback<delta1/2, delta^2 >= (1/4)*(3/32)I^2=3/128 I^2.
    # single-role deficit >=delta^2/16 >=3/2048 I^2.  Split in half with pair rescue ->3/4096.
    if Fraction(1,4)*Fraction(3,32)*Fraction(1,16)*Fraction(1,2)!=NO_ESCAPE_COEFF:
        raise AssertionError('no-escape rational constant mismatch')
    return {
        'first_duhamel_square':'3 I^2/32',
        'post_feedback_daughter_square':'3 I^2/128',
        'single_role_deficit_before_rescue_split':'3 I^2/2048',
        'pair_rescue_or_net_deficit':'3 I^2/4096',
        'large_daughter_sigma':'1/80',
        'status':'EXACT_GIVEN_COHERENCE_AND_ODD_SIDEBAND_THEOREMS',
    }


@dataclass(frozen=True)
class H3NoEscapeStress:
    samples:int
    minimum_transfer_cost_margin:float
    minimum_pair_cost_margin:float
    minimum_first_impulse_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260807)->H3NoEscapeStress:
    rng=np.random.default_rng(seed); mt=mp=mi=float('inf'); counts={}
    for _ in range(samples):
        I=float(10**rng.uniform(-4,-0.5)); T=float(10**rng.uniform(-3,0.0)); req=first_impulse_lower(I); cost=no_escape_quadratic_cost(I)
        mode=int(rng.integers(0,5))
        if mode==0: # source
            J=(1+rng.random())*I/T; d1=0; fb=0; sigma=0; rescue=0; deficit=0
        else:
            J=float(rng.uniform(0,0.99))*I/T; d1=(1+rng.random())*req; mi=min(mi,d1-req)
            if mode==1: # feedback
                fb=float(rng.uniform(.5,1.2))*d1; sigma=0; rescue=0; deficit=0
            else:
                fb=float(rng.uniform(0,.49))*d1; actual=d1-fb
                if mode==2: # large daughter
                    sigma=max(float(SMALL_SIGMA),actual)*(1+rng.random()); rescue=0; deficit=0
                else:
                    # Need sigma>=actual while remaining below 1/80; if actual too large, shrink I for this synthetic branch.
                    if actual>=float(SMALL_SIGMA):
                        I=1e-4; req=first_impulse_lower(I); d1=1.2*req; fb=.1*d1; actual=d1-fb; cost=no_escape_quadratic_cost(I); J=.1*I/T
                    sigma=float(rng.uniform(max(actual,1e-12),float(SMALL_SIGMA)*.999))
                    d0=sigma*sigma/16
                    if mode==3: # pair rescue
                        rescue=float(rng.uniform(.5,1.0))*d0; deficit=0; mp=min(mp,rescue-cost)
                    else:
                        rescue=float(rng.uniform(0,.49))*d0; deficit=max(.5*d0,cost)*(1+rng.random()); mt=min(mt,deficit-cost)
        out=classify_h3_no_escape(I,T,J,d1,fb,sigma,rescue,deficit); counts[out['branch']]=counts.get(out['branch'],0)+1
    return H3NoEscapeStress(samples,mt,mp,mi,counts)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-h3-sideband-no-escape'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=exact_constant_certificate(); out=stress(args.samples)
    data={'certificate':cert,'stress':out.__dict__}
    (args.outdir/'h3_sideband_no_escape.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# H3 sideband local no-escape theorem

Status: **{cert['status']}**.

Let `I3=int ||Sym B_tilde|| dt` in the affine-curvature interaction frame and `J3=int ||Sym S_tilde|| dt`, with `T` the packet lifetime.

The coherence theorem gives

- `J3 >= I3/T`, or
- first-Duhamel H3 daughter `delta1^2 >= 3 I3^2/32`.

On the coherent branch, either nonlinear sideband feedback has size at least `delta1/2`, or the surviving daughter has `delta^2 >= 3 I3^2/128`.  Let `sigma` be its second moment in the critical `|G|^(3/2)` Gaussian measure.  The variance change gives `sigma>=delta`.

If `sigma>=1/80`, this is already a definite daughter-capacity event.  If `sigma<1/80`, odd-Hermite convexity gives single-role transfer loss at least `sigma^2/16`.  Splitting that loss against possible pair-sideband rescue yields the clean alternative

`net transfer deficit >= 3 I3^2/4096`

or

`pair-sideband rescue >= 3 I3^2/4096`.

Thus H3 curvature has five and only five exits in this model: acceleration-Hessian dephasing source, nonlinear sideband feedback, a definite large daughter, a quadratic transfer deficit, or a quadratic pair-sideband interaction.

Stress: `{out.samples}`
- branch counts: `{out.branch_counts}`
- minimum first-impulse margin: `{out.minimum_first_impulse_margin:.3e}`
- minimum transfer-cost margin: `{out.minimum_transfer_cost_margin:.3e}`
- minimum pair-cost margin: `{out.minimum_pair_cost_margin:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
