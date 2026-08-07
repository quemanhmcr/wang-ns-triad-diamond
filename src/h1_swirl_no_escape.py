from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict,dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

H1_COST=Fraction(1,25600)
FULL_MILD_COST=Fraction(1,102400)


def h1_first_impulse_lower(hook_impulse:float)->float:
    if hook_impulse<0: raise ValueError('nonnegative impulse required')
    return hook_impulse/(2*math.sqrt(50.0))


def h1_post_feedback_daughter_lower(hook_impulse:float)->float:
    return .5*h1_first_impulse_lower(hook_impulse)


def h1_quadratic_cost(hook_impulse:float)->float:
    if hook_impulse<0: raise ValueError('nonnegative impulse required')
    return float(H1_COST)*hook_impulse*hook_impulse


def h1_dephasing_threshold(hook_impulse:float,lifetime:float)->float:
    if hook_impulse<0 or lifetime<=0: raise ValueError('invalid data')
    return hook_impulse/(math.sqrt(50.0)*lifetime)


def classify_h1_mild_no_escape(
    hook_impulse:float,lifetime:float,covariant_forcing_variation:float,
    first_duhamel_norm:float,nonlinear_feedback_norm:float,
    critical_sigma:float,pair_rescue:float,net_transfer_deficit:float,
)->dict[str,float|str]:
    vals=[hook_impulse,covariant_forcing_variation,first_duhamel_norm,nonlinear_feedback_norm,critical_sigma,pair_rescue,net_transfer_deficit]
    if lifetime<=0 or any(v<0 for v in vals): raise ValueError('invalid H1 no-escape data')
    I=hook_impulse; varth=h1_dephasing_threshold(I,lifetime); req=h1_first_impulse_lower(I); cost=h1_quadratic_cost(I)
    if covariant_forcing_variation>=varth-1e-12*max(1.,varth): branch='H1_covariant_dephasing'
    else:
        if first_duhamel_norm+1e-11*max(1.,req)<req: raise ValueError('coherent H1 branch violates first-Duhamel lower')
        if nonlinear_feedback_norm>=.5*first_duhamel_norm-1e-12*max(1.,first_duhamel_norm): branch='nonlinear_sideband_feedback'
        else:
            actual=max(0.,first_duhamel_norm-nonlinear_feedback_norm)
            if critical_sigma+1e-11<actual: raise ValueError('critical H1 sigma must dominate L2 daughter; variance factor is 4/3')
            if critical_sigma>=1/80-1e-14: branch='large_daughter_capacity'
            else:
                d0=critical_sigma**2/16
                if pair_rescue>=.5*d0-1e-13*max(1.,d0):
                    branch='pair_sideband_rescue'
                    if pair_rescue+2e-12<cost: raise AssertionError('H1 pair cost failed')
                else:
                    branch='transfer_deficit'
                    if net_transfer_deficit+2e-12<.5*d0 or net_transfer_deficit+2e-12<cost: raise AssertionError('H1 deficit cost failed')
    return {'branch':branch,'dephasing_threshold':varth,'first_duhamel_required':req,'clean_quadratic_cost':cost}


def full_curvature_channel(I_full:float,I_h3:float,I_hook:float)->str:
    """From I_full<=sqrt(6) I_h3+I_hook, either H3>=I/(2sqrt6) or hook>=I/2."""
    if min(I_full,I_h3,I_hook)<0: raise ValueError('nonnegative impulses required')
    if I_full>math.sqrt(6)*I_h3+I_hook+1e-11*max(1.,I_full): raise ValueError('impulses violate curvature split upper bound')
    if I_h3>=I_full/(2*math.sqrt(6))-1e-13: return 'H3'
    if I_hook+1e-12<I_full/2: raise AssertionError('full curvature channel dichotomy failed')
    return 'H1_hook'


def full_mild_aspect_quadratic_cost(I_full:float)->float:
    if I_full<0: raise ValueError('nonnegative impulse required')
    return float(FULL_MILD_COST)*I_full*I_full


def exact_constant_certificate()->dict[str,str]:
    # H1: forcing^2 >=1/50 B_hook^2.  First Duhamel half -> 1/(4*50)=1/200.
    # Half feedback ->1/800; single-role deficit /16 ->1/12800; rescue split /2 ->1/25600.
    if Fraction(1,200)*Fraction(1,4)*Fraction(1,16)*Fraction(1,2)!=H1_COST: raise AssertionError('H1 constant mismatch')
    # full curvature: H1 branch I_hook>=I_full/2 -> cost >=(1/25600)*(1/4)=1/102400.
    if H1_COST*Fraction(1,4)!=FULL_MILD_COST: raise AssertionError('full mild cost mismatch')
    # H3 branch is stronger: (3/4096)*(1/24)=1/32768 >1/102400.
    if not Fraction(1,32768)>FULL_MILD_COST: raise AssertionError('H3 branch should dominate clean full cost')
    return {
        'H1_first_duhamel_square':'I_1^2/200','H1_post_feedback_square':'I_1^2/800',
        'H1_pair_or_deficit':'I_1^2/25600','full_mild_pair_or_deficit':'I_B^2/102400',
        'aspect_condition':'cond(L)<=21/20 throughout lifetime','status':'EXACT_GIVEN_MILD_ASPECT_H1_BRIDGE_AND_SIDEBAND_THEOREMS'
    }


@dataclass(frozen=True)
class H1NoEscapeStress:
    samples:int
    minimum_H1_pair_margin:float
    minimum_H1_deficit_margin:float
    minimum_full_channel_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260807)->H1NoEscapeStress:
    rng=np.random.default_rng(seed); mp=md=mf=float('inf'); counts={}
    for _ in range(samples):
        I=float(10**rng.uniform(-4,-.5)); T=float(10**rng.uniform(-3,0)); req=h1_first_impulse_lower(I); cost=h1_quadratic_cost(I)
        mode=int(rng.integers(0,5))
        if mode==0: J=(1+rng.random())*h1_dephasing_threshold(I,T); d1=fb=sigma=rescue=deficit=0.
        else:
            J=float(rng.uniform(0,.99))*h1_dephasing_threshold(I,T); d1=(1+rng.random())*req
            if mode==1: fb=float(rng.uniform(.5,1.2))*d1; sigma=rescue=deficit=0.
            else:
                fb=float(rng.uniform(0,.49))*d1; actual=d1-fb
                if mode==2: sigma=max(1/80,actual)*(1+rng.random()); rescue=deficit=0.
                else:
                    upper=.99/80
                    if actual>=upper:
                        I=1e-4; req=h1_first_impulse_lower(I); d1=1.2*req; fb=.1*d1; actual=d1-fb; cost=h1_quadratic_cost(I); J=.1*h1_dephasing_threshold(I,T)
                    sigma=float(rng.uniform(max(actual,1e-12),upper)); d0=sigma*sigma/16
                    if mode==3: rescue=float(rng.uniform(.5,1))*d0; deficit=0.; mp=min(mp,rescue-cost)
                    else: rescue=float(rng.uniform(0,.49))*d0; deficit=max(.5*d0,cost)*(1+rng.random()); md=min(md,deficit-cost)
        out=classify_h1_mild_no_escape(I,T,J,d1,fb,sigma,rescue,deficit); counts[out['branch']]=counts.get(out['branch'],0)+1
        # independent full-curvature channel regression
        Ifull=float(10**rng.uniform(-4,0)); Ih3=float(rng.uniform(0,Ifull/math.sqrt(6))); Ihook=max(0.,Ifull-math.sqrt(6)*Ih3)+float(rng.uniform(0,.2))*Ifull
        ch=full_curvature_channel(Ifull,Ih3,Ihook)
        if ch=='H3': margin=Ih3-Ifull/(2*math.sqrt(6))
        else: margin=Ihook-Ifull/2
        mf=min(mf,margin)
        if margin<-2e-12: raise AssertionError('full channel stress failed')
    return H1NoEscapeStress(samples,mp,md,mf,counts)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-h1-swirl-no-escape'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=exact_constant_certificate(); out=stress(args.samples)
    (args.outdir/'h1_swirl_no_escape.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f"""# Mild-aspect H1/swirl local no-escape theorem

Status: **{cert['status']}**.

Assume `cond(L(t))<=21/20` throughout one packet lifetime and let `I1=int||B_hook||dt`.  The physical relative-parent/child bridge gives `||f_H1^rel||_2 >= ||B_hook||/sqrt(50)`.  The Banach/Hilbert variation theorem therefore yields

`H1 covariant forcing variation >= I1/(sqrt(50) T)`

or a first-Duhamel daughter with `delta1^2>=I1^2/200`.

After the same half-feedback, odd-Hermite capacity and pair-rescue split used in the H3 theorem, the clean alternative is

`net transfer deficit >= I1^2/25600`

or

`pair-sideband rescue >= I1^2/25600`,

unless there is a large daughter (`sigma>=1/80`) or nonlinear sideband feedback.

The intrinsic curvature split also gives `I_B <= sqrt(6) I3 + I1`.  Hence either `I3>=I_B/(2sqrt(6))` or `I1>=I_B/2`.  Combining the H3 constant `3/4096` with the H1 constant gives the clean mild-aspect full-curvature cost

`pair rescue or transfer deficit >= I_B^2/102400`

outside the source/dephasing, nonlinear-feedback and large-daughter branches.

The H1 **covariant dephasing** branch is not yet expanded into pressure/SGS/viscous derivative sources with explicit constants.  That is the remaining physical source-calculus frontier.  High-aspect grains are not charged here; they remain in affine ancestry/reuse.

Stress: `{out.samples}`
- branch counts: `{out.branch_counts}`
- minimum H1 pair-cost margin: `{out.minimum_H1_pair_margin:.3e}`
- minimum H1 deficit-cost margin: `{out.minimum_H1_deficit_margin:.3e}`
- minimum full-channel margin: `{out.minimum_full_channel_margin:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
