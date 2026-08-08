from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

CURVATURE_DENOM = 1_036_800_000  # 36 * 28,800,000, from extended H1 channel
PHYSICAL_HODGE_FACTOR = Fraction(106,25)
EXTENDED_ASPECT = Fraction(567,500)
LOW_STRAIN_ACTION = Fraction(1,30)
PHASE_HOLONOMY_FLAT = Fraction(1,5)


def clean_thresholds(tau: float) -> dict[str,float]:
    if not (0 < tau <= 0.1):
        raise ValueError('certificate uses 0<tau<=1/10')
    delta=tau*tau/CURVATURE_DENOM
    return {
        'flatness_target':tau,
        'block_transfer_deficit':delta,
        'sideband_transfer_deficit':delta,
        'pair_rescue':delta,
        'h3_source_impulse':tau/(6.0*math.sqrt(6.0)), # T J3
        'h1_source_impulse':tau/792.0,               # T J1; 132*(T J1)<=tau/6
        'objective_strain_variation_action':tau/60.0, # T int ||Ddot_obj||dt
        'aspect_threshold':float(EXTENDED_ASPECT),
        'low_strain_action':float(LOW_STRAIN_ACTION),
        'phase_holonomy':float(PHASE_HOLONOMY_FLAT),
    }


def hodge_flatness_upper(avg_transfer_deficit: float) -> float:
    """sqrt(E_H^phys) <= sqrt((106/25) epsilon) on the physical low-cost core."""
    if avg_transfer_deficit < 0: raise ValueError('nonnegative deficit required')
    if avg_transfer_deficit >= 1/20_000:
        return math.inf
    return math.sqrt(float(PHYSICAL_HODGE_FACTOR)*avg_transfer_deficit)


def affine_strain_flatness_upper(avg_transfer_deficit: float, objective_variation_action: float) -> float:
    """Upper on d*T if the low-strain theorem applies.

    Coherent strain: avg Def >=(dT)^2/24.
    Incoherent strain: T int||Ddot_obj||dt >=dT/20.
    Hence without knowing which branch occurs, dT <=max(sqrt(24 eps),20 A_obj).
    """
    if min(avg_transfer_deficit,objective_variation_action)<0: raise ValueError('nonnegative inputs required')
    return max(math.sqrt(24.0*avg_transfer_deficit),20.0*objective_variation_action)


def curvature_flatness_upper(
    sideband_transfer_deficit: float,
    pair_rescue: float,
    h3_source_impulse: float,
    h1_source_impulse: float,
) -> float:
    """Full curvature impulse upper after excluding feedback/large-daughter branches.

    Use I_B<=sqrt(6)I3+I1 and the exact channel dichotomy:
      H3 channel: I3>=I_B/(2sqrt6), source if TJ3>=I3, otherwise cost>=3I3^2/4096.
      H1 channel: I1>=I_B/2, source if TJ1>=I1/132, otherwise cost>=I1^2/28,800,000.
    Thus I_B is bounded by the worse of the two channel-specific inversions below.
    """
    vals=[sideband_transfer_deficit,pair_rescue,h3_source_impulse,h1_source_impulse]
    if min(vals)<0: raise ValueError('nonnegative inputs required')
    C=max(sideband_transfer_deficit,pair_rescue)
    h3_cost_inverse=math.sqrt((4096.0/3.0)*C)
    h3=2.0*math.sqrt(6.0)*max(h3_source_impulse,h3_cost_inverse)
    h1_cost_inverse=math.sqrt(28_800_000.0*C)
    h1=2.0*max(132.0*h1_source_impulse,h1_cost_inverse)
    return max(h3,h1)


def connection_flatness_upper(
    avg_transfer_deficit: float,
    sideband_transfer_deficit: float,
    pair_rescue: float,
    h3_source_impulse: float,
    h1_source_impulse: float,
    objective_variation_action: float,
) -> dict[str,float]:
    h=hodge_flatness_upper(avg_transfer_deficit)
    s=affine_strain_flatness_upper(avg_transfer_deficit,objective_variation_action)
    c=curvature_flatness_upper(sideband_transfer_deficit,pair_rescue,h3_source_impulse,h1_source_impulse)
    return {'hodge_rms':h,'nonconformal_strain_number':s,'curvature_impulse':c,'kelvin_connection_flatness':h+s+c}


def classify_service_or_flat(
    tau: float,
    avg_transfer_deficit: float,
    sideband_transfer_deficit: float,
    pair_rescue: float,
    h3_source_impulse: float,
    h1_source_impulse: float,
    objective_variation_action: float,
    total_strain_action: float,
    aspect: float,
    has_predecessor: bool,
    nonlinear_sideband_feedback: bool=False,
    large_daughter: bool=False,
    weighted_phase_holonomy: float=0.0,
) -> dict[str,float|str|bool]:
    """Uniform service/action-or-Kelvin-flat packet theorem.

    The function only assembles existing theorem-level branches.  It does not
    assert that a continuum PDE extraction has already supplied these inputs.
    """
    th=clean_thresholds(tau)
    vals=[avg_transfer_deficit,sideband_transfer_deficit,pair_rescue,h3_source_impulse,h1_source_impulse,objective_variation_action,total_strain_action,weighted_phase_holonomy]
    if min(vals)<0 or aspect<1: raise ValueError('invalid block data')
    delta=th['block_transfer_deficit']
    if avg_transfer_deficit>=delta: return {'branch':'physical_transfer_cost','threshold':delta,'value':avg_transfer_deficit}
    if sideband_transfer_deficit>=delta: return {'branch':'sideband_transfer_cost','threshold':delta,'value':sideband_transfer_deficit}
    if pair_rescue>=delta: return {'branch':'pair_rescue_ancestry','threshold':delta,'value':pair_rescue}
    if nonlinear_sideband_feedback: return {'branch':'nonlinear_sideband_feedback'}
    if large_daughter: return {'branch':'large_daughter_capacity'}
    if h3_source_impulse>=th['h3_source_impulse']: return {'branch':'H3_physical_source','threshold':th['h3_source_impulse'],'value':h3_source_impulse}
    if h1_source_impulse>=th['h1_source_impulse']: return {'branch':'H1_physical_source','threshold':th['h1_source_impulse'],'value':h1_source_impulse}
    if objective_variation_action>=th['objective_strain_variation_action']: return {'branch':'objective_strain_source_action','threshold':th['objective_strain_variation_action'],'value':objective_variation_action}
    if total_strain_action>th['low_strain_action']:
        return {'branch':'high_strain_lifetime','threshold':th['low_strain_action'],'value':total_strain_action}
    if aspect>th['aspect_threshold']:
        return {'branch':'inherited_high_aspect' if has_predecessor else 'fresh_high_aspect','threshold':th['aspect_threshold'],'value':aspect}
    if weighted_phase_holonomy>=th['phase_holonomy']:
        return {'branch':'helical_phase_holonomy','threshold':th['phase_holonomy'],'value':weighted_phase_holonomy}
    flat=connection_flatness_upper(avg_transfer_deficit,sideband_transfer_deficit,pair_rescue,h3_source_impulse,h1_source_impulse,objective_variation_action)
    # The clean threshold design gives each of the three pieces <=tau/3.
    if flat['hodge_rms']>tau/3+2e-13 or flat['nonconformal_strain_number']>tau/3+2e-13 or flat['curvature_impulse']>tau/3+2e-13:
        raise AssertionError(('flatness threshold design failed',flat,th))
    return {'branch':'kelvin_extremal_flat','tau':tau,'phase_flat':True,**flat}


def arb_threshold_certificate(tau_num:int=1,tau_den:int=100)->dict[str,str]:
    try:
        from flint import arb,ctx
    except ImportError as exc: raise RuntimeError('python-flint required') from exc
    ctx.prec=180
    tau=arb(tau_num)/tau_den
    if not (tau>0 and tau<=arb(1)/10): raise ValueError('tau outside certificate range')
    delta=tau*tau/arb(CURVATURE_DENOM)
    # Hodge low-cost and tau/3 target.
    if not (((arb(106)/25)*delta).sqrt() < tau/3): raise AssertionError('Hodge flat threshold failed')
    if not ((24*delta).sqrt() < tau/3): raise AssertionError('coherent strain threshold failed')
    # H3 cost inversion is strictly stronger than target because extended H1 controls delta.
    h3=2*arb(6).sqrt()*((arb(4096)/3)*delta).sqrt()
    if not (h3 < tau/3): raise AssertionError('H3 curvature cost threshold failed')
    # H1 equality is exact algebraically: 2 sqrt(28.8m delta)=tau/3.
    if Fraction(4*28_800_000, CURVATURE_DENOM) != Fraction(1,9): raise AssertionError('H1 exact target identity failed')
    return {
        'tau':f'{tau_num}/{tau_den}',
        'uniform_transfer_pair_threshold':f'tau^2/{CURVATURE_DENOM}',
        'tau_1_100_threshold':str(Fraction(tau_num*tau_num,CURVATURE_DENOM*tau_den*tau_den)),
        'hodge_flat':'sqrt(E_H^phys)<=tau/3',
        'affine_nonconformal_strain_flat':'dT<=tau/3',
        'full_curvature_flat':'I_B<=tau/3',
        'H3_source_impulse':'T J3 < tau/(6 sqrt6)',
        'H1_source_impulse':'T J1 < tau/792',
        'objective_variation_action':'T int||Ddot_obj|| < tau/60',
        'status':'CERTIFIED_ASSEMBLY_GIVEN_EXISTING_PACKET_THEOREMS',
    }


@dataclass(frozen=True)
class ServiceFlatStress:
    samples:int
    minimum_flat_margin:float
    maximum_flat_ratio:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260808)->ServiceFlatStress:
    rng=np.random.default_rng(seed); mm=float('inf'); mr=0.0; counts={}; tau=.01; th=clean_thresholds(tau); delta=th['block_transfer_deficit']
    for _ in range(samples):
        mode=int(rng.integers(0,11))
        # start safely inside flat thresholds
        eps=float(rng.uniform(0,.8))*delta; sd=float(rng.uniform(0,.8))*delta; pr=float(rng.uniform(0,.8))*delta
        s3=float(rng.uniform(0,.8))*th['h3_source_impulse']; s1=float(rng.uniform(0,.8))*th['h1_source_impulse']; ao=float(rng.uniform(0,.8))*th['objective_strain_variation_action']; strain=float(rng.uniform(0,.8))*th['low_strain_action']; aspect=float(rng.uniform(1,th['aspect_threshold'])); pred=bool(rng.integers(0,2)); fb=ld=False; ph=float(rng.uniform(0,.8))*th['phase_holonomy']
        if mode==0: eps=1.1*delta
        elif mode==1: sd=1.1*delta
        elif mode==2: pr=1.1*delta
        elif mode==3: s3=1.1*th['h3_source_impulse']
        elif mode==4: s1=1.1*th['h1_source_impulse']
        elif mode==5: ao=1.1*th['objective_strain_variation_action']
        elif mode==6: strain=1.1*th['low_strain_action']
        elif mode==7: aspect=1.1*th['aspect_threshold']; pred=True
        elif mode==8: aspect=1.1*th['aspect_threshold']; pred=False
        elif mode==9: fb=True
        elif mode==10: pass
        out=classify_service_or_flat(tau,eps,sd,pr,s3,s1,ao,strain,aspect,pred,fb,ld,ph)
        b=str(out['branch']); counts[b]=counts.get(b,0)+1
        if b=='kelvin_extremal_flat':
            f=float(out['kelvin_connection_flatness']); mr=max(mr,f/tau); mm=min(mm,tau-f)
            if f>tau+3e-13: raise AssertionError('Kelvin flatness total exceeded tau')
    return ServiceFlatStress(samples,mm,mr,counts)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-service-or-flat-rigidity'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=arb_threshold_certificate(); out=stress(args.samples)
    data={'certificate':cert,'stress':asdict(out),'tau_1_100':clean_thresholds(.01)}
    (args.outdir/'service_or_flat_rigidity.json').write_text(json.dumps(data,indent=2))
    md=f'''# Uniform service/action-or-Kelvin-flat rigidity

Status: **{cert['status']}**.

For any `0<tau<=1/10`, set

`delta_tau = tau^2/{CURVATURE_DENOM}`.

On the signed-good physical low-cost core, the transition aspect strip `cond(L)<=567/500`, and the low-strain packet branch, suppose there is no nonlinear sideband-feedback or large-daughter event.  If all named currencies remain below

- block and sideband transfer deficit `<delta_tau`;
- pair rescue `<delta_tau`;
- H3 source impulse `T J3 < tau/(6 sqrt(6))`;
- H1 source impulse `T J1 < tau/792`;
- objective strain-variation action `<tau/60`;

then three independent rigidity estimates hold:

`sqrt(E_H^phys) <= tau/3`,

`(dT)_nonconformal <= tau/3`,

`I_B <= tau/3`.

Consequently the gauge-quotiented Kelvin connection flatness

`F_K = sqrt(E_H^phys) + (dT)_nonconformal + I_B`

satisfies `F_K<=tau`.

The curvature step uses the irreducible channel dichotomy rather than adding H1 and H3 constants: if full curvature is non-flat, either the H3 channel carries at least `I_B/(2 sqrt6)` or the hook H1 channel at least `I_B/2`; the latter fixes the clean denominator `{CURVATURE_DENOM}`.  Common affine/Kelvin motion is not charged.

Contrapositively, any block which is not `tau`-flat must enter a uniformly positive named branch: transfer/pair cost, H1/H3 physical source, objective strain source action, high-strain lifetime, inherited/fresh high aspect, nonlinear daughter/feedback, or helical phase holonomy.  The source branches feed the already-certified SGS/pressure/viscous coherent-service ledgers.

At the concrete value `tau=1/100`, the common transfer/pair threshold is `{Fraction(1,CURVATURE_DENOM*10000)}`.

Stress: `{out.samples}`
- maximum flatness / tau ratio: `{out.maximum_flat_ratio:.9f}`
- minimum flatness margin: `{out.minimum_flat_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
