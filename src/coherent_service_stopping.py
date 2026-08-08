from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.ancestor_reservoir_sync import physical_energy_service_ratio_upper
from src.reservoir_pool_erosion import old_pool_service_capacity_upper


def clean_old_pool_ratio() -> Fraction:
    return physical_energy_service_ratio_upper()


def first_forced_cost_generation(service_threshold_y: float, initial_old_capacity: float) -> int:
    """First q with C_old(q)<=Y/8, where C_old(q)=C0 r^q and r<1/2.

    On an efficient coherent-service block with d_high<Y/4, the actual low service
    is >=Y/2.  Once old capacity is <=Y/8, the coherent increment theorem forces
    interface Xi, new coherent mass, entropy/cycle, or dissipation.
    """
    if service_threshold_y <= 0 or initial_old_capacity < 0:
        raise ValueError('positive Y and nonnegative old capacity required')
    target=service_threshold_y/8.0
    if initial_old_capacity<=target:
        return 0
    r=float(clean_old_pool_ratio())
    q=int(math.ceil(math.log(target/initial_old_capacity)/math.log(r)))
    # Correct possible floating endpoint ambiguity by exact monotone checks.
    q=max(0,q)
    while initial_old_capacity*r**q>target:
        q+=1
    while q>0 and initial_old_capacity*r**(q-1)<=target:
        q-=1
    return q


def forced_cost_generation_from_physical_pool(
    service_threshold_y: float,
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    global_energy: float,
    beta_filter_radius: float=1.0,
    frame_energy_bound: float=1.0,
)->int:
    c0=old_pool_service_capacity_upper(
        0,initial_low_cut_ratio,initial_block_frequency,frame_energy_bound,global_energy,beta_filter_radius
    )
    return first_forced_cost_generation(service_threshold_y,c0)


def epoch_certificate(service_threshold_y: float, initial_old_capacity: float) -> dict[str,float|int|str]:
    q=first_forced_cost_generation(service_threshold_y,initial_old_capacity)
    r=float(clean_old_pool_ratio()); target=service_threshold_y/8.0
    at=initial_old_capacity*r**q
    prev=initial_old_capacity*r**(q-1) if q>0 else initial_old_capacity
    return {
        'first_forced_generation':q,
        'maximum_cost_free_old_pool_blocks':q,
        'old_capacity_at_forced_generation':at,
        'target_old_capacity':target,
        'previous_old_capacity':prev,
        'one_step_ratio':r,
        'forced_alternatives':'d_high>=Y/4 OR Xi_cell>=Y/8 OR new coherent mass>=Y/32 OR H_anc>=log2 OR cycle mass>=1/4',
    }


def exact_certificate() -> dict[str,str]:
    r=clean_old_pool_ratio()
    return {
        'old_pool_ratio':f'{r.numerator}/{r.denominator}',
        'half_life':f'{r.numerator}/{r.denominator}<1/2',
        'stopping_condition':'C_old(q)<=Y0/8',
        'epoch_statement':'a cost-free old-pool coherent-service epoch has length < first q with C0 r^q<=Y0/8',
        'post_stop':'d_high>=Y0/4 OR Xi_cell>=Y0/8 OR new coherent mass>=Y0/32 OR log2 entropy OR 1/4 cycle mass',
        'status':'EXACT_GEOMETRIC_STOPPING_GIVEN_UNIFORM_COHERENT_SERVICE_THRESHOLD',
    }


@dataclass(frozen=True)
class CoherentStoppingStress:
    samples:int
    minimum_forced_margin:float
    minimum_previous_margin:float
    maximum_cost_free_blocks:int


def stress(samples:int=50_000,seed:int=20260808)->CoherentStoppingStress:
    rng=np.random.default_rng(seed); mf=mp=float('inf'); mq=0; r=float(clean_old_pool_ratio())
    for _ in range(samples):
        Y=float(rng.lognormal(-2,1.2)); C0=float(rng.lognormal(-1,1.5)); q=first_forced_cost_generation(Y,C0); target=Y/8
        cq=C0*r**q; mf=min(mf,target-cq); mq=max(mq,q)
        if cq>target+2e-13*max(1.0,target): raise AssertionError('forced generation still has too much old capacity')
        if q>0:
            cp=C0*r**(q-1); mp=min(mp,cp-target)
            if cp<=target-2e-13*max(1.0,target): raise AssertionError('stopping generation not minimal')
        else:
            mp=min(mp,0.0)
    return CoherentStoppingStress(samples,mf,mp,mq)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-coherent-service-stopping'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=exact_certificate(); out=stress(args.samples)
    (args.outdir/'coherent_service_stopping.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f'''# Coherent service stopping: a sticky old-pool epoch has finite length

Status: **{cert['status']}**.

The whole-old-pool service ratio is exactly

`r = {cert['old_pool_ratio']} < 1/2`.

Suppose every block in one sticky ancestry epoch carries a uniform coherent increment square-service threshold `Y>=Y0`.  If the high-frequency dissipation, selected-interface, fresh-coherent, entropy and cycle exits all fail, the coherent increment theorem requires old-pool capacity `>Y0/8`.

But `C_old(q)<=C0 r^q`.  Hence the epoch must stop by the first generation `q_*` satisfying

`C0 r^(q_*) <= Y0/8`.

At and after this stopping generation at least one named physical currency occurs:

- `d_high >= Y0/4`;
- `Xi_cell >= Y0/8`;
- new coherent critical mass `>=Y0/32`;
- ancestry Bellman entropy `>=log 2`; or
- same-ancestry cycle mass `>=1/4`.

A relinking event may start a new old-pool epoch, but that relinking is itself one of the charged exits.  Thus a uniform service threshold converts the geometric reservoir half-life into a genuine finite stopping-time statement suitable for the master episode ledger.

Stress: `{out.samples}`
- minimum forced-generation margin: `{out.minimum_forced_margin:.3e}`
- minimum previous-generation minimality margin: `{out.minimum_previous_margin:.3e}`
- largest sampled cost-free epoch length: `{out.maximum_cost_free_blocks}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
