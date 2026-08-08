from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.coherent_increment_service import cubic_to_square_threshold
from src.sgs_source_collision import cubic_increment_from_sgs_source_lower, fresh_radius_mass_lower
from src.source_episode_collision import h1_channel_normalized_integral_lower


def coherent_square_service_per_sgs_source(
    scale_radius_cap: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> float:
    """c_Y in Y_coh(rho)>=c_Y rho, exact by 3/2 -> 2/3 homogeneity."""
    q1=cubic_increment_from_sgs_source_lower(1.0,scale_radius_cap,filter_l1)
    return cubic_to_square_threshold(q1,filter_l1,lp_constant,bernstein_constant)


def coherent_sgs_episode_costs(
    I1: float,
    lifetime_c: float,
    scale_radius_cap: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    old_pool_capacity: float,
    source_divisor: float=132.0,
) -> dict[str,float|str]:
    """Source-weighted coherent alternatives with no persistence hypothesis.

    Sigma is the H1 source weight of one SGS channel.  Large-radius vs scale-matched
    costs one factor 1/2.  On scale-matched times, y=cY*rho.  High-enstrophy vs
    low coherent service costs another factor 1/2 in y-weight.  Old/interface/new
    then uses the coherent service-edge routing constants.
    """
    if old_pool_capacity<0: raise ValueError('nonnegative old capacity required')
    sigma=h1_channel_normalized_integral_lower(I1,lifetime_c,source_divisor)
    cy=coherent_square_service_per_sgs_source(scale_radius_cap,filter_l1,lp_constant,bernstein_constant)
    y_scale=0.5*cy*sigma
    high_diss=y_scale/8.0
    low_service=y_scale/4.0
    old_integrated=lifetime_c*old_pool_capacity
    old_target=cy*sigma/32.0
    xi_target=cy*sigma/32.0
    new_service=cy*sigma/16.0
    integrated_new_mass=cy*sigma/128.0
    peak_new_mass=integrated_new_mass/lifetime_c
    return {
        'total_source_weight':sigma,
        'coherent_square_service_per_source':cy,
        'large_radius_source_weight':0.5*sigma,
        'large_radius_mass':fresh_radius_mass_lower(scale_radius_cap),
        'scale_matched_coherent_y_weight':y_scale,
        'high_frequency_dissipation':high_diss,
        'integrated_low_coherent_service':low_service,
        'old_pool_integrated_capacity':old_integrated,
        'old_pool_erosion_target':old_target,
        'selected_interface_Xi':xi_target,
        'new_new_coherent_service':new_service,
        'integrated_new_coherent_critical_mass':integrated_new_mass,
        'peak_new_coherent_critical_mass':peak_new_mass,
        'ancestry_entropy':math.log(2.0),
        'same_ancestry_pair_mass':0.25,
        'status':'SOURCE_WEIGHTED_COHERENT_ROUTING_NO_PERSISTENCE',
    }


def source_homogeneity_residual(
    rho: float,
    scale_radius_cap: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
)->float:
    q=cubic_increment_from_sgs_source_lower(rho,scale_radius_cap,filter_l1)
    y=cubic_to_square_threshold(q,filter_l1,lp_constant,bernstein_constant)
    cy=coherent_square_service_per_sgs_source(scale_radius_cap,filter_l1,lp_constant,bernstein_constant)
    return y-cy*rho


@dataclass(frozen=True)
class CoherentSgsEpisodeStress:
    samples:int
    worst_relative_homogeneity_residual:float
    minimum_partition_margin:float
    minimum_peak_occupation_margin:float


def stress(samples:int=50_000,seed:int=20260808)->CoherentSgsEpisodeStress:
    rng=np.random.default_rng(seed); wh=0.0; mp=mo=float('inf')
    for _ in range(samples):
        s0=float(rng.uniform(.4,4)); g1=float(rng.uniform(1,2)); clp=float(rng.uniform(1,3)); cb=float(rng.uniform(1,2)); rho=float(rng.lognormal(-5,1.3))
        res=source_homogeneity_residual(rho,s0,g1,clp,cb); cy=coherent_square_service_per_sgs_source(s0,g1,clp,cb); scale=max(1e-300,cy*rho); wh=max(wh,abs(res)/scale)
        if abs(res)>2e-12*max(1.0,scale): raise AssertionError('SGS source -> coherent service homogeneity failed')
        I1=float(rng.uniform(1e-4,.3)); c=float(rng.uniform(.05,1)); old=float(rng.lognormal(-3,1)); out=coherent_sgs_episode_costs(I1,c,s0,g1,clp,cb,old,source_divisor=float(rng.choice([132.,1800.])))
        y=float(out['scale_matched_coherent_y_weight']); hd=float(out['high_frequency_dissipation']); low=float(out['integrated_low_coherent_service'])
        mp=min(mp,hd-y/8.0,low-y/4.0)
        occ=float(out['integrated_new_coherent_critical_mass']); peak=float(out['peak_new_coherent_critical_mass']); mo=min(mo,peak*c-occ)
        if peak*c+2e-13<occ: raise AssertionError('peak coherent occupation lower failed')
    return CoherentSgsEpisodeStress(samples,wh,mp,mo)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-coherent-sgs-episode'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); out=stress(args.samples)
    data={'status':'SOURCE_WEIGHTED_COHERENT_ROUTING_NO_PERSISTENCE','stress':asdict(out),'theorem':{
        'homogeneity':'rho_R -> Q>=c_Q rho_R^(3/2) -> Y_coh>=c_Y rho_R',
        'scale_matched_y_weight':'int Y >= c_Y Sigma_R/2',
        'dissipation':'D_high >= c_Y Sigma_R/16',
        'low_service':'int S_low >= c_Y Sigma_R/8',
        'old_target':'c C_old <= c_Y Sigma_R/32',
        'Xi':'int Xi_cell >= c_Y Sigma_R/32',
        'new_service':'int S_new >= c_Y Sigma_R/16',
        'new_mass_occupation':'int mu_new >= c_Y Sigma_R/128',
        'peak_new_mass':'sup mu_new >= c_Y Sigma_R/(128c)',
    }}
    (args.outdir/'coherent_sgs_episode.json').write_text(json.dumps(data,indent=2))
    md=f'''# Source-weighted SGS to coherent ancestry: no persistence hypothesis

Status: **SOURCE_WEIGHTED_COHERENT_ROUTING_NO_PERSISTENCE**.

The differentiated-SGS source density `rho_R` first forces cubic increment charge `Q>=c_Q rho_R^(3/2)`.  The coherent increment square threshold takes the `2/3` power, so the composition is exactly linear:

`Y_coh >= c_Y rho_R`.

Therefore integrated H1 SGS source weight `Sigma_R` produces integrated coherent service weight without assuming that the source persists above a pointwise threshold.  After the existing large-radius / scale-matched split, the scale-matched branch has `int Y >= c_Y Sigma_R/2`.  Hence either

`D_high >= c_Y Sigma_R/16`,

or

`int S_low >= c_Y Sigma_R/8`.

If the old material pool has integrated capacity `c C_old <= c_Y Sigma_R/32`, then either

`int Xi_cell >= c_Y Sigma_R/32`,

or new--new coherent service is at least `c_Y Sigma_R/16`.  A quarter-dominant new service edge gives

`int mu_coh,new >= c_Y Sigma_R/128`,

and therefore on a scaled lifetime of length at most `c`,

`sup_tau mu_coh,new >= c_Y Sigma_R/(128c)`.

If no dominant new edge exists, the same service-edge law pays `log 2` ancestry entropy or `1/4` same-ancestry pair/cycle mass.

Stress: `{out.samples}`
- worst source-homogeneity relative residual: `{out.worst_relative_homogeneity_residual:.3e}`
- minimum source-weight partition margin: `{out.minimum_partition_margin:.3e}`
- minimum peak-occupation margin: `{out.minimum_peak_occupation_margin:.3e}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
