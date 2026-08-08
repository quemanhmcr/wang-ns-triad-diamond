from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict,dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.spherical_erosion import C_STAR

BETA=Fraction(99,100)
TAU=Fraction(1,100)
CAP_CHORD=Fraction(3,10)
ENTROPY_CLEAN=Fraction(1,200)
CAP_MASS_CLEAN=Fraction(7,9)
CAP_GAP_CLEAN=Fraction(1,3)


def geometric_child_barycenter_error_upper(hodge_energy:float)->float:
    if hodge_energy<0: raise ValueError('nonnegative Hodge energy required')
    return 2.0*math.sqrt(hodge_energy)+0.5*hodge_energy


def collision_entropy_lower_from_barycenter(bary_norm:float)->float:
    if not 0<=bary_norm<=1: raise ValueError('barycenter norm in [0,1] required')
    return math.log(2.0/(1.0+bary_norm))


def barycenter_direction_separation_lower(hodge_energy:float,beta:float=float(BETA))->float:
    """Angle lower when both parent barycenter norms >=beta and |b_child|<=1."""
    if not 0<beta<=1: raise ValueError('beta in (0,1] required')
    e=geometric_child_barycenter_error_upper(hodge_energy)
    ratio=C_STAR*(1.0+e)/beta
    if ratio>=1: return 0.0
    return 2.0*math.acos(max(-1.0,min(1.0,ratio)))


def cap_mass_lower_from_barycenter(beta:float,chord_radius:float)->float:
    """Markov from E|X-bhat|^2=2(1-|b|)."""
    if not (0<beta<=1 and chord_radius>0): raise ValueError('bad cap data')
    return max(0.0,1.0-2.0*(1.0-beta)/(chord_radius*chord_radius))


def separated_cap_gap_chord_lower(direction_angle:float,chord_radius:float)->float:
    if direction_angle<0 or chord_radius<0: raise ValueError('nonnegative geometry required')
    center_chord=2.0*math.sin(direction_angle/2.0)
    return max(0.0,center_chord-2.0*chord_radius)


def classify_flat_companion(hodge_energy:float,old_bary_norm:float,companion_bary_norm:float)->dict[str,float|str]:
    vals=[hodge_energy,old_bary_norm,companion_bary_norm]
    if min(vals)<0 or max(old_bary_norm,companion_bary_norm)>1: raise ValueError('invalid data')
    beta=float(BETA)
    if old_bary_norm<=beta:
        h=collision_entropy_lower_from_barycenter(old_bary_norm)
        return {'branch':'old_parent_collision_entropy','entropy_lower':h,'clean_entropy_lower':float(ENTROPY_CLEAN)}
    if companion_bary_norm<=beta:
        h=collision_entropy_lower_from_barycenter(companion_bary_norm)
        return {'branch':'companion_collision_entropy','entropy_lower':h,'clean_entropy_lower':float(ENTROPY_CLEAN)}
    phi=barycenter_direction_separation_lower(hodge_energy,beta)
    mass=cap_mass_lower_from_barycenter(beta,float(CAP_CHORD))
    gap=separated_cap_gap_chord_lower(phi,float(CAP_CHORD))
    return {
        'branch':'two_trackable_parent_cores',
        'barycenter_direction_angle_lower':phi,
        'cap_chord_radius':float(CAP_CHORD),
        'mass_each_cap_lower':mass,
        'cap_gap_chord_lower':gap,
        'clean_angle_lower':1.0,
        'clean_mass_lower':float(CAP_MASS_CLEAN),
        'clean_gap_lower':float(CAP_GAP_CLEAN),
    }


def arb_flat_companion_certificate()->dict[str,str]:
    try:
        from flint import arb,ctx
    except ImportError as exc: raise RuntimeError('python-flint required') from exc
    ctx.prec=180
    # use certified r* bracket to enclose c*=1/(2r*)
    from src.single_edge_certificate import RSTAR_LO,RSTAR_HI
    def aq(q): return arb(q.numerator)/q.denominator
    r=aq(RSTAR_LO).union(aq(RSTAR_HI)); c=1/(2*r); tau=arb(1)/100; H=tau*tau/9; eg=2*H.sqrt()+H/2; beta=arb(99)/100
    ratio=c*(1+eg)/beta
    if not (ratio < (arb(1)/2).cos()): raise AssertionError(f'parent barycenter angle >1 failed: {ratio}')
    if not ((arb(200)/199).log()>arb(1)/200): raise AssertionError('clean entropy 1/200 failed')
    # Markov cap mass: 1-(1/50)/(9/100)=7/9 exactly.
    if Fraction(1,1)-Fraction(1,50)/Fraction(9,100)!=CAP_MASS_CLEAN: raise AssertionError('cap mass identity failed')
    # center chord at angle >1 is >2 sin(1/2), subtract two radius 3/10.
    gap=2*(arb(1)/2).sin()-arb(3)/5
    if not (gap>arb(1)/3): raise AssertionError(f'clean separated-cap gap failed {gap}')
    return {
        'tau':'1/100','Hodge_flat':'sqrt(H)<=tau/3','barycenter_threshold':'99/100',
        'low_barycenter_entropy':'log(200/199)>1/200',
        'both_concentrated_angle':'angle(b1,b2)>1 radian',
        'cap_chord_radius':'3/10','mass_each_cap':'at least 7/9','cap_chord_gap':'greater than 1/3',
        'status':'CERTIFIED_FLAT_ENTROPY_OR_TRACKABLE_COMPANION',
    }


@dataclass(frozen=True)
class FlatCompanionStress:
    samples:int
    minimum_entropy_margin:float
    minimum_angle_margin:float
    minimum_cap_gap_margin:float
    branch_counts:dict[str,int]


def stress(samples:int=50_000,seed:int=20260808)->FlatCompanionStress:
    rng=np.random.default_rng(seed); me=ma=mg=float('inf'); counts={}; Hmax=(.01/3)**2
    for _ in range(samples):
        H=float(rng.uniform(0,Hmax)); a=float(rng.uniform(.75,1)); b=float(rng.uniform(.75,1)); out=classify_flat_companion(H,a,b); br=str(out['branch']); counts[br]=counts.get(br,0)+1
        if 'entropy_lower' in out:
            me=min(me,float(out['entropy_lower'])-1/200)
            if float(out['entropy_lower'])<1/200-1e-13: raise AssertionError('flat companion entropy cost failed')
        else:
            ma=min(ma,float(out['barycenter_direction_angle_lower'])-1.0); mg=min(mg,float(out['cap_gap_chord_lower'])-1/3)
            if float(out['barycenter_direction_angle_lower'])<=1-1e-13 or float(out['mass_each_cap_lower'])<7/9-1e-13 or float(out['cap_gap_chord_lower'])<=1/3-1e-13:
                raise AssertionError('trackable companion core constants failed')
    return FlatCompanionStress(samples,me,ma,mg,counts)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-flat-companion-gate'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=arb_flat_companion_certificate(); out=stress(args.samples)
    (args.outdir/'flat_companion_gate.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f'''# Flat companion gate: entropy or two trackable parent cores

Status: **{cert['status']}**.

On a physical `tau=1/100` Kelvin-flat signed-good block, `H=E_H^phys<=tau^2/9`.  Without assuming equal parent marginals, triad geometry gives

`||b_child-(b_1+b_2)/(2c_*)|| <= 2 sqrt(H)+H/2`.

Choose `beta=99/100`.

If either parent barycenter has norm at most `beta`, the atomic barycenter--collision theorem gives

`H_2 >= log(200/199) >1/200`.

If both parent barycenter norms exceed `0.99`, physical `|b_child|<=1` forces their barycenter directions to be separated by more than `1` radian.  Each such marginal obeys `E|X-bhat|^2<1/50`; Markov at chord radius `3/10` therefore puts at least `7/9` of its transfer mass in a directional cap.  The two caps have chord gap greater than `1/3`.

Thus a `1%` Kelvin-flat block does not require an equal-marginal hypothesis: it either pays a fixed collision entropy `1/200`, or exposes two quantitatively separated, high-mass directional parent cores.  Relative to a distinguished old lineage, the second core is a trackable companion and must enter the existing fresh/reuse ancestry classification.

Stress: `{out.samples}`
- minimum entropy margin over `1/200`: `{out.minimum_entropy_margin:.3e}`
- minimum barycenter-angle margin over `1` radian: `{out.minimum_angle_margin:.3e}`
- minimum cap-gap margin over `1/3`: `{out.minimum_cap_gap_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
