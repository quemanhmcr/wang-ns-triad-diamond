from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.curvature_sideband_irrep import hook_from_M
from src.full_strain_observability import extremal_parent_directions, tracefree_2x2, transverse_frame

RSTAR_LO=Fraction(61090410158,100_000_000_000)
RSTAR_HI=Fraction(61090410160,100_000_000_000)


def physical_hook_strain_slices(B:np.ndarray,L:np.ndarray)->np.ndarray:
    """C_c=Sym(L B_c L^-1), the physical strain gradient per grain coordinate."""
    B=np.asarray(B,float); L=np.asarray(L,float); Li=np.linalg.inv(L)
    C=np.zeros((3,3,3))
    for c in range(3):
        G=L@B[:,:,c]@Li
        C[:,:,c]=0.5*(G+G.T)
    return C


def polarization_only_observable_slice(S:np.ndarray,rstar:float=.610904101586766)->float:
    """||D1-D2||_F^2+||D3||_F^2; deliberately excludes scalar D_Pi."""
    S=np.asarray(S,float); k1,k2=extremal_parent_directions(rstar)
    D1=tracefree_2x2(transverse_frame(k1).T@S@transverse_frame(k1))
    D2=tracefree_2x2(transverse_frame(k2).T@S@transverse_frame(k2))
    D3=tracefree_2x2(transverse_frame(np.array([1.,0.,0.])).T@S@transverse_frame(np.array([1.,0.,0.])))
    return float(np.sum((D1-D2)**2)+np.sum(D3**2))


def polarization_only_observable(C:np.ndarray,rstar:float=.610904101586766)->float:
    C=np.asarray(C,float)
    return sum(polarization_only_observable_slice(C[:,:,c],rstar) for c in range(3))


def arb_isotropic_hook_certificate()->dict[str,str]:
    """Certify Q_pol(C_hook)>=1/10||B_hook||^2 on the r* bracket.

    Write M=[[a,b,x],[b,d,y],[x,y,-a-d]], B=hook(M), C_c=Sym B_c,
    C=cos(phi)^2=1/(4r^2). Direct expansion gives Q_pol-lambda||B||^2
    block diagonal in (a,d),b,x,y.  Here ||B||^2=12(a^2+ad+b^2+d^2+x^2+y^2).
    """
    try:
        from flint import arb,ctx
    except ImportError as exc: raise RuntimeError('python-flint required') from exc
    ctx.prec=160
    def aq(q:Fraction): return arb(q.numerator)/q.denominator
    r=aq(RSTAR_LO).union(aq(RSTAR_HI)); C=1/(4*r*r); lam=arb(1)/10
    # Q coefficients from exact hook contraction.
    A=-2*C*C-6*C+arb(17)/2-12*lam
    D=-2*C*C+4-12*lam
    L=4*C*C-12*C+10-12*lam
    Bc=3-2*C-12*lam
    X=-2*C*C+2*C+1-12*lam
    Y=-2*C*C+4-12*lam
    det=A*D-(L/2)*(L/2)
    for name,val in {'A':A,'D':D,'det_ad':det,'b':Bc,'x':X,'y':Y}.items():
        if not(val>arb(0)): raise AssertionError(f'isotropic hook polarization certificate failed {name}: {val}')
    # Mild-aspect perturbation arithmetic: (1/sqrt(10)-sqrt(5)/20)^2 >1/25.
    mild=(arb(1)/10).sqrt()-arb(5).sqrt()/20
    mild2=mild*mild
    if not(mild2>arb(1)/25): raise AssertionError(f'mild-aspect arithmetic failed: {mild2}')
    return {
        'rstar_ball':str(r),'C_ball':str(C),'isotropic_Qpol_lower':'1/10',
        'aspect_threshold':'21/20','mild_raw_square_ball':str(mild2),
        'mild_Qpol_lower':'1/25','spinor_H1_energy_lower':'1/50','status':'CERTIFIED'
    }


def mild_aspect_lower(condL:float)->float:
    if condL<1 or condL>21/20+1e-12: raise ValueError('mild-aspect branch requires 1<=cond(L)<=21/20')
    return max(0.,1/math.sqrt(10)-math.sqrt(5)*(condL-1))**2


def spinor_action_energy_from_qpol(qpol:float)->float:
    if qpol<0: raise ValueError('nonnegative qpol required')
    return .5*qpol


def random_stf_M(rng:np.random.Generator)->np.ndarray:
    X=rng.normal(size=(3,3)); M=.5*(X+X.T); M-=np.trace(M)/3*np.eye(3); return M


def random_mild_spd(rng:np.random.Generator,kmax:float=21/20)->np.ndarray:
    Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    # log axes centered so global scalar is removed; enforce exact condition <=kmax.
    vals=np.sort(rng.uniform(-.5,.5,size=3)); span=vals[-1]-vals[0]
    if span>0: vals*=math.log(kmax)/span
    axes=np.exp(vals); axes/=math.sqrt(axes.max()*axes.min())
    return Q@np.diag(axes)@Q.T


@dataclass(frozen=True)
class H1MildAspectStress:
    samples:int
    worst_isotropic_ratio:float
    worst_mild_ratio:float
    worst_mild_bound_margin:float
    worst_spinor_identity_residual:float
    maximum_condition_number:float


def stress(samples:int=50_000,seed:int=20260807)->H1MildAspectStress:
    rng=np.random.default_rng(seed); wi=wm=float('inf'); mb=float('inf'); ws=0.; kc=1.
    for _ in range(samples):
        M=random_stf_M(rng); B=hook_from_M(M); B2=float(np.sum(B*B))
        if B2<1e-16: continue
        Ci=physical_hook_strain_slices(B,np.eye(3)); qi=polarization_only_observable(Ci); wi=min(wi,qi/B2)
        if qi+2e-12<.1*B2: raise AssertionError('isotropic hook Qpol lower failed')
        L=random_mild_spd(rng); cond=float(np.linalg.cond(L)); kc=max(kc,cond)
        C=physical_hook_strain_slices(B,L); q=polarization_only_observable(C); wm=min(wm,q/B2)
        raw=mild_aspect_lower(cond); mb=min(mb,q/B2-raw)
        if q+3e-11<raw*B2 or q+3e-11<.04*B2: raise AssertionError('mild-aspect Qpol lower failed')
        # D^2=(||D||_F^2/2)I gives spinor-action identity for arbitrary complex unit vectors.
        for cc in range(3):
            k1,k2=extremal_parent_directions(.610904101586766)
            D1=tracefree_2x2(transverse_frame(k1).T@C[:,:,cc]@transverse_frame(k1))
            D2=tracefree_2x2(transverse_frame(k2).T@C[:,:,cc]@transverse_frame(k2))
            D3=tracefree_2x2(transverse_frame(np.array([1.,0.,0.])).T@C[:,:,cc]@transverse_frame(np.array([1.,0.,0.])))
            for D in (D1-D2,D3):
                v=rng.normal(size=2)+1j*rng.normal(size=2); v/=np.linalg.norm(v)
                lhs=float(np.vdot(D@v,D@v).real); rhs=.5*float(np.sum(D*D)); ws=max(ws,abs(lhs-rhs))
        if ws>2e-11: raise AssertionError('tracefree spinor action identity failed')
    return H1MildAspectStress(samples,wi,wm,mb,ws,kc)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-h1-swirl-mild-aspect'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_isotropic_hook_certificate(); out=stress(args.samples)
    (args.outdir/'h1_swirl_mild_aspect.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f"""# Physical H1/swirl bridge on mild-aspect affine grains

Status: **{cert['status']}**.

For the five-dimensional hook/swirl curvature sector, the physical polarization-only observable deliberately excludes the scalar shape term:

`Q_pol=sum_c (||D1(C_c)-D2(C_c)||_F^2+||D3(C_c)||_F^2)`.

At an isotropic grain, Arb certifies

`Q_pol >= (1/10)||B_hook||_F^2`.

For a general grain factor, after polar/global-scalar normalization its spectrum lies in `[cond(L)^(-1/2),cond(L)^(1/2)]`.  The physical slices obey `C_c=Sym(L B_c L^-1)` and

`||C(L)-C(I)|| <= (cond(L)-1)||B||`.

Since `sqrt(Q_pol(S))<=sqrt(5)||S||`, the observable triangle inequality gives

`sqrt(Q_pol(L)) >= [1/sqrt(10)-sqrt(5)(cond(L)-1)] ||B_hook||`.

Hence on the mild-aspect branch `cond(L)<=21/20`, Arb certifies the clean bound

`Q_pol(L) >= (1/25)||B_hook||^2`.

Every real symmetric trace-free 2x2 generator satisfies `D^2=(||D||_F^2/2)I`.  Thus for arbitrary unit helicity spinors the combined relative-parent/child H1 sideband forcing energy is exactly half of `Q_pol`, giving

`E ||F_H1^rel||^2 >= (1/50)||B_hook||^2`.

This is a transfer-facing **relative-parent/child polarization** statement; no `D_Pi` scalar-shape term is used.  It is intentionally only a mild-aspect theorem.  Grains with larger condition number are not charged by aspect: they remain in the affine fresh/reuse ancestry branch.

Stress: `{out.samples}`
- worst isotropic `Q_pol/||B||^2`: `{out.worst_isotropic_ratio:.9f}`
- worst mild-aspect ratio: `{out.worst_mild_ratio:.9f}`
- minimum margin above condition-dependent perturbation bound: `{out.worst_mild_bound_margin:.3e}`
- worst spinor-action identity residual: `{out.worst_spinor_identity_residual:.3e}`
- maximum condition number tested: `{out.maximum_condition_number:.9f}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
