from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.affine_coherent_moyal import periodic_discrete_stft


def gaussian_covariance_overlap(Sigma: np.ndarray, Theta: np.ndarray) -> float:
    """Exact L2 overlap of centered normalized real Gaussian windows with |g|^2 covariances Sigma,Theta."""
    S=np.asarray(Sigma,float); T=np.asarray(Theta,float)
    if S.shape!=(3,3) or T.shape!=(3,3): raise ValueError('3x3 covariances required')
    es=np.linalg.eigvalsh(S); et=np.linalg.eigvalsh(T)
    if es[0]<=0 or et[0]<=0: raise ValueError('positive definite covariances required')
    sign,ld=np.linalg.slogdet(S+T)
    if sign<=0: raise ValueError('bad covariance sum')
    logov=1.5*math.log(2.0)+0.25*(np.linalg.slogdet(S)[1]+np.linalg.slogdet(T)[1])-0.5*ld
    return float(math.exp(logov))


def generalized_log_covariance_distance(Sigma: np.ndarray, Theta: np.ndarray) -> float:
    """Affine-invariant SPD log distance from generalized eigenvalues."""
    S=np.asarray(Sigma,float); T=np.asarray(Theta,float)
    vals,vecs=np.linalg.eigh(S)
    if vals[0]<=0: raise ValueError('Sigma must be SPD')
    Sinvhalf=(vecs*(vals**-0.5))@vecs.T
    R=Sinvhalf@T@Sinvhalf
    lam=np.linalg.eigvalsh(0.5*(R+R.T))
    if lam[0]<=0: raise ValueError('Theta must be SPD')
    return float(np.linalg.norm(np.log(lam)))


def gaussian_window_distance(Sigma: np.ndarray, Theta: np.ndarray) -> float:
    ov=gaussian_covariance_overlap(Sigma,Theta)
    return math.sqrt(max(0.0,2.0-2.0*ov))


def clean_window_distance_upper(Sigma: np.ndarray, Theta: np.ndarray) -> float:
    return generalized_log_covariance_distance(Sigma,Theta)/(2.0*math.sqrt(2.0))


def covariance_energy_tv_upper(norm_f_sq: float, Sigma: np.ndarray, Theta: np.ndarray) -> float:
    """L1 total variation upper for |V_g f|^2 under a covariance-window change."""
    if norm_f_sq<0: raise ValueError('nonnegative energy required')
    return math.sqrt(2.0)*0.5*generalized_log_covariance_distance(Sigma,Theta)*norm_f_sq


def discrete_window_change_residual(f: np.ndarray, g: np.ndarray, h: np.ndarray) -> float:
    """Polarized Moyal on the window slot: ||V_g f-V_h f||^2=||f||^2||g-h||^2."""
    f=np.asarray(f,complex); g=np.asarray(g,complex); h=np.asarray(h,complex)
    if f.ndim!=1 or g.shape!=f.shape or h.shape!=f.shape: raise ValueError('same 1D shape required')
    A=periodic_discrete_stft(f,g); B=periodic_discrete_stft(f,h); n=len(f)
    lhs=float(np.sum(np.abs(A-B)**2)/n)
    rhs=float(np.vdot(f,f).real*np.vdot(g-h,g-h).real)
    return lhs-rhs


def discrete_energy_density_tv(f: np.ndarray, g: np.ndarray, h: np.ndarray) -> tuple[float,float]:
    f=np.asarray(f,complex); g=np.asarray(g,complex); h=np.asarray(h,complex); n=len(f)
    A=periodic_discrete_stft(f,g); B=periodic_discrete_stft(f,h)
    tv=float(np.sum(np.abs(np.abs(A)**2-np.abs(B)**2))/n)
    bound=float(np.linalg.norm(f)**2*np.linalg.norm(g-h)*(np.linalg.norm(g)+np.linalg.norm(h)))
    return tv,bound


def partition_energy_change(cell_old: Sequence[float], cell_new: Sequence[float]) -> float:
    a=np.asarray(cell_old,float); b=np.asarray(cell_new,float)
    if a.shape!=b.shape or np.any(a<0) or np.any(b<0): raise ValueError('matching nonnegative cell laws required')
    return float(np.sum(np.abs(a-b)))


def exact_certificate() -> dict[str,str]:
    return {
        'window_moyal':'||V_g f-V_h f||_2=||f||_2 ||g-h||_2',
        'gaussian_fidelity':'<g_Sigma,g_Theta>=2^(3/2)(det Sigma det Theta)^(1/4)/det(Sigma+Theta)^(1/2)',
        'eigenvalue_form':'overlap=product_i cosh(a_i/2)^(-1/2), exp(a_i)=generalized covariance eigenvalues',
        'clean_distance':'||g_Sigma-g_Theta||_2 <= ||a||_2/(2 sqrt(2))',
        'cell_tv':'sum_C |E_Sigma(C)-E_Theta(C)| <= d_log ||f||_2^2/sqrt(2)',
        'status':'EXACT_GAUSSIAN_FIDELITY_AND_MOYAL_WINDOW_STABILITY',
    }


@dataclass(frozen=True)
class CovarianceInterfaceStress:
    samples:int
    maximum_distance_ratio:float
    worst_overlap_formula_residual:float
    worst_window_moyal_relative_residual:float
    maximum_energy_tv_ratio:float


def random_spd(rng:np.random.Generator)->np.ndarray:
    Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    vals=np.exp(rng.uniform(-2,2,size=3))
    return Q@np.diag(vals)@Q.T


def direct_generalized_overlap(S:np.ndarray,T:np.ndarray)->float:
    vals,vecs=np.linalg.eigh(S); H=(vecs*(vals**-0.5))@vecs.T; lam=np.linalg.eigvalsh(H@T@H)
    return float(np.prod(np.cosh(0.5*np.log(lam))**-0.5))


def stress(samples:int=50_000,seed:int=20260808)->CovarianceInterfaceStress:
    rng=np.random.default_rng(seed); md=mo=mw=mt=0.0
    for _ in range(samples):
        S=random_spd(rng); T=random_spd(rng)
        ov=gaussian_covariance_overlap(S,T); ov2=direct_generalized_overlap(S,T)
        mo=max(mo,abs(ov-ov2))
        if abs(ov-ov2)>2e-12: raise AssertionError('Gaussian covariance overlap identity failed')
        d=gaussian_window_distance(S,T); ub=clean_window_distance_upper(S,T)
        if ub>1e-15: md=max(md,d/ub)
        if d>ub+2e-12: raise AssertionError('clean log-covariance window bound failed')
    for _ in range(min(samples,5000)):
        n=int(rng.integers(8,42)); f=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); h=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g); h/=np.linalg.norm(h)
        res=discrete_window_change_residual(f,g,h); scale=max(1.0,np.linalg.norm(f)**2*np.linalg.norm(g-h)**2); mw=max(mw,abs(res)/scale)
        if abs(res)>3e-10*scale: raise AssertionError('window-slot Moyal failed')
        tv,b=discrete_energy_density_tv(f,g,h); mt=max(mt,tv/max(b,1e-300))
        if tv>b+3e-10*max(1.0,b): raise AssertionError('coherent energy TV bound failed')
    return CovarianceInterfaceStress(samples,md,mo,mw,mt)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-coherent-covariance-interface'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True); cert=exact_certificate(); out=stress(args.samples)
    (args.outdir/'coherent_covariance_interface.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2))
    md=f'''# Coherent covariance interface: changing Gaussian analysis cells has a quantified Xi cost

Status: **{cert['status']}**.

For centered normalized Gaussian windows with physical covariances `Sigma,Theta`,

`<g_Sigma,g_Theta> = 2^(3/2)(det Sigma det Theta)^(1/4)/det(Sigma+Theta)^(1/2)`.

If `exp(a_i)` are the generalized covariance eigenvalues, this is `prod_i cosh(a_i/2)^(-1/2)`.  Since `log cosh x <= x^2/2`,

`||g_Sigma-g_Theta||_2 <= ||a||_2/(2 sqrt(2))`.

Polarized Moyal in the **window slot** is exact:

`||V_(g_Sigma) f - V_(g_Theta) f||_(L2 phase)^2 = ||f||_2^2 ||g_Sigma-g_Theta||_2^2`.

Hence the total variation of the positive coherent energy density obeys

`int ||V_Sigma f|^2-|V_Theta f|^2| dmu <= d_log(Sigma,Theta) ||f||_2^2/sqrt(2)`.

The same bound holds after any common phase-space partition.  Thus changing the representative Gaussian covariance inside a small covariance cell creates a deterministic transfer-interface error `Xi_cov`, while common affine center/carrier transport remains exact gauge.  Large covariance jumps are not hidden in this estimate and remain a genuine relink/strain/source branch.

Stress: `{out.samples}`
- maximum exact window distance / clean log bound: `{out.maximum_distance_ratio:.9f}`
- worst overlap-formula residual: `{out.worst_overlap_formula_residual:.3e}`
- worst window-Moyal relative residual: `{out.worst_window_moyal_relative_residual:.3e}`
- maximum coherent energy-TV / bound: `{out.maximum_energy_tv_ratio:.9f}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
