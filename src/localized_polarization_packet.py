from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.full_strain_observability import tracefree_2x2


def unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    n = float(np.linalg.norm(v))
    if n <= 1e-14:
        raise ValueError("zero vector")
    return v / n


def triad_normal(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    w = np.cross(unit(a), unit(b))
    n = float(np.linalg.norm(w))
    if n <= 1e-14:
        raise ValueError("degenerate triad")
    return w / n


def real_triad_frame(direction: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Columns (n x khat, n), matching the triad-normal helical gauge."""
    kh = unit(direction)
    n = triad_normal(a, b)
    e1 = np.cross(n, kh)
    e1 = unit(e1)
    return np.column_stack([e1, n])


def transverse_tracefree(S: np.ndarray, E: np.ndarray) -> np.ndarray:
    return tracefree_2x2(np.asarray(E, float).T @ np.asarray(S, float) @ np.asarray(E, float))


def polarization_generators(S: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    D1 = transverse_tracefree(S, real_triad_frame(a, a, b))
    D2 = transverse_tracefree(S, real_triad_frame(b, a, b))
    D3 = transverse_tracefree(S, real_triad_frame(c, a, b))
    return D1 - D2, D3


def kelvin_direction_rhs(khat: np.ndarray, A: np.ndarray) -> np.ndarray:
    kh = unit(khat)
    A = np.asarray(A, float)
    return -(np.eye(3) - np.outer(kh, kh)) @ A.T @ kh


def kelvin_direction_lipschitz_bound(A: np.ndarray, B: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """4 L |a-b| + |A-B| for L=max(||A||op,||B||op)."""
    L = max(float(np.linalg.norm(A, 2)), float(np.linalg.norm(B, 2)))
    return 4.0 * L * float(np.linalg.norm(unit(a) - unit(b))) + float(np.linalg.norm(np.asarray(A)-np.asarray(B), 2))


def kelvin_packet_direction_bound(c: float, sigma0: float, kappa: float, M: float, h: float) -> float:
    """Scale-free Gronwall bound at T=c N^-2.

    Assumes ||A||op <= sigma0 N^2 and ||A_x-A_0||op <= kappa M N^2.
    """
    return math.exp(4.0 * c * sigma0) * (h + c * kappa * M)


def generator_freezing_bound(deltaS_fro: float, strain_op: float, direction_delta: float) -> float:
    """Bound the pair norm sqrt(||d(D1-D2)||^2+||dD3||^2).

    Both triads are assumed to have sin(parent angle)>=9/10 and each of the
    three carrier directions differs by at most direction_delta.
    """
    return math.sqrt(5.0) * (deltaS_fro + 16.0 * strain_op * direction_delta)


def actual_generator_freezing_error(S0: np.ndarray, S1: np.ndarray, dirs0: tuple[np.ndarray,np.ndarray,np.ndarray], dirs1: tuple[np.ndarray,np.ndarray,np.ndarray]) -> float:
    G0, C0 = polarization_generators(S0, *dirs0)
    G1, C1 = polarization_generators(S1, *dirs1)
    return math.sqrt(float(np.sum((G1-G0)**2) + np.sum((C1-C0)**2)))


def low_strain_packet_error(c: float, sigma0: float, kappa: float, M: float, h: float) -> float:
    """Exact stated bound before rational simplification.

    This is the normalized additional polarization forcing coefficient after
    integrating the generator-freezing residual over T=c N^-2 and using the
    exact relative-polarization capacity bound.
    """
    delta = kelvin_packet_direction_bound(c, sigma0, kappa, M, h)
    eps_gen = math.sqrt(5.0) * (kappa * M + 16.0 * sigma0 * delta)
    return 2.0 * c * eps_gen


def simplified_low_strain_bound(c: float, kappa: float, M: float, h: float) -> float:
    """For c*sigma0<=1/30: E_pol <= 3h +(15/2)c kappa M."""
    return 3.0 * h + 7.5 * c * kappa * M


def combined_localization_optimum(a: float, b: float, c: float, kappa: float, h: float) -> tuple[float,float]:
    """Old a/M+b kappa M plus polarization correction."""
    if min(a, b, c, kappa) <= 0:
        raise ValueError("positive a,b,c,kappa required")
    beff = b + 7.5 * c
    M = math.sqrt(a / (beff * kappa))
    E = 3.0 * h + 2.0 * math.sqrt(a * beff * kappa)
    return M, E


def arb_constant_certificate() -> dict[str,str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    if not ((arb(2)/15).exp() < arb(6)/5):
        raise AssertionError("exp(2/15)<6/5 failed")
    if not (arb(5).sqrt() < arb(9)/4):
        raise AssertionError("sqrt(5)<9/4 failed")
    # Starting from 2 c sqrt(5)[kM+16 sigma exp(4c sigma)(h+c kM)]
    # and c sigma <=1/30 gives coefficients <72/25 for h and <369/50 for c kM.
    if not (arb(72)/25 < arb(3)):
        raise AssertionError("frequency simplification failed")
    if not (arb(369)/50 < arb(15)/2):
        raise AssertionError("curvature simplification failed")
    return {
        "exp_bound": "exp(2/15)<6/5",
        "sqrt5_bound": "sqrt(5)<9/4",
        "raw_h_coefficient_upper": "72/25",
        "simplified_h_coefficient": "3",
        "raw_ckappaM_coefficient_upper": "369/50",
        "simplified_ckappaM_coefficient": "15/2",
        "status": "CERTIFIED",
    }


def random_good_pair(rng: np.random.Generator) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
    # Rotate a good pair with angle in the certified interval cos in (1/4,2/5).
    cosang = float(rng.uniform(0.26, 0.39))
    sinang = math.sqrt(1.0-cosang*cosang)
    a = np.array([1.0,0.0,0.0])
    b = np.array([cosang,sinang,0.0])
    # random SO(3) by QR
    Q,_ = np.linalg.qr(rng.normal(size=(3,3)))
    if np.linalg.det(Q)<0: Q[:,0]*=-1
    a=Q@a; b=Q@b
    c=unit(a+b)
    return a,b,c


def perturb_direction(v: np.ndarray, rng: np.random.Generator, scale: float) -> np.ndarray:
    w = rng.normal(size=3)
    w -= np.dot(w,v)*v
    nw=float(np.linalg.norm(w))
    if nw<1e-14: return v.copy()
    w/=nw
    ang=float(rng.uniform(-scale,scale))
    return math.cos(ang)*v+math.sin(ang)*w


@dataclass(frozen=True)
class PacketStress:
    samples: int
    worst_kelvin_lipschitz_ratio: float
    worst_generator_bound_ratio: float
    worst_simplified_margin: float
    worst_optimality_residual: float


def stress(samples: int=50_000, seed: int=20260807) -> PacketStress:
    rng=np.random.default_rng(seed)
    wk=wg=wo=0.0
    wm=float("inf")
    for _ in range(samples):
        A=rng.normal(size=(3,3)); B=A+0.2*rng.normal(size=(3,3))
        a=unit(rng.normal(size=3)); b=unit(a+0.05*rng.normal(size=3))
        lhs=float(np.linalg.norm(kelvin_direction_rhs(a,A)-kelvin_direction_rhs(b,B)))
        rhs=kelvin_direction_lipschitz_bound(A,B,a,b)
        if lhs>rhs+2e-12: raise AssertionError("Kelvin direction Lipschitz bound violated")
        if rhs>1e-14: wk=max(wk,lhs/rhs)

        d0=random_good_pair(rng)
        scale=float(rng.uniform(1e-5,0.015))
        d1=tuple(perturb_direction(v,rng,scale) for v in d0)
        # Keep the perturbed parent pair in the certified angular moat.
        if np.linalg.norm(np.cross(d1[0],d1[1])) <=0.9:
            continue
        X=rng.normal(size=(3,3)); S0=0.5*(X+X.T); S0-=np.trace(S0)/3*np.eye(3)
        dS=0.05*rng.normal(size=(3,3)); dS=0.5*(dS+dS.T); dS-=np.trace(dS)/3*np.eye(3)
        S1=S0+dS
        delta=max(float(np.linalg.norm(d0[i]-d1[i])) for i in range(3))
        bound=generator_freezing_bound(float(np.linalg.norm(dS,'fro')), max(float(np.linalg.norm(S0,2)),float(np.linalg.norm(S1,2))), delta)
        actual=actual_generator_freezing_error(S0,S1,d0,d1)
        if actual>bound+3e-10: raise AssertionError("generator freezing bound violated")
        if bound>1e-14: wg=max(wg,actual/bound)

        c=float(rng.uniform(0.05,1.0)); sigma0=1.0/(30.0*c)*float(rng.uniform(0.05,1.0))
        kappa=10**float(rng.uniform(-6,-2)); M=float(rng.uniform(1.0,30.0)); h=10**float(rng.uniform(-6,-2))
        exact=low_strain_packet_error(c,sigma0,kappa,M,h)
        simp=simplified_low_strain_bound(c,kappa,M,h)
        wm=min(wm,simp-exact)
        if exact>simp+2e-12: raise AssertionError("simplified low-strain bound violated")

        aa=float(rng.uniform(0.1,3)); bb=float(rng.uniform(0.1,3))
        Mopt,Eopt=combined_localization_optimum(aa,bb,c,kappa,h)
        direct=aa/Mopt+(bb+7.5*c)*kappa*Mopt+3*h
        wo=max(wo,abs(direct-Eopt)/max(1.0,Eopt))
    return PacketStress(samples,wk,wg,wm,wo)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--outdir",type=Path,default=Path("results-localized-polarization"))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_constant_certificate(); out=stress(args.samples)
    (args.outdir/"localized_polarization_packet.json").write_text(json.dumps({"certificate":cert,"stress":out.__dict__},indent=2))
    md=f"""# Localized relative-polarization packet bridge

Status: **{cert['status']}** for the clean low-strain constants.

- Kelvin direction stability: `delta(T)<=exp(4 c sigma0)(h+c kappa M)`
- generator-freezing pair bound: `eps_pol<=sqrt(5)(deltaS_F+16 sigma delta_dir)`
- if `c sigma0<=1/30`, integrated additional polarization forcing:
  `E_pol <= 3 h +(15/2)c kappa M`
- combined with the old localization ledger:
  `E_total <= a/M +(b+15c/2) kappa M +3h`
- optimized spatial width: `M*=sqrt(a/((b+15c/2)kappa))`
- optimized error: `3h+2 sqrt(a(b+15c/2)kappa)`
- random checks: `{out.samples}`
- worst Kelvin Lipschitz ratio: `{out.worst_kelvin_lipschitz_ratio:.9f}`
- worst generator-bound ratio: `{out.worst_generator_bound_ratio:.9f}`
- minimum simplified-bound margin: `{out.worst_simplified_margin:.3e}`
- worst optimizer residual: `{out.worst_optimality_residual:.3e}`

Thus localized helical polarization introduces no new spatial error currency on
the low-strain branch: frequency-cell variation contributes a summable `O(h)`
term and spatial frame variation is absorbed into the same `kappa M` curvature
term already balanced against the `a/M` commutator.
"""
    (args.outdir/"summary.md").write_text(md); print(md)


if __name__=="__main__": main()
