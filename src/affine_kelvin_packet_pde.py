from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def unit(v: np.ndarray)->np.ndarray:
    v=np.asarray(v,float); n=float(np.linalg.norm(v))
    if n<1e-14: raise ValueError('zero vector')
    return v/n


def kelvin_direction_rate(n: np.ndarray,A: np.ndarray)->np.ndarray:
    n=unit(n); A=np.asarray(A,float)
    return -(np.eye(3)-np.outer(n,n))@A.T@n


def kelvin_amplitude_generator(A: np.ndarray,k: np.ndarray)->np.ndarray:
    """G_A(k)=-A+2 khat khat^T A in the affine linearized NS equation."""
    A=np.asarray(A,float); n=unit(k)
    return -A+2*np.outer(n,n)@A


def kelvin_amplitude_rate(A: np.ndarray,k: np.ndarray,a: np.ndarray,nu:float=0.0)->np.ndarray:
    k=np.asarray(k,float); a=np.asarray(a,complex)
    return kelvin_amplitude_generator(A,k)@a - nu*float(np.dot(k,k))*a


def transverse_frame(n: np.ndarray)->np.ndarray:
    n=unit(n)
    ref=np.array([1.,0.,0.]) if abs(n[0])<.8 else np.array([0.,1.,0.])
    e1=unit(np.cross(n,ref)); e2=np.cross(n,e1)
    return np.column_stack([e1,e2])


def objective_frame_rate(A: np.ndarray,n: np.ndarray,E: np.ndarray)->np.ndarray:
    """Frame derivative preserving E^T E=I, n^T E=0 and cancelling transverse spin.

    Let Omega=E^T dot E=-skew(E^T A E).  The normal-motion term is forced by
    d(n^T E)/dt=0.
    """
    A=np.asarray(A,float); n=unit(n); E=np.asarray(E,float)
    dn=kelvin_direction_rate(n,A)
    A2=E.T@A@E
    skew=.5*(A2-A2.T)
    Omega=-skew
    return -np.outer(n,dn)@E + E@Omega


def objective_coordinate_generator(A: np.ndarray,E: np.ndarray,nu_k2:float=0.0)->np.ndarray:
    A2=np.asarray(E,float).T@np.asarray(A,float)@np.asarray(E,float)
    S=.5*(A2+A2.T)
    return -S-nu_k2*np.eye(2)


def coordinate_rate_from_3d(A: np.ndarray,k: np.ndarray,E: np.ndarray,c: np.ndarray,nu:float=0.0)->np.ndarray:
    """Differentiate c=E^T a using full 3D Kelvin amplitude and objective frame."""
    n=unit(k); E=np.asarray(E,float); c=np.asarray(c,complex); a=E@c
    dE=objective_frame_rate(A,n,E)
    da=kelvin_amplitude_rate(A,k,a,nu)
    return dE.T@a+E.T@da


def affine_fourier_pressure_i(A: np.ndarray,k: np.ndarray,a: np.ndarray)->complex:
    """Return i*pi for the affine linearized pressure Fourier amplitude.

    i*pi = -2 (k.Aa)/|k|^2, so -i k pi gives the Kelvin pressure correction.
    """
    k=np.asarray(k,float); a=np.asarray(a,complex)
    return -2.0*np.dot(k,np.asarray(A,float)@a)/float(np.dot(k,k))


def affine_fourier_rhs_from_pressure(A: np.ndarray,k: np.ndarray,a: np.ndarray,nu:float=0.0)->np.ndarray:
    """Characteristic RHS from -A a - i k pi - nu |k|^2 a."""
    ipi=affine_fourier_pressure_i(A,k,a)
    # Since ipi=i*pi, the term -i k*pi is -k*(i*pi).
    return -np.asarray(A,float)@a - np.asarray(k,float)*ipi - nu*float(np.dot(k,k))*a


def scalar_cell_transport_residual(A: np.ndarray,xi: np.ndarray,grad_m:np.ndarray,dt_m:float)->float:
    return float(dt_m-np.dot(np.asarray(A,float).T@np.asarray(xi,float),np.asarray(grad_m,float)))


@dataclass(frozen=True)
class AffineKelvinPDEStress:
    samples:int
    worst_pressure_kelvin_residual:float
    worst_transversality_rate_residual:float
    worst_frame_orthogonality_rate_residual:float
    worst_frame_transverse_rate_residual:float
    worst_objective_coordinate_residual:float


def stress(samples:int=50_000,seed:int=20260807)->AffineKelvinPDEStress:
    rng=np.random.default_rng(seed); wp=wt=wo=wx=wc=0.0
    for _ in range(samples):
        A=rng.normal(size=(3,3)); A-=np.trace(A)/3*np.eye(3)
        k=rng.normal(size=3); n=unit(k); E=transverse_frame(n)
        # random rotate frame in its plane
        th=float(rng.uniform(-math.pi,math.pi)); R=np.array([[math.cos(th),-math.sin(th)],[math.sin(th),math.cos(th)]])
        E=E@R
        c=rng.normal(size=2)+1j*rng.normal(size=2); a=E@c
        nu=float(rng.uniform(0,1.5))
        da1=kelvin_amplitude_rate(A,k,a,nu); da2=affine_fourier_rhs_from_pressure(A,k,a,nu)
        wp=max(wp,float(np.linalg.norm(da1-da2)))
        if np.linalg.norm(da1-da2)>3e-12*max(1.,np.linalg.norm(da1)): raise AssertionError('pressure elimination/Kelvin generator mismatch')
        dk=-A.T@k
        trans=abs(np.dot(dk,a)+np.dot(k,da1))
        wt=max(wt,float(trans))
        if trans>5e-12*max(1.,np.linalg.norm(k)*np.linalg.norm(a)): raise AssertionError('k.a transversality not preserved')
        dn=kelvin_direction_rate(n,A); dE=objective_frame_rate(A,n,E)
        ort=float(np.linalg.norm(dE.T@E+E.T@dE)); wo=max(wo,ort)
        tr=float(np.linalg.norm(dn@E+n@dE)); wx=max(wx,tr)
        if ort>3e-12 or tr>3e-12: raise AssertionError('objective frame constraints failed')
        dc1=coordinate_rate_from_3d(A,k,E,c,nu)
        dc2=objective_coordinate_generator(A,E,nu*float(np.dot(k,k)))@c
        wc=max(wc,float(np.linalg.norm(dc1-dc2)))
        if np.linalg.norm(dc1-dc2)>5e-12*max(1.,np.linalg.norm(dc2)): raise AssertionError('objective 2x2 generator mismatch')
    return AffineKelvinPDEStress(samples,wp,wt,wo,wx,wc)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-affine-kelvin-packet-pde'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    data={'stress':out.__dict__,'equations':{'fourier_characteristic':'(partial_t-(A^T xi).grad_xi) w_hat = G_A(xi)w_hat-nu|xi|^2w_hat','G_A':'-A+2 khat khat^T A','kelvin_k':'dot k=-A^T k','objective_spinor':'dot c=-sym(E^T A E)c-nu|k|^2c'}}
    (args.outdir/'affine_kelvin_packet_pde.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Affine low-high Navier--Stokes packet equation

For the exact affine low-frequency jet `V=A(x-X)` in co-moving coordinates, the divergence-free linearized role obeys in Fourier space

`(partial_t-(A^T xi).grad_xi) w_hat = (-A+2 khat khat^T A) w_hat - nu|xi|^2 w_hat`.

Thus the characteristic laws are exactly
`dot k=-A^T k` and `dot a=(-A+2 khat khat^T A)a-nu|k|^2a`.
The pressure correction is exactly the `2 khat khat^T A` term; there is no additional packet pressure force after Leray projection.

In an objective transverse frame with `E^T dot E=-skew(E^T A E)`, the exact two-component coefficient equation is

`dot c=-sym(E^T A E)c-nu|k|^2c`,

the same generator used by the helical polarization ledger.  Hence the PDE low-high linearization, Kelvin carrier dynamics and objective helicity spinor are one identity, not separate modeling assumptions.

Stress checks: `{out.samples}`
- worst pressure/Kelvin residual: `{out.worst_pressure_kelvin_residual:.3e}`
- worst transversality-rate residual: `{out.worst_transversality_rate_residual:.3e}`
- worst frame orthogonality-rate residual: `{out.worst_frame_orthogonality_rate_residual:.3e}`
- worst frame transverse-constraint residual: `{out.worst_frame_transverse_rate_residual:.3e}`
- worst objective 2x2-coordinate residual: `{out.worst_objective_coordinate_residual:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
