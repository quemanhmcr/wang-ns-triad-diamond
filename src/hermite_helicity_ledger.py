from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.polynomial.hermite_e import hermegauss


def sym3(T: np.ndarray) -> np.ndarray:
    T=np.asarray(T,float)
    out=np.zeros_like(T)
    import itertools
    perms=list(itertools.permutations(range(3)))
    for p in perms:
        out += np.transpose(T,p)
    return out/6.0


def h3_tensor(z: np.ndarray) -> np.ndarray:
    z=np.asarray(z,float); I=np.eye(3)
    H=np.einsum('a,b,c->abc',z,z,z)
    H-=np.einsum('ab,c->abc',I,z)
    H-=np.einsum('ac,b->abc',I,z)
    H-=np.einsum('bc,a->abc',I,z)
    return H


def third_hermite_polynomial(T: np.ndarray, z: np.ndarray) -> float:
    return float(np.einsum('abc,abc',np.asarray(T,float),h3_tensor(z)))


def third_hermite_norm_sq(T: np.ndarray) -> float:
    """E[(T:H3(Z))^2]=3! ||T||_F^2 for fully symmetric T."""
    Ts=sym3(T)
    return 6.0*float(np.sum(Ts*Ts))


def affine_envelope_sideband_norm_sq(B: np.ndarray) -> float:
    """Norm^2 of -(1/4) Sym(B):H3 times normalized Gaussian."""
    T=sym3(B)
    return 3.0/8.0*float(np.sum(T*T))


def first_hermite_sideband_norm_sq(C: np.ndarray) -> float:
    """For vector/tensor coefficients C[...,c], E|sum_c C_c Z_c|^2=||C||_F^2."""
    C=np.asarray(C,float)
    return float(np.sum(C*C))


def gauss_hermite_standard_normal(order: int=6) -> tuple[np.ndarray,np.ndarray]:
    """Tensor-product nodes/weights for standard normal expectation.

    hermegauss integrates exp(-x^2/2); divide weights by sqrt(2pi).
    """
    x,w=hermegauss(order); w=w/math.sqrt(2*math.pi)
    pts=[]; ws=[]
    for i in range(order):
        for j in range(order):
            for k in range(order):
                pts.append((x[i],x[j],x[k])); ws.append(w[i]*w[j]*w[k])
    return np.asarray(pts,float),np.asarray(ws,float)


def h3_projection_onto_degree2(T: np.ndarray, coeffs: tuple[float,np.ndarray,np.ndarray], order:int=6) -> float:
    """Gaussian inner product of T:H3 with arbitrary scalar polynomial degree<=2."""
    c0,c1,C2=coeffs; c1=np.asarray(c1,float); C2=np.asarray(C2,float)
    pts,ws=gauss_hermite_standard_normal(order)
    vals=[]
    for z in pts:
        p=float(c0+np.dot(c1,z)+z@C2@z)
        vals.append(third_hermite_polynomial(T,z)*p)
    return float(np.dot(ws,np.asarray(vals)))


def h1_projection_onto_base(C: np.ndarray, base: np.ndarray, order:int=4) -> np.ndarray:
    """Mean of (sum_c C[...,c] Z_c) dotted/projected against a constant base.

    For any constant base vector the n=0 projection vanishes.
    Returns the componentwise mean sideband before the final base contraction.
    """
    C=np.asarray(C,float); pts,ws=gauss_hermite_standard_normal(order)
    acc=np.zeros(C.shape[:-1])
    for z,w in zip(pts,ws):
        acc += w*np.tensordot(C,z,axes=([-1],[0]))
    return acc


def field_residual_ledger() -> dict[str,str]:
    return {
        "n0_base_helicity":"only this degree-zero projection enters F_i in the forced symplectic spinor identity",
        "n1_polarization_sideband":"spatial polarization curvature; orthogonal to the base Gaussian role and routed to coherence/branching",
        "n_le_2_scalar_tangent":"center/carrier/covariance/chirp and bulk viscosity; quotient as Gaussian tangent motion",
        "n3_envelope_sideband":"Sym B third-Hermite forcing; orthogonal to all scalar Gaussian tangent modes",
        "higher":"unresolved/nonquadratic packet forcing; route to cross-error/fresh-grain ledger",
    }


@dataclass(frozen=True)
class HermiteHelicityStress:
    samples: int
    worst_h3_degree2_projection: float
    worst_h3_norm_relative_residual: float
    worst_h1_base_projection: float
    worst_envelope_formula_residual: float


def stress(samples:int=20_000,seed:int=20260807)->HermiteHelicityStress:
    rng=np.random.default_rng(seed)
    wp=wn=wh1=we=0.0
    pts6,ws6=gauss_hermite_standard_normal(6)
    pts4,ws4=gauss_hermite_standard_normal(4)
    for _ in range(samples):
        T=sym3(rng.normal(size=(3,3,3)))
        c0=float(rng.normal()); c1=rng.normal(size=3); C2=rng.normal(size=(3,3)); C2=.5*(C2+C2.T)
        # Reuse precomputed quadrature for speed.
        proj=0.0; qnorm=0.0
        for z,w in zip(pts6,ws6):
            h=third_hermite_polynomial(T,z)
            p=c0+np.dot(c1,z)+z@C2@z
            proj += w*h*p; qnorm += w*h*h
        wp=max(wp,abs(proj))
        exact=third_hermite_norm_sq(T)
        wn=max(wn,abs(qnorm-exact)/max(1.0,exact))
        if abs(proj)>2e-10: raise AssertionError("H3 not orthogonal to degree<=2 Gaussian tangent")
        if abs(qnorm-exact)>3e-10*max(1.0,exact): raise AssertionError("H3 norm isometry failed")
        C=rng.normal(size=(2,3))
        acc=np.zeros(2)
        for z,w in zip(pts4,ws4): acc += w*(C@z)
        wh1=max(wh1,float(np.linalg.norm(acc)))
        if np.linalg.norm(acc)>2e-12: raise AssertionError("H1 sideband has nonzero base projection")
        B=rng.normal(size=(3,3,3)); Ts=sym3(B)
        lhs=affine_envelope_sideband_norm_sq(B); rhs=third_hermite_norm_sq(Ts)/16.0
        we=max(we,abs(lhs-rhs))
        if abs(lhs-rhs)>2e-12*max(1.0,lhs,rhs): raise AssertionError("envelope H3 norm coefficient failed")
    return HermiteHelicityStress(samples,wp,wn,wh1,we)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=20_000); ap.add_argument('--outdir',type=Path,default=Path('results-hermite-helicity-ledger'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    data={'ledger':field_residual_ledger(),'stress':out.__dict__,'identities':{'H3_norm':'E[(T:H3)^2]=6||T||_F^2','envelope_sideband':'||F_H3||^2/||psi||^2=(3/8)||Sym B||_F^2','H1_norm':'E|C.Z|^2=||C||_F^2','base_projection':'H1 and H3 have zero n=0 Gaussian projection'}}
    (args.outdir/'hermite_helicity_ledger.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Hermite-helicity forcing ledger

The field residual of a moving affine role must be projected by Hermite degree before it is inserted into the polarization equations.

- `n=0`: base Gaussian helicity spinor; only this component is the `F_i` in the forced symplectic identity.
- scalar degree `<=2`: center/carrier/covariance/chirp/bulk-viscosity tangent motion; quotient it.
- `n=1` vector polarization sideband: spatial polarization curvature; its base Gaussian projection is exactly zero.
- `n=3` scalar envelope sideband: `Sym B`; it is exactly orthogonal to every scalar Gaussian tangent polynomial of degree `<=2`.

Exact Gaussian-chaos identities:
- `E[(T:H3)^2]=6 ||T||_F^2`;
- `||F_H3||^2/||psi||^2=(3/8)||Sym B||_F^2`;
- `E|C.Z|^2=||C||_F^2`;
- both H1 and H3 have zero `n=0` projection.

Therefore third-Hermite leakage must not be double-counted as a direct spinor/phase force.  It is an orthogonal daughter-mode/coherence event.

Quadrature stress: `{out.samples}` random tensors
- worst H3/degree<=2 projection: `{out.worst_h3_degree2_projection:.3e}`
- worst H3 norm relative residual: `{out.worst_h3_norm_relative_residual:.3e}`
- worst H1/base projection: `{out.worst_h1_base_projection:.3e}`
- worst envelope coefficient residual: `{out.worst_envelope_formula_residual:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
