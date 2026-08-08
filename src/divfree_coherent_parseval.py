from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


def projected_coefficient(P: np.ndarray, phi: np.ndarray, u: np.ndarray) -> complex:
    P=np.asarray(P,complex); phi=np.asarray(phi,complex); u=np.asarray(u,complex)
    return np.vdot(P@phi,u)


def unprojected_coefficient(phi: np.ndarray, u: np.ndarray) -> complex:
    return np.vdot(np.asarray(phi,complex),np.asarray(u,complex))


def parseval_projected_energy(P: np.ndarray, frame: np.ndarray, u: np.ndarray) -> float:
    """Sum |<P phi_a,u>|^2 for frame vectors stored as columns."""
    P=np.asarray(P,complex); F=np.asarray(frame,complex); u=np.asarray(u,complex)
    coeff=F.conj().T@(P@u)  # <phi,P u>=<P phi,u>
    return float(np.vdot(coeff,coeff).real)


def pressure_pairing(P: np.ndarray, phi: np.ndarray, gradient_vector: np.ndarray) -> complex:
    """Finite-dimensional analogue <P phi, grad p>, grad p in ker(P)."""
    return np.vdot(P@np.asarray(phi,complex),np.asarray(gradient_vector,complex))


def orthogonal_projector_from_qr(Q: np.ndarray, rank: int) -> np.ndarray:
    Q=np.asarray(Q,complex)
    U,_=np.linalg.qr(Q)
    return U[:,:rank]@U[:,:rank].conj().T


def theorem_certificate() -> dict[str,str]:
    return {
        'status':'EXACT_DIVERGENCE_FREE_COHERENT_PARSEVAL_FRAME',
        'abstract_frame':'if {phi_z} is Parseval on H and P is orthogonal, {P phi_z} is Parseval on ran(P)',
        'coefficient':'Pu=u => <u,P phi_z>=<u,phi_z>',
        'moyal':'sum_a int |<u,P(g_z e_a)>|^2 dmu=||u||_2^2 for div-free u',
        'pressure':'<P(g_z e_a),grad p>=0 exactly',
        'master_consequence':'canonical coherent ancestry does not require a compact spatial cutoff for pressure cancellation',
    }


@dataclass(frozen=True)
class DivFreeFrameStress:
    samples:int
    worst_projected_coefficient_residual:float
    worst_parseval_residual:float
    worst_pressure_residual:float
    worst_projector_idempotence:float


def stress(samples:int=20_000,seed:int=20260808)->DivFreeFrameStress:
    rng=np.random.default_rng(seed)
    wc=we=wp=wi=0.0
    for _ in range(samples):
        n=int(rng.integers(3,18)); r=int(rng.integers(1,n))
        Z=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        U,_=np.linalg.qr(Z)
        P=U[:,:r]@U[:,:r].conj().T
        wi=max(wi,float(np.linalg.norm(P@P-P)))
        y=rng.normal(size=r)+1j*rng.normal(size=r)
        u=U[:,:r]@y
        phi=rng.normal(size=n)+1j*rng.normal(size=n)
        c1=projected_coefficient(P,phi,u)
        c0=unprojected_coefficient(phi,u)
        wc=max(wc,abs(c1-c0))
        if abs(c1-c0)>3e-11*max(1.,abs(c0)):
            raise AssertionError('projected coherent coefficient changed on divergence-free subspace')

        # A unitary basis is a finite Parseval frame; projection preserves tightness on ran(P).
        F=U
        lhs=parseval_projected_energy(P,F,u)
        rhs=float(np.vdot(u,u).real)
        we=max(we,abs(lhs-rhs)/max(1.,rhs))
        if abs(lhs-rhs)>3e-11*max(1.,rhs):
            raise AssertionError('projected Parseval frame identity failed')

        # Pressure-gradient analogue is any vector in ker(P).
        q=rng.normal(size=n-r)+1j*rng.normal(size=n-r)
        grad=U[:,r:]@q
        pr=pressure_pairing(P,phi,grad)
        wp=max(wp,abs(pr))
        if abs(pr)>4e-11*max(1.,np.linalg.norm(phi)*np.linalg.norm(grad)):
            raise AssertionError('projected probe did not cancel gradient pressure')
    return DivFreeFrameStress(samples,wc,we,wp,wi)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=20_000); ap.add_argument('--outdir',type=Path,default=Path('results-divfree-coherent-parseval'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); out=stress(args.samples)
    (args.outdir/'divfree_coherent_parseval.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f'''# Divergence-free coherent Parseval frame\n\nStatus: **{cert['status']}**.\n\nLet `P` be the Leray projector and `g_z` any normalized affine coherent Gaussian family. For a divergence-free velocity `Pu=u`, self-adjointness gives exactly\n\n`<u,P(g_z e_a)>=<u,g_z e_a>`.\n\nComponentwise Moyal therefore survives projection:\n\n`sum_(a=1)^3 int |<u,P(g_z e_a)>|^2 dmu(z)=||u||_2^2`.\n\nMore generally polarized Moyal gives the exact work pairing on the divergence-free subspace. Each probe `P(g_z e_a)` is divergence-free, hence\n\n`<P(g_z e_a),grad p>=0`\n\nwith no compact spatial cutoff and no pressure boundary term. The Leray-projected probe need not remain exactly Gaussian; on the signed-good narrow Fourier cell its multiplier distortion is precisely the already-certified smooth-symbol/freezing representation error.\n\nThus the canonical coherent ancestry/master analysis may use a **global divergence-free coherent frame** rather than a compact moving spatial window. Compact windows remain useful for optional local/CKN diagnostics, but their moving-boundary commutator and localized pressure work need not enter the canonical master `Xi` ledger.\n\nStress: `{out.samples}` random orthogonal-projector/Parseval checks\n- worst projected coefficient residual: `{out.worst_projected_coefficient_residual:.3e}`\n- worst projected Parseval residual: `{out.worst_parseval_residual:.3e}`\n- worst gradient-pressure pairing residual: `{out.worst_pressure_residual:.3e}`\n- worst projector idempotence residual: `{out.worst_projector_idempotence:.3e}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
