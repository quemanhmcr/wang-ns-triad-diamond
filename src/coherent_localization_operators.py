from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


def cell_localization_operator(frame: np.ndarray, indices: Sequence[int]) -> np.ndarray:
    """A_C=sum_{a in C} phi_a phi_a^* for a finite Parseval frame.

    Frame vectors are columns. The continuous coherent theorem is the identical
    weak-operator formula with the sum replaced by the Moyal integral.
    """
    F=np.asarray(frame,complex)
    idx=np.asarray(list(indices),dtype=int)
    if F.ndim!=2:
        raise ValueError('frame matrix required')
    if idx.size==0:
        return np.zeros((F.shape[0],F.shape[0]),complex)
    Fc=F[:,idx]
    return Fc@Fc.conj().T


def partition_operators(frame: np.ndarray, cells: Sequence[Sequence[int]]) -> list[np.ndarray]:
    return [cell_localization_operator(frame,c) for c in cells]


def cell_energy(A: np.ndarray, f: np.ndarray) -> float:
    A=np.asarray(A,complex); f=np.asarray(f,complex)
    return float(np.vdot(f,A@f).real)


def synthesis_piece_energy(A: np.ndarray, f: np.ndarray) -> float:
    y=np.asarray(A,complex)@np.asarray(f,complex)
    return float(np.vdot(y,y).real)


def trilinear_tensor_value(T: np.ndarray, f: np.ndarray, g: np.ndarray, h: np.ndarray) -> complex:
    T=np.asarray(T,complex); f=np.asarray(f,complex); g=np.asarray(g,complex); h=np.asarray(h,complex)
    return np.einsum('ijk,i,j,k->',T,f,g,h,optimize=True)


def trilinear_partition_sum(T: np.ndarray, f: np.ndarray, g: np.ndarray, h: np.ndarray, Af: Sequence[np.ndarray], Ag: Sequence[np.ndarray], Ah: Sequence[np.ndarray]) -> complex:
    total=0j
    for A in Af:
        fa=A@f
        for B in Ag:
            gb=B@g
            for C in Ah:
                hc=C@h
                total += trilinear_tensor_value(T,fa,gb,hc)
    return total


def theorem_certificate() -> dict[str,str]:
    return {
        'status':'EXACT_COHERENT_LOCALIZATION_OPERATOR_DECOMPOSITION',
        'operator':'A_C=sum_a int_C |P(g_z e_a)><P(g_z e_a)| dmu',
        'partition':'sum_C A_C=I on L2_sigma',
        'positivity':'0<=A_C<=I',
        'cell_energy':'E_C=<f,A_C f>',
        'piece_budget':'||A_C f||_2^2<=E_C and sum_C||A_C f||_2^2<=||f||_2^2',
        'trilinear':'T(f,g,h)=sum_CDE T(A_C f,A_D g,A_E h) exactly for finite partitions',
        'interpretation':'no arbitrary synthesis/reconstruction Xi; only omitted physical cross-cell interactions remain',
    }


@dataclass(frozen=True)
class LocalizationStress:
    samples:int
    worst_partition_identity_residual:float
    minimum_positive_eigenvalue:float
    minimum_contraction_margin:float
    minimum_piece_energy_margin:float
    minimum_total_piece_budget_margin:float
    worst_trilinear_reconstruction_residual:float


def random_parseval_frame(rng: np.random.Generator,n:int,m:int)->np.ndarray:
    # Rows orthonormal => F F^*=I_n, so columns form a Parseval frame in C^n.
    Z=rng.normal(size=(m,n))+1j*rng.normal(size=(m,n))
    Q,_=np.linalg.qr(Z)
    return Q.conj().T  # n x m


def random_partition(rng: np.random.Generator,m:int,k:int)->list[list[int]]:
    labels=rng.integers(0,k,size=m)
    return [np.where(labels==j)[0].tolist() for j in range(k)]


def stress(samples:int=20_000,seed:int=20260808)->LocalizationStress:
    rng=np.random.default_rng(seed)
    wi=wt=0.; mp=mc=me=mb=float('inf')
    for _ in range(samples):
        n=int(rng.integers(2,8)); m=int(rng.integers(n,3*n+5)); k=int(rng.integers(1,min(m,7)+1))
        F=random_parseval_frame(rng,n,m); cells=random_partition(rng,m,k); ops=partition_operators(F,cells)
        S=sum(ops,np.zeros((n,n),complex)); ires=float(np.linalg.norm(S-np.eye(n)))
        wi=max(wi,ires)
        if ires>3e-11: raise AssertionError('localization operators do not resolve identity')
        for A in ops:
            ev=np.linalg.eigvalsh((A+A.conj().T)/2)
            mp=min(mp,float(ev.min(initial=0.0)))
            mc=min(mc,float(1.0-ev.max(initial=0.0)))
            if ev.min(initial=0.0)<-3e-11 or ev.max(initial=0.0)>1+3e-11:
                raise AssertionError('cell localization operator left [0,I]')
        f=rng.normal(size=n)+1j*rng.normal(size=n)
        energies=[cell_energy(A,f) for A in ops]
        pieces=[synthesis_piece_energy(A,f) for A in ops]
        for E,P in zip(energies,pieces):
            me=min(me,E-P)
            if P>E+3e-11*max(1.,E): raise AssertionError('canonical synthesis piece exceeded Moyal cell energy')
        nf=float(np.vdot(f,f).real); mb=min(mb,nf-sum(pieces))
        if sum(pieces)>nf+3e-11*max(1.,nf): raise AssertionError('canonical synthesis pieces lost Bessel P=1 budget')

        g=rng.normal(size=n)+1j*rng.normal(size=n); h=rng.normal(size=n)+1j*rng.normal(size=n)
        T=rng.normal(size=(n,n,n))+1j*rng.normal(size=(n,n,n))
        direct=trilinear_tensor_value(T,f,g,h)
        expanded=trilinear_partition_sum(T,f,g,h,ops,ops,ops)
        tres=abs(direct-expanded)/max(1.,abs(direct))
        wt=max(wt,tres)
        if tres>2e-10: raise AssertionError('trilinear coherent-cell reconstruction failed')
    return LocalizationStress(samples,wi,mp,mc,me,mb,wt)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=20_000); ap.add_argument('--outdir',type=Path,default=Path('results-coherent-localization-operators'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); out=stress(args.samples)
    (args.outdir/'coherent_localization_operators.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f'''# Coherent localization operators: exact synthesis without arbitrary packet coefficients\n\nStatus: **{cert['status']}**.\n\nFor the divergence-free coherent Parseval frame define the positive localization operator\n\n`A_C = sum_a int_C |P(g_z e_a)><P(g_z e_a)| dmu(z)`.\n\nFor every measurable partition, weakly on `L2_sigma`,\n\n`sum_C A_C=I`,  `0<=A_C<=I`.\n\nThe Moyal cell energy is `E_C=<f,A_C f>`. Since `A_C^2<=A_C`,\n\n`||A_C f||_2^2 <= E_C`,\n\nand summing gives the canonical synthesis budget\n\n`sum_C ||A_C f||_2^2 <= ||f||_2^2`.\n\nThus continuous coherent cells already provide synthesis pieces with frame constant `P=1`; no 5-separated representative family or arbitrary redundant coefficients are required for the canonical analysis.\n\nFor a finite cell partition and any continuous trilinear form, multilinearity plus `sum A_C=I` gives exactly\n\n`T(f,g,h)=sum_(C,D,E) T(A_C f,A_D g,A_E h)`.\n\nTherefore there is **zero reconstruction/synthesis Xi**. The remaining interface problem is only the actual physical cross-cell interaction mass discarded when one selects a dominant coherent ancestry component. Estimating that mass requires the Gaussian triad/resonance kernel; it is not a frame-reconstruction issue.\n\nStress: `{out.samples}` random finite Parseval-frame/localization/trilinear checks\n- worst partition identity residual: `{out.worst_partition_identity_residual:.3e}`\n- minimum positive eigenvalue: `{out.minimum_positive_eigenvalue:.3e}`\n- minimum contraction margin: `{out.minimum_contraction_margin:.3e}`\n- minimum cell-energy/piece margin: `{out.minimum_piece_energy_margin:.3e}`\n- minimum total piece-budget margin: `{out.minimum_total_piece_budget_margin:.3e}`\n- worst trilinear reconstruction residual: `{out.worst_trilinear_reconstruction_residual:.3e}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
