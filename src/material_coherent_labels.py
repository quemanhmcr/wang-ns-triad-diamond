from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.smooth_symbol_freezing import sharp_young_constant_3d


def intrinsic_zeta(L: np.ndarray, X: np.ndarray, k: np.ndarray) -> np.ndarray:
    L=np.asarray(L,float); X=np.asarray(X,float); k=np.asarray(k,float)
    if L.shape!=(3,3) or X.shape!=(3,) or k.shape!=(3,):
        raise ValueError('expected L 3x3 and X,k in R3')
    return np.concatenate((0.5*np.linalg.solve(L,X), L.T@k))


def affine_transport_zeta_residual(M: np.ndarray, L: np.ndarray, X: np.ndarray, k: np.ndarray) -> float:
    M=np.asarray(M,float); L=np.asarray(L,float); X=np.asarray(X,float); k=np.asarray(k,float)
    if abs(np.linalg.det(M))<1e-12:
        raise ValueError('invertible affine map required')
    z0=intrinsic_zeta(L,X,k)
    z1=intrinsic_zeta(M@L, M@X, np.linalg.solve(M.T,k))
    return float(np.linalg.norm(z1-z0))


def dyadic_address(zeta: Sequence[float], depth: int, base_width: float=1.0) -> tuple[int,...]:
    z=np.asarray(zeta,float)
    if z.ndim!=1 or depth<0 or base_width<=0:
        raise ValueError('invalid dyadic cell data')
    width=base_width/(2**depth)
    return tuple(np.floor(z/width).astype(np.int64).tolist())


def parent_address(address: Sequence[int]) -> tuple[int,...]:
    return tuple(int(a)//2 for a in address)


def nested_address_chain(zeta: Sequence[float], max_depth: int, base_width: float=1.0) -> list[tuple[int,...]]:
    if max_depth<0:
        raise ValueError('nonnegative depth required')
    return [dyadic_address(zeta,j,base_width) for j in range(max_depth+1)]


def refinement_energy_residual(coarse_energy: float, child_energies: Sequence[float]) -> float:
    c=np.asarray(child_energies,float)
    if coarse_energy<0 or np.any(c<0):
        raise ValueError('nonnegative Moyal energies required')
    return float(np.sum(c)-coarse_energy)


def geometric_schedule(initial: float, depth: int, ratio: float=0.5) -> list[float]:
    if initial<0 or depth<0 or not (0<ratio<1):
        raise ValueError('invalid geometric schedule')
    return [initial*ratio**j for j in range(depth)]


def geometric_schedule_sum_upper(initial: float, ratio: float=0.5) -> float:
    if initial<0 or not (0<ratio<1):
        raise ValueError('invalid geometric schedule')
    return initial/(1.0-ratio)


def frequency_representation_xi_upper(
    initial_cell_diameter: float,
    symbol_lipschitz: float,
    role_norm_product_upper: float=1.0,
    ratio: float=0.5,
) -> float:
    if min(initial_cell_diameter,symbol_lipschitz,role_norm_product_upper)<0:
        raise ValueError('nonnegative frequency data required')
    return (
        sharp_young_constant_3d()
        * symbol_lipschitz
        * role_norm_product_upper
        * geometric_schedule_sum_upper(initial_cell_diameter,ratio)
    )


def covariance_representation_xi_upper(
    initial_log_mesh: float,
    global_energy: float,
    ratio: float=0.5,
) -> float:
    if min(initial_log_mesh,global_energy)<0:
        raise ValueError('nonnegative covariance data required')
    return global_energy*geometric_schedule_sum_upper(initial_log_mesh,ratio)/math.sqrt(2.0)


def total_representation_xi_upper(
    initial_cell_diameter: float,
    symbol_lipschitz: float,
    role_norm_product_upper: float,
    initial_log_mesh: float,
    global_energy: float,
    ratio: float=0.5,
) -> float:
    return frequency_representation_xi_upper(
        initial_cell_diameter,symbol_lipschitz,role_norm_product_upper,ratio
    ) + covariance_representation_xi_upper(initial_log_mesh,global_energy,ratio)


@dataclass(frozen=True)
class LabelStress:
    samples: int
    worst_affine_zeta_residual: float
    maximum_nested_address_failure: float
    worst_refinement_energy_residual: float
    minimum_frequency_schedule_margin: float
    minimum_covariance_schedule_margin: float


def stress(samples: int=50_000, seed: int=20260808) -> LabelStress:
    rng=np.random.default_rng(seed)
    wz=wr=0.0
    nested_fail=0.0
    mf=mc=float('inf')
    for _ in range(samples):
        # Moderately conditioned random affine geometry; identity is algebraic for any invertible M.
        Q,_=np.linalg.qr(rng.normal(size=(3,3)))
        s=np.exp(rng.uniform(-0.6,0.6,size=3))
        L=Q@np.diag(s)
        M=np.eye(3)+0.12*rng.normal(size=(3,3))
        if abs(np.linalg.det(M))<0.2:
            M += np.eye(3)
        X=rng.normal(size=3); k=rng.normal(size=3)
        res=affine_transport_zeta_residual(M,L,X,k)
        wz=max(wz,res)
        if res>2e-11:
            raise AssertionError('intrinsic material label changed under common affine transport')

        z=intrinsic_zeta(L,X,k)
        D=int(rng.integers(1,12))
        chain=nested_address_chain(z,D,base_width=float(rng.uniform(.4,2.0)))
        for j in range(D):
            ok=parent_address(chain[j+1])==chain[j]
            nested_fail=max(nested_fail,0.0 if ok else 1.0)
            if not ok:
                raise AssertionError('dyadic coherent label refinement is not nested')

        # Positive Moyal refinement is exactly additive by definition of a partition.
        n=64
        children=rng.random(n)
        coarse=float(children.sum())
        rr=refinement_energy_residual(coarse,children)
        wr=max(wr,abs(rr))
        if abs(rr)>3e-13*max(1.0,coarse):
            raise AssertionError('Moyal cell refinement created artificial switch mass')

        h0=float(rng.uniform(1e-5,.2)); delta0=float(rng.uniform(1e-5,.2))
        depth=int(rng.integers(1,40)); r=.5
        hs=sum(geometric_schedule(h0,depth,r))
        ds=sum(geometric_schedule(delta0,depth,r))
        mf=min(mf,geometric_schedule_sum_upper(h0,r)-hs)
        mc=min(mc,geometric_schedule_sum_upper(delta0,r)-ds)
        if hs>2*h0+1e-14 or ds>2*delta0+1e-14:
            raise AssertionError('geometric representation schedule lost summability')
    return LabelStress(samples,wz,nested_fail,wr,mf,mc)


def theorem_certificate() -> dict[str,object]:
    return {
        'status':'EXACT_CANONICAL_MATERIAL_LABELS_AND_SUMMABLE_REPRESENTATION_XI',
        'label':'dyadic address of intrinsic zeta=(L^-1 X/2,L^T k)',
        'common_transport':'exactly label preserving',
        'nested_refinement':'zero Moyal switch charge',
        'frequency_schedule':'h_j=h_0 2^-j, sum h_j<=2h_0',
        'frequency_xi':'Xi_sym<=2 A_3 L_* B_* h_0',
        'covariance_schedule':'delta_j=delta_0 2^-j, sum delta_j<=2delta_0',
        'covariance_xi':'Xi_cov<=sqrt(2) delta_0 E_global',
        'physical_switch':'not representation error; remains Moyal relink/backflow/fresh currency',
    }


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--samples',type=int,default=50_000)
    ap.add_argument('--outdir',type=Path,default=Path('results-material-coherent-labels'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples); cert=theorem_certificate()
    (args.outdir/'material_coherent_labels.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f'''# Canonical material coherent labels and summable representation Xi\n\nStatus: **{cert['status']}**.\n\nUse the intrinsic material coordinate `zeta=(L^-1 X/2,L^T k)` and a nested dyadic grid. Common affine/Kelvin transport leaves `zeta` and hence every cell address exactly unchanged. Refining a selected cell into its dyadic children is a partition identity, so Moyal energy is additive and refinement itself carries **zero** switch/relink charge.\n\nAt causal depth `j`, choose normalized frequency-cell diameter `h_j=h_0 2^-j`. Smooth SGS symbol freezing then has total representation error\n\n`Xi_sym <= 2 A_3 L_* B_* h_0`.\n\nChoose covariance representative mesh `delta_j=delta_0 2^-j`. The exact coherent covariance-interface theorem gives\n\n`Xi_cov <= sqrt(2) delta_0 E_global`.\n\nThus the causal Duhamel pushforward, Shannon/Renyi reuse and Hodge/resistance/holonomy graph may use the **same material dyadic cell address by construction**. Frequency/covariance representatives are auxiliary and have a summable, tunably small error. A genuine physical change of selected material cell is not hidden in this representation theorem; it remains the existing Moyal switch/fresh/relink/backflow currency.\n\nStress: `{out.samples}` affine/nested-grid/schedule checks\n- worst affine zeta residual: `{out.worst_affine_zeta_residual:.3e}`\n- maximum nested-address failure: `{out.maximum_nested_address_failure:.3e}`\n- worst refinement-energy residual: `{out.worst_refinement_energy_residual:.3e}`\n- minimum frequency schedule margin: `{out.minimum_frequency_schedule_margin:.3e}`\n- minimum covariance schedule margin: `{out.minimum_covariance_schedule_margin:.3e}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__':
    main()
