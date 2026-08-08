from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

LOW_STRAIN_THRESHOLD = Fraction(1, 30)
TRANSPORTER_RADIUS = Fraction(1, 4)


def gradient_linf_bernstein_constant(radius_factor: float = 0.25) -> float:
    """C in ||grad V||_infty <= C N^(3/2)||grad V||_2.

    Uses unitary Fourier normalization and supp(Vhat) subset B_(radius_factor N).
    C=(2pi)^(-3/2) |B_1|^(1/2) radius_factor^(3/2).
    """
    if radius_factor <= 0:
        raise ValueError('positive support radius required')
    return (2.0*math.pi)**(-1.5)*math.sqrt(4.0*math.pi/3.0)*radius_factor**1.5


def normalized_dissipation_lower(
    strain_action: float,
    scaled_lifetime: float,
    radius_factor: float = 0.25,
) -> float:
    """D_V=N int||grad V||_2^2 dt forced by int||S||_op dt >= K.

    T=cN^-2. Since ||S||_op<=||grad V||_F, Bernstein plus time Cauchy gives
    K <= C sqrt(c D_V).
    """
    if strain_action < 0 or scaled_lifetime <= 0:
        raise ValueError('nonnegative action and positive scaled lifetime required')
    C = gradient_linf_bernstein_constant(radius_factor)
    return strain_action**2/(C*C*scaled_lifetime)


def clean_high_strain_dissipation_lower(scaled_lifetime: float) -> float:
    return normalized_dissipation_lower(float(LOW_STRAIN_THRESHOLD), scaled_lifetime, float(TRANSPORTER_RADIUS))


def geometric_physical_cost_sum(
    normalized_cost: float,
    base_frequency: float,
    scale_ratio: float,
    viscosity: float = 1.0,
) -> float:
    """Physical viscous cost for an infinite chain with D_j=normalized_cost.

    N_j=N0 q^j and int||grad V_j||_2^2 dt=D_j/N_j, so the infinite sum is
    nu D0/N0 * q/(q-1), finite for q>1.
    """
    if normalized_cost < 0 or base_frequency <= 0 or scale_ratio <= 1 or viscosity < 0:
        raise ValueError('invalid geometric-chain data')
    return viscosity*normalized_cost/base_frequency*scale_ratio/(scale_ratio-1.0)


def geometric_fresh_energy_sum(
    critical_mass: float,
    base_frequency: float,
    scale_ratio: float,
) -> float:
    """Energy of infinitely many critical fresh packets E_j=mu/N_j."""
    if critical_mass < 0 or base_frequency <= 0 or scale_ratio <= 1:
        raise ValueError('invalid critical chain data')
    return critical_mass/base_frequency*scale_ratio/(scale_ratio-1.0)


def arb_clean_certificate() -> dict[str,str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:
        raise RuntimeError('python-flint required') from exc
    ctx.prec=180
    pi=arb.pi()
    # C^2 = (2pi)^-3 * (4pi/3) * (1/4)^3 = 1/(384 pi^2).
    C2=(2*pi)**(-3)*(4*pi/3)*(arb(1)/4)**3
    target=arb(1)/(384*pi*pi)
    if not (abs(C2-target) < arb('1e-50')):
        raise AssertionError('clean Bernstein constant identity failed')
    d=(arb(1)/30)**2/(C2)  # c=1
    clean=arb(32)*pi*pi/75
    if not (abs(d-clean) < arb('1e-45')):
        raise AssertionError('high-strain clean dissipation identity failed')
    return {
        'unitary_fourier_C':'1/(8 sqrt(6) pi)',
        'C_squared':'1/(384 pi^2)',
        'high_strain':'K>1/30 => D_V>32 pi^2/(75 c)',
        'status':'ARB_CERTIFIED_HIGH_STRAIN_TO_CRITICAL_DISSIPATION',
    }


@dataclass(frozen=True)
class HighStrainStress:
    samples:int
    minimum_collision_margin:float
    minimum_geometric_sum_margin:float
    maximum_finite_chain_fraction:float


def stress(samples:int=50_000,seed:int=20260808)->HighStrainStress:
    rng=np.random.default_rng(seed)
    mc=mg=float('inf'); mf=0.0
    C=gradient_linf_bernstein_constant(.25)
    for _ in range(samples):
        c=float(math.exp(rng.uniform(-3,2)))
        D=float(math.exp(rng.uniform(-4,5)))
        # Saturate the analytic upper action then choose an admissible smaller action.
        Kup=C*math.sqrt(c*D)
        K=float(rng.uniform(0,1))*Kup
        low=normalized_dissipation_lower(K,c,.25)
        mc=min(mc,D-low)
        if low>D+2e-12*max(1.0,D):
            raise AssertionError('high-strain/dissipation collision failed')

        q=float(rng.uniform(1.05,3.0)); N0=float(math.exp(rng.uniform(-2,5)))
        d0=float(math.exp(rng.uniform(-4,3))); nu=float(rng.uniform(0,2))
        exact=geometric_physical_cost_sum(d0,N0,q,nu)
        L=int(rng.integers(1,80))
        finite=sum(nu*d0/(N0*q**j) for j in range(L))
        mg=min(mg,exact-finite)
        mf=max(mf,finite/exact if exact>0 else 0.0)
        if finite>exact+2e-12*max(1.0,exact):
            raise AssertionError('geometric physical-cost sum exceeded infinite budget')
    return HighStrainStress(samples,mc,mg,mf)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-high-strain-dissipation-collision'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_clean_certificate(); out=stress(args.samples)
    data={'certificate':cert,'stress':asdict(out)}
    (args.outdir/'high_strain_dissipation_collision.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f'''# High-strain lifetime to critical dissipation collision\n\nStatus: **{cert['status']}**.\n\nFor the strict transporter `V=S_(N/4)u`, unitary Fourier Cauchy--Schwarz gives\n\n`||grad V||_infty <= [N^(3/2)/(8 sqrt(6) pi)] ||grad V||_2`.\n\nOn a natural packet lifetime `T=c N^-2`, time Cauchy therefore yields\n\n`K:=int ||S||_op dt <= int||grad V||_infty dt <= sqrt(c D_V)/(8 sqrt(6) pi)`,\n\nwhere `D_V=N int||grad V||_2^2 dt`. Hence\n\n`D_V >= 384 pi^2 K^2/c`,\n\nand the old low-strain threshold gives the clean branch\n\n`K>1/30 => D_V>32 pi^2/(75 c)`.\n\nThis is a genuine physical collision but **not a global reset-count budget**. If `N_j=N_0 q^j` and every generation pays the same normalized `D_V=d_0`, the actual viscous energy cost is `nu d_0/N_j`, so\n\n`sum_(j>=0) nu d_0/N_j = (nu d_0/N_0) q/(q-1) < infinity`.\n\nLikewise infinitely many critical fresh packets with `N_j E_j=mu` have finite total energy `mu q/[N_0(q-1)]`. Thus critical energy/dissipation currencies cannot be inserted into the multi-currency master as if each event consumed one scale-independent global amount. They need the existing branching/reuse/entropy structure or a genuinely weighted telescope.\n\nStress: `{out.samples}` collision/geometric-chain checks\n- minimum collision margin: `{out.minimum_collision_margin:.3e}`\n- minimum geometric-sum margin: `{out.minimum_geometric_sum_margin:.3e}`\n- maximum finite/infinite chain fraction: `{out.maximum_finite_chain_fraction:.9f}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
