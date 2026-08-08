from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping

import numpy as np


def additive_reset_count_upper(budgets: Mapping[str, tuple[float,float]]) -> float:
    """Upper bound sum_r B_r/b_r for single-charged additive reset events.

    Every reset event is assigned one primary resource r, consumes at least b_r,
    and total consumption of resource r is at most B_r.
    """
    total=0.0
    for name,(B,b) in budgets.items():
        if B<0 or b<=0:
            raise ValueError(f'invalid additive budget {name}')
        total += B/b
    return total


def transfer_costly_count_lower(
    depth: float,
    additive_resets: float,
    kappa0: float,
    potential_reset: float,
    potential_error: float=0.0,
) -> float:
    """Master lower bound for multiplicative transfer-cost blocks.

    Flat blocks obey N_F kappa0 <= (N_T+N_A+1)Pmax+Z and
    L=N_F+N_T+N_A. Solving gives the displayed lower bound.
    """
    if min(depth,additive_resets,potential_error)<0 or kappa0<=0 or potential_reset<0:
        raise ValueError('invalid master data')
    raw=(kappa0*depth-potential_reset-potential_error)/(kappa0+potential_reset)-additive_resets
    return max(0.0,raw)


def multicurrency_log_efficiency_lower(
    *,
    depth: float,
    budgets: Mapping[str,tuple[float,float]],
    c0: float,
    kappa0: float,
    potential_reset: float,
    potential_error: float=0.0,
    xi: float=0.0,
) -> dict[str,float]:
    if c0<=0 or xi<0:
        raise ValueError('positive c0 and nonnegative Xi required')
    na=additive_reset_count_upper(budgets)
    nt=transfer_costly_count_lower(depth,na,kappa0,potential_reset,potential_error)
    lower=c0*nt-xi
    rate=c0*kappa0/(kappa0+potential_reset)
    finite_offset=(
        c0*(potential_reset+potential_error)/(kappa0+potential_reset)
        + c0*na + xi
    )
    return {
        'depth':float(depth),
        'additive_reset_upper':na,
        'transfer_costly_lower':nt,
        'log_efficiency_lower':lower,
        'efficiency_upper':math.exp(-max(-700.0,min(700.0,lower))),
        'asymptotic_rate':rate,
        'finite_offset_upper':finite_offset,
    }


def direct_flat_ledger_margin(
    depth: int,
    transfer_costly: int,
    additive_resets: int,
    flat_steps: int,
    kappa0: float,
    potential_reset: float,
    potential_error: float,
) -> float:
    if transfer_costly+additive_resets+flat_steps!=depth:
        raise ValueError('counts must sum to depth')
    return (transfer_costly+additive_resets+1)*potential_reset + potential_error-flat_steps*kappa0


@dataclass(frozen=True)
class MultiCurrencyStress:
    samples: int
    minimum_transfer_count_margin: float
    minimum_log_efficiency_margin: float
    minimum_resource_count_margin: float
    minimum_asymptotic_rate: float


def stress(samples: int=50_000, seed: int=20260808) -> MultiCurrencyStress:
    rng=np.random.default_rng(seed)
    mt=ml=mr=float('inf'); rate_min=float('inf')
    for _ in range(samples):
        kappa=float(rng.uniform(.02,.5)); P=float(rng.uniform(.01,1.0)); Z=float(rng.uniform(0,.3))
        c0=float(10**rng.uniform(-6,-1))
        # Generate actual additive event counts and budgets with slack.
        classes=int(rng.integers(1,6)); budgets={}; actual_na=0
        for r in range(classes):
            b=float(10**rng.uniform(-4,-1)); n=int(rng.integers(0,20)); actual_na+=n
            B=b*(n+float(rng.uniform(0,2)))
            budgets[f'r{r}']=(B,b)
        NT=int(rng.integers(0,30))
        # Choose the largest flat count compatible with the exact episode ledger plus random slack.
        max_flat=max(0,int(math.floor(((NT+actual_na+1)*P+Z)/kappa)))
        NF=int(rng.integers(0,max_flat+1)) if max_flat>0 else 0
        L=NT+actual_na+NF
        if direct_flat_ledger_margin(L,NT,actual_na,NF,kappa,P,Z)<-1e-12:
            raise AssertionError('synthetic flat ledger is inadmissible')
        na_up=additive_reset_count_upper(budgets)
        mr=min(mr,na_up-actual_na)
        if na_up+1e-12<actual_na:
            raise AssertionError('resource budget did not dominate actual reset count')
        nt_low=transfer_costly_count_lower(L,na_up,kappa,P,Z)
        mt=min(mt,NT-nt_low)
        if nt_low>NT+2e-10:
            raise AssertionError('multicurrency master overcounted transfer-cost blocks')
        xi=float(rng.uniform(0,.2))
        out=multicurrency_log_efficiency_lower(depth=L,budgets=budgets,c0=c0,kappa0=kappa,potential_reset=P,potential_error=Z,xi=xi)
        true_log=c0*NT-xi
        ml=min(ml,true_log-out['log_efficiency_lower'])
        if out['log_efficiency_lower']>true_log+2e-10:
            raise AssertionError('log-efficiency lower bound exceeded actual synthetic cost')
        rate_min=min(rate_min,out['asymptotic_rate'])
    return MultiCurrencyStress(samples,mt,ml,mr,rate_min)


def theorem_certificate() -> dict[str,object]:
    return {
        'status':'EXACT_MULTICURRENCY_EPISODE_TELESCOPE',
        'single_charge_rule':'each additive reset selects exactly one primary globally bounded resource',
        'additive_reset_bound':'N_A <= sum_r B_r/b_r',
        'transfer_count':'N_T >= [kappa0 L-Pmax-Z]/[kappa0+Pmax]-N_A',
        'efficiency':'-log prod R_j >= c0 N_T-Xi',
        'rate':'c_eff=c0 kappa0/(kappa0+Pmax)',
        'physics':'only genuinely uniform globally bounded resource classes affect the finite offset; critical NE or normalized D_V are not automatically such classes',
    }


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-physical-multicurrency-master'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples); cert=theorem_certificate()
    (args.outdir/'physical_multicurrency_master.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f'''# Physical multi-currency master telescope\n\nStatus: **{cert['status']}**.\n\nSeparate non-flat resets into two kinds. A multiplicative transfer-cost block pays `C_j>=c0`. An additive physical-resource reset is assigned **one** primary currency `r`, consumes at least `b_r`, and that currency has total global budget `B_r`. Hence\n\n`N_A <= sum_r B_r/b_r`.\n\nIf flat blocks erode barycentric potential by `kappa0` up to total perturbation `Z`, while every transfer/additive reset may restart the potential below `Pmax`, then\n\n`N_F kappa0 <= (N_T+N_A+1)Pmax+Z`.\n\nSolving with `L=N_F+N_T+N_A` gives\n\n`N_T >= [kappa0 L-Pmax-Z]/[kappa0+Pmax]-N_A`.\n\nTherefore, with one global cross/interface penalty `Xi`,\n\n`-log prod_(j<L) R_j >= c0 N_T-Xi`\n\nand the asymptotic depth rate remains\n\n`c_eff=c0 kappa0/(kappa0+Pmax)>0`.\n\nAny reset class with a proved scale-independent threshold `b_r` in a globally bounded resource changes only the finite offset through `sum B_r/b_r`. Critical fresh mass `N E` and normalized dissipation `D_V` are **not** automatically eligible: their physical costs decay like `1/N` on a geometric scale chain. Entropy/Hodge/resistance/Renyi events that already pay multiplicative Bellman cost remain in `N_T`, not in the additive resource count.\n\nStress: `{out.samples}` synthetic multi-ledger episodes\n- minimum transfer-count margin: `{out.minimum_transfer_count_margin:.3e}`\n- minimum log-efficiency margin: `{out.minimum_log_efficiency_margin:.3e}`\n- minimum resource-count margin: `{out.minimum_resource_count_margin:.3e}`\n- minimum sampled asymptotic rate: `{out.minimum_asymptotic_rate:.3e}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__':
    main()
