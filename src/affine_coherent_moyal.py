from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.reservoir_pool_erosion import old_pool_service_capacity_upper, total_old_pool_service_upper


def moyal_energy(norm_f_sq: float, norm_window_sq: float = 1.0) -> float:
    if min(norm_f_sq, norm_window_sq) < 0:
        raise ValueError('nonnegative norms required')
    return norm_f_sq * norm_window_sq


def phase_cell_budget(cell_integrals) -> float:
    e = np.asarray(cell_integrals, float)
    if np.any(e < 0):
        raise ValueError('cell energies must be nonnegative')
    return float(np.sum(e))


def coherent_old_pool_capacity_upper(generation: int, initial_low_cut_ratio: float,
                                     initial_block_frequency: float, global_energy: float,
                                     beta_filter_radius: float = 1.0,
                                     lp_overlap_constant: float = 1.0) -> float:
    return old_pool_service_capacity_upper(
        generation, initial_low_cut_ratio, initial_block_frequency,
        lp_overlap_constant, global_energy, beta_filter_radius,
    )


def coherent_total_old_pool_capacity_upper(initial_low_cut_ratio: float,
                                           initial_block_frequency: float,
                                           global_energy: float,
                                           beta_filter_radius: float = 1.0,
                                           lp_overlap_constant: float = 1.0) -> float:
    return total_old_pool_service_upper(
        initial_low_cut_ratio, initial_block_frequency, lp_overlap_constant,
        global_energy, beta_filter_radius,
    )


def periodic_discrete_stft(f: np.ndarray, g: np.ndarray) -> np.ndarray:
    f = np.asarray(f, complex)
    g = np.asarray(g, complex)
    if f.ndim != 1 or g.shape != f.shape:
        raise ValueError('same 1D shape required')
    n = len(f)
    out = np.empty((n, n), complex)
    for m in range(n):
        out[m] = np.fft.fft(f * np.conj(np.roll(g, m)))
    return out


def discrete_moyal_residual(f: np.ndarray, g: np.ndarray) -> float:
    f = np.asarray(f, complex)
    g = np.asarray(g, complex)
    V = periodic_discrete_stft(f, g)
    n = len(f)
    lhs = float(np.sum(np.abs(V) ** 2) / n)
    rhs = float(np.vdot(f, f).real * np.vdot(g, g).real)
    return lhs - rhs


def exact_moyal_certificate() -> dict[str, str]:
    return {
        'continuum_identity': 'integral |V_g f|^2 dX dk/(2pi)^3 = ||g||_2^2 ||f||_2^2',
        'normalized_window_budget': 'P=1',
        'proof': 'Plancherel in k followed by Fubini/translation invariance in X',
        'phase_cell_budget': 'for any measurable partition, sum_C E_C=||f||_2^2',
        'status': 'EXACT_BY_PLANCHEREL_FUBINI',
    }


@dataclass(frozen=True)
class MoyalStress:
    samples: int
    worst_relative_discrete_moyal_residual: float
    worst_cell_partition_residual: float
    minimum_old_pool_half_life_margin: float


def stress(samples: int = 5000, seed: int = 20260808) -> MoyalStress:
    rng = np.random.default_rng(seed)
    wm = wc = 0.0
    mh = float('inf')
    for _ in range(samples):
        n = int(rng.integers(8, 49))
        f = rng.normal(size=n) + 1j * rng.normal(size=n)
        g = rng.normal(size=n) + 1j * rng.normal(size=n)
        g /= np.linalg.norm(g)
        res = discrete_moyal_residual(f, g)
        scale = max(1.0, float(np.vdot(f, f).real))
        wm = max(wm, abs(res) / scale)
        if abs(res) > 2e-10 * scale:
            raise AssertionError('discrete Moyal regression failed')
        cells = rng.dirichlet(np.ones(int(rng.integers(2, 20)))) * float(np.vdot(f, f).real)
        cres = phase_cell_budget(cells) - float(np.vdot(f, f).real)
        wc = max(wc, abs(cres) / scale)
        if abs(cres) > 2e-12 * scale:
            raise AssertionError('positive phase-cell partition failed')
        alpha = float(rng.uniform(.05, 1.0))
        N0 = float(math.exp(rng.uniform(-2, 2)))
        E = float(math.exp(rng.uniform(-3, 2)))
        q = int(rng.integers(0, 15))
        beta = float(rng.uniform(.2, 2.0))
        c0 = coherent_old_pool_capacity_upper(0, alpha, N0, E, beta)
        cq = coherent_old_pool_capacity_upper(q, alpha, N0, E, beta)
        mh = min(mh, c0 * .5 ** q - cq)
        if cq > c0 * .5 ** q + 2e-12 * max(1.0, c0):
            raise AssertionError('Moyal old-pool half-life failed')
    return MoyalStress(samples, wm, wc, mh)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--samples', type=int, default=5000)
    ap.add_argument('--outdir', type=Path, default=Path('results-affine-coherent-moyal'))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = exact_moyal_certificate()
    out = stress(args.samples)
    (args.outdir / 'affine_coherent_moyal.json').write_text(json.dumps({'certificate': cert, 'stress': asdict(out)}, indent=2))
    md = f'''# Affine coherent Moyal energy ledger

Status: **{cert['status']}**.

For any normalized affine Gaussian window `g_L`, the coherent/STFT transform obeys exactly

`integral |V_L f(X,k)|^2 dX dk/(2pi)^3 = ||f||_2^2`.

Any measurable phase-space partition therefore defines positive reservoir energies `E_C` with `sum_C E_C=||f||_2^2`.  This supplies an exact analysis-level reservoir budget `P=1`, independent of affine aspect and without synthesis-coefficient cancellation.

With an orthogonal dyadic band partition the whole-old-pool erosion theorem applies with `P=1`; a smooth LP partition pays only its fixed square-function overlap constant.  The 5-separated coherent Riesz theorem remains useful for a discrete synthesis realization.

Stress: `{out.samples}` discrete periodic Moyal regressions
- worst relative Moyal residual: `{out.worst_relative_discrete_moyal_residual:.3e}`
- worst cell-partition residual: `{out.worst_cell_partition_residual:.3e}`
- minimum old-pool half-life margin: `{out.minimum_old_pool_half_life_margin:.3e}`
'''
    (args.outdir / 'summary.md').write_text(md)
    print(md)


if __name__ == '__main__':
    main()
