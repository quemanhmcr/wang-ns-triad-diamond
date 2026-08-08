from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.affine_coherent_moyal import periodic_discrete_stft


def discrete_cross_moyal(f: np.ndarray, h: np.ndarray, g: np.ndarray) -> complex:
    """Periodic polarized Moyal integral with the repository FFT convention."""
    f = np.asarray(f, complex); h = np.asarray(h, complex); g = np.asarray(g, complex)
    if f.ndim != 1 or h.shape != f.shape or g.shape != f.shape:
        raise ValueError('same 1D shape required')
    vf = periodic_discrete_stft(f, g)
    vh = periodic_discrete_stft(h, g)
    return complex(np.sum(vf * np.conj(vh)) / len(f))


def cell_work(Vf: np.ndarray, VF: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Signed coherent-cell work 2 Re int_C V(f) conjugate(V(F)).

    Discrete normalization is inherited from periodic_discrete_stft: divide by n.
    Labels may be any nonnegative integer partition of the discrete phase plane.
    """
    Vf = np.asarray(Vf, complex); VF = np.asarray(VF, complex); labels = np.asarray(labels)
    if Vf.shape != VF.shape or labels.shape != Vf.shape or Vf.ndim != 2:
        raise ValueError('matching 2D STFT arrays and labels required')
    if np.any(labels < 0):
        raise ValueError('nonnegative cell labels required')
    n = Vf.shape[1]
    m = int(labels.max()) + 1 if labels.size else 0
    out = np.zeros(m, float)
    density = 2.0 * np.real(Vf * np.conj(VF)) / n
    for a in range(m):
        out[a] = float(np.sum(density[labels == a]))
    return out


def cell_energy(Vf: np.ndarray, labels: np.ndarray) -> np.ndarray:
    Vf = np.asarray(Vf, complex); labels = np.asarray(labels)
    if labels.shape != Vf.shape or Vf.ndim != 2:
        raise ValueError('matching 2D STFT and labels required')
    n = Vf.shape[1]
    m = int(labels.max()) + 1 if labels.size else 0
    out = np.zeros(m, float)
    density = np.abs(Vf) ** 2 / n
    for a in range(m):
        out[a] = float(np.sum(density[labels == a]))
    return out


def symmetric_difference_energy(energies: Sequence[float], old_selected: Sequence[int], new_selected: Sequence[int]) -> float:
    e = np.asarray(energies, float)
    if np.any(e < 0):
        raise ValueError('nonnegative cell energies required')
    a, b = set(map(int, old_selected)), set(map(int, new_selected))
    if any(i < 0 or i >= len(e) for i in a | b):
        raise ValueError('cell index out of range')
    return float(sum(e[i] for i in a.symmetric_difference(b)))


def selection_jump(energies: Sequence[float], old_selected: Sequence[int], new_selected: Sequence[int]) -> float:
    e = np.asarray(energies, float)
    a, b = set(map(int, old_selected)), set(map(int, new_selected))
    return float(sum(e[i] for i in b) - sum(e[i] for i in a))


def service_no_escape(positive_work: float, negative_work: float, final_energy: float, relink_energy: float) -> dict[str, float | str]:
    """Finite-time selected-cell work identity after discarding nonnegative initial energy.

    For piecewise material selected sets,
      E_T = E_0 + P_+ - P_- + sum_j J_j,
    and |J_j| is bounded by the Moyal energy in the symmetric difference.
    Hence P_+ <= E_T + P_- + R_switch.
    """
    vals = [positive_work, negative_work, final_energy, relink_energy]
    if min(vals) < 0:
        raise ValueError('nonnegative ledger values required')
    rhs = negative_work + final_energy + relink_energy
    if positive_work > rhs + 1e-12 * max(1.0, positive_work, rhs):
        raise ValueError('inputs violate the coherent selected-cell energy balance')
    if positive_work == 0:
        return {'branch': 'zero_service', 'threshold': 0.0, 'margin': rhs}
    th = positive_work / 3.0
    candidates = {
        'terminal_coherent_energy': final_energy,
        'backflow_or_cancellation': negative_work,
        'relink_symmetric_difference': relink_energy,
    }
    branch, val = max(candidates.items(), key=lambda kv: kv[1])
    if val + 1e-14 < th:
        raise AssertionError('one-third coherent service no-escape failed')
    return {'branch': branch, 'threshold': th, 'branch_value': val, 'margin': val - th}


def affine_phase_point(L: np.ndarray, X: np.ndarray, k: np.ndarray) -> np.ndarray:
    L = np.asarray(L, float); X = np.asarray(X, float); k = np.asarray(k, float)
    return np.concatenate((0.5 * np.linalg.solve(L, X), L.T @ k))


def affine_phase_covariance_residual(M: np.ndarray, L: np.ndarray, X: np.ndarray, k: np.ndarray) -> float:
    """Common affine map L->ML, X->MX, k->M^-T k leaves intrinsic phase point fixed."""
    M = np.asarray(M, float); L = np.asarray(L, float)
    z0 = affine_phase_point(L, X, k)
    z1 = affine_phase_point(M @ L, M @ X, np.linalg.solve(M.T, k))
    return float(np.linalg.norm(z1 - z0))


def exact_certificate() -> dict[str, str]:
    return {
        'polarized_moyal': 'int V_g f conjugate(V_g F) dmu = <f,F> for ||g||_2=1',
        'cell_work': 'W_C=2 Re int_C V_g f conjugate(V_g F), sum_C W_C=2 Re<f,F>',
        'material_affine_cells': 'zeta=(L^-1 X/2,L^T k) invariant under L->ML,X->MX,k->M^-T k',
        'switch_cost': '|E(S_new)-E(S_old)| <= E(S_new symmetric_difference S_old)',
        'finite_time_no_escape': 'P_plus <= E_final + P_minus + R_switch; hence one branch >=P_plus/3',
        'status': 'EXACT_BY_POLARIZED_MOYAL_AND_POSITIVE_CELL_ENERGY',
    }


@dataclass(frozen=True)
class CoherentTransferStress:
    samples: int
    worst_cross_moyal_relative_residual: float
    worst_cell_work_relative_residual: float
    worst_affine_phase_residual: float
    worst_switch_ratio: float
    minimum_no_escape_margin: float


def stress(samples: int = 20_000, seed: int = 20260808) -> CoherentTransferStress:
    rng = np.random.default_rng(seed)
    wm = ww = wa = ws = 0.0
    mn = float('inf')
    for _ in range(samples):
        n = int(rng.integers(8, 40))
        f = rng.normal(size=n) + 1j * rng.normal(size=n)
        F = rng.normal(size=n) + 1j * rng.normal(size=n)
        g = rng.normal(size=n) + 1j * rng.normal(size=n)
        g /= np.linalg.norm(g)
        lhs = discrete_cross_moyal(f, F, g)
        rhs = np.vdot(F, f)  # sum f conjugate(F)
        scale = max(1.0, abs(rhs), np.linalg.norm(f) * np.linalg.norm(F))
        wm = max(wm, abs(lhs - rhs) / scale)
        if abs(lhs - rhs) > 3e-10 * scale:
            raise AssertionError('polarized Moyal regression failed')

        Vf = periodic_discrete_stft(f, g); VF = periodic_discrete_stft(F, g)
        labels = rng.integers(0, int(rng.integers(2, 12)), size=Vf.shape)
        work = cell_work(Vf, VF, labels)
        exact = 2.0 * float(np.vdot(F, f).real)
        wr = abs(float(work.sum()) - exact) / max(1.0, abs(exact), np.linalg.norm(f) * np.linalg.norm(F))
        ww = max(ww, wr)
        if wr > 3e-10:
            raise AssertionError('cell work partition failed')

        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        if np.linalg.det(Q) < 0: Q[:, 0] *= -1
        scales = np.exp(rng.uniform(-1.5, 1.5, size=3))
        M = Q @ np.diag(scales) @ Q.T
        L = np.linalg.qr(rng.normal(size=(3, 3)))[0] @ np.diag(np.exp(rng.uniform(-1,1,size=3)))
        X = rng.normal(size=3); k = rng.normal(size=3)
        ar = affine_phase_covariance_residual(M, L, X, k)
        wa = max(wa, ar)
        if ar > 5e-10 * max(1.0, np.linalg.norm(X), np.linalg.norm(k)):
            raise AssertionError('affine phase covariance failed')

        e = rng.lognormal(mean=-2, sigma=1, size=int(rng.integers(4, 30)))
        ids = np.arange(len(e)); rng.shuffle(ids)
        a = ids[:len(e)//3].tolist(); b = ids[len(e)//4:2*len(e)//3].tolist()
        jump = abs(selection_jump(e, a, b)); sym = symmetric_difference_energy(e, a, b)
        ws = max(ws, jump / max(sym, 1e-300))
        if jump > sym + 2e-13 * max(1.0, sym):
            raise AssertionError('selection jump exceeded symmetric-difference energy')

        Pm = float(rng.lognormal(-1, .8)); N = float(rng.lognormal(-1, .8)); R = float(rng.lognormal(-1, .8)); ET = float(rng.lognormal(-1, .8))
        Pp = float(rng.uniform(0.0, N + R + ET))
        route = service_no_escape(Pp, N, ET, R)
        mn = min(mn, float(route['margin']))
    return CoherentTransferStress(samples, wm, ww, wa, ws, mn)


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument('--samples', type=int, default=20_000); ap.add_argument('--outdir', type=Path, default=Path('results-coherent-transfer-cells'))
    args = ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    cert = exact_certificate(); out = stress(args.samples)
    (args.outdir/'coherent_transfer_cells.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)}, indent=2))
    md = f'''# Coherent transfer cells: exact nonlinear work and relinking balance

Status: **{cert['status']}**.

Polarized Moyal gives an exact cellwise work decomposition

`W_C = 2 Re int_C V_g f conjugate(V_g F) dmu`,  `sum_C W_C = 2 Re <f,F>`.

If the coherent cells are transported by the common affine phase map `L->ML, X->MX, k->M^-T k`, their intrinsic coordinate `zeta=(L^-1X/2,L^Tk)` is unchanged.  Thus common affine motion creates no coherent-cell interface forcing.

For a piecewise material selected family, switching from one selected cell set to another has jump bounded by the positive Moyal energy in the symmetric difference.  Integrating the exact cell work balance gives

`P_plus <= E_final + P_minus + R_switch`.

Hence positive nonlinear service cannot disappear: at least one of terminal coherent energy, backflow/cancellation, or relinking symmetric-difference energy is at least `P_plus/3`.

Stress: `{out.samples}`
- worst polarized-Moyal relative residual: `{out.worst_cross_moyal_relative_residual:.3e}`
- worst cell-work relative residual: `{out.worst_cell_work_relative_residual:.3e}`
- worst affine phase residual: `{out.worst_affine_phase_residual:.3e}`
- worst switch jump / symmetric-difference ratio: `{out.worst_switch_ratio:.9f}`
- minimum one-third routing margin: `{out.minimum_no_escape_margin:.3e}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__ == '__main__': main()
