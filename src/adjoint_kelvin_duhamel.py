from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.linalg import expm

from src.affine_kelvin_packet_pde import objective_coordinate_generator, transverse_frame
from src.single_edge_certificate import RSTAR_LO, RSTAR_HI

INHERIT_FRACTION = Fraction(1, 4)
RESIDUAL_FRACTION = Fraction(1, 4)
GENERATE_FRACTION = Fraction(1, 2)


def forcing_propagator(G: np.ndarray, dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Return U=exp(G dt) and B=int_0^dt exp(G(dt-s)) ds exactly via block exponential."""
    G = np.asarray(G, complex)
    if G.ndim != 2 or G.shape[0] != G.shape[1] or dt < 0:
        raise ValueError('square generator and nonnegative dt required')
    n = G.shape[0]
    M = np.zeros((2*n, 2*n), complex)
    M[:n, :n] = G
    M[:n, n:] = np.eye(n)
    E = expm(M * float(dt))
    return E[:n, :n], E[:n, n:]


def adjoint_pairing_derivative_residual(
    G: np.ndarray,
    c: np.ndarray,
    psi: np.ndarray,
    forcing: np.ndarray,
) -> complex:
    """Exact algebra: d<psi,c>/dt=<psi,F> for c_dot=Gc+F, psi_dot=-G^*psi."""
    G = np.asarray(G, complex); c = np.asarray(c, complex); psi = np.asarray(psi, complex); F = np.asarray(forcing, complex)
    dc = G @ c + F
    dpsi = -G.conj().T @ psi
    lhs = np.vdot(dpsi, c) + np.vdot(psi, dc)
    rhs = np.vdot(psi, F)
    return lhs - rhs


def exact_piecewise_duhamel(
    c0: np.ndarray,
    psi_terminal: np.ndarray,
    generators: Sequence[np.ndarray],
    dts: Sequence[float],
    high_high_forcings: Sequence[np.ndarray],
    residual_forcings: Sequence[np.ndarray],
) -> dict[str, complex | np.ndarray | list[complex]]:
    """Exact piecewise-constant adjoint Duhamel ledger.

    Bulk Kelvin/strain/viscosity is inside G.  Only high-high generation and already
    classified residual forcing appear in the Duhamel impulses.
    """
    m = len(generators)
    if not (len(dts) == len(high_high_forcings) == len(residual_forcings) == m):
        raise ValueError('slab arrays must have equal length')
    c = np.asarray(c0, complex).copy()
    n = len(c)
    Us: list[np.ndarray] = []
    Bs: list[np.ndarray] = []
    states = [c.copy()]
    for G, dt, F, R in zip(generators, dts, high_high_forcings, residual_forcings):
        G = np.asarray(G, complex); F = np.asarray(F, complex); R = np.asarray(R, complex)
        if G.shape != (n, n) or F.shape != (n,) or R.shape != (n,):
            raise ValueError('dimension mismatch')
        U, B = forcing_propagator(G, float(dt))
        c = U @ c + B @ (F + R)
        Us.append(U); Bs.append(B); states.append(c.copy())
    psi = np.asarray(psi_terminal, complex).copy()
    if psi.shape != (n,): raise ValueError('terminal dual dimension mismatch')
    psis = [None] * (m + 1); psis[m] = psi.copy()
    hh_atoms = [0j] * m; rr_atoms = [0j] * m
    for i in range(m - 1, -1, -1):
        hh_atoms[i] = complex(np.vdot(psis[i+1], Bs[i] @ np.asarray(high_high_forcings[i], complex)))
        rr_atoms[i] = complex(np.vdot(psis[i+1], Bs[i] @ np.asarray(residual_forcings[i], complex)))
        psis[i] = Us[i].conj().T @ psis[i+1]
    z0 = complex(np.vdot(psis[0], states[0]))
    z1 = complex(np.vdot(psis[m], states[m]))
    ihh = sum(hh_atoms, 0j); ir = sum(rr_atoms, 0j)
    residual = z1 - z0 - ihh - ir
    return {
        'initial_pairing': z0,
        'terminal_pairing': z1,
        'high_high_impulse': ihh,
        'residual_impulse': ir,
        'duhamel_residual': residual,
        'high_high_slab_atoms': hh_atoms,
        'residual_slab_atoms': rr_atoms,
        'terminal_state': states[-1],
        'initial_dual': psis[0],
    }


def terminal_unit_dual(c_terminal: np.ndarray) -> np.ndarray:
    c = np.asarray(c_terminal, complex)
    n = float(np.linalg.norm(c))
    if n <= 0: raise ValueError('nonzero terminal coefficient required')
    return c / n


def inherit_or_generate_route(
    terminal_pairing: complex,
    initial_pairing: complex,
    residual_impulse: complex,
    high_high_impulse: complex | None = None,
) -> dict[str, float | str]:
    """A terminal coefficient is inherited, classified residual, or high-high generated.

    If |z1|=A and both |z0|<A/4 and |IR|<A/4, exact Duhamel forces |IHH|>A/2.
    """
    A = abs(terminal_pairing)
    if A <= 0: return {'branch': 'zero_terminal', 'terminal_amplitude': 0.0}
    if high_high_impulse is None:
        high_high_impulse = terminal_pairing - initial_pairing - residual_impulse
    if abs(initial_pairing) >= A * float(INHERIT_FRACTION):
        return {'branch': 'material_inheritance', 'terminal_amplitude': A, 'threshold': A/4, 'value': abs(initial_pairing)}
    if abs(residual_impulse) >= A * float(RESIDUAL_FRACTION):
        return {'branch': 'classified_residual', 'terminal_amplitude': A, 'threshold': A/4, 'value': abs(residual_impulse)}
    lower = A * float(GENERATE_FRACTION)
    val = abs(high_high_impulse)
    if val + 2e-13 * max(1.0, A) < lower:
        raise AssertionError('inherit-or-generate triangle inequality failed')
    return {'branch': 'high_high_generation', 'terminal_amplitude': A, 'threshold': lower, 'value': val}


def phase_aligned_positive_generation(atoms: Sequence[complex]) -> dict[str, object]:
    """Turn an integrated complex high-high impulse into a positive causal atom law.

    This is an amplitude-generation law.  It is NOT automatically the physical
    child-energy transfer law; the latter requires the existing near-extremal
    phase-lock / physical-weight bridge.
    """
    z = np.asarray(atoms, complex)
    total = complex(np.sum(z))
    A = abs(total)
    if A == 0:
        return {'total_impulse': total, 'positive_mass': 0.0, 'negative_mass': 0.0, 'weights': np.zeros(len(z)), 'aligned_real': np.zeros(len(z))}
    phase = total / A
    aligned = np.real(np.conj(phase) * z)
    pos = np.maximum(aligned, 0.0)
    neg = np.maximum(-aligned, 0.0)
    p = float(pos.sum()); n = float(neg.sum())
    if abs((p - n) - A) > 2e-12 * max(1.0, A, p+n):
        raise AssertionError('phase-aligned positive/negative decomposition failed')
    if p + 2e-13 * max(1.0, A) < A:
        raise AssertionError('positive generation mass below total impulse')
    weights = pos / p if p > 0 else np.zeros(len(z))
    return {'total_impulse': total, 'phase': phase, 'positive_mass': p, 'negative_mass': n, 'weights': weights, 'aligned_real': aligned}


def half_slab_generation(atoms: Sequence[complex], times: Sequence[float], t0: float, t1: float) -> dict[str, object]:
    """One half of a child slab carries at least half the positive aligned generation mass."""
    if not t1 > t0: raise ValueError('nontrivial time interval required')
    z = np.asarray(atoms, complex); t = np.asarray(times, float)
    if len(z) != len(t) or np.any(t < t0) or np.any(t > t1): raise ValueError('atoms/times outside slab')
    out = phase_aligned_positive_generation(z)
    aligned = np.asarray(out['aligned_real'], float); pos = np.maximum(aligned, 0.0)
    mid = 0.5 * (t0 + t1)
    left = float(pos[t <= mid].sum()); right = float(pos[t > mid].sum())
    mask = (t <= mid) if left >= right else (t > mid)
    chosen = left if left >= right else right
    total = float(pos.sum())
    if chosen + 2e-13 * max(1.0, total) < 0.5 * total:
        raise AssertionError('half-slab pigeonhole failed')
    return {'half': 'left' if left >= right else 'right', 'chosen_positive_mass': chosen, 'total_positive_mass': total, 'mask': mask, 'midpoint': mid}


def signed_good_lifetime_ratio_bounds() -> tuple[Fraction, Fraction]:
    """T_parent/T_child=(N_child/N_parent)^2 from 3/5<parent/child<5/8."""
    return Fraction(64, 25), Fraction(25, 9)


def common_parent_natural_window_lower(child_lifetime: float, interaction_cluster_width_fraction: float = 0.5) -> float:
    """Common overlap of parent *natural* backward windows for interactions in a short cluster.

    This is geometric only: it does not assert packet persistence on the entire window.
    """
    if child_lifetime <= 0 or not (0 <= interaction_cluster_width_fraction < float(Fraction(64,25))):
        raise ValueError('bad lifetime/window data')
    return (float(Fraction(64,25)) - interaction_cluster_width_fraction) * child_lifetime


def arb_causal_time_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:
        raise RuntimeError('python-flint required') from exc
    ctx.prec = 180
    def aq(q): return arb(q.numerator) / q.denominator
    r = aq(RSTAR_LO).union(aq(RSTAR_HI))
    lo = r * (-arb(1)/80).exp(); hi = r * (arb(1)/80).exp()
    if not (lo > arb(3)/5 and hi < arb(5)/8):
        raise AssertionError('signed-good parent/child frequency window failed')
    if Fraction(1,1) / Fraction(5,8) ** 2 != Fraction(64,25):
        raise AssertionError('lower lifetime ratio identity failed')
    if Fraction(1,1) / Fraction(3,5) ** 2 != Fraction(25,9):
        raise AssertionError('upper lifetime ratio identity failed')
    if Fraction(64,25) - Fraction(1,2) != Fraction(103,50):
        raise AssertionError('half-slab common-window identity failed')
    return {
        'signed_good_frequency': '3/5 < N_parent/N_child < 5/8',
        'natural_lifetime_ratio': '64/25 < T_parent/T_child < 25/9',
        'half_child_slab_parent_window_overlap': '>103/50 T_child',
        'scope': 'natural parabolic windows only; actual packet persistence is not claimed',
        'status': 'CERTIFIED_ADJOINT_DUHAMEL_AND_CAUSAL_TIME_GEOMETRY',
    }


@dataclass(frozen=True)
class AdjointDuhamelStress:
    samples: int
    worst_instantaneous_pairing_residual: float
    worst_piecewise_duhamel_residual: float
    minimum_generation_margin: float
    minimum_positive_mass_margin: float
    minimum_half_slab_margin: float
    branch_counts: dict[str, int]


def random_objective_generator(rng: np.random.Generator) -> np.ndarray:
    A = rng.normal(size=(3,3)); A -= np.trace(A)/3 * np.eye(3)
    k = rng.normal(size=3); E = transverse_frame(k)
    nu = float(rng.uniform(0,1.5)); k2 = float(np.dot(k,k))
    return objective_coordinate_generator(A, E, nu*k2).astype(complex)


def stress(samples: int = 50_000, seed: int = 20260808) -> AdjointDuhamelStress:
    rng = np.random.default_rng(seed)
    wi = wp = 0.0; mg = mp = mh = float('inf'); branches: dict[str,int] = {}
    # Full 50k: instantaneous adjoint algebra and cheap phase/time routing.
    for _ in range(samples):
        G = random_objective_generator(rng)
        c = rng.normal(size=2) + 1j*rng.normal(size=2)
        psi = rng.normal(size=2) + 1j*rng.normal(size=2)
        F = rng.normal(size=2) + 1j*rng.normal(size=2)
        rr = adjoint_pairing_derivative_residual(G,c,psi,F)
        wi = max(wi, abs(rr))
        if abs(rr) > 3e-12 * max(1.0, np.linalg.norm(c)*np.linalg.norm(psi), np.linalg.norm(F)*np.linalg.norm(psi)):
            raise AssertionError('instantaneous adjoint pairing cancellation failed')
        n=int(rng.integers(1,10)); atoms=rng.normal(size=n)+1j*rng.normal(size=n)
        phase=phase_aligned_positive_generation(atoms); mp=min(mp,float(phase['positive_mass'])-abs(complex(phase['total_impulse'])))
        times=np.sort(rng.uniform(0,1,size=n)); hs=half_slab_generation(atoms,times,0.0,1.0); mh=min(mh,float(hs['chosen_positive_mass'])-.5*float(hs['total_positive_mass']))
    # Exact block-exponential Duhamel histories on a representative 5k subset.
    for _ in range(min(samples, 5_000)):
        m = int(rng.integers(1,7)); Gs=[random_objective_generator(rng) for _ in range(m)]; dts=rng.uniform(.002,.08,size=m)
        HH=[rng.normal(size=2)+1j*rng.normal(size=2) for _ in range(m)]; R=[.1*(rng.normal(size=2)+1j*rng.normal(size=2)) for _ in range(m)]
        c0=rng.normal(size=2)+1j*rng.normal(size=2)
        trial=exact_piecewise_duhamel(c0,np.array([1.+0j,0j]),Gs,dts,HH,R)
        psi1=terminal_unit_dual(np.asarray(trial['terminal_state']))
        out=exact_piecewise_duhamel(c0,psi1,Gs,dts,HH,R)
        dr=abs(complex(out['duhamel_residual'])); wp=max(wp,dr)
        if dr>5e-11*max(1.0,abs(complex(out['terminal_pairing']))): raise AssertionError('piecewise exact Duhamel telescope failed')
        route=inherit_or_generate_route(complex(out['terminal_pairing']),complex(out['initial_pairing']),complex(out['residual_impulse']),complex(out['high_high_impulse']))
        b=str(route['branch']); branches[b]=branches.get(b,0)+1
        if b=='high_high_generation': mg=min(mg,float(route['value'])-float(route['threshold']))
    # Deterministic branch probes guarantee every clean triangle branch is exercised.
    probes=[
        inherit_or_generate_route(1+0j,.3+0j,0j,.7+0j),
        inherit_or_generate_route(1+0j,.1+0j,.3j,.9-.3j),
        inherit_or_generate_route(1+0j,.1+0j,.1j,.9-.1j),
    ]
    for route in probes:
        b=str(route['branch']); branches[b]=branches.get(b,0)+1
        if b=='high_high_generation': mg=min(mg,float(route['value'])-float(route['threshold']))
    if not math.isfinite(mg): mg=0.0
    return AdjointDuhamelStress(samples,wi,wp,mg,mp,mh,branches)


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-adjoint-kelvin-duhamel'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_causal_time_certificate(); out=stress(args.samples)
    data={'certificate':cert,'stress':asdict(out),'routing':{
        'terminal':'|z1|=A',
        'inherit':'|z0|>=A/4',
        'classified_residual':'|IR|>=A/4',
        'generate':'|IHH|>=A/2',
        'positive_generation':'phase-align total IHH; positive aligned atoms sum >=|IHH|',
        'important_scope':'positive Duhamel generation weights are amplitude-generation weights, not automatically physical child-energy transfer weights',
    }}
    (args.outdir/'adjoint_kelvin_duhamel.json').write_text(json.dumps(data,indent=2))
    md=f'''# Adjoint Kelvin--Duhamel causal gate

Status: **{cert['status']}**.

Put the common low--high affine/Kelvin transport **and bulk viscosity** into the linear objective generator `G(t)`.  Along a Kelvin characteristic let the selected transverse coefficient obey

`c_dot = G c + F_HH + R_class`,

where `F_HH` is actual high--high generation and `R_class` contains only already classified cross-cell / moving-projector / H1-H3 / window / profile residuals.  Pressure is absent because the role equation has already been Leray projected.

For the backward adjoint dual

`psi_dot = -G^* psi`,

one has exactly

`d <psi,c>/dt = <psi,F_HH+R_class>`.

Thus on a child slab

`z_1 = z_0 + I_HH + I_R`.

If the terminal dual is chosen along the terminal coefficient, `|z_1|=||c(t_1)||=:A`, and the exact triangle inequality gives the clean causal alternative

- inherited material coefficient `|z_0|>=A/4`; or
- classified residual `|I_R|>=A/4`; or
- genuine high--high generation `|I_HH|>=A/2`.

No common affine strain, pressure or bulk viscosity is counted again as generation.

For the generated branch, decompose the high--high Duhamel impulse into quadratic parent-pair atoms `z_alpha`.  A **single** phase aligned with the total impulse gives

`sum [Re(conj(phase) z_alpha)]_+ >= |I_HH|`.

After normalization this is a positive causal-generation law on same-time quadratic parent-pair events.  No pointwise persistence is required.  This is an amplitude-generation law; identification with the physical positive child-energy transfer law is made only on the already-certified near-extremal phase-locked core.

Signed-good scale geometry also gives

`64/25 < T_parent/T_child < 25/9`

for natural parabolic lifetimes.  A half-child-slab carrying at least half of the positive generation mass therefore has parent **natural** backward windows with common overlap longer than `103/50 T_child`.  This is geometric synchronization only; actual packet persistence on that whole common window is not claimed.

Stress: `{out.samples}` instantaneous/phase checks plus `min(samples,5000)` exact block-exponential Duhamel histories
- worst instantaneous adjoint-pairing residual: `{out.worst_instantaneous_pairing_residual:.3e}`
- worst exact piecewise Duhamel residual: `{out.worst_piecewise_duhamel_residual:.3e}`
- minimum high-high generation margin: `{out.minimum_generation_margin:.3e}`
- minimum positive aligned-mass margin: `{out.minimum_positive_mass_margin:.3e}`
- minimum half-slab pigeonhole margin: `{out.minimum_half_slab_margin:.3e}`
- branches: `{out.branch_counts}`
'''
    (args.outdir/'summary.md').write_text(md); print(md)

if __name__=='__main__': main()
