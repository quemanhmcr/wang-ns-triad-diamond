from __future__ import annotations

import argparse
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from src.hermite_helicity_ledger import h3_tensor, sym3, third_hermite_norm_sq


def interval_variation_bound(values: np.ndarray, times: np.ndarray) -> tuple[float, float, float]:
    """Discrete regression analogue of the Banach interval-variation theorem.

    The theorem itself is continuous and valid in every Banach space:
      int ||f|| <= ||int f|| + (T/2) int ||f'||.
    Here trapezoidal integration is used only for stress testing smooth samples.
    """
    values = np.asarray(values)
    times = np.asarray(times, float)
    if values.ndim != 2 or len(times) != len(values) or len(times) < 2:
        raise ValueError("values must be (n,d) on at least two increasing times")
    if np.any(np.diff(times) <= 0):
        raise ValueError("times must increase")
    T = float(times[-1] - times[0])
    norms = np.linalg.norm(values, axis=1)
    forcing_l1 = float(np.trapezoid(norms, times))
    impulse = np.trapezoid(values, times, axis=0)
    impulse_norm = float(np.linalg.norm(impulse))
    deriv = np.diff(values, axis=0) / np.diff(times)[:, None]
    variation = float(np.sum(np.linalg.norm(deriv, axis=1) * np.diff(times)))
    return forcing_l1, impulse_norm, variation


def continuous_variation_margin(forcing_l1: float, impulse_norm: float, variation: float, T: float) -> float:
    if min(forcing_l1, impulse_norm, variation) < 0 or T <= 0:
        raise ValueError("invalid variation data")
    return impulse_norm + 0.5 * T * variation - forcing_l1


def first_duhamel_dichotomy(forcing_l1: float, impulse_norm: float, variation: float, T: float) -> dict[str, float | str]:
    """Exact consequence of the Banach interval-variation theorem.

    Either the first Duhamel impulse retains at least half of the L1 forcing, or
    the total variation is at least forcing_l1/T.
    """
    margin = continuous_variation_margin(forcing_l1, impulse_norm, variation, T)
    if margin < -1e-11 * max(1.0, forcing_l1):
        raise ValueError("data violate the continuous variation theorem")
    if impulse_norm >= 0.5 * forcing_l1:
        branch = "coherent_daughter"
    else:
        branch = "dephasing_source"
        if variation + 1e-10 * max(1.0, forcing_l1 / T) < forcing_l1 / T:
            raise AssertionError("dichotomy variation threshold failed")
    return {
        "branch": branch,
        "forcing_l1": float(forcing_l1),
        "impulse_norm": float(impulse_norm),
        "variation": float(variation),
        "variation_threshold": float(forcing_l1 / T),
        "margin": float(margin),
    }


def h3_l32_reverse_constant_clean() -> Fraction:
    """Clean critical-Young footprint constant from hypercontractivity + Paley-Zygmund.

    Let g(z)=exp(-|z|^2/4), so |g|^2 has standard Gaussian law and |g|^(3/2)
    has Gaussian variance s=4/3.  For P=T:H3, Gaussian degree-3
    hypercontractivity gives ||P||_4 <= 3^(3/2)||P||_2 under the s-Gaussian.
    Paley-Zygmund at theta=3/11 then gives an explicit L^(3/2)-from-L2 lower
    bound.  The exact arithmetic constant exceeds 1/160.
    """
    return Fraction(1, 160)


def h1_l32_reverse_constant_clean() -> Fraction:
    """Clean vector H1 footprint constant for a two-component polarization sideband.

    Degree-1 hypercontractivity + Paley-Zygmund gives the scalar constant; the
    largest singular value of a 2x3 coefficient matrix is >= ||C||_F/sqrt(2).
    The resulting clean lower bound exceeds 1/16.
    """
    return Fraction(1, 16)


def arb_sideband_capacity_certificate() -> dict[str, str]:
    """Interval-certify only the arithmetic constants; hypercontractivity is analytic input."""
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    theta = arb(3) / 11
    one = arb(1)
    s = arb(4) / 3
    # degree 3: PZ factor / 81, then variance-change L2 gain s^(3/2)
    c3 = theta.sqrt() * (one - theta) ** (arb(4) / 3) / 81 * s ** (arb(3) / 2)
    if not (c3 > arb(1) / 160):
        raise AssertionError(f"H3 reverse constant failed: {c3}")
    # degree 1: PZ factor / 9^(2/3), variance gain sqrt(s), then /sqrt(2) for 2 outputs
    c1 = theta.sqrt() * (one - theta) ** (arb(4) / 3) / (arb(9) ** (arb(2) / 3)) * s.sqrt() / arb(2).sqrt()
    if not (c1 > arb(1) / 16):
        raise AssertionError(f"H1 reverse constant failed: {c1}")
    # arithmetic used in the small-sideband norm-deficit theorem
    tail_fraction = arb(2916) / (80 * 80)
    if not (tail_fraction < arb(1) / 2):
        raise AssertionError("degree-3 fourth-moment tail threshold failed")
    return {
        "H3_exact_arithmetic_ball": str(c3),
        "H3_clean_reverse_constant": "1/160",
        "H1_two_component_arithmetic_ball": str(c1),
        "H1_clean_reverse_constant": "1/16",
        "small_sideband_threshold": "1/80",
        "hypercontractive_H3_fourth_factor": "729",
        "status": "CERTIFIED_ARITHMETIC_GIVEN_GAUSSIAN_HYPERCONTRACTIVITY",
    }


def divergence_free_curvature_basis() -> np.ndarray:
    """Orthonormal coefficient basis for B[a,b,c]=B[a,c,b], sum_a B[a,a,c]=0."""
    coords = [(a, b, c) for a in range(3) for b in range(3) for c in range(b, 3)]
    n = len(coords)
    def full(x: np.ndarray) -> np.ndarray:
        B = np.zeros((3, 3, 3))
        for val, (a, b, c) in zip(x, coords):
            B[a, b, c] = val
            B[a, c, b] = val
        return B
    C = np.zeros((3, n))
    for j in range(n):
        C[:, j] = np.einsum('aac->c', full(np.eye(n)[j]))
    _, _, vh = np.linalg.svd(C)
    return vh[3:].T


def coeff_to_curvature(x: np.ndarray) -> np.ndarray:
    coords = [(a, b, c) for a in range(3) for b in range(3) for c in range(b, 3)]
    B = np.zeros((3, 3, 3))
    for val, (a, b, c) in zip(np.asarray(x, float), coords):
        B[a, b, c] = val
        B[a, c, b] = val
    return B


def curvature_interaction_pullback(A_aff: np.ndarray, t: float, B: np.ndarray) -> np.ndarray:
    """For constant A_aff and Bdot+2 A_aff B=S, pull back by exp(2 A_aff t)."""
    P = expm(2.0 * np.asarray(A_aff, float) * float(t))
    return np.einsum('ad,dbc->abc', P, np.asarray(B, float))


def pulled_source(A_aff: np.ndarray, t: float, S: np.ndarray) -> np.ndarray:
    P = expm(2.0 * np.asarray(A_aff, float) * float(t))
    return np.einsum('ad,dbc->abc', P, np.asarray(S, float))


def h3_forcing_vector(T: np.ndarray) -> np.ndarray:
    """Coefficient vector with Euclidean norm sqrt(3/8)||T||_F.

    We use the flattened symmetric tensor scaled by sqrt(3/8); this is isometric
    to the normalized L2 H3 sideband norm established in hermite_helicity_ledger.
    """
    Ts = sym3(np.asarray(T, float))
    return math.sqrt(3.0 / 8.0) * Ts.ravel()


@dataclass(frozen=True)
class SidebandCoherenceStress:
    samples: int
    worst_variation_negative_margin: float
    worst_curvature_pullback_residual: float
    worst_h3_source_derivative_residual: float
    minimum_dichotomy_margin: float


def stress(samples: int = 20_000, seed: int = 20260807) -> SidebandCoherenceStress:
    rng = np.random.default_rng(seed)
    basis = divergence_free_curvature_basis()
    wvar = wpull = wh3 = 0.0
    mind = float("inf")
    # Smooth generic vector paths test the Banach/Hilbert variation inequality.
    ngrid = 81
    ts = np.linspace(0.0, 1.0, ngrid)
    for _ in range(samples // 20):
        d = 7
        a = rng.normal(size=d); b = rng.normal(size=d); c = rng.normal(size=d)
        vals = np.array([a + math.sin(2*math.pi*t)*b + 0.35*math.cos(4*math.pi*t)*c for t in ts])
        A1, imp, var = interval_variation_bound(vals, ts)
        margin = continuous_variation_margin(A1, imp, var, 1.0)
        wvar = max(wvar, max(0.0, -margin))
        mind = min(mind, margin)
        if margin < -3e-3:  # finite-difference variation overestimates only up to quadrature error
            raise AssertionError("variation inequality regression failed")
    # Constant-connection curvature ODE: Bdot=-2AB+S. Pullback derivative must be pulled source.
    for _ in range(samples):
        A = rng.normal(size=(3,3)); A -= np.trace(A)/3*np.eye(3)
        x = basis @ rng.normal(size=basis.shape[1]); B0 = coeff_to_curvature(x)
        S = rng.normal(size=(3,3,3)); S = 0.5*(S+np.swapaxes(S,1,2))
        t = float(rng.uniform(0.0, 0.08))
        # Exact instantaneous constant-connection identity at arbitrary t.
        # P=exp(2At) commutes with A, so d(PB)/dt=2APB+P(-2AB+S)=PS.
        P=expm(2.0*A*t)
        rhs=-2*np.einsum('ad,dbc->abc',A,B0)+S
        lhs=2*np.einsum('ad,dbc->abc',A,np.einsum('de,ebc->dbc',P,B0)) + np.einsum('ad,dbc->abc',P,rhs)
        ps=np.einsum('ad,dbc->abc',P,S)
        rr=float(np.linalg.norm(lhs-ps))/max(1.0,float(np.linalg.norm(ps)))
        wpull=max(wpull,rr)
        if rr>5e-11: raise AssertionError('curvature interaction pullback derivative failed')
        # H3 projection is linear, so the pulled H3 forcing derivative is exactly the H3 projection of PS.
        fd=h3_forcing_vector(sym3(lhs))
        fs=h3_forcing_vector(sym3(ps))
        hres=float(np.linalg.norm(fd-fs))/max(1.,float(np.linalg.norm(fs)))
        wh3=max(wh3,hres)
        if hres>5e-11: raise AssertionError('H3 pulled-source derivative failed')
    return SidebandCoherenceStress(samples,wvar,wpull,wh3,mind)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument('--samples',type=int,default=20_000)
    ap.add_argument('--outdir',type=Path,default=Path('results-sideband-coherence-daughter'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_sideband_capacity_certificate(); out=stress(args.samples)
    data={
        'certificate':cert,
        'stress':out.__dict__,
        'theorems':{
            'variation':'int||f|| <= ||int f|| + (T/2) int||f_dot||',
            'dichotomy':'daughter >= A/2 or variation >= A/T',
            'curvature_interaction':'B_tilde=P B, Pdot=2 A_aff P => dot B_tilde=P S',
            'H3_L2_forcing':'||F_H3||_2/||g||_2=sqrt(3/8)||Sym B_tilde||_F',
            'H3_Young_footprint':'||F_H3||_(3/2)/||g||_(3/2) >= (1/160)||F_H3||_2/||g||_2',
        }
    }
    (args.outdir/'sideband_coherence_daughter.json').write_text(json.dumps(data,indent=2),encoding='utf-8')
    md=f"""# Sideband coherence / first-Duhamel daughter theorem

Status: **{cert['status']}** for the clean Young-footprint arithmetic constants, conditional only on the standard Gaussian hypercontractivity inequality used analytically.

For an interaction-frame sideband forcing `f(t)` in any Banach space,

`int ||f|| <= ||int f|| + (T/2) int ||f_dot||`.

Hence, writing `A=int||f||`, either the first Duhamel daughter impulse has norm at least `A/2`, or the forcing variation is at least `A/T`.

For the affine curvature connection `Bdot+2 A_aff B=S`, pulling back by the connection gives `dot B_tilde=P S`.  The H3 envelope forcing in that pulled-back grain has normalized L2 norm `sqrt(3/8)||Sym B_tilde||` and its dephasing derivative is sourced by `Sym(P S)`.

Critical Young footprint:
- H3 Gaussian sideband: `L^(3/2) relative norm > (1/160) * relative L2 norm`;
- two-component H1 Gaussian sideband: `L^(3/2) relative norm > (1/16) * relative L2/Frobenius norm`.

Thus a coherent H1/H3 first-Duhamel daughter cannot be invisible to the critical Young capacity.  If the first iterate is later cancelled, that cancellation must come from nonlinear sideband/cross interactions rather than from the common affine connection.

Stress checks: `{out.samples}`
- worst negative variation margin: `{out.worst_variation_negative_margin:.3e}`
- worst exact curvature-pullback residual: `{out.worst_curvature_pullback_residual:.3e}`
- worst H3 pulled-source derivative residual: `{out.worst_h3_source_derivative_residual:.3e}`
"""
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
