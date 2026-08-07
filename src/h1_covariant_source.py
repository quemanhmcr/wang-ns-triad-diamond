from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.curvature_sideband_irrep import hook_component, random_divfree_curvature

KAPPA_MAX = 21.0 / 20.0
BASE_ACTION_MAX = 1.0 / 30.0
FRAME_RATE_COEFF = 15.0 / math.sqrt(2.0)  # real-frame Frobenius rate from child helical 15/2 bound


def hook_rate(Bdot: np.ndarray) -> np.ndarray:
    """Derivative of the fixed grain-coordinate hook projection."""
    return hook_component(np.asarray(Bdot, float))


def physical_hook_matrix(Bh_slice: np.ndarray, L: np.ndarray) -> np.ndarray:
    L = np.asarray(L, float)
    return L @ np.asarray(Bh_slice, float) @ np.linalg.inv(L)


def physical_hook_matrix_rate(
    A: np.ndarray, L: np.ndarray, Bh_slice: np.ndarray, Bh_dot_slice: np.ndarray
) -> np.ndarray:
    """Exact derivative of G=L Bh L^-1 under Ldot=A L."""
    A = np.asarray(A, float)
    G = physical_hook_matrix(Bh_slice, L)
    return A @ G + physical_hook_matrix(Bh_dot_slice, L) - G @ A


def per_role_covariant_D_bound(source_norm: float, strain_norm: float, curvature_norm: float, kappa: float = KAPPA_MAX) -> float:
    """Clean bound before three-role direct sum and base-propagator conditioning.

    C=Sym(L Bh L^-1).  Since Bh_dot=P_H Bdot and Bdot=-2 A_aff B+S,
      ||Bh_dot|| <= 2 kappa ||A|| ||B|| + ||S||.
    The physical matrix derivative contributes (2 kappa+2 kappa^2)||A||||B||+kappa||S||.
    Moving triad-normal frames contribute 2*(15/sqrt2)*kappa ||A||||B||.
    The interaction commutator [D0,D] contributes 2 kappa ||A||||B||.
    """
    if min(source_norm, strain_norm, curvature_norm) < 0 or kappa < 1:
        raise ValueError("invalid source bound data")
    conn = 4 * kappa + 2 * kappa * kappa + 15 * math.sqrt(2.0) * kappa
    return kappa * source_norm + conn * strain_norm * curvature_norm


def interaction_forcing_variation_density_upper(
    source_norm: float,
    strain_norm: float,
    curvature_norm: float,
    kappa: float = KAPPA_MAX,
    action_budget: float = BASE_ACTION_MAX,
) -> float:
    """Three-role interaction-picture H1 forcing derivative upper bound."""
    if action_budget < 0:
        raise ValueError("nonnegative action budget required")
    return math.exp(2 * action_budget) * math.sqrt(3.0) * per_role_covariant_D_bound(
        source_norm, strain_norm, curvature_norm, kappa
    )


def clean_interaction_variation_upper(source_l1: float, strain_curvature_l1: float) -> float:
    """Arb-certified clean integrated bound J1 <= 2 S1 + 54 AB1."""
    if source_l1 < 0 or strain_curvature_l1 < 0:
        raise ValueError("nonnegative ledgers required")
    return 2.0 * source_l1 + 54.0 * strain_curvature_l1


def dephasing_source_route(I1: float, T: float, J1: float, source_l1: float, strain_curvature_l1: float) -> dict[str, float | str]:
    """Route a physical H1 dephasing branch through the clean source bound.

    If J1>=I1/(11T) and J1<=2 S+54 AB, then either
      S>=I1/(44T) or AB>=I1/(1188T).
    """
    vals = [I1, J1, source_l1, strain_curvature_l1]
    if T <= 0 or any(v < 0 for v in vals):
        raise ValueError("invalid dephasing route data")
    lower = I1 / (11.0 * T)
    if J1 + 1e-12 * max(1.0, lower) < lower:
        raise ValueError("not on the H1 dephasing branch")
    upper = clean_interaction_variation_upper(source_l1, strain_curvature_l1)
    if J1 > upper + 1e-11 * max(1.0, J1):
        raise ValueError("data violate the H1 source upper bound")
    source_thresh = I1 / (44.0 * T)
    coupling_thresh = I1 / (1188.0 * T)
    if source_l1 >= source_thresh - 1e-13 * max(1.0, source_thresh):
        branch = "curvature_source"
    else:
        branch = "strain_curvature_coupling"
        if strain_curvature_l1 + 1e-12 * max(1.0, coupling_thresh) < coupling_thresh:
            raise AssertionError("source/coupling dichotomy failed")
    return {
        "branch": branch,
        "source_threshold": source_thresh,
        "strain_curvature_threshold": coupling_thresh,
        "per_source_channel_threshold": source_thresh / 3.0,
    }


def h1_dominant_strain_or_source(I1: float, IB: float, T: float, source_l1: float, sigma_max: float) -> dict[str, float | str]:
    """On I1>=IB/2, coupling can be bounded by sigma_max*IB.

    If sigma_max*T <1/2376, the coupling alternative AB>=I1/(1188T)
    is impossible, hence the curvature source must be >=I1/(44T).
    """
    if min(I1, IB, source_l1, sigma_max) < 0 or T <= 0:
        raise ValueError("invalid H1-dominant data")
    if I1 + 1e-12 < IB / 2:
        raise ValueError("H1-dominant branch requires I1>=IB/2")
    strain_action = sigma_max * T
    if strain_action >= 1 / 2376 - 1e-15:
        branch = "base_strain_action"
    else:
        branch = "curvature_source"
        if source_l1 + 1e-12 * max(1.0, I1 / (44 * T)) < I1 / (44 * T):
            raise ValueError("small-strain H1 dephasing requires the source lower bound")
    return {
        "branch": branch,
        "strain_action_threshold": 1 / 2376,
        "source_threshold": I1 / (44 * T),
        "per_source_channel_threshold": I1 / (132 * T),
    }


def arb_source_constants_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    k = arb(21) / 20
    K = arb(1) / 30
    src = (2 * K).exp() * arb(3).sqrt() * k
    conn0 = 4 * k + 2 * k * k + 15 * arb(2).sqrt() * k
    conn = (2 * K).exp() * arb(3).sqrt() * conn0
    if not (src < arb(2)):
        raise AssertionError(f"source coefficient not below 2: {src}")
    if not (conn < arb(54)):
        raise AssertionError(f"connection coefficient not below 54: {conn}")
    # Exact threshold arithmetic belongs to rationals, not Arb touching-ball comparisons.
    from fractions import Fraction
    if Fraction(54, 1188) != Fraction(1, 22):
        raise AssertionError("dephasing coupling split arithmetic failed")
    if Fraction(108, 2376) != Fraction(1, 22):
        raise AssertionError("H1-dominant strain threshold arithmetic failed")
    return {
        "mild_aspect_kappa": "21/20",
        "base_action_budget": "1/30",
        "raw_source_coefficient_ball": str(src),
        "clean_source_coefficient": "2",
        "raw_connection_coefficient_ball": str(conn),
        "clean_connection_coefficient": "54",
        "dephasing_source_threshold": "I_1/(44 T)",
        "dephasing_strain_curvature_threshold": "I_1/(1188 T)",
        "H1_dominant_strain_action_threshold": "1/2376",
        "per_NS_source_channel_threshold": "I_1/(132 T)",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class H1SourceStress:
    samples: int
    worst_hook_rate_projection_residual: float
    worst_physical_matrix_rate_residual: float
    worst_clean_density_margin: float
    minimum_route_margin: float


def random_mild_L(rng: np.random.Generator) -> np.ndarray:
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    vals = np.sort(rng.uniform(-0.5, 0.5, size=3))
    span = vals[-1] - vals[0]
    if span > 0:
        vals *= math.log(KAPPA_MAX) / span
    axes = np.exp(vals)
    axes /= math.sqrt(axes.max() * axes.min())
    return Q @ np.diag(axes) @ Q.T


def stress(samples: int = 50_000, seed: int = 20260807) -> H1SourceStress:
    rng = np.random.default_rng(seed)
    wh = wg = 0.0
    wmargin = float("inf")
    rmargin = float("inf")
    for _ in range(samples):
        B = random_divfree_curvature(rng)
        Bdot = random_divfree_curvature(rng)
        Bh = hook_component(B)
        Bhd = hook_rate(Bdot)
        # Exact linearity of the fixed hook projector.
        proj_res = hook_component(B + Bdot) - Bh - Bhd
        wh = max(wh, float(np.linalg.norm(proj_res)))
        A = rng.normal(size=(3, 3)); A -= np.trace(A) / 3 * np.eye(3)
        sigma = float(np.linalg.norm(A, 2))
        if sigma > 0:
            A /= sigma
            sigma = float(rng.uniform(0, 2.0))
            A *= sigma
        L = random_mild_L(rng)
        c = int(rng.integers(0, 3))
        G = physical_hook_matrix(Bh[:, :, c], L)
        Gd = physical_hook_matrix_rate(A, L, Bh[:, :, c], Bhd[:, :, c])
        # Independent exact product-rule expansion of d(L Bh L^-1)/dt.
        Li = np.linalg.inv(L)
        direct = (A @ L) @ Bh[:, :, c] @ Li + L @ Bhd[:, :, c] @ Li \
            - L @ Bh[:, :, c] @ Li @ (A @ L) @ Li
        wg = max(wg, float(np.linalg.norm(direct - Gd)))
        # Create a compatible curvature source S = Bdot + 2 A_aff B.
        Aa = np.linalg.inv(L) @ A @ L
        S = Bdot + 2 * np.einsum('ad,dbc->abc', Aa, B)
        source = float(np.linalg.norm(S)); curv = float(np.linalg.norm(B))
        raw = interaction_forcing_variation_density_upper(source, sigma, curv, float(np.linalg.cond(L)), BASE_ACTION_MAX)
        clean = 2 * source + 54 * sigma * curv
        wmargin = min(wmargin, clean - raw)
        if raw > clean + 2e-10 * max(1.0, clean):
            raise AssertionError("clean H1 variation density bound failed")
        # Synthetic exact route on the boundary algebra.
        I1 = float(10 ** rng.uniform(-4, -0.2)); T = float(10 ** rng.uniform(-3, 0))
        low = I1 / (11 * T)
        if rng.random() < 0.5:
            src = (1 + rng.random()) * I1 / (44 * T)
            ab = 0.0
        else:
            src = 0.99 * I1 / (44 * T)
            ab = (1 + rng.random()) * I1 / (1188 * T)
        J = min(2 * src + 54 * ab, (1 + rng.random()) * low)
        J = max(J, low)
        # Ensure upper bound if max pushed above it.
        if J > 2 * src + 54 * ab:
            J = 2 * src + 54 * ab
        if J + 1e-13 >= low:
            out = dephasing_source_route(I1, T, J, src, ab)
            margin = src - out['source_threshold'] if out['branch'] == 'curvature_source' else ab - out['strain_curvature_threshold']
            rmargin = min(rmargin, float(margin))
    if wh > 5e-11 or wg > 5e-11:
        raise AssertionError("H1 source exact-identity regression failed")
    return H1SourceStress(samples, wh, wg, wmargin, rmargin)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-h1-covariant-source"))
    args = ap.parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_source_constants_certificate(); out = stress(args.samples)
    (args.outdir / "h1_covariant_source.json").write_text(json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8")
    md = f"""# Physical H1 covariant source calculus

Status: **{cert['status']}**.

Let `B_H=P_H B` be the five-dimensional hook curvature.  The full affine curvature equation is `Bdot=-2 A_aff B+S`, hence the fixed hook projector gives the exact identity `B_Hdot=P_H(-2 A_aff B+S)`.  With `G_c=L B_H,c L^-1`, `Ldot=A L`, the physical hook matrix obeys `Gdot=A G+L B_Hdot L^-1-G A`.

Using `cond(L)<=21/20`, the good-core triad-normal frame rate, the objective base generator bound and low-strain action `K<=1/30`, the three-role interaction-picture H1 forcing satisfies the clean pointwise/integrated estimate

`J1 <= 2 int||S|| + 54 int ||A|| ||B||`.

Therefore the H1 dephasing branch `J1>=I1/(11T)` forces

`int||S|| >= I1/(44T)`

or

`int||A||||B|| >= I1/(1188T)`.

On the H1-dominant full-curvature branch `I1>=I_B/2`, if `||A||_infty T<1/2376`, the strain-curvature alternative is impossible and the curvature source is mandatory.  Since `S=S_P+S_R+S_nu` with `S_P` from pressure third derivatives, `S_R` from differentiated SGS stress and `S_nu` from viscous fourth derivatives, one source channel has integrated norm at least `I1/(132T)`.

This source theorem does not charge base strain: if the `1/2376` action threshold is crossed, it is handed to the existing objective-strain/source branch.  Pressure-third far-field locality retains the previous `6-3=3` summable exponent.

Stress: `{out.samples}`
- worst exact hook-projector linearity residual: `{out.worst_hook_rate_projection_residual:.3e}`
- worst exact physical-hook product-rule residual: `{out.worst_physical_matrix_rate_residual:.3e}`
- minimum clean-density margin: `{out.worst_clean_density_margin:.3e}`
- minimum routing margin: `{out.minimum_route_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8"); print(md)


if __name__ == "__main__":
    main()
