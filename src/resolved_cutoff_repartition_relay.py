from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.outer_moving_role_extraction import bilinear_apply, linearized_resolved

PARENT_RATIO_LO = 3.0 / 5.0
PARENT_RATIO_HI = 5.0 / 8.0
SMOOTH_ENVELOPE_LOWER = 11.0 / 20.0
LOW_STRAIN_ACTION = 1.0 / 30.0
TRANSPORT_CUT = 1.0 / 4.0


def cutoff_repartition_expression(
    tensor: np.ndarray,
    u: np.ndarray,
    V: np.ndarray,
    Q: np.ndarray,
) -> np.ndarray:
    """Nonlinear selected-role expression after moving L_V to the RHS.

    For w=Q u, h=u-V and L_V f=B(V,f)+B(f,V), define

      G_V = -L_V(Q u) + Q B(V,V) - Q B(h,h)
            + (L_V Q-Q L_V)u.

    Direct bilinear algebra gives G_V=-Q B(u,u), independently of V.  Thus a
    change of resolved cutoff is a repartition of the same NS interaction, not a
    new forcing term or a boundary currency.
    """
    T = np.asarray(tensor, complex)
    u = np.asarray(u, complex)
    V = np.asarray(V, complex)
    Q = np.asarray(Q, complex)
    n = len(u)
    if V.shape != (n,) or Q.shape != (n, n):
        raise ValueError("cutoff repartition dimensions do not match")
    h = u - V
    w = Q @ u
    LVw = linearized_resolved(T, V, w)
    LVu = linearized_resolved(T, V, u)
    heisenberg = linearized_resolved(T, V, Q @ u) - Q @ LVu
    return -LVw + Q @ bilinear_apply(T, V, V) - Q @ bilinear_apply(T, h, h) + heisenberg


def cutoff_repartition_residual(
    tensor: np.ndarray,
    u: np.ndarray,
    V_old: np.ndarray,
    V_new: np.ndarray,
    Q: np.ndarray,
) -> dict[str, float]:
    old = cutoff_repartition_expression(tensor, u, V_old, Q)
    new = cutoff_repartition_expression(tensor, u, V_new, Q)
    exact = -np.asarray(Q, complex) @ bilinear_apply(np.asarray(tensor, complex), np.asarray(u, complex), np.asarray(u, complex))
    return {
        "old_to_new": float(np.linalg.norm(old - new)),
        "old_to_full_ns": float(np.linalg.norm(old - exact)),
        "new_to_full_ns": float(np.linalg.norm(new - exact)),
    }


def renewed_parent_role_lower_relative_to_parent(
    parent_child_ratio: float,
    old_slab_strain: float = LOW_STRAIN_ACTION,
    renewed_slab_strain: float = LOW_STRAIN_ACTION,
) -> float:
    """Worst lower Fourier edge during the renewed parent-scale slab.

    The event envelope starts at 11N/20.  Backward transport to the common slice
    costs exp(-K_old).  Re-anchor the same carrier with parent scale Np=rN and
    transport through the renewed slab, costing exp(-K_new).  Hence

      |xi|/Np >= (11/20) exp(-(K_old+K_new)) / r.
    """
    r = float(parent_child_ratio)
    ko = float(old_slab_strain)
    kn = float(renewed_slab_strain)
    if r <= 0 or min(ko, kn) < 0 or not math.isfinite(r + ko + kn):
        raise ValueError("positive finite scale ratio and nonnegative strain actions required")
    return SMOOTH_ENVELOPE_LOWER * math.exp(-(ko + kn)) / r


def renewed_low_low_gap(
    parent_child_ratio: float = PARENT_RATIO_HI,
    old_slab_strain: float = LOW_STRAIN_ACTION,
    renewed_slab_strain: float = LOW_STRAIN_ACTION,
    transport_cut: float = TRANSPORT_CUT,
) -> float:
    return renewed_parent_role_lower_relative_to_parent(parent_child_ratio, old_slab_strain, renewed_slab_strain) - 2.0 * float(transport_cut)


def renewed_lifetime_ratio(parent_child_ratio: float) -> float:
    """Parabolic T~N^-2 lifetime ratio T_parent/T_child."""
    r = float(parent_child_ratio)
    if r <= 0 or not math.isfinite(r):
        raise ValueError("positive finite parent/child scale ratio required")
    return 1.0 / (r * r)


def arb_scale_renewal_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 180
    lower = (arb(11) / 20) * (-arb(1) / 15).exp() / (arb(5) / 8)
    if not (lower > arb(1) / 2):
        raise AssertionError(f"renewed parent low-low gap failed: {lower}")
    tmin = Fraction(1, 1) / Fraction(5, 8) ** 2
    tmax = Fraction(1, 1) / Fraction(3, 5) ** 2
    if tmin != Fraction(64, 25) or tmax != Fraction(25, 9):
        raise AssertionError("parent lifetime window changed")
    return {
        "renewed_relative_lower": str(lower),
        "renewed_low_low_target": "1/2",
        "parent_child_scale": "3/5 < Np/N < 5/8",
        "parent_child_lifetime": "64/25 < Tp/T < 25/9",
        "status": "CERTIFIED_SCALE_RENEWAL_SUPPORT_AND_PARABOLIC_LIFETIME_WINDOW",
    }


def theorem_certificate() -> dict[str, object]:
    gap = renewed_low_low_gap()
    return {
        "status": "EXACT_RESOLVED_CUTOFF_REPARTITION_GAUGE_AND_PARENT_SCALE_RENEWAL__GENERATED_SURVIVOR_REENTRY_NO_NEW_INTERFACE",
        "cutoff_identity": "-L_V(Q u)+Q B(V,V)-Q B(u-V,u-V)+(L_V Q-Q L_V)u = -Q B(u,u), independently of the resolved cutoff V",
        "meaning": "changing S_(N/4)u to S_(Np/4)u redistributes the same interaction between transport, HH source and Heisenberg interface; it creates no boundary forcing/currency",
        "renewed_support": f"from envelope 11N/20, two low-strain transports and Np/N<=5/8 give |xi|/Np >= (22/25)e^(-1/15)>1/2; clean gap {gap:.12g}",
        "renewed_low_low": "with Vp=S_(Np/4)u, B(Vp,Vp) is supported below Np/2 and is still excluded by the relayed carrier",
        "lifetime": "3/5<Np/N<5/8 implies 64/25<Tp/T<25/9 for T~N^-2; the existing asynchronous parabolic geometry is unchanged",
        "relay": "the common-slice material carrier may therefore be re-anchored at parent scale and the resolved cutoff changed without manufacturing a hard packet or a cutoff-switch stop",
        "scope": "this closes the algebra/support part of generated-survivor scale renewal; actual renewed slabs must still route failures of efficiency/moat/service/strain/material hypotheses through the existing physical first-stop system",
    }


@dataclass(frozen=True)
class CutoffRelayStress:
    samples: int
    worst_old_new_residual: float
    worst_full_ns_residual: float
    minimum_renewed_low_low_gap: float
    minimum_lifetime_margin: float


def stress(samples: int = 50_000, seed: int = 20260809) -> CutoffRelayStress:
    rng = np.random.default_rng(seed)
    won = wfull = 0.0
    mgap = mlife = float("inf")
    for _ in range(samples):
        n = int(rng.integers(2, 9))
        T = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
        u = rng.normal(size=n) + 1j * rng.normal(size=n)
        Vo = rng.normal(size=n) + 1j * rng.normal(size=n)
        Vn = rng.normal(size=n) + 1j * rng.normal(size=n)
        Q = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        rr = cutoff_repartition_residual(T, u, Vo, Vn, Q)
        scale = max(1.0, np.linalg.norm(T) * np.linalg.norm(u) ** 2 * np.linalg.norm(Q))
        won = max(won, rr["old_to_new"] / scale)
        wfull = max(wfull, rr["old_to_full_ns"] / scale, rr["new_to_full_ns"] / scale)
        if rr["old_to_new"] > 8e-11 * scale or max(rr["old_to_full_ns"], rr["new_to_full_ns"]) > 8e-11 * scale:
            raise AssertionError("resolved cutoff created a spurious nonlinear term")

        r = float(rng.uniform(PARENT_RATIO_LO, PARENT_RATIO_HI))
        ko = float(rng.uniform(0.0, LOW_STRAIN_ACTION))
        kn = float(rng.uniform(0.0, LOW_STRAIN_ACTION))
        gap = renewed_low_low_gap(r, ko, kn)
        mgap = min(mgap, gap)
        if gap <= 0:
            raise AssertionError("renewed parent carrier lost low-low exclusion")
        tr = renewed_lifetime_ratio(r)
        margin = min(tr - 64.0 / 25.0, 25.0 / 9.0 - tr)
        mlife = min(mlife, margin)
        if margin < -2e-14:
            raise AssertionError("renewed parabolic lifetime left signed-good window")
    return CutoffRelayStress(samples, won, wfull, mgap, mlife)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-cutoff-repartition-relay"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    arb = arb_scale_renewal_certificate()
    out = stress(args.samples)
    (args.outdir / "resolved_cutoff_repartition_relay.json").write_text(
        json.dumps({"certificate": cert, "arb": arb, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Resolved-cutoff repartition and parent-scale renewal\n\nStatus: **{cert['status']}**.\n\nLet `w=Q u`, choose any resolved field `V`, put `h=u-V`, and write `L_V f=B(V,f)+B(f,V)`.  Move the resolved transporter to the nonlinear side of the exact outer-role equation.  The entire cutoff-dependent nonlinear expression is\n\n`G_V=-L_V(Q u)+Q B(V,V)-Q B(h,h)+(L_V Q-Q L_V)u`.\n\nBilinearity gives the exact cancellation\n\n`G_V = Q[B(V,V)-B(h,h)-L_V u] = -Q B(u,u)`.\n\nTherefore `G_V` is **independent of V**.  Replacing `V=S_(N/4)u` by `V_p=S_(N_p/4)u` at a common-slice relay does not create a new physical source.  It only repartitions the same Navier--Stokes interaction between resolved transport, HH source and Heisenberg interface.  No cutoff-switch currency is permitted or needed.\n\nThe support geometry also renews with room to spare.  The smooth parent envelope starts at `11N/20`.  Transport to the common slice costs at most `e^(-1/30)`, and a renewed parent slab costs another `e^(-1/30)`.  Since a signed-good parent has `N_p/N<5/8`, throughout the renewed slab\n\n`|xi|/N_p >= (11/20)e^(-1/15)/(5/8) = (22/25)e^(-1/15) > 1/2`.\n\nArb certifies this strict inequality.  Meanwhile `V_p tensor V_p` is supported below `N_p/2`, so the exact low-low exclusion survives the scale change.  The parabolic lifetime also renews without modifying constants:\n\n`64/25 < T_p/T < 25/9`.\n\nThus a registered smooth material carrier can be re-anchored at the parent scale, switch to the parent resolved cutoff, and enter the next local role equation without a hard reselection and without a cutoff interface charge.\n\nStress: `{out.samples}` arbitrary bilinear/cutoff/scale states\n- worst old/new cutoff identity residual: `{out.worst_old_new_residual:.3e}`\n- worst residual against full `-Q B(u,u)`: `{out.worst_full_ns_residual:.3e}`\n- minimum renewed low-low support gap: `{out.minimum_renewed_low_low_gap:.6e}`\n- minimum sampled lifetime-window margin: `{out.minimum_lifetime_margin:.3e}`\n\nThis closes the **algebraic and support-geometric scale renewal of a generated survivor**.  What remains is not a cutoff problem: on each renewed physical slab, any failure of efficient transfer, moat, service, strain/coherent deformation or material coherence must be shown to enter the already named first-stop causes, and the non-generated critical/source/relink recursive routes must be attached to the same renewal mechanism.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
