from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.resolved_role_egorov import exact_affine_subtracted_commutator
from src.smooth_sgs_packet_equation import moving_multiplier_residual, transported_gaussian_symbol

TRANSPORT_CUT = 1.0 / 4.0
ROLE_LOWER = 3.0 / 5.0
LOW_STRAIN_ACTION = 1.0 / 30.0


def transported_role_lower_radius(
    role_lower: float = ROLE_LOWER,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Worst radial lower edge after trace-free affine transport.

    If Fdot=A F and K=int ||sym A||_op, singular values of F and F^{-1}
    lie in [exp(-K),exp(K)].  A dual-transported role which at its anchor
    time is supported in |xi|>=r N therefore stays in

        |xi| >= r exp(-K) N.
    """
    r = float(role_lower)
    K = float(strain_action)
    if r <= 0 or K < 0 or not math.isfinite(r + K):
        raise ValueError("positive finite role lower edge and nonnegative strain action required")
    return r * math.exp(-K)


def low_low_output_upper(transport_cut: float = TRANSPORT_CUT) -> float:
    """Minkowski support radius of V tensor V in units of N."""
    c = float(transport_cut)
    if c < 0 or not math.isfinite(c):
        raise ValueError("finite nonnegative transport cut required")
    return 2.0 * c


def persistent_low_low_gap(
    role_lower: float = ROLE_LOWER,
    strain_action: float = LOW_STRAIN_ACTION,
    transport_cut: float = TRANSPORT_CUT,
) -> float:
    """Strict support gap between the moving selected role and V tensor V."""
    return transported_role_lower_radius(role_lower, strain_action) - low_low_output_upper(transport_cut)


def arb_persistent_support_certificate() -> dict[str, str]:
    """Arb certificate for (3/5)exp(-1/30)>1/2."""
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 180
    moving_lower = (arb(3) / 5) * (-arb(1) / 30).exp()
    low_low = arb(1) / 2
    if not (moving_lower > low_low):
        raise AssertionError(f"moving-role low-low separation failed: {moving_lower}")
    return {
        "anchored_role_lower": "3/5",
        "strain_action_upper": "1/30",
        "transport_cut": "1/4",
        "low_low_output_upper": "1/2",
        "moving_role_lower": str(moving_lower),
        "status": "CERTIFIED_PERSISTENT_LOWLOW_EXCLUSION",
    }


def bilinear_apply(tensor: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Finite-dimensional bilinear model B(a,b), used only for exact algebra stress."""
    T = np.asarray(tensor, complex)
    a = np.asarray(a, complex)
    b = np.asarray(b, complex)
    if T.ndim != 3 or T.shape[0] != T.shape[1] or T.shape[1] != T.shape[2]:
        raise ValueError("cubic output/input tensor required")
    n = T.shape[0]
    if a.shape != (n,) or b.shape != (n,):
        raise ValueError("bilinear input dimension mismatch")
    return np.einsum("ijk,j,k->i", T, a, b)


def linearized_resolved(tensor: np.ndarray, V: np.ndarray, f: np.ndarray) -> np.ndarray:
    """L_V f=B(V,f)+B(f,V), the algebraic resolved low-high transporter."""
    return bilinear_apply(tensor, V, f) + bilinear_apply(tensor, f, V)


def outer_role_identity_residual(
    *,
    tensor: np.ndarray,
    u: np.ndarray,
    V: np.ndarray,
    Q: np.ndarray,
    dQ: np.ndarray,
    viscosity_operator: np.ndarray | None = None,
    viscosity: float = 0.0,
) -> np.ndarray:
    r"""Residual of the exact moving outer-role identity.

    Abstract Navier--Stokes algebra:

      u_dot = -B(u,u) + nu D u,
      L_V f = B(V,f)+B(f,V),
      h = u-V,
      w = Q u.

    Direct differentiation gives exactly

      (d_t+L_V-nu D)w
       = Q B(V,V) - Q B(h,h)
         + (dQ + L_V Q - Q L_V)u
         + nu (Q D-D Q)u.

    For scalar Fourier Q, D=Delta commutes with Q.  On the persistent support
    moat Q B(V,V)=0, leaving one genuine quadratic h-h source and one
    Heisenberg interface defect.
    """
    u = np.asarray(u, complex)
    V = np.asarray(V, complex)
    Q = np.asarray(Q, complex)
    dQ = np.asarray(dQ, complex)
    n = len(u)
    if V.shape != (n,) or Q.shape != (n, n) or dQ.shape != (n, n):
        raise ValueError("outer-role dimensions do not match")
    if viscosity_operator is None:
        D = np.zeros((n, n), complex)
    else:
        D = np.asarray(viscosity_operator, complex)
        if D.shape != (n, n):
            raise ValueError("viscosity operator dimension mismatch")
    nu = float(viscosity)
    if nu < 0:
        raise ValueError("nonnegative viscosity required")

    h = u - V
    du = -bilinear_apply(tensor, u, u) + nu * (D @ u)
    w = Q @ u
    dw = dQ @ u + Q @ du
    lhs = dw + linearized_resolved(tensor, V, w) - nu * (D @ w)

    LVu = linearized_resolved(tensor, V, u)
    comm_u = dQ @ u + linearized_resolved(tensor, V, Q @ u) - Q @ LVu
    visc_comm_u = nu * ((Q @ D - D @ Q) @ u)
    rhs = (
        Q @ bilinear_apply(tensor, V, V)
        - Q @ bilinear_apply(tensor, h, h)
        + comm_u
        + visc_comm_u
    )
    return lhs - rhs


def selected_outer_role_sources() -> tuple[str, ...]:
    """Only terms left after exact low-low exclusion on the selected role."""
    return (
        "physical_high_high_quadratic_source_-Q_P_div_h_tensor_h",
        "nonaffine_heisenberg_role_interface_(dtQ+[L_V,Q])u",
        "bulk_viscosity_inside_propagator",
    )


def theorem_certificate() -> dict[str, object]:
    gap = persistent_low_low_gap()
    return {
        "status": "EXACT_OUTER_MOVING_ROLE_IDENTITY_AND_PERSISTENT_LOWLOW_EXCLUSION__INTERFACE_PROVENANCE_SUPPLIED",
        "identity": "(dt+L_V-nuDelta)(Q u)=Q B(V,V)-Q B(h,h)+(dtQ+[L_V,Q])u, with h=u-V and [Q,Delta]=0 for scalar Fourier roles",
        "support": f"on K<=1/30 an anchor |xi|>=3N/5 stays above (3/5)e^(-1/30)N>N/2; clean dimensionless gap {gap:.12g}",
        "low_low": "supp Vhat subset B_(N/4) implies supp B(V,V) subset B_(N/2), hence Q B(V,V)=0 for the entire low-strain role interval",
        "high_high": "the sole nonlinear source after low-low exclusion is -Q P div(h tensor h); its coherent Hahn disintegration is the existing physical binary causal work law",
        "affine_gauge": "choose the scalar role symbol by dual affine transport; the affine part of dtQ+[V.grad,Q] is exactly zero, while constant affine stretching commutes with scalar Q",
        "egorov": "the remaining Heisenberg defect is exactly nonaffine resolved low-high role-interface work and vanishes for a genuinely affine transporter",
        "pressure": "Leray is applied before role extraction; no pressure forcing is introduced",
        "viscosity": "scalar Fourier Q commutes exactly with Delta; viscosity remains in the propagator and creates no role-interface source",
        "normalization": "Q is anchored to the already selected frozen physical transfer cell at the event time, so evolving the role does not alter that event's physical transfer normalization",
        "continuum_status": "the outer moving-role PDE and persistent low-low exclusion are exact; the companion interface theorem splits nonaffine work into conservative role flux plus existing strain provenance, leaving only recursive first-stop assembly",
    }


@dataclass(frozen=True)
class OuterRoleStress:
    samples: int
    worst_role_identity_residual: float
    worst_viscosity_commutator: float
    worst_affine_heisenberg_residual: float
    worst_egorov_identity_residual: float
    minimum_support_gap: float


def stress(samples: int = 50_000, seed: int = 20260809) -> OuterRoleStress:
    rng = np.random.default_rng(seed)
    wr = wv = wh = we = 0.0
    mg = float("inf")
    n = 4

    for _ in range(samples):
        T = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
        u = rng.normal(size=n) + 1j * rng.normal(size=n)
        V = rng.normal(size=n) + 1j * rng.normal(size=n)
        qdiag = rng.uniform(0.0, 1.0, size=n)
        Q = np.diag(qdiag).astype(complex)
        dQ = np.diag(rng.normal(size=n)).astype(complex)
        D = np.diag(-rng.uniform(0.0, 4.0, size=n)).astype(complex)
        nu = float(rng.uniform(0.0, 2.0))
        res = outer_role_identity_residual(
            tensor=T, u=u, V=V, Q=Q, dQ=dQ,
            viscosity_operator=D, viscosity=nu,
        )
        scale = max(1.0, np.linalg.norm(u), np.linalg.norm(V), np.linalg.norm(T))
        rel = float(np.linalg.norm(res)) / scale
        wr = max(wr, rel)
        if rel > 5e-12:
            raise AssertionError("exact outer-role algebra failed")
        vc = float(np.linalg.norm(Q @ D - D @ Q))
        wv = max(wv, vc)
        if vc > 2e-14:
            raise AssertionError("scalar Fourier/viscosity commutation model failed")

        A = rng.normal(size=(3, 3))
        A -= np.trace(A) / 3.0 * np.eye(3)
        xi = rng.normal(size=3)
        C = rng.normal(size=(3, 3))
        C = C.T @ C + 0.3 * np.eye(3)
        t = float(rng.uniform(-0.05, 0.05))
        _, grad, dtm = transported_gaussian_symbol(A, t, xi, C)
        hres = moving_multiplier_residual(A, xi, grad, dtm)
        wh = max(wh, abs(hres))
        if abs(hres) > 3e-11 * max(1.0, abs(dtm), np.linalg.norm(grad)):
            raise AssertionError("affine role did not have zero Heisenberg residual")

        m = 9
        kernel = rng.normal(size=m) + 1j * rng.normal(size=m)
        offsets = np.arange(m, dtype=float) - m // 2
        vel = rng.normal(size=m)
        df = rng.normal(size=m) + 1j * rng.normal(size=m)
        rate = float(rng.normal())
        lhs, rhs = exact_affine_subtracted_commutator(kernel, offsets, vel, df, rate)
        eres = float(np.linalg.norm(lhs - rhs)) / max(1.0, np.linalg.norm(lhs), np.linalg.norm(rhs))
        we = max(we, eres)
        if eres > 5e-12:
            raise AssertionError("exact affine-subtracted Egorov identity failed")

        K = float(rng.uniform(0.0, LOW_STRAIN_ACTION))
        role_lower = float(rng.uniform(ROLE_LOWER, 0.9))
        cut = float(rng.uniform(0.02, TRANSPORT_CUT))
        gap = persistent_low_low_gap(role_lower, K, cut)
        mg = min(mg, gap)
        if gap <= 0:
            raise AssertionError("transported selected role met the low-low output")

    return OuterRoleStress(samples, wr, wv, wh, we, mg)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-outer-moving-role-extraction"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    arb = arb_persistent_support_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "arb": arb, "stress": asdict(out), "sources": selected_outer_role_sources()}
    (args.outdir / "outer_moving_role_extraction.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Exact outer moving-role extraction\n\nStatus: **{cert['status']}**.\n\nLet `V=S_(N/4)u`, `h=u-V`, and `L_V f=P div(V tensor f+f tensor V)`.  For a time-dependent scalar divergence-free Fourier role `w=Q(t,D)u`, direct differentiation of Leray Navier--Stokes gives exactly\n\n`(dt+L_V-nu Delta)w = Q B(V,V) - Q B(h,h) + (dtQ+[L_V,Q])u`.\n\nThere is no pressure source, and scalar Fourier `Q` commutes with `Delta`.  Anchor `Q` at the already selected frozen physical transfer cell and transport its scalar symbol by the coherent affine dual flow.  On the low-strain branch `K<=1/30`, an anchored lower edge `3N/5` stays above `(3/5)e^(-1/30)N>N/2`, whereas `V tensor V` is supported in `B_(N/2)`.  Hence `Q B(V,V)=0` **throughout the whole role interval**, not merely at the anchor slice.\n\nThe selected outer-role equation is therefore\n\n`(dt+L_V-nu Delta)w = -Q P div(h tensor h) + R_Q`,\n`R_Q=(dtQ+[L_V,Q])u`.\n\nThe first term is the unique genuine quadratic high--high source; its positive coherent Hahn disintegration is already the physical binary causal work law.  For the second term, affine dual transport cancels the common affine advection Heisenberg term exactly, and constant affine stretching commutes with scalar `Q`.  The exact Egorov identity shows that `R_Q` is purely **nonaffine resolved low--high role-interface work** and vanishes for a genuinely affine transporter.\n\nStress: `{out.samples}` exact algebra/support/Egorov states\n- worst outer-role identity residual: `{out.worst_role_identity_residual:.3e}`\n- worst scalar-Q/viscosity commutator: `{out.worst_viscosity_commutator:.3e}`\n- worst affine Heisenberg residual: `{out.worst_affine_heisenberg_residual:.3e}`\n- worst affine-subtracted Egorov residual: `{out.worst_egorov_identity_residual:.3e}`\n- minimum sampled persistent low-low gap: `{out.minimum_support_gap:.6e}`\n\nThis changes the continuum frontier: constructing the moving outer role itself is no longer the missing PDE step.  The only unclosed outer-role term is the **work-level routing of the nonaffine Heisenberg interface**.  It must be shown to enter exactly once as either coherent deformation/critical `D_V` or physical role-relink/transfer loss; it must not be promoted to a new currency or silently absorbed as representation `Xi`.  No Navier--Stokes global-regularity conclusion is asserted.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
