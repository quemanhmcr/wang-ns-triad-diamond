from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


def periodic_convolution(kernel: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Qf[x]=sum_j K[j] f[x-j] on a periodic discrete line."""
    K = np.asarray(kernel, complex)
    f = np.asarray(f, complex)
    if K.ndim != 1 or f.ndim != 1 or K.shape != f.shape:
        raise ValueError("matching periodic vectors required")
    out = np.zeros_like(f, dtype=complex)
    for j, kj in enumerate(K):
        out += kj * np.roll(f, j)
    return out


def affine_dt_action(kernel: np.ndarray, offsets: np.ndarray, affine_rate: float, df: np.ndarray) -> np.ndarray:
    """Instantaneous multiplier motion cancelling the affine commutator.

    This is the kernel-side form of partial_t m-(A^T xi).grad_xi m=0 in one
    scalar dimension: dot Q f = -sum_j K_j A y_j (df)(x-y_j).
    """
    K = np.asarray(kernel, complex)
    y = np.asarray(offsets, float)
    df = np.asarray(df, complex)
    if K.shape != y.shape or K.shape != df.shape or K.ndim != 1:
        raise ValueError("matching one-dimensional data required")
    out = np.zeros_like(df, dtype=complex)
    for j, kj in enumerate(K):
        out -= kj * affine_rate * y[j] * np.roll(df, j)
    return out


def exact_affine_subtracted_commutator(
    kernel: np.ndarray,
    offsets: np.ndarray,
    velocity: np.ndarray,
    df: np.ndarray,
    affine_rate: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Exact kernel identity for dot Q+[V d_x,Q] after affine cancellation.

    lhs = dot Q f + V Q(df) - Q(V df)
    rhs = sum_j K_j [V(x)-V(x-y_j)-A y_j] df(x-y_j).

    No Taylor expansion is used in this identity.
    """
    K = np.asarray(kernel, complex)
    y = np.asarray(offsets, float)
    V = np.asarray(velocity, complex)
    df = np.asarray(df, complex)
    if not (K.shape == y.shape == V.shape == df.shape) or K.ndim != 1:
        raise ValueError("matching periodic data required")
    lhs = affine_dt_action(K, y, affine_rate, df) + V * periodic_convolution(K, df) - periodic_convolution(K, V * df)
    rhs = np.zeros_like(df, dtype=complex)
    for j, kj in enumerate(K):
        rhs += kj * (V - np.roll(V, j) - affine_rate * y[j]) * np.roll(df, j)
    return lhs, rhs


def affine_subtracted_increment_upper(hessian_sup: float, x_distance: float, y_distance: float) -> float:
    """Integral-Taylor bound |V(x)-V(x-y)-A y|.

    A=grad V(X).  If ||Hess V||<=H on the relevant segments, then
      <= H (|x-X||y| + |y|^2/2).
    """
    if min(hessian_sup, x_distance, y_distance) < 0:
        raise ValueError("nonnegative geometric data required")
    return hessian_sup * (x_distance * y_distance + 0.5 * y_distance**2)


def scalar_cell_egorov_l2_upper(
    *,
    hessian_sup: float,
    frequency_scale: float,
    kernel_l1: float,
    kernel_m1: float,
    kernel_m2: float,
    packet_x_moment: float,
    packet_grad_moment: float,
    packet_xgrad_moment: float,
) -> float:
    """L2 upper for advection commutator plus affine-amplitude remainder.

    Dimensionless kernel moments are defined by
      ||K_N||_1<=k0,
      int |y||K_N(y)|dy <= m1/N,
      int |y|^2|K_N(y)|dy <= m2/N^2.
    Packet moments are
      Mx=N |||x-X|f||/||f||,
      Mg=||grad f||/(N||f||),
      Mxg=|||x-X|grad f||/||f||.

    Then ||R_Eg f||/||f|| <= (H/N)[k0 Mx+m1 Mxg+(3/2)m2 Mg].
    """
    nonnegative = (hessian_sup, kernel_l1, kernel_m1, kernel_m2, packet_x_moment, packet_grad_moment, packet_xgrad_moment)
    if not math.isfinite(frequency_scale) or frequency_scale <= 0 or any((not math.isfinite(v) or v < 0) for v in nonnegative):
        raise ValueError("finite nonnegative moments and positive finite frequency required")
    C = kernel_l1 * packet_x_moment + kernel_m1 * packet_xgrad_moment + 1.5 * kernel_m2 * packet_grad_moment
    return hessian_sup / frequency_scale * C


def isotropic_gaussian_packet_moments(carrier_ratio: float) -> dict[str, float]:
    """Exact moments for |g|^2 with z=N(x-X) standard Gaussian.

    g is proportional to exp(-N^2|x-X|^2/4) exp(i k.x), q=|k|/N.
    """
    if carrier_ratio < 0 or not math.isfinite(carrier_ratio):
        raise ValueError("nonnegative finite carrier ratio required")
    q = carrier_ratio
    return {
        "Mx": math.sqrt(3.0),
        "Mg": math.sqrt(q * q + 3.0 / 4.0),
        "Mxg": math.sqrt(15.0 / 4.0 + 3.0 * q * q),
    }


def center_jet_shear_countermodel(amplitude: float = 1.0, frequency_ratio: float = 1.0 / 8.0) -> dict[str, float]:
    """Divergence-free B_(N/4) shear with zero center affine/quadratic jet.

    V_2(x1)=a[sin(r x1)-1/2 sin(2r x1)], r=frequency_ratio*N.
    At x1=0: V=A=H=0 but d_1^3 V_2=3 a r^3.
    Taking r/N=1/8 keeps both Fourier modes inside N/4.
    Values are normalized by setting N=1.
    """
    if amplitude == 0 or frequency_ratio <= 0 or 2.0 * frequency_ratio > 0.25 + 1e-15:
        raise ValueError("nonzero amplitude and both modes inside the N/4 transporter required")
    a = float(amplitude)
    r = float(frequency_ratio)
    return {
        "V_center": 0.0,
        "gradient_center": 0.0,
        "hessian_center": 0.0,
        "third_derivative_center": 3.0 * a * r**3,
        "first_mode_over_N": r,
        "second_mode_over_N": 2.0 * r,
        "lowpass_cut_over_N": 0.25,
        "cubic_taylor_coefficient": 0.5 * a * r**3,
    }


def coherent_deformation_variance_from_hermite(degrees: Sequence[int], squared_coefficients: Sequence[float]) -> tuple[float, float]:
    """Gaussian Poincare identity/inequality in Hermite coordinates.

    For F(z)-E F=sum_{alpha!=0} c_alpha H_alpha, orthonormal Hermites give
      Var(F)=sum |c_alpha|^2,
      E|grad F|^2=sum |alpha| |c_alpha|^2 >= Var(F).
    Matrix-valued F follows componentwise, so only degrees and squared norms are
    needed here.
    """
    d = np.asarray(degrees, int)
    w = np.asarray(squared_coefficients, float)
    if d.ndim != 1 or w.shape != d.shape or np.any(d < 1) or np.any(w < 0):
        raise ValueError("positive Hermite degrees and nonnegative squared coefficients required")
    var = float(w.sum())
    grad = float(np.dot(d, w))
    if grad + 1e-14 < var:
        raise AssertionError("Gaussian Poincare spectral gap failed")
    return var, grad


def affine_invariant_gradient_matrix(M: np.ndarray, L: np.ndarray, A: np.ndarray) -> np.ndarray:
    """Intrinsic resolved-gradient matrix L^-1 A L and its common-affine covariance."""
    M = np.asarray(M, float)
    L = np.asarray(L, float)
    A = np.asarray(A, float)
    if M.shape != (3, 3) or L.shape != (3, 3) or A.shape != (3, 3):
        raise ValueError("3x3 matrices required")
    return np.linalg.solve(L, A @ L)


def affine_invariant_gradient_residual(M: np.ndarray, L: np.ndarray, A: np.ndarray) -> float:
    """Under L->ML, A->M A M^-1 the intrinsic gradient is unchanged."""
    M = np.asarray(M, float)
    if abs(np.linalg.det(M)) < 1e-12:
        raise ValueError("invertible affine map required")
    old = affine_invariant_gradient_matrix(M, L, A)
    newL = M @ np.asarray(L, float)
    newA = M @ np.asarray(A, float) @ np.linalg.inv(M)
    new = np.linalg.solve(newL, newA @ newL)
    return float(np.linalg.norm(old - new))


@dataclass(frozen=True)
class ResolvedEgorovStress:
    samples: int
    worst_exact_commutator_relative_residual: float
    minimum_taylor_bound_margin: float
    minimum_gaussian_poincare_margin: float
    worst_affine_invariant_gradient_residual: float
    countermodel_third_derivative: float


def stress(samples: int = 50_000, seed: int = 20260808) -> ResolvedEgorovStress:
    rng = np.random.default_rng(seed)
    wc = wa = 0.0
    mt = mp = float("inf")
    for _ in range(samples):
        n = int(rng.integers(5, 35))
        K = rng.normal(size=n) + 1j * rng.normal(size=n)
        K /= max(np.sum(np.abs(K)), 1e-30)
        V = rng.normal(size=n)
        df = rng.normal(size=n) + 1j * rng.normal(size=n)
        # Offsets need not be the torus geometry for the algebraic cancellation test;
        # they are the displacements carried by the convolution slots.
        y = rng.normal(size=n)
        A = float(rng.normal())
        lhs, rhs = exact_affine_subtracted_commutator(K, y, V, df, A)
        scale = max(1.0, float(np.linalg.norm(lhs)), float(np.linalg.norm(rhs)))
        r = float(np.linalg.norm(lhs - rhs)) / scale
        wc = max(wc, r)
        if r > 3e-11:
            raise AssertionError("affine-subtracted commutator identity failed")

        # Constant-Hessian quadratic velocity gives an exact test of the integral Taylor bound.
        H = rng.normal(size=(3, 3, 3))
        H = 0.5 * (H + H.swapaxes(1, 2))
        X = rng.normal(size=3)
        x = X + rng.normal(size=3)
        yy = rng.normal(size=3)
        A3 = rng.normal(size=(3, 3))
        V0 = rng.normal(size=3)
        def qvel(z: np.ndarray) -> np.ndarray:
            zz = z - X
            return V0 + A3 @ zz + 0.5 * np.einsum("ijk,j,k->i", H, zz, zz)
        inc = qvel(x) - qvel(x - yy) - A3 @ yy
        Hnorm = float(np.linalg.norm(H))
        ub = affine_subtracted_increment_upper(Hnorm, float(np.linalg.norm(x - X)), float(np.linalg.norm(yy)))
        margin = ub - float(np.linalg.norm(inc))
        mt = min(mt, margin)
        if margin < -3e-12 * max(1.0, ub):
            raise AssertionError("curvature Taylor bound failed")

        m = int(rng.integers(1, 50))
        deg = rng.integers(1, 9, size=m)
        weights = rng.random(m)
        var, grad = coherent_deformation_variance_from_hermite(deg, weights)
        mp = min(mp, grad - var)

        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        scales = np.exp(rng.uniform(-1.0, 1.0, size=3))
        L = Q @ np.diag(scales)
        MM = np.eye(3) + 0.1 * rng.normal(size=(3, 3))
        if abs(np.linalg.det(MM)) < 0.2:
            MM += np.eye(3)
        AA = rng.normal(size=(3, 3))
        ar = affine_invariant_gradient_residual(MM, L, AA)
        wa = max(wa, ar)
        if ar > 2e-11 * max(1.0, np.linalg.norm(AA)):
            raise AssertionError("intrinsic gradient lost common-affine invariance")

    cm = center_jet_shear_countermodel()
    if cm["hessian_center"] != 0.0 or abs(cm["third_derivative_center"]) <= 0.0:
        raise AssertionError("center-jet countermodel failed")
    return ResolvedEgorovStress(samples, wc, mt, mp, wa, float(cm["third_derivative_center"]))


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_AFFINE_SUBTRACTED_EGOROV_IDENTITY_AND_CENTER_JET_COUNTERMODEL__COHERENT_VARIANCE_REDESIGN",
        "egorov_identity": "(dot Q+[V.grad,Q])f = int K_N(y)[V(x)-V(x-y)-A y].grad f(x-y)dy",
        "curvature_bound": "|V(x)-V(x-y)-A y| <= H(|x-X||y|+|y|^2/2)",
        "l2_bound": "||R_Eg f||/||f|| <= (H/N)[k0 Mx+m1 Mxg+(3/2)m2 Mg] including affine-amplitude mismatch",
        "center_countermodel": "V=(0,a[sin(rx1)-sin(2rx1)/2],0), r=N/8: A(0)=H(0)=0 but d1^3 V2(0)=3ar^3",
        "false_bridge": "full non-affine moving-role residual cannot be bounded solely by the center Hessian B(X)",
        "coherent_repair": "use the affine-invariant coherent gradient variance Var_gamma[L^-1 grad V(X+Lz)L]",
        "gaussian_poincare": "coherent gradient variance <= E_gamma ||L^-1 Hess V(X+Lz)[L,L]||^2",
        "physical_interpretation": "measures resolved deformation variation across the coherent eddy; common affine deformation has zero variance",
        "continuum_status": "must still collide the coherent/cellwise curvature variance with the existing H1/H3/source/transfer currencies or prove it is uniformly absorbable on the flat branch",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-role-egorov"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    payload = {"certificate": cert, "stress": asdict(out), "countermodel": center_jet_shear_countermodel(), "gaussian_moments_q1": isotropic_gaussian_packet_moments(1.0)}
    (args.outdir / "resolved_role_egorov.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Resolved-role Egorov identity and coherent-deformation redesign\n\nStatus: **{cert['status']}**.\n\nFor a smooth moving Fourier cell with convolution kernel `K_N`, after transporting the multiplier by the affine jet `A=grad V(X)`, the complete scalar advection commutator is exactly\n\n`(dot Q+[V.grad,Q])f = int K_N(y)[V(x)-V(x-y)-A y].grad f(x-y) dy`.\n\nNo higher Taylor expansion is hidden here.  If `||Hess V||<=H`, the exact integral-Taylor bound is\n\n`|V(x)-V(x-y)-A y| <= H(|x-X||y|+|y|^2/2)`.\n\nIncluding the vector amplitude mismatch `Q[(grad V-A)f]`, fixed dimensionless kernel moments and packet moments give\n\n`||R_Eg f||/||f|| <= (H/N)[k0 Mx+m1 Mxg+(3/2)m2 Mg]`.\n\nHowever **center Hessian alone is not enough**.  The divergence-free strict-lowpass shear\n\n`V=(0,a[sin(r x1)-sin(2r x1)/2],0)`, `r=N/8`,\n\nhas `V(0)=0`, `grad V(0)=0`, `Hess V(0)=0`, but `d1^3 V2(0)=3 a r^3 !=0`.  Thus a bridge claiming that every non-affine moving-role residual is controlled only by the center tensor `B(X)` is false even for a smooth `B_(N/4)` transporter.\n\nThe natural replacement is a coherent, affine-invariant deformation observable.  For standard Gaussian coherent coordinate `z`, put\n\n`F(z)=L^-1 grad V(X+Lz)L`,  `Abar=E_gamma F`,\n`K_coh^2=E_gamma ||F-Abar||^2`.\n\nCommon affine flow makes `K_coh=0`.  Gaussian Poincare, componentwise, gives\n\n`K_coh^2 <= E_gamma ||L^-1 Hess V(X+Lz)[L,L]||^2`.\n\nSo the unresolved outer-role defect is not an arbitrary packet error: it is **resolved deformation variance across the physical coherent eddy**.  This catches spatial non-affinity missed by point sampling while preserving affine covariance.\n\nStress: `{out.samples}` exact-commutator/curvature/Hermite/covariance states\n- worst exact commutator relative residual: `{out.worst_exact_commutator_relative_residual:.3e}`\n- minimum Taylor-bound margin: `{out.minimum_taylor_bound_margin:.3e}`\n- minimum Gaussian-Poincare margin: `{out.minimum_gaussian_poincare_margin:.3e}`\n- worst affine-invariant gradient residual: `{out.worst_affine_invariant_gradient_residual:.3e}`\n- countermodel normalized third derivative: `{out.countermodel_third_derivative:.6e}`\n\nThe next theorem must either route this coherent/cellwise curvature variance into the existing H1/H3/source/transfer currencies, or prove that on the certified low-strain/near-extremal branch its contribution is uniformly absorbable by flat erosion.  No such closure is asserted here.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
