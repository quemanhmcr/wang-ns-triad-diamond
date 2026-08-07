from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from .single_edge_certificate import RSTAR_HI, RSTAR_LO, float_rstar

# Conservative theorem constants.
LOCAL_LOG_BOX = Fraction(2, 25)          # existing |u|,|v| local box
FROZEN_STRAIN_TIME = Fraction(1, 25)     # d*T <= 0.04
STRAIN_COERCIVITY = Fraction(43, 100)    # c_dyn > 0.43
POINTWISE_HODGE = Fraction(3, 5)         # H(t) >= 0.6 d^2 t^2
AVERAGE_DEFICIT = Fraction(1, 10)        # avg Def >= 0.1 (dT)^2


def sym(A: np.ndarray) -> np.ndarray:
    return 0.5 * (A + A.T)


def tracefree_2x2(S: np.ndarray) -> np.ndarray:
    S = np.asarray(S, dtype=float)
    if S.shape != (2, 2):
        raise ValueError("expected a 2x2 matrix")
    return S - 0.5 * np.trace(S) * np.eye(2)


def extremal_parent_directions(r: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit parent directions and child direction at the symmetric extremizer.

    The two parents have magnitude r and sum to a unit child.  In the child-plane
    basis e1=child, e2=transverse, cos(phi)=1/(2r).
    """
    if r is None:
        r = float_rstar()
    c = 1.0 / (2.0 * r)
    if not (0.0 < c < 1.0):
        raise ValueError("invalid extremal radius")
    s = math.sqrt(1.0 - c * c)
    a = np.array([c, s], dtype=float)
    b = np.array([c, -s], dtype=float)
    child = np.array([1.0, 0.0], dtype=float)
    return a, b, child


def extremal_shape_rates(S_plane: np.ndarray, r: float | None = None) -> tuple[float, float, float]:
    """Signed log-shape rates (u_dot, v_dot, H_dot_speed^2).

    All three carrier wavevectors obey k_dot=-B k with symmetric driver
    S_plane=sym(B)|_plane.  Common scalar strain drops out.  The signed Hodge
    coordinates are u=log|b|-log|a| and
    v=log|c|-(log|a|+log|b|)/2-gamma_*.
    """
    D = tracefree_2x2(S_plane)
    na, nb, nc = extremal_parent_directions(r)
    qa = float(na @ D @ na)
    qb = float(nb @ D @ nb)
    qc = float(nc @ D @ nc)
    # d/dt log|k| = - n^T D n.
    udot = -(qb - qa)
    vdot = -qc + 0.5 * (qa + qb)
    speed2 = 0.5 * udot * udot + 2.0 * vdot * vdot
    return udot, vdot, speed2


def planar_strain_coercivity_ratio(S_plane: np.ndarray, r: float | None = None) -> float:
    D = tracefree_2x2(S_plane)
    n2 = float(np.sum(D * D))
    if n2 == 0.0:
        return math.inf
    return extremal_shape_rates(D, r)[2] / n2


def frozen_strain_shape_coords(D_plane: np.ndarray, t: float, r: float | None = None) -> tuple[float, float, float]:
    """Exact shape coordinates under a frozen symmetric planar strain.

    The common scalar planar strain is irrelevant, so D is reduced to its
    trace-free symmetric part and the carrier vectors evolve by exp(-D t).
    Returns signed u, v, and H=u^2/2+2v^2.
    """
    if t < 0:
        raise ValueError("t must be nonnegative")
    if r is None:
        r = float_rstar()
    D = tracefree_2x2(np.asarray(D_plane, dtype=float))
    if not np.allclose(D, D.T, atol=1e-13):
        raise ValueError("D must be symmetric")
    na, nb, _ = extremal_parent_directions(r)
    ka0 = r * na
    kb0 = r * nb
    kc0 = ka0 + kb0
    E = expm(-D * t)
    ka = E @ ka0
    kb = E @ kb0
    kc = E @ kc0
    la, lb, lc = math.log(np.linalg.norm(ka)), math.log(np.linalg.norm(kb)), math.log(np.linalg.norm(kc))
    gamma = -math.log(r)
    u = lb - la
    v = lc - 0.5 * (la + lb) - gamma
    H = 0.5 * u * u + 2.0 * v * v
    return u, v, H


def frozen_strain_average_deficit_lower(d: float, T: float) -> float:
    """Theorem lower bound for a frozen principal-strain episode.

    Here +/-d are the eigenvalues of the trace-free planar strain.  Starting
    from the exact extremal triad, if d*T<=1/25 then the entire episode stays in
    the certified single-edge local box and

        (1/T) int Def(t) dt >= (1/10) (dT)^2.
    """
    if d < 0 or T < 0:
        raise ValueError("d and T must be nonnegative")
    if d * T > float(FROZEN_STRAIN_TIME) + 1e-15:
        raise ValueError("the theorem requires d*T <= 1/25")
    return float(AVERAGE_DEFICIT) * (d * T) ** 2


def fourier_gaussian_rhs(P: np.ndarray, kappa: np.ndarray, A: np.ndarray, nu: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact scalar-envelope ODE for affine advection-diffusion in Fourier space.

    For
      f_t + (A x).grad f = nu Delta f,   tr A=0,
    and a Gaussian Fourier exponent with precision P and peak kappa,

      P_dot = A P + P A^T + 2 nu I,
      kappa_dot = -A^T kappa - 2 nu P^{-1} kappa.

    The vector stretching/Leray dynamics of a velocity packet acts on its
    polarization and is recorded separately by kelvin_amplitude_rhs.
    """
    P = np.asarray(P, dtype=float)
    kappa = np.asarray(kappa, dtype=float)
    A = np.asarray(A, dtype=float)
    n = P.shape[0]
    if P.shape != (n, n) or A.shape != (n, n) or kappa.shape != (n,):
        raise ValueError("dimension mismatch")
    if nu < 0:
        raise ValueError("nu must be nonnegative")
    Pinv_k = np.linalg.solve(P, kappa)
    Pdot = A @ P + P @ A.T + 2.0 * nu * np.eye(n)
    kdot = -A.T @ kappa - 2.0 * nu * Pinv_k
    return Pdot, kdot


def dual_center_derivative(P: np.ndarray, kappa: np.ndarray, A: np.ndarray, nu: float) -> np.ndarray:
    """Derivative of b=P kappa; exactly A b under the Gaussian ODE."""
    Pdot, kdot = fourier_gaussian_rhs(P, kappa, A, nu)
    return Pdot @ kappa + P @ kdot


def logdet_precision_rate(P: np.ndarray, A: np.ndarray, nu: float) -> float:
    """d/dt log det P = 2 tr A + 2 nu tr(P^{-1})."""
    P = np.asarray(P, dtype=float)
    A = np.asarray(A, dtype=float)
    return 2.0 * float(np.trace(A)) + 2.0 * nu * float(np.trace(np.linalg.inv(P)))


def effective_shape_driver(P: np.ndarray, A: np.ndarray, nu: float) -> np.ndarray:
    """Symmetric carrier-length driver for the Gaussian spectral peak.

    kappa_dot = -B kappa with B=A^T+2 nu P^{-1}; only sym(B) changes lengths.
    """
    return sym(A.T + 2.0 * nu * np.linalg.inv(P))


def kelvin_amplitude_rhs(a: np.ndarray, k: np.ndarray, A: np.ndarray, nu: float) -> np.ndarray:
    """Exact Kelvin-mode amplitude ODE for linearized incompressible NS.

    Around U=A x, a transverse plane wave has
      k_dot=-A^T k,
      a_dot=-A a + 2 k (k.Aa)/|k|^2 - nu |k|^2 a.
    The pressure/Leray term preserves k.a=0 and does no direct amplitude work.
    """
    a = np.asarray(a, dtype=float)
    k = np.asarray(k, dtype=float)
    A = np.asarray(A, dtype=float)
    k2 = float(k @ k)
    if k2 <= 0:
        raise ValueError("k must be nonzero")
    Aa = A @ a
    return -Aa + 2.0 * k * float(k @ Aa) / k2 - nu * k2 * a


def kelvin_energy_rate(a: np.ndarray, k: np.ndarray, A: np.ndarray, nu: float) -> float:
    """d/dt |a|^2 for a transverse Kelvin mode."""
    a = np.asarray(a, dtype=float)
    k = np.asarray(k, dtype=float)
    S = sym(np.asarray(A, dtype=float))
    return -2.0 * float(a @ S @ a) - 2.0 * nu * float(k @ k) * float(a @ a)


def affine_window_material_error_bound(
    hessian_u: float,
    material_radius: float,
    deformation_norm: float,
    inverse_deformation_norm: float,
    grad_chi0: float,
) -> float:
    """Pointwise bound for the moving affine-window material derivative.

    Let X_dot=U(X), F_dot=(grad U)(X) F and
      chi(x,t)=chi0(F^{-1}(x-X)/R).
    The affine Taylor part cancels exactly.  If ||grad^2 U||<=H on the window,
    then
      |(partial_t+U.grad)chi|
      <= (H/2) R ||F||^2 ||F^{-1}|| ||grad chi0||.
    """
    vals = (hessian_u, material_radius, deformation_norm, inverse_deformation_norm, grad_chi0)
    if any(v < 0 for v in vals):
        raise ValueError("bounds must be nonnegative")
    return 0.5 * hessian_u * material_radius * deformation_norm ** 2 * inverse_deformation_norm * grad_chi0


def arb_affine_grain_certificate() -> dict[str, str]:
    """Rigorous constants for extremal strain rigidity and frozen-strain cost."""
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover - exercised in Actions
        raise RuntimeError("python-flint is required for the rigorous affine-grain certificate") from exc
    ctx.prec = 160

    def aq(q: Fraction):
        return arb(q.numerator) / q.denominator

    rlo, rhi = aq(RSTAR_LO), aq(RSTAR_HI)
    r = rlo.union(rhi)
    # cos(phi)=1/(2r); c_dyn=4 sin(phi)^4 is the minimum singular-value
    # coefficient relating Hodge shape speed to planar trace-free strain.
    sin2 = 1 - 1 / (4 * r * r)
    cdyn = 4 * sin2 * sin2
    if not (cdyn > aq(STRAIN_COERCIVITY)):
        raise AssertionError(f"strain coercivity failed: {cdyn}")

    # If D has eigenvalues +/-d then ||D||_F^2=2d^2, so initial Hodge-speed
    # squared is >= 2*c_dyn*d^2 > 0.86 d^2.  Each log-length acceleration is
    # <=2d^2.  In Hodge coordinates z=(u/sqrt2,sqrt2 v),
    # ||z_dot(t)-z_dot(0)||<=sqrt(40)d^2 t.  Integrating gives
    # ||z(t)|| >= d t [sqrt(0.86)-sqrt(10) d t].  At d t<=1/25 this square is
    # >3/5, yielding H>=3/5 d^2t^2 and avg Def>=1/10(dT)^2.
    speed2 = arb(86) / 100
    bracket = speed2.sqrt() - arb(10).sqrt() / 25
    bracket2 = bracket * bracket
    if not (bracket2 > aq(POINTWISE_HODGE)):
        raise AssertionError(f"frozen strain pointwise bound failed: {bracket2}")

    # Local-box containment: |u|,|v| <= 2 d t <=2/25.
    if not (2 * aq(FROZEN_STRAIN_TIME) <= aq(LOCAL_LOG_BOX)):
        raise AssertionError("frozen strain episode does not stay in local box")

    # 1/2 from Def>=H/2, 1/3 from averaging t^2, 3/5 pointwise H coefficient.
    avg = arb(1) / 2 * aq(POINTWISE_HODGE) / 3
    if not (avg >= aq(AVERAGE_DEFICIT)):
        raise AssertionError(f"average deficit coefficient failed: {avg}")

    return {
        "rstar_ball": str(r),
        "strain_coercivity_ball": str(cdyn),
        "strain_coercivity_lower": "43/100",
        "frozen_strain_time": "1/25",
        "pointwise_hodge_bracket_ball": str(bracket2),
        "pointwise_hodge_lower": "3/5",
        "average_deficit_lower": "1/10",
        "status": "CERTIFIED",
    }


def _random_tracefree_matrix(rng: np.random.Generator, n: int = 3) -> np.ndarray:
    A = rng.normal(size=(n, n))
    A -= np.trace(A) / n * np.eye(n)
    return A


def stress(samples: int = 50_000, seed: int = 20260807) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    r = float_rstar()
    worst_coercivity = float("inf")
    worst_hodge_ratio = float("inf")
    worst_dual_residual = 0.0
    worst_logdet_residual = 0.0
    worst_kelvin_constraint = 0.0
    worst_kelvin_energy = 0.0
    worst_viscous_neutrality = 0.0

    for _ in range(samples):
        # Random orientation of a fixed planar trace-free principal strain.
        psi = float(rng.uniform(0.0, math.pi))
        cp, sp = math.cos(psi), math.sin(psi)
        R = np.array([[cp, -sp], [sp, cp]])
        D = R @ np.diag([1.0, -1.0]) @ R.T
        ratio = planar_strain_coercivity_ratio(D, r)
        worst_coercivity = min(worst_coercivity, ratio)

        t = float(rng.uniform(1e-8, float(FROZEN_STRAIN_TIME)))
        u, v, H = frozen_strain_shape_coords(D, t, r)
        if abs(u) > float(LOCAL_LOG_BOX) + 1e-10 or abs(v) > float(LOCAL_LOG_BOX) + 1e-10:
            raise AssertionError("frozen strain left certified local box")
        hratio = H / (t * t)
        worst_hodge_ratio = min(worst_hodge_ratio, hratio)
        if hratio + 2e-10 < float(POINTWISE_HODGE):
            raise AssertionError(("frozen strain Hodge lower bound failed", hratio, psi, t))

        # Gaussian affine-advection/diffusion identities.
        A = _random_tracefree_matrix(rng, 3)
        X = rng.normal(size=(3, 3))
        P = X @ X.T + 0.3 * np.eye(3)
        k = rng.normal(size=3)
        nu = float(rng.uniform(0.0, 2.0))
        bdot = dual_center_derivative(P, k, A, nu)
        target_bdot = A @ (P @ k)
        scale = max(1.0, np.linalg.norm(target_bdot))
        worst_dual_residual = max(worst_dual_residual, np.linalg.norm(bdot - target_bdot) / scale)
        Pdot, _ = fourier_gaussian_rhs(P, k, A, nu)
        direct_logdet = float(np.trace(np.linalg.solve(P, Pdot)))
        target_logdet = logdet_precision_rate(P, A, nu)
        worst_logdet_residual = max(worst_logdet_residual, abs(direct_logdet - target_logdet) / max(1.0, abs(target_logdet)))

        # At isotropic width viscosity adds only scalar length drift, hence no
        # trace-free planar shape driver at first order.
        p = float(rng.uniform(0.2, 4.0))
        Piso = p * np.eye(3)
        Mnu = effective_shape_driver(Piso, A, nu)
        M0 = sym(A)
        Dnu = tracefree_2x2(Mnu[:2, :2])
        D0 = tracefree_2x2(M0[:2, :2])
        worst_viscous_neutrality = max(worst_viscous_neutrality, float(np.linalg.norm(Dnu - D0)))

        # Kelvin polarization: preserve k.a=0 and satisfy exact energy work.
        kK = rng.normal(size=3)
        kK /= np.linalg.norm(kK)
        a = rng.normal(size=3)
        a -= kK * float(kK @ a)
        kdot = -A.T @ kK
        adot = kelvin_amplitude_rhs(a, kK, A, nu)
        constraint_derivative = float(kdot @ a + kK @ adot)
        worst_kelvin_constraint = max(worst_kelvin_constraint, abs(constraint_derivative))
        energy_direct = 2.0 * float(a @ adot)
        energy_target = kelvin_energy_rate(a, kK, A, nu)
        worst_kelvin_energy = max(worst_kelvin_energy, abs(energy_direct - energy_target) / max(1.0, abs(energy_target)))

    return {
        "samples": samples,
        "worst_planar_strain_coercivity": worst_coercivity,
        "worst_frozen_hodge_ratio": worst_hodge_ratio,
        "worst_dual_center_relative_residual": worst_dual_residual,
        "worst_logdet_relative_residual": worst_logdet_residual,
        "worst_kelvin_constraint_residual": worst_kelvin_constraint,
        "worst_kelvin_energy_relative_residual": worst_kelvin_energy,
        "worst_isotropic_viscosity_shape_residual": worst_viscous_neutrality,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-affine-grain"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_affine_grain_certificate()
    result = stress(args.samples)
    payload = {"certificate": cert, "stress": result}
    (args.outdir / "affine_grain_dynamics.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Affine Gaussian grain dynamics\n\nStatus: **{cert['status']}** for the extremal strain-rigidity and frozen-strain cost.\n\n- extremal planar strain coercivity: `> {cert['strain_coercivity_lower']}`\n- certified coercivity enclosure: `{cert['strain_coercivity_ball']}`\n- frozen dimensionless strain-time: `d T <= {cert['frozen_strain_time']}`\n- pointwise Hodge coefficient: `> {cert['pointwise_hodge_lower']}`\n- bracket enclosure: `{cert['pointwise_hodge_bracket_ball']}`\n- time-averaged single-edge deficit: `>= {cert['average_deficit_lower']} (dT)^2`\n- adversarial affine/Gaussian checks: `{result['samples']}`\n- worst planar coercivity seen: `{result['worst_planar_strain_coercivity']:.9f}`\n- worst frozen Hodge ratio `H/(d t)^2`: `{result['worst_frozen_hodge_ratio']:.9f}`\n- worst Gaussian dual-center residual: `{result['worst_dual_center_relative_residual']:.3e}`\n- worst log-det residual: `{result['worst_logdet_relative_residual']:.3e}`\n- worst Kelvin transversality residual: `{result['worst_kelvin_constraint_residual']:.3e}`\n- worst Kelvin energy-work residual: `{result['worst_kelvin_energy_relative_residual']:.3e}`\n- worst isotropic-viscosity shape residual: `{result['worst_isotropic_viscosity_shape_residual']:.3e}`\n\nThe theorem is a local affine/Kelvin-Gaussian dynamics statement, not a full\nNavier--Stokes packet-lifetime theorem.  A common rotation and common planar\nscalar strain are gauge directions; the charged quantity is the non-conformal\ntrace-free strain seen by the extremal triad plane.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
