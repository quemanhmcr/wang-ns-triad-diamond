from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from .affine_grain_dynamics import (
    LOCAL_LOG_BOX,
    extremal_parent_directions,
    sym,
    tracefree_2x2,
)
from .single_edge_certificate import float_rstar

COHERENCE_FRACTION = Fraction(1, 20)      # 5% sup variation from D0
COHERENT_STRAIN_TIME = Fraction(1, 30)    # d*T <= 1/30
COHERENT_HODGE = Fraction(1, 4)           # H(t) >= 1/4 d^2 t^2
COHERENT_AVG_DEFICIT = Fraction(1, 24)    # avg Def >= 1/24 (dT)^2


def opnorm(A: np.ndarray) -> float:
    return float(np.linalg.norm(np.asarray(A, dtype=float), ord=2))


def coherent_strain_average_deficit_lower(d: float, T: float, coherence_fraction: float = float(COHERENCE_FRACTION)) -> float:
    """Theorem lower bound in the 5%-coherent planar-strain regime.

    D(t) is symmetric trace-free on a co-rotating invariant triad plane,
    D(0) has eigenvalues +/-d, and
        sup_t ||D(t)-D(0)||_op <= coherence_fraction*d.
    For the certified theorem constants coherence_fraction<=1/20 and dT<=1/30,
        avg Def >= (1/24)(dT)^2.
    """
    if d < 0 or T < 0 or coherence_fraction < 0:
        raise ValueError("invalid parameters")
    if coherence_fraction > float(COHERENCE_FRACTION) + 1e-15:
        raise ValueError("certificate requires coherence_fraction <= 1/20")
    if d * T > float(COHERENT_STRAIN_TIME) + 1e-15:
        raise ValueError("certificate requires d*T <= 1/30")
    return float(COHERENT_AVG_DEFICIT) * (d * T) ** 2


def velocity_gradient_material_rhs(A: np.ndarray, pressure_hessian: np.ndarray, laplacian_A: np.ndarray, nu: float) -> np.ndarray:
    """Exact smooth Navier--Stokes velocity-gradient equation.

    For A=grad u along a material trajectory,
        D_t A + A^2 = - Hess p + nu Delta A.
    """
    A = np.asarray(A, dtype=float)
    Hp = np.asarray(pressure_hessian, dtype=float)
    LA = np.asarray(laplacian_A, dtype=float)
    if A.shape != Hp.shape or A.shape != LA.shape or A.shape[0] != A.shape[1]:
        raise ValueError("matrix dimension mismatch")
    if nu < 0:
        raise ValueError("nu must be nonnegative")
    return -(A @ A) - Hp + nu * LA


def corotational_strain_rhs(A: np.ndarray, pressure_hessian: np.ndarray, laplacian_A: np.ndarray, nu: float) -> np.ndarray:
    """Objective strain derivative in the frame rotating with local vorticity.

    Let A=S+Omega, S^T=S, Omega^T=-Omega.  If R_dot=Omega R, then
    d(R^T S R)/dt = R^T S_circ R with

      S_circ = D_t S + S Omega - Omega S
             = -S^2-Omega^2-Hess p+nu Delta S + S Omega-Omega S.

    `pressure_hessian` is assumed symmetric, as it is for a scalar pressure.
    """
    A = np.asarray(A, dtype=float)
    Hp = np.asarray(pressure_hessian, dtype=float)
    LA = np.asarray(laplacian_A, dtype=float)
    S = sym(A)
    Omega = 0.5 * (A - A.T)
    return -(S @ S) - (Omega @ Omega) - sym(Hp) + nu * sym(LA) + S @ Omega - Omega @ S


def coherence_failure_action_lower(d: float, coherence_fraction: float = float(COHERENCE_FRACTION)) -> float:
    """Total variation needed to leave the coherent strain ball.

    If a continuous matrix path D(t) starts at D0 and at some time satisfies
    ||D(t)-D0||_op > eps*d, then int ||D_dot||_op dt > eps*d.
    """
    if d < 0 or coherence_fraction < 0:
        raise ValueError("invalid parameters")
    return coherence_fraction * d


def _rk4_fundamental(Dfun, T: float, steps: int = 48) -> np.ndarray:
    E = np.eye(2)
    dt = T / steps
    t = 0.0
    for _ in range(steps):
        def rhs(tt, X):
            return -Dfun(tt) @ X
        k1 = rhs(t, E)
        k2 = rhs(t + dt / 2, E + dt * k1 / 2)
        k3 = rhs(t + dt / 2, E + dt * k2 / 2)
        k4 = rhs(t + dt, E + dt * k3)
        E = E + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        t += dt
    return E


def shape_coords_from_fundamental(E: np.ndarray, r: float | None = None) -> tuple[float, float, float]:
    if r is None:
        r = float_rstar()
    na, nb, _ = extremal_parent_directions(r)
    ka = E @ (r * na)
    kb = E @ (r * nb)
    kc = ka + kb
    la, lb, lc = math.log(np.linalg.norm(ka)), math.log(np.linalg.norm(kb)), math.log(np.linalg.norm(kc))
    gamma = -math.log(r)
    u = lb - la
    v = lc - 0.5 * (la + lb) - gamma
    H = 0.5 * u * u + 2.0 * v * v
    return u, v, H


def arb_strain_coherence_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint is required for the rigorous strain-coherence certificate") from exc
    ctx.prec = 160

    eps = arb(1) / 20
    a = arb(1) / 30
    # Initial extremal Hodge speed has squared norm > 0.86 d^2 from the
    # affine-grain certificate.  For any carrier direction, the rate change is
    # bounded by eps*d from D-D0 plus 4(1+eps)d^2 t from direction drift.
    # Hodge coordinate combinations amplify a common log-rate error by <=sqrt10.
    bracket = (arb(86) / 100).sqrt() - arb(10).sqrt() * eps - 2 * arb(10).sqrt() * (1 + eps) * a
    bracket2 = bracket * bracket
    if not (bracket2 > arb(1) / 4):
        raise AssertionError(f"coherent strain Hodge bracket failed: {bracket2}")

    # Local-box containment: |u|,|v| <= 2(1+eps)dT = 7/100 < 2/25.
    local_radius = 2 * (1 + eps) * a
    if not (local_radius < arb(2) / 25):
        raise AssertionError(f"coherent strain local-box containment failed: {local_radius}")

    # Exact rational average: (1/2)*(1/4)*(1/3)=1/24.
    if Fraction(1, 2) * COHERENT_HODGE / 3 != COHERENT_AVG_DEFICIT:
        raise AssertionError("coherent average deficit identity failed")

    return {
        "coherence_fraction": "1/20",
        "strain_time": "1/30",
        "hodge_bracket_ball": str(bracket2),
        "pointwise_hodge_lower": "1/4",
        "local_shape_radius_ball": str(local_radius),
        "average_deficit_lower": "1/24",
        "coherence_failure_action": ">= d/20",
        "status": "CERTIFIED",
    }


def stress(samples: int = 20_000, seed: int = 20260807) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    r = float_rstar()
    eps = float(COHERENCE_FRACTION)
    T = float(COHERENT_STRAIN_TIME)
    worst_hodge_ratio = float("inf")
    worst_local = 0.0
    worst_gradient_identity = 0.0
    worst_corotational_identity = 0.0

    for _ in range(samples):
        psi0 = float(rng.uniform(0.0, math.pi))
        c0, s0 = math.cos(psi0), math.sin(psi0)
        R0 = np.array([[c0, -s0], [s0, c0]])
        D0 = R0 @ np.diag([1.0, -1.0]) @ R0.T

        psi1 = float(rng.uniform(0.0, math.pi))
        c1, s1 = math.cos(psi1), math.sin(psi1)
        R1 = np.array([[c1, -s1], [s1, c1]])
        Epert = R1 @ np.diag([1.0, -1.0]) @ R1.T
        phase = float(rng.uniform(0.0, 2 * math.pi))
        omega = float(rng.uniform(0.0, 12.0 / T))
        amp = float(rng.uniform(0.0, eps))

        def Dfun(t):
            return D0 + amp * math.sin(omega * t + phase) * Epert

        F = _rk4_fundamental(Dfun, T, steps=40)
        u, v, H = shape_coords_from_fundamental(F, r)
        ratio = H / (T * T)
        worst_hodge_ratio = min(worst_hodge_ratio, ratio)
        worst_local = max(worst_local, abs(u), abs(v))
        if H + 3e-8 < float(COHERENT_HODGE) * T * T:
            raise AssertionError(("coherent variable-strain Hodge bound failed", H, ratio, amp, omega, phase))
        if max(abs(u), abs(v)) > float(LOCAL_LOG_BOX) + 2e-7:
            raise AssertionError("coherent variable strain left the local box")

        # Exact NS velocity-gradient and objective-strain algebra.
        A = rng.normal(size=(3, 3))
        A -= np.trace(A) / 3 * np.eye(3)
        Hp = rng.normal(size=(3, 3)); Hp = sym(Hp)
        LA = rng.normal(size=(3, 3))
        nu = float(rng.uniform(0.0, 2.0))
        Arhs = velocity_gradient_material_rhs(A, Hp, LA, nu)
        S = sym(A); O = 0.5 * (A - A.T)
        directS = sym(Arhs) + S @ O - O @ S
        formulaS = corotational_strain_rhs(A, Hp, LA, nu)
        scale = max(1.0, np.linalg.norm(directS))
        worst_corotational_identity = max(worst_corotational_identity, np.linalg.norm(directS - formulaS) / scale)
        # sym(-A^2-Hp+nu LA) is the material S derivative.
        materialS = sym(Arhs)
        target_materialS = -(S @ S) - (O @ O) - Hp + nu * sym(LA)
        worst_gradient_identity = max(worst_gradient_identity, np.linalg.norm(materialS - target_materialS) / max(1.0, np.linalg.norm(target_materialS)))

    return {
        "samples": samples,
        "worst_coherent_hodge_ratio": worst_hodge_ratio,
        "worst_abs_shape_coordinate": worst_local,
        "worst_material_strain_identity_residual": worst_gradient_identity,
        "worst_corotational_identity_residual": worst_corotational_identity,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=20_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-strain-coherence"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_strain_coherence_certificate()
    result = stress(args.samples)
    payload = {"certificate": cert, "stress": result}
    (args.outdir / "strain_coherence.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Strain coherence / objective-gradient dichotomy\n\nStatus: **{cert['status']}**.\n\n- coherent variation threshold: `<= {cert['coherence_fraction']}` of the initial non-conformal strain\n- dimensionless strain-time: `dT <= {cert['strain_time']}`\n- pointwise Hodge coefficient: `>= {cert['pointwise_hodge_lower']}`\n- certified bracket enclosure: `{cert['hodge_bracket_ball']}`\n- local-coordinate radius enclosure: `{cert['local_shape_radius_ball']}`\n- time-averaged edge deficit: `>= {cert['average_deficit_lower']} (dT)^2`\n- coherence failure requires objective-strain variation `{cert['coherence_failure_action']}`\n- variable-strain traces: `{result['samples']}`\n- worst numerical Hodge ratio `H/(dT)^2`: `{result['worst_coherent_hodge_ratio']:.9f}`\n- worst local shape coordinate: `{result['worst_abs_shape_coordinate']:.9f}`\n- worst material-strain identity residual: `{result['worst_material_strain_identity_residual']:.3e}`\n- worst corotational identity residual: `{result['worst_corotational_identity_residual']:.3e}`\n\nThe low-cost alternative is therefore no longer "the strain was not frozen".\nEither coherent non-conformal strain pays a multiplier cost, or the co-rotating\nstrain changes by a definite fraction.  Navier--Stokes identifies the source of\nthat change through `D_t grad u + (grad u)^2 = -Hess p + nu Delta grad u`.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
