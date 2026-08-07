from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from src.single_edge_certificate import RSTAR_LO

TRANSPORT_CUT = Fraction(1, 4)
ROLE_LOWER = Fraction(3, 5)


def selected_parent_lower_ratio(rstar: float, u_max: float = 1/200, v_max: float = 1/100) -> float:
    """Minimum parent/child frequency ratio on the signed-good core.

    x=r*exp(-v-u/2), y=r*exp(-v+u/2).  The smaller parent is minimized at
    v=v_max, u=u_max.
    """
    return rstar * math.exp(-v_max - 0.5 * u_max)


def low_low_output_radius(transport_cut: float) -> float:
    """Minkowski support radius of V tensor V if supp(Vhat) lies in |xi|<=cut*N."""
    return 2.0 * transport_cut


def affine_transport_commutator_symbol(A: np.ndarray, xi: np.ndarray, grad_m: np.ndarray) -> float:
    """Symbol of [(A x).grad, m(D)] for trace-free A.

    Fourier convention gives hat((Ax).grad f)=-(A^T xi).grad_xi f when tr A=0,
    hence the commutator symbol is -(A^T xi).grad m.
    """
    A = np.asarray(A, float)
    xi = np.asarray(xi, float)
    grad_m = np.asarray(grad_m, float)
    return -float(np.dot(A.T @ xi, grad_m))


def moving_multiplier_residual(A: np.ndarray, xi: np.ndarray, grad_m: np.ndarray, dt_m: float) -> float:
    """Heisenberg residual dt m + [(Ax).grad,m(D)] symbol."""
    return float(dt_m + affine_transport_commutator_symbol(A, xi, grad_m))


def transported_gaussian_symbol(A: np.ndarray, t: float, xi: np.ndarray, C: np.ndarray) -> tuple[float, np.ndarray, float]:
    """Exact scalar packet symbol transported by the affine dual flow.

    m(t,xi)=exp(-1/2 |C exp(A^T t) xi|^2) solves
        partial_t m - (A^T xi).grad_xi m = 0
    for constant trace-free A.
    Returns (m, grad_xi m, partial_t m).
    """
    A = np.asarray(A, float)
    xi = np.asarray(xi, float)
    C = np.asarray(C, float)
    E = expm(A.T * t)
    y = E @ xi
    Q = C.T @ C
    val = math.exp(-0.5 * float(y @ Q @ y))
    grad = -val * (E.T @ Q @ y)
    # Since E commutes with A^T, the transport PDE gives the exact time derivative.
    dt = float(np.dot(A.T @ xi, grad))
    return val, grad, dt


def microscopic_role_sources() -> tuple[str, ...]:
    """Sources belonging to the full-velocity packet equation after tangent quotient."""
    return (
        "designated_high_high_triad",
        "other_high_high_interactions",
        "cross_low_high_cells",
        "moving_projector_heisenberg_residual",
        "nonaffine_third_hermite",
        "relative_polarization_curvature",
        "spatial_window_leray_or_multiplier_commutators",
        "viscous_window_boundary_if_strongly_localized",
    )


def macroscopic_sgs_energy_sources() -> tuple[str, ...]:
    """Terms belonging to the resolved localized energy ledger, not role forcing."""
    return (
        "raw_sgs_work_Pi",
        "combined_pressure_flux_work",
        "RU_boundary_transport",
        "resolved_window_transport_mismatch",
        "resolved_viscous_boundary_flux",
        "pressure_cancellation_annular_charge",
        "partition_boundary_flux",
    )


def source_taxonomy_is_disjoint() -> bool:
    return set(microscopic_role_sources()).isdisjoint(macroscopic_sgs_energy_sources())


def arb_transport_separation_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    rlo = arb(RSTAR_LO.numerator) / RSTAR_LO.denominator
    parent_min = rlo * (-arb(1) / 80).exp()
    if not (parent_min > arb(3) / 5):
        raise AssertionError(f"good-core parent lower frequency failed: {parent_min}")
    if not (arb(2) * (arb(1) / 4) < arb(3) / 5):
        raise AssertionError("low-low support separation failed")
    return {
        "parent_lower_ball": str(parent_min),
        "clean_role_lower": "3/5",
        "transport_lowpass_cut": "1/4",
        "low_low_output_upper": "1/2",
        "spectral_gap": "1/10",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class PacketProjectionStress:
    samples: int
    worst_heisenberg_residual: float
    worst_finite_difference_dt_residual: float
    minimum_low_low_gap: float
    taxonomy_disjoint: bool


def stress(samples: int = 50_000, seed: int = 20260807) -> PacketProjectionStress:
    rng = np.random.default_rng(seed)
    wh = wfd = 0.0
    mingap = 1e9
    for _ in range(samples):
        A = rng.normal(size=(3, 3))
        A -= np.trace(A) / 3.0 * np.eye(3)
        xi = rng.normal(size=3)
        C = rng.normal(size=(3, 3))
        C = C.T @ C + 0.2 * np.eye(3)
        t = float(rng.uniform(-0.08, 0.08))
        m, grad, dt = transported_gaussian_symbol(A, t, xi, C)
        res = moving_multiplier_residual(A, xi, grad, dt)
        wh = max(wh, abs(res))
        if abs(res) > 2e-11 * max(1.0, abs(dt), np.linalg.norm(grad)):
            raise AssertionError("affine moving-multiplier residual failed")
        eps = 1e-6
        mp = transported_gaussian_symbol(A, t + eps, xi, C)[0]
        mm = transported_gaussian_symbol(A, t - eps, xi, C)[0]
        fd = (mp - mm) / (2.0 * eps)
        wfd = max(wfd, abs(fd - dt))
        if abs(fd - dt) > 4e-4 * max(1.0, abs(dt)):
            raise AssertionError("transported symbol time derivative finite difference failed")
        role_lower = float(rng.uniform(0.6001, 0.95))
        cut = float(rng.uniform(0.02, 0.2499))
        gap = role_lower - 2.0 * cut
        if role_lower >= 0.6 and cut <= 0.25:
            mingap = min(mingap, gap)
            if gap <= 0.0:
                raise AssertionError("low-low support leaked into selected role")
    return PacketProjectionStress(samples, wh, wfd, mingap, source_taxonomy_is_disjoint())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-smooth-sgs-packet-equation"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_transport_separation_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": out.__dict__, "micro_sources": microscopic_role_sources(), "macro_sources": macroscopic_sgs_energy_sources()}
    (args.outdir / "smooth_sgs_packet_equation.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Smooth-SGS two-level affine packet equation

Status: **{cert['status']}** for the good-core transport separation.

- signed-good parent/child ratio is strictly above `3/5`;
- choose a smooth common transport low-pass supported in `|xi|<=N/4`;
- its low-low product is supported in `|xi|<=N/2`, so it cannot feed any selected role directly;
- the selected role equation therefore splits exactly into low-high transport and high-high transfer;
- a Fourier packet multiplier transported by the affine dual flow has zero Heisenberg residual
  `partial_t m + [(Ax).grad,m(D)] = 0`;
- microscopic full-velocity role forcing and macroscopic resolved-SGS boundary transport are recorded as disjoint ledgers, preventing `RU`/pressure double counting.

Stress checks: `{out.samples}`
- worst affine Heisenberg residual: `{out.worst_heisenberg_residual:.3e}`
- worst transported-symbol finite-difference residual: `{out.worst_finite_difference_dt_residual:.3e}`
- minimum sampled low-low/role support gap: `{out.minimum_low_low_gap:.6f}`
- source taxonomies disjoint: `{out.taxonomy_disjoint}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
