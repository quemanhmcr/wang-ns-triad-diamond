from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from src.helical import edge_metrics, helical_basis
from src.helical_spin_transport import forward_normal_coupling, helical_with_normal, rotation_matrix, triad_normal

Array = np.ndarray


def wrap_angle(x: float) -> float:
    return float((x + math.pi) % (2.0 * math.pi) - math.pi)


def phase_line_ratio(reference: Array, target: Array) -> float:
    """Phase chi when target=e^{i chi} reference for unit vectors on one line."""
    z = np.vdot(reference, target)
    if abs(z) < 1e-12:
        raise ValueError("vectors are not on the same complex line")
    return float(np.angle(z))


def forward_edge_data(x: Array, y: Array, z: Array, sx: int, sy: int, sz: int) -> dict[str, float]:
    """Global and triad-normal phases for a forward edge x+y=z."""
    e = edge_metrics(x, y, z, sx, sy, sz)
    n = triad_normal(x, y)
    hx_edge = helical_with_normal(x, sx, n)
    hy_edge = helical_with_normal(y, sy, n)
    hz_edge = helical_with_normal(z, sz, n)
    chi_x = phase_line_ratio(helical_basis(x, sx), hx_edge)
    chi_y = phase_line_ratio(helical_basis(y, sy), hy_edge)
    chi_z = phase_line_ratio(helical_basis(z, sz), hz_edge)
    gamma = float(np.angle(forward_normal_coupling(x, y, z, sx, sy, sz)))
    reconstructed = wrap_angle(gamma + chi_x + chi_y - chi_z)
    return {
        "global_phase": e.g_phase,
        "normal_phase": gamma,
        "chi_x": chi_x,
        "chi_y": chi_y,
        "chi_z": chi_z,
        "reconstructed_phase": reconstructed,
        "target_phase": e.target_phase,
        "g_abs": e.g_abs,
    }


def diamond_edge_data(a: Array, b: Array, c: Array, signs: Sequence[int]) -> dict[str, dict[str, float]]:
    sa, sb, sc, sm, sn, sd = (int(x) for x in signs)
    m = np.asarray(a) + np.asarray(b)
    n = np.asarray(b) + np.asarray(c)
    d = np.asarray(a) + np.asarray(b) + np.asarray(c)
    return {
        "ab_m": forward_edge_data(a, b, m, sa, sb, sm),
        "mc_d": forward_edge_data(m, c, d, sm, sc, sd),
        "bc_n": forward_edge_data(b, c, n, sb, sc, sn),
        "an_d": forward_edge_data(a, n, d, sa, sn, sd),
    }


def diamond_holonomy_from_edges(edges: Mapping[str, Mapping[str, float]]) -> tuple[float, float, float]:
    geom = wrap_angle(
        edges["ab_m"]["global_phase"] + edges["mc_d"]["global_phase"]
        - edges["bc_n"]["global_phase"] - edges["an_d"]["global_phase"]
    )
    target = wrap_angle(
        edges["ab_m"]["target_phase"] + edges["mc_d"]["target_phase"]
        - edges["bc_n"]["target_phase"] - edges["an_d"]["target_phase"]
    )
    return geom, target, wrap_angle(geom - target)


def diamond_phase_residuals(
    a: Array,
    b: Array,
    c: Array,
    signs: Sequence[int],
    modal_phases: Mapping[str, float],
) -> dict[str, float]:
    """Gauge-invariant forward-transfer residuals in a four-edge diamond.

    delta_e=arg g_e - theta_parent1 - theta_parent2 + theta_child - target_e.
    The signed sum of the four deltas is independent of all six modal phases.
    """
    edges = diamond_edge_data(a, b, c, signs)
    th = modal_phases
    return {
        "ab_m": wrap_angle(edges["ab_m"]["global_phase"] - th["a"] - th["b"] + th["m"] - edges["ab_m"]["target_phase"]),
        "mc_d": wrap_angle(edges["mc_d"]["global_phase"] - th["m"] - th["c"] + th["d"] - edges["mc_d"]["target_phase"]),
        "bc_n": wrap_angle(edges["bc_n"]["global_phase"] - th["b"] - th["c"] + th["n"] - edges["bc_n"]["target_phase"]),
        "an_d": wrap_angle(edges["an_d"]["global_phase"] - th["a"] - th["n"] + th["d"] - edges["an_d"]["target_phase"]),
    }


def residual_holonomy(residuals: Mapping[str, float]) -> float:
    return wrap_angle(residuals["ab_m"] + residuals["mc_d"] - residuals["bc_n"] - residuals["an_d"])


def sharp_four_phase_cost(holonomy: float) -> float:
    """Sharp minimum of sum_i (1-cos delta_i) at four-edge holonomy H.

    For principal |H|<=pi the global maximizer of sum cos(delta_i), subject to
    delta1+delta2-delta3-delta4=H mod 2pi, is the equal lift
    (+H/4,+H/4,-H/4,-H/4).  Thus the exact minimum cost is
    4(1-cos(|H|/4)).
    """
    H = abs(wrap_angle(holonomy))
    return 4.0 * (1.0 - math.cos(H / 4.0))


def weighted_phase_cost_lower(holonomy: float, edge_weight_multiplier_floor: float) -> float:
    """Lower bound sum w_i m_i(1-cos delta_i) if every w_i m_i>=floor."""
    if edge_weight_multiplier_floor < 0:
        raise ValueError("floor must be nonnegative")
    return edge_weight_multiplier_floor * sharp_four_phase_cost(holonomy)


def arb_phase_holonomy_certificate() -> dict[str, str]:
    """Clean finite-packet branch: |H|>=1/5 and w_i>=beta, m_i>=1-1e-4.

    Then phase polarization deficit >= beta/250.
    """
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required for the phase certificate") from exc
    ctx.prec = 160
    eta = arb(1) / 10000
    H0 = arb(1) / 5
    coefficient = 4 * (1 - eta) * (1 - (H0 / 4).cos())
    if not (coefficient > arb(1) / 250):
        raise AssertionError(f"phase holonomy coefficient failed: {coefficient}")
    return {
        "holonomy_threshold": "1/5",
        "good_multiplier_floor": "1-1/10000",
        "phase_deficit_per_edge_weight_floor_ball": str(coefficient),
        "clean_lower_bound": "1/250",
        "statement": "if each diamond edge has capacity weight >= beta, phase deficit >= beta/250",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class HolonomyStress:
    samples: int
    worst_edge_normal_reconstruction: float
    worst_modal_cancellation_residual: float
    worst_rotation_holonomy_residual: float
    minimum_sharp_cost_margin: float


def _good_random_diamond(rng: np.random.Generator) -> tuple[Array, Array, Array]:
    for _ in range(100):
        a = rng.normal(size=3)
        b = rng.normal(size=3)
        c = rng.normal(size=3)
        m, n, d = a + b, b + c, a + b + c
        pairs = [(a, b), (m, c), (b, c), (a, n)]
        if min(np.linalg.norm(np.cross(x, y)) / max(1e-12, np.linalg.norm(x) * np.linalg.norm(y)) for x, y in pairs) < 0.08:
            continue
        if min(np.linalg.norm(x) for x in (a, b, c, m, n, d)) < 0.1:
            continue
        return a, b, c
    raise RuntimeError("failed to sample nondegenerate diamond")


def stress(samples: int = 50_000, seed: int = 20260807) -> HolonomyStress:
    rng = np.random.default_rng(seed)
    worst_recon = worst_cancel = worst_rot = 0.0
    min_margin = float("inf")
    accepted = 0
    while accepted < samples:
        a, b, c = _good_random_diamond(rng)
        signs = tuple(int(x) for x in rng.choice([-1, 1], size=6))
        try:
            edges = diamond_edge_data(a, b, c, signs)
        except ValueError:
            continue
        if min(e["g_abs"] for e in edges.values()) < 1e-7:
            continue
        accepted += 1
        for e in edges.values():
            worst_recon = max(worst_recon, abs(wrap_angle(e["global_phase"] - e["reconstructed_phase"])))
        geom, target, H = diamond_holonomy_from_edges(edges)

        phases = {name: float(rng.uniform(-math.pi, math.pi)) for name in ("a", "b", "c", "m", "n", "d")}
        residuals = diamond_phase_residuals(a, b, c, signs, phases)
        worst_cancel = max(worst_cancel, abs(wrap_angle(residual_holonomy(residuals) - H)))

        # Rigid rotation may cross deterministic basis charts, but the diamond
        # holonomy is per-mode gauge invariant and must not change.
        R = rotation_matrix(rng.normal(size=3), float(rng.uniform(-math.pi, math.pi)))
        edges_R = diamond_edge_data(R @ a, R @ b, R @ c, signs)
        geom_R, target_R, H_R = diamond_holonomy_from_edges(edges_R)
        worst_rot = max(worst_rot, abs(wrap_angle(H_R - H)), abs(wrap_angle(target_R - target)))

        # Adversarial phase residuals: exact sharp cost inequality.
        ds = rng.uniform(-math.pi, math.pi, size=4)
        HH = wrap_angle(float(ds[0] + ds[1] - ds[2] - ds[3]))
        actual = float(np.sum(1.0 - np.cos(ds)))
        lower = sharp_four_phase_cost(HH)
        margin = actual - lower
        min_margin = min(min_margin, margin)
        if margin < -2e-11:
            raise AssertionError(f"sharp four-phase cost violated: {margin}")

    return HolonomyStress(accepted, worst_recon, worst_cancel, worst_rot, min_margin)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-helical-holonomy"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_phase_holonomy_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": out.__dict__}
    (args.outdir / "helical_phase_holonomy.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Helical phase holonomy / diamond obstruction

Status: **{cert['status']}** for the clean finite-packet phase branch.

- exact diamond identity: modal phases cancel from the signed four-edge residual sum
- sharp phase cost: `sum(1-cos delta_i) >= 4(1-cos(|H|/4))`
- clean branch: `|H|>=1/5`, each edge capacity weight `>= beta`, each multiplier `>=1-1e-4`
  implies total polarization deficit `>= beta/250`
- certified coefficient before multiplying by beta: `{cert['phase_deficit_per_edge_weight_floor_ball']}`
- random nondegenerate diamonds: `{out.samples}`
- worst triad-normal/global coupling reconstruction residual: `{out.worst_edge_normal_reconstruction:.3e}`
- worst modal-phase cancellation residual: `{out.worst_modal_cancellation_residual:.3e}`
- worst rigid-rotation holonomy residual: `{out.worst_rotation_holonomy_residual:.3e}`
- minimum numerical sharp-cost margin: `{out.minimum_sharp_cost_margin:.3e}`

The observable is relative incidence holonomy, not the Berry phase of one mode or
one triad.  Rigid rotation is exactly free.  A nonzero diamond holonomy is a
phase-lock obstruction and therefore feeds the existing positive polarization
deficit.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
