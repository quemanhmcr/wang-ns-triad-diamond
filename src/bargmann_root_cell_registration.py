from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.dual_gaussian_root_registration import (
    SCALE_COLORS,
    covariance_cover_number_upper,
    dual_probe_critical_mass_lower,
    same_color_scale_shells_are_disjoint,
)

PHASE_REAL_DIMENSION = 6
BARGMANN_COMPLEX_DIMENSION = 3


def bargmann_ball_energy_fraction(radius: float, complex_dimension: int = BARGMANN_COMPLEX_DIMENSION) -> float:
    """Local Husimi/Moyal mass forced by one coherent coefficient.

    In normalized coherent coordinates with |<g_z,g_w>|=exp(-|z-w|^2/2),
      c(z)=exp(-|z|^2/2) F(z)
    after Weyl translation to z0=0, with F entire holomorphic in C^d.
    Mean-value/subharmonicity gives
      int_{B_R}|F|^2 pi^-d dz >= R^(2d)/d! |F(0)|^2.
    Since exp(-|z|^2)>=exp(-R^2) on the ball,
      E(B_R)>=exp(-R^2)R^(2d)/d! |c(0)|^2.
    """
    R = float(radius)
    d = int(complex_dimension)
    if R < 0 or d <= 0 or not math.isfinite(R):
        raise ValueError("finite nonnegative radius and positive dimension required")
    return math.exp(-R * R) * R ** (2 * d) / math.factorial(d)


def optimal_bargmann_radius(complex_dimension: int = BARGMANN_COMPLEX_DIMENSION) -> float:
    if complex_dimension <= 0:
        raise ValueError("positive complex dimension required")
    return math.sqrt(float(complex_dimension))


def optimal_bargmann_fraction(complex_dimension: int = BARGMANN_COMPLEX_DIMENSION) -> float:
    d = int(complex_dimension)
    if d <= 0:
        raise ValueError("positive complex dimension required")
    return math.exp(-d) * d**d / math.factorial(d)


def unit_grid_cells_intersecting_ball_upper(radius: float, real_dimension: int = PHASE_REAL_DIMENSION) -> int:
    """Axis-aligned unit cells intersecting a Euclidean ball, by box enclosure."""
    R = float(radius)
    d = int(real_dimension)
    if R < 0 or d <= 0 or not math.isfinite(R):
        raise ValueError("finite nonnegative radius and positive dimension required")
    per_axis = math.ceil(2.0 * R) + 1
    return per_axis**d


def canonical_cell_critical_mass_lower(
    probe_critical_mass: float,
    radius: float | None = None,
) -> float:
    """N E_C lower for one canonical unit coherent cell near the probe mark."""
    eta = float(probe_critical_mass)
    if eta < 0 or not math.isfinite(eta):
        raise ValueError("finite nonnegative probe critical mass required")
    R = optimal_bargmann_radius() if radius is None else float(radius)
    frac = bargmann_ball_energy_fraction(R)
    cells = unit_grid_cells_intersecting_ball_upper(R)
    return frac * eta / cells


def default_canonical_cell_quantum() -> float:
    return canonical_cell_critical_mass_lower(dual_probe_critical_mass_lower())


def canonical_cell_frame_budget() -> int:
    """Uniform energy budget for cell witnesses across scale/covariance strata.

    For one exact outer frequency/helicity subrole and one covariance
    representative, Moyal cell energies sum exactly to that role's L2 energy.
    Within one log-scale bin the selected frequency subroles form a disjoint
    partition.  Same-color scale bins are orthogonal.  Therefore only the four
    scale colors and finite covariance representatives multiply the global energy
    budget; no phase-space coloring or Gaussian Riesz constant is needed here.
    """
    if not same_color_scale_shells_are_disjoint():
        raise AssertionError("scale coloring must give orthogonal outer role subspaces")
    return SCALE_COLORS * covariance_cover_number_upper()


def registered_material_root_count_upper(
    global_energy: float,
    root_scale_upper: float,
    cell_quantum: float | None = None,
) -> float:
    """Bound distinct canonical energy-anchor root cells on a common root slice."""
    if global_energy < 0 or root_scale_upper <= 0:
        raise ValueError("valid energy/root scale required")
    eta = default_canonical_cell_quantum() if cell_quantum is None else float(cell_quantum)
    if eta <= 0:
        raise ValueError("positive cell quantum required")
    return canonical_cell_frame_budget() * global_energy * root_scale_upper / eta


def renyi_reuse_action_lower_from_material_cells(
    depth: int,
    global_energy: float,
    base_scale: float,
    cell_quantum: float | None = None,
) -> float:
    """Existing binary root action after replacing packet roots by energetic cells."""
    if depth < 0 or global_energy <= 0 or base_scale <= 0:
        raise ValueError("valid depth/energy/base scale required")
    root_scale_upper = base_scale * (25.0 / 24.0) ** depth
    n0 = registered_material_root_count_upper(global_energy, root_scale_upper, cell_quantum)
    return depth * math.log(2.0) - math.log(max(n0, 1.0))


def deterministic_energy_anchor(cell_energies: dict[tuple[int, ...], float]) -> tuple[int, ...]:
    """Choose an actual positive-energy canonical cell without arbitrary jitter.

    Maximize Moyal cell energy; break exact ties lexicographically.  The rule is
    deterministic and depends only on the physical analysis measure.
    """
    if not cell_energies:
        raise ValueError("nonempty cell-energy family required")
    for key, value in cell_energies.items():
        if not key or value < 0 or not math.isfinite(value):
            raise ValueError("valid canonical cell energies required")
    best_value = max(cell_energies.values())
    return min(key for key, value in cell_energies.items() if value == best_value)


def pushforward_parent_slot_weights(
    slot_weights: dict[str, float],
    parent_to_anchor: dict[str, tuple[int, ...]],
) -> dict[tuple[int, ...], float]:
    """Positive causal parent law pushed to energetic material-cell anchors."""
    out: dict[tuple[int, ...], float] = {}
    total = 0.0
    for parent, weight in slot_weights.items():
        if weight < 0 or not math.isfinite(weight):
            raise ValueError("finite nonnegative slot weights required")
        if parent not in parent_to_anchor:
            raise ValueError("every parent slot needs one energy anchor")
        key = parent_to_anchor[parent]
        out[key] = out.get(key, 0.0) + weight
        total += weight
    if abs(sum(out.values()) - total) > 2e-13 * max(1.0, total):
        raise AssertionError("parent-label pushforward changed positive causal mass")
    return out


@dataclass(frozen=True)
class BargmannRootCellStress:
    samples: int
    minimum_optimality_margin: float
    minimum_cell_quantum: float
    worst_pushforward_mass_residual: float
    canonical_ball_cells: int
    cell_frame_budget: int


def stress(samples: int = 50_000, seed: int = 20260808) -> BargmannRootCellStress:
    rng = np.random.default_rng(seed)
    mo = float("inf")
    wp = 0.0
    Ropt = optimal_bargmann_radius()
    fopt = optimal_bargmann_fraction()
    eta = default_canonical_cell_quantum()
    if eta <= 0:
        raise AssertionError("canonical coherent cell quantum is not positive")
    for _ in range(samples):
        R = float(rng.uniform(0.0, 5.0))
        margin = fopt - bargmann_ball_energy_fraction(R)
        mo = min(mo, margin)
        if margin < -3e-14:
            raise AssertionError("claimed Bargmann radius is not optimal")

        n = int(rng.integers(1, 30))
        weights = {f"p{j}": float(rng.random()) for j in range(n)}
        anchors = {f"p{j}": (int(rng.integers(-4, 5)), int(rng.integers(-4, 5))) for j in range(n)}
        pushed = pushforward_parent_slot_weights(weights, anchors)
        res = abs(sum(pushed.values()) - sum(weights.values()))
        wp = max(wp, res)
        if res > 3e-13 * max(1.0, sum(weights.values())):
            raise AssertionError("positive causal pushforward lost mass")

        # Deterministic maximum-energy anchor really carries at least average mass.
        m = int(rng.integers(1, 50))
        cells = {(j,): float(rng.random()) for j in range(m)}
        anchor = deterministic_energy_anchor(cells)
        if cells[anchor] + 1e-15 < sum(cells.values()) / m:
            raise AssertionError("max-energy anchor fell below the average")

    return BargmannRootCellStress(
        samples,
        mo,
        eta,
        wp,
        unit_grid_cells_intersecting_ball_upper(Ropt),
        canonical_cell_frame_budget(),
    )


def theorem_certificate() -> dict[str, object]:
    eta_probe = dual_probe_critical_mass_lower()
    R = optimal_bargmann_radius()
    frac = optimal_bargmann_fraction()
    cells = unit_grid_cells_intersecting_ball_upper(R)
    eta_cell = canonical_cell_critical_mass_lower(eta_probe, R)
    budget = canonical_cell_frame_budget()
    return {
        "status": "EXACT_BARGMANN_LOCAL_MOYAL_ROOT_CELL_QUANTUM__SELECTED_OUTER_ROLE_PARTITION_ASSUMED",
        "local_submean": f"E(B_R)>=exp(-R^2)R^6/3! |<f,g_z0>|^2; optimum R=sqrt3 gives fraction {frac:.12g}",
        "canonical_cell": f"B_sqrt3 intersects at most {cells} unit six-dimensional material cells, so one carries N E_C>={eta_cell:.12g}",
        "parent_label": "choose the maximum-energy canonical cell near the complex Gaussian mark; exact ties use a fixed lexicographic rule",
        "causal_pushforward": "positive parent-slot weights may be pushed to these energy anchors without changing total physical causal mass; collisions only increase reuse",
        "energy_budget": f"Moyal P=1 inside each exact outer role/covariance frame; 4 scale colors times finite covariance bins give uniform cell budget {budget}",
        "no_work_alignment_needed": "the anchor labels physical parent identity through actual energy; it is not asserted that nonlinear work is spatially localized in the same cell",
        "important_scope": "the outer frozen Fourier/helical roles must form the exact disjoint/orthogonal selected partition used in the global energy budget",
        "continuum_status": "with the companion complex-Young and dual-probe theorems, remaining assembly is to make this parent energy-anchor map on every recursive selected physical event and use it consistently as the canonical material label",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-bargmann-root-cell-registration"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "bargmann_root_cell_registration.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Bargmann root-cell registration: a large Gaussian coefficient forces actual Moyal cell energy\n\nStatus: **{cert['status']}**.\n\nThe remaining registration does not require the nonlinear work atom to sit in the same coherent cell as the Christ mark.  A causal parent needs a physical **identity anchor** with a global energy budget.  The Gaussian coherent transform supplies one.\n\nFor normalized coherent states with overlap `|<g_z,g_w>|=exp(-|z-w|^2/2)`, Weyl translation to a mark `z0` writes the coefficient function as\n\n`c(z)=exp(-|z|^2/2) F(z)`,\n\nwhere `F` is Bargmann holomorphic on `C^3`.  Mean-value/subharmonicity gives for every `R`\n\n`E(B_R(z0)) >= exp(-R^2) R^6/3! |c(z0)|^2`.\n\nThe optimum is `R=sqrt(3)`, with fraction\n\n`9/(2e^3) = {optimal_bargmann_fraction():.12g}`.\n\nThe radius-`sqrt3` ball intersects at most `{out.canonical_ball_cells}` unit cells of the fixed six-dimensional intrinsic grid.  Therefore one actual canonical coherent cell has\n\n`N E_C >= {out.minimum_cell_quantum:.12g}`\n\nwhen the dual-Gaussian parent coefficient uses the default one-percent/`delta=0.4` quantum.  The number is small but **strictly scale independent**; causal reuse only needs a fixed positive root quantum.\n\nChoose the maximum-Moyal-energy canonical cell in that finite neighborhood, with deterministic tie breaking, and call it the parent energy anchor.  Push the positive physical parent-slot law to these anchors.  Total positive transfer mass is unchanged.  If several slots choose the same anchor they are physically merged/reused; if they choose distinct anchors, each distinct root has a positive Moyal energy quantum.\n\nThere is no need to claim that the nonlinear work itself is localized in the anchor cell.  Work weights come from the exact physical child-transfer measure; the anchor answers a different physical question: **which material coherent reservoir carries this parent role?**\n\nThe global root budget is also depth independent.  Inside one exact outer frequency/helicity subrole and one covariance representative, Moyal cell energies sum with `P=1`.  The selected relative frequency roles within a scale bin are disjoint; same-color logarithmic scale bins are orthogonal; only four scale colors and the finite covariance representatives multiply the energy budget.  Thus the effective cell budget is `{out.cell_frame_budget}` with no packet-count or causal-depth factor.\n\nStress: `{out.samples}` radius/pushforward/anchor states\n- minimum optimal-radius margin: `{out.minimum_optimality_margin:.3e}`\n- canonical root-cell critical mass: `{out.minimum_cell_quantum:.12g}`\n- worst positive pushforward mass residual: `{out.worst_pushforward_mass_residual:.3e}`\n\nCombined with complex Young parent marking and dual-Gaussian analysis, this replaces the old transfer-cell-alignment demand by a simpler physical registration: each selected parent role carries its own nearby energetic material cell, and the causal law is pushed to that cell.  The remaining continuum step is to verify that the recursive physical selector really is an exact disjoint/orthogonal outer role partition and to install this energy-anchor rule as the one canonical material label used by Duhamel/Renyi/Hodge/service bookkeeping.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
