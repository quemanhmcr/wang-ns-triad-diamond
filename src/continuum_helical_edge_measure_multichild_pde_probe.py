from __future__ import annotations

import argparse
import json
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from src.continuum_helical_edge_measure_pde_probe import (
    _deterministic_smooth_initial_state,
    _divergence_norm,
    _gradient_energy,
    _nonlinear_term,
    _relative_scalar,
    _rk4_step,
    _snapshot,
    _spectral_average_inner,
    _spectral_geometry,
    _trapezoid,
)

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "MULTIPLE_CHILD_MODES_ON_ONE_PDE_ORBIT__"
    "OUTER_CHILD_AGGREGATION_OF_UNORDERED_HELICAL_EDGE_WORK_BEFORE_HAHN"
)

MULTI_CHILD_MODES: tuple[tuple[int, int, int], ...] = (
    (5, 1, 0),
    (3, 4, 1),
    (1, 5, 2),
)


def _joint_snapshot(rows: Sequence[dict[str, float]]) -> dict[str, float]:
    rs = tuple(rows)
    if not rs:
        raise ValueError("at least one child registration row required")
    actual = math.fsum(float(row["actual_child_work"]) for row in rs)
    direct = math.fsum(float(row["ledger_signed_direct_work"]) for row in rs)
    modal = math.fsum(float(row["ledger_signed_modal_work"]) for row in rs)
    positive = math.fsum(float(row["positive_edge_work"]) for row in rs)
    negative = math.fsum(float(row["negative_edge_work"]) for row in rs)
    direct_progress = math.fsum(float(row["direct_progress"]) for row in rs)
    registered_progress = math.fsum(float(row["registered_progress"]) for row in rs)
    work_scale = max(
        abs(actual), abs(direct), abs(modal), positive + negative,
        math.fsum(abs(float(row["actual_child_work"])) for row in rs),
    )
    progress_scale = max(
        abs(direct_progress), abs(registered_progress),
        math.fsum(abs(float(row["direct_progress"])) for row in rs),
    )
    aggregate_positive = math.fsum(max(float(row["actual_child_work"]), 0.0) for row in rs)
    dominance_scale = max(positive, aggregate_positive)
    dominance_margin = positive - aggregate_positive
    if dominance_margin < -3.0e-10 * dominance_scale:
        raise AssertionError("multi-child positive edge Hahn law failed outer-child work dominance")
    return {
        "work_residual": _relative_scalar(actual, direct, work_scale),
        "modal_work_residual": _relative_scalar(actual, modal, work_scale),
        "hahn_residual": _relative_scalar(positive - negative, actual, work_scale),
        "progress_residual": _relative_scalar(direct_progress, registered_progress, progress_scale),
        "positive_dominance_margin": dominance_margin,
    }


@dataclass(frozen=True)
class MultiChildGalerkinEdgeRun:
    resolution: int
    spectral_cutoff: int
    child_modes: tuple[tuple[int, int, int], ...]
    steps: int
    snapshots: int
    child_snapshot_registrations: int
    minimum_unordered_pairs_per_child: int
    minimum_modal_edges_per_child: int
    worst_child_source_residual: float
    worst_child_work_residual: float
    worst_child_progress_residual: float
    worst_child_hahn_residual: float
    worst_joint_work_residual: float
    worst_joint_modal_work_residual: float
    worst_joint_hahn_residual: float
    worst_joint_progress_residual: float
    minimum_joint_positive_dominance_margin: float
    positive_nonforward_child_snapshots: int
    maximum_divergence_relative_to_initial_l2: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float


def simulate_multichild_edge_measure_on_galerkin_ns(
    *,
    resolution: int = 24,
    spectral_cutoff: int = 5,
    child_modes: Sequence[tuple[int, int, int]] = MULTI_CHILD_MODES,
    steps: int = 32,
    viscosity: float = 0.03,
    amplitude: float = 4.0,
    duration: float = 0.006,
    snapshot_count: int = 3,
) -> MultiChildGalerkinEdgeRun:
    n = int(resolution)
    cutoff = int(spectral_cutoff)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    snaps = int(snapshot_count)
    children = tuple(tuple(int(v) for v in child) for child in child_modes)
    if count < 12 or snaps < 2 or snaps > count + 1:
        raise ValueError("at least twelve RK4 steps and two snapshots required")
    if not children or len(set(children)) != len(children):
        raise ValueError("distinct nonempty child-mode family required")
    if any(child == (0, 0, 0) or max(abs(v) for v in child) > cutoff for child in children):
        raise ValueError("every child mode must lie inside the common Galerkin cutoff")
    if any(not any(v % 2 for v in child) for child in children):
        raise ValueError("audit children must avoid the discrete fixed parent orbit x=y")
    if not all(math.isfinite(x) and x > 0.0 for x in (nu, amp, horizon)):
        raise ValueError("positive finite physical PDE parameters required")

    k, k2, dealias, actual_cutoff = _spectral_geometry(n, cutoff)
    if actual_cutoff != cutoff:
        raise AssertionError("requested multi-child cutoff changed under FFT representation")
    state = _deterministic_smooth_initial_state(n, k, k2, dealias, amp)
    dt = horizon / count
    sample_indices = tuple(sorted({round(j * count / (snaps - 1)) for j in range(snaps)}))
    all_energy: list[float] = []
    all_gradient: list[float] = []
    all_nonlinear: list[float] = []
    all_divergence: list[float] = []
    child_rows: list[dict[str, float]] = []
    joint_rows: list[dict[str, float]] = []

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        all_energy.append(_spectral_average_inner(state, state, n))
        all_gradient.append(_gradient_energy(state, k2, n))
        all_nonlinear.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        all_divergence.append(_divergence_norm(state, k, n))
        if step in sample_indices:
            rows = tuple(
                _snapshot(
                    state, k, k2, dealias, cutoff,
                    child_mode=child, nonlinear_hat=nonlinear,
                )
                for child in children
            )
            if any(float(row["actual_source_norm"]) == 0.0 for row in rows):
                raise AssertionError("one child lost its actual nonlinear source on the common PDE orbit")
            child_rows.extend(rows)
            joint_rows.append(_joint_snapshot(rows))
        if step < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    names = ("unordered_source_residual", "signed_work_residual", "progress_residual", "hahn_residual")
    worst_child = {name: max(float(row[name]) for row in child_rows) for name in names}
    worst_joint_work = max(float(row["work_residual"]) for row in joint_rows)
    worst_joint_modal = max(float(row["modal_work_residual"]) for row in joint_rows)
    worst_joint_hahn = max(float(row["hahn_residual"]) for row in joint_rows)
    worst_joint_progress = max(float(row["progress_residual"]) for row in joint_rows)
    if max((*worst_child.values(), worst_joint_work, worst_joint_modal, worst_joint_hahn, worst_joint_progress)) > 3.0e-8:
        raise AssertionError("multi-child edge law lost actual Galerkin NS work/progress/Hahn identity")

    grid_times = tuple(i * dt for i in range(count + 1))
    initial_energy = all_energy[0]
    gradient_action = _trapezoid(all_gradient, grid_times)
    balance = abs(all_energy[-1] - initial_energy + 2.0 * nu * gradient_action) / initial_energy
    nonlinear_scale = initial_energy / horizon
    max_nonlinear = max(abs(x) for x in all_nonlinear) / nonlinear_scale
    max_divergence = max(all_divergence) / math.sqrt(initial_energy)
    if balance > 5.0e-5 or max_nonlinear > 5.0e-10 or max_divergence > 5.0e-11:
        raise AssertionError("multi-child audit trajectory lost a native Navier-Stokes invariant")

    return MultiChildGalerkinEdgeRun(
        resolution=n,
        spectral_cutoff=cutoff,
        child_modes=children,
        steps=count,
        snapshots=len(sample_indices),
        child_snapshot_registrations=len(child_rows),
        minimum_unordered_pairs_per_child=min(int(row["unordered_pairs"]) for row in child_rows),
        minimum_modal_edges_per_child=min(int(row["modal_edges"]) for row in child_rows),
        worst_child_source_residual=worst_child["unordered_source_residual"],
        worst_child_work_residual=worst_child["signed_work_residual"],
        worst_child_progress_residual=worst_child["progress_residual"],
        worst_child_hahn_residual=worst_child["hahn_residual"],
        worst_joint_work_residual=worst_joint_work,
        worst_joint_modal_work_residual=worst_joint_modal,
        worst_joint_hahn_residual=worst_joint_hahn,
        worst_joint_progress_residual=worst_joint_progress,
        minimum_joint_positive_dominance_margin=min(float(row["positive_dominance_margin"]) for row in joint_rows),
        positive_nonforward_child_snapshots=sum(float(row["positive_nonforward_work"]) > 0.0 for row in child_rows),
        maximum_divergence_relative_to_initial_l2=max_divergence,
        global_energy_balance_relative_residual=balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolution", type=int, default=24)
    ap.add_argument("--cutoff", type=int, default=5)
    ap.add_argument("--steps", type=int, default=32)
    ap.add_argument("--viscosity", type=float, default=0.03)
    ap.add_argument("--amplitude", type=float, default=4.0)
    ap.add_argument("--duration", type=float, default=0.006)
    ap.add_argument("--snapshots", type=int, default=3)
    ap.add_argument("--outdir", type=Path, default=Path("results-continuum-helical-edge-measure-multichild-pde"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = simulate_multichild_edge_measure_on_galerkin_ns(
        resolution=args.resolution,
        spectral_cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        snapshot_count=args.snapshots,
    )
    (args.outdir / "continuum_helical_edge_measure_multichild_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    summary = f"""# Multi-child continuum edge measure on one actual Galerkin Navier--Stokes orbit

Resolution `{result.resolution}`, common cutoff `{result.spectral_cutoff}`, child modes `{result.child_modes}`, `{result.steps}` RK4 steps, `{result.snapshots}` sampled PDE times.

Every child is reconstructed from all retained ordered parents, then unordered parent orbits and all eight helicity sectors. The outer child family is aggregated **before** the Hahn identity is checked on the same evolved PDE state.

- child registrations: `{result.child_snapshot_registrations}`;
- minimum unordered pairs/helical edges per child: `{result.minimum_unordered_pairs_per_child}` / `{result.minimum_modal_edges_per_child}`;
- worst child source/work/progress/Hahn residuals: `{result.worst_child_source_residual:.3e}`, `{result.worst_child_work_residual:.3e}`, `{result.worst_child_progress_residual:.3e}`, `{result.worst_child_hahn_residual:.3e}`;
- worst outer-child work/modal/Hahn/progress residuals: `{result.worst_joint_work_residual:.3e}`, `{result.worst_joint_modal_work_residual:.3e}`, `{result.worst_joint_hahn_residual:.3e}`, `{result.worst_joint_progress_residual:.3e}`;
- minimum outer-child positive-Hahn dominance margin: `{result.minimum_joint_positive_dominance_margin:.3e}`;
- child snapshots with positive nonforward physical work: `{result.positive_nonforward_child_snapshots}`;
- NS energy-balance residual: `{result.global_energy_balance_relative_residual:.3e}`;
- maximum normalized divergence: `{result.maximum_divergence_relative_to_initial_l2:.3e}`.

This is a finite Galerkin outer-child aggregation falsifier of the companion analytic joint Radon theorem `dLambda_unord=(1/16) dz d(q_#dr)`.  The PDE probe does not replace that proof; it attacks its normalization, parent quotient, helicity reconstruction, and Hahn ordering on one evolved Navier--Stokes orbit.
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
