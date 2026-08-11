from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.canonical_positive_edge_work_routing import (
    compress_signed_edge_work_to_hard_cells,
    exact_mode_role_map,
    route_canonical_positive_edge_work,
    single_hard_role_map,
)
from src.continuum_helical_edge_measure_pde_probe import (
    CHILD_MODE,
    _deterministic_smooth_initial_state,
    _nonlinear_term,
    _pair_orbits_for_child,
    _rk4_step,
    _series_coefficient,
    _snapshot,
    _spectral_geometry,
)
from src.continuum_helical_edge_measure_registration import (
    continuum_edge_measure_ledger,
    register_continuum_triad_fiber,
    unitary_fourier_convolution_factor,
)

STATUS = (
    "ACTUAL_GALERKIN_NS_CANONICAL_POSITIVE_EDGE_ROUTING_AUDIT__"
    "SIGNED_EDGE_RECONSTRUCTION__DW_PLUS_FATE__HARD_PUSHFORWARD_NO_REHAHN"
)


def _ledger_from_actual_state(
    state_hat: np.ndarray,
    *,
    cutoff: int,
    child_mode: tuple[int, int, int] = CHILD_MODE,
):
    n = int(state_hat.shape[1])
    child = tuple(int(v) for v in child_mode)
    z = np.asarray(child, dtype=float)
    uz = _series_coefficient(state_hat, child)
    qmass = 1.0 / unitary_fourier_convolution_factor()
    fibers = []
    for x, y in _pair_orbits_for_child(child, cutoff):
        fibers.append(
            register_continuum_triad_fiber(
                x=np.asarray(x, dtype=float),
                y=np.asarray(y, dtype=float),
                z=z,
                ux=_series_coefficient(state_hat, x),
                uy=_series_coefficient(state_hat, y),
                uz=uz,
                quotient_measure_mass=qmass,
            )
        )
    if not fibers:
        raise AssertionError("actual NS routing audit found no retained parent orbit")
    return continuum_edge_measure_ledger(tuple(fibers))


@dataclass(frozen=True)
class CanonicalPositiveEdgeRoutingPDEProbe:
    status: str
    resolution: int
    cutoff: int
    steps: int
    snapshots: int
    positive_work_snapshots: int
    bad_work_snapshots: int
    positive_nonforward_snapshots: int
    worst_signed_ns_reconstruction_relative: float
    worst_positive_mass_reconstruction_relative: float
    worst_hard_pushforward_relative: float
    maximum_coarsened_cancellation_fraction: float
    minimum_bad_deficit_margin: float | None
    minimum_bad_fixed_transfer_margin: float | None
    stage_zero_first_time_failures: int
    geometry_good_marking_promotions: int


def run_probe(
    *,
    resolution: int = 24,
    steps: int = 48,
    viscosity: float = 0.03,
    amplitude: float = 4.0,
    duration: float = 0.008,
    snapshot_count: int = 5,
    tau: float = 0.1,
) -> CanonicalPositiveEdgeRoutingPDEProbe:
    n = int(resolution)
    count = int(steps)
    snaps = int(snapshot_count)
    if count < 16 or snaps < 3 or snaps > count + 1:
        raise ValueError("actual NS routing audit needs at least sixteen steps and three snapshots")
    k, k2, dealias, cutoff = _spectral_geometry(n, None)
    state = _deterministic_smooth_initial_state(n, k, k2, dealias, float(amplitude))
    dt = float(duration) / count
    sample_indices = tuple(sorted({round(j * count / (snaps - 1)) for j in range(snaps)}))

    positive_snapshots = 0
    bad_snapshots = 0
    nonforward_snapshots = 0
    worst_signed = 0.0
    worst_mass = 0.0
    worst_push = 0.0
    max_cancel = 0.0
    min_bad = math.inf
    min_fixed = math.inf
    first_time_failures = 0
    marking_promotions = 0

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        if step in sample_indices:
            row = _snapshot(state, k, k2, dealias, cutoff, nonlinear_hat=nonlinear)
            ledger = _ledger_from_actual_state(state, cutoff=cutoff)
            exact_roles = exact_mode_role_map(ledger)
            route = route_canonical_positive_edge_work(ledger, tau=tau, mode_roles=exact_roles)
            coarse = compress_signed_edge_work_to_hard_cells(ledger, single_hard_role_map(ledger))

            signed_scale = max(
                abs(float(row["ledger_signed_direct_work"])),
                abs(ledger.signed_direct_work),
                ledger.positive_edge_work + ledger.negative_edge_work,
                1.0e-300,
            )
            worst_signed = max(
                worst_signed,
                abs(float(row["ledger_signed_direct_work"]) - ledger.signed_direct_work) / signed_scale,
            )
            mass_scale = max(route.total_positive_work, 1.0e-300)
            worst_mass = max(worst_mass, abs(route.mass_reconstruction_residual) / mass_scale)
            worst_push = max(
                worst_push,
                abs(route.hard_cell_compression.inherited_positive_work - route.total_positive_work) / mass_scale,
            )
            if coarse.inherited_positive_work > 0.0:
                max_cancel = max(max_cancel, coarse.cancellation_gap / coarse.inherited_positive_work)

            positive_snapshots += route.total_positive_work > 0.0
            bad_snapshots += route.bad_positive_work > 0.0
            nonforward_snapshots += ledger.positive_nonforward_work > 0.0
            marking_promotions += int(route.young_eligible.marking_good or route.young_eligible.young_certified)
            if route.bad_route is not None:
                min_bad = min(min_bad, route.bad_route.deficit - 1.0e-4)
                min_fixed = min(
                    min_fixed,
                    route.bad_route.deficit - route.bad_route.fixed_transfer_gate.threshold,
                )
                first_time_failures += int(route.bad_route.transfer_partition.first_time is not None)
                if ledger.positive_nonforward_work > 0.0 and route.bad_positive_work <= 0.0:
                    raise AssertionError("actual positive nonforward NS work was dropped instead of routed bad")

            if abs(route.total_positive_work - ledger.positive_edge_work) > 8.0e-10 * max(
                route.total_positive_work, ledger.positive_edge_work, 1.0e-300
            ):
                raise AssertionError("actual NS routing changed canonical dW+ mass")
            if abs(coarse.inherited_positive_work - ledger.positive_edge_work) > 8.0e-10 * max(
                coarse.inherited_positive_work, ledger.positive_edge_work, 1.0e-300
            ):
                raise AssertionError("actual NS hard coarsening changed inherited causal mass")
        if step < count:
            state = _rk4_step(state, dt, float(viscosity), k, k2, dealias)

    if positive_snapshots == 0:
        raise AssertionError("actual NS routing audit saw no positive canonical child work")
    if bad_snapshots == 0:
        raise AssertionError("actual NS routing audit saw no geometry-bad positive work")
    if first_time_failures or marking_promotions:
        raise AssertionError("actual NS audit exposed a stage-zero clock or geometry-to-marking promotion")

    return CanonicalPositiveEdgeRoutingPDEProbe(
        status=STATUS,
        resolution=n,
        cutoff=cutoff,
        steps=count,
        snapshots=len(sample_indices),
        positive_work_snapshots=positive_snapshots,
        bad_work_snapshots=bad_snapshots,
        positive_nonforward_snapshots=nonforward_snapshots,
        worst_signed_ns_reconstruction_relative=worst_signed,
        worst_positive_mass_reconstruction_relative=worst_mass,
        worst_hard_pushforward_relative=worst_push,
        maximum_coarsened_cancellation_fraction=max_cancel,
        minimum_bad_deficit_margin=None if math.isinf(min_bad) else min_bad,
        minimum_bad_fixed_transfer_margin=None if math.isinf(min_fixed) else min_fixed,
        stage_zero_first_time_failures=first_time_failures,
        geometry_good_marking_promotions=marking_promotions,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=24)
    parser.add_argument("--steps", type=int, default=48)
    parser.add_argument("--viscosity", type=float, default=0.03)
    parser.add_argument("--amplitude", type=float, default=4.0)
    parser.add_argument("--duration", type=float, default=0.008)
    parser.add_argument("--snapshots", type=int, default=5)
    parser.add_argument("--tau", type=float, default=0.1)
    parser.add_argument("--outdir", type=Path, default=Path("results-canonical-positive-edge-work-routing-pde"))
    args = parser.parse_args()
    out = run_probe(
        resolution=args.resolution,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        snapshot_count=args.snapshots,
        tau=args.tau,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "probe.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    summary = f"""# Actual Galerkin NS canonical positive edge routing audit\n\nStatus: **{STATUS}**.\n\n- resolution/cutoff: `{out.resolution}` / `{out.cutoff}`\n- steps/snapshots: `{out.steps}` / `{out.snapshots}`\n- positive-work snapshots: `{out.positive_work_snapshots}`\n- bad-work snapshots: `{out.bad_work_snapshots}`\n- positive-nonforward snapshots: `{out.positive_nonforward_snapshots}`\n- worst signed NS reconstruction relative residual: `{out.worst_signed_ns_reconstruction_relative:.12g}`\n- worst dW+ fate reconstruction relative residual: `{out.worst_positive_mass_reconstruction_relative:.12g}`\n- worst hard pushforward relative residual: `{out.worst_hard_pushforward_relative:.12g}`\n- maximum coarse cancellation fraction: `{out.maximum_coarsened_cancellation_fraction:.12g}`\n- minimum bad deficit margin above eta0: `{out.minimum_bad_deficit_margin}`\n- minimum bad deficit margin above fixed-transfer threshold: `{out.minimum_bad_fixed_transfer_margin}`\n- stage-zero first-time failures: `{out.stage_zero_first_time_failures}`\n- geometry-good marking promotions: `{out.geometry_good_marking_promotions}`\n\nThis audit evolves dealiased Fourier--Galerkin incompressible Navier--Stokes and routes the replayed physical edge law itself.  It does not infer causality from capacity or from a later coherent Hahn split, and it makes no global-regularity claim.\n"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
