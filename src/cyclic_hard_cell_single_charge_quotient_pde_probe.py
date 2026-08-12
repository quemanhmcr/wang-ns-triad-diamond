from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.continuum_helical_edge_measure_pde_probe import (
    _divergence_norm,
    _gradient_energy,
    _nonlinear_term,
    _rk4_step,
    _spectral_average_inner,
    _spectral_geometry,
    _trapezoid,
)
from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_hard_cell_single_charge_quotient import (
    fine_hard_role_map,
    hard_cell_single_charge_quotient,
    pushforward_restricted_hard_cell_donor_work,
    single_hard_role_map,
)
from src.cyclic_helical_triad_donor_kernel_pde_probe import _selected_closed_triad
from src.mixed_fate_reserved_young_handoff_pde_probe import adversarial_mixed_fate_initial_state

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_CYCLIC_HARD_CELL_SINGLE_CHARGE_AUDIT__"
    "ACTUAL_DW_MINUS_ROWS_TO_CANONICAL_DW_PLUS_FATE_COLUMNS__COARSE_SELF_LOOPS_ZERO_DEPTH__"
    "ONE_DONOR_TWO_RECIPIENT_AND_PHASE_REVERSED_TWO_DONOR_ONE_RECIPIENT"
)


@dataclass(frozen=True)
class HardCellSingleChargeNSRun:
    resolution: int
    cutoff: int
    steps: int
    snapshots: int
    duration: float
    viscosity: float
    amplitude: float
    phase_sign: int
    worst_balance_native_residual: float
    worst_donor_marginal_native_residual: float
    worst_recipient_marginal_native_residual: float
    worst_fate_partition_native_residual: float
    worst_restricted_pushforward_native_residual: float
    worst_coarse_total_native_residual: float
    worst_coarse_fate_native_residual: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float
    initial_donor_count: int
    initial_recipient_count: int
    initial_good_recipient_mass: float
    initial_bad_recipient_mass: float
    initial_overlapping_recipient_charge_count: int
    initial_coarse_self_loop_fraction: float
    selected_total_positive_work_snapshots: tuple[float, ...]
    selected_good_recipient_work_snapshots: tuple[float, ...]
    selected_bad_recipient_work_snapshots: tuple[float, ...]
    selected_native_work_mass_scale_snapshots: tuple[float, ...]


@dataclass(frozen=True)
class HardCellSingleChargePDEProbe:
    status: str
    runs: tuple[HardCellSingleChargeNSRun, ...]
    common_cutoff: int
    phase_sign: int
    maximum_total_positive_work_representation_native_residual: float
    maximum_good_recipient_work_representation_native_residual: float
    maximum_bad_recipient_work_representation_native_residual: float


def _run_one(
    *,
    resolution: int,
    cutoff: int,
    steps: int,
    viscosity: float,
    amplitude: float,
    duration: float,
    snapshot_count: int,
    phase_sign: int,
) -> HardCellSingleChargeNSRun:
    n = int(resolution)
    count = int(steps)
    snaps = int(snapshot_count)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    sign = int(phase_sign)
    if n < 24 or n % 2:
        raise ValueError("hard-cell single-charge NS probe requires an even FFT grid at least 24")
    if count < 16 or snaps < 3 or snaps > count + 1:
        raise ValueError("hard-cell single-charge NS probe requires at least sixteen RK4 steps and three snapshots")
    if not all(math.isfinite(v) and v > 0.0 for v in (nu, amp, horizon)):
        raise ValueError("positive finite NS audit parameters required")
    if sign not in (-1, 1):
        raise ValueError("phase_sign must be +1 or -1")

    k, k2, dealias, actual_cutoff = _spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("requested common Galerkin cutoff changed")
    state = adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=amp)
    state = sign * state
    dt = horizon / count
    sample_indices = tuple(sorted({round(j * count / (snaps - 1)) for j in range(snaps)}))
    grid_times = tuple(i * dt for i in range(count + 1))

    energy: list[float] = []
    gradient: list[float] = []
    nonlinear_work: list[float] = []
    divergence: list[float] = []
    total_snapshots: list[float] = []
    good_snapshots: list[float] = []
    bad_snapshots: list[float] = []
    native_snapshots: list[float] = []
    wb = wd = wr = wf = wrestrict = wcoarse_total = wcoarse_fate = 0.0
    initial_donors = initial_recipients = initial_overlap = 0
    initial_good = initial_bad = initial_self_fraction = 0.0
    discrete_triad_qmass = 1.0 / unitary_fourier_convolution_factor()

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        energy.append(_spectral_average_inner(state, state, n))
        gradient.append(_gradient_energy(state, k2, n))
        nonlinear_work.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        divergence.append(_divergence_norm(state, k, n))
        if step in sample_indices:
            triad = _selected_closed_triad(state)
            fine = hard_cell_single_charge_quotient(
                triad,
                quotient_measure_mass=discrete_triad_qmass,
                mode_roles=fine_hard_role_map(triad),
            )
            coarse = hard_cell_single_charge_quotient(
                triad,
                quotient_measure_mass=discrete_triad_qmass,
                mode_roles=single_hard_role_map(triad),
            )
            donor_cell = fine.donor_charges[0].cell
            restricted = pushforward_restricted_hard_cell_donor_work(
                fine, donor_cells=(donor_cell,)
            )
            native = fine.native_work_mass_scale
            wb = max(wb, fine.total_balance_native_residual)
            wd = max(wd, fine.worst_donor_marginal_native_residual)
            wr = max(wr, fine.worst_recipient_marginal_native_residual)
            wf = max(wf, fine.recipient_fate_partition_native_residual)
            wrestrict = max(wrestrict, restricted.mass_conservation_native_residual)
            wcoarse_total = max(
                wcoarse_total,
                abs(coarse.total_positive_work_mass - fine.total_positive_work_mass) / native,
            )
            wcoarse_fate = max(
                wcoarse_fate,
                abs(coarse.good_recipient_mass - fine.good_recipient_mass) / native,
                abs(coarse.bad_recipient_mass - fine.bad_recipient_mass) / native,
            )
            if coarse.self_loop_atom_count != len(coarse.atoms):
                raise AssertionError("maximal hard coarsening lost its same-time self-loop structure on evolved NS")
            if abs(coarse.self_loop_mass - coarse.total_positive_work_mass) > 5.0e-10 * native:
                raise AssertionError("coarse self-loop did not retain all same-time transported work")
            total_snapshots.append(fine.total_positive_work_mass)
            good_snapshots.append(fine.good_recipient_mass)
            bad_snapshots.append(fine.bad_recipient_mass)
            native_snapshots.append(native)
            if step == 0:
                initial_donors = triad.donor_kernel.donor_count
                initial_recipients = triad.donor_kernel.recipient_count
                initial_good = fine.good_recipient_mass
                initial_bad = fine.bad_recipient_mass
                initial_overlap = fine.overlapping_recipient_charge_count
                initial_self_fraction = coarse.self_loop_mass / coarse.total_positive_work_mass
                if sign == 1:
                    if not (initial_donors == 1 and initial_recipients == 2):
                        raise AssertionError("physical + phase NS initial triad lost one-donor/two-recipient pattern")
                    if not (initial_good > 0.0 and initial_bad > 0.0):
                        raise AssertionError("physical + phase NS initial triad lost good/bad recipient split")
                else:
                    if not (initial_donors == 2 and initial_recipients == 1):
                        raise AssertionError("phase-reversed physical NS initial triad lost two-donor/one-recipient pattern")
                    if initial_overlap != 1:
                        raise AssertionError("phase-reversed two-donor NS triad did not recombine to one recipient charge")
        if step < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    initial_energy = energy[0]
    balance = abs(energy[-1] - initial_energy + 2.0 * nu * _trapezoid(gradient, grid_times)) / initial_energy
    nonlinear_scale = initial_energy / horizon
    max_nonlinear = max(abs(v) for v in nonlinear_work) / nonlinear_scale
    max_divergence = max(divergence) / math.sqrt(initial_energy)
    if balance > 5.0e-5 or max_nonlinear > 5.0e-10 or max_divergence > 5.0e-11:
        raise AssertionError("hard-cell single-charge probe trajectory lost a native Navier-Stokes invariant")
    if max(wb, wd, wr, wf, wrestrict, wcoarse_total, wcoarse_fate) > 5.0e-8:
        raise AssertionError("hard-cell single-charge quotient left its native physical work scale on evolved NS")
    return HardCellSingleChargeNSRun(
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        snapshots=len(sample_indices),
        duration=horizon,
        viscosity=nu,
        amplitude=amp,
        phase_sign=sign,
        worst_balance_native_residual=wb,
        worst_donor_marginal_native_residual=wd,
        worst_recipient_marginal_native_residual=wr,
        worst_fate_partition_native_residual=wf,
        worst_restricted_pushforward_native_residual=wrestrict,
        worst_coarse_total_native_residual=wcoarse_total,
        worst_coarse_fate_native_residual=wcoarse_fate,
        global_energy_balance_relative_residual=balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
        maximum_divergence_relative_to_initial_l2=max_divergence,
        initial_donor_count=initial_donors,
        initial_recipient_count=initial_recipients,
        initial_good_recipient_mass=initial_good,
        initial_bad_recipient_mass=initial_bad,
        initial_overlapping_recipient_charge_count=initial_overlap,
        initial_coarse_self_loop_fraction=initial_self_fraction,
        selected_total_positive_work_snapshots=tuple(total_snapshots),
        selected_good_recipient_work_snapshots=tuple(good_snapshots),
        selected_bad_recipient_work_snapshots=tuple(bad_snapshots),
        selected_native_work_mass_scale_snapshots=tuple(native_snapshots),
    )


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    cutoff: int = 7,
    steps: int = 32,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    snapshot_count: int = 5,
    phase_sign: int = 1,
) -> HardCellSingleChargePDEProbe:
    resolved = tuple(int(n) for n in resolutions)
    if not resolved:
        raise ValueError("at least one FFT representation required")
    runs = tuple(
        _run_one(
            resolution=n,
            cutoff=int(cutoff),
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            duration=float(duration),
            snapshot_count=int(snapshot_count),
            phase_sign=int(phase_sign),
        )
        for n in resolved
    )
    if len({run.snapshots for run in runs}) != 1:
        raise AssertionError("representation runs used different physical snapshot counts")
    worst_total = worst_good = worst_bad = 0.0
    for j in range(runs[0].snapshots):
        native = max(run.selected_native_work_mass_scale_snapshots[j] for run in runs)
        if native <= 0.0:
            raise AssertionError("selected evolved triad lost positive native work scale")
        totals = tuple(run.selected_total_positive_work_snapshots[j] for run in runs)
        goods = tuple(run.selected_good_recipient_work_snapshots[j] for run in runs)
        bads = tuple(run.selected_bad_recipient_work_snapshots[j] for run in runs)
        worst_total = max(worst_total, (max(totals) - min(totals)) / native)
        worst_good = max(worst_good, (max(goods) - min(goods)) / native)
        worst_bad = max(worst_bad, (max(bads) - min(bads)) / native)
    if max(worst_total, worst_good, worst_bad) > 5.0e-8:
        raise AssertionError("same cutoff Galerkin NS hard-cell charges changed under FFT representation")
    return HardCellSingleChargePDEProbe(
        status=STATUS,
        runs=runs,
        common_cutoff=int(cutoff),
        phase_sign=int(phase_sign),
        maximum_total_positive_work_representation_native_residual=worst_total,
        maximum_good_recipient_work_representation_native_residual=worst_good,
        maximum_bad_recipient_work_representation_native_residual=worst_bad,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", type=int, nargs="+", default=(24, 28))
    parser.add_argument("--cutoff", type=int, default=7)
    parser.add_argument("--steps", type=int, default=32)
    parser.add_argument("--viscosity", type=float, default=0.03)
    parser.add_argument("--amplitude", type=float, default=1.0)
    parser.add_argument("--duration", type=float, default=0.001)
    parser.add_argument("--snapshots", type=int, default=5)
    parser.add_argument("--phase-sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--outdir", type=Path, default=Path("results-cyclic-hard-cell-single-charge-pde"))
    args = parser.parse_args()
    out = run_probe(
        resolutions=args.resolutions,
        cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        snapshot_count=args.snapshots,
        phase_sign=args.phase_sign,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "probe.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    lines = [
        "# Actual Galerkin NS cyclic hard-cell single-charge audit",
        "",
        f"Status: **{STATUS}**.",
        "",
        "The probe evolves the repository's real 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system. At each physical snapshot it reads the actual three-root signed helical triad work, forms the certified cyclic donor kernel, and only then pushes donor/recipient roots through deterministic hard-cell maps.",
        "",
        f"- common cutoff: `{out.common_cutoff}`",
        f"- phase sign: `{out.phase_sign:+d}`",
        f"- FFT representations: `{', '.join(str(run.resolution) for run in out.runs)}`",
        f"- maximum total-work representation native residual: `{out.maximum_total_positive_work_representation_native_residual:.3e}`",
        f"- maximum good-recipient representation native residual: `{out.maximum_good_recipient_work_representation_native_residual:.3e}`",
        f"- maximum bad-recipient representation native residual: `{out.maximum_bad_recipient_work_representation_native_residual:.3e}`",
    ]
    for run in out.runs:
        lines.extend([
            "",
            f"## resolution {run.resolution}",
            f"- steps/snapshots: `{run.steps}` / `{run.snapshots}`",
            f"- worst row/column native residuals: `{run.worst_donor_marginal_native_residual:.3e}` / `{run.worst_recipient_marginal_native_residual:.3e}`",
            f"- worst fate-partition native residual: `{run.worst_fate_partition_native_residual:.3e}`",
            f"- worst restricted-donor pushforward native residual: `{run.worst_restricted_pushforward_native_residual:.3e}`",
            f"- worst fine/coarse total/fate native residuals: `{run.worst_coarse_total_native_residual:.3e}` / `{run.worst_coarse_fate_native_residual:.3e}`",
            f"- NS energy-balance residual: `{run.global_energy_balance_relative_residual:.3e}`",
            f"- initial donor/recipient counts: `{run.initial_donor_count}` / `{run.initial_recipient_count}`",
            f"- initial good/bad recipient masses: `{run.initial_good_recipient_mass:.12g}` / `{run.initial_bad_recipient_mass:.12g}`",
            f"- initial overlapping recipient charges: `{run.initial_overlapping_recipient_charge_count}`",
            f"- initial maximal-coarsening self-loop fraction: `{run.initial_coarse_self_loop_fraction:.12g}`",
        ])
    lines.extend([
        "",
        "Changing the sign of the whole divergence-free initial field is used only to expose the opposite physical cubic-work sign pattern at t=0: the + branch begins one-donor/two-recipient, while the - branch begins two-donor/one-recipient. Both are then evolved by the same Navier--Stokes equations.",
        "",
        "Coarse self-loops remain real same-time redistribution. They create no event depth and no scale progress. Geometry-bad recipient work remains the existing transfer-loss causal sublaw, not vanished PDE energy. No between-time deposit/withdrawal matching and no global-regularity claim are made.",
        "",
    ])
    summary = "\n".join(lines)
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
