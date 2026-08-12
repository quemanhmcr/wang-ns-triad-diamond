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
    _series_coefficient,
    _spectral_average_inner,
    _spectral_geometry,
    _trapezoid,
)
from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    register_closed_helical_triad,
    signed_good_side_recipient_certificate,
)
from src.helical import helical_basis
from src.mixed_fate_reserved_young_handoff_pde_probe import (
    CHILD,
    GOOD_PARENTS,
    HELICITIES,
    adversarial_mixed_fate_initial_state,
)

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_CYCLIC_HELICAL_TRIAD_DONOR_KERNEL_AUDIT__"
    "ACTUAL_MODAL_COEFFICIENTS__THREE_ROOT_ENERGY_CONSERVATION__DW_MINUS_TO_DW_PLUS_MARGINALS__"
    "SIGNED_GOOD_NONFORWARD_SIDE_RECIPIENT"
)

CLOSED_MODES = (
    tuple(-v for v in CHILD),
    GOOD_PARENTS[0],
    GOOD_PARENTS[1],
)
CLOSED_HELICITIES = (HELICITIES[2], HELICITIES[0], HELICITIES[1])


def _helical_amplitude(
    state_hat: np.ndarray, wavevector: tuple[int, int, int], helicity: int
) -> complex:
    coeff = _series_coefficient(state_hat, wavevector)
    h = helical_basis(np.asarray(wavevector, dtype=float), int(helicity))
    return complex(np.vdot(h, coeff))


def _selected_closed_triad(state_hat: np.ndarray):
    amps = tuple(
        _helical_amplitude(state_hat, k, s)
        for k, s in zip(CLOSED_MODES, CLOSED_HELICITIES)
    )
    return register_closed_helical_triad(
        wavevectors=tuple(np.asarray(k, dtype=float) for k in CLOSED_MODES),
        helicities=CLOSED_HELICITIES,
        amplitudes=amps,
    )


@dataclass(frozen=True)
class CyclicDonorKernelNSRun:
    resolution: int
    cutoff: int
    steps: int
    snapshots: int
    duration: float
    viscosity: float
    amplitude: float
    worst_cyclic_energy_conservation_relative: float
    worst_cyclic_coupling_native_residual: float
    worst_measure_donor_marginal_relative: float
    worst_measure_recipient_marginal_relative: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float
    initial_signed_good_efficiency: float
    initial_side_to_child_ratio: float
    initial_child_to_donor_ratio: float
    initial_side_to_donor_ratio: float
    initial_side_forward_ratio: float
    initial_side_geometric_multiplier: float
    initial_donor_count: int
    initial_recipient_count: int
    selected_triad_positive_work_snapshots: tuple[float, ...]
    selected_root_work_snapshots: tuple[tuple[float, float, float], ...]


@dataclass(frozen=True)
class CyclicDonorKernelPDEProbe:
    status: str
    runs: tuple[CyclicDonorKernelNSRun, ...]
    common_cutoff: int
    maximum_selected_triad_positive_work_representation_spread: float
    maximum_selected_root_work_representation_spread: float


def _run_one(
    *,
    resolution: int,
    cutoff: int,
    steps: int,
    viscosity: float,
    amplitude: float,
    duration: float,
    snapshot_count: int,
) -> CyclicDonorKernelNSRun:
    n = int(resolution)
    count = int(steps)
    snaps = int(snapshot_count)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    if n < 24 or n % 2:
        raise ValueError("cyclic donor NS probe requires an even FFT grid at least 24")
    if count < 16 or snaps < 3 or snaps > count + 1:
        raise ValueError("cyclic donor NS probe requires at least sixteen RK4 steps and three snapshots")
    if not all(math.isfinite(v) and v > 0.0 for v in (nu, amp, horizon)):
        raise ValueError("positive finite NS audit parameters required")
    k, k2, dealias, actual_cutoff = _spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("requested common Galerkin cutoff changed")
    state = adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=amp)
    dt = horizon / count
    sample_indices = tuple(sorted({round(j * count / (snaps - 1)) for j in range(snaps)}))
    grid_times = tuple(i * dt for i in range(count + 1))

    energy: list[float] = []
    gradient: list[float] = []
    nonlinear_work: list[float] = []
    divergence: list[float] = []
    total_positive_snapshots: list[float] = []
    root_work_snapshots: list[tuple[float, float, float]] = []
    worst_energy = worst_coupling = worst_donor = worst_recipient = 0.0
    initial_efficiency = initial_side = initial_child_donor = initial_side_donor = 0.0
    initial_side_forward = initial_side_j = 0.0
    initial_donors = initial_recipients = 0
    discrete_triad_qmass = 1.0 / unitary_fourier_convolution_factor()

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        energy.append(_spectral_average_inner(state, state, n))
        gradient.append(_gradient_energy(state, k2, n))
        nonlinear_work.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        divergence.append(_divergence_norm(state, k, n))
        if step in sample_indices:
            triad = _selected_closed_triad(state)
            measure = cyclic_triad_measure_kernel(
                triad, quotient_measure_mass=discrete_triad_qmass
            )
            worst_energy = max(worst_energy, triad.signed_energy_conservation_residual)
            worst_coupling = max(worst_coupling, triad.cyclic_coupling_native_residual)
            worst_donor = max(worst_donor, measure.donor_marginal_residual)
            worst_recipient = max(worst_recipient, measure.recipient_marginal_residual)
            total_positive_snapshots.append(measure.total_mass)
            root_work_snapshots.append(tuple(slot.signed_work for slot in triad.slots))
            if step == 0:
                recipient = triad.slot_for_edge_child_wavevector(CHILD)
                side = signed_good_side_recipient_certificate(
                    triad, recipient_closed_mode_index=recipient.closed_mode_index
                )
                initial_efficiency = side.recipient_signed_efficiency
                initial_side = side.side_to_recipient_ratio
                initial_child_donor = side.recipient_to_donor_ratio
                initial_side_donor = side.side_to_donor_ratio
                initial_side_forward = side.side_forward_ratio
                initial_side_j = side.side_geometric_multiplier
                initial_donors = triad.donor_kernel.donor_count
                initial_recipients = triad.donor_kernel.recipient_count
                if not (initial_donors == 1 and initial_recipients == 2):
                    raise AssertionError("engineered signed-good NS triad lost its one-donor/two-recipient pattern")
        if step < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    initial_energy = energy[0]
    balance = abs(energy[-1] - initial_energy + 2.0 * nu * _trapezoid(gradient, grid_times)) / initial_energy
    nonlinear_scale = initial_energy / horizon
    max_nonlinear = max(abs(v) for v in nonlinear_work) / nonlinear_scale
    max_divergence = max(divergence) / math.sqrt(initial_energy)
    if balance > 5.0e-5 or max_nonlinear > 5.0e-10 or max_divergence > 5.0e-11:
        raise AssertionError("cyclic donor probe trajectory lost a native Navier-Stokes invariant")
    if max(worst_energy, worst_coupling, worst_donor, worst_recipient) > 4.0e-8:
        raise AssertionError("cyclic donor kernel lost physical work provenance on evolved NS")
    return CyclicDonorKernelNSRun(
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        snapshots=len(sample_indices),
        duration=horizon,
        viscosity=nu,
        amplitude=amp,
        worst_cyclic_energy_conservation_relative=worst_energy,
        worst_cyclic_coupling_native_residual=worst_coupling,
        worst_measure_donor_marginal_relative=worst_donor,
        worst_measure_recipient_marginal_relative=worst_recipient,
        global_energy_balance_relative_residual=balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
        maximum_divergence_relative_to_initial_l2=max_divergence,
        initial_signed_good_efficiency=initial_efficiency,
        initial_side_to_child_ratio=initial_side,
        initial_child_to_donor_ratio=initial_child_donor,
        initial_side_to_donor_ratio=initial_side_donor,
        initial_side_forward_ratio=initial_side_forward,
        initial_side_geometric_multiplier=initial_side_j,
        initial_donor_count=initial_donors,
        initial_recipient_count=initial_recipients,
        selected_triad_positive_work_snapshots=tuple(total_positive_snapshots),
        selected_root_work_snapshots=tuple(root_work_snapshots),
    )


def _relative_spread(values: Sequence[float]) -> float:
    if not values:
        return 0.0
    scale = max(abs(float(v)) for v in values)
    return 0.0 if scale == 0.0 else (max(values) - min(values)) / scale


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    cutoff: int = 7,
    steps: int = 32,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    snapshot_count: int = 5,
) -> CyclicDonorKernelPDEProbe:
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
        )
        for n in resolved
    )
    if len({run.snapshots for run in runs}) != 1:
        raise AssertionError("representation runs used different physical snapshot counts")
    worst_total = 0.0
    worst_root = 0.0
    for j in range(runs[0].snapshots):
        worst_total = max(
            worst_total,
            _relative_spread(tuple(run.selected_triad_positive_work_snapshots[j] for run in runs)),
        )
        for root in range(3):
            values = tuple(run.selected_root_work_snapshots[j][root] for run in runs)
            scale = max(abs(v) for v in values)
            spread = 0.0 if scale == 0.0 else (max(values) - min(values)) / scale
            worst_root = max(worst_root, abs(spread))
    if worst_total > 5.0e-8 or worst_root > 5.0e-8:
        raise AssertionError("the same cutoff-7 Galerkin NS system changed under FFT representation")
    return CyclicDonorKernelPDEProbe(
        status=STATUS,
        runs=runs,
        common_cutoff=int(cutoff),
        maximum_selected_triad_positive_work_representation_spread=worst_total,
        maximum_selected_root_work_representation_spread=worst_root,
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
    parser.add_argument("--outdir", type=Path, default=Path("results-cyclic-helical-triad-donor-kernel-pde"))
    args = parser.parse_args()
    out = run_probe(
        resolutions=args.resolutions,
        cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        snapshot_count=args.snapshots,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "probe.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    lines = [
        "# Actual Galerkin NS cyclic helical-triad donor-kernel audit",
        "",
        f"Status: **{STATUS}**.",
        "",
        "The probe evolves the repository's real 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system.  At each physical snapshot it reads the actual evolving helical coefficients of one closed triad and registers all three cyclic child-work roots before Hahn.",
        "",
        f"- common cutoff: `{out.common_cutoff}`",
        f"- FFT representations: `{', '.join(str(run.resolution) for run in out.runs)}`",
        f"- maximum selected-triad positive-work representation spread: `{out.maximum_selected_triad_positive_work_representation_spread:.3e}`",
        f"- maximum selected root-work representation spread: `{out.maximum_selected_root_work_representation_spread:.3e}`",
    ]
    for run in out.runs:
        lines.extend([
            "",
            f"## resolution {run.resolution}",
            f"- steps/snapshots: `{run.steps}` / `{run.snapshots}`",
            f"- worst cyclic energy-conservation residual: `{run.worst_cyclic_energy_conservation_relative:.3e}`",
            f"- worst donor/recipient measure marginal residuals: `{run.worst_measure_donor_marginal_relative:.3e}` / `{run.worst_measure_recipient_marginal_relative:.3e}`",
            f"- NS energy-balance residual: `{run.global_energy_balance_relative_residual:.3e}`",
            f"- initial signed-good efficiency: `{run.initial_signed_good_efficiency:.12g}`",
            f"- initial side/child work ratio: `{run.initial_side_to_child_ratio:.12g}`",
            f"- initial child/donor work ratio: `{run.initial_child_to_donor_ratio:.12g}`",
            f"- initial side/donor work ratio: `{run.initial_side_to_donor_ratio:.12g}`",
            f"- initial side forward ratio/J: `{run.initial_side_forward_ratio:.12g}` / `{run.initial_side_geometric_multiplier:.12g}`",
            f"- initial donor/recipient counts: `{run.initial_donor_count}` / `{run.initial_recipient_count}`",
        ])
    lines.extend([
        "",
        "The side recipient is actual positive nonforward NS work at the same triad/time.  The probe does not reinterpret it as dissipation or a reset.  The donor kernel adds provenance to canonical dW+ and creates no new event.  No global-regularity claim is made.",
        "",
    ])
    summary = "\n".join(lines)
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
