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
    _nonlinear_term,
    _rk4_step,
    _spectral_average_inner,
    _spectral_geometry,
    _trapezoid,
)
from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel
from src.cyclic_helical_triad_donor_kernel_pde_probe import _selected_closed_triad
from src.helical_mode_set_energy_continuity import flow_atoms_from_cyclic_kernel
from src.mixed_fate_reserved_young_handoff_pde_probe import adversarial_mixed_fate_initial_state
from src.radial_spectral_crossing_layer_cake import (
    finite_radial_log_action,
    mode_radius,
    radial_exterior_balance,
    truncated_radial_layer_cake,
)

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_RADIAL_SPECTRAL_CROSSING__"
    "ACTUAL_TAIL_STOCK_SIGNED_NONLINEAR_WORK_VISCOSITY__"
    "SELECTED_CYCLIC_UPWARD_DOWNWARD_CROSSING__SAME_CUTOFF_CROSS_FFT"
)


def _relative(actual: float, expected: float, scale: float) -> float:
    gap = abs(float(actual) - float(expected))
    native = abs(float(scale))
    if native == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / native


def _native_spread(values: Sequence[float], scale: float) -> float:
    if not values:
        return 0.0
    native = abs(float(scale))
    if not math.isfinite(native) or native <= 0.0:
        raise ValueError("positive finite representation error envelope required")
    return (max(values) - min(values)) / native


@dataclass(frozen=True)
class RadialSpectralCrossingNSRun:
    resolution: int
    cutoff: int
    steps: int
    duration: float
    viscosity: float
    amplitude: float
    phase_sign: int
    radial_boundary: float
    initial_tail_energy: float
    final_tail_energy: float
    integrated_signed_tail_nonlinear_work: float
    tail_viscous_dissipation: float
    tail_interval_balance_native_residual: float
    native_tail_energy_throughput_scale: float
    integrated_selected_upward_crossing: float
    integrated_selected_downward_crossing: float
    integrated_selected_total_flow: float
    integrated_selected_signed_crossing: float
    selected_crossing_time_action_scale: float
    initial_selected_upward_crossing: float
    initial_selected_downward_crossing: float
    worst_selected_radial_divergence_native_residual: float
    worst_selected_truncated_layer_cake_native_residual: float
    worst_selected_full_log_marginal_native_residual: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float


@dataclass(frozen=True)
class RadialSpectralCrossingPDEProbe:
    status: str
    runs: tuple[RadialSpectralCrossingNSRun, ...]
    common_cutoff: int
    radial_boundary: float
    phase_sign: int
    tail_representation_native_scale: float
    selected_crossing_representation_native_scale: float
    maximum_initial_tail_energy_representation_native_residual: float
    maximum_final_tail_energy_representation_native_residual: float
    maximum_integrated_tail_work_representation_native_residual: float
    maximum_tail_viscosity_representation_native_residual: float
    maximum_integrated_selected_upward_representation_native_residual: float
    maximum_integrated_selected_downward_representation_native_residual: float


def _run_one(
    *,
    resolution: int,
    cutoff: int,
    steps: int,
    viscosity: float,
    amplitude: float,
    duration: float,
    radial_boundary: float,
    phase_sign: int,
) -> RadialSpectralCrossingNSRun:
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    R = float(radial_boundary)
    sign = int(phase_sign)
    if n < 24 or n % 2:
        raise ValueError("radial spectral crossing NS audit requires an even FFT grid at least 24")
    if count < 16:
        raise ValueError("radial spectral crossing NS audit requires at least sixteen RK4 steps")
    if not all(math.isfinite(v) and v > 0.0 for v in (nu, amp, horizon, R)):
        raise ValueError("positive finite radial NS audit parameters required")
    if sign not in (-1, 1):
        raise ValueError("phase_sign must be +1 or -1")

    k, k2, dealias, actual_cutoff = _spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("common Galerkin cutoff changed")
    state = sign * adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=amp)
    dt = horizon / count
    times = tuple(j * dt for j in range(count + 1))
    radius_grid = np.sqrt(k2)
    tail = radius_grid > R
    tail_mask = tail[None, ...]
    discrete_qmass = 1.0 / unitary_fourier_convolution_factor()

    tail_energy: list[float] = []
    tail_gradient: list[float] = []
    tail_signed_work: list[float] = []
    selected_up: list[float] = []
    selected_down: list[float] = []
    selected_total: list[float] = []
    selected_signed: list[float] = []
    global_energy: list[float] = []
    global_gradient: list[float] = []
    global_nonlinear: list[float] = []
    divergence: list[float] = []
    worst_div = worst_layer = worst_full = 0.0
    initial_up = initial_down = 0.0

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        tail_state = state * tail_mask
        tail_nonlinear = nonlinear * tail_mask
        tail_energy.append(_spectral_average_inner(tail_state, tail_state, n))
        weighted_tail = np.sqrt(k2)[None, ...] * tail_state
        tail_gradient.append(_spectral_average_inner(weighted_tail, weighted_tail, n))
        tail_signed_work.append(-2.0 * _spectral_average_inner(tail_state, tail_nonlinear, n))

        triad = _selected_closed_triad(state)
        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=discrete_qmass)
        if not kernel.numerically_resolved_transport:
            raise AssertionError("selected evolved physical triad fell below donor-sign resolution")
        atoms = flow_atoms_from_cyclic_kernel(kernel)
        balance = radial_exterior_balance(atoms, radius=R)
        all_radii = tuple(
            mode_radius(mode)
            for atom in atoms
            for mode in (atom.donor_mode, atom.recipient_mode)
        )
        layer = truncated_radial_layer_cake(
            atoms,
            lower_radius=0.5 * min(all_radii),
            upper_radius=1.5 * max(all_radii),
        )
        full = finite_radial_log_action(atoms)
        selected_up.append(balance.upward_crossing_flow)
        selected_down.append(balance.downward_crossing_flow)
        selected_total.append(kernel.total_mass)
        selected_signed.append(balance.tail_signed_work)
        worst_div = max(worst_div, balance.tail_divergence_native_residual)
        worst_layer = max(
            worst_layer,
            layer.upward_atomwise_identity_native_residual,
            layer.downward_atomwise_identity_native_residual,
            layer.signed_marginal_identity_native_residual,
        )
        worst_full = max(worst_full, full.signed_marginal_identity_native_residual)
        direct_selected_tail = math.fsum(
            slot.signed_work
            for slot in triad.slots
            if mode_radius(slot.edge_identity.child) > R
        )
        direct_native = max(kernel.native_work_mass_scale, 1.0e-300)
        if _relative(balance.tail_signed_work, direct_selected_tail, direct_native) > 5.0e-9:
            raise AssertionError("selected actual NS triad radial crossing did not reconstruct its high-mode signed work")
        if step == 0:
            initial_up = balance.upward_crossing_flow
            initial_down = balance.downward_crossing_flow
        global_energy.append(_spectral_average_inner(state, state, n))
        weighted = np.sqrt(k2)[None, ...] * state
        global_gradient.append(_spectral_average_inner(weighted, weighted, n))
        global_nonlinear.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        divergence.append(_divergence_norm(state, k, n))
        if step < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    int_tail_signed = _trapezoid(tail_signed_work, times)
    tail_visc = 2.0 * nu * _trapezoid(tail_gradient, times)
    tail_lhs = tail_energy[-1] + tail_visc
    tail_rhs = tail_energy[0] + int_tail_signed
    tail_native = max(
        tail_energy[0] + abs(int_tail_signed),
        tail_energy[-1] + tail_visc,
        abs(int_tail_signed) + tail_visc,
        1.0e-300,
    )
    tail_res = _relative(tail_lhs, tail_rhs, tail_native)

    int_up = _trapezoid(selected_up, times)
    int_down = _trapezoid(selected_down, times)
    int_total = _trapezoid(selected_total, times)
    int_selected_signed = _trapezoid(selected_signed, times)
    crossing_native = max(int_total, int_up + int_down, 1.0e-300)
    if _relative(int_selected_signed, int_up - int_down, crossing_native) > 5.0e-9:
        raise AssertionError("integrated selected radial crossing lost upward-minus-downward identity")

    e0 = global_energy[0]
    global_balance = abs(global_energy[-1] - e0 + 2.0 * nu * _trapezoid(global_gradient, times)) / e0
    nonlinear_scale = e0 / horizon
    max_nonlinear = max(abs(v) for v in global_nonlinear) / nonlinear_scale
    max_div = max(divergence) / math.sqrt(e0)
    if tail_res > 5.0e-5:
        raise AssertionError("actual radial tail stock/work/viscosity balance left its finite-step native scale")
    if global_balance > 5.0e-5 or max_nonlinear > 5.0e-10 or max_div > 5.0e-11:
        raise AssertionError("radial crossing probe trajectory lost a native global Navier-Stokes invariant")
    if max(worst_div, worst_layer, worst_full) > 5.0e-9:
        raise AssertionError("selected evolved cyclic radial law left its native physical scale")

    direction_scale = max(initial_up + initial_down, 1.0e-300)
    if sign > 0:
        if not (initial_up > 0.0 and initial_down <= 5.0e-10 * direction_scale):
            raise AssertionError("positive-phase selected physical triad lost low-to-high radial crossing")
    else:
        if not (initial_down > 0.0 and initial_up <= 5.0e-10 * direction_scale):
            raise AssertionError("sign-reversed selected physical triad lost high-to-low radial crossing")

    return RadialSpectralCrossingNSRun(
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        duration=horizon,
        viscosity=nu,
        amplitude=amp,
        phase_sign=sign,
        radial_boundary=R,
        initial_tail_energy=tail_energy[0],
        final_tail_energy=tail_energy[-1],
        integrated_signed_tail_nonlinear_work=int_tail_signed,
        tail_viscous_dissipation=tail_visc,
        tail_interval_balance_native_residual=tail_res,
        native_tail_energy_throughput_scale=tail_native,
        integrated_selected_upward_crossing=int_up,
        integrated_selected_downward_crossing=int_down,
        integrated_selected_total_flow=int_total,
        integrated_selected_signed_crossing=int_selected_signed,
        selected_crossing_time_action_scale=crossing_native,
        initial_selected_upward_crossing=initial_up,
        initial_selected_downward_crossing=initial_down,
        worst_selected_radial_divergence_native_residual=worst_div,
        worst_selected_truncated_layer_cake_native_residual=worst_layer,
        worst_selected_full_log_marginal_native_residual=worst_full,
        global_energy_balance_relative_residual=global_balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
        maximum_divergence_relative_to_initial_l2=max_div,
    )


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    cutoff: int = 7,
    steps: int = 48,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    radial_boundary: float = 8.0,
    phase_sign: int = 1,
) -> RadialSpectralCrossingPDEProbe:
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
            radial_boundary=float(radial_boundary),
            phase_sign=int(phase_sign),
        )
        for n in resolved
    )
    tail_native = max(r.native_tail_energy_throughput_scale for r in runs)
    crossing_native = max(r.selected_crossing_time_action_scale for r in runs)
    metrics = (
        _native_spread([r.initial_tail_energy for r in runs], tail_native),
        _native_spread([r.final_tail_energy for r in runs], tail_native),
        _native_spread([r.integrated_signed_tail_nonlinear_work for r in runs], tail_native),
        _native_spread([r.tail_viscous_dissipation for r in runs], tail_native),
        _native_spread([r.integrated_selected_upward_crossing for r in runs], crossing_native),
        _native_spread([r.integrated_selected_downward_crossing for r in runs], crossing_native),
    )
    if max(metrics) > 5.0e-7:
        raise AssertionError("same finite cutoff radial NS system changed under FFT representation on its native physical envelope")
    return RadialSpectralCrossingPDEProbe(
        status=STATUS,
        runs=runs,
        common_cutoff=int(cutoff),
        radial_boundary=float(radial_boundary),
        phase_sign=int(phase_sign),
        tail_representation_native_scale=tail_native,
        selected_crossing_representation_native_scale=crossing_native,
        maximum_initial_tail_energy_representation_native_residual=metrics[0],
        maximum_final_tail_energy_representation_native_residual=metrics[1],
        maximum_integrated_tail_work_representation_native_residual=metrics[2],
        maximum_tail_viscosity_representation_native_residual=metrics[3],
        maximum_integrated_selected_upward_representation_native_residual=metrics[4],
        maximum_integrated_selected_downward_representation_native_residual=metrics[5],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=STATUS)
    parser.add_argument("--resolutions", type=int, nargs="+", default=(24, 28))
    parser.add_argument("--cutoff", type=int, default=7)
    parser.add_argument("--steps", type=int, default=48)
    parser.add_argument("--viscosity", type=float, default=0.03)
    parser.add_argument("--amplitude", type=float, default=1.0)
    parser.add_argument("--duration", type=float, default=0.001)
    parser.add_argument("--radius", type=float, default=8.0)
    parser.add_argument("--phase-sign", type=int, default=1)
    parser.add_argument("--outdir", type=Path, default=Path("results-radial-spectral-crossing-ns"))
    args = parser.parse_args()
    out = run_probe(
        resolutions=args.resolutions,
        cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        radial_boundary=args.radius,
        phase_sign=args.phase_sign,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "probe.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    lines = [
        "# Actual Galerkin NS radial spectral crossing audit",
        "",
        f"Status: **{STATUS}**.",
        "",
        "The probe evolves the actual 2/3-dealiased incompressible Fourier--Galerkin Navier--Stokes system.  It reads full radial-tail stock, signed nonlinear tail work and viscosity at every RK4 output time.  Independently, it reads one actual evolving closed helical triad from the same state and restricts its already-certified cyclic donor flow by the same Fourier radius; the selected triad is not substituted for the full tail law.",
        "",
        f"- common cutoff: `{out.common_cutoff}`",
        f"- radial boundary: `{out.radial_boundary:.12g}`",
        f"- phase sign: `{out.phase_sign}`",
        f"- FFT representations: `{', '.join(str(r.resolution) for r in out.runs)}`",
        f"- tail representation native scale: `{out.tail_representation_native_scale:.12g}`",
        f"- selected crossing representation native scale: `{out.selected_crossing_representation_native_scale:.12g}`",
        f"- max tail energy/work/viscosity representation residual: `{max(out.maximum_initial_tail_energy_representation_native_residual, out.maximum_final_tail_energy_representation_native_residual, out.maximum_integrated_tail_work_representation_native_residual, out.maximum_tail_viscosity_representation_native_residual):.3e}`",
        f"- max selected up/down representation residual: `{max(out.maximum_integrated_selected_upward_representation_native_residual, out.maximum_integrated_selected_downward_representation_native_residual):.3e}`",
    ]
    for run in out.runs:
        lines.extend([
            "",
            f"## resolution {run.resolution}",
            f"- steps: `{run.steps}`",
            f"- initial/final tail energy: `{run.initial_tail_energy:.12g}` / `{run.final_tail_energy:.12g}`",
            f"- integrated signed tail nonlinear work: `{run.integrated_signed_tail_nonlinear_work:.12g}`",
            f"- tail viscous dissipation: `{run.tail_viscous_dissipation:.12g}`",
            f"- tail interval balance native residual: `{run.tail_interval_balance_native_residual:.3e}`",
            f"- integrated selected upward/downward crossing: `{run.integrated_selected_upward_crossing:.12g}` / `{run.integrated_selected_downward_crossing:.12g}`",
            f"- initial selected upward/downward crossing: `{run.initial_selected_upward_crossing:.12g}` / `{run.initial_selected_downward_crossing:.12g}`",
            f"- worst selected radial divergence residual: `{run.worst_selected_radial_divergence_native_residual:.3e}`",
            f"- worst selected layer-cake residual: `{run.worst_selected_truncated_layer_cake_native_residual:.3e}`",
            f"- global NS energy-balance residual: `{run.global_energy_balance_relative_residual:.3e}`",
        ])
    lines.extend([
        "",
        "The full-tail reading and selected-triad reading are separate physical observables on the same evolved state.  No net-tail Hahn law is minted, no selected triad is claimed to exhaust the full tail, and no radial crossing is promoted to a finite traffic budget or recursive event count.",
    ])
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
