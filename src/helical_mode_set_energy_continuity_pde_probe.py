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
    _index,
    _nonlinear_term,
    _rk4_step,
    _series_coefficient,
    _snapshot_with_ledger,
    _spectral_average_inner,
    _spectral_geometry,
    _trapezoid,
)
from src.continuum_helical_edge_measure_registration import (
    HelicalModeIdentity,
    helical_coefficients,
)
from src.helical_mode_set_energy_continuity import interval_continuity_certificate
from src.mixed_fate_reserved_young_handoff_pde_probe import (
    CHILD,
    HELICITIES,
    adversarial_mixed_fate_initial_state,
)

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_HELICAL_MODE_ENERGY_CONTINUITY__"
    "GROSS_CANONICAL_EDGE_DW_PLUS_DW_MINUS__MODE_STOCK_PLUS_VISCOSITY__"
    "SAME_CUTOFF_CROSS_FFT_REPRESENTATION__NO_REHAHN_OR_TEMPORAL_MATCHING"
)


def _relative(actual: float, expected: float, scale: float) -> float:
    gap = abs(float(actual) - float(expected))
    s = abs(float(scale))
    if s == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / s


@dataclass(frozen=True)
class HelicalModeSnapshot:
    mode: HelicalModeIdentity
    energy: float
    positive_edge_work: float
    negative_edge_work: float
    signed_edge_work: float
    actual_modal_work: float
    signed_reconstruction_native_residual: float
    gross_work_scale: float


def _helical_mode_snapshot(
    state_hat: np.ndarray,
    nonlinear_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
    cutoff: int,
    *,
    child_mode: tuple[int, int, int],
    helicity: int,
) -> HelicalModeSnapshot:
    n = int(state_hat.shape[1])
    child = tuple(int(v) for v in child_mode)
    s = int(helicity)
    if s not in (-1, 1):
        raise ValueError("helicity must be plus or minus one")
    row, ledger = _snapshot_with_ledger(
        state_hat,
        k,
        k2,
        dealias,
        cutoff,
        child_mode=child,
        nonlinear_hat=nonlinear_hat,
    )
    uz = _series_coefficient(state_hat, child)
    source = -np.asarray(
        nonlinear_hat[(slice(None),) + _index(child, n)], dtype=complex
    ) / float(n**3)
    z = np.asarray(child, dtype=float)
    a = helical_coefficients(z, uz)[s]
    f = helical_coefficients(z, source)[s]
    energy = float(abs(a) ** 2)
    actual = 2.0 * float(np.real(np.conjugate(a) * f))

    signed_atoms: list[float] = []
    for fiber in ledger.physical_fibers:
        for atom in fiber.modal_atoms:
            identity = atom.physical_edge_identity
            if identity.child.wavevector == tuple(float(v) for v in child) and identity.child.helicity == s:
                signed_atoms.append(float(atom.signed_work_mass))
    if not signed_atoms:
        raise AssertionError("actual PDE helical mode lost its canonical child-edge atoms")
    positive = math.fsum(max(0.0, value) for value in signed_atoms)
    negative = math.fsum(max(0.0, -value) for value in signed_atoms)
    signed = positive - negative
    gross = max(positive + negative, abs(actual), abs(signed), 1.0e-300)
    residual = _relative(signed, actual, gross)
    if residual > 5.0e-9:
        raise AssertionError("gross canonical edge Hahn law did not reconstruct actual helical modal work")
    # The wavevector-level row is independently checked upstream.  Ensure the two
    # helicity rows remain a true refinement, not a second law.
    if positive > float(row["positive_edge_work"]) + 5.0e-9 * max(gross, float(row["positive_edge_work"])):
        raise AssertionError("helical positive edge restriction exceeded the canonical wavevector child law")
    if negative > float(row["negative_edge_work"]) + 5.0e-9 * max(gross, float(row["negative_edge_work"])):
        raise AssertionError("helical negative edge restriction exceeded the canonical wavevector child law")
    return HelicalModeSnapshot(
        mode=HelicalModeIdentity(tuple(float(v) for v in child), s),
        energy=energy,
        positive_edge_work=positive,
        negative_edge_work=negative,
        signed_edge_work=signed,
        actual_modal_work=actual,
        signed_reconstruction_native_residual=residual,
        gross_work_scale=gross,
    )


@dataclass(frozen=True)
class HelicalModeContinuityNSRun:
    resolution: int
    cutoff: int
    steps: int
    duration: float
    viscosity: float
    amplitude: float
    phase_sign: int
    child_mode: tuple[int, int, int]
    helicity: int
    initial_energy: float
    final_energy: float
    integrated_positive_work: float
    integrated_negative_work: float
    viscous_dissipation: float
    interval_continuity_native_residual: float
    worst_instantaneous_signed_reconstruction_native_residual: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float
    positive_work_steps: int
    negative_work_steps: int


@dataclass(frozen=True)
class HelicalModeContinuityPDEProbe:
    status: str
    runs: tuple[HelicalModeContinuityNSRun, ...]
    common_cutoff: int
    child_mode: tuple[int, int, int]
    helicity: int
    phase_sign: int
    maximum_initial_energy_representation_relative_residual: float
    maximum_final_energy_representation_relative_residual: float
    maximum_integrated_positive_work_representation_relative_residual: float
    maximum_integrated_negative_work_representation_relative_residual: float
    maximum_viscous_dissipation_representation_relative_residual: float


def _run_one(
    *,
    resolution: int,
    cutoff: int,
    steps: int,
    viscosity: float,
    amplitude: float,
    duration: float,
    child_mode: tuple[int, int, int],
    helicity: int,
    phase_sign: int,
) -> HelicalModeContinuityNSRun:
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    sign = int(phase_sign)
    if n < 24 or n % 2:
        raise ValueError("helical mode continuity NS audit requires an even FFT grid at least 24")
    if count < 16:
        raise ValueError("helical mode continuity NS audit requires at least sixteen RK4 steps")
    if not all(math.isfinite(v) and v > 0.0 for v in (nu, amp, horizon)):
        raise ValueError("positive finite NS audit parameters required")
    if sign not in (-1, 1):
        raise ValueError("phase_sign must be +1 or -1")

    k, k2, dealias, actual_cutoff = _spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("common Galerkin cutoff changed")
    state = adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=amp)
    state = sign * state
    dt = horizon / count
    times = tuple(j * dt for j in range(count + 1))
    energies: list[float] = []
    positives: list[float] = []
    negatives: list[float] = []
    signed_residuals: list[float] = []
    global_energy: list[float] = []
    gradient: list[float] = []
    nonlinear_work: list[float] = []
    divergence: list[float] = []

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        snap = _helical_mode_snapshot(
            state,
            nonlinear,
            k,
            k2,
            dealias,
            int(cutoff),
            child_mode=child_mode,
            helicity=int(helicity),
        )
        energies.append(snap.energy)
        positives.append(snap.positive_edge_work)
        negatives.append(snap.negative_edge_work)
        signed_residuals.append(snap.signed_reconstruction_native_residual)
        global_energy.append(_spectral_average_inner(state, state, n))
        gradient.append(_gradient_energy(state, k2, n))
        nonlinear_work.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        divergence.append(_divergence_norm(state, k, n))
        if step < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    p_int = _trapezoid(positives, times)
    n_int = _trapezoid(negatives, times)
    e_int = _trapezoid(energies, times)
    wave2 = math.fsum(float(v) * float(v) for v in child_mode)
    visc = 2.0 * nu * wave2 * e_int
    throughput = max(energies[0] + p_int, energies[-1] + n_int + visc, 1.0e-300)
    cert = interval_continuity_certificate(
        modes=(HelicalModeIdentity(tuple(float(v) for v in child_mode), int(helicity)),),
        initial_energy=energies[0],
        final_energy=energies[-1],
        integrated_inward_flow=p_int,
        integrated_outward_flow=n_int,
        viscous_dissipation=visc,
        native_energy_throughput_scale=throughput,
    )
    # Here singleton inward/outward are the gross canonical positive/negative edge
    # marginals.  The flow theorem separately proves these are boundary fluxes.
    if cert.balance_native_residual > 5.0e-5:
        raise AssertionError("actual helical modal stock/work/viscosity continuity failed")

    initial_global = global_energy[0]
    global_balance = abs(
        global_energy[-1] - initial_global + 2.0 * nu * _trapezoid(gradient, times)
    ) / initial_global
    nonlinear_scale = initial_global / horizon
    max_nonlin = max(abs(v) for v in nonlinear_work) / nonlinear_scale
    max_div = max(divergence) / math.sqrt(initial_global)
    if global_balance > 5.0e-5 or max_nonlin > 5.0e-10 or max_div > 5.0e-11:
        raise AssertionError("actual NS trajectory lost a native global invariant")

    return HelicalModeContinuityNSRun(
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        duration=horizon,
        viscosity=nu,
        amplitude=amp,
        phase_sign=sign,
        child_mode=tuple(int(v) for v in child_mode),
        helicity=int(helicity),
        initial_energy=energies[0],
        final_energy=energies[-1],
        integrated_positive_work=p_int,
        integrated_negative_work=n_int,
        viscous_dissipation=visc,
        interval_continuity_native_residual=cert.balance_native_residual,
        worst_instantaneous_signed_reconstruction_native_residual=max(signed_residuals),
        global_energy_balance_relative_residual=global_balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlin,
        maximum_divergence_relative_to_initial_l2=max_div,
        positive_work_steps=sum(v > 0.0 for v in positives),
        negative_work_steps=sum(v > 0.0 for v in negatives),
    )


def _spread(values: Sequence[float]) -> float:
    scale = max(max(abs(v) for v in values), 1.0e-300)
    return (max(values) - min(values)) / scale


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    cutoff: int = 7,
    steps: int = 32,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    child_mode: tuple[int, int, int] = CHILD,
    helicity: int = HELICITIES[2],
    phase_sign: int = 1,
) -> HelicalModeContinuityPDEProbe:
    grids = tuple(int(n) for n in resolutions)
    if not grids:
        raise ValueError("at least one FFT representation required")
    runs = tuple(
        _run_one(
            resolution=n,
            cutoff=int(cutoff),
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            duration=float(duration),
            child_mode=tuple(int(v) for v in child_mode),
            helicity=int(helicity),
            phase_sign=int(phase_sign),
        )
        for n in grids
    )
    metrics = (
        _spread([r.initial_energy for r in runs]),
        _spread([r.final_energy for r in runs]),
        _spread([r.integrated_positive_work for r in runs]),
        _spread([r.integrated_negative_work for r in runs]),
        _spread([r.viscous_dissipation for r in runs]),
    )
    if max(metrics) > 5.0e-7:
        raise AssertionError("same cutoff physical helical-mode continuity changed under FFT representation")
    return HelicalModeContinuityPDEProbe(
        status=STATUS,
        runs=runs,
        common_cutoff=int(cutoff),
        child_mode=tuple(int(v) for v in child_mode),
        helicity=int(helicity),
        phase_sign=int(phase_sign),
        maximum_initial_energy_representation_relative_residual=metrics[0],
        maximum_final_energy_representation_relative_residual=metrics[1],
        maximum_integrated_positive_work_representation_relative_residual=metrics[2],
        maximum_integrated_negative_work_representation_relative_residual=metrics[3],
        maximum_viscous_dissipation_representation_relative_residual=metrics[4],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", type=int, nargs="+", default=(24, 28))
    parser.add_argument("--cutoff", type=int, default=7)
    parser.add_argument("--steps", type=int, default=32)
    parser.add_argument("--viscosity", type=float, default=0.03)
    parser.add_argument("--amplitude", type=float, default=1.0)
    parser.add_argument("--duration", type=float, default=0.001)
    parser.add_argument("--helicity", type=int, choices=(-1, 1), default=HELICITIES[2])
    parser.add_argument("--phase-sign", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--outdir", type=Path, default=Path("results-helical-mode-set-energy-continuity-pde"))
    args = parser.parse_args()
    out = run_probe(
        resolutions=args.resolutions,
        cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        helicity=args.helicity,
        phase_sign=args.phase_sign,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "probe.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    lines = [
        "# Actual Galerkin NS helical-mode energy continuity audit",
        "",
        f"Status: **{STATUS}**.",
        "",
        "At every RK4 output time the probe reconstructs the full unordered canonical child-edge ledger from the actual evolving Fourier coefficients, restricts it to one physical child helicity, and integrates gross dW+ and dW- against the actual helical modal energy and viscous dissipation.",
        "",
        f"- common cutoff: `{out.common_cutoff}`",
        f"- child/helicity: `{out.child_mode}` / `{out.helicity:+d}`",
        f"- phase sign: `{out.phase_sign:+d}`",
        f"- FFT representations: `{', '.join(str(r.resolution) for r in out.runs)}`",
        f"- initial/final energy representation residuals: `{out.maximum_initial_energy_representation_relative_residual:.3e}` / `{out.maximum_final_energy_representation_relative_residual:.3e}`",
        f"- integrated positive/negative work representation residuals: `{out.maximum_integrated_positive_work_representation_relative_residual:.3e}` / `{out.maximum_integrated_negative_work_representation_relative_residual:.3e}`",
        f"- viscous dissipation representation residual: `{out.maximum_viscous_dissipation_representation_relative_residual:.3e}`",
    ]
    for run in out.runs:
        lines.extend([
            "",
            f"## resolution {run.resolution}",
            f"- steps: `{run.steps}`",
            f"- initial/final helical modal energy: `{run.initial_energy:.12g}` / `{run.final_energy:.12g}`",
            f"- integrated gross positive/negative work: `{run.integrated_positive_work:.12g}` / `{run.integrated_negative_work:.12g}`",
            f"- viscous dissipation: `{run.viscous_dissipation:.12g}`",
            f"- interval continuity native residual: `{run.interval_continuity_native_residual:.3e}`",
            f"- worst instantaneous signed reconstruction native residual: `{run.worst_instantaneous_signed_reconstruction_native_residual:.3e}`",
            f"- global NS energy-balance residual: `{run.global_energy_balance_relative_residual:.3e}`",
            f"- steps with positive/negative gross edge work: `{run.positive_work_steps}` / `{run.negative_work_steps}`",
        ])
    lines.extend([
        "",
        "The positive and negative terms are the gross canonical edge Hahn marginals for this physical helical child mode, not a fresh Hahn split of its summed modal net work.  The identity is aggregate stock/flow continuity only: it creates no FIFO/LIFO matching between earlier deposits and later withdrawals and does not bound total gross nonlinear transfer variation.",
        "",
    ])
    summary = "\n".join(lines)
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
