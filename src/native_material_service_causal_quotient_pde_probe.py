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
from src.mixed_fate_reserved_young_handoff_pde_probe import adversarial_mixed_fate_initial_state
from src.native_material_service_causal_quotient import (
    PositiveMaterialServiceQuotientCertificate,
    material_ownership_rereading_anti_theorem,
)


STATUS = (
    "ACTUAL_DEALIASED_FOURIER_GALERKIN_NS_INCREMENT_SERVICE__"
    "SAME_PHYSICAL_LAW_DIFFERENT_MATERIAL_READINGS__"
    "MATERIAL_INTERFACE_MASS_IS_NOT_GENERATION"
)


@dataclass(frozen=True)
class ActualNSMaterialServiceObservation:
    resolution: int
    time: float
    physical_displacement_fraction: float
    total_increment_service: float
    x_material_old_old_service: float
    x_material_old_new_service: float
    x_material_new_new_service: float
    y_material_old_old_service: float
    y_material_old_new_service: float
    y_material_new_new_service: float
    ownership_partition_l1_change: float
    service_rereading_residual: float
    global_nonlinear_work_rate: float
    divergence_relative_to_initial_l2: float

    def __post_init__(self) -> None:
        vals = (
            self.time,
            self.total_increment_service,
            self.x_material_old_old_service,
            self.x_material_old_new_service,
            self.x_material_new_new_service,
            self.y_material_old_old_service,
            self.y_material_old_new_service,
            self.y_material_new_new_service,
            self.ownership_partition_l1_change,
            self.service_rereading_residual,
            abs(self.global_nonlinear_work_rate),
            self.divergence_relative_to_initial_l2,
        )
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in vals):
            raise ValueError("finite nonnegative actual NS material-service diagnostics required")
        if self.total_increment_service <= 0.0:
            raise AssertionError("actual NS snapshot carried no positive velocity-increment service")
        if self.service_rereading_residual > 2e-12 * max(1.0, self.total_increment_service):
            raise AssertionError("changing only material reading changed actual NS increment service")
        if self.y_material_old_new_service > 2e-12 * max(1.0, self.total_increment_service):
            raise AssertionError("y-slab material set acquired x-displacement interface service")


@dataclass(frozen=True)
class ActualNSMaterialServiceRun:
    resolution: int
    cutoff: int
    steps: int
    duration: float
    viscosity: float
    amplitude: float
    phase_sign: int
    observations: tuple[ActualNSMaterialServiceObservation, ...]
    global_energy_balance_relative_residual: float
    minimum_x_interface_service: float
    maximum_service_rereading_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float

    def __post_init__(self) -> None:
        if not self.observations:
            raise ValueError("actual NS material-service run requires observations")
        if self.global_energy_balance_relative_residual > 5e-5:
            raise AssertionError("actual NS global energy balance left certified tolerance")
        if self.minimum_x_interface_service <= 0.0:
            raise AssertionError("actual NS x-material boundary never carried positive increment service")
        if self.maximum_service_rereading_relative_residual > 5e-12:
            raise AssertionError("material rereading altered the underlying actual NS service law")
        if self.maximum_global_nonlinear_work_relative_rate > 5e-10:
            raise AssertionError("actual Galerkin NS nonlinearity lost global energy conservation")
        if self.maximum_divergence_relative_to_initial_l2 > 5e-11:
            raise AssertionError("actual Galerkin NS trajectory lost incompressibility")


@dataclass(frozen=True)
class NativeMaterialServicePDEProbe:
    status: str
    runs: tuple[ActualNSMaterialServiceRun, ...]
    common_cutoff: int
    phase_sign: int
    maximum_cross_representation_service_relative_residual: float
    minimum_material_partition_change_fraction: float

    def __post_init__(self) -> None:
        if self.status != STATUS or not self.runs:
            raise ValueError("actual NS material-service referee provenance mismatch")
        if self.maximum_cross_representation_service_relative_residual > 3e-8:
            raise AssertionError("same-cutoff physical increment service changed across FFT representations")
        if self.minimum_material_partition_change_fraction <= 1e-6:
            raise AssertionError("actual NS anti-theorem failed to change material partition at fixed service")


def _physical_velocity(state_hat: np.ndarray) -> np.ndarray:
    return np.fft.ifftn(np.asarray(state_hat, dtype=complex), axes=(1, 2, 3))


def _material_slabs(n: int) -> tuple[np.ndarray, np.ndarray]:
    idx = np.indices((n, n, n), dtype=int)
    return idx[0] < n // 2, idx[1] < n // 2


def _snapshot_material_service(
    state_hat: np.ndarray,
    nonlinear_hat: np.ndarray,
    *,
    time: float,
    k: np.ndarray,
    initial_l2: float,
) -> ActualNSMaterialServiceObservation:
    n = int(state_hat.shape[1])
    if n % 4:
        raise ValueError("actual NS material-service referee requires grid divisible by four")
    # A quarter-period x translation is the same physical displacement pi/2 on
    # every representation.  The service density is the actual positive square
    # velocity increment |u(x-r)-u(x)|^2 of the evolving NS field.
    shift = n // 4
    u = _physical_velocity(state_hat)
    nbr = np.roll(u, shift=shift, axis=1)
    density = np.sum(np.abs(nbr - u) ** 2, axis=0).real
    weights = density.reshape(-1)
    x_old, y_old = _material_slabs(n)
    x_neighbor = np.roll(x_old, shift=shift, axis=0)
    y_neighbor = np.roll(y_old, shift=shift, axis=0)
    anti = material_ownership_rereading_anti_theorem(
        service_measure="actual_NS_quarter_period_velocity_increment_square_service",
        native_owner=None,
        edge_weights=weights,
        first_old_here=x_old.reshape(-1),
        first_old_neighbor=x_neighbor.reshape(-1),
        second_old_here=y_old.reshape(-1),
        second_old_neighbor=y_neighbor.reshape(-1),
    )
    qa = anti["first"]
    qb = anti["second"]
    if not isinstance(qa, PositiveMaterialServiceQuotientCertificate) or not isinstance(qb, PositiveMaterialServiceQuotientCertificate):
        raise AssertionError("actual NS material-service referee lost typed quotient")
    norm = float(n**3)
    total = qa.total_service / norm
    residual = abs(qa.total_service - qb.total_service) / norm
    nonlinear_work = -2.0 * _spectral_average_inner(state_hat, nonlinear_hat, n)
    divergence = _divergence_norm(state_hat, k, n) / math.sqrt(initial_l2)
    return ActualNSMaterialServiceObservation(
        resolution=n,
        time=float(time),
        physical_displacement_fraction=0.25,
        total_increment_service=total,
        x_material_old_old_service=qa.old_old_service / norm,
        x_material_old_new_service=qa.old_new_interface_service / norm,
        x_material_new_new_service=qa.new_new_service / norm,
        y_material_old_old_service=qb.old_old_service / norm,
        y_material_old_new_service=qb.old_new_interface_service / norm,
        y_material_new_new_service=qb.new_new_service / norm,
        ownership_partition_l1_change=float(anti["ownership_partition_l1_change"]) / norm,
        service_rereading_residual=residual,
        global_nonlinear_work_rate=nonlinear_work,
        divergence_relative_to_initial_l2=divergence,
    )


def _run_one(
    *,
    resolution: int,
    cutoff: int,
    steps: int,
    viscosity: float,
    amplitude: float,
    duration: float,
    phase_sign: int,
) -> ActualNSMaterialServiceRun:
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    sign = int(phase_sign)
    if n < 24 or n % 4:
        raise ValueError("actual NS material-service referee requires FFT grid >=24 divisible by four")
    if count < 16 or nu <= 0.0 or amp <= 0.0 or horizon <= 0.0 or sign not in (-1, 1):
        raise ValueError("valid positive Galerkin NS referee parameters required")
    k, k2, dealias, actual_cutoff = _spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("common Galerkin cutoff changed")
    state = sign * adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=amp)
    initial_l2 = _spectral_average_inner(state, state, n)
    dt = horizon / count
    times = [j * dt for j in range(count + 1)]
    energy: list[float] = []
    gradient: list[float] = []
    nonlinear_work: list[float] = []
    divergence: list[float] = []
    obs: list[ActualNSMaterialServiceObservation] = []
    for j, t in enumerate(times):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        energy.append(_spectral_average_inner(state, state, n))
        gradient.append(_gradient_energy(state, k2, n))
        nonlinear_work.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        divergence.append(_divergence_norm(state, k, n))
        # Sample several genuine trajectory times, including both endpoints.
        if j in {0, count // 3, (2 * count) // 3, count}:
            obs.append(_snapshot_material_service(state, nonlinear, time=t, k=k, initial_l2=initial_l2))
        if j < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)
    global_balance = abs(energy[-1] - energy[0] + 2.0 * nu * _trapezoid(gradient, times)) / energy[0]
    nonlinear_scale = energy[0] / horizon
    max_nonlinear = max(abs(x) for x in nonlinear_work) / nonlinear_scale
    max_divergence = max(divergence) / math.sqrt(energy[0])
    return ActualNSMaterialServiceRun(
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        duration=horizon,
        viscosity=nu,
        amplitude=amp,
        phase_sign=sign,
        observations=tuple(obs),
        global_energy_balance_relative_residual=global_balance,
        minimum_x_interface_service=min(o.x_material_old_new_service for o in obs),
        maximum_service_rereading_relative_residual=max(
            o.service_rereading_residual / max(o.total_increment_service, 1e-300) for o in obs
        ),
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
        maximum_divergence_relative_to_initial_l2=max_divergence,
    )


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    cutoff: int = 7,
    steps: int = 32,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    phase_sign: int = 1,
) -> NativeMaterialServicePDEProbe:
    runs = tuple(
        _run_one(
            resolution=int(n),
            cutoff=int(cutoff),
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            duration=float(duration),
            phase_sign=int(phase_sign),
        )
        for n in resolutions
    )
    if not runs:
        raise ValueError("at least one actual NS representation required")
    if len({r.cutoff for r in runs}) != 1:
        raise ValueError("actual NS representations must share one Galerkin cutoff")
    sample_count = len(runs[0].observations)
    if any(len(r.observations) != sample_count for r in runs):
        raise AssertionError("cross-representation NS referee sampled different physical times")
    worst = 0.0
    for q in range(sample_count):
        values = [r.observations[q].total_increment_service for r in runs]
        scale = max(max(abs(v) for v in values), 1e-300)
        worst = max(worst, (max(values) - min(values)) / scale)
    change_fraction = min(
        o.ownership_partition_l1_change / max(o.total_increment_service, 1e-300)
        for r in runs
        for o in r.observations
    )
    return NativeMaterialServicePDEProbe(
        status=STATUS,
        runs=runs,
        common_cutoff=runs[0].cutoff,
        phase_sign=int(phase_sign),
        maximum_cross_representation_service_relative_residual=worst,
        minimum_material_partition_change_fraction=change_fraction,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--resolutions", type=int, nargs="+", default=(24, 28))
    ap.add_argument("--cutoff", type=int, default=7)
    ap.add_argument("--steps", type=int, default=32)
    ap.add_argument("--viscosity", type=float, default=0.03)
    ap.add_argument("--amplitude", type=float, default=1.0)
    ap.add_argument("--duration", type=float, default=0.001)
    ap.add_argument("--phase-sign", type=int, default=1)
    ap.add_argument("--outdir", type=Path, default=Path("results-native-material-service-causal-quotient-ns"))
    args = ap.parse_args()
    out = run_probe(
        resolutions=args.resolutions,
        cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        phase_sign=args.phase_sign,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "native_material_service_causal_quotient_pde.json").write_text(
        json.dumps(asdict(out), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "# Actual Fourier--Galerkin Navier--Stokes referee: material ownership is downstream of service",
        "",
        f"Status: **{STATUS}**.",
        "",
        "The referee evolves the same dealiased incompressible Galerkin Navier--Stokes state used by the physical helical audits. At genuine trajectory times it forms one actual positive velocity-increment square law at the fixed physical displacement r=(pi/2,0,0). The state and service law are then held literally fixed while two geometric old-material readings are applied: an x-slab, whose endpoints cross the displacement boundary, and a y-slab, whose endpoints do not. The OO/ON/NN masses change, but the underlying service and Navier--Stokes trajectory do not.",
        "",
        f"- maximum cross-FFT representation service residual: `{out.maximum_cross_representation_service_relative_residual:.3e}`",
        f"- minimum OO/ON/NN partition-change fraction at fixed service: `{out.minimum_material_partition_change_fraction:.12g}`",
    ]
    for run in out.runs:
        lines.extend([
            "",
            f"## resolution {run.resolution}",
            f"- global NS energy-balance residual: `{run.global_energy_balance_relative_residual:.3e}`",
            f"- minimum positive x-boundary interface service: `{run.minimum_x_interface_service:.12g}`",
            f"- maximum material-rereading service residual: `{run.maximum_service_rereading_relative_residual:.3e}`",
            f"- maximum global nonlinear-work relative rate: `{run.maximum_global_nonlinear_work_relative_rate:.3e}`",
            f"- maximum divergence / initial L2: `{run.maximum_divergence_relative_to_initial_l2:.3e}`",
        ])
        for o in run.observations:
            lines.append(
                f"- t={o.time:.9g}: total service `{o.total_increment_service:.12g}`, x-ON `{o.x_material_old_new_service:.12g}`, y-ON `{o.y_material_old_new_service:.12g}`, partition L1 change `{o.ownership_partition_l1_change:.12g}`"
            )
    lines.extend([
        "",
        "This referee does not invent a material PDE, does not infer a causal supplier for the snapshot increment law, and does not claim that the chosen geometric slabs are persistent material packets. Its point is narrower and exact: for an actual NS field, changing only the material reading can change interface/fresh provenance while the positive physical service law is unchanged. Therefore material ownership alone cannot be its generator. No global-regularity claim is made.",
    ])
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
