from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.full_natural_service_corridor_pde_probe import (
    _nonlinear_term,
    _rk4_step,
    _spectral_geometry,
)
from src.helical import helical_basis
from src.helical_physical_edge_registration import register_helical_physical_edge
from src.high_strain_descending_epoch_pde_probe import (
    TORUS_VOLUME,
    _divergence_norm,
    _gradient_energy,
    _physical_inner,
)
from src.signed_good_generated_epoch_pde_probe import (
    CHILD_K,
    PARENT_P,
    PARENT_Q,
    _mode_mask,
    _signed_good_triad_initial_state,
    _trapezoid,
)


STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "ACTUAL_P_PLUS_Q_SOURCE_EQUALS_SUM_OF_EIGHT_HELICAL_EDGES__"
    "SIGNED_CHILD_WORK_AND_UPPER_PROGRESS_REGISTER_ON_THE_SAME_PDE_ORBIT"
)


def _index(wavevector: tuple[int, int, int], resolution: int) -> tuple[int, int, int]:
    return tuple(int(value) % int(resolution) for value in wavevector)


def _complex_norm3(value: np.ndarray) -> float:
    q = np.asarray(value, dtype=complex)
    if q.shape != (3,) or np.any(~np.isfinite(q.real)) or np.any(~np.isfinite(q.imag)):
        raise ValueError("finite complex three-vector required")
    return float(math.hypot(*(abs(complex(x)) for x in q)))


def _relative_scalar(a: float, b: float, scale: float) -> float:
    native_scale = abs(float(scale))
    gap = abs(float(a) - float(b))
    if native_scale == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / native_scale


def _relative_vector(a: np.ndarray, b: np.ndarray, scale: float) -> float:
    native_scale = abs(float(scale))
    gap = _complex_norm3(np.asarray(a, complex) - np.asarray(b, complex))
    if native_scale == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / native_scale


def _helical_components(
    wavevector: tuple[int, int, int], value: np.ndarray
) -> dict[int, complex]:
    k = np.asarray(wavevector, dtype=float)
    v = np.asarray(value, dtype=complex)
    return {s: complex(np.vdot(helical_basis(k, s), v)) for s in (-1, 1)}


def _helical_reconstruction(
    wavevector: tuple[int, int, int], components: dict[int, complex]
) -> np.ndarray:
    k = np.asarray(wavevector, dtype=float)
    return sum(
        (components[s] * helical_basis(k, s) for s in (-1, 1)),
        start=np.zeros(3, dtype=complex),
    )


def _snapshot(
    state_hat: np.ndarray,
    viscosity: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
    parent_mask: np.ndarray,
    child_mask: np.ndarray,
) -> dict[str, float]:
    n = int(state_hat.shape[1])
    series_normalization = float(n**3)
    p_index = _index(PARENT_P, n)
    q_index = _index(PARENT_Q, n)
    z_index = _index(CHILD_K, n)

    parent_hat = parent_mask[None, ...] * state_hat
    child_hat = child_mask[None, ...] * state_hat
    designated_hat = -(
        child_mask[None, ...] * _nonlinear_term(parent_hat, k, k2, dealias)
    )
    full_nonlinear = _nonlinear_term(state_hat, k, k2, dealias)

    up = state_hat[(slice(None),) + p_index] / series_normalization
    uq = state_hat[(slice(None),) + q_index] / series_normalization
    uz = state_hat[(slice(None),) + z_index] / series_normalization
    actual_source = designated_hat[(slice(None),) + z_index] / series_normalization
    p_components = _helical_components(PARENT_P, up)
    q_components = _helical_components(PARENT_Q, uq)
    z_components = _helical_components(CHILD_K, uz)

    velocity_scale = sum(abs(a) for a in (*p_components.values(), *q_components.values()))
    parent_reconstruction_residual = max(
        _relative_vector(
            up,
            _helical_reconstruction(PARENT_P, p_components),
            sum(abs(a) for a in p_components.values()),
        ),
        _relative_vector(
            uq,
            _helical_reconstruction(PARENT_Q, q_components),
            sum(abs(a) for a in q_components.values()),
        ),
    )
    child_reconstruction_residual = _relative_vector(
        uz,
        _helical_reconstruction(CHILD_K, z_components),
        sum(abs(a) for a in z_components.values()),
    )

    source_components = {s: 0.0j for s in (-1, 1)}
    signed_mode_work = 0.0
    registered_upper = 0.0
    direct_upper = 0.0
    source_scale = 0.0
    work_scale = 0.0
    upper_scale = 0.0
    maximum_multiplier = 0.0
    for sp in (-1, 1):
        for sq in (-1, 1):
            for sz in (-1, 1):
                row = register_helical_physical_edge(
                    x=np.asarray(PARENT_P, dtype=float),
                    y=np.asarray(PARENT_Q, dtype=float),
                    z=np.asarray(CHILD_K, dtype=float),
                    sx=sp,
                    sy=sq,
                    sz=sz,
                    ax=p_components[sp],
                    ay=q_components[sq],
                    az=z_components[sz],
                )
                source_components[sz] += row.direct_child_source_coefficient
                signed_mode_work += row.signed_child_energy_work
                direct_upper += row.signed_upper_progress_work
                registered_upper += row.registered_upper_progress_work
                source_scale += row.native_source_coefficient_scale
                work_scale += abs(row.signed_child_energy_work)
                upper_scale += row.native_modal_capacity * row.geometric_multiplier_J
                maximum_multiplier = max(maximum_multiplier, row.normalized_multiplier)

    registered_source = sum(
        (
            source_components[s] * helical_basis(np.asarray(CHILD_K, float), s)
            for s in (-1, 1)
        ),
        start=np.zeros(3, dtype=complex),
    )
    actual_mode_work = 2.0 * float(np.real(np.vdot(uz, actual_source)))
    actual_full_pair_work = (
        2.0 * _physical_inner(child_hat, designated_hat, n) / TORUS_VOLUME
    )
    progress = max(
        0.0,
        math.log(
            math.sqrt(sum(value * value for value in CHILD_K))
            / math.sqrt(sum(value * value for value in PARENT_P))
        ),
    )
    actual_upper = actual_mode_work * progress

    return {
        "parent_reconstruction_residual": parent_reconstruction_residual,
        "child_reconstruction_residual": child_reconstruction_residual,
        "source_registration_residual": _relative_vector(
            actual_source, registered_source, source_scale
        ),
        "mode_work_registration_residual": _relative_scalar(
            actual_mode_work, signed_mode_work, work_scale
        ),
        "full_pair_work_registration_residual": _relative_scalar(
            actual_full_pair_work, 2.0 * signed_mode_work, 2.0 * work_scale
        ),
        "direct_upper_registration_residual": _relative_scalar(
            actual_upper, direct_upper, upper_scale
        ),
        "registered_upper_residual": _relative_scalar(
            actual_upper, registered_upper, upper_scale
        ),
        "actual_full_pair_work": actual_full_pair_work,
        "actual_source_norm": _complex_norm3(actual_source),
        "maximum_normalized_multiplier": maximum_multiplier,
        "global_energy": _physical_inner(state_hat, state_hat, n),
        "global_gradient": _gradient_energy(
            state_hat, np.ones_like(k2), k2, n
        ),
        "global_nonlinear_work": -2.0
        * _physical_inner(state_hat, full_nonlinear, n),
        "divergence_norm": _divergence_norm(state_hat, k, n),
        "child_energy": _physical_inner(child_hat, child_hat, n),
        "velocity_scale": velocity_scale,
    }


@dataclass(frozen=True)
class GalerkinHelicalEdgeRun:
    resolution: int
    steps: int
    dt: float
    duration: float
    viscosity: float
    amplitude: float
    samples: int
    positive_physical_work_samples: int
    integrated_actual_pair_work: float
    final_child_energy: float
    minimum_actual_source_norm: float
    worst_parent_reconstruction_residual: float
    worst_child_reconstruction_residual: float
    worst_source_registration_residual: float
    worst_mode_work_registration_residual: float
    worst_full_pair_work_registration_residual: float
    worst_direct_upper_registration_residual: float
    worst_registered_upper_residual: float
    maximum_normalized_multiplier: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float


def simulate_helical_edge_on_galerkin_ns(
    *,
    resolution: int,
    steps: int,
    viscosity: float = 0.02,
    amplitude: float = 64.0,
    scaled_lifetime: float = 0.05,
) -> GalerkinHelicalEdgeRun:
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    lifetime = float(scaled_lifetime)
    if n < 24 or n % 2 or count < 16:
        raise ValueError("even resolution >=24 and at least sixteen RK4 steps required")
    if not all(math.isfinite(x) and x > 0.0 for x in (nu, amp, lifetime)):
        raise ValueError("positive finite Galerkin parameters required")

    child_frequency = math.sqrt(sum(value * value for value in CHILD_K))
    duration = lifetime / (child_frequency * child_frequency)
    dt = duration / count
    k, k2, dealias = _spectral_geometry(n)
    parent_mask = _mode_mask(
        k,
        (
            PARENT_P,
            tuple(-x for x in PARENT_P),
            PARENT_Q,
            tuple(-x for x in PARENT_Q),
        ),
    )
    child_mask = _mode_mask(k, (CHILD_K, tuple(-x for x in CHILD_K)))
    state = _signed_good_triad_initial_state(n, amp, k, k2, dealias)
    observations: list[dict[str, float]] = []
    for index in range(count + 1):
        observations.append(
            _snapshot(state, nu, k, k2, dealias, parent_mask, child_mask)
        )
        if index < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    residual_names = (
        "parent_reconstruction_residual",
        "child_reconstruction_residual",
        "source_registration_residual",
        "mode_work_registration_residual",
        "full_pair_work_registration_residual",
        "direct_upper_registration_residual",
        "registered_upper_residual",
    )
    worst = {
        name: max(float(row[name]) for row in observations) for name in residual_names
    }
    if max(worst.values()) > 2.0e-9:
        raise AssertionError(f"helical edge lost the actual Galerkin NS law: {worst}")

    initial_energy = observations[0]["global_energy"]
    final_energy = observations[-1]["global_energy"]
    gradient_action = _trapezoid(
        tuple(row["global_gradient"] for row in observations), dt
    )
    balance_residual = abs(
        final_energy - initial_energy + 2.0 * nu * gradient_action
    ) / initial_energy
    nonlinear_rate_scale = initial_energy / duration
    maximum_nonlinear = max(
        abs(row["global_nonlinear_work"]) for row in observations
    ) / nonlinear_rate_scale
    maximum_divergence = max(
        row["divergence_norm"] for row in observations
    ) / math.sqrt(initial_energy)
    if balance_residual > 2.0e-5 or maximum_nonlinear > 2.0e-11 or maximum_divergence > 2.0e-11:
        raise AssertionError("helical audit trajectory lost a native Navier-Stokes invariant")

    positive_samples = sum(row["actual_full_pair_work"] > 0.0 for row in observations)
    minimum_source = min(row["actual_source_norm"] for row in observations)
    if positive_samples == 0 or minimum_source <= 0.0:
        raise AssertionError("evolved PDE supplied no nonzero positive physical edge work")

    return GalerkinHelicalEdgeRun(
        resolution=n,
        steps=count,
        dt=dt,
        duration=duration,
        viscosity=nu,
        amplitude=amp,
        samples=len(observations),
        positive_physical_work_samples=positive_samples,
        integrated_actual_pair_work=_trapezoid(
            tuple(row["actual_full_pair_work"] for row in observations), dt
        ),
        final_child_energy=observations[-1]["child_energy"],
        minimum_actual_source_norm=minimum_source,
        worst_parent_reconstruction_residual=worst["parent_reconstruction_residual"],
        worst_child_reconstruction_residual=worst["child_reconstruction_residual"],
        worst_source_registration_residual=worst["source_registration_residual"],
        worst_mode_work_registration_residual=worst["mode_work_registration_residual"],
        worst_full_pair_work_registration_residual=worst[
            "full_pair_work_registration_residual"
        ],
        worst_direct_upper_registration_residual=worst[
            "direct_upper_registration_residual"
        ],
        worst_registered_upper_residual=worst["registered_upper_residual"],
        maximum_normalized_multiplier=max(
            row["maximum_normalized_multiplier"] for row in observations
        ),
        global_energy_balance_relative_residual=balance_residual,
        maximum_global_nonlinear_work_relative_rate=maximum_nonlinear,
        maximum_divergence_relative_to_initial_l2=maximum_divergence,
    )


@dataclass(frozen=True)
class HelicalPhysicalEdgePDEProbe:
    status: str
    runs: tuple[GalerkinHelicalEdgeRun, ...]
    child_energy_resolution_spread: float
    integrated_work_resolution_spread: float


def run_probe(
    resolutions: Sequence[int] = (24, 28, 32),
    *,
    steps: int = 96,
    viscosity: float = 0.02,
    amplitude: float = 64.0,
    scaled_lifetime: float = 0.05,
) -> HelicalPhysicalEdgePDEProbe:
    runs = tuple(
        simulate_helical_edge_on_galerkin_ns(
            resolution=int(n),
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            scaled_lifetime=float(scaled_lifetime),
        )
        for n in resolutions
    )
    if not runs:
        raise ValueError("at least one Galerkin resolution required")
    energies = tuple(row.final_child_energy for row in runs)
    works = tuple(row.integrated_actual_pair_work for row in runs)
    energy_spread = (max(energies) - min(energies)) / max(energies)
    work_spread = (max(works) - min(works)) / max(abs(x) for x in works)
    if energy_spread > 1.0e-2 or work_spread > 1.0e-2:
        raise AssertionError("helical physical-edge observables did not refine stably")
    return HelicalPhysicalEdgePDEProbe(
        status=STATUS,
        runs=runs,
        child_energy_resolution_spread=energy_spread,
        integrated_work_resolution_spread=work_spread,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", type=int, nargs="+", default=(24, 28, 32))
    ap.add_argument("--steps", type=int, default=96)
    ap.add_argument("--viscosity", type=float, default=0.02)
    ap.add_argument("--amplitude", type=float, default=64.0)
    ap.add_argument("--scaled-lifetime", type=float, default=0.05)
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("results-helical-physical-edge-pde-probe"),
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = run_probe(
        args.resolutions,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        scaled_lifetime=args.scaled_lifetime,
    )
    (args.outdir / "helical_physical_edge_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    table = "\n".join(
        f"| {row.resolution} | {row.steps} | {row.worst_source_registration_residual:.3e} | "
        f"{row.worst_full_pair_work_registration_residual:.3e} | "
        f"{row.worst_registered_upper_residual:.3e} | {row.global_energy_balance_relative_residual:.3e} |"
        for row in result.runs
    )
    summary = f"""# Helical edge on actual NS

`p+q=z` is read from one evolved dealiased Fourier--Galerkin NS orbit and split
into all eight helical sign edges. Their sum reconstructs the PDE source, signed
child work, and `T log_+=A J c` on that same orbit.

| n | steps | source | work | upper | NS energy balance |
|---:|---:|---:|---:|---:|---:|
{table}

This audits one resolved physical triad; it is not a continuum or generic-HH proof.
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
