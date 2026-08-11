from __future__ import annotations

import argparse
import json
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.continuum_helical_edge_measure_registration import (
    continuum_edge_measure_ledger,
    ordered_parent_curl_source,
    register_continuum_triad_fiber,
    unitary_fourier_convolution_factor,
    unordered_parent_curl_source_vector,
)
from src.helical_physical_edge_registration import leray_project

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "FULL_RETAINED_PARENT_CONVOLUTION_QUOTIENTED_INTO_UNORDERED_EDGES__"
    "EIGHT_HELICITIES_RECONSTRUCT_SIGNED_WORK_HAHN_AND_PROGRESS_ON_ONE_PDE_ORBIT"
)

TORUS_VOLUME = (2.0 * math.pi) ** 3
CHILD_MODE = (5, 1, 0)


def _index(wavevector: tuple[int, int, int], resolution: int) -> tuple[int, int, int]:
    return tuple(int(value) % int(resolution) for value in wavevector)


def _norm3(value: np.ndarray) -> float:
    q = np.asarray(value, dtype=complex)
    if q.shape != (3,) or np.any(~np.isfinite(q.real)) or np.any(~np.isfinite(q.imag)):
        raise ValueError("finite complex three-vector required")
    return float(math.hypot(*(abs(complex(x)) for x in q)))


def _relative_scalar(actual: float, expected: float, scale: float) -> float:
    native = abs(float(scale))
    gap = abs(float(actual) - float(expected))
    if native == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / native


def _relative_vector(actual: np.ndarray, expected: np.ndarray, scale: float) -> float:
    native = abs(float(scale))
    gap = _norm3(np.asarray(actual, complex) - np.asarray(expected, complex))
    if native == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / native


def _spectral_geometry(
    resolution: int, cutoff_override: int | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    n = int(resolution)
    if n < 20 or n % 2:
        raise ValueError("an even Fourier resolution at least 20 is required")
    one = np.fft.fftfreq(n, d=1.0 / n)
    k = np.asarray(np.meshgrid(one, one, one, indexing="ij"), dtype=float)
    k2 = np.sum(k * k, axis=0)
    native_cutoff = n // 3 - 1
    cutoff = native_cutoff if cutoff_override is None else int(cutoff_override)
    if cutoff <= 0 or cutoff > native_cutoff:
        raise ValueError("Galerkin cutoff must be positive and no larger than the native dealiased cutoff")
    dealias = np.max(np.abs(k), axis=0) <= cutoff
    if max(abs(v) for v in CHILD_MODE) > cutoff:
        raise ValueError("selected child mode is outside the dealiased Galerkin set")
    return k, k2, dealias, cutoff


def _leray_dealias(
    field_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    out = np.asarray(field_hat, dtype=complex) * dealias[None, ...]
    inverse = np.zeros_like(k2)
    nonzero = k2 > 0.0
    inverse[nonzero] = 1.0 / k2[nonzero]
    longitudinal = np.sum(k * out, axis=0) * inverse
    out = out - k * longitudinal[None, ...]
    out[:, 0, 0, 0] = 0.0
    return out * dealias[None, ...]


def _nonlinear_term(
    state_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Leray-projected 2/3-dealiased Fourier-Galerkin ``(u.grad)u``."""
    velocity = np.fft.ifftn(state_hat, axes=(1, 2, 3)).real
    convection = np.zeros_like(velocity)
    for component in range(3):
        for direction in range(3):
            derivative = np.fft.ifftn(
                1j * k[direction] * state_hat[component],
                axes=(0, 1, 2),
            ).real
            convection[component] += velocity[direction] * derivative
    return _leray_dealias(
        np.fft.fftn(convection, axes=(1, 2, 3)),
        k,
        k2,
        dealias,
    )


def _rhs(
    state_hat: np.ndarray,
    viscosity: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    nonlinear = _nonlinear_term(state_hat, k, k2, dealias)
    return _leray_dealias(
        -nonlinear - viscosity * k2[None, ...] * state_hat,
        k,
        k2,
        dealias,
    )


def _rk4_step(
    state_hat: np.ndarray,
    dt: float,
    viscosity: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    k1 = _rhs(state_hat, viscosity, k, k2, dealias)
    k2s = _rhs(state_hat + 0.5 * dt * k1, viscosity, k, k2, dealias)
    k3 = _rhs(state_hat + 0.5 * dt * k2s, viscosity, k, k2, dealias)
    k4 = _rhs(state_hat + dt * k3, viscosity, k, k2, dealias)
    return _leray_dealias(
        state_hat + (dt / 6.0) * (k1 + 2.0 * k2s + 2.0 * k3 + k4),
        k,
        k2,
        dealias,
    )


def _canonical_half(wavevector: tuple[int, int, int]) -> bool:
    for value in wavevector:
        if value:
            return value > 0
    return False


def _deterministic_smooth_initial_state(
    resolution: int,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
    amplitude: float,
    band: int = 5,
) -> np.ndarray:
    """One fixed real smooth divergence-free Fourier polynomial on every grid."""
    n = int(resolution)
    coeff = np.zeros((3, n, n, n), dtype=complex)
    for a in range(-band, band + 1):
        for b in range(-band, band + 1):
            for c in range(-band, band + 1):
                q = (a, b, c)
                if q == (0, 0, 0) or not _canonical_half(q):
                    continue
                idx = _index(q, n)
                if not bool(dealias[idx]):
                    continue
                qv = np.asarray(q, dtype=float)
                phase = 0.173 * a - 0.119 * b + 0.071 * c
                raw = np.asarray(
                    (
                        math.sin(0.37 * a + 0.11 * b) + 1j * math.cos(0.23 * c + phase),
                        math.cos(0.29 * b - 0.07 * c) + 1j * math.sin(0.31 * a - phase),
                        math.sin(0.19 * c + 0.13 * a) + 1j * math.cos(0.17 * b + phase),
                    ),
                    dtype=complex,
                )
                projected = leray_project(qv, raw)
                decay = math.exp(-0.18 * float(np.dot(qv, qv)))
                value = decay * projected
                coeff[(slice(None),) + idx] = value
                coeff[(slice(None),) + _index(tuple(-v for v in q), n)] = np.conjugate(value)

    energy = float(np.vdot(coeff, coeff).real)
    if not math.isfinite(energy) or energy <= 0.0:
        raise AssertionError("deterministic Galerkin initial data lost positive energy")
    coeff *= float(amplitude) / math.sqrt(energy)
    # numpy's inverse FFT divides by n^3.  Store FFT coefficients so ``coeff`` is
    # exactly the common Fourier-series data on every audit resolution.
    state_hat = coeff * float(n**3)
    return _leray_dealias(state_hat, k, k2, dealias)


def _series_coefficient(state_hat: np.ndarray, wavevector: tuple[int, int, int]) -> np.ndarray:
    n = int(state_hat.shape[1])
    return np.asarray(state_hat[(slice(None),) + _index(wavevector, n)], dtype=complex) / float(n**3)


def _retained_modes(cutoff: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (a, b, c)
        for a in range(-cutoff, cutoff + 1)
        for b in range(-cutoff, cutoff + 1)
        for c in range(-cutoff, cutoff + 1)
        if (a, b, c) != (0, 0, 0)
    )


def _pair_orbits_for_child(
    child: tuple[int, int, int], cutoff: int
) -> tuple[tuple[tuple[int, int, int], tuple[int, int, int]], ...]:
    retained = set(_retained_modes(cutoff))
    seen: set[frozenset[tuple[int, int, int]]] = set()
    out: list[tuple[tuple[int, int, int], tuple[int, int, int]]] = []
    for x in retained:
        y = tuple(child[i] - x[i] for i in range(3))
        if y == (0, 0, 0) or y not in retained:
            continue
        if x == y:
            raise AssertionError("PDE audit child unexpectedly contains a discrete fixed parent orbit")
        orbit = frozenset((x, y))
        if orbit in seen:
            continue
        seen.add(orbit)
        out.append((x, y))
    if not out:
        raise AssertionError("selected child has no retained unordered parent orbits")
    return tuple(out)


def _spectral_average_inner(left: np.ndarray, right: np.ndarray, resolution: int) -> float:
    return float(np.vdot(left, right).real / float(resolution**6))


def _gradient_energy(state_hat: np.ndarray, k2: np.ndarray, resolution: int) -> float:
    weighted = np.sqrt(k2)[None, ...] * state_hat
    return _spectral_average_inner(weighted, weighted, resolution)


def _divergence_norm(state_hat: np.ndarray, k: np.ndarray, resolution: int) -> float:
    divergence_hat = 1j * np.sum(k * state_hat, axis=0)
    return math.sqrt(max(0.0, float(np.vdot(divergence_hat, divergence_hat).real / float(resolution**6))))


def _snapshot(
    state_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
    cutoff: int,
) -> dict[str, float]:
    n = int(state_hat.shape[1])
    child = CHILD_MODE
    z = np.asarray(child, dtype=float)
    uz = _series_coefficient(state_hat, child)
    nonlinear = _nonlinear_term(state_hat, k, k2, dealias)
    actual_source = -np.asarray(
        nonlinear[(slice(None),) + _index(child, n)], dtype=complex
    ) / float(n**3)

    ordered_source = np.zeros(3, dtype=complex)
    ordered_source_scale = 0.0
    unordered_source = np.zeros(3, dtype=complex)
    unordered_source_scale = 0.0
    fibers = []
    direct_progress = 0.0
    direct_progress_scale = 0.0
    direct_pair_work_scale = 0.0
    discrete_qmass = 1.0 / unitary_fourier_convolution_factor()

    # Ordered convolution is checked independently of the unordered quotient.
    retained = set(_retained_modes(cutoff))
    for x in retained:
        y = tuple(child[i] - x[i] for i in range(3))
        if y == (0, 0, 0) or y not in retained:
            continue
        ux = _series_coefficient(state_hat, x)
        uy = _series_coefficient(state_hat, y)
        term = ordered_parent_curl_source(
            np.asarray(x, float), np.asarray(y, float), z, ux, uy
        )
        ordered_source += term
        ordered_source_scale += _norm3(term)

    for x, y in _pair_orbits_for_child(child, cutoff):
        ux = _series_coefficient(state_hat, x)
        uy = _series_coefficient(state_hat, y)
        pair_source = unordered_parent_curl_source_vector(
            np.asarray(x, float), np.asarray(y, float), z, ux, uy
        )
        unordered_source += pair_source
        unordered_source_scale += _norm3(pair_source)
        fiber = register_continuum_triad_fiber(
            x=np.asarray(x, float),
            y=np.asarray(y, float),
            z=z,
            ux=ux,
            uy=uy,
            uz=uz,
            quotient_measure_mass=discrete_qmass,
        )
        fibers.append(fiber)
        direct_pair_work_scale += abs(fiber.direct_signed_work_density)
        direct_progress += fiber.direct_signed_progress_density
        direct_progress_scale += abs(fiber.direct_signed_progress_density)

    ledger = continuum_edge_measure_ledger(tuple(fibers))
    actual_work = 2.0 * float(np.real(np.vdot(uz, actual_source)))
    source_scale = max(ordered_source_scale, unordered_source_scale, _norm3(actual_source))
    work_scale = max(direct_pair_work_scale, abs(actual_work), ledger.positive_edge_work + ledger.negative_edge_work)
    progress_scale = max(direct_progress_scale, abs(direct_progress), abs(ledger.signed_registered_progress))

    global_energy = _spectral_average_inner(state_hat, state_hat, n)
    global_gradient = _gradient_energy(state_hat, k2, n)
    global_nonlinear_work = -2.0 * _spectral_average_inner(state_hat, nonlinear, n)

    return {
        "unordered_pairs": float(len(fibers)),
        "modal_edges": float(ledger.modal_edges),
        "actual_source_norm": _norm3(actual_source),
        "ordered_source_residual": _relative_vector(actual_source, ordered_source, source_scale),
        "unordered_source_residual": _relative_vector(actual_source, unordered_source, source_scale),
        "ordered_unordered_residual": _relative_vector(ordered_source, unordered_source, source_scale),
        "signed_work_residual": _relative_scalar(actual_work, ledger.signed_direct_work, work_scale),
        "signed_modal_work_residual": _relative_scalar(actual_work, ledger.signed_modal_work, work_scale),
        "progress_residual": _relative_scalar(direct_progress, ledger.signed_registered_progress, progress_scale),
        "hahn_residual": _relative_scalar(
            ledger.positive_edge_work - ledger.negative_edge_work,
            ledger.signed_direct_work,
            ledger.positive_edge_work + ledger.negative_edge_work,
        ),
        "positive_edge_work": ledger.positive_edge_work,
        "negative_edge_work": ledger.negative_edge_work,
        "aggregate_positive_work": ledger.aggregate_positive_work,
        "fiber_positive_work": ledger.fiber_positive_work,
        "positive_forward_work": ledger.positive_forward_work,
        "positive_nonforward_work": ledger.positive_nonforward_work,
        "block_transfer_deficit": ledger.block_transfer_deficit,
        "capacity_mass": ledger.capacity_mass,
        "actual_child_work": actual_work,
        "global_energy": global_energy,
        "global_gradient": global_gradient,
        "global_nonlinear_work": global_nonlinear_work,
        "divergence_norm": _divergence_norm(state_hat, k, n),
        "child_energy": float(np.vdot(uz, uz).real),
    }


@dataclass(frozen=True)
class GalerkinContinuumEdgeRun:
    resolution: int
    spectral_cutoff: int
    steps: int
    snapshots: int
    dt: float
    duration: float
    viscosity: float
    amplitude: float
    unordered_pairs: int
    modal_edges: int
    nonzero_source_snapshots: int
    positive_child_work_snapshots: int
    positive_nonforward_snapshots: int
    worst_ordered_source_residual: float
    worst_unordered_source_residual: float
    worst_ordered_unordered_residual: float
    worst_signed_work_residual: float
    worst_signed_modal_work_residual: float
    worst_progress_residual: float
    worst_hahn_residual: float
    maximum_block_transfer_deficit: float
    minimum_block_transfer_deficit: float
    integrated_child_work: float
    final_child_energy: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float


def _trapezoid(values: Sequence[float], times: Sequence[float]) -> float:
    if len(values) != len(times) or not values:
        raise ValueError("matching nonempty trapezoid data required")
    return float(sum(0.5 * (values[i] + values[i + 1]) * (times[i + 1] - times[i]) for i in range(len(values) - 1)))


def simulate_continuum_edge_measure_on_galerkin_ns(
    *,
    resolution: int,
    steps: int,
    viscosity: float = 0.03,
    amplitude: float = 4.0,
    duration: float = 0.012,
    snapshot_count: int = 5,
    cutoff_override: int | None = None,
) -> GalerkinContinuumEdgeRun:
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    horizon = float(duration)
    snaps = int(snapshot_count)
    if count < 16 or snaps < 3 or snaps > count + 1:
        raise ValueError("at least sixteen RK4 steps and three valid snapshots required")
    if not all(math.isfinite(x) and x > 0.0 for x in (nu, amp, horizon)):
        raise ValueError("positive finite Galerkin audit parameters required")

    k, k2, dealias, cutoff = _spectral_geometry(n, cutoff_override)
    state = _deterministic_smooth_initial_state(n, k, k2, dealias, amp)
    dt = horizon / count
    sample_indices = tuple(sorted({round(j * count / (snaps - 1)) for j in range(snaps)}))
    observations: list[dict[str, float]] = []
    times: list[float] = []
    all_energy: list[float] = []
    all_gradient: list[float] = []
    all_nonlinear: list[float] = []
    all_divergence: list[float] = []

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        all_energy.append(_spectral_average_inner(state, state, n))
        all_gradient.append(_gradient_energy(state, k2, n))
        all_nonlinear.append(-2.0 * _spectral_average_inner(state, nonlinear, n))
        all_divergence.append(_divergence_norm(state, k, n))
        if step in sample_indices:
            observations.append(_snapshot(state, k, k2, dealias, cutoff))
            times.append(step * dt)
        if step < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    residual_fields = (
        "ordered_source_residual",
        "unordered_source_residual",
        "ordered_unordered_residual",
        "signed_work_residual",
        "signed_modal_work_residual",
        "progress_residual",
        "hahn_residual",
    )
    worst = {name: max(float(row[name]) for row in observations) for name in residual_fields}
    if max(worst.values()) > 3.0e-8:
        raise AssertionError(f"continuum edge measure lost the actual Galerkin NS convolution: {worst}")

    initial_energy = all_energy[0]
    final_energy = all_energy[-1]
    grid_times = tuple(i * dt for i in range(count + 1))
    gradient_action = _trapezoid(all_gradient, grid_times)
    balance = abs(final_energy - initial_energy + 2.0 * nu * gradient_action) / initial_energy
    nonlinear_scale = initial_energy / horizon
    max_nonlinear = max(abs(x) for x in all_nonlinear) / nonlinear_scale
    max_divergence = max(all_divergence) / math.sqrt(initial_energy)
    if balance > 5.0e-5 or max_nonlinear > 5.0e-10 or max_divergence > 5.0e-11:
        raise AssertionError("continuum edge audit trajectory lost a native Navier-Stokes invariant")

    nonzero_sources = sum(float(row["actual_source_norm"]) > 0.0 for row in observations)
    if nonzero_sources != len(observations):
        raise AssertionError("selected child lost its actual nonlinear source on the audit orbit")

    child_works = tuple(float(row["actual_child_work"]) for row in observations)
    return GalerkinContinuumEdgeRun(
        resolution=n,
        spectral_cutoff=cutoff,
        steps=count,
        snapshots=len(observations),
        dt=dt,
        duration=horizon,
        viscosity=nu,
        amplitude=amp,
        unordered_pairs=int(observations[0]["unordered_pairs"]),
        modal_edges=int(observations[0]["modal_edges"]),
        nonzero_source_snapshots=nonzero_sources,
        positive_child_work_snapshots=sum(x > 0.0 for x in child_works),
        positive_nonforward_snapshots=sum(float(row["positive_nonforward_work"]) > 0.0 for row in observations),
        worst_ordered_source_residual=worst["ordered_source_residual"],
        worst_unordered_source_residual=worst["unordered_source_residual"],
        worst_ordered_unordered_residual=worst["ordered_unordered_residual"],
        worst_signed_work_residual=worst["signed_work_residual"],
        worst_signed_modal_work_residual=worst["signed_modal_work_residual"],
        worst_progress_residual=worst["progress_residual"],
        worst_hahn_residual=worst["hahn_residual"],
        maximum_block_transfer_deficit=max(float(row["block_transfer_deficit"]) for row in observations),
        minimum_block_transfer_deficit=min(float(row["block_transfer_deficit"]) for row in observations),
        integrated_child_work=_trapezoid(child_works, times),
        final_child_energy=float(observations[-1]["child_energy"]),
        global_energy_balance_relative_residual=balance,
        maximum_global_nonlinear_work_relative_rate=max_nonlinear,
        maximum_divergence_relative_to_initial_l2=max_divergence,
    )


@dataclass(frozen=True)
class ContinuumEdgeMeasurePDEProbe:
    status: str
    child_mode: tuple[int, int, int]
    runs: tuple[GalerkinContinuumEdgeRun, ...]
    native_final_child_energy_resolution_spread: float
    native_integrated_child_work_resolution_spread: float
    common_cutoff: int
    common_cutoff_runs: tuple[GalerkinContinuumEdgeRun, ...]
    common_final_child_energy_resolution_spread: float
    common_integrated_child_work_resolution_spread: float


def _resolution_spreads(
    runs: Sequence[GalerkinContinuumEdgeRun],
) -> tuple[float, float]:
    if not runs:
        raise ValueError("at least one Galerkin run required")
    energies = tuple(row.final_child_energy for row in runs)
    works = tuple(row.integrated_child_work for row in runs)
    energy_scale = max(abs(x) for x in energies)
    work_scale = max(abs(x) for x in works)
    energy_spread = 0.0 if energy_scale == 0.0 else (max(energies) - min(energies)) / energy_scale
    work_spread = 0.0 if work_scale == 0.0 else (max(works) - min(works)) / work_scale
    return energy_spread, work_spread


def run_probe(
    resolutions: Sequence[int] = (20, 24, 28),
    *,
    steps: int = 64,
    viscosity: float = 0.03,
    amplitude: float = 4.0,
    duration: float = 0.012,
    snapshot_count: int = 5,
) -> ContinuumEdgeMeasurePDEProbe:
    resolved = tuple(int(n) for n in resolutions)
    if not resolved:
        raise ValueError("at least one Galerkin resolution required")

    # These are genuinely different Galerkin truncations.  Every run must obey
    # the NS convolution/work identities, but cross-resolution convergence is a
    # separate numerical-analysis statement and is therefore reported only.
    runs = tuple(
        simulate_continuum_edge_measure_on_galerkin_ns(
            resolution=n,
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            duration=float(duration),
            snapshot_count=int(snapshot_count),
        )
        for n in resolved
    )
    native_energy_spread, native_work_spread = _resolution_spreads(runs)

    # Now embed exactly the same finite Galerkin system on every FFT grid.  This
    # is a representation-invariance test, not a convergence assumption.  The
    # common cutoff is the largest cutoff admissible on every requested grid.
    common_cutoff = min(n // 3 - 1 for n in resolved)
    if common_cutoff < max(abs(v) for v in CHILD_MODE):
        raise ValueError("requested resolutions do not share a Galerkin cutoff containing the audit child")
    common_runs = tuple(
        simulate_continuum_edge_measure_on_galerkin_ns(
            resolution=n,
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            duration=float(duration),
            snapshot_count=int(snapshot_count),
            cutoff_override=common_cutoff,
        )
        for n in resolved
    )
    common_energy_spread, common_work_spread = _resolution_spreads(common_runs)
    if common_energy_spread > 5.0e-8 or common_work_spread > 5.0e-8:
        raise AssertionError(
            "the same Galerkin Navier-Stokes system changed under FFT-grid representation"
        )

    return ContinuumEdgeMeasurePDEProbe(
        status=STATUS,
        child_mode=CHILD_MODE,
        runs=runs,
        native_final_child_energy_resolution_spread=native_energy_spread,
        native_integrated_child_work_resolution_spread=native_work_spread,
        common_cutoff=common_cutoff,
        common_cutoff_runs=common_runs,
        common_final_child_energy_resolution_spread=common_energy_spread,
        common_integrated_child_work_resolution_spread=common_work_spread,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", nargs="+", type=int, default=(20, 24, 28))
    ap.add_argument("--steps", type=int, default=64)
    ap.add_argument("--viscosity", type=float, default=0.03)
    ap.add_argument("--amplitude", type=float, default=4.0)
    ap.add_argument("--duration", type=float, default=0.012)
    ap.add_argument("--snapshots", type=int, default=5)
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("results-continuum-helical-edge-measure-pde-probe"),
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = run_probe(
        args.resolutions,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        snapshot_count=args.snapshots,
    )
    (args.outdir / "continuum_helical_edge_measure_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    rows = "\n".join(
        f"| {row.resolution} | {row.spectral_cutoff} | {row.unordered_pairs} | {row.modal_edges} | "
        f"{row.worst_unordered_source_residual:.3e} | {row.worst_signed_work_residual:.3e} | "
        f"{row.worst_progress_residual:.3e} | {row.global_energy_balance_relative_residual:.3e} |"
        for row in result.runs
    )
    common_rows = "\n".join(
        f"| {row.resolution} | {row.spectral_cutoff} | {row.unordered_pairs} | {row.modal_edges} | "
        f"{row.worst_unordered_source_residual:.3e} | {row.worst_signed_work_residual:.3e} | "
        f"{row.worst_progress_residual:.3e} | {row.global_energy_balance_relative_residual:.3e} |"
        for row in result.common_cutoff_runs
    )
    summary = f"""# Continuum helical edge measure on actual Galerkin Navier--Stokes

A single real divergence-free smooth Fourier polynomial is evolved by the
unforced three-dimensional incompressible Navier--Stokes Galerkin system with
Leray projection, viscosity, 2/3 dealiasing and RK4.  For child mode
`{result.child_mode}`, every retained ordered convolution parent is read from the
same evolved state, quotiented into unordered parent orbits, and then resolved
into all eight helicity edges.

The torus Fourier-series convolution has coefficient one.  To test the theorem's
unitary-R3 measure normalization without changing the physical Galerkin source,
each discrete unordered orbit is assigned quotient mass `1/C_F`, so the theorem
factor `C_F` cancels exactly against the discrete counting-measure embedding.

| n | cutoff | unordered pairs | helical edges | source residual | work residual | progress residual | NS energy balance |
|---:|---:|---:|---:|---:|---:|---:|---:|
{rows}

The native cutoffs above define different Galerkin PDEs. Their final-child-energy spread `{result.native_final_child_energy_resolution_spread:.3e}` and integrated-child-work spread `{result.native_integrated_child_work_resolution_spread:.3e}` are **diagnostics only**; no unproved convergence threshold is imposed.

The following table embeds the **same** Galerkin cutoff `{result.common_cutoff}` on every FFT grid:

| n | cutoff | unordered pairs | helical edges | source residual | work residual | progress residual | NS energy balance |
|---:|---:|---:|---:|---:|---:|---:|---:|
{common_rows}

Same-system final-child-energy representation spread: `{result.common_final_child_energy_resolution_spread:.3e}`.
Same-system integrated-child-work representation spread: `{result.common_integrated_child_work_resolution_spread:.3e}`.

This is direct falsification evidence on finite Fourier--Galerkin Navier--Stokes,
not a continuum PDE proof.  It checks that the proposed signed edge measure is
actually the same nonlinear convolution/work law on evolved PDE states before
Hahn splitting.
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
