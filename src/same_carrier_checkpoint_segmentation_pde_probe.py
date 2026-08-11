from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.full_natural_service_corridor_pde_probe import (
    _analytic_divergence_free_initial_state,
    _carrier_symbol,
    _leray_dealias,
    _nonlinear_term,
    _observables,
    _rk4_step,
    _spectral_geometry,
    _spectral_inner,
)
from src.same_carrier_checkpoint_segmentation_quotient import (
    SameCarrierProvenance,
    fixed_carrier_natural_window_capacity,
    segmentation_invariance,
)


STATUS = (
    "DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "ONE_EVOLVED_TRAJECTORY_FIXED_Q_AND_TERMINAL_DUAL__"
    "ACTUAL_COMPLEX_HH_RESIDUAL_IMPULSES__"
    "CHECKPOINT_PARTITION_INVARIANCE"
)


def _complex_inner(left: np.ndarray, right: np.ndarray, resolution: int) -> complex:
    return complex(np.vdot(left, right) / float(resolution**6))


def _bilinear_convection(
    advecting_hat: np.ndarray,
    advected_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Leray-projected dealiased B(a,b)=(a.grad)b on the actual grid."""
    velocity = np.fft.ifftn(advecting_hat, axes=(1, 2, 3)).real
    convection = np.zeros_like(velocity)
    for component in range(3):
        for direction in range(3):
            derivative = np.fft.ifftn(
                1j * k[direction] * advected_hat[component],
                axes=(0, 1, 2),
            ).real
            convection[component] += velocity[direction] * derivative
    return _leray_dealias(
        np.fft.fftn(convection, axes=(1, 2, 3)),
        k,
        k2,
        dealias,
    )


def _resolved_strain_rate(
    state_hat: np.ndarray,
    carrier_frequency: float,
    k: np.ndarray,
    k2: np.ndarray,
) -> float:
    """Spatial sup of the actual resolved symmetric-gradient operator norm."""
    n = int(state_hat.shape[1])
    resolved = (np.sqrt(k2) <= 0.25 * carrier_frequency)[None, ...] * state_hat
    gradient = np.empty((n, n, n, 3, 3), dtype=float)
    for component in range(3):
        for direction in range(3):
            gradient[..., component, direction] = np.fft.ifftn(
                1j * k[direction] * resolved[component],
                axes=(0, 1, 2),
            ).real
    symmetric = 0.5 * (gradient + np.swapaxes(gradient, -1, -2))
    eigenvalues = np.linalg.eigvalsh(symmetric)
    return float(np.max(np.abs(eigenvalues)))


def _coefficient_parts(
    state_hat: np.ndarray,
    dual_hat: np.ndarray,
    carrier_symbol: np.ndarray,
    carrier_frequency: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> dict[str, float | complex]:
    """Read HH, resolved-interface, and low-low pieces from the same NS RHS."""
    n = int(state_hat.shape[1])
    resolved_mask = np.sqrt(k2) <= 0.25 * carrier_frequency
    resolved = resolved_mask[None, ...] * state_hat
    high = state_hat - resolved
    low_low = _bilinear_convection(resolved, resolved, k, k2, dealias)
    low_high = _bilinear_convection(resolved, high, k, k2, dealias)
    high_low = _bilinear_convection(high, resolved, k, k2, dealias)
    high_high = _bilinear_convection(high, high, k, k2, dealias)
    reconstructed = low_low + low_high + high_low + high_high
    actual = _nonlinear_term(state_hat, k, k2, dealias)
    difference = actual - reconstructed
    actual_norm = math.sqrt(max(0.0, _spectral_inner(actual, actual, n)))
    difference_norm = math.sqrt(max(0.0, _spectral_inner(difference, difference, n)))
    decomposition_relative = difference_norm / max(actual_norm, 1.0e-300)

    q = carrier_symbol[None, ...]
    forcing_low_low = -q * low_low
    forcing_residual = -q * (low_high + high_low)
    forcing_hh = -q * high_high
    forcing_actual = -q * actual
    coefficient_low_low = _complex_inner(dual_hat, forcing_low_low, n)
    coefficient_residual = _complex_inner(dual_hat, forcing_residual, n)
    coefficient_hh = _complex_inner(dual_hat, forcing_hh, n)
    coefficient_actual = _complex_inner(dual_hat, forcing_actual, n)
    coefficient_scale = max(abs(coefficient_actual), 1.0e-300)
    coefficient_reconstruction = abs(
        coefficient_actual
        - coefficient_low_low
        - coefficient_residual
        - coefficient_hh
    ) / coefficient_scale
    return {
        "hh": coefficient_hh,
        "residual": coefficient_residual,
        "low_low": coefficient_low_low,
        "actual": coefficient_actual,
        "decomposition_relative_residual": decomposition_relative,
        "coefficient_reconstruction_relative_residual": coefficient_reconstruction,
        "low_low_to_actual_coefficient_ratio": abs(coefficient_low_low) / coefficient_scale,
    }


def _reverse_cumulative(step_values: Sequence[complex | float]) -> tuple[complex, ...]:
    total = 0.0j
    cumulative: list[complex] = [0.0j]
    for value in reversed(tuple(step_values)):
        total += complex(value)
        cumulative.append(total)
    return tuple(cumulative)


@dataclass(frozen=True)
class GalerkinSameCarrierRun:
    resolution: int
    steps: int
    dt: float
    terminal_carrier_amplitude: float
    first_elapsed: float | None
    joint_first_stops: tuple[str, ...]
    maximum_divergence_norm: float
    maximum_global_nonlinear_work: float
    global_energy_balance_relative_residual: float
    carrier_energy_balance_relative_residual: float
    maximum_direct_q2_identity_relative_residual: float
    maximum_nonlinear_split_relative_residual: float
    maximum_coefficient_split_relative_residual: float
    maximum_low_low_moat_ratio: float
    maximum_complex_duhamel_relative_residual: float
    maximum_segmentation_first_time_residual: float
    maximum_absolute_hh_impulse: float
    maximum_absolute_residual_impulse: float
    maximum_absolute_imaginary_impulse: float
    fixed_natural_windows_before_t0: int
    fixed_window_interior_zeno_possible: bool


def simulate_same_carrier_galerkin(
    *,
    resolution: int,
    steps: int,
    duration: float,
    viscosity: float,
    carrier_frequency: float,
    natural_window_count: int,
) -> GalerkinSameCarrierRun:
    n = int(resolution)
    count = int(steps)
    T = float(duration)
    nu = float(viscosity)
    A = float(carrier_frequency)
    windows = int(natural_window_count)
    if count < 8 or count % windows or windows < 1:
        raise ValueError("at least eight RK4 steps divisible by the positive natural-window count required")
    if min(T, nu, A) <= 0 or not all(math.isfinite(x) for x in (T, nu, A)):
        raise ValueError("positive finite physical PDE parameters required")
    dt = T / count
    k, k2, dealias = _spectral_geometry(n)
    state = _analytic_divergence_free_initial_state(n, k, k2, dealias)
    states = [state]
    for _ in range(count):
        state = _rk4_step(state, dt, nu, k, k2, dealias)
        states.append(state)

    q = _carrier_symbol(k2, A)
    terminal_carrier = q[None, ...] * states[-1]
    carrier_norm = math.sqrt(max(0.0, _spectral_inner(terminal_carrier, terminal_carrier, n)))
    if carrier_norm <= 0:
        raise AssertionError("evolved NS trajectory supplied no terminal Q carrier")
    # A fixed normalized complex Fourier dual.  The even phase keeps the terminal
    # coefficient large while exposing genuine complex impulse phase on the real
    # Navier--Stokes trajectory.
    phase = np.exp(1j * 0.01 * k2)
    terminal_dual = phase[None, ...] * terminal_carrier / carrier_norm
    terminal_coefficient = _complex_inner(terminal_dual, terminal_carrier, n)
    if abs(terminal_coefficient) <= 0:
        raise AssertionError("terminal complex dual lost the carrier coefficient")

    observables: list[dict[str, float]] = []
    strain_rates: list[float] = []
    hh_rates: list[complex] = []
    residual_rates: list[complex] = []
    low_low_rates: list[complex] = []
    actual_rates: list[complex] = []
    coefficient_values: list[complex] = []
    maximum_split = 0.0
    maximum_coefficient_split = 0.0
    maximum_low_low = 0.0

    for index, current in enumerate(states):
        physical_time = index * dt
        dual = np.exp(nu * k2 * (physical_time - T))[None, ...] * terminal_dual
        carrier = q[None, ...] * current
        coefficient_values.append(_complex_inner(dual, carrier, n))
        parts = _coefficient_parts(current, dual, q, A, k, k2, dealias)
        hh_rates.append(complex(parts["hh"]))
        residual_rates.append(complex(parts["residual"]))
        low_low_rates.append(complex(parts["low_low"]))
        actual_rates.append(complex(parts["actual"]))
        maximum_split = max(maximum_split, float(parts["decomposition_relative_residual"]))
        maximum_coefficient_split = max(
            maximum_coefficient_split,
            float(parts["coefficient_reconstruction_relative_residual"]),
        )
        maximum_low_low = max(maximum_low_low, float(parts["low_low_to_actual_coefficient_ratio"]))
        observables.append(_observables(current, nu, A, k, k2, dealias))
        strain_rates.append(_resolved_strain_rate(current, A, k, k2))

    def trapezoid_steps(values: Sequence[complex | float]) -> tuple[complex, ...]:
        return tuple(0.5 * dt * (complex(values[j]) + complex(values[j + 1])) for j in range(count))

    cumulative_hh = _reverse_cumulative(trapezoid_steps(hh_rates))
    cumulative_residual = _reverse_cumulative(trapezoid_steps(residual_rates))
    cumulative_low_low = _reverse_cumulative(trapezoid_steps(low_low_rates))
    cumulative_actual = _reverse_cumulative(trapezoid_steps(actual_rates))
    cumulative_strain_complex = _reverse_cumulative(trapezoid_steps(strain_rates))
    cumulative_strain = tuple(float(x.real) for x in cumulative_strain_complex)
    elapsed = tuple(j * dt for j in range(count + 1))
    state_tokens = tuple(f"N{n}-state-{count - j}" for j in range(count + 1))
    scaled_lifetime = (T / windows) * A * A
    provenance = SameCarrierProvenance(
        event_id=f"NS-Galerkin-event-N{n}",
        carrier_id=f"fixed-Q-{A:g}-N{n}",
        terminal_dual_id=f"complex-terminal-dual-N{n}",
        trajectory_id=f"unforced-dealiased-NS-N{n}",
        terminal_state_token=state_tokens[0],
        terminal_time=T,
        carrier_frequency=A,
        scaled_lifetime=scaled_lifetime,
        terminal_coefficient=terminal_coefficient,
    )
    cuts = tuple(sorted({count // 4, count // 2, 3 * count // 4}))
    invariant = segmentation_invariance(
        provenance=provenance,
        state_tokens=state_tokens,
        elapsed_times=elapsed,
        strain_action=cumulative_strain,
        residual_impulse=cumulative_residual,
        hh_impulse=cumulative_hh,
        checkpoint_indices=cuts,
        tie_tolerance=0.0,
    )

    maximum_duhamel = 0.0
    for backward_index in range(count + 1):
        forward_index = count - backward_index
        coefficient_delta = terminal_coefficient - coefficient_values[forward_index]
        reconstructed = (
            cumulative_hh[backward_index]
            + cumulative_residual[backward_index]
            + cumulative_low_low[backward_index]
        )
        # The independently accumulated full nonlinear rate is checked as well;
        # both use the same evolved states, not a proxy ODE.
        scale = max(
            abs(terminal_coefficient),
            abs(coefficient_delta),
            abs(reconstructed),
            abs(cumulative_actual[backward_index]),
            1.0e-300,
        )
        maximum_duhamel = max(
            maximum_duhamel,
            abs(coefficient_delta - reconstructed) / scale,
            abs(coefficient_delta - cumulative_actual[backward_index]) / scale,
        )

    carrier_action = sum(
        0.5 * dt * (observables[j]["carrier_power"] + observables[j + 1]["carrier_power"])
        for j in range(count)
    )
    carrier_delta = observables[-1]["carrier_energy"] - observables[0]["carrier_energy"]
    carrier_balance = abs(carrier_delta - carrier_action) / max(
        abs(carrier_delta), abs(carrier_action), observables[0]["carrier_energy"], 1.0e-300
    )
    global_action = sum(
        0.5 * dt * (observables[j]["global_power"] + observables[j + 1]["global_power"])
        for j in range(count)
    )
    global_delta = observables[-1]["global_energy"] - observables[0]["global_energy"]
    global_balance = abs(global_delta - global_action) / max(
        abs(global_delta), abs(global_action), observables[0]["global_energy"], 1.0e-300
    )
    capacity = fixed_carrier_natural_window_capacity(
        event_time=T,
        carrier_frequency=A,
        scaled_lifetime=scaled_lifetime,
    )
    all_impulses = cumulative_hh + cumulative_residual
    return GalerkinSameCarrierRun(
        resolution=n,
        steps=count,
        dt=dt,
        terminal_carrier_amplitude=abs(terminal_coefficient),
        first_elapsed=None if invariant["first_elapsed"] is None else float(invariant["first_elapsed"]),
        joint_first_stops=tuple(invariant["joint_first_stops"]),
        maximum_divergence_norm=max(float(x["divergence_norm"]) for x in observables),
        maximum_global_nonlinear_work=max(abs(float(x["global_nonlinear_work"])) for x in observables),
        global_energy_balance_relative_residual=global_balance,
        carrier_energy_balance_relative_residual=carrier_balance,
        maximum_direct_q2_identity_relative_residual=max(
            float(x["direct_identity_relative_residual"]) for x in observables
        ),
        maximum_nonlinear_split_relative_residual=maximum_split,
        maximum_coefficient_split_relative_residual=maximum_coefficient_split,
        maximum_low_low_moat_ratio=maximum_low_low,
        maximum_complex_duhamel_relative_residual=maximum_duhamel,
        maximum_segmentation_first_time_residual=float(invariant["first_time_residual"]),
        maximum_absolute_hh_impulse=max(abs(x) for x in cumulative_hh),
        maximum_absolute_residual_impulse=max(abs(x) for x in cumulative_residual),
        maximum_absolute_imaginary_impulse=max(abs(x.imag) for x in all_impulses),
        fixed_natural_windows_before_t0=int(capacity["maximum_complete_windows_before_t0"]),
        fixed_window_interior_zeno_possible=bool(capacity["interior_zeno_possible"]),
    )


@dataclass(frozen=True)
class SameCarrierPhysicalPDEProbe:
    status: str
    duration: float
    viscosity: float
    carrier_frequency: float
    natural_window_count: int
    runs: tuple[GalerkinSameCarrierRun, ...]
    terminal_amplitude_resolution_spread: float


def run_same_carrier_physical_pde_probe(
    *,
    resolutions: Sequence[int] = (20, 24, 28),
    steps: int = 80,
    duration: float = 0.015625,
    viscosity: float = 0.05,
    carrier_frequency: float = 4.0,
    natural_window_count: int = 4,
) -> SameCarrierPhysicalPDEProbe:
    runs = tuple(
        simulate_same_carrier_galerkin(
            resolution=int(resolution),
            steps=int(steps),
            duration=float(duration),
            viscosity=float(viscosity),
            carrier_frequency=float(carrier_frequency),
            natural_window_count=int(natural_window_count),
        )
        for resolution in resolutions
    )
    if not runs:
        raise ValueError("at least one Fourier-Galerkin resolution required")
    for run in runs:
        if run.maximum_divergence_norm > 2.0e-11:
            raise AssertionError("same-carrier NS trajectory lost incompressibility")
        if run.maximum_global_nonlinear_work > 2.0e-10:
            raise AssertionError("dealiased NS nonlinearity lost global energy conservation")
        if run.global_energy_balance_relative_residual > 2.0e-5:
            raise AssertionError("global NS energy balance failed on the segmented trajectory")
        if run.carrier_energy_balance_relative_residual > 2.0e-5:
            raise AssertionError("Q^2 carrier energy balance failed on the segmented trajectory")
        if run.maximum_direct_q2_identity_relative_residual > 2.0e-11:
            raise AssertionError("direct Q^2 identity failed on the actual NS state")
        if run.maximum_nonlinear_split_relative_residual > 2.0e-11:
            raise AssertionError("resolved/HH split failed to reconstruct the actual NS nonlinearity")
        if run.maximum_coefficient_split_relative_residual > 2.0e-11:
            raise AssertionError("complex coefficient pieces failed to reconstruct the actual nonlinear rate")
        if run.maximum_low_low_moat_ratio > 2.0e-10:
            raise AssertionError("resolved low-low output crossed the fixed Q support moat")
        if run.maximum_complex_duhamel_relative_residual > 2.0e-4:
            raise AssertionError("actual complex coefficient/Duhamel budget failed")
        if run.maximum_segmentation_first_time_residual != 0.0:
            raise AssertionError("checkpoint cuts changed the actual NS first-stop path")
        if run.maximum_absolute_hh_impulse <= 1.0e-10:
            raise AssertionError("same-carrier PDE probe degenerated to zero HH interaction")
        if run.maximum_absolute_residual_impulse <= 1.0e-12:
            raise AssertionError("same-carrier PDE probe degenerated to zero resolved-interface interaction")
        if run.maximum_absolute_imaginary_impulse <= 1.0e-12:
            raise AssertionError("complex terminal dual lost all physical impulse phase")
        if run.fixed_natural_windows_before_t0 != natural_window_count:
            raise AssertionError("fixed positive natural-window count lost its native time")
        if run.fixed_window_interior_zeno_possible:
            raise AssertionError("fixed-carrier natural windows formed an interior Zeno sequence")
    amplitudes = [run.terminal_carrier_amplitude for run in runs]
    spread = (max(amplitudes) - min(amplitudes)) / max(max(amplitudes), 1.0e-300)
    if spread > 2.0e-2:
        raise AssertionError("terminal same-carrier coefficient did not stabilize under resolution refinement")
    return SameCarrierPhysicalPDEProbe(
        status=STATUS,
        duration=float(duration),
        viscosity=float(viscosity),
        carrier_frequency=float(carrier_frequency),
        natural_window_count=int(natural_window_count),
        runs=runs,
        terminal_amplitude_resolution_spread=spread,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", type=int, nargs="+", default=[20, 24, 28])
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--duration", type=float, default=0.015625)
    parser.add_argument("--viscosity", type=float, default=0.05)
    parser.add_argument("--carrier-frequency", type=float, default=4.0)
    parser.add_argument("--natural-window-count", type=int, default=4)
    parser.add_argument("--outdir", type=Path, default=Path("results-same-carrier-checkpoint-segmentation-pde-probe"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = run_same_carrier_physical_pde_probe(
        resolutions=args.resolutions,
        steps=args.steps,
        duration=args.duration,
        viscosity=args.viscosity,
        carrier_frequency=args.carrier_frequency,
        natural_window_count=args.natural_window_count,
    )
    (args.outdir / "same_carrier_checkpoint_segmentation_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    lines = [
        "# Same-carrier checkpoint segmentation: physical PDE probe",
        "",
        f"Status: **{result.status}**.",
        "",
        "The trajectory is the unforced 3D incompressible Navier--Stokes Fourier-Galerkin system with Leray projection, viscosity, 2/3 dealiasing and RK4.  One fixed Q and one terminal complex dual are read on the same evolved trajectory.  This is numerical falsification evidence, not a continuum proof.",
        "",
        f"Physical interval: `T={result.duration:.12g}`, `A={result.carrier_frequency:.12g}`, `nu={result.viscosity:.12g}`, fixed natural windows `{result.natural_window_count}`.",
        "",
        "| N | steps | div | global balance | Q2 balance | nonlinear split | Duhamel | low-low moat | HH impulse | residual impulse | segmentation |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for run in result.runs:
        lines.append(
            f"| {run.resolution} | {run.steps} | {run.maximum_divergence_norm:.3e} | "
            f"{run.global_energy_balance_relative_residual:.3e} | {run.carrier_energy_balance_relative_residual:.3e} | "
            f"{run.maximum_nonlinear_split_relative_residual:.3e} | {run.maximum_complex_duhamel_relative_residual:.3e} | "
            f"{run.maximum_low_low_moat_ratio:.3e} | {run.maximum_absolute_hh_impulse:.3e} | "
            f"{run.maximum_absolute_residual_impulse:.3e} | {run.maximum_segmentation_first_time_residual:.3e} |"
        )
    lines.extend(
        (
            "",
            f"Terminal coefficient resolution spread: `{result.terminal_amplitude_resolution_spread:.3e}`.",
            "",
            "All cumulative complex impulses, resolved strain, carrier/global balances and checkpoint partitions are read from the same evolved PDE states. No proxy evolution, segment reset, diagnostic-scale duration or synthetic resource is introduced.",
        )
    )
    summary = "\n".join(lines) + "\n"
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
