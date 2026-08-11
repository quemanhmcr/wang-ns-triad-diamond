from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.critical_annular_carrier_service_reentry import (
    BOUNDED_HEAT_RADIUS,
    bounded_heat_defect_fraction_lower,
    gaussian_3d_tail_probability,
    heat_defect_fraction_lower,
    transported_annular_support_ratios,
)


STATUS = (
    "DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "Q2_CARRIER_ENERGY_AND_ENDPOINT_HARD_SHELL_COVER__"
    "SAME_INTERVAL_BOUNDED_INCREMENT_SERVICE"
)


def _spectral_geometry(resolution: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = int(resolution)
    if n < 12 or n % 2:
        raise ValueError("an even Fourier resolution at least 12 is required")
    one = np.fft.fftfreq(n, d=1.0 / n)
    k = np.asarray(np.meshgrid(one, one, one, indexing="ij"), dtype=float)
    k2 = np.sum(k * k, axis=0)
    cutoff = n // 3 - 1
    dealias = np.max(np.abs(k), axis=0) <= cutoff
    return k, k2, dealias


def _leray_dealias(field_hat: np.ndarray, k: np.ndarray, k2: np.ndarray, dealias: np.ndarray) -> np.ndarray:
    out = np.asarray(field_hat, dtype=complex) * dealias[None, ...]
    inverse = np.zeros_like(k2)
    nonzero = k2 > 0
    inverse[nonzero] = 1.0 / k2[nonzero]
    longitudinal = np.sum(k * out, axis=0) * inverse
    out = out - k * longitudinal[None, ...]
    out[:, 0, 0, 0] = 0.0
    return out * dealias[None, ...]


def _spectral_inner(left: np.ndarray, right: np.ndarray, resolution: int) -> float:
    return float(np.vdot(left, right).real / float(resolution**6))


def _analytic_divergence_free_initial_state(
    resolution: int,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """A fixed multiscale field, sampled exactly on every audit resolution."""
    n = int(resolution)
    grid = 2.0 * math.pi * np.arange(n, dtype=float) / n
    x, y, z = np.meshgrid(grid, grid, grid, indexing="ij")

    # A stream-function mode, an ABC mode, and two transverse plane waves.
    # Each summand is analytically divergence free; their interactions make the
    # retained carrier exchange energy nontrivially through the NS nonlinearity.
    u = np.zeros((3, n, n, n), dtype=float)
    u[0] += 2.0 * np.sin(x) * np.cos(2.0 * y) * np.cos(z)
    u[1] += -np.cos(x) * np.sin(2.0 * y) * np.cos(z)

    u[0] += 0.35 * (np.sin(z) + 0.7 * np.cos(y))
    u[1] += 0.35 * (0.8 * np.sin(x) + np.cos(z))
    u[2] += 0.35 * (0.7 * np.sin(y) + 0.8 * np.cos(x))

    phase_one = x + 2.0 * y + 0.37
    u[0] += 0.17 * 2.0 * np.cos(phase_one)
    u[1] += -0.17 * np.cos(phase_one)
    u[2] += 0.17 * 0.5 * np.cos(phase_one)

    phase_two = 2.0 * x - y + z + 0.23
    u[0] += 0.11 * np.sin(phase_two)
    u[1] += 0.11 * np.sin(phase_two)
    u[2] += -0.11 * np.sin(phase_two)

    state = np.fft.fftn(u, axes=(1, 2, 3))
    state = _leray_dealias(state, k, k2, dealias)
    energy = _spectral_inner(state, state, n)
    if energy <= 0 or not math.isfinite(energy):
        raise AssertionError("analytic NS initial state lost positive finite energy")
    return state / math.sqrt(energy)


def _nonlinear_term(
    state_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Leray-projected, 2/3-dealiased Fourier-Galerkin (u.grad)u."""
    velocity = np.fft.ifftn(state_hat, axes=(1, 2, 3)).real
    convection = np.zeros_like(velocity)
    for component in range(3):
        for direction in range(3):
            derivative = np.fft.ifftn(
                1j * k[direction] * state_hat[component],
                axes=(0, 1, 2),
            ).real
            convection[component] += velocity[direction] * derivative
    convection_hat = np.fft.fftn(convection, axes=(1, 2, 3))
    return _leray_dealias(convection_hat, k, k2, dealias)


def _rhs(
    state_hat: np.ndarray,
    viscosity: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    nonlinear = _nonlinear_term(state_hat, k, k2, dealias)
    return _leray_dealias(-nonlinear - viscosity * k2[None, ...] * state_hat, k, k2, dealias)


def _rk4_step(
    state_hat: np.ndarray,
    dt: float,
    viscosity: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    k1 = _rhs(state_hat, viscosity, k, k2, dealias)
    k2_stage = _rhs(state_hat + 0.5 * dt * k1, viscosity, k, k2, dealias)
    k3 = _rhs(state_hat + 0.5 * dt * k2_stage, viscosity, k, k2, dealias)
    k4 = _rhs(state_hat + dt * k3, viscosity, k, k2, dealias)
    advanced = state_hat + (dt / 6.0) * (k1 + 2.0 * k2_stage + 2.0 * k3 + k4)
    return _leray_dealias(advanced, k, k2, dealias)


def _carrier_symbol(k2: np.ndarray, carrier_frequency: float) -> np.ndarray:
    lo, hi = transported_annular_support_ratios()
    rho = np.sqrt(k2) / float(carrier_frequency)
    coordinate = (rho - lo) / (hi - lo)
    inside = (coordinate > 0.0) & (coordinate < 1.0)
    symbol = np.zeros_like(rho)
    symbol[inside] = np.sin(math.pi * coordinate[inside]) ** 2
    return symbol


def _intrinsic_heat_increment_energies(
    carrier_hat: np.ndarray,
    k2: np.ndarray,
    carrier_frequency: float,
    resolution: int,
) -> tuple[float, float]:
    """Exact heat-law service and its radius-3/A retained lower.

    For the intrinsic displacement r~N(0,A^-2 I), Plancherel gives the exact
    multiplier 2(1-exp(-|k|^2/(2A^2))).  Removing |r|>3/A costs at most four
    times the Gaussian tail probability times ||Q_A u||_2^2.  This is the same
    heat-semigroup law and the same truncation used by the continuum theorem;
    no observer-chosen list of displacement directions is substituted.
    """
    A = float(carrier_frequency)
    heat_multiplier = 2.0 * (1.0 - np.exp(-0.5 * k2 / (A * A)))
    full_heat = _spectral_inner(
        np.sqrt(heat_multiplier)[None, ...] * carrier_hat,
        np.sqrt(heat_multiplier)[None, ...] * carrier_hat,
        resolution,
    )
    carrier_energy = _spectral_inner(carrier_hat, carrier_hat, resolution)
    tail = gaussian_3d_tail_probability(BOUNDED_HEAT_RADIUS)
    bounded_heat_lower = full_heat - 4.0 * tail * carrier_energy
    return full_heat, bounded_heat_lower


def _observables(
    state_hat: np.ndarray,
    viscosity: float,
    carrier_frequency: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> dict[str, float]:
    n = int(state_hat.shape[1])
    q = _carrier_symbol(k2, carrier_frequency)
    nonlinear = _nonlinear_term(state_hat, k, k2, dealias)
    derivative = -nonlinear - viscosity * k2[None, ...] * state_hat
    carrier = q[None, ...] * state_hat
    quadratic_effect_state = (q * q)[None, ...] * state_hat

    carrier_energy = _spectral_inner(carrier, carrier, n)
    carrier_critical_mass = carrier_frequency * carrier_energy
    carrier_gradient = _spectral_inner(
        np.sqrt(k2)[None, ...] * carrier,
        np.sqrt(k2)[None, ...] * carrier,
        n,
    )
    nonlinear_work = -2.0 * _spectral_inner(quadratic_effect_state, nonlinear, n)
    carrier_power = nonlinear_work - 2.0 * viscosity * carrier_gradient
    direct_carrier_derivative = 2.0 * _spectral_inner(carrier, q[None, ...] * derivative, n)
    identity_scale = max(abs(carrier_power), abs(direct_carrier_derivative), 1.0e-300)

    radius = np.sqrt(k2)
    shell_a = (radius > 0.5 * carrier_frequency) & (radius <= carrier_frequency)
    shell_2a = (radius > carrier_frequency) & (radius <= 2.0 * carrier_frequency)
    shell_a_energy = _spectral_inner(shell_a[None, ...] * state_hat, shell_a[None, ...] * state_hat, n)
    shell_2a_energy = _spectral_inner(shell_2a[None, ...] * state_hat, shell_2a[None, ...] * state_hat, n)
    mu_a = carrier_frequency * shell_a_energy
    mu_2a = 2.0 * carrier_frequency * shell_2a_energy
    cover_sum = mu_a + 0.5 * mu_2a
    cover_scale = max(carrier_critical_mass, 1.0e-300)

    global_energy = _spectral_inner(state_hat, state_hat, n)
    global_gradient = _spectral_inner(
        np.sqrt(k2)[None, ...] * state_hat,
        np.sqrt(k2)[None, ...] * state_hat,
        n,
    )
    global_nonlinear_work = -2.0 * _spectral_inner(state_hat, nonlinear, n)
    global_power = global_nonlinear_work - 2.0 * viscosity * global_gradient
    divergence_hat = 1j * np.sum(k * state_hat, axis=0)
    divergence_norm = math.sqrt(
        max(0.0, float(np.vdot(divergence_hat, divergence_hat).real / float(n**6)))
    )

    full_heat_energy, bounded_heat_energy_lower = _intrinsic_heat_increment_energies(
        carrier,
        k2,
        carrier_frequency,
        n,
    )
    full_heat_service = carrier_frequency * full_heat_energy
    bounded_heat_service_lower = carrier_frequency * bounded_heat_energy_lower
    return {
        "carrier_energy": carrier_energy,
        "carrier_critical_mass": carrier_critical_mass,
        "carrier_power": carrier_power,
        "carrier_nonlinear_work": nonlinear_work,
        "direct_identity_relative_residual": abs(direct_carrier_derivative - carrier_power) / identity_scale,
        "hard_shell_cover_relative_margin": (cover_sum - carrier_critical_mass) / cover_scale,
        "hard_shell_max_relative_margin": (max(mu_a, mu_2a) - (2.0 / 3.0) * carrier_critical_mass) / cover_scale,
        "full_heat_service": full_heat_service,
        "full_heat_service_to_carrier_mass_ratio": full_heat_service / cover_scale,
        "bounded_heat_service_lower": bounded_heat_service_lower,
        "bounded_heat_lower_to_carrier_mass_ratio": bounded_heat_service_lower / cover_scale,
        "global_energy": global_energy,
        "global_power": global_power,
        "global_nonlinear_work": global_nonlinear_work,
        "divergence_norm": divergence_norm,
    }


@dataclass(frozen=True)
class GalerkinCorridorRun:
    resolution: int
    steps: int
    dt: float
    initial_carrier_energy: float
    final_carrier_energy: float
    carrier_energy_balance_relative_residual: float
    global_energy_balance_relative_residual: float
    maximum_direct_q2_identity_relative_residual: float
    maximum_global_nonlinear_work: float
    maximum_divergence_norm: float
    minimum_hard_shell_cover_relative_margin: float
    minimum_hard_shell_max_relative_margin: float
    minimum_full_heat_service_to_carrier_mass_ratio: float
    minimum_bounded_heat_lower_to_carrier_mass_ratio: float
    integrated_same_interval_bounded_heat_service_lower: float
    maximum_absolute_carrier_nonlinear_work: float


def simulate_galerkin_corridor(
    *,
    resolution: int,
    steps: int,
    duration: float,
    viscosity: float,
    carrier_frequency: float,
) -> GalerkinCorridorRun:
    n = int(resolution)
    count = int(steps)
    T = float(duration)
    nu = float(viscosity)
    A = float(carrier_frequency)
    if count < 2 or min(T, nu, A) <= 0 or not all(math.isfinite(x) for x in (T, nu, A)):
        raise ValueError("positive finite PDE corridor parameters and at least two RK4 steps required")
    dt = T / count
    k, k2, dealias = _spectral_geometry(n)
    state = _analytic_divergence_free_initial_state(n, k, k2, dealias)
    before = _observables(state, nu, A, k, k2, dealias)
    initial = before

    carrier_action = 0.0
    global_action = 0.0
    integrated_service = 0.0
    max_identity = float(before["direct_identity_relative_residual"])
    max_global_nonlinear = abs(float(before["global_nonlinear_work"]))
    max_divergence = float(before["divergence_norm"])
    min_cover = float(before["hard_shell_cover_relative_margin"])
    min_max_cover = float(before["hard_shell_max_relative_margin"])
    min_full_heat = float(before["full_heat_service_to_carrier_mass_ratio"])
    min_bounded_heat = float(before["bounded_heat_lower_to_carrier_mass_ratio"])
    max_carrier_nonlinear = abs(float(before["carrier_nonlinear_work"]))

    for _ in range(count):
        state = _rk4_step(state, dt, nu, k, k2, dealias)
        after = _observables(state, nu, A, k, k2, dealias)
        carrier_action += 0.5 * dt * (float(before["carrier_power"]) + float(after["carrier_power"]))
        global_action += 0.5 * dt * (float(before["global_power"]) + float(after["global_power"]))
        integrated_service += 0.5 * dt * A**2 * (
            float(before["bounded_heat_service_lower"]) + float(after["bounded_heat_service_lower"])
        )
        max_identity = max(max_identity, float(after["direct_identity_relative_residual"]))
        max_global_nonlinear = max(max_global_nonlinear, abs(float(after["global_nonlinear_work"])))
        max_divergence = max(max_divergence, float(after["divergence_norm"]))
        min_cover = min(min_cover, float(after["hard_shell_cover_relative_margin"]))
        min_max_cover = min(min_max_cover, float(after["hard_shell_max_relative_margin"]))
        min_full_heat = min(min_full_heat, float(after["full_heat_service_to_carrier_mass_ratio"]))
        min_bounded_heat = min(
            min_bounded_heat,
            float(after["bounded_heat_lower_to_carrier_mass_ratio"]),
        )
        max_carrier_nonlinear = max(max_carrier_nonlinear, abs(float(after["carrier_nonlinear_work"])))
        before = after

    carrier_delta = float(before["carrier_energy"]) - float(initial["carrier_energy"])
    carrier_scale = max(abs(carrier_delta), abs(carrier_action), float(initial["carrier_energy"]), 1.0e-300)
    global_delta = float(before["global_energy"]) - float(initial["global_energy"])
    global_scale = max(abs(global_delta), abs(global_action), float(initial["global_energy"]), 1.0e-300)
    return GalerkinCorridorRun(
        resolution=n,
        steps=count,
        dt=dt,
        initial_carrier_energy=float(initial["carrier_energy"]),
        final_carrier_energy=float(before["carrier_energy"]),
        carrier_energy_balance_relative_residual=abs(carrier_delta - carrier_action) / carrier_scale,
        global_energy_balance_relative_residual=abs(global_delta - global_action) / global_scale,
        maximum_direct_q2_identity_relative_residual=max_identity,
        maximum_global_nonlinear_work=max_global_nonlinear,
        maximum_divergence_norm=max_divergence,
        minimum_hard_shell_cover_relative_margin=min_cover,
        minimum_hard_shell_max_relative_margin=min_max_cover,
        minimum_full_heat_service_to_carrier_mass_ratio=min_full_heat,
        minimum_bounded_heat_lower_to_carrier_mass_ratio=min_bounded_heat,
        integrated_same_interval_bounded_heat_service_lower=integrated_service,
        maximum_absolute_carrier_nonlinear_work=max_carrier_nonlinear,
    )


@dataclass(frozen=True)
class PhysicalPDEProbe:
    status: str
    duration: float
    scaled_lifetime: float
    viscosity: float
    carrier_frequency: float
    runs: tuple[GalerkinCorridorRun, ...]
    final_carrier_energy_resolution_spread: float


def run_physical_pde_probe(
    *,
    resolutions: Sequence[int] = (12, 16, 20),
    steps: int = 80,
    duration: float = 0.0025,
    viscosity: float = 0.05,
    carrier_frequency: float = 2.0,
) -> PhysicalPDEProbe:
    A = float(carrier_frequency)
    T = float(duration)
    runs = tuple(
        simulate_galerkin_corridor(
            resolution=int(n),
            steps=int(steps),
            duration=T,
            viscosity=float(viscosity),
            carrier_frequency=A,
        )
        for n in resolutions
    )
    if not runs:
        raise ValueError("at least one Galerkin resolution is required")

    for run in runs:
        if run.maximum_direct_q2_identity_relative_residual > 2.0e-11:
            raise AssertionError("direct Q^2 carrier identity failed on the evolved NS state")
        if run.carrier_energy_balance_relative_residual > 2.0e-6:
            raise AssertionError("time-integrated carrier energy balance failed")
        if run.global_energy_balance_relative_residual > 2.0e-6:
            raise AssertionError("time-integrated global NS energy balance failed")
        if run.maximum_global_nonlinear_work > 2.0e-10:
            raise AssertionError("dealiased incompressible nonlinearity failed global energy conservation")
        if run.maximum_divergence_norm > 2.0e-11:
            raise AssertionError("Fourier-Galerkin trajectory lost incompressibility")
        if run.minimum_hard_shell_cover_relative_margin < -2.0e-12:
            raise AssertionError("actual evolved Q carrier escaped its two hard shells")
        if run.minimum_hard_shell_max_relative_margin < -2.0e-12:
            raise AssertionError("actual evolved endpoint shell witness lost the two-thirds lower")
        if run.minimum_full_heat_service_to_carrier_mass_ratio + 2.0e-12 < heat_defect_fraction_lower():
            raise AssertionError("evolved annular carrier violated the intrinsic heat-service lower")
        if (
            run.minimum_bounded_heat_lower_to_carrier_mass_ratio + 2.0e-12
            < bounded_heat_defect_fraction_lower()
        ):
            raise AssertionError("evolved annular carrier violated the radius-3 bounded heat lower")
        if run.integrated_same_interval_bounded_heat_service_lower <= 0:
            raise AssertionError("same physical corridor carried no integrated bounded service")
        if run.maximum_absolute_carrier_nonlinear_work <= 1.0e-8:
            raise AssertionError("PDE probe degenerated to a nonlinear-null carrier")

    finals = [run.final_carrier_energy for run in runs]
    spread = (max(finals) - min(finals)) / max(max(finals), 1.0e-300)
    if spread > 2.0e-2:
        raise AssertionError("carrier observable did not stabilize under Galerkin resolution refinement")
    return PhysicalPDEProbe(
        status=STATUS,
        duration=T,
        scaled_lifetime=T * A * A,
        viscosity=float(viscosity),
        carrier_frequency=A,
        runs=runs,
        final_carrier_energy_resolution_spread=spread,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", type=int, nargs="+", default=[12, 16, 20])
    parser.add_argument("--steps", type=int, default=80)
    parser.add_argument("--duration", type=float, default=0.0025)
    parser.add_argument("--viscosity", type=float, default=0.05)
    parser.add_argument("--carrier-frequency", type=float, default=2.0)
    parser.add_argument("--outdir", type=Path, default=Path("results-full-natural-service-corridor-pde-probe"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = run_physical_pde_probe(
        resolutions=args.resolutions,
        steps=args.steps,
        duration=args.duration,
        viscosity=args.viscosity,
        carrier_frequency=args.carrier_frequency,
    )
    payload = asdict(result)
    (args.outdir / "full_natural_service_corridor_pde_probe.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    lines = [
        "# Full-natural corridor physical PDE probe",
        "",
        f"Status: **{result.status}**.",
        "",
        "The trajectory is the unforced 3D incompressible Navier--Stokes Fourier-Galerkin system with Leray projection, viscosity, 2/3 dealiasing, and RK4 time integration.  It is not a proxy evolution.  The experiment is numerical falsification evidence, not a continuum proof.",
        "",
        f"Physical corridor: `T={result.duration:.12g}=c A^-2`, `c={result.scaled_lifetime:.12g}`, `A={result.carrier_frequency:.12g}`, `nu={result.viscosity:.12g}`.",
        "",
        f"Intrinsic heat-law thresholds: full annular fraction `q={heat_defect_fraction_lower():.12g}`; radius-{BOUNDED_HEAT_RADIUS:g}/A retained fraction `q_b={bounded_heat_defect_fraction_lower():.12g}`.",
        "",
        "| N | steps | Q2 identity | carrier balance | global balance | shell-cover margin | full heat/carrier | bounded heat lower/carrier | nonlinear work |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for run in result.runs:
        lines.append(
            f"| {run.resolution} | {run.steps} | {run.maximum_direct_q2_identity_relative_residual:.3e} | "
            f"{run.carrier_energy_balance_relative_residual:.3e} | {run.global_energy_balance_relative_residual:.3e} | "
            f"{run.minimum_hard_shell_cover_relative_margin:.3e} | "
            f"{run.minimum_full_heat_service_to_carrier_mass_ratio:.3e} | "
            f"{run.minimum_bounded_heat_lower_to_carrier_mass_ratio:.3e} | "
            f"{run.maximum_absolute_carrier_nonlinear_work:.3e} |"
        )
    lines.extend(
        (
            "",
            f"Final carrier-energy resolution spread: `{result.final_carrier_energy_resolution_spread:.3e}`.",
            "",
            "The `Q^2` carrier balance, exact two-hard-shell cover, positive bounded increment service, and physical time integration are all read from the same evolved PDE corridor; none is introduced as an extra event or synthetic resource.",
        )
    )
    summary = "\n".join(lines) + "\n"
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
