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
    _nonlinear_term,
    _rk4_step,
    _spectral_geometry,
    _spectral_inner,
)
from src.high_strain_descending_epoch_telescope import (
    HighStrainRenewalStep,
    high_strain_epoch_telescope,
    kinetic_energy_gradient_dissipation_upper,
)
from src.high_strain_dissipation_collision import (
    LOW_STRAIN_THRESHOLD,
    clean_high_strain_dissipation_lower,
    normalized_dissipation_lower,
)
from src.high_strain_resolved_ancestor import high_strain_ancestor_mass_threshold
from src.nn_critical_heat_carrier_seed import renewal_scale


STATUS = (
    "DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "ACTUAL_LOW_PASS_STRAIN_AND_NORMALIZED_DISSIPATION__"
    "PHYSICAL_DV_CRITICAL_SHELL_DISINTEGRATION__"
    "DESCENDING_HIGH_STRAIN_EPOCH_TELESCOPE"
)

TORUS_VOLUME = (2.0 * math.pi) ** 3


def _strict_low_pass_symbol(k2: np.ndarray, radius: float) -> np.ndarray:
    """One fixed compact strict low pass supported in ``|k|<radius``."""
    R = float(radius)
    if R <= 0.0 or not math.isfinite(R):
        raise ValueError("positive finite low-pass radius required")
    rho = np.sqrt(k2) / R
    out = np.zeros_like(rho)
    out[rho <= 0.5] = 1.0
    transition = (rho > 0.5) & (rho < 1.0)
    out[transition] = np.cos(math.pi * (rho[transition] - 0.5)) ** 2
    return out


def _physical_inner(left: np.ndarray, right: np.ndarray, resolution: int) -> float:
    """Integral L2 pairing on the physical ``(2 pi)^3`` Galerkin torus."""
    return TORUS_VOLUME * _spectral_inner(left, right, resolution)


def _field_energy(state_hat: np.ndarray, resolution: int) -> float:
    return _physical_inner(state_hat, state_hat, resolution)


def _gradient_energy(
    state_hat: np.ndarray,
    symbol: np.ndarray,
    k2: np.ndarray,
    resolution: int,
) -> float:
    weighted = np.sqrt(k2)[None, ...] * symbol[None, ...] * state_hat
    return _physical_inner(weighted, weighted, resolution)


def _strain_operator_supremum(
    state_hat: np.ndarray,
    symbol: np.ndarray,
    k: np.ndarray,
) -> float:
    """Grid lower reading of ``||sym grad V||_{L-infinity,op}``.

    The evolved state and low-pass multiplier are the ones used in the
    dissipation integral.  Sampling can underestimate the continuum supremum,
    which is safe for falsifying its lower dissipation consequence.
    """
    resolved = symbol[None, ...] * state_hat
    n = int(state_hat.shape[1])
    gradient = np.empty((3, 3, n, n, n), dtype=float)
    for component in range(3):
        for direction in range(3):
            gradient[component, direction] = np.fft.ifftn(
                1j * k[direction] * resolved[component], axes=(0, 1, 2)
            ).real
    strain = 0.5 * (gradient + np.swapaxes(gradient, 0, 1))
    matrices = np.moveaxis(strain, (0, 1), (-2, -1))
    eigenvalues = np.linalg.eigvalsh(matrices)
    return float(np.max(np.abs(eigenvalues)))


def _divergence_norm(
    state_hat: np.ndarray,
    k: np.ndarray,
    resolution: int,
) -> float:
    divergence_hat = 1j * np.sum(k * state_hat, axis=0)
    mean_square = float(
        np.vdot(divergence_hat, divergence_hat).real / float(resolution**6)
    )
    return math.sqrt(max(0.0, TORUS_VOLUME * mean_square))


def _trapezoid(values: Sequence[float], dt: float) -> float:
    rows = tuple(float(x) for x in values)
    if len(rows) < 2 or dt <= 0.0:
        raise ValueError("at least two samples and a positive time step required")
    return dt * (0.5 * rows[0] + math.fsum(rows[1:-1]) + 0.5 * rows[-1])


def _dyadic_shells(child_frequency: float) -> tuple[float, ...]:
    N = float(child_frequency)
    rows: list[float] = []
    for j in range(32):
        upper = 0.25 * N * (2.0 ** (-j))
        if upper < 0.5:
            break
        rows.append(upper)
    if not rows:
        raise ValueError("child frequency exposes no nonzero Galerkin ancestor shell")
    return tuple(rows)


@dataclass(frozen=True)
class GalerkinHighStrainEpochRun:
    resolution: int
    steps: int
    dt: float
    duration: float
    viscosity: float
    amplitude: float
    child_frequency: float
    selected_ancestor_frequency: float
    renewal_frequency: float
    maximum_divergence_relative_to_initial_l2: float
    maximum_global_nonlinear_work_relative_rate: float
    global_energy_balance_relative_residual: float
    initial_energy: float
    final_energy: float
    global_gradient_reservoir: float
    energy_inequality_gradient_upper: float
    root_strain_action: float
    high_strain_action_threshold: float
    root_normalized_resolved_dissipation: float
    high_strain_dissipation_lower: float
    collision_dissipation_lower_from_measured_action: float
    collision_relative_margin: float
    global_reservoir_relative_margin: float
    retained_critical_dissipation: float
    retained_critical_fraction: float
    retained_half_law_relative_margin: float
    selected_shell_critical_dissipation: float
    maximum_selected_shell_critical_mass: float
    critical_shell_mass_threshold: float
    maximum_descendant_resolved_gradient_relative_to_root: float
    descendant_high_strain_excluded_by_spectral_gap: bool
    observed_epoch_steps: int
    certified_epoch_count_upper: int


def simulate_high_strain_epoch_galerkin(
    *,
    resolution: int,
    steps: int,
    viscosity: float = 0.05,
    amplitude: float = 256.0,
    child_frequency: float = 16.0,
    scaled_lifetime: float = 1.0,
) -> GalerkinHighStrainEpochRun:
    """Run one actual Galerkin NS high-strain event and its physical shell route."""
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    N = float(child_frequency)
    c = float(scaled_lifetime)
    if count < 16 or min(nu, amp, N, c) <= 0.0 or not all(
        math.isfinite(x) for x in (nu, amp, N, c)
    ):
        raise ValueError("positive finite PDE data and at least sixteen steps required")

    duration = c / (N * N)
    dt = duration / count
    k, k2, dealias = _spectral_geometry(n)
    state = amp * _analytic_divergence_free_initial_state(n, k, k2, dealias)
    root_symbol = _strict_low_pass_symbol(k2, 0.25 * N)
    shell_frequencies = _dyadic_shells(N)
    selected_M = shell_frequencies[0]
    selected_A = renewal_scale(selected_M)
    descendant_radius = 0.25 * selected_A
    if descendant_radius >= 1.0:
        raise ValueError(
            "the falsification fixture requires a renewed cutoff below the first nonzero torus mode"
        )
    descendant_symbol = _strict_low_pass_symbol(k2, descendant_radius)
    Dstar = clean_high_strain_dissipation_lower(c)
    mustar = high_strain_ancestor_mass_threshold(c)

    energy: list[float] = []
    global_gradient: list[float] = []
    resolved_gradient: list[float] = []
    strain_sup: list[float] = []
    divergence: list[float] = []
    nonlinear_work: list[float] = []
    descendant_gradient: list[float] = []
    critical_density: list[float] = []
    selected_critical_density: list[float] = []
    selected_mu: list[float] = []

    radius = np.sqrt(k2)
    shell_masks = tuple(
        (radius > 0.5 * M) & (radius <= M) for M in shell_frequencies
    )

    for index in range(count + 1):
        energy.append(_field_energy(state, n))
        global_gradient.append(_gradient_energy(state, np.ones_like(k2), k2, n))
        resolved_gradient.append(_gradient_energy(state, root_symbol, k2, n))
        strain_sup.append(_strain_operator_supremum(state, root_symbol, k))
        divergence.append(_divergence_norm(state, k, n))
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        nonlinear_work.append(abs(-2.0 * _physical_inner(state, nonlinear, n)))
        descendant_gradient.append(
            _gradient_energy(state, descendant_symbol, k2, n)
        )

        good_density = 0.0
        selected_density = 0.0
        selected_mass = 0.0
        for M, mask in zip(shell_frequencies, shell_masks):
            shell_state = mask[None, ...] * state
            shell_energy = _field_energy(shell_state, n)
            critical_mass = M * shell_energy
            shell_resolved_gradient = _gradient_energy(
                state, mask.astype(float) * root_symbol, k2, n
            )
            density = N * shell_resolved_gradient
            if critical_mass >= mustar * (1.0 - 2.0e-12):
                good_density += density
                if M == selected_M:
                    selected_density = density
            if M == selected_M:
                selected_mass = critical_mass
        critical_density.append(good_density)
        selected_critical_density.append(selected_density)
        selected_mu.append(selected_mass)

        if index < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    E0 = energy[0]
    E1 = energy[-1]
    G = _trapezoid(global_gradient, dt)
    D = N * _trapezoid(resolved_gradient, dt)
    K = _trapezoid(strain_sup, dt)
    Dgood = _trapezoid(critical_density, dt)
    Dselected = _trapezoid(selected_critical_density, dt)
    descendant_max = max(descendant_gradient)
    collision_lower = normalized_dissipation_lower(K, c, 0.25)
    energy_upper = kinetic_energy_gradient_dissipation_upper(E0, nu)

    energy_balance = abs(E1 - E0 + 2.0 * nu * G) / E0
    divergence_relative = max(divergence) / math.sqrt(E0)
    nonlinear_rate_scale = E0 / duration
    nonlinear_relative = max(nonlinear_work) / nonlinear_rate_scale
    collision_margin = (D - collision_lower) / max(D, collision_lower)
    reservoir_capacity = N * G
    reservoir_margin = (reservoir_capacity - D) / reservoir_capacity
    retained_fraction = Dgood / D
    retained_half_margin = (Dgood - 0.5 * D) / D
    descendant_relative = descendant_max / max(resolved_gradient)

    if divergence_relative > 2.0e-11:
        raise AssertionError("Galerkin trajectory lost incompressibility")
    if nonlinear_relative > 2.0e-11:
        raise AssertionError("Galerkin convection lost global energy skewness")
    if energy_balance > 5.0e-3:
        raise AssertionError("Galerkin trajectory lost the viscous energy balance")
    if K < float(LOW_STRAIN_THRESHOLD):
        raise AssertionError("physical Galerkin fixture did not reach high strain")
    if D < Dstar * (1.0 - 2.0e-10):
        raise AssertionError("measured high-strain event fell below D_*")
    if collision_margin < -5.0e-3:
        raise AssertionError("measured strain/dissipation collision failed")
    if reservoir_margin < -2.0e-10:
        raise AssertionError("resolved dissipation exceeded the actual global gradient reservoir")
    if retained_half_margin < -5.0e-3:
        raise AssertionError("actual critical-shell D_V law lost its retained half")
    if Dselected <= 0.0 or max(selected_mu) < mustar:
        raise AssertionError("selected physical M=N/4 shell carried no critical D_V law")
    if descendant_relative > 1.0e-24:
        raise AssertionError("renewed cutoff below the torus spectral gap retained a false mode")
    if G > energy_upper * (1.0 + 2.0e-10):
        raise AssertionError("finite trajectory exceeded the NS energy-inequality reservoir")

    step = HighStrainRenewalStep(N, selected_M, selected_A, D)
    epoch = high_strain_epoch_telescope(
        (step,), total_gradient_dissipation=G, scaled_lifetime=c
    )

    return GalerkinHighStrainEpochRun(
        resolution=n,
        steps=count,
        dt=dt,
        duration=duration,
        viscosity=nu,
        amplitude=amp,
        child_frequency=N,
        selected_ancestor_frequency=selected_M,
        renewal_frequency=selected_A,
        maximum_divergence_relative_to_initial_l2=divergence_relative,
        maximum_global_nonlinear_work_relative_rate=nonlinear_relative,
        global_energy_balance_relative_residual=energy_balance,
        initial_energy=E0,
        final_energy=E1,
        global_gradient_reservoir=G,
        energy_inequality_gradient_upper=energy_upper,
        root_strain_action=K,
        high_strain_action_threshold=float(LOW_STRAIN_THRESHOLD),
        root_normalized_resolved_dissipation=D,
        high_strain_dissipation_lower=Dstar,
        collision_dissipation_lower_from_measured_action=collision_lower,
        collision_relative_margin=collision_margin,
        global_reservoir_relative_margin=reservoir_margin,
        retained_critical_dissipation=Dgood,
        retained_critical_fraction=retained_fraction,
        retained_half_law_relative_margin=retained_half_margin,
        selected_shell_critical_dissipation=Dselected,
        maximum_selected_shell_critical_mass=max(selected_mu),
        critical_shell_mass_threshold=mustar,
        maximum_descendant_resolved_gradient_relative_to_root=descendant_relative,
        descendant_high_strain_excluded_by_spectral_gap=True,
        observed_epoch_steps=epoch.step_count,
        certified_epoch_count_upper=epoch.certified_count_upper,
    )


def run_probe(
    resolutions: Sequence[int],
    *,
    steps: int,
    viscosity: float,
    amplitude: float,
    child_frequency: float,
    scaled_lifetime: float,
) -> tuple[GalerkinHighStrainEpochRun, ...]:
    rows = tuple(
        simulate_high_strain_epoch_galerkin(
            resolution=int(n),
            steps=steps,
            viscosity=viscosity,
            amplitude=amplitude,
            child_frequency=child_frequency,
            scaled_lifetime=scaled_lifetime,
        )
        for n in resolutions
    )
    if not rows:
        raise ValueError("at least one Galerkin resolution required")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", type=int, nargs="+", default=(12, 16, 20))
    ap.add_argument("--steps", type=int, default=128)
    ap.add_argument("--viscosity", type=float, default=0.05)
    ap.add_argument("--amplitude", type=float, default=256.0)
    ap.add_argument("--child-frequency", type=float, default=16.0)
    ap.add_argument("--scaled-lifetime", type=float, default=1.0)
    ap.add_argument(
        "--outdir", type=Path, default=Path("results-high-strain-descending-epoch-pde-probe")
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = run_probe(
        args.resolutions,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        child_frequency=args.child_frequency,
        scaled_lifetime=args.scaled_lifetime,
    )
    payload = {"status": STATUS, "runs": [asdict(x) for x in rows]}
    (args.outdir / "high_strain_descending_epoch_pde_probe.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    scale = max(x.root_normalized_resolved_dissipation for x in rows)
    resolution_spread = (
        max(x.root_normalized_resolved_dissipation for x in rows)
        - min(x.root_normalized_resolved_dissipation for x in rows)
    ) / scale
    table = "\n".join(
        "| {n} | {steps} | {div:.3e} | {bal:.3e} | {K:.3e} | {D:.3e} | "
        "{coll:.3e} | {reservoir:.3e} | {retained:.6f} | {half:.3e} | "
        "{desc:.3e} | {bound} |".format(
            n=x.resolution,
            steps=x.steps,
            div=x.maximum_divergence_relative_to_initial_l2,
            bal=x.global_energy_balance_relative_residual,
            K=x.root_strain_action,
            D=x.root_normalized_resolved_dissipation,
            coll=x.collision_relative_margin,
            reservoir=x.global_reservoir_relative_margin,
            retained=x.retained_critical_fraction,
            half=x.retained_half_law_relative_margin,
            desc=x.maximum_descendant_resolved_gradient_relative_to_root,
            bound=x.certified_epoch_count_upper,
        )
        for x in rows
    )
    md = f"""# High-strain descending epoch: physical PDE probe

Status: **{STATUS}**.

The probe evolves the unforced three-dimensional incompressible Fourier--Galerkin
Navier--Stokes system with Leray projection, viscosity, `2/3` dealiasing and RK4.
On the same evolved states it reads the strict low pass `S_(N/4)u`, integrates
the actual strain action `K_N`, normalized resolved dissipation `D_N`, and global
gradient reservoir `G_*`, and disintegrates the positive `D_V` density over the
physical dyadic ancestor shells.  No random operator or recurrence proxy supplies
the reported epoch.

The fixture uses `N={args.child_frequency:g}`, `c={args.scaled_lifetime:g}`,
`nu={args.viscosity:g}`, amplitude `{args.amplitude:g}` and the actual natural
duration `cN^-2`.  Its selected critical shell is `M=N/4` and renewal is
`A=3M/4`.  On this periodic falsification fixture the renewed cutoff `A/4<1`
lies below the first nonzero Fourier mode, so the next high-strain step is
physically absent on the evolved mean-zero Galerkin trajectory.

| resolution | steps | div/||u0|| | energy balance | K_N | D_N | collision margin | reservoir margin | retained D_V fraction | half-law margin | descendant/root grad | certified count upper |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

Root-dissipation resolution spread: `{resolution_spread:.3e}`.

This is a numerical falsification test on an actual finite Galerkin NS system,
not a continuum proof and not evidence that every continuum high-strain epoch
terminates after one step.  The exact continuum conclusion still comes from
`D_*<=D_j<=N_jG_*` and the physical `N_(j+1)/N_j<=3/16` renewal law.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
