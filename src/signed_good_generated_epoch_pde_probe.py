from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.full_natural_service_corridor_pde_probe import (
    _leray_dealias,
    _nonlinear_term,
    _rk4_step,
    _spectral_geometry,
)
from src.high_strain_descending_epoch_pde_probe import (
    _divergence_norm,
    _gradient_energy,
    _physical_inner,
    _strict_low_pass_symbol,
)
from src.physical_energy_causal_bridge import (
    heavy_half_physical_transfer,
    route_physical_energy_causality,
)
from src.signed_good_generated_epoch_time_telescope import (
    ACTUAL_HH_GENERATION_BRANCH,
    SignedGoodGeneratedWorkProvenance,
    signed_good_generated_epoch_telescope,
    signed_good_step_from_energy_reentry,
)


STATUS = (
    "DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "EXACT_P_PLUS_Q_SIGNED_GOOD_TRIAD__"
    "ACTUAL_POSITIVE_CHILD_ENERGY_HH_WORK__"
    "TYPED_PHYSICAL_PROVENANCE_TO_T0"
)

PARENT_P = (3, 2, 0)
PARENT_Q = (3, -2, 0)
CHILD_K = (6, 0, 0)


def _mode_mask(k: np.ndarray, wavevectors: Sequence[tuple[int, int, int]]) -> np.ndarray:
    mask = np.zeros(k.shape[1:], dtype=bool)
    for wavevector in wavevectors:
        row = np.ones(k.shape[1:], dtype=bool)
        for direction, value in enumerate(wavevector):
            row &= k[direction] == float(value)
        mask |= row
    return mask


def _signed_good_triad_initial_state(
    resolution: int,
    amplitude: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Real, divergence-free parent pair whose exact sum generates ``k=(6,0,0)``."""
    n = int(resolution)
    state = np.zeros((3, n, n, n), dtype=complex)
    p_polarization = np.asarray((2.0, -3.0, 0.0)) / math.sqrt(13.0)
    q_polarization = np.asarray((0.0, 0.0, 1.0))
    for wavevector, polarization in (
        (PARENT_P, p_polarization),
        (PARENT_Q, q_polarization),
    ):
        positive = tuple(value % n for value in wavevector)
        negative = tuple((-value) % n for value in wavevector)
        state[(slice(None),) + positive] = polarization
        state[(slice(None),) + negative] = polarization

    state = _leray_dealias(state, k, k2, dealias)
    energy = _physical_inner(state, state, n)
    if not math.isfinite(energy) or energy <= 0.0:
        raise AssertionError("signed-good parent triad lost positive physical energy")
    return float(amplitude) * state / math.sqrt(energy)


def _fourier_strain_linf_upper(
    state_hat: np.ndarray,
    symbol: np.ndarray,
    k: np.ndarray,
) -> float:
    """Rigorous triangle-inequality upper for the represented low-pass strain.

    Unlike a grid maximum, this cannot underestimate the continuum supremum of
    the finite Fourier polynomial used by the Galerkin energy gate.
    """
    n = int(state_hat.shape[1])
    resolved = symbol[None, ...] * state_hat
    derivative_bounds = np.empty((3, 3), dtype=float)
    normalization = float(n**3)
    for component in range(3):
        for direction in range(3):
            derivative_bounds[component, direction] = float(
                np.sum(np.abs(k[direction] * resolved[component])) / normalization
            )
    symmetric_bounds = 0.5 * (derivative_bounds + derivative_bounds.T)
    return float(np.linalg.norm(symmetric_bounds))


def _trapezoid(values: Sequence[float], dt: float) -> float:
    rows = tuple(float(value) for value in values)
    if len(rows) < 2 or dt <= 0.0:
        raise ValueError("at least two PDE samples and a positive time step required")
    return dt * (0.5 * rows[0] + math.fsum(rows[1:-1]) + 0.5 * rows[-1])


def _trapezoid_atom_weights(values: Sequence[float], dt: float) -> tuple[float, ...]:
    rows = [float(value) * dt for value in values]
    if len(rows) < 2 or dt <= 0.0:
        raise ValueError("at least two PDE samples and a positive time step required")
    rows[0] *= 0.5
    rows[-1] *= 0.5
    return tuple(rows)


def _observables(
    state_hat: np.ndarray,
    viscosity: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
    parent_mask: np.ndarray,
    child_mask: np.ndarray,
    low_pass_symbol: np.ndarray,
) -> dict[str, float]:
    n = int(state_hat.shape[1])
    nonlinear = _nonlinear_term(state_hat, k, k2, dealias)
    child = child_mask[None, ...] * state_hat
    parent_pair = parent_mask[None, ...] * state_hat

    # A single divergence-free plane wave has zero self-advection.  Therefore
    # N(P_p+P_q), restricted to +/-k, is exactly the two ordered p-q atoms.
    parent_pair_nonlinear = _nonlinear_term(parent_pair, k, k2, dealias)
    designated_hh_forcing = -child_mask[None, ...] * parent_pair_nonlinear
    full_child_nonlinear_forcing = -child_mask[None, ...] * nonlinear
    residual_forcing = full_child_nonlinear_forcing - designated_hh_forcing

    child_energy = _physical_inner(child, child, n)
    child_gradient = _gradient_energy(child, np.ones_like(k2), k2, n)
    designated_signed_work = 2.0 * _physical_inner(child, designated_hh_forcing, n)
    residual_signed_work = 2.0 * _physical_inner(child, residual_forcing, n)
    full_signed_work = 2.0 * _physical_inner(child, full_child_nonlinear_forcing, n)
    decomposition_scale = max(
        abs(full_signed_work),
        abs(designated_signed_work),
        abs(residual_signed_work),
    )
    decomposition_residual = abs(
        full_signed_work - designated_signed_work - residual_signed_work
    )
    if decomposition_scale > 0.0:
        decomposition_residual /= decomposition_scale

    global_energy = _physical_inner(state_hat, state_hat, n)
    global_gradient = _gradient_energy(state_hat, np.ones_like(k2), k2, n)
    global_nonlinear_work = -2.0 * _physical_inner(state_hat, nonlinear, n)
    return {
        "child_energy": child_energy,
        "child_gradient": child_gradient,
        "designated_signed_work": designated_signed_work,
        "designated_positive_work": max(0.0, designated_signed_work),
        "residual_signed_work": residual_signed_work,
        "residual_positive_work": max(0.0, residual_signed_work),
        "full_signed_work": full_signed_work,
        "forcing_decomposition_relative_residual": decomposition_residual,
        "designated_forcing_energy": _physical_inner(
            designated_hh_forcing, designated_hh_forcing, n
        ),
        "global_energy": global_energy,
        "global_gradient": global_gradient,
        "global_nonlinear_work": global_nonlinear_work,
        "divergence_norm": _divergence_norm(state_hat, k, n),
        "low_pass_strain_linf_upper": _fourier_strain_linf_upper(
            state_hat, low_pass_symbol, k
        ),
    }


@dataclass(frozen=True)
class GalerkinSignedGoodRun:
    resolution: int
    steps: int
    dt: float
    duration: float
    viscosity: float
    amplitude: float
    child_frequency: float
    parent_frequency: float
    signed_good_parent_child_ratio: float
    initial_global_energy: float
    final_global_energy: float
    initial_child_energy: float
    final_child_energy: float
    final_child_to_initial_global_energy_ratio: float
    actual_positive_hh_work: float
    energy_gate_hh_work_lower: float
    selected_heavy_half_work: float
    selected_heavy_half_fraction: float
    residual_positive_work: float
    residual_positive_work_to_final_child_energy: float
    low_pass_strain_action_upper: float
    child_energy_balance_relative_residual: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float
    maximum_forcing_decomposition_relative_residual: float
    minimum_designated_forcing_energy: float
    selected_support_start: float
    selected_support_end: float
    normalized_selected_parent_span: float
    route_branch: str
    epoch_hits_initial_boundary: bool


def simulate_signed_good_triad_galerkin(
    *,
    resolution: int,
    steps: int,
    viscosity: float = 0.02,
    amplitude: float = 64.0,
    scaled_lifetime: float = 0.05,
) -> GalerkinSignedGoodRun:
    """Measure one signed-good generated-HH step on an actual Galerkin NS orbit."""
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    c = float(scaled_lifetime)
    if n < 24 or n % 2 or count < 16:
        raise ValueError("an even resolution at least 24 and at least sixteen RK4 steps are required")
    if min(nu, amp, c) <= 0.0 or not all(math.isfinite(x) for x in (nu, amp, c)):
        raise ValueError("positive finite physical Galerkin parameters required")

    child_frequency = math.sqrt(sum(value * value for value in CHILD_K))
    parent_frequency = math.sqrt(sum(value * value for value in PARENT_P))
    duration = c / (child_frequency * child_frequency)
    dt = duration / count
    k, k2, dealias = _spectral_geometry(n)
    parent_mask = _mode_mask(
        k,
        (PARENT_P, tuple(-x for x in PARENT_P), PARENT_Q, tuple(-x for x in PARENT_Q)),
    )
    child_mask = _mode_mask(k, (CHILD_K, tuple(-x for x in CHILD_K)))
    if int(np.sum(parent_mask)) != 4 or int(np.sum(child_mask)) != 2:
        raise AssertionError("Galerkin lattice did not retain the exact signed-good triad")
    low_pass_symbol = _strict_low_pass_symbol(k2, 0.25 * child_frequency)
    state = _signed_good_triad_initial_state(n, amp, k, k2, dealias)

    times = tuple(duration * index / count for index in range(count + 1))
    observations: list[dict[str, float]] = []
    for index in range(count + 1):
        observations.append(
            _observables(
                state,
                nu,
                k,
                k2,
                dealias,
                parent_mask,
                child_mask,
                low_pass_symbol,
            )
        )
        if index < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    initial = observations[0]
    final = observations[-1]
    designated_positive = tuple(row["designated_positive_work"] for row in observations)
    residual_positive = tuple(row["residual_positive_work"] for row in observations)
    actual_positive_hh_work = _trapezoid(designated_positive, dt)
    residual_positive_work = _trapezoid(residual_positive, dt)
    strain_action = _trapezoid(
        tuple(row["low_pass_strain_linf_upper"] for row in observations), dt
    )
    designated_signed = _trapezoid(
        tuple(row["designated_signed_work"] for row in observations), dt
    )
    residual_signed = _trapezoid(
        tuple(row["residual_signed_work"] for row in observations), dt
    )
    child_gradient = _trapezoid(
        tuple(row["child_gradient"] for row in observations), dt
    )
    global_gradient = _trapezoid(
        tuple(row["global_gradient"] for row in observations), dt
    )

    child_delta = final["child_energy"] - initial["child_energy"]
    child_balance_rhs = designated_signed + residual_signed - 2.0 * nu * child_gradient
    child_balance_scale = max(
        abs(child_delta),
        abs(child_balance_rhs),
        final["child_energy"],
    )
    child_balance_residual = abs(child_delta - child_balance_rhs) / child_balance_scale
    global_delta = final["global_energy"] - initial["global_energy"]
    global_balance_rhs = -2.0 * nu * global_gradient
    global_balance_residual = abs(global_delta - global_balance_rhs) / initial["global_energy"]
    nonlinear_rate_scale = initial["global_energy"] / duration
    maximum_global_nonlinear = max(
        abs(row["global_nonlinear_work"]) for row in observations
    ) / nonlinear_rate_scale
    maximum_divergence = max(row["divergence_norm"] for row in observations) / math.sqrt(
        initial["global_energy"]
    )
    maximum_decomposition_residual = max(
        row["forcing_decomposition_relative_residual"] for row in observations
    )

    gate = dict(
        route_physical_energy_causality(
            terminal_energy=final["child_energy"],
            initial_energy=initial["child_energy"],
            residual_positive_work=residual_positive_work,
            strain_action=strain_action,
        )
    )
    if gate.get("branch") != ACTUAL_HH_GENERATION_BRANCH:
        raise AssertionError("actual signed-good Galerkin triad did not enter the HH generation route")

    half = dict(
        heavy_half_physical_transfer(
            times=times,
            positive_work_weights=_trapezoid_atom_weights(designated_positive, dt),
            slab_start=0.0,
            slab_end=duration,
        )
    )
    provenance = SignedGoodGeneratedWorkProvenance(
        event_id=f"galerkin-signed-good-event-n{n}",
        trajectory_id=f"unforced-dealiased-NS-triad-n{n}",
        child_carrier_id="fourier-carrier-(+/-6,0,0)",
        generated_parent_carrier_id="fourier-carrier-|p|=|q|=sqrt(13)",
        work_law_id=f"actual-positive-pq-child-energy-work-n{n}",
        child_frequency=child_frequency,
        parent_frequency=parent_frequency,
        scaled_lifetime=c,
        slab_start=0.0,
        slab_end=duration,
    )
    gate["provenance"] = provenance
    half["provenance"] = provenance
    step = signed_good_step_from_energy_reentry(
        reentry=gate,
        selected_physical_half_slab=half,
        child_frequency=child_frequency,
        parent_frequency=parent_frequency,
        scaled_lifetime=c,
    )
    epoch = signed_good_generated_epoch_telescope((step,))

    work_lower = float(gate["physical_hh_work_lower"])
    work_scale = max(actual_positive_hh_work, work_lower)
    if actual_positive_hh_work + 2.0e-6 * work_scale < work_lower:
        raise AssertionError("measured actual HH work fell below its physical energy-gate lower")
    if child_balance_residual > 2.0e-3:
        raise AssertionError("selected child lost its measured NS energy balance")
    if global_balance_residual > 2.0e-5:
        raise AssertionError("Galerkin orbit lost its global viscous energy balance")
    if maximum_global_nonlinear > 2.0e-11:
        raise AssertionError("dealiased incompressible nonlinearity lost global energy skewness")
    if maximum_divergence > 2.0e-11:
        raise AssertionError("Galerkin orbit lost incompressibility")
    if maximum_decomposition_residual > 2.0e-12:
        raise AssertionError("designated p-q work and residual no longer decompose child work")
    if strain_action > 1.0e-10:
        raise AssertionError("triad fixture acquired a spurious low-pass strain route")
    if final["child_energy"] / initial["global_energy"] <= 1.0e-10:
        raise AssertionError("signed-good PDE fixture generated no resolved child energy")
    if min(row["designated_forcing_energy"] for row in observations) <= 0.0:
        raise AssertionError("exact p-q parent pair lost its physical child forcing")
    if not epoch.hits_initial_boundary:
        raise AssertionError("one natural slab did not register the generated carrier at t=0")

    return GalerkinSignedGoodRun(
        resolution=n,
        steps=count,
        dt=dt,
        duration=duration,
        viscosity=nu,
        amplitude=amp,
        child_frequency=child_frequency,
        parent_frequency=parent_frequency,
        signed_good_parent_child_ratio=parent_frequency / child_frequency,
        initial_global_energy=initial["global_energy"],
        final_global_energy=final["global_energy"],
        initial_child_energy=initial["child_energy"],
        final_child_energy=final["child_energy"],
        final_child_to_initial_global_energy_ratio=final["child_energy"]
        / initial["global_energy"],
        actual_positive_hh_work=actual_positive_hh_work,
        energy_gate_hh_work_lower=work_lower,
        selected_heavy_half_work=float(half["mass"]),
        selected_heavy_half_fraction=float(half["mass"]) / float(half["total"]),
        residual_positive_work=residual_positive_work,
        residual_positive_work_to_final_child_energy=residual_positive_work
        / final["child_energy"],
        low_pass_strain_action_upper=strain_action,
        child_energy_balance_relative_residual=child_balance_residual,
        global_energy_balance_relative_residual=global_balance_residual,
        maximum_global_nonlinear_work_relative_rate=maximum_global_nonlinear,
        maximum_divergence_relative_to_initial_l2=maximum_divergence,
        maximum_forcing_decomposition_relative_residual=maximum_decomposition_residual,
        minimum_designated_forcing_energy=min(
            row["designated_forcing_energy"] for row in observations
        ),
        selected_support_start=float(half["start"]),
        selected_support_end=float(half["end"]),
        normalized_selected_parent_span=step.normalized_parent_span,
        route_branch=str(gate["branch"]),
        epoch_hits_initial_boundary=epoch.hits_initial_boundary,
    )


@dataclass(frozen=True)
class SignedGoodPDEProbe:
    status: str
    runs: tuple[GalerkinSignedGoodRun, ...]
    final_child_energy_resolution_spread: float
    actual_hh_work_resolution_spread: float


def run_probe(
    resolutions: Sequence[int] = (24, 28, 32),
    *,
    steps: int = 128,
    viscosity: float = 0.02,
    amplitude: float = 64.0,
    scaled_lifetime: float = 0.05,
) -> SignedGoodPDEProbe:
    runs = tuple(
        simulate_signed_good_triad_galerkin(
            resolution=int(resolution),
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            scaled_lifetime=float(scaled_lifetime),
        )
        for resolution in resolutions
    )
    if not runs:
        raise ValueError("at least one Galerkin resolution is required")

    child_energies = tuple(row.final_child_energy for row in runs)
    work = tuple(row.actual_positive_hh_work for row in runs)
    child_spread = (max(child_energies) - min(child_energies)) / max(child_energies)
    work_spread = (max(work) - min(work)) / max(work)
    if child_spread > 1.0e-2 or work_spread > 1.0e-2:
        raise AssertionError("signed-good physical observables did not stabilize under resolution refinement")
    return SignedGoodPDEProbe(
        status=STATUS,
        runs=runs,
        final_child_energy_resolution_spread=child_spread,
        actual_hh_work_resolution_spread=work_spread,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", type=int, nargs="+", default=(24, 28, 32))
    ap.add_argument("--steps", type=int, default=128)
    ap.add_argument("--viscosity", type=float, default=0.02)
    ap.add_argument("--amplitude", type=float, default=64.0)
    ap.add_argument("--scaled-lifetime", type=float, default=0.05)
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("results-signed-good-generated-epoch-pde-probe"),
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
    (args.outdir / "signed_good_generated_epoch_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    table = "\n".join(
        "| {n} | {steps} | {child:.3e} | {work:.3e} | {lower:.3e} | "
        "{half:.6f} | {residual:.3e} | {K:.3e} | {cbal:.3e} | {gbal:.3e} |".format(
            n=row.resolution,
            steps=row.steps,
            child=row.final_child_energy,
            work=row.actual_positive_hh_work,
            lower=row.energy_gate_hh_work_lower,
            half=row.selected_heavy_half_fraction,
            residual=row.residual_positive_work_to_final_child_energy,
            K=row.low_pass_strain_action_upper,
            cbal=row.child_energy_balance_relative_residual,
            gbal=row.global_energy_balance_relative_residual,
        )
        for row in result.runs
    )
    md = f"""# Signed-good generated epoch: physical PDE probe

Status: **{STATUS}**.

This evolves the unforced 3D incompressible Fourier--Galerkin Navier--Stokes
system with Leray projection, viscosity, 2/3 dealiasing and RK4.  The real
divergence-free modes `p=(3,2,0)` and `q=(3,-2,0)` generate `k=(6,0,0)` through
the actual quadratic PDE term.  Thus `|p|/|k|=sqrt(13)/6` lies strictly in
`(3/5,5/8)`.  Every reported causal weight is the measured physical law
`2[Re <u_k,F_(p,q)>]_+ dt`; no coefficient impulse or recurrence proxy is used.

| resolution | steps | final child E | actual HH work | gate lower | heavy-half fraction | residual/E_child | low-strain K upper | child balance | global balance |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{table}

Child-energy resolution spread: `{result.final_child_energy_resolution_spread:.3e}`.
Actual-HH-work resolution spread: `{result.actual_hh_work_resolution_spread:.3e}`.

The finite Fourier triangle bound supplies an upper, not a sampled lower, for
the low-pass strain.  The resulting typed work provenance passes through the
physical energy gate, heavy-half selection and one-layer signed-good telescope
on the same orbit.  This is strong numerical falsification evidence for that
executable bridge; it is not a continuum proof or a generic-HH claim.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
