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
    _rk4_step,
    _spectral_geometry,
)
from src.high_strain_descending_epoch_pde_probe import _strict_low_pass_symbol
from src.physical_energy_causal_bridge import (
    heavy_half_physical_transfer,
    route_physical_energy_causality,
)
from src.signed_good_generated_epoch_pde_probe import (
    CHILD_K,
    PARENT_P,
    PARENT_Q,
    _mode_mask,
    _observables,
    _physical_inner,
    _trapezoid,
    _trapezoid_atom_weights,
)
from src.signed_good_generated_epoch_time_telescope import (
    ACTUAL_HH_GENERATION_BRANCH,
    SignedGoodGeneratedWorkProvenance,
    signed_good_generated_epoch_telescope,
    signed_good_step_from_energy_reentry,
)


STATUS = (
    "DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "TWO_CONSECUTIVE_EXACT_SIGNED_GOOD_TRIADS__"
    "ONE_PDE_TRAJECTORY_AND_NATIVE_PHYSICAL_WORK__"
    "TWO_LAYER_TELESCOPE_TO_T0"
)

MIDDLE_CHILD = CHILD_K
MIDDLE_PARTNER = (2, 4, 4)
TOP_CHILD = (8, 4, 4)


def _negative(wavevector: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(-value for value in wavevector)


def _two_layer_initial_state(
    resolution: int,
    amplitude: float,
    middle_seed_weight: float,
    partner_polarization_angle: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Real divergence-free data for ``A+B=P`` and ``P+C=K``.

    The tiny phase-aligned ``P`` seed is permitted by the physical-energy gate:
    its energy must remain below one fifth of the later generated ``P`` energy.
    It makes the selected ``P+C`` work visible before competing two-interaction
    paths can masquerade as the designated parent-pair law.
    """
    n = int(resolution)
    state = np.zeros((3, n, n, n), dtype=complex)
    p_polarization = np.asarray((2.0, -3.0, 0.0)) / math.sqrt(13.0)
    q_polarization = np.asarray((0.0, 0.0, 1.0))
    partner_basis_one = np.asarray((0.0, 1.0, -1.0)) / math.sqrt(2.0)
    partner_basis_two = np.asarray((-4.0, 1.0, 1.0)) / math.sqrt(18.0)
    partner_polarization = (
        math.cos(partner_polarization_angle) * partner_basis_one
        + math.sin(partner_polarization_angle) * partner_basis_two
    )

    for wavevector, polarization in (
        (PARENT_P, p_polarization),
        (PARENT_Q, q_polarization),
        (MIDDLE_PARTNER, partner_polarization),
    ):
        positive = tuple(value % n for value in wavevector)
        negative = tuple((-value) % n for value in wavevector)
        state[(slice(None),) + positive] = polarization
        state[(slice(None),) + negative] = polarization

    seed = -1j * float(middle_seed_weight) * q_polarization
    positive = tuple(value % n for value in MIDDLE_CHILD)
    negative = tuple((-value) % n for value in MIDDLE_CHILD)
    state[(slice(None),) + positive] = seed
    state[(slice(None),) + negative] = np.conjugate(seed)

    state = _leray_dealias(state, k, k2, dealias)
    energy = _physical_inner(state, state, n)
    if not math.isfinite(energy) or energy <= 0.0:
        raise AssertionError("two-layer signed-good data lost positive physical energy")
    return float(amplitude) * state / math.sqrt(energy)


def _window_metrics(
    observations: Sequence[dict[str, float]],
    start_index: int,
    end_index: int,
    dt: float,
    viscosity: float,
) -> dict[str, object]:
    rows = tuple(observations[start_index : end_index + 1])
    designated_positive = tuple(row["designated_positive_work"] for row in rows)
    residual_positive = tuple(row["residual_positive_work"] for row in rows)
    designated_signed = _trapezoid(
        tuple(row["designated_signed_work"] for row in rows), dt
    )
    residual_signed = _trapezoid(
        tuple(row["residual_signed_work"] for row in rows), dt
    )
    gradient = _trapezoid(tuple(row["child_gradient"] for row in rows), dt)
    initial_energy = float(rows[0]["child_energy"])
    final_energy = float(rows[-1]["child_energy"])
    delta = final_energy - initial_energy
    balance_rhs = designated_signed + residual_signed - 2.0 * viscosity * gradient
    balance_scale = max(abs(delta), abs(balance_rhs), final_energy)
    if balance_scale <= 0.0:
        raise AssertionError("two-layer child carried no physical energy scale")
    return {
        "rows": rows,
        "initial_energy": initial_energy,
        "final_energy": final_energy,
        "designated_positive": designated_positive,
        "actual_positive_work": _trapezoid(designated_positive, dt),
        "residual_positive_work": _trapezoid(residual_positive, dt),
        "strain_action": _trapezoid(
            tuple(row["low_pass_strain_linf_upper"] for row in rows), dt
        ),
        "balance_relative_residual": abs(delta - balance_rhs) / balance_scale,
        "minimum_designated_forcing_energy": min(
            row["designated_forcing_energy"] for row in rows
        ),
    }


@dataclass(frozen=True)
class GalerkinSignedGoodTwoLayerRun:
    resolution: int
    steps: int
    dt: float
    duration: float
    viscosity: float
    amplitude: float
    middle_seed_weight: float
    partner_polarization_angle: float
    top_parent_child_ratio: float
    middle_parent_child_ratio: float
    initial_global_energy: float
    final_global_energy: float
    middle_initial_energy_fraction: float
    top_initial_energy_fraction: float
    middle_actual_positive_hh_work: float
    top_actual_positive_hh_work: float
    middle_energy_gate_lower: float
    top_energy_gate_lower: float
    middle_heavy_half_fraction: float
    top_heavy_half_fraction: float
    middle_residual_work_fraction: float
    top_residual_work_fraction: float
    middle_child_balance_relative_residual: float
    top_child_balance_relative_residual: float
    global_energy_balance_relative_residual: float
    maximum_global_nonlinear_work_relative_rate: float
    maximum_divergence_relative_to_initial_l2: float
    maximum_forcing_decomposition_relative_residual: float
    middle_low_pass_strain_action_upper: float
    top_low_pass_strain_action_upper: float
    middle_support_start: float
    middle_support_end: float
    top_support_start: float
    top_support_end: float
    first_common_reference_time: float
    last_common_reference_time: float
    epoch_layer_count: int
    epoch_hits_initial_boundary: bool


def simulate_signed_good_two_layer_galerkin(
    *,
    resolution: int,
    steps: int,
    viscosity: float = 0.02,
    amplitude: float = 96.0,
    scaled_lifetime: float = 0.05,
    middle_seed_weight: float = 5.0e-4,
    partner_polarization_angle: float = math.pi / 4.0,
) -> GalerkinSignedGoodTwoLayerRun:
    """Falsify two consecutive signed-good steps on one actual NS orbit."""
    n = int(resolution)
    count = int(steps)
    nu = float(viscosity)
    amp = float(amplitude)
    c = float(scaled_lifetime)
    seed_weight = float(middle_seed_weight)
    partner_angle = float(partner_polarization_angle)
    if n < 28 or n % 2 or count < 32 or count % 32:
        raise ValueError("resolution >=28 and an RK4 step count divisible by 32 required")
    if not all(math.isfinite(x) and x > 0.0 for x in (nu, amp, c, seed_weight)) or not math.isfinite(
        partner_angle
    ):
        raise ValueError("positive finite two-layer Galerkin parameters required")

    top_frequency = math.sqrt(sum(value * value for value in TOP_CHILD))
    middle_frequency = math.sqrt(sum(value * value for value in MIDDLE_CHILD))
    base_frequency = math.sqrt(sum(value * value for value in PARENT_P))
    top_lifetime = c / (top_frequency * top_frequency)
    duration = 1.6 * top_lifetime
    dt = duration / count
    top_start_index = 3 * count // 8
    top_mid_index = 11 * count // 16
    middle_mid_index = 15 * count // 32
    middle_end_index = 15 * count // 16

    k, k2, dealias = _spectral_geometry(n)
    base_pair_mask = _mode_mask(
        k,
        (PARENT_P, _negative(PARENT_P), PARENT_Q, _negative(PARENT_Q)),
    )
    middle_mask = _mode_mask(k, (MIDDLE_CHILD, _negative(MIDDLE_CHILD)))
    top_pair_mask = _mode_mask(
        k,
        (
            MIDDLE_CHILD,
            _negative(MIDDLE_CHILD),
            MIDDLE_PARTNER,
            _negative(MIDDLE_PARTNER),
        ),
    )
    top_mask = _mode_mask(k, (TOP_CHILD, _negative(TOP_CHILD)))
    if tuple(int(np.sum(mask)) for mask in (base_pair_mask, middle_mask, top_pair_mask, top_mask)) != (
        4,
        2,
        4,
        2,
    ):
        raise AssertionError("Galerkin lattice did not retain the exact two-layer triad chain")

    middle_low_pass = _strict_low_pass_symbol(k2, 0.25 * middle_frequency)
    top_low_pass = _strict_low_pass_symbol(k2, 0.25 * top_frequency)
    state = _two_layer_initial_state(
        n, amp, seed_weight, partner_angle, k, k2, dealias
    )
    times = tuple(dt * index for index in range(count + 1))
    middle_observations: list[dict[str, float]] = []
    top_observations: list[dict[str, float]] = []
    for index in range(count + 1):
        middle_observations.append(
            _observables(
                state,
                nu,
                k,
                k2,
                dealias,
                base_pair_mask,
                middle_mask,
                middle_low_pass,
            )
        )
        top_observations.append(
            _observables(
                state,
                nu,
                k,
                k2,
                dealias,
                top_pair_mask,
                top_mask,
                top_low_pass,
            )
        )
        if index < count:
            state = _rk4_step(state, dt, nu, k, k2, dealias)

    middle = _window_metrics(
        middle_observations, 0, middle_end_index, dt, nu
    )
    top = _window_metrics(top_observations, top_start_index, count, dt, nu)
    middle_times = times[: middle_end_index + 1]
    top_times = times[top_start_index:]

    middle_gate = dict(
        route_physical_energy_causality(
            terminal_energy=float(middle["final_energy"]),
            initial_energy=float(middle["initial_energy"]),
            residual_positive_work=float(middle["residual_positive_work"]),
            strain_action=float(middle["strain_action"]),
        )
    )
    top_gate = dict(
        route_physical_energy_causality(
            terminal_energy=float(top["final_energy"]),
            initial_energy=float(top["initial_energy"]),
            residual_positive_work=float(top["residual_positive_work"]),
            strain_action=float(top["strain_action"]),
        )
    )
    if middle_gate.get("branch") != ACTUAL_HH_GENERATION_BRANCH:
        raise AssertionError(
            "middle exact triad missed the physical HH gate: "
            f"branch={middle_gate.get('branch')}, "
            f"initial/final={float(middle['initial_energy']) / float(middle['final_energy']):.6e}, "
            f"residual/final={float(middle['residual_positive_work']) / float(middle['final_energy']):.6e}, "
            f"strain={float(middle['strain_action']):.6e}"
        )
    if top_gate.get("branch") != ACTUAL_HH_GENERATION_BRANCH:
        raise AssertionError(
            "top exact triad missed the physical HH gate: "
            f"branch={top_gate.get('branch')}, "
            f"initial/final={float(top['initial_energy']) / float(top['final_energy']):.6e}, "
            f"residual/final={float(top['residual_positive_work']) / float(top['final_energy']):.6e}, "
            f"strain={float(top['strain_action']):.6e}"
        )

    middle_half = dict(
        heavy_half_physical_transfer(
            times=middle_times,
            positive_work_weights=_trapezoid_atom_weights(
                middle["designated_positive"], dt
            ),
            slab_start=middle_times[0],
            slab_end=middle_times[-1],
        )
    )
    top_half = dict(
        heavy_half_physical_transfer(
            times=top_times,
            positive_work_weights=_trapezoid_atom_weights(
                top["designated_positive"], dt
            ),
            slab_start=top_times[0],
            slab_end=top_times[-1],
        )
    )
    if int(middle_half["half"]) != 1 or int(top_half["half"]) != 1:
        raise AssertionError(
            "actual two-layer work missed the later half-slabs: "
            f"middle={middle_half['half']}, top={top_half['half']}"
        )
    if abs(float(middle_half["start"]) - times[middle_mid_index]) > 8.0e-12 * (
        c / (middle_frequency * middle_frequency)
    ):
        raise AssertionError("middle physical half-slab lost its native time boundary")
    if abs(float(top_half["start"]) - times[top_mid_index]) > 8.0e-12 * top_lifetime:
        raise AssertionError("top physical half-slab lost its native time boundary")

    trajectory_id = f"unforced-dealiased-NS-two-layer-n{n}"
    top_provenance = SignedGoodGeneratedWorkProvenance(
        event_id=f"two-layer-top-event-n{n}",
        trajectory_id=trajectory_id,
        child_carrier_id="fourier-carrier-(+/-8,+/-4,+/-4)",
        generated_parent_carrier_id="fourier-carrier-(+/-6,0,0)",
        work_law_id=f"actual-positive-middle-partner-top-work-n{n}",
        child_frequency=top_frequency,
        parent_frequency=middle_frequency,
        scaled_lifetime=c,
        slab_start=top_times[0],
        slab_end=top_times[-1],
    )
    middle_provenance = SignedGoodGeneratedWorkProvenance(
        event_id=f"two-layer-middle-event-n{n}",
        trajectory_id=trajectory_id,
        child_carrier_id="fourier-carrier-(+/-6,0,0)",
        generated_parent_carrier_id="fourier-carrier-|p|=|q|=sqrt(13)",
        work_law_id=f"actual-positive-pq-middle-work-n{n}",
        child_frequency=middle_frequency,
        parent_frequency=base_frequency,
        scaled_lifetime=c,
        slab_start=middle_times[0],
        slab_end=middle_times[-1],
    )
    for gate, half, provenance in (
        (top_gate, top_half, top_provenance),
        (middle_gate, middle_half, middle_provenance),
    ):
        gate["provenance"] = provenance
        half["provenance"] = provenance

    top_step = signed_good_step_from_energy_reentry(
        reentry=top_gate,
        selected_physical_half_slab=top_half,
        child_frequency=top_frequency,
        parent_frequency=middle_frequency,
        scaled_lifetime=c,
    )
    middle_step = signed_good_step_from_energy_reentry(
        reentry=middle_gate,
        selected_physical_half_slab=middle_half,
        child_frequency=middle_frequency,
        parent_frequency=base_frequency,
        scaled_lifetime=c,
    )
    epoch = signed_good_generated_epoch_telescope((top_step, middle_step))

    for label, window, gate in (
        ("middle", middle, middle_gate),
        ("top", top, top_gate),
    ):
        work = float(window["actual_positive_work"])
        lower = float(gate["physical_hh_work_lower"])
        if work + 2.0e-5 * max(work, lower) < lower:
            raise AssertionError(
                f"{label} measured HH work {work:.6e} fell below gate lower {lower:.6e}"
            )
        if float(window["balance_relative_residual"]) > 5.0e-3:
            raise AssertionError(f"{label} child lost its measured NS energy balance")
        if float(window["minimum_designated_forcing_energy"]) <= 0.0:
            raise AssertionError(f"{label} exact parent pair lost its physical forcing")

    initial_global = float(top_observations[0]["global_energy"])
    final_global = float(top_observations[-1]["global_energy"])
    global_gradient = _trapezoid(
        tuple(row["global_gradient"] for row in top_observations), dt
    )
    global_balance_rhs = -2.0 * nu * global_gradient
    global_balance_residual = abs(
        final_global - initial_global - global_balance_rhs
    ) / initial_global
    nonlinear_rate_scale = initial_global / duration
    maximum_global_nonlinear = max(
        abs(row["global_nonlinear_work"]) for row in top_observations
    ) / nonlinear_rate_scale
    maximum_divergence = max(
        row["divergence_norm"] for row in top_observations
    ) / math.sqrt(initial_global)
    maximum_decomposition = max(
        row["forcing_decomposition_relative_residual"]
        for row in (*middle_observations, *top_observations)
    )
    if global_balance_residual > 2.0e-5:
        raise AssertionError("two-layer Galerkin orbit lost global viscous energy balance")
    if maximum_global_nonlinear > 2.0e-11:
        raise AssertionError("two-layer Galerkin nonlinearity lost global energy skewness")
    if maximum_divergence > 2.0e-11:
        raise AssertionError("two-layer Galerkin orbit lost incompressibility")
    if maximum_decomposition > 2.0e-12:
        raise AssertionError("two-layer designated and residual work lost exact decomposition")
    if epoch.layer_count != 2 or not epoch.hits_initial_boundary:
        raise AssertionError("actual two-layer work did not telescope from the interior to t=0")

    return GalerkinSignedGoodTwoLayerRun(
        resolution=n,
        steps=count,
        dt=dt,
        duration=duration,
        viscosity=nu,
        amplitude=amp,
        middle_seed_weight=seed_weight,
        partner_polarization_angle=partner_angle,
        top_parent_child_ratio=middle_frequency / top_frequency,
        middle_parent_child_ratio=base_frequency / middle_frequency,
        initial_global_energy=initial_global,
        final_global_energy=final_global,
        middle_initial_energy_fraction=float(middle["initial_energy"])
        / float(middle["final_energy"]),
        top_initial_energy_fraction=float(top["initial_energy"])
        / float(top["final_energy"]),
        middle_actual_positive_hh_work=float(middle["actual_positive_work"]),
        top_actual_positive_hh_work=float(top["actual_positive_work"]),
        middle_energy_gate_lower=float(middle_gate["physical_hh_work_lower"]),
        top_energy_gate_lower=float(top_gate["physical_hh_work_lower"]),
        middle_heavy_half_fraction=float(middle_half["mass"])
        / float(middle_half["total"]),
        top_heavy_half_fraction=float(top_half["mass"]) / float(top_half["total"]),
        middle_residual_work_fraction=float(middle["residual_positive_work"])
        / float(middle["final_energy"]),
        top_residual_work_fraction=float(top["residual_positive_work"])
        / float(top["final_energy"]),
        middle_child_balance_relative_residual=float(middle["balance_relative_residual"]),
        top_child_balance_relative_residual=float(top["balance_relative_residual"]),
        global_energy_balance_relative_residual=global_balance_residual,
        maximum_global_nonlinear_work_relative_rate=maximum_global_nonlinear,
        maximum_divergence_relative_to_initial_l2=maximum_divergence,
        maximum_forcing_decomposition_relative_residual=maximum_decomposition,
        middle_low_pass_strain_action_upper=float(middle["strain_action"]),
        top_low_pass_strain_action_upper=float(top["strain_action"]),
        middle_support_start=float(middle_half["start"]),
        middle_support_end=float(middle_half["end"]),
        top_support_start=float(top_half["start"]),
        top_support_end=float(top_half["end"]),
        first_common_reference_time=epoch.first_common_reference_time,
        last_common_reference_time=epoch.last_common_reference_time,
        epoch_layer_count=epoch.layer_count,
        epoch_hits_initial_boundary=epoch.hits_initial_boundary,
    )


@dataclass(frozen=True)
class SignedGoodTwoLayerPDEProbe:
    status: str
    runs: tuple[GalerkinSignedGoodTwoLayerRun, ...]
    top_work_resolution_spread: float
    middle_work_resolution_spread: float


def run_probe(
    resolutions: Sequence[int] = (28, 32),
    *,
    steps: int = 128,
    viscosity: float = 0.02,
    amplitude: float = 96.0,
    scaled_lifetime: float = 0.05,
    middle_seed_weight: float = 5.0e-4,
    partner_polarization_angle: float = math.pi / 4.0,
) -> SignedGoodTwoLayerPDEProbe:
    runs = tuple(
        simulate_signed_good_two_layer_galerkin(
            resolution=int(resolution),
            steps=int(steps),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            scaled_lifetime=float(scaled_lifetime),
            middle_seed_weight=float(middle_seed_weight),
            partner_polarization_angle=float(partner_polarization_angle),
        )
        for resolution in resolutions
    )
    if not runs:
        raise ValueError("at least one two-layer Galerkin resolution required")
    top_work = tuple(row.top_actual_positive_hh_work for row in runs)
    middle_work = tuple(row.middle_actual_positive_hh_work for row in runs)
    top_spread = (max(top_work) - min(top_work)) / max(top_work)
    middle_spread = (max(middle_work) - min(middle_work)) / max(middle_work)
    if top_spread > 5.0e-2 or middle_spread > 5.0e-2:
        raise AssertionError("two-layer physical work did not stabilize under refinement")
    return SignedGoodTwoLayerPDEProbe(
        status=STATUS,
        runs=runs,
        top_work_resolution_spread=top_spread,
        middle_work_resolution_spread=middle_spread,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", type=int, nargs="+", default=(28, 32))
    ap.add_argument("--steps", type=int, default=128)
    ap.add_argument("--viscosity", type=float, default=0.02)
    ap.add_argument("--amplitude", type=float, default=96.0)
    ap.add_argument("--scaled-lifetime", type=float, default=0.05)
    ap.add_argument("--middle-seed-weight", type=float, default=5.0e-4)
    ap.add_argument("--partner-polarization-angle", type=float, default=math.pi / 4.0)
    ap.add_argument(
        "--outdir",
        type=Path,
        default=Path("results-signed-good-generated-epoch-two-layer-pde-probe"),
    )
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = run_probe(
        args.resolutions,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        scaled_lifetime=args.scaled_lifetime,
        middle_seed_weight=args.middle_seed_weight,
        partner_polarization_angle=args.partner_polarization_angle,
    )
    (args.outdir / "signed_good_generated_epoch_two_layer_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    table = "\n".join(
        "| {n} | {middle:.3e} | {top:.3e} | {mi:.3e} | {ti:.3e} | {mr:.3e} | {tr:.3e} |".format(
            n=row.resolution,
            middle=row.middle_actual_positive_hh_work,
            top=row.top_actual_positive_hh_work,
            mi=row.middle_initial_energy_fraction,
            ti=row.top_initial_energy_fraction,
            mr=row.middle_residual_work_fraction,
            tr=row.top_residual_work_fraction,
        )
        for row in result.runs
    )
    md = f"""# Signed-good two-layer physical PDE probe

Status: **{STATUS}**.

`(3,2,0)+(3,-2,0)=(6,0,0)` and
`(6,0,0)+(2,4,4)=(8,4,4)` are followed on one unforced, dealiased,
incompressible Fourier--Galerkin Navier--Stokes trajectory.  The exact ratios
are `sqrt(13)/6` and `6/sqrt(96)`, both strictly in `(3/5,5/8)`.

| resolution | middle work | top work | middle initial/final E | top initial/final E | middle residual/final E | top residual/final E |
|---:|---:|---:|---:|---:|---:|---:|
{table}

Work spreads: middle `{result.middle_work_resolution_spread:.3e}`, top
`{result.top_work_resolution_spread:.3e}`.  Both physical heavy-half supports
form one certified two-layer epoch on the same PDE history and reach `t=0`.
This is numerical falsification evidence, not a continuum or generic-HH proof.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
