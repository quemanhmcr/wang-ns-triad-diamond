from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.continuum_master_event_quotient import energy_reentry_master_route
from src.full_natural_service_corridor_pde_probe import (
    _analytic_divergence_free_initial_state,
    _nonlinear_term,
    _observables,
    _rk4_step,
    _spectral_geometry,
    _spectral_inner,
)
from src.same_carrier_checkpoint_segmentation_pde_probe import _bilinear_convection
from src.smooth_quadratic_carrier_interface import (
    RELINK_OWNER,
    STRAIN_OWNER,
    GaugeQuotientedInterfaceWork,
    positive_smooth_interface_split,
)
from src.smooth_relink_donor_quotient import (
    SMOOTH_RELINK_SAME_EVENT_RELAY,
    smooth_relink_donor_quotient,
)


STATUS = (
    "DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES__"
    "ACTUAL_RESOLVED_LINEARIZED_OPERATOR__"
    "SMOOTH_KPHYS_PAIR_FLUX_AND_DONOR_CLOSURE__"
    "NATIVE_WORK_MASTER_REPLAY"
)


def _smooth_quadratic_role_effects(
    k2: np.ndarray,
    carrier_frequency: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A fixed smooth positive Fourier partition eta_0+eta_1+eta_2=1.

    Put s=|k|^2/(|k|^2+(A/2)^2) and use the quadratic Bernstein
    partition ((1-s)^2,2s(1-s),s^2).  The analysis multipliers are the
    positive square roots A_a=sqrt(eta_a), so these eta_a are precisely the
    Q^2 energy effects required by the smooth-interface theorem.
    """
    A = float(carrier_frequency)
    transition2 = (0.5 * A) ** 2
    s = k2 / (k2 + transition2)
    return (1.0 - s) ** 2, 2.0 * s * (1.0 - s), s**2


def _linearized_resolved(
    resolved_hat: np.ndarray,
    field_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """The actual Galerkin L_V f=B(V,f)+B(f,V), including Leray projection."""
    return _bilinear_convection(resolved_hat, field_hat, k, k2, dealias) + _bilinear_convection(
        field_hat,
        resolved_hat,
        k,
        k2,
        dealias,
    )


def _field_norm(field_hat: np.ndarray, resolution: int) -> float:
    return math.sqrt(max(0.0, _spectral_inner(field_hat, field_hat, resolution)))


def _relative(residual: float, scale: float) -> float:
    return abs(float(residual)) / max(abs(float(scale)), 1.0e-300)


@dataclass(frozen=True)
class GalerkinSmoothRelinkRun:
    resolution: int
    steps: int
    dt: float
    snapshots: int
    positive_native_work_snapshots: int
    relink_owner_snapshots: int
    mixed_relink_strain_snapshots: int
    maximum_divergence_norm: float
    maximum_global_nonlinear_work: float
    global_energy_balance_relative_residual: float
    maximum_role_partition_residual: float
    maximum_linearized_partition_relative_residual: float
    maximum_pair_antisymmetry_relative_residual: float
    maximum_strain_pair_symmetry_relative_residual: float
    maximum_relink_row_binding_relative_residual: float
    maximum_strain_row_binding_relative_residual: float
    maximum_native_split_relative_residual: float
    maximum_total_relink_relative_residual: float
    minimum_incoming_minus_gain_relative_margin: float
    maximum_shortest_donor_path_length: int
    maximum_positive_relink_work: float
    final_positive_relink_work: float
    maximum_positive_native_interface_work: float
    master_route_failures: int


def _snapshot_interface_diagnostics(
    state_hat: np.ndarray,
    *,
    carrier_frequency: float,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> dict[str, object]:
    """Read the K/S pair laws from one evolved incompressible NS state.

    No dense surrogate operator is assembled.  Every pairing is evaluated by
    applying the actual dealiased Galerkin bilinear form.  For the real Hilbert
    pairing and L=L_V,

      -2<x,K y> = -( <x,Ly>-<Lx,y> ),
      -2<x,S y> = -( <x,Ly>+<Lx,y> ).

    This exposes the adjoint split through its defining physical work form.
    """
    n = int(state_hat.shape[1])
    A = float(carrier_frequency)
    radius = np.sqrt(k2)
    resolved_hat = (radius <= 0.25 * A)[None, ...] * state_hat
    etas = _smooth_quadratic_role_effects(k2, A)
    role_partition_residual = float(np.max(np.abs(sum(etas) - 1.0)))
    roles = tuple(eta[None, ...] * state_hat for eta in etas)
    role_images = tuple(
        _linearized_resolved(resolved_hat, role, k, k2, dealias) for role in roles
    )
    full_image = _linearized_resolved(resolved_hat, state_hat, k, k2, dealias)
    image_partition_defect = full_image - sum(role_images)
    image_scale = max(_field_norm(full_image, n), 1.0e-300)
    linearized_partition_relative = _field_norm(image_partition_defect, n) / image_scale

    role_count = len(roles)
    relink_pair = np.empty((role_count, role_count), dtype=float)
    strain_pair = np.empty((role_count, role_count), dtype=float)
    for a in range(role_count):
        for b in range(role_count):
            left = _spectral_inner(roles[a], role_images[b], n)
            adjoint_left = _spectral_inner(role_images[a], roles[b], n)
            relink_pair[a, b] = -(left - adjoint_left)
            strain_pair[a, b] = -(left + adjoint_left)

    direct_relink = np.empty(role_count, dtype=float)
    direct_strain = np.empty(role_count, dtype=float)
    direct_native = np.empty(role_count, dtype=float)
    for a in range(role_count):
        left = _spectral_inner(roles[a], full_image, n)
        adjoint_left = _spectral_inner(role_images[a], state_hat, n)
        direct_relink[a] = -(left - adjoint_left)
        direct_strain[a] = -(left + adjoint_left)
        direct_native[a] = -2.0 * left

    native_scale = max(
        float(np.max(np.abs(relink_pair))),
        float(np.max(np.abs(strain_pair))),
        float(np.max(np.abs(direct_native))),
        1.0e-300,
    )
    relink_row_residual = float(np.max(np.abs(relink_pair.sum(axis=1) - direct_relink)))
    strain_row_residual = float(np.max(np.abs(strain_pair.sum(axis=1) - direct_strain)))
    native_split_residual = float(np.max(np.abs(direct_native - direct_relink - direct_strain)))
    pair_antisymmetry_residual = float(np.max(np.abs(relink_pair + relink_pair.T)))
    strain_pair_symmetry_residual = float(np.max(np.abs(strain_pair - strain_pair.T)))
    total_relink_residual = abs(float(direct_relink.sum()))

    work = GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=tuple(float(x) for x in direct_native),
        signed_physical_relink_atoms=tuple(float(x) for x in direct_relink),
        signed_existing_strain_atoms=tuple(float(x) for x in direct_strain),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=tuple(
            tuple(float(x) for x in row) for row in relink_pair
        ),
    )
    split = positive_smooth_interface_split(work)
    owners = tuple(str(x) for x in split["joint_physical_owners"])
    route = None
    donor = None
    positive_native = float(split["positive_native_interface_work"])
    if owners:
        route = energy_reentry_master_route(
            "positive native smooth interface work",
            positive_native,
            {
                "branch": "smooth_interface_physical_work",
                "joint_interface_owners": owners,
                "coefficient_impulse_used_as_physical_work": False,
                "observer_partition_motion_charged_as_physics": False,
                "gauge_quotiented_interface_work_certificate": work,
            },
        )
    if RELINK_OWNER in owners:
        donor = smooth_relink_donor_quotient(work)

    return {
        "work": work,
        "split": split,
        "owners": owners,
        "route": route,
        "donor": donor,
        "role_partition_residual": role_partition_residual,
        "linearized_partition_relative_residual": linearized_partition_relative,
        "pair_antisymmetry_relative_residual": _relative(
            pair_antisymmetry_residual, native_scale
        ),
        "strain_pair_symmetry_relative_residual": _relative(
            strain_pair_symmetry_residual, native_scale
        ),
        "relink_row_binding_relative_residual": _relative(
            relink_row_residual, native_scale
        ),
        "strain_row_binding_relative_residual": _relative(
            strain_row_residual, native_scale
        ),
        "native_split_relative_residual": _relative(native_split_residual, native_scale),
        "total_relink_relative_residual": _relative(total_relink_residual, native_scale),
        "positive_relink_work": float(split["positive_conservative_relink_work"]),
        "positive_native_interface_work": positive_native,
        "native_scale": native_scale,
    }


def simulate_smooth_relink_galerkin(
    *,
    resolution: int,
    steps: int,
    duration: float,
    viscosity: float,
    carrier_frequency: float,
) -> GalerkinSmoothRelinkRun:
    n = int(resolution)
    count = int(steps)
    T = float(duration)
    nu = float(viscosity)
    A = float(carrier_frequency)
    if count < 4 or min(T, nu, A) <= 0.0 or not all(
        math.isfinite(x) for x in (T, nu, A)
    ):
        raise ValueError("positive finite PDE parameters and at least four RK4 steps required")
    dt = T / count
    k, k2, dealias = _spectral_geometry(n)
    state = _analytic_divergence_free_initial_state(n, k, k2, dealias)
    before_observable = _observables(state, nu, A, k, k2, dealias)
    initial_observable = before_observable

    positive_snapshots = 0
    relink_snapshots = 0
    mixed_snapshots = 0
    route_failures = 0
    maximum_divergence = 0.0
    maximum_global_nonlinear = 0.0
    maximum_partition = 0.0
    maximum_linearized_partition = 0.0
    maximum_antisymmetry = 0.0
    maximum_strain_symmetry = 0.0
    maximum_relink_row = 0.0
    maximum_strain_row = 0.0
    maximum_native_split = 0.0
    maximum_total_relink = 0.0
    minimum_donor_margin = float("inf")
    maximum_path = 0
    maximum_positive_relink = 0.0
    final_positive_relink = 0.0
    maximum_positive_native = 0.0
    global_action = 0.0

    for index in range(count + 1):
        observable = before_observable
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        maximum_divergence = max(maximum_divergence, float(observable["divergence_norm"]))
        maximum_global_nonlinear = max(
            maximum_global_nonlinear,
            abs(float(observable["global_nonlinear_work"])),
            abs(2.0 * _spectral_inner(state, nonlinear, n)),
        )
        diagnostics = _snapshot_interface_diagnostics(
            state,
            carrier_frequency=A,
            k=k,
            k2=k2,
            dealias=dealias,
        )
        owners = tuple(diagnostics["owners"])
        split = diagnostics["split"]
        route = diagnostics["route"]
        donor = diagnostics["donor"]
        positive_native = float(diagnostics["positive_native_interface_work"])
        positive_relink = float(diagnostics["positive_relink_work"])
        if positive_native > 0.0:
            positive_snapshots += 1
        if RELINK_OWNER in owners:
            relink_snapshots += 1
            if donor is None or route is None:
                route_failures += 1
            else:
                certificate = donor["certificate"]
                scale = float(diagnostics["native_scale"])
                margin = (
                    float(donor["recipient_positive_incoming_flux"])
                    - float(donor["positive_relink_work"])
                ) / scale
                minimum_donor_margin = min(minimum_donor_margin, margin)
                maximum_path = max(
                    maximum_path,
                    int(donor["maximum_shortest_donor_path_length"]),
                )
                if route.same_event_relays != (SMOOTH_RELINK_SAME_EVENT_RELAY,):
                    route_failures += 1
                if route.smooth_relink_donor_certificate != certificate:
                    route_failures += 1
        if RELINK_OWNER in owners and STRAIN_OWNER in owners:
            mixed_snapshots += 1
        if route is not None:
            if STRAIN_OWNER in owners:
                if route.owner_bundle is None:
                    route_failures += 1
                else:
                    expected_strain = float(split["positive_existing_strain_work"])
                    if route.owner_bundle.physical_measure != "positive existing smooth strain work":
                        route_failures += 1
                    if abs(route.owner_bundle.mass - expected_strain) > 8.0e-12 * max(
                        expected_strain, 1.0e-300
                    ):
                        route_failures += 1
            elif route.owner_bundle is not None:
                route_failures += 1
            if abs(route.mass - positive_native) > 8.0e-12 * max(
                positive_native, 1.0e-300
            ):
                route_failures += 1

        maximum_partition = max(
            maximum_partition, float(diagnostics["role_partition_residual"])
        )
        maximum_linearized_partition = max(
            maximum_linearized_partition,
            float(diagnostics["linearized_partition_relative_residual"]),
        )
        maximum_antisymmetry = max(
            maximum_antisymmetry,
            float(diagnostics["pair_antisymmetry_relative_residual"]),
        )
        maximum_strain_symmetry = max(
            maximum_strain_symmetry,
            float(diagnostics["strain_pair_symmetry_relative_residual"]),
        )
        maximum_relink_row = max(
            maximum_relink_row,
            float(diagnostics["relink_row_binding_relative_residual"]),
        )
        maximum_strain_row = max(
            maximum_strain_row,
            float(diagnostics["strain_row_binding_relative_residual"]),
        )
        maximum_native_split = max(
            maximum_native_split,
            float(diagnostics["native_split_relative_residual"]),
        )
        maximum_total_relink = max(
            maximum_total_relink,
            float(diagnostics["total_relink_relative_residual"]),
        )
        maximum_positive_relink = max(maximum_positive_relink, positive_relink)
        maximum_positive_native = max(maximum_positive_native, positive_native)
        if index == count:
            final_positive_relink = positive_relink
            break

        state = _rk4_step(state, dt, nu, k, k2, dealias)
        after_observable = _observables(state, nu, A, k, k2, dealias)
        global_action += 0.5 * dt * (
            float(before_observable["global_power"])
            + float(after_observable["global_power"])
        )
        before_observable = after_observable

    global_delta = float(before_observable["global_energy"]) - float(
        initial_observable["global_energy"]
    )
    global_scale = max(
        abs(global_delta),
        abs(global_action),
        float(initial_observable["global_energy"]),
        1.0e-300,
    )
    if not math.isfinite(minimum_donor_margin):
        minimum_donor_margin = -float("inf")
    return GalerkinSmoothRelinkRun(
        resolution=n,
        steps=count,
        dt=dt,
        snapshots=count + 1,
        positive_native_work_snapshots=positive_snapshots,
        relink_owner_snapshots=relink_snapshots,
        mixed_relink_strain_snapshots=mixed_snapshots,
        maximum_divergence_norm=maximum_divergence,
        maximum_global_nonlinear_work=maximum_global_nonlinear,
        global_energy_balance_relative_residual=abs(global_delta - global_action) / global_scale,
        maximum_role_partition_residual=maximum_partition,
        maximum_linearized_partition_relative_residual=maximum_linearized_partition,
        maximum_pair_antisymmetry_relative_residual=maximum_antisymmetry,
        maximum_strain_pair_symmetry_relative_residual=maximum_strain_symmetry,
        maximum_relink_row_binding_relative_residual=maximum_relink_row,
        maximum_strain_row_binding_relative_residual=maximum_strain_row,
        maximum_native_split_relative_residual=maximum_native_split,
        maximum_total_relink_relative_residual=maximum_total_relink,
        minimum_incoming_minus_gain_relative_margin=minimum_donor_margin,
        maximum_shortest_donor_path_length=maximum_path,
        maximum_positive_relink_work=maximum_positive_relink,
        final_positive_relink_work=final_positive_relink,
        maximum_positive_native_interface_work=maximum_positive_native,
        master_route_failures=route_failures,
    )


@dataclass(frozen=True)
class SmoothRelinkPhysicalPDEProbe:
    status: str
    duration: float
    viscosity: float
    carrier_frequency: float
    runs: tuple[GalerkinSmoothRelinkRun, ...]
    final_positive_relink_resolution_spread: float


def run_smooth_relink_physical_pde_probe(
    *,
    resolutions: Sequence[int] = (12, 16, 20),
    steps: int = 48,
    duration: float = 0.01,
    viscosity: float = 0.05,
    carrier_frequency: float = 4.0,
) -> SmoothRelinkPhysicalPDEProbe:
    runs = tuple(
        simulate_smooth_relink_galerkin(
            resolution=int(resolution),
            steps=int(steps),
            duration=float(duration),
            viscosity=float(viscosity),
            carrier_frequency=float(carrier_frequency),
        )
        for resolution in resolutions
    )
    if not runs:
        raise ValueError("at least one Fourier-Galerkin resolution required")
    for run in runs:
        if run.maximum_divergence_norm > 2.0e-11:
            raise AssertionError("smooth-relink NS trajectory lost incompressibility")
        if run.maximum_global_nonlinear_work > 2.0e-10:
            raise AssertionError("dealiased NS nonlinearity lost global energy conservation")
        if run.global_energy_balance_relative_residual > 2.0e-5:
            raise AssertionError("global NS energy balance failed on the smooth-relink trajectory")
        if run.maximum_role_partition_residual > 2.0e-15:
            raise AssertionError("smooth Q^2 role effects lost their exact partition")
        if run.maximum_linearized_partition_relative_residual > 3.0e-12:
            raise AssertionError("actual L_V failed to respect the smooth role partition")
        if run.maximum_pair_antisymmetry_relative_residual > 3.0e-12:
            raise AssertionError("actual K_phys pair work lost antisymmetry")
        if run.maximum_strain_pair_symmetry_relative_residual > 3.0e-12:
            raise AssertionError("actual S pair work lost symmetry")
        if run.maximum_relink_row_binding_relative_residual > 3.0e-12:
            raise AssertionError("actual K_phys pair rows failed to reconstruct relink work")
        if run.maximum_strain_row_binding_relative_residual > 3.0e-12:
            raise AssertionError("actual S pair rows failed to reconstruct strain work")
        if run.maximum_native_split_relative_residual > 3.0e-12:
            raise AssertionError("actual native interface work failed its K_phys/S split")
        if run.maximum_total_relink_relative_residual > 3.0e-12:
            raise AssertionError("actual smooth K_phys relink created net energy")
        if run.positive_native_work_snapshots != run.snapshots:
            raise AssertionError("evolved trajectory lost positive native interface work")
        if run.relink_owner_snapshots <= 0:
            raise AssertionError("PDE probe never entered the smooth relink owner branch")
        if run.maximum_positive_relink_work <= 1.0e-10:
            raise AssertionError("PDE probe degenerated to numerically zero K_phys relink")
        if run.minimum_incoming_minus_gain_relative_margin < -1.0e-10:
            raise AssertionError("actual relink recipient gain exceeded incoming same-event flux")
        if run.maximum_shortest_donor_path_length > 2:
            raise AssertionError("three-role donor path exceeded the finite simple-path bound")
        if run.master_route_failures:
            raise AssertionError("actual native work failed master replay/component routing")

    final_relink = [run.final_positive_relink_work for run in runs]
    spread = (max(final_relink) - min(final_relink)) / max(max(final_relink), 1.0e-300)
    if spread > 5.0e-2:
        raise AssertionError("actual K_phys relink work did not stabilize under resolution refinement")
    return SmoothRelinkPhysicalPDEProbe(
        status=STATUS,
        duration=float(duration),
        viscosity=float(viscosity),
        carrier_frequency=float(carrier_frequency),
        runs=runs,
        final_positive_relink_resolution_spread=spread,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolutions", type=int, nargs="+", default=[12, 16, 20])
    parser.add_argument("--steps", type=int, default=48)
    parser.add_argument("--duration", type=float, default=0.01)
    parser.add_argument("--viscosity", type=float, default=0.05)
    parser.add_argument("--carrier-frequency", type=float, default=4.0)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("results-smooth-relink-donor-pde-probe"),
    )
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = run_smooth_relink_physical_pde_probe(
        resolutions=args.resolutions,
        steps=args.steps,
        duration=args.duration,
        viscosity=args.viscosity,
        carrier_frequency=args.carrier_frequency,
    )
    (args.outdir / "smooth_relink_donor_pde_probe.json").write_text(
        json.dumps(asdict(result), indent=2), encoding="utf-8"
    )
    lines = [
        "# Smooth relink donor: physical PDE probe",
        "",
        f"Status: **{result.status}**.",
        "",
        "The trajectory is the unforced 3D incompressible Navier--Stokes Fourier-Galerkin system with Leray projection, viscosity, 2/3 dealiasing and RK4.  At every snapshot the probe applies the actual resolved linearized operator `L_V f=B(V,f)+B(f,V)` and reads its adjoint `K_phys/S` work split through physical pairings.  No proxy evolution or random matrix supplies the reported donor law.  This is numerical falsification evidence, not a continuum proof.",
        "",
        f"Physical interval: `T={result.duration:.12g}`, `A={result.carrier_frequency:.12g}`, `nu={result.viscosity:.12g}`.",
        "",
        "| N | steps | relink snapshots | mixed snapshots | div | global balance | K antisym | K rows | S rows | donor margin | max path | master failures |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for run in result.runs:
        lines.append(
            f"| {run.resolution} | {run.steps} | {run.relink_owner_snapshots}/{run.snapshots} | "
            f"{run.mixed_relink_strain_snapshots}/{run.snapshots} | {run.maximum_divergence_norm:.3e} | "
            f"{run.global_energy_balance_relative_residual:.3e} | "
            f"{run.maximum_pair_antisymmetry_relative_residual:.3e} | "
            f"{run.maximum_relink_row_binding_relative_residual:.3e} | "
            f"{run.maximum_strain_row_binding_relative_residual:.3e} | "
            f"{run.minimum_incoming_minus_gain_relative_margin:.3e} | "
            f"{run.maximum_shortest_donor_path_length} | {run.master_route_failures} |"
        )
    lines.extend(
        (
            "",
            f"Final positive relink-work resolution spread: `{result.final_positive_relink_resolution_spread:.3e}`.",
            "",
            "The smooth roles are one fixed positive quadratic Fourier partition `eta_a=A_a^2`, the resolved field and every work pairing come from the same evolved PDE state, and the master receives the replayed positive native work while any surviving strain bundle receives only its own positive component.",
        )
    )
    summary = "\n".join(lines) + "\n"
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
