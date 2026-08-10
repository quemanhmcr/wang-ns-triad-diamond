from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    INHERIT_FRACTION,
    RESIDUAL_FRACTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
    exact_adjoint_residual,
)
from src.critical_annular_carrier_service_reentry import (
    BOUNDED_HEAT_RADIUS,
    bounded_heat_defect_fraction_lower,
    persistent_carrier_critical_mass_lower as high_strain_persistent_carrier_mass_lower,
    renewed_analysis_probe_growth_upper,
    uniform_bounded_square_service_lower as high_strain_uniform_service_lower,
)
from src.high_strain_resolved_ancestor import (
    high_strain_ancestor_mass_threshold,
    retained_dissipation_lower,
)
from src.nn_critical_heat_carrier_seed import (
    LOW_STRAIN_ACTION,
    persistent_seed_low_low_gap,
    renewal_scale,
)
from src.nn_seed_temporal_first_stop import backward_natural_endpoint
from src.smooth_sgs_first_hit_extraction import (
    PhysicalPathMonitor,
    ThresholdTopology,
    first_physical_corridor_exit,
    rescale_monitor_units,
)

SHELL_TO_RENEWAL_RATIO = 3.0 / 4.0
DEFAULT_DOMINANT_FRESH_FRACTION = 1.0 / 4.0


def critical_shell_terminal_mass_lower(shell_critical_mass_lower: float) -> float:
    """Lower bound for A|z(t)|^2 after whole-shell registration.

    If mu0 <= M||P_M u(t)||_2^2 and A=3M/4, choose Q_A=1 on
    {M/2<|xi|<=M} and psi=P_Mu/||P_Mu||.  Then

        A|z(t)|^2 = A||P_Mu||_2^2 >= (3/4) mu0.

    No coherent cell, packet maximizer, or material label appears.
    """
    mu0 = float(shell_critical_mass_lower)
    if mu0 <= 0 or not math.isfinite(mu0):
        raise ValueError("positive finite critical shell mass lower required")
    return SHELL_TO_RENEWAL_RATIO * mu0


def critical_shell_survivor_coefficient_mass_lower(shell_critical_mass_lower: float) -> float:
    """Clean A|z(s)|^2 lower on every prefix of a full no-hit corridor."""
    return (INHERIT_FRACTION**2) * critical_shell_terminal_mass_lower(shell_critical_mass_lower)


def critical_shell_persistent_carrier_mass_lower(
    shell_critical_mass_lower: float,
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Uniform A||Q_Au(s)||_2^2 lower throughout a full no-hit natural slab.

    The coefficient survives by 1/4 on every prefix.  The registered
    affine/Kelvin/viscous analysis dual has norm-growth factor J, while all
    nonaffine/full-transport mismatch remains in the separately monitored
    role-interface coefficient obstruction.  Cauchy therefore gives

        A||Q_Au(s)||_2^2 >= (3 mu0)/(64 J^2).
    """
    J = renewed_analysis_probe_growth_upper(scaled_lifetime, viscosity, strain_action)
    return critical_shell_survivor_coefficient_mass_lower(shell_critical_mass_lower) / (J * J)


def critical_shell_bounded_service_lower(
    shell_critical_mass_lower: float,
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Uniform own-scale bounded-displacement square-service threshold.

    The existing Arb-certified radius-3 truncation retains the clean fraction
    q_b of the annular heat defect.  Hence some actual |r_s|<=3/A satisfies

        A||delta_{r_s}Q_Au(s)||_2^2 >= q_b * 3 mu0/(64 J^2).
    """
    q_b = bounded_heat_defect_fraction_lower(strain_action)
    return q_b * critical_shell_persistent_carrier_mass_lower(
        shell_critical_mass_lower,
        scaled_lifetime,
        viscosity,
        strain_action,
    )


def critical_shell_integrated_service_lower(
    shell_critical_mass_lower: float,
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Normalized bounded heat service on one full A-natural slab: c*Y_shell."""
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    return c * critical_shell_bounded_service_lower(
        shell_critical_mass_lower,
        c,
        viscosity,
        strain_action,
    )


def dissipation_supplier_shell_mass_threshold(
    resolved_dissipation_lower: float,
    scaled_lifetime: float,
) -> float:
    """Critical-shell threshold supplied by any certified D_V>=D0>0 block.

    Set mu0=D0/c.  The standard dyadic low-mass estimate then costs exactly
    c*mu0/2=D0/2, independent of the parent scale.
    """
    D0 = float(resolved_dissipation_lower)
    c = float(scaled_lifetime)
    if D0 <= 0 or c <= 0 or not math.isfinite(D0 + c):
        raise ValueError("positive finite dissipation lower and lifetime required")
    return D0 / c


def dissipation_supplier_retained_lower(
    total_resolved_dissipation: float,
    resolved_dissipation_lower: float,
    scaled_lifetime: float,
) -> float:
    """Actual D_V mass retained on shells with M||P_Mu||^2>=D0/c.

    This is a positive diagnostic restriction of resolved dissipation.  It is
    not a child-energy causal probability and must not be used for Shannon/Renyi.
    """
    D = float(total_resolved_dissipation)
    D0 = float(resolved_dissipation_lower)
    c = float(scaled_lifetime)
    if D < D0 or D0 <= 0 or c <= 0 or not all(math.isfinite(x) for x in (D, D0, c)):
        raise ValueError("require finite D_V>=D0>0 and positive lifetime")
    mu0 = dissipation_supplier_shell_mass_threshold(D0, c)
    return retained_dissipation_lower(D, c, mu0)


def dissipation_supplier_retained_fraction_lower(
    total_resolved_dissipation: float,
    resolved_dissipation_lower: float,
    scaled_lifetime: float,
) -> float:
    D = float(total_resolved_dissipation)
    return dissipation_supplier_retained_lower(D, resolved_dissipation_lower, scaled_lifetime) / D


def two_cell_cluster_to_whole_shell_mass_lower(two_cell_critical_mass_lower: float) -> float:
    """Convert M(E_C+E_C-r) lower to whole-shell M||P_Mu||_2^2 lower.

    Each coherent cell energy is bounded by the full shell energy, so
    E_C+E_C-r <= 2 E_shell.  No cell is promoted to a whole-shell identity.
    """
    m = float(two_cell_critical_mass_lower)
    if m < 0 or not math.isfinite(m):
        raise ValueError("finite nonnegative two-cell mass lower required")
    return 0.5 * m


def fresh_dominant_service_shell_mass_lower(
    square_service_threshold: float,
    dominant_fraction: float = DEFAULT_DOMINANT_FRESH_FRACTION,
) -> float:
    """Whole-shell supplier from the existing dominant fresh service branch.

    On the certified route: fresh service >=Y/4; a dominant edge carries at
    least theta of it; local capacity gives M(E_C+E_C-r)>=edge/2; the full shell
    has at least half that two-cell mass.  Therefore

        M E_shell >= theta Y/16,

    equal to Y/64 for theta=1/4.
    """
    Y = float(square_service_threshold)
    theta = float(dominant_fraction)
    if Y <= 0 or not math.isfinite(Y) or not (0 < theta < 1):
        raise ValueError("positive finite service threshold and 0<theta<1 required")
    return theta * Y / 16.0


def _finite_path(name: str, values: Sequence[float], n: int) -> tuple[float, ...]:
    x = tuple(float(v) for v in values)
    if len(x) != n or any(not math.isfinite(v) for v in x):
        raise ValueError(f"{name} must be a matching finite path")
    return x


def critical_shell_backward_first_hit(
    elapsed_times: Sequence[float],
    *,
    terminal_amplitude: float,
    strain_action: Sequence[float],
    residual_impulse_abs: Sequence[float],
    hh_impulse_abs: Sequence[float],
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """Three-monitor first hit for a shell carrier before materiality is assigned.

    The wrapper records the actually observed backward horizon.  A downstream
    outcome is not allowed to claim t=0 or a full natural survivor unless these
    monitors cover the whole required interval.
    """
    ell = np.asarray(elapsed_times, float)
    if ell.ndim != 1 or len(ell) < 2 or abs(float(ell[0])) > 1e-14 or np.any(np.diff(ell) <= 0):
        raise ValueError("elapsed times must start at zero and increase strictly")
    if np.any(~np.isfinite(ell)):
        raise ValueError("finite elapsed times required")
    amp = float(terminal_amplitude)
    if amp <= 0 or not math.isfinite(amp):
        raise ValueError("positive finite terminal amplitude required")
    n = len(ell)
    K = _finite_path("strain action", strain_action, n)
    IR = _finite_path("residual impulse", residual_impulse_abs, n)
    IH = _finite_path("HH impulse", hh_impulse_abs, n)
    if min(K) < 0 or min(IR) < 0 or min(IH) < 0:
        raise ValueError("nonnegative native monitor paths required")
    monitors = (
        PhysicalPathMonitor("high_strain_critical_dissipation", LOW_STRAIN_ACTION, K, ThresholdTopology.CLOSED),
        PhysicalPathMonitor(ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, RESIDUAL_FRACTION * amp, IR, ThresholdTopology.CLOSED),
        PhysicalPathMonitor(HH_COEFFICIENT_OBSTRUCTION, GENERATED_FRACTION * amp, IH, ThresholdTopology.CLOSED),
    )
    out = first_physical_corridor_exit(ell, monitors, tie_tolerance=tie_tolerance)
    needs_energy_reentry = any(
        label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for label in out.joint_first_stops
    )
    return {
        "first_elapsed": out.first_time,
        "joint_causes": out.joint_first_stops,
        "joint_first_stops": out.joint_first_stops,
        "individual_debuts": out.individual_debuts,
        "terminal_amplitude": amp,
        "observed_elapsed_end": float(ell[-1]),
        "requires_physical_energy_reentry": needs_energy_reentry,
        "coefficient_impulses_used_as_work": False,
    }


def critical_shell_natural_outcome(
    *,
    event_time: float,
    renewal_frequency: float,
    shell_critical_mass_lower: float,
    scaled_lifetime: float,
    viscosity: float,
    terminal_coefficient: complex,
    endpoint_coefficient: complex,
    hh_impulse: complex,
    residual_interface_impulse: complex,
    first_hit: dict[str, object],
) -> dict[str, object]:
    """Critical shell -> named recursive stop, t=0, or own-scale service.

    Material identity is intentionally absent before service.  A strain hit is
    the already named critical-dissipation recursion at the renewed scale;
    interface and HH coefficient hits only locate physical-energy reentry.  If no
    hit occurs, the monitor horizon must actually cover the entire backward
    interval required by the natural window or initial boundary before a survivor
    is certified.
    """
    geom = backward_natural_endpoint(event_time, renewal_frequency, scaled_lifetime)
    A = float(renewal_frequency)
    mu0 = float(shell_critical_mass_lower)
    c = float(scaled_lifetime)
    nu = float(viscosity)
    zt = complex(terminal_coefficient)
    zs = complex(endpoint_coefficient)
    ih = complex(hh_impulse)
    ir = complex(residual_interface_impulse)
    amp = abs(zt)
    if A <= 0 or mu0 <= 0 or c <= 0 or nu < 0 or amp <= 0:
        raise ValueError("positive renewal scale/shell mass/lifetime/coefficient and nonnegative viscosity required")
    if not all(math.isfinite(x) for x in (A, mu0, c, nu, amp)):
        raise ValueError("finite critical-shell outcome data required")
    terminal_lower = critical_shell_terminal_mass_lower(mu0)
    terminal_mass = A * amp * amp
    if terminal_mass < terminal_lower - 4e-12 * max(1.0, terminal_lower):
        raise ValueError("terminal coefficient does not realize the certified critical shell lower")
    monitor_amp = float(first_hit.get("terminal_amplitude", -1.0))
    amp_tol = 4e-12 * max(1.0, amp, abs(monitor_amp))
    if monitor_amp <= 0 or not math.isfinite(monitor_amp) or abs(monitor_amp - amp) > amp_tol:
        raise ValueError("first-hit monitor thresholds do not match the terminal coefficient amplitude")

    required = float(geom["elapsed_available"])
    horizon = float(first_hit.get("observed_elapsed_end", -1.0))
    if horizon + 2e-12 * max(1.0, required) < required:
        raise ValueError("first-hit monitors do not cover the required backward corridor")

    hit_time = first_hit.get("first_elapsed")
    causes = tuple(str(x) for x in first_hit.get("joint_first_stops", first_hit.get("joint_causes", ())))
    needs_energy_reentry = any(
        label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for label in causes
    )
    if hit_time is not None and float(hit_time) <= required + 2e-12 * max(1.0, required):
        return {
            "classification": "named_first_stop",
            "joint_causes": causes,
            "joint_first_stops": causes,
            "first_elapsed": float(hit_time),
            "required_elapsed": required,
            "observed_elapsed_end": horizon,
            "materiality_assigned": False,
            "primary_selected": False,
            "requires_physical_energy_reentry": needs_energy_reentry,
            "coefficient_impulses_used_as_work": False,
        }

    if bool(geom["hits_initial_boundary"]):
        return {
            "classification": "initial_boundary_root",
            "joint_causes": ("t=0",),
            "joint_first_stops": ("t=0",),
            "first_elapsed": required,
            "required_elapsed": required,
            "observed_elapsed_end": horizon,
            "materiality_assigned": False,
            "requires_physical_energy_reentry": False,
            "coefficient_impulses_used_as_work": False,
        }

    res = abs(exact_adjoint_residual(zt, zs, ih, ir))
    tol = 4e-12 * max(1.0, amp, abs(zs), abs(ih), abs(ir))
    if res > tol:
        raise ValueError("critical-shell Duhamel decomposition is not exact")
    if abs(ir) >= RESIDUAL_FRACTION * amp - tol:
        raise ValueError("endpoint interface impulse contradicts full no-hit corridor")
    if abs(ih) >= GENERATED_FRACTION * amp - tol:
        raise ValueError("endpoint HH impulse contradicts full no-hit corridor")
    inherited = abs(zs)
    if inherited < INHERIT_FRACTION * amp - tol:
        raise AssertionError("full natural critical-shell corridor lost the quarter coefficient")
    retained_mass = A * inherited * inherited
    clean_retained = critical_shell_survivor_coefficient_mass_lower(mu0)
    if retained_mass < clean_retained - 5e-12 * max(1.0, clean_retained):
        raise AssertionError("critical-shell survivor lost its clean coefficient mass")

    carrier = critical_shell_persistent_carrier_mass_lower(mu0, c, nu)
    Y = critical_shell_bounded_service_lower(mu0, c, nu)
    Sint = critical_shell_integrated_service_lower(mu0, c, nu)
    if min(carrier, Y, Sint) <= 0:
        raise AssertionError("positive critical shell failed to create positive own-scale service")
    return {
        "classification": "full_natural_own_scale_service",
        "joint_causes": (),
        "joint_first_stops": (),
        "required_elapsed": required,
        "observed_elapsed_end": horizon,
        "terminal_coefficient_mass": terminal_mass,
        "terminal_coefficient_mass_lower": terminal_lower,
        "retained_coefficient_mass": retained_mass,
        "clean_retained_coefficient_mass_lower": clean_retained,
        "uniform_carrier_mass_lower": carrier,
        "bounded_displacement_radius_over_A": BOUNDED_HEAT_RADIUS,
        "uniform_square_service_lower": Y,
        "integrated_bounded_heat_service_lower": Sint,
        "duhamel_residual": res,
        "materiality_assigned": "only_after_service_via_actual_Moyal_endpoints",
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def theorem_certificate(scaled_lifetime: float = 1.0, viscosity: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    nu = float(viscosity)
    if c <= 0 or nu < 0 or not math.isfinite(c + nu):
        raise ValueError("positive finite lifetime and nonnegative viscosity required")
    mu_hs = high_strain_ancestor_mass_threshold(c)
    terminal_hs = critical_shell_terminal_mass_lower(mu_hs)
    survivor_hs = critical_shell_survivor_coefficient_mass_lower(mu_hs)
    carrier_hs = critical_shell_persistent_carrier_mass_lower(mu_hs, c, nu)
    service_hs = critical_shell_bounded_service_lower(mu_hs, c, nu)
    existing_carrier = high_strain_persistent_carrier_mass_lower(c, nu)
    existing_service = high_strain_uniform_service_lower(c, nu)
    if not math.isclose(carrier_hs, existing_carrier, rel_tol=2e-14, abs_tol=2e-14):
        raise AssertionError("generic shell theorem did not specialize to the existing high-strain carrier lower")
    if not math.isclose(service_hs, existing_service, rel_tol=2e-14, abs_tol=2e-14):
        raise AssertionError("generic shell theorem did not specialize to the existing high-strain service lower")
    if persistent_seed_low_low_gap() <= 0:
        raise AssertionError("canonical shell envelope lost the low-low moat")
    return {
        "status": "EXACT_CRITICAL_SHELL_MASS_TO_OWN_SCALE_SERVICE_REENTRY__GENERIC_DV_AND_FRESH_CLUSTER_SUPPLIERS__MATERIALITY_DEFERRED",
        "local_input": "one actual shell-time event with M||P_Mu(t)||_2^2>=mu0>0; no packet, material label or probability law is part of the core theorem",
        "registration": "A=3M/4 and Q_A=1 on {M/2<|xi|<=M} give A|z(t)|^2>=(3/4)mu0 exactly with the shell's normalized state as terminal analysis probe",
        "first_stop": "before materiality is assigned, use renewed strain plus role-interface and HH coefficient obstructions; coefficient hits only locate physical-energy reentry, and the observed monitor horizon must cover the claimed natural/boundary interval",
        "survivor": "a full no-hit corridor keeps A|z(s)|^2>=(3/64)mu0 on every prefix; the registered analysis-dual cost J gives A||Q_Au(s)||^2>=(3mu0)/(64J^2)",
        "service": "the Arb-certified radius-3 annular heat truncation gives some actual |r|<=3/A with A||delta_r Q_Au||^2 >= (q_b 3mu0)/(64J^2), and normalized full-slab service at least c times this threshold",
        "material_order": "material OO/ON/NN is read only after the renewed positive service law exists; pure shell criticality carries no synthetic old/new label",
        "high_strain_specialization": f"mu0=32pi^2/(75c^2) gives terminal mass={terminal_hs:.12g}, quarter-survivor mass={survivor_hs:.12g}, and exactly the previously certified high-strain carrier/service constants",
        "dissipation_supplier": "any certified resolved D_V>=D0>0 supplies mu0=D0/c on at least D_V-D0/2>=D_V/2 actual dissipation mass; normalized restriction weights remain diagnostic, never HH causal probabilities",
        "fresh_service_supplier": "the existing dominant fresh service branch gives M(E_C+E_C-r)>=theta Y/8; because each cell energy is <=whole-shell energy, it supplies M E_shell>=theta Y/16=Y/64 at theta=1/4",
        "scale_scope": "this theorem is shell-local and proves own-scale service re-entry; supplier-specific signed-good progress relative to a previous block scale is a separate geometric statement and is not fabricated here",
        "scope": "this unifies high-strain critical shells, generic resolved-dissipation critical shells and dominant fresh coherent clusters at the carrier/service level; source/pressure routing, pure material-label transparency, and final continuum master assembly remain separate tasks",
    }


@dataclass(frozen=True)
class CriticalShellServiceStress:
    samples: int
    minimum_high_strain_specialization_margin: float
    minimum_dissipation_retained_fraction: float
    minimum_fresh_shell_mass_margin: float
    minimum_survivor_mass_margin: float
    minimum_service_threshold: float
    worst_duhamel_residual: float
    order_invariance_failures: int
    unit_invariance_failures: int
    horizon_guard_failures: int
    maximum_joint_first_stop_count: int
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260809) -> CriticalShellServiceStress:
    rng = np.random.default_rng(seed)
    mh = md = mf = ms = my = float("inf")
    wr = 0.0
    order_fail = unit_fail = horizon_fail = 0
    max_joint = 0
    branches: dict[str, int] = {}
    for _ in range(samples):
        c = float(math.exp(rng.uniform(math.log(0.3), math.log(2.5))))
        nu = float(rng.uniform(0.0, 2.0))
        mu0 = float(math.exp(rng.uniform(-10.0, 5.0)))

        # High-strain specialization identity is exact up to floating arithmetic.
        mu_hs = high_strain_ancestor_mass_threshold(c)
        spec = critical_shell_bounded_service_lower(mu_hs, c, nu)
        ref = high_strain_uniform_service_lower(c, nu)
        mh = min(mh, spec - ref)
        if abs(spec - ref) > 3e-13 * max(1.0, abs(ref)):
            raise AssertionError("generic critical-shell service lost high-strain specialization")

        # Generic resolved-D_V supplier: low-mass part costs exactly D0/2.
        D0 = float(math.exp(rng.uniform(-9.0, 4.0)))
        D = float(rng.uniform(1.0, 5.0)) * D0
        frac = dissipation_supplier_retained_fraction_lower(D, D0, c)
        md = min(md, frac)
        if frac < 0.5 - 3e-13:
            raise AssertionError("generic D0 supplier retained less than half actual D_V")

        # Existing dominant fresh branch -> two-cell mass -> whole-shell mass.
        Y = float(math.exp(rng.uniform(-10.0, 4.0)))
        theta = DEFAULT_DOMINANT_FRESH_FRACTION
        clean_pair = theta * Y / 8.0
        shell = two_cell_cluster_to_whole_shell_mass_lower(clean_pair)
        expected_shell = fresh_dominant_service_shell_mass_lower(Y, theta)
        mf = min(mf, shell - expected_shell)
        if shell + 2e-14 * max(1.0, Y) < expected_shell:
            raise AssertionError("fresh coherent cluster lost whole-shell supplier mass")

        # Generic three-monitor corridor.
        M = float(math.exp(rng.uniform(-2.0, 5.0)))
        A = renewal_scale(M)
        T = c / (A * A)
        terminal_lower = critical_shell_terminal_mass_lower(mu0)
        terminal_mass = float(rng.uniform(1.0, 2.0)) * terminal_lower
        amp = math.sqrt(terminal_mass / A)
        phase = float(rng.uniform(-math.pi, math.pi))
        zt = amp * complex(math.cos(phase), math.sin(phase))
        event_time = float(rng.uniform(1.1, 2.5)) * T
        mode = int(rng.integers(0, 6))
        if mode == 4:
            event_time = float(rng.uniform(0.2, 0.9)) * T
        required = min(T, event_time)
        ell = np.linspace(0.0, required, 6)
        Kend = float(rng.uniform(0.05, 0.70)) * LOW_STRAIN_ACTION
        IRend = float(rng.uniform(0.05, 0.65)) * RESIDUAL_FRACTION * amp
        IHend = float(rng.uniform(0.05, 0.65)) * GENERATED_FRACTION * amp
        if mode == 0:
            Kend = float(rng.uniform(1.05, 1.7)) * LOW_STRAIN_ACTION
        elif mode == 1:
            IRend = float(rng.uniform(1.05, 1.7)) * RESIDUAL_FRACTION * amp
        elif mode == 2:
            IHend = float(rng.uniform(1.05, 1.7)) * GENERATED_FRACTION * amp
        elif mode == 5:
            # Exact two-face tie at the final observed time.
            Kend = LOW_STRAIN_ACTION
            IRend = RESIDUAL_FRACTION * amp
        Kpath = np.linspace(0.0, Kend, 6)
        IRpath = np.linspace(0.0, IRend, 6)
        IHpath = np.linspace(0.0, IHend, 6)
        time_tol = 2e-10 * max(1.0, abs(required))
        hit = critical_shell_backward_first_hit(
            ell,
            terminal_amplitude=amp,
            strain_action=Kpath,
            residual_impulse_abs=IRpath,
            hh_impulse_abs=IHpath,
            tie_tolerance=time_tol,
        )
        max_joint = max(max_joint, len(tuple(hit["joint_first_stops"])))

        mons = [
            PhysicalPathMonitor("high_strain_critical_dissipation", LOW_STRAIN_ACTION, tuple(Kpath), ThresholdTopology.CLOSED),
            PhysicalPathMonitor(ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, RESIDUAL_FRACTION * amp, tuple(IRpath), ThresholdTopology.CLOSED),
            PhysicalPathMonitor(HH_COEFFICIENT_OBSTRUCTION, GENERATED_FRACTION * amp, tuple(IHpath), ThresholdTopology.CLOSED),
        ]
        base = first_physical_corridor_exit(ell, mons, tie_tolerance=time_tol)
        perm = rng.permutation(3)
        alt = first_physical_corridor_exit(ell, [mons[int(i)] for i in perm], tie_tolerance=time_tol)
        if base.first_time != alt.first_time or set(base.joint_first_stops) != set(alt.joint_first_stops):
            order_fail += 1
            raise AssertionError("generic critical-shell first stop depended on monitor order")
        scaled = [rescale_monitor_units(m, float(math.exp(rng.uniform(-8.0, 8.0)))) for m in mons]
        altu = first_physical_corridor_exit(ell, scaled, tie_tolerance=time_tol)
        if (base.first_time is None) != (altu.first_time is None):
            unit_fail += 1
            raise AssertionError("generic critical-shell first stop depended on monitor units")
        if base.first_time is not None and (
            abs(float(base.first_time) - float(altu.first_time)) > time_tol
            or set(base.joint_first_stops) != set(altu.joint_first_stops)
        ):
            unit_fail += 1
            raise AssertionError("generic critical-shell first stop depended on monitor units")

        rp = float(rng.uniform(-math.pi, math.pi))
        hp = float(rng.uniform(-math.pi, math.pi))
        ir = IRend * complex(math.cos(rp), math.sin(rp))
        ih = IHend * complex(math.cos(hp), math.sin(hp))
        zs = zt - ir - ih
        out = critical_shell_natural_outcome(
            event_time=event_time,
            renewal_frequency=A,
            shell_critical_mass_lower=mu0,
            scaled_lifetime=c,
            viscosity=nu,
            terminal_coefficient=zt,
            endpoint_coefficient=zs,
            hh_impulse=ih,
            residual_interface_impulse=ir,
            first_hit=hit,
        )
        b = str(out["classification"])
        branches[b] = branches.get(b, 0) + 1
        if b == "full_natural_own_scale_service":
            clean = critical_shell_survivor_coefficient_mass_lower(mu0)
            ms = min(ms, float(out["retained_coefficient_mass"]) - clean)
            my = min(my, float(out["uniform_square_service_lower"]))
            wr = max(wr, float(out["duhamel_residual"]))

        # Horizon guard: a no-hit claim with only half the required history must fail.
        if mode == 3 and required > 0:
            short_ell = np.linspace(0.0, 0.5 * required, 4)
            short_hit = critical_shell_backward_first_hit(
                short_ell,
                terminal_amplitude=amp,
                strain_action=np.linspace(0.0, 0.2 * LOW_STRAIN_ACTION, 4),
                residual_impulse_abs=np.linspace(0.0, 0.1 * RESIDUAL_FRACTION * amp, 4),
                hh_impulse_abs=np.linspace(0.0, 0.1 * GENERATED_FRACTION * amp, 4),
            )
            try:
                critical_shell_natural_outcome(
                    event_time=event_time,
                    renewal_frequency=A,
                    shell_critical_mass_lower=mu0,
                    scaled_lifetime=c,
                    viscosity=nu,
                    terminal_coefficient=zt,
                    endpoint_coefficient=zt,
                    hh_impulse=0j,
                    residual_interface_impulse=0j,
                    first_hit=short_hit,
                )
            except ValueError as exc:
                if "do not cover" not in str(exc):
                    raise
            else:
                horizon_fail += 1
                raise AssertionError("incomplete critical-shell corridor was certified as complete")

    if not math.isfinite(ms):
        ms = 0.0
    if not math.isfinite(my):
        my = critical_shell_bounded_service_lower(1.0, 1.0, 1.0)
    return CriticalShellServiceStress(
        samples,
        mh,
        md,
        mf,
        ms,
        my,
        wr,
        order_fail,
        unit_fail,
        horizon_fail,
        max_joint,
        branches,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-critical-shell-service-reentry"))
    ap.add_argument("--scaled-lifetime", type=float, default=1.0)
    ap.add_argument("--viscosity", type=float, default=1.0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate(args.scaled_lifetime, args.viscosity)
    out = stress(args.samples)
    (args.outdir / "critical_shell_service_reentry.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    mu_hs = high_strain_ancestor_mass_threshold(args.scaled_lifetime)
    md = f"""# Generic critical shell -> own-scale coherent-service re-entry

Status: **{cert['status']}**.

The core input is deliberately local and deterministic.  At one physical shell-time event assume only

`M ||P_M u(t)||_2^2 >= mu0 > 0`.

Set `A=3M/4`.  Choose the same smooth scalar role `Q_A` which equals one on the hard shell and whose support remains separated from the `A/2` low--low output on the low-strain branch.  With the shell's normalized state as terminal analysis probe,

`A|z(t)|^2 >= (3/4)mu0`.

Inspect backward through the required `A`-natural interval using exactly three native monitors: renewed strain, role-interface coefficient obstruction, and HH-regeneration coefficient obstruction.  The two coefficient monitors are interval locators: at a hit, the same smooth carrier must reenter the physical-energy gate before any work owner is named.  Materiality is not assigned yet.  Exact ties stay unsplit.  The first-hit record also stores the actually observed backward horizon; neither a `t=0` root nor a full-natural survivor may be certified from a shorter monitor path.

If renewed strain hits, it keeps its critical-dissipation owner.  If a coefficient monitor hits, it only locates reentry of the same carrier into the physical-energy gate.  If the required interval reaches `t=0`, the initial boundary absorbs it.  Otherwise every prefix remains below both coefficient-obstruction faces and

`|z(s)|>|z(t)|/4`.

Hence throughout the full no-hit slab

`A|z(s)|^2 >= 3mu0/64`.

For the already registered affine/Kelvin/viscous analysis dual,

`||psi(s)|| <= J||psi(t)||`,

so

`A||Q_Au(s)||_2^2 >= 3mu0/(64J^2)`.

The annular heat multiplier and the Arb-certified radius-3 Gaussian truncation then force an actual bounded displacement

`|r_s|<=3/A`,

`A||delta_(r_s)Q_Au(s)||_2^2 >= Y_shell(mu0)`,

where

`Y_shell(mu0)=q_b 3mu0/(64J^2)>0`.

The normalized full-natural-slab bounded heat service is at least `c Y_shell(mu0)`.  Only after this positive renewed service exists do we apply exact Moyal and read material OO/ON/NN from its actual two endpoints.

Two physically different suppliers now enter the same local theorem without being conflated with causal probability:

1. **Resolved dissipation.**  If a block already has a certified `D_V>=D0>0`, choose `mu0=D0/c`.  The standard dyadic low-mass part costs exactly `D0/2`, so at least half the actual `D_V` mass lies on qualifying shells.  Normalizing that restriction is diagnostic sampling only; Shannon/Renyi still use actual positive HH child-energy work.
2. **Dominant fresh coherent service.**  The existing fresh-service route gives `M(E_C+E_(C-r))>=theta Y/8`.  Since each coherent cell energy is bounded by the whole shell energy, this supplies `M E_shell>=theta Y/16`, i.e. `Y/64` at `theta=1/4`, with no packet persistence claim.

The high-strain theorem is an exact specialization.  With

`mu0=32pi^2/(75c^2)`,

the generic terminal mass is `8pi^2/(25c^2)`, the quarter-survivor mass is `pi^2/(50c^2)`, and the carrier/service constants coincide with the already certified high-strain annular-service theorem.  At the requested `(c,nu)=({args.scaled_lifetime:g},{args.viscosity:g})`, the specialized bounded-service threshold is `{critical_shell_bounded_service_lower(mu_hs,args.scaled_lifetime,args.viscosity):.12g}`.

Stress: `{out.samples}` generic shell/supplier/corridor/service states
- minimum high-strain specialization signed margin: `{out.minimum_high_strain_specialization_margin:.3e}`
- minimum generic D0 retained fraction: `{out.minimum_dissipation_retained_fraction:.9f}`
- minimum fresh-cluster whole-shell identity margin: `{out.minimum_fresh_shell_mass_margin:.3e}`
- minimum full-survivor coefficient-mass margin: `{out.minimum_survivor_mass_margin:.3e}`
- minimum sampled own-scale service threshold: `{out.minimum_service_threshold:.3e}`
- worst exact survivor Duhamel residual: `{out.worst_duhamel_residual:.3e}`
- monitor-order failures: `{out.order_invariance_failures}`
- monitor-unit failures: `{out.unit_invariance_failures}`
- incomplete-horizon certification failures: `{out.horizon_guard_failures}`
- maximum sampled exact joint first-stop count: `{out.maximum_joint_first_stop_count}`
- outcomes: `{out.branch_counts}`

This theorem is intentionally **shell-local**.  It proves carrier/service re-entry but does not manufacture signed-good scale progress relative to whichever previous block supplied the shell.  Pressure/source routing and pure material-label transparency remain separate continuum tasks.  No global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
