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
    integrated_bounded_heat_service_lower,
    persistent_carrier_critical_mass_lower,
    uniform_bounded_square_service_lower,
)
from src.high_strain_resolved_ancestor import (
    high_strain_ancestor_mass_threshold,
    retained_fraction_lower,
)
from src.nn_critical_heat_carrier_seed import (
    LOW_STRAIN_ACTION,
    RENEWAL_SCALE_FACTOR,
    persistent_seed_low_low_gap,
    renewal_carrier_critical_mass_lower,
    renewal_natural_lifetime_ratio,
    renewal_scale,
)
from src.nn_seed_temporal_first_stop import backward_natural_endpoint, renewed_natural_duration
from src.smooth_sgs_first_hit_extraction import (
    PhysicalPathMonitor,
    ThresholdTopology,
    first_physical_corridor_exit,
    rescale_monitor_units,
)


@dataclass(frozen=True)
class CriticalDissipationAtom:
    """One shell-time atom of the retained positive D_V|_G law.

    ``mass`` is actual positive normalized resolved dissipation mass.  The law is
    already restricted to G, so the shell-time mark independently satisfies
    M||P_j u(t)||_2^2>=mu_*.  No coherent/material mark is assumed here.
    """

    mass: float
    child_frequency: float
    shell_upper_frequency: float
    shell_energy_u: float
    time: float


@dataclass(frozen=True)
class CriticalCarrierSeed:
    normalized_dissipation_weight: float
    dissipation_mass: float
    time: float
    scaled_lifetime: float
    child_frequency: float
    shell_upper_frequency: float
    renewal_frequency: float
    shell_critical_mass: float
    renewal_critical_mass: float
    natural_lifetime_ratio: float


def pushforward_critical_dissipation_law(
    atoms: Sequence[CriticalDissipationAtom],
    *,
    scaled_lifetime: float,
) -> tuple[CriticalCarrierSeed, ...]:
    """Push actual positive D_V|_G mass to whole-shell smooth carrier seeds.

    This is the same shell registration geometry previously used after the
    NN-intersect-critical heat law, but its mathematical input is only the
    critical shell-time mark.  The positive law is now D_V restricted to G,
    which the high-strain ancestor theorem already proves has at least half of
    total D_V mass.  Material ownership is intentionally absent at this stage.
    """
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    rows = tuple(atoms)
    if not rows:
        raise ValueError("nonempty retained critical dissipation law required")
    mu_star = high_strain_ancestor_mass_threshold(c)
    total = 0.0
    for a in rows:
        vals = (a.mass, a.child_frequency, a.shell_upper_frequency, a.shell_energy_u, a.time)
        if a.mass <= 0 or a.child_frequency <= 0 or a.shell_upper_frequency <= 0 or a.shell_energy_u < 0 or a.time <= 0:
            raise ValueError("positive mass/frequencies/time and nonnegative shell energy required")
        if not all(math.isfinite(x) for x in vals):
            raise ValueError("finite critical dissipation atom data required")
        if a.shell_upper_frequency > a.child_frequency / 4.0 * (1.0 + 1e-13):
            raise ValueError("critical ancestor shell must satisfy M<=N/4")
        shell_mass = a.shell_upper_frequency * a.shell_energy_u
        shell_tol = 3e-13 * max(shell_mass, mu_star)
        if shell_mass + shell_tol < mu_star:
            raise ValueError("dissipation atom is not on the critical shell-time set G")
        total += a.mass
    seeds: list[CriticalCarrierSeed] = []
    for a in rows:
        M = a.shell_upper_frequency
        A = renewal_scale(M)
        seeds.append(
            CriticalCarrierSeed(
                normalized_dissipation_weight=a.mass / total,
                dissipation_mass=a.mass,
                time=a.time,
                scaled_lifetime=c,
                child_frequency=a.child_frequency,
                shell_upper_frequency=M,
                renewal_frequency=A,
                shell_critical_mass=M * a.shell_energy_u,
                renewal_critical_mass=A * a.shell_energy_u,
                natural_lifetime_ratio=renewal_natural_lifetime_ratio(a.child_frequency, M),
            )
        )
    if not math.isclose(sum(x.normalized_dissipation_weight for x in seeds), 1.0, rel_tol=2e-14, abs_tol=2e-14):
        raise AssertionError("critical dissipation-law pushforward failed to normalize")
    lower = renewal_carrier_critical_mass_lower(c)
    for x in seeds:
        renewed_tol = 4e-13 * max(x.renewal_critical_mass, lower)
        if x.renewal_critical_mass + renewed_tol < lower:
            raise AssertionError("critical D_V atom lost renewed whole-shell mass")
    return tuple(seeds)


def _finite_path(name: str, values: Sequence[float], n: int) -> tuple[float, ...]:
    x = tuple(float(v) for v in values)
    if len(x) != n or any(not math.isfinite(v) for v in x):
        raise ValueError(f"{name} must be a matching finite path")
    return x


def critical_seed_backward_first_hit(
    elapsed_times: Sequence[float],
    *,
    terminal_amplitude: float,
    strain_action: Sequence[float],
    residual_impulse_abs: Sequence[float],
    hh_impulse_abs: Sequence[float],
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """First obstruction for a critical shell seed before materiality is assigned.

    There are exactly three native monitors: renewed strain, role-interface
    coefficient obstruction, and HH coefficient obstruction.
    No artificial material-distance observable is inserted because this route has
    not yet assigned the carrier to old/new material.  Materiality will be read
    only after an actual renewed positive service law exists.
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
        raise ValueError("nonnegative physical monitor paths required")
    monitors = (
        PhysicalPathMonitor("high_strain_critical_dissipation", float(LOW_STRAIN_ACTION), K, ThresholdTopology.CLOSED),
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


def critical_seed_natural_outcome(
    *,
    source_seed: CriticalCarrierSeed,
    event_time: float,
    renewal_frequency: float,
    scaled_lifetime: float,
    viscosity: float,
    terminal_coefficient: complex,
    endpoint_coefficient: complex,
    hh_impulse: complex,
    residual_interface_impulse: complex,
    first_hit: dict[str, object],
) -> dict[str, object]:
    """High-strain critical seed -> named stop, t=0, or own-scale service.

    The service conclusion is available only for a full natural no-hit corridor.
    Large interface/HH coefficient impulses retain their obstruction provenance,
    locate physical-energy reentry, and are never promoted to work.  A renewed
    high-strain hit is another critical dissipation recursion, not an additive
    reset.
    """
    seed = source_seed
    t = float(event_time)
    A = float(renewal_frequency)
    zt = complex(terminal_coefficient)
    zs = complex(endpoint_coefficient)
    ih = complex(hh_impulse)
    ir = complex(residual_interface_impulse)
    amp = abs(zt)
    c = float(scaled_lifetime)
    nu = float(viscosity)
    source_values = (
        seed.normalized_dissipation_weight,
        seed.dissipation_mass,
        seed.time,
        seed.scaled_lifetime,
        seed.child_frequency,
        seed.shell_upper_frequency,
        seed.renewal_frequency,
        seed.shell_critical_mass,
        seed.renewal_critical_mass,
        seed.natural_lifetime_ratio,
    )
    if (
        A <= 0
        or amp <= 0
        or c <= 0
        or nu < 0
        or t <= 0
        or not all(
            math.isfinite(x)
            for x in (
                A,
                amp,
                c,
                nu,
                t,
                zt.real,
                zt.imag,
                zs.real,
                zs.imag,
                ih.real,
                ih.imag,
                ir.real,
                ir.imag,
            )
        )
        or not all(
        math.isfinite(x) and x > 0 for x in source_values
        )
    ):
        raise ValueError("positive finite source seed, event/lifetime and nonnegative viscosity required")
    expected_A = renewal_scale(seed.shell_upper_frequency)
    expected_renewed_mass = RENEWAL_SCALE_FACTOR * seed.shell_critical_mass
    expected_lifetime_ratio = renewal_natural_lifetime_ratio(seed.child_frequency, seed.shell_upper_frequency)
    provenance = (
        (t, seed.time, "event time"),
        (A, seed.renewal_frequency, "renewal frequency"),
        (A, expected_A, "parent-shell renewal scale"),
        (c, seed.scaled_lifetime, "scaled lifetime"),
        (seed.renewal_critical_mass, expected_renewed_mass, "renewed critical mass"),
        (seed.natural_lifetime_ratio, expected_lifetime_ratio, "natural lifetime ratio"),
    )
    for supplied, carried, name in provenance:
        tol = 4e-12 * max(abs(supplied), abs(carried))
        if abs(supplied - carried) > tol:
            raise ValueError(f"source seed {name} provenance does not match the requested corridor")
    if seed.shell_upper_frequency > seed.child_frequency / 4.0 * (1.0 + 1e-13):
        raise ValueError("source seed parent shell escaped its certified child-scale range")
    geom = backward_natural_endpoint(t, A, c)
    seed_lower = renewal_carrier_critical_mass_lower(c)
    terminal_mass = A * amp * amp
    terminal_tol = 4e-12 * max(terminal_mass, seed.renewal_critical_mass)
    if abs(terminal_mass - seed.renewal_critical_mass) > terminal_tol:
        raise ValueError("terminal coefficient does not realize the carried source seed mass")
    lower_tol = 4e-12 * max(terminal_mass, seed_lower)
    if terminal_mass + lower_tol < seed_lower:
        raise ValueError("terminal coefficient is not a certified critical shell seed")
    monitor_amp = float(first_hit.get("terminal_amplitude", -1.0))
    monitor_amp_tol = 4e-12 * max(amp, abs(monitor_amp))
    if monitor_amp <= 0 or not math.isfinite(monitor_amp) or abs(monitor_amp - amp) > monitor_amp_tol:
        raise ValueError("first-hit monitor thresholds do not match the terminal coefficient amplitude")
    res = abs(exact_adjoint_residual(zt, zs, ih, ir))
    duhamel_tol = 4e-12 * max(amp, abs(zs), abs(ih), abs(ir))
    if res > duhamel_tol:
        raise ValueError("critical-seed Duhamel decomposition is not exact")

    hit_time = first_hit.get("first_elapsed")
    causes = tuple(str(x) for x in first_hit.get("joint_first_stops", first_hit.get("joint_causes", ())))
    needs_energy_reentry = any(
        label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for label in causes
    )
    elapsed = float(geom["elapsed_available"])
    horizon = float(first_hit.get("observed_elapsed_end", -1.0))
    horizon_tol = 2e-12 * max(abs(horizon), elapsed)
    if horizon <= 0 or horizon + horizon_tol < elapsed:
        raise ValueError("first-hit monitors do not cover the required backward corridor")
    hit = None if hit_time is None else float(hit_time)
    allowed_causes = {
        "high_strain_critical_dissipation",
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
        HH_COEFFICIENT_OBSTRUCTION,
    }
    if any(label not in allowed_causes for label in causes):
        raise ValueError("first-hit certificate contains an unknown physical monitor")
    if hit is None and causes:
        raise ValueError("first-hit causes require an actual finite debut time")
    if hit is not None and (not math.isfinite(hit) or hit < 0 or not causes):
        raise ValueError("first-hit debut and joint cause set are inconsistent")
    if hit is not None and hit > horizon + 2e-12 * max(abs(hit), horizon):
        raise ValueError("first-hit debut lies beyond the observed monitor horizon")
    hit_tol = 0.0 if hit is None else 2e-12 * max(abs(hit), elapsed)
    if hit is not None and hit <= elapsed + hit_tol:
        return {
            "classification": "named_first_stop",
            "joint_causes": causes,
            "joint_first_stops": causes,
            "first_elapsed": hit,
            "primary_selected": False,
            "duhamel_residual": res,
            "materiality_assigned": False,
            "requires_physical_energy_reentry": needs_energy_reentry,
            "coefficient_impulses_used_as_work": False,
        }
    if bool(geom["hits_initial_boundary"]):
        return {
            "classification": "initial_boundary_root",
            "joint_causes": ("t=0",),
            "joint_first_stops": ("t=0",),
            "first_elapsed": elapsed,
            "duhamel_residual": res,
            "materiality_assigned": False,
            "requires_physical_energy_reentry": False,
            "coefficient_impulses_used_as_work": False,
        }
    coefficient_tol = 4e-12 * amp
    if abs(ir) >= RESIDUAL_FRACTION * amp - coefficient_tol:
        raise ValueError("endpoint interface impulse contradicts no-hit corridor")
    if abs(ih) >= GENERATED_FRACTION * amp - coefficient_tol:
        raise ValueError("endpoint HH impulse contradicts no-hit corridor")
    inherited = abs(zs)
    clean = INHERIT_FRACTION * amp
    if inherited < clean - coefficient_tol:
        raise AssertionError("critical seed full natural corridor lost quarter coefficient")
    retained_mass = A * inherited * inherited
    clean_retained = INHERIT_FRACTION**2 * seed_lower
    retained_tol = 5e-12 * max(retained_mass, clean_retained)
    if retained_mass + retained_tol < clean_retained:
        raise AssertionError("critical seed survivor lost clean critical coefficient mass")
    carrier = persistent_carrier_critical_mass_lower(c, nu)
    Y0 = uniform_bounded_square_service_lower(c, nu)
    Sint = integrated_bounded_heat_service_lower(c, nu)
    if Y0 <= 0 or Sint <= 0:
        raise AssertionError("full natural critical carrier did not create positive own-scale service")
    return {
        "classification": "full_natural_own_scale_service",
        "joint_causes": (),
        "joint_first_stops": (),
        "natural_duration": float(geom["natural_duration"]),
        "backward_endpoint": float(geom["backward_endpoint"]),
        "required_elapsed": elapsed,
        "observed_elapsed_end": horizon,
        "corridor_terminal_time": t,
        "corridor_endpoint_time": float(geom["backward_endpoint"]),
        "corridor_endpoint_elapsed_from_terminal": float(geom["elapsed_available"]),
        "physical_time_drop": float(geom["natural_duration"]),
        "renewal_frequency": A,
        "scaled_lifetime": c,
        "parent_shell_frequency": seed.shell_upper_frequency,
        "source_child_frequency": seed.child_frequency,
        "source_dissipation_mass": seed.dissipation_mass,
        "source_event_time": seed.time,
        "service_same_corridor_witness": True,
        "service_adds_recursion_depth": False,
        "terminal_critical_mass": terminal_mass,
        "retained_coefficient_critical_mass": retained_mass,
        "clean_retained_coefficient_mass_lower": clean_retained,
        "endpoint_carrier_critical_mass_lower": carrier,
        "bounded_displacement_radius_over_A": BOUNDED_HEAT_RADIUS,
        "uniform_square_service_lower": Y0,
        "integrated_bounded_heat_service_lower": Sint,
        "duhamel_residual": res,
        "materiality_assigned": "only_after_service_via_exact_Moyal_OO_ON_NN",
        "nn_seed_required": False,
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def theorem_certificate(scaled_lifetime: float = 1.0, viscosity: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    nu = float(viscosity)
    if c <= 0 or nu < 0 or not math.isfinite(c + nu):
        raise ValueError("positive finite lifetime and nonnegative viscosity required")
    mu = high_strain_ancestor_mass_threshold(c)
    seed = renewal_carrier_critical_mass_lower(c)
    retained = INHERIT_FRACTION**2 * seed
    gap = persistent_seed_low_low_gap()
    Y0 = uniform_bounded_square_service_lower(c, nu)
    Sint = integrated_bounded_heat_service_lower(c, nu)
    # On the high-strain branch D>=D_*=c*mu, the existing theorem retains >=1/2.
    frac = retained_fraction_lower(c * mu, c, mu)
    if frac < 0.5 - 2e-14:
        raise AssertionError("critical dissipation law lost its clean half")
    if gap <= 0 or Y0 <= 0:
        raise AssertionError("renewed role lost moat or service threshold")
    return {
        "status": "EXACT_HIGH_STRAIN_CRITICAL_DISSIPATION_LAW_TO_OWN_SCALE_SERVICE_REENTRY__NN_NOT_REQUIRED_FOR_RENEWAL_ENTRANCE__UNIVERSAL_SOURCE_RELINK_REMAINS",
        "positive_input_law": f"high strain gives D_V|_G >= D_V/2 on the actual positive resolved-dissipation law; clean endpoint fraction={frac:.12g}",
        "pushforward": "normalize D_V restricted to G and push its deterministic shell-time mark (j,t) to A=3M/4; the typed seed carries t, c, child N, parent M, A and actual renewed mass, all of which the corridor must reuse; no heat edge, coherent-cell argmax, or material label is required for this entrance",
        "critical_seed": f"M||P_j u||^2>=mu_*={mu:.12g} gives A||P_j u||^2>=(3/4)mu_*={seed:.12g}; Q_A=1 on the whole shell registers that coefficient exactly",
        "first_stop": "backward over one A-natural window use renewed strain plus role-interface and HH coefficient obstructions; exact ties remain unsplit and coefficient hits only locate physical-energy reentry",
        "corridor": f"if no monitor hits and t=0 is not reached, |z(s)|>=|z(t)|/4 and the clean retained coefficient mass is >={retained:.12g}",
        "service": f"the completed full-natural corridor carries its own bounded A-scale increment service: some |r|<=3/A has A||delta_r Q_Au||^2>=Y0={Y0:.12g}, and integrated normalized bounded heat service is >={Sint:.12g}; this is a same-corridor witness and adds no second recursion edge",
        "material_order": "material ownership is deliberately deferred: after the renewed positive service exists, exact Moyal assigns its actual endpoints and OO/ON/NN is read from that new law; no child-scale NN witness is propagated as whole-carrier ownership",
        "architectural_shortcut": "the shortest high-strain renewal entrance no longer requires child-scale heat ownership, old-incident erosion, or NN-intersect-critical selection before a carrier can renew; those theorems remain valid refinements for material capacity/provenance",
        "causal_scope": "normalized D_V|_G weights are only a positive diagnostic sampling law for the high-strain recursive route; they never replace actual positive HH child-energy work in Shannon/Renyi or transfer causality",
        "currency": "a renewed high-strain first hit recursively re-enters critical dissipation; interface/HH coefficient hits only locate physical-energy reentry, where actual gauge-quotiented work receives its native owner; none is promoted to a scale-independent additive reset",
        "scope": "this closes the high-strain route to either an already named recursive stop, the absorbing initial boundary, or a completed full-natural corridor carrying its own coherent-service witness; the service witness adds no event depth. Universal source/SGS and genuine relink slab renewal, and global master closure, remain open",
    }


@dataclass(frozen=True)
class HighStrainCriticalReentryStress:
    samples: int
    minimum_seed_mass_margin: float
    minimum_normalized_weight: float
    minimum_survivor_mass_margin: float
    minimum_service_threshold: float
    worst_duhamel_residual: float
    order_invariance_failures: int
    unit_invariance_failures: int
    maximum_joint_first_stop_count: int
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260809) -> HighStrainCriticalReentryStress:
    rng = np.random.default_rng(seed)
    ms = mp = mr = my = float("inf")
    wr = 0.0
    order_fail = unit_fail = 0
    max_joint = 0
    branches: dict[str, int] = {}
    for _ in range(samples):
        c = float(math.exp(rng.uniform(math.log(0.3), math.log(2.5))))
        nu = float(rng.uniform(0.0, 2.0))
        N = float(math.exp(rng.uniform(-1.0, 7.0)))
        mu = high_strain_ancestor_mass_threshold(c)
        n_atoms = int(rng.integers(1, 10))
        atoms: list[CriticalDissipationAtom] = []
        for _a in range(n_atoms):
            j = int(rng.integers(0, 8))
            M = (N / 4.0) * (2.0 ** (-j))
            A_atom = renewal_scale(M)
            T_atom = renewed_natural_duration(A_atom, c)
            Eu = float(rng.uniform(1.0, 4.0)) * mu / M
            atoms.append(
                CriticalDissipationAtom(
                    mass=float(math.exp(rng.uniform(-7.0, 4.0))),
                    child_frequency=N,
                    shell_upper_frequency=M,
                    shell_energy_u=Eu,
                    time=float(rng.uniform(0.15, 2.5)) * T_atom,
                )
            )
        seeds = pushforward_critical_dissipation_law(atoms, scaled_lifetime=c)
        lower = renewal_carrier_critical_mass_lower(c)
        ms = min(ms, min(x.renewal_critical_mass for x in seeds) - lower)
        mp = min(mp, min(x.normalized_dissipation_weight for x in seeds))
        if abs(sum(x.normalized_dissipation_weight for x in seeds) - 1.0) > 3e-13:
            raise AssertionError("critical D_V seed law lost normalized mass")
        chosen = seeds[int(rng.integers(0, len(seeds)))]
        A = chosen.renewal_frequency
        T = renewed_natural_duration(A, c)
        amp = math.sqrt(chosen.renewal_critical_mass / A)
        theta = float(rng.uniform(-math.pi, math.pi))
        zt = amp * complex(math.cos(theta), math.sin(theta))
        ell = np.linspace(0.0, T, 5)
        mode = int(rng.integers(0, 5))
        Kend = float(rng.uniform(0.10, 0.80)) * LOW_STRAIN_ACTION
        IRend = float(rng.uniform(0.05, 0.70)) * RESIDUAL_FRACTION * amp
        IHend = float(rng.uniform(0.05, 0.70)) * GENERATED_FRACTION * amp
        if mode == 0:
            Kend = float(rng.uniform(1.02, 1.8)) * LOW_STRAIN_ACTION
        elif mode == 1:
            IRend = float(rng.uniform(1.02, 1.8)) * RESIDUAL_FRACTION * amp
        elif mode == 2:
            IHend = float(rng.uniform(1.02, 1.8)) * GENERATED_FRACTION * amp
        Kpath = np.linspace(0.0, Kend, 5)
        IRpath = np.linspace(0.0, IRend, 5)
        IHpath = np.linspace(0.0, IHend, 5)
        hit = critical_seed_backward_first_hit(
            ell,
            terminal_amplitude=amp,
            strain_action=Kpath,
            residual_impulse_abs=IRpath,
            hh_impulse_abs=IHpath,
        )
        max_joint = max(max_joint, len(tuple(hit["joint_first_stops"])))

        # Native-unit/order invariance of the first-stop set.
        mons = [
            PhysicalPathMonitor("high_strain_critical_dissipation", LOW_STRAIN_ACTION, tuple(Kpath), ThresholdTopology.CLOSED),
            PhysicalPathMonitor(ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, RESIDUAL_FRACTION * amp, tuple(IRpath), ThresholdTopology.CLOSED),
            PhysicalPathMonitor(HH_COEFFICIENT_OBSTRUCTION, GENERATED_FRACTION * amp, tuple(IHpath), ThresholdTopology.CLOSED),
        ]
        time_tol = 2e-10 * abs(T)
        base = first_physical_corridor_exit(ell, mons, tie_tolerance=time_tol)
        perm = rng.permutation(3)
        alt = first_physical_corridor_exit(ell, [mons[int(i)] for i in perm], tie_tolerance=time_tol)
        if base.joint_first_stops != alt.joint_first_stops or base.first_time != alt.first_time:
            order_fail += 1
            raise AssertionError("critical-seed first stop depended on monitor order")
        scaled = [rescale_monitor_units(m, float(math.exp(rng.uniform(-8.0, 8.0)))) for m in mons]
        altu = first_physical_corridor_exit(ell, scaled, tie_tolerance=time_tol)
        if (base.first_time is None) != (altu.first_time is None):
            unit_fail += 1
            raise AssertionError("critical-seed first stop depended on monitor units")
        if base.first_time is not None and (
            abs(float(base.first_time) - float(altu.first_time)) > time_tol
            or set(base.joint_first_stops) != set(altu.joint_first_stops)
        ):
            unit_fail += 1
            raise AssertionError("critical-seed first stop depended on monitor units")

        # Endpoint impulses are the actual cumulative complex impulses.  For stop
        # modes they may cross thresholds; for no-hit modes they stay safely inside.
        rp = float(rng.uniform(-math.pi, math.pi))
        hp = float(rng.uniform(-math.pi, math.pi))
        ir = IRend * complex(math.cos(rp), math.sin(rp))
        ih = IHend * complex(math.cos(hp), math.sin(hp))
        zs = zt - ir - ih
        out = critical_seed_natural_outcome(
            source_seed=chosen,
            event_time=chosen.time,
            renewal_frequency=A,
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
        wr = max(wr, float(out["duhamel_residual"]))
        if b == "full_natural_own_scale_service":
            clean = INHERIT_FRACTION**2 * lower
            mr = min(mr, float(out["retained_coefficient_critical_mass"]) - clean)
            my = min(my, float(out["uniform_square_service_lower"]))
            if bool(out["nn_seed_required"]):
                raise AssertionError("generic high-strain service route unexpectedly required NN")
    if not math.isfinite(mr):
        mr = 0.0
    if not math.isfinite(my):
        my = uniform_bounded_square_service_lower(1.0, 1.0)
    return HighStrainCriticalReentryStress(samples, ms, mp, mr, my, wr, order_fail, unit_fail, max_joint, branches)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-strain-critical-carrier-reentry"))
    ap.add_argument("--scaled-lifetime", type=float, default=1.0)
    ap.add_argument("--viscosity", type=float, default=1.0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate(args.scaled_lifetime, args.viscosity)
    out = stress(args.samples)
    (args.outdir / "high_strain_critical_carrier_reentry.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# High strain -> critical dissipation law -> own-scale service re-entry

Status: **{cert['status']}**.

The high-strain ancestor theorem already supplies the correct positive selector.  On the critical set

`G={{(j,t): M_j||P_j u(t)||_2^2>=mu_*}}`,

the actual normalized resolved-dissipation law satisfies

`D_V(G)>=D_V/2`.

Normalize `D_V|_G` and push only its deterministic shell-time mark `(j,t)`.  Every atom has `M<=N/4`.  Set `A=3M/4`.  Then

`A||P_j u(t)||_2^2 >= 8 pi^2/(25c^2)`,

and the shell's own normalized state registers this coefficient exactly into the smooth scalar envelope `Q_A=1` on the shell.  No child-scale heat edge, coherent-cell maximizer, or material label is needed to create the renewed carrier seed.

Inspect backward through one `A`-natural window with exactly three native first-stop monitors: renewed strain, role-interface coefficient obstruction, and HH-regeneration coefficient obstruction.  Coefficient hits only locate reentry of the same smooth carrier into the physical-energy gate; they do not supply work weights.  There is intentionally no material-boundary monitor because material ownership has not yet been assigned.  Exact ties remain unsplit.

If strain fires, its already named critical-dissipation cause owns the stop.  If a coefficient monitor fires, it only locates physical-energy reentry; the raw impulse owns no work.  If the interval reaches `t=0`, the initial boundary absorbs it.  Otherwise the exact Duhamel triangle gives

`|z(s)|>=|z(t)|/4`.

The surviving smooth annular carrier therefore has a scale-critical coefficient throughout the full natural slab.  The companion annular-service theorem, including inverse-heat growth of the registered analysis dual and the Arb-certified radius-3 Gaussian truncation, then gives an actual bounded own-scale displacement service

`|r|<=3/A`,  `A||delta_r Q_Au||_2^2 >= Y0`,

with

`Y0={uniform_bounded_square_service_lower(args.scaled_lifetime,args.viscosity):.12g}`

at the default requested `(c,nu)`, and normalized integrated bounded heat service at least

`{integrated_bounded_heat_service_lower(args.scaled_lifetime,args.viscosity):.12g}`.

Only **after this renewed positive service exists** do we apply Moyal and read OO/ON/NN from its actual two intrinsic endpoints.  The short high-strain renewal route therefore does not need an NN witness to be transported into the carrier.  Child-scale heat ownership, old-incident erosion, and NN-intersect-critical extraction remain valid and useful material-capacity refinements, but they are no longer prerequisites for the entrance from high strain to a renewed service carrier.

Stress: `{out.samples}` critical-dissipation-law/first-stop/service states
- minimum renewed seed critical-mass margin: `{out.minimum_seed_mass_margin:.3e}`
- minimum sampled normalized dissipation weight: `{out.minimum_normalized_weight:.3e}`
- minimum full-survivor critical-mass margin: `{out.minimum_survivor_mass_margin:.3e}`
- minimum sampled own-scale service threshold: `{out.minimum_service_threshold:.3e}`
- worst exact Duhamel residual: `{out.worst_duhamel_residual:.3e}`
- monitor-order failures: `{out.order_invariance_failures}`
- monitor-unit failures: `{out.unit_invariance_failures}`
- maximum sampled exact joint first-stop count: `{out.maximum_joint_first_stop_count}`
- outcomes: `{out.branch_counts}`

This closes the **high-strain corridor alternative** to named recursive stop / initial boundary / a completed full-natural corridor carrying own-scale coherent service.  That service is attached to the corridor and does not add a second recursion edge; `D_V` is not a reset and no whole-carrier NN assertion is made.  Universal source/SGS and genuine material-relink slab renewal remain the master-facing continuum frontier.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
