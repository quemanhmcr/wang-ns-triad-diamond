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
from src.nn_critical_heat_carrier_seed import (
    LOW_STRAIN_ACTION,
    persistent_seed_low_low_gap,
    renewal_carrier_critical_mass_lower,
)
from src.smooth_sgs_first_hit_extraction import (
    PhysicalPathMonitor,
    ThresholdTopology,
    first_physical_corridor_exit,
)


def renewed_natural_duration(renewal_frequency: float, scaled_lifetime: float) -> float:
    A = float(renewal_frequency)
    c = float(scaled_lifetime)
    if A <= 0 or c <= 0 or not all(math.isfinite(x) for x in (A, c)):
        raise ValueError("positive finite renewed frequency and scaled lifetime required")
    # Form sqrt(c)/A first.  The mathematically equivalent c/(A*A) can
    # overflow its denominator and return zero even when the natural duration
    # itself is a positive representable subnormal number.
    root_duration = math.sqrt(c) / A
    duration = root_duration * root_duration
    if duration <= 0 or not math.isfinite(duration):
        raise ValueError("renewed natural duration is outside the finite floating certificate range")
    return duration


def backward_natural_endpoint(event_time: float, renewal_frequency: float, scaled_lifetime: float) -> dict[str, float | bool]:
    t = float(event_time)
    if t < 0 or not math.isfinite(t):
        raise ValueError("finite nonnegative event time required")
    T = renewed_natural_duration(renewal_frequency, scaled_lifetime)
    elapsed = min(t, T)
    hits_boundary = t <= T
    s = 0.0 if hits_boundary else t - T
    return {
        "event_time": t,
        "natural_duration": T,
        "backward_endpoint": s,
        "endpoint_elapsed_from_event": elapsed,
        "elapsed_available": elapsed,
        "hits_initial_boundary": hits_boundary,
        "full_natural_interval_before_boundary": t >= T,
    }


def inherited_seed_critical_mass_lower(scaled_lifetime: float) -> float:
    """After the exact 1/4 coefficient gate, A|z(s)|^2 retains 1/16 of seed mass."""
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    return (INHERIT_FRACTION**2) * renewal_carrier_critical_mass_lower(c)


def _validated_path(name: str, values: Sequence[float], n: int) -> tuple[float, ...]:
    x = tuple(float(v) for v in values)
    if len(x) != n or any(not math.isfinite(v) for v in x):
        raise ValueError(f"{name} must be a matching finite path")
    return x


def seed_backward_first_hit(
    elapsed_times: Sequence[float],
    *,
    terminal_amplitude: float,
    strain_action: Sequence[float],
    residual_impulse_abs: Sequence[float],
    hh_impulse_abs: Sequence[float],
    material_boundary_distance: Sequence[float],
    tie_tolerance: float | None = None,
) -> dict[str, object]:
    """First physical obstruction while elapsed backward time grows from zero.

    ``elapsed_times`` is ell=t_event-s.  All monitors remain in native units:
    strain action, coefficient impulses, and geometric material distance are not
    normalized against one another.  Exact simultaneous hits remain unsplit.

    The material observable is ``-dist(endpoint witnesses, boundary(old pool))``.
    It starts negative for an interior NN witness and reaches zero exactly at
    first material-boundary contact.  For the two endpoints use their minimum
    boundary distance before calling this helper.
    """
    ell = np.asarray(elapsed_times, float)
    if ell.ndim != 1 or len(ell) < 2 or abs(float(ell[0])) > 1e-14 or np.any(np.diff(ell) <= 0):
        raise ValueError("elapsed backward times must start at zero and increase strictly")
    if np.any(~np.isfinite(ell)):
        raise ValueError("finite elapsed times required")
    amp = float(terminal_amplitude)
    if amp <= 0 or not math.isfinite(amp):
        raise ValueError("positive finite terminal amplitude required")
    n = len(ell)
    K = _validated_path("strain action", strain_action, n)
    IR = _validated_path("residual impulse", residual_impulse_abs, n)
    IH = _validated_path("HH impulse", hh_impulse_abs, n)
    dmat = _validated_path("material boundary distance", material_boundary_distance, n)
    if min(K) < 0 or min(IR) < 0 or min(IH) < 0 or min(dmat) < 0:
        raise ValueError("nonnegative physical monitor paths required")

    monitors = (
        PhysicalPathMonitor("high_strain_critical_dissipation", float(LOW_STRAIN_ACTION), K, ThresholdTopology.CLOSED),
        PhysicalPathMonitor(ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, RESIDUAL_FRACTION * amp, IR, ThresholdTopology.CLOSED),
        PhysicalPathMonitor(HH_COEFFICIENT_OBSTRUCTION, GENERATED_FRACTION * amp, IH, ThresholdTopology.CLOSED),
        PhysicalPathMonitor("material_boundary_contact", 0.0, tuple(-x for x in dmat), ThresholdTopology.CLOSED),
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
        "requires_physical_energy_reentry": needs_energy_reentry,
        "coefficient_impulses_used_as_work": False,
    }


def natural_corridor_outcome(
    *,
    event_time: float,
    renewal_frequency: float,
    scaled_lifetime: float,
    terminal_coefficient: complex,
    endpoint_coefficient: complex,
    hh_impulse: complex,
    residual_interface_impulse: complex,
    first_hit: dict[str, object],
) -> dict[str, object]:
    """Full backward natural-window outcome after the measurable first-hit pass.

    If a monitor hits before or at the requested endpoint, return its unsplit
    joint first-stop set.  If t=0 truncates the natural window, the initial
    boundary is absorbing.  Otherwise the exact Duhamel identity plus strict
    subthreshold impulses forces |z(s)|>=|z(t)|/4.

    A large residual coefficient impulse is only a role-interface obstruction
    locator: this function does not convert it into physical work.  A large HH
    coefficient impulse is likewise an obstruction and must re-enter the existing
    physical-energy gate before HH generation or any causal weighting is named.
    """
    geom = backward_natural_endpoint(event_time, renewal_frequency, scaled_lifetime)
    zt = complex(terminal_coefficient)
    zs = complex(endpoint_coefficient)
    ih = complex(hh_impulse)
    ir = complex(residual_interface_impulse)
    amp = abs(zt)
    if amp <= 0:
        raise ValueError("nonzero terminal coefficient required")
    res = abs(exact_adjoint_residual(zt, zs, ih, ir))
    tol = 3e-12 * max(1.0, amp, abs(zs), abs(ih), abs(ir))
    if res > tol:
        raise ValueError("seed natural-window Duhamel decomposition is not exact")

    hit_time = first_hit.get("first_elapsed")
    causes = tuple(str(x) for x in first_hit.get("joint_first_stops", first_hit.get("joint_causes", ())))
    needs_energy_reentry = any(
        label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for label in causes
    )
    elapsed_available = float(geom["elapsed_available"])
    if hit_time is not None and float(hit_time) <= elapsed_available + 2e-12 * max(1.0, elapsed_available):
        return {
            "classification": "named_first_stop",
            "joint_causes": causes,
            "joint_first_stops": causes,
            "first_elapsed": float(hit_time),
            "primary_selected": False,
            "terminal_amplitude": amp,
            "duhamel_residual": res,
            "requires_physical_energy_reentry": needs_energy_reentry,
            "coefficient_impulses_used_as_work": False,
        }

    if bool(geom["hits_initial_boundary"]):
        return {
            "classification": "initial_boundary_root",
            "joint_causes": ("t=0",),
            "joint_first_stops": ("t=0",),
            "first_elapsed": elapsed_available,
            "primary_selected": False,
            "terminal_amplitude": amp,
            "duhamel_residual": res,
            "requires_physical_energy_reentry": False,
            "coefficient_impulses_used_as_work": False,
        }

    if abs(ir) >= RESIDUAL_FRACTION * amp - 3e-12 * max(1.0, amp):
        raise ValueError("endpoint residual impulse contradicts no-hit corridor")
    if abs(ih) >= GENERATED_FRACTION * amp - 3e-12 * max(1.0, amp):
        raise ValueError("endpoint HH impulse contradicts no-hit corridor")
    inherited = abs(zs)
    clean = INHERIT_FRACTION * amp
    if inherited < clean - 4e-12 * max(1.0, amp):
        raise AssertionError("full natural corridor lost the exact 1/4 carrier coefficient")
    A = float(renewal_frequency)
    retained_mass = A * inherited * inherited
    terminal_mass = A * amp * amp
    return {
        "classification": "full_natural_corridor_survivor",
        "joint_causes": (),
        "joint_first_stops": (),
        "natural_duration": float(geom["natural_duration"]),
        "backward_endpoint": float(geom["backward_endpoint"]),
        "terminal_amplitude": amp,
        "endpoint_amplitude": inherited,
        "registered_amplitude_lower": clean,
        "terminal_critical_mass": terminal_mass,
        "endpoint_critical_mass": retained_mass,
        "clean_endpoint_critical_mass_lower": (INHERIT_FRACTION**2) * terminal_mass,
        "duhamel_residual": res,
        "nn_endpoint_witness_survives": True,
        "whole_carrier_declared_nn": False,
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def theorem_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    seed = renewal_carrier_critical_mass_lower(c)
    inherited = inherited_seed_critical_mass_lower(c)
    expected = math.pi * math.pi / (50.0 * c * c)
    if abs(inherited - expected) > 2e-13 * max(1.0, expected):
        raise AssertionError("clean inherited seed critical mass identity changed")
    gap = persistent_seed_low_low_gap(LOW_STRAIN_ACTION)
    if gap <= 0:
        raise AssertionError("low-low moat unavailable inside seed strain corridor")
    return {
        "status": "EXACT_NN_SEED_BACKWARD_NATURAL_FIRST_STOP_CORRIDOR__ONE_QUARTER_CRITICAL_SURVIVOR_OR_NAMED_STOP__MATERIAL_ATTACHMENT_REMAINS",
        "interval": "for a seed at time t and scale A, inspect the backward natural interval [max(0,t-cA^-2),t] using elapsed time ell=t-s",
        "first_hit": "native-unit AC monitors are K_A, |I_role-interface|, |I_HH| and minus the minimum distance of the two NN endpoint witnesses to the old-material boundary; exact ties remain an unsplit first-stop set",
        "outer_role": "while K_A<1/30, the dual-transported smooth Q_A carrier stays strictly above the A/2 low-low output; the exact outer equation has only HH source plus the nonaffine role-interface residual",
        "duhamel": "z(t)=z(s)+I_HH[s,t]+I_interface[s,t]; if neither impulse hits A_z/2 or A_z/4 and no material/strain stop occurs, |z(s)|>=|z(t)|/4",
        "critical_survival": f"the seed starts with A|z(t)|^2>={seed:.12g}; a full no-hit natural interval retains A|z(s)|^2>={inherited:.12g}=pi^2/(50c^2) at c={c:.12g}",
        "material": "the two NN heat-edge endpoint witnesses remain outside the old pool until their continuous boundary distance first reaches zero; this preserves witness provenance only and never declares the whole Q_Au carrier NN material",
        "boundary": "if the requested backward natural interval reaches t=0 first, the initial boundary is absorbing rather than a free survivor",
        "interface_scope": "a large role-interface coefficient impulse is not promoted to physical work or a new currency; it only locates smooth-carrier physical-energy reentry, whose actual gauge-quotiented native interface work may then route to conservative relink or existing strain",
        "hh_scope": "a large HH coefficient impulse is an obstruction locator and must pass through the existing physical-energy causal gate before HH generation or any work-based probability is named",
        "scope": "this closes temporal coefficient survival-or-first-stop for the high-strain carrier seed itself; attaching the surviving coefficient energy to its NN material witness and proving full efficiency/service renewal, plus universal source/relink re-entry, remain open",
    }


@dataclass(frozen=True)
class TemporalSeedStress:
    samples: int
    minimum_clean_critical_mass_margin: float
    minimum_survivor_coefficient_margin: float
    minimum_low_low_gap: float
    worst_duhamel_residual: float
    order_invariance_failures: int
    unit_rescaling_failures: int
    maximum_joint_first_stop_count: int
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260809) -> TemporalSeedStress:
    rng = np.random.default_rng(seed)
    mm = mc = mg = float("inf")
    wr = 0.0
    order_fail = unit_fail = 0
    max_joint = 0
    branches: dict[str, int] = {}

    for _ in range(samples):
        c = float(math.exp(rng.uniform(-2.0, 2.0)))
        A = float(math.exp(rng.uniform(-1.0, 6.0)))
        clean_seed = renewal_carrier_critical_mass_lower(c)
        amp = math.sqrt(float(rng.uniform(1.0, 4.0)) * clean_seed / A)
        phase = float(rng.uniform(-math.pi, math.pi))
        zt = amp * complex(math.cos(phase), math.sin(phase))
        T = renewed_natural_duration(A, c)
        event_time = float(rng.uniform(1.05, 3.0)) * T
        n = 5
        ell = np.linspace(0.0, T, n)

        mode = int(rng.integers(0, 5))
        K_end = float(rng.uniform(0.2, 0.9)) * LOW_STRAIN_ACTION
        IR_end = float(rng.uniform(0.05, 0.8)) * RESIDUAL_FRACTION * amp
        IH_end = float(rng.uniform(0.05, 0.8)) * GENERATED_FRACTION * amp
        d0 = float(math.exp(rng.uniform(-4.0, 1.0)))
        d_end = float(rng.uniform(0.2, 1.0)) * d0
        if mode == 1:
            K_end = float(rng.uniform(1.0, 1.5)) * LOW_STRAIN_ACTION
        elif mode == 2:
            IR_end = float(rng.uniform(1.0, 1.5)) * RESIDUAL_FRACTION * amp
        elif mode == 3:
            IH_end = float(rng.uniform(1.0, 1.5)) * GENERATED_FRACTION * amp
        elif mode == 4:
            d_end = 0.0

        K = tuple(np.linspace(0.0, K_end, n))
        IR = tuple(np.linspace(0.0, IR_end, n))
        IH = tuple(np.linspace(0.0, IH_end, n))
        dmat = tuple(np.linspace(d0, d_end, n))
        hit = seed_backward_first_hit(
            ell,
            terminal_amplitude=amp,
            strain_action=K,
            residual_impulse_abs=IR,
            hh_impulse_abs=IH,
            material_boundary_distance=dmat,
            tie_tolerance=2e-12,
        )
        max_joint = max(max_joint, len(tuple(hit["joint_first_stops"])))

        ir_phase = float(rng.uniform(-math.pi, math.pi))
        ih_phase = float(rng.uniform(-math.pi, math.pi))
        ir = IR_end * complex(math.cos(ir_phase), math.sin(ir_phase))
        ih = IH_end * complex(math.cos(ih_phase), math.sin(ih_phase))
        zs = zt - ir - ih
        out = natural_corridor_outcome(
            event_time=event_time,
            renewal_frequency=A,
            scaled_lifetime=c,
            terminal_coefficient=zt,
            endpoint_coefficient=zs,
            hh_impulse=ih,
            residual_interface_impulse=ir,
            first_hit=hit,
        )
        b = str(out["classification"])
        branches[b] = branches.get(b, 0) + 1
        wr = max(wr, float(out["duhamel_residual"]))
        if b == "full_natural_corridor_survivor":
            cm = float(out["endpoint_critical_mass"]) - inherited_seed_critical_mass_lower(c)
            mm = min(mm, cm)
            cc = float(out["endpoint_amplitude"]) - INHERIT_FRACTION * amp
            mc = min(mc, cc)
            if cm < -4e-12 * max(1.0, inherited_seed_critical_mass_lower(c)):
                raise AssertionError("full natural survivor lost clean critical mass")
            if cc < -4e-12 * max(1.0, amp):
                raise AssertionError("full natural survivor lost 1/4 coefficient")

        gap = persistent_seed_low_low_gap(float(min(K_end, LOW_STRAIN_ACTION)))
        mg = min(mg, gap)
        if gap <= 0:
            raise AssertionError("safe seed corridor lost low-low support moat")

        # Unit/order invariance of the native monitor first-hit set.  We call the
        # shared theorem directly with independently rescaled monitors.
        mons = (
            PhysicalPathMonitor("strain", LOW_STRAIN_ACTION, K),
            PhysicalPathMonitor("interface", RESIDUAL_FRACTION * amp, IR),
            PhysicalPathMonitor("hh", GENERATED_FRACTION * amp, IH),
            PhysicalPathMonitor("material", 0.0, tuple(-x for x in dmat)),
        )
        base = first_physical_corridor_exit(ell, mons, tie_tolerance=2e-12)
        perm = rng.permutation(len(mons))
        reordered = tuple(mons[int(i)] for i in perm)
        chk = first_physical_corridor_exit(ell, reordered, tie_tolerance=2e-12)
        if base.first_time != chk.first_time or base.joint_first_stops != chk.joint_first_stops:
            order_fail += 1
            raise AssertionError("seed first stop depended on monitor order")
        factors = [float(math.exp(rng.uniform(-8.0, 8.0))) for _ in mons]
        scaled = tuple(
            PhysicalPathMonitor(m.label, m.threshold * f, tuple(f * x for x in m.values), m.topology)
            for m, f in zip(mons, factors)
        )
        chk2 = first_physical_corridor_exit(ell, scaled, tie_tolerance=2e-10)
        if (base.first_time is None) != (chk2.first_time is None):
            unit_fail += 1
            raise AssertionError("seed first stop changed under independent units")
        if base.first_time is not None and (
            abs(float(base.first_time) - float(chk2.first_time)) > 2e-10
            or set(base.joint_first_stops) != set(chk2.joint_first_stops)
        ):
            unit_fail += 1
            raise AssertionError("seed first stop changed under independent units")

    # Explicit absorbing-boundary sample.
    c = 1.0
    A = 2.0
    T = renewed_natural_duration(A, c)
    amp = math.sqrt(renewal_carrier_critical_mass_lower(c) / A)
    ell = np.linspace(0.0, 0.5 * T, 5)
    hit = seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0, 0.1 * LOW_STRAIN_ACTION, 5),
        residual_impulse_abs=np.linspace(0, 0.1 * RESIDUAL_FRACTION * amp, 5),
        hh_impulse_abs=np.linspace(0, 0.1 * GENERATED_FRACTION * amp, 5),
        material_boundary_distance=np.ones(5),
    )
    ir = 0.1 * RESIDUAL_FRACTION * amp
    ih = 0.1j * GENERATED_FRACTION * amp
    out = natural_corridor_outcome(
        event_time=0.5 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
        terminal_coefficient=amp,
        endpoint_coefficient=amp - ir - ih,
        hh_impulse=ih,
        residual_interface_impulse=ir,
        first_hit=hit,
    )
    branches[str(out["classification"])] = branches.get(str(out["classification"]), 0) + 1
    if out["classification"] != "initial_boundary_root":
        raise AssertionError("t=0 failed to absorb truncated seed corridor")

    if not math.isfinite(mm):
        mm = 0.0
    if not math.isfinite(mc):
        mc = 0.0
    return TemporalSeedStress(samples, mm, mc, mg, wr, order_fail, unit_fail, max_joint, branches)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-nn-seed-temporal-first-stop"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "nn_seed_temporal_first_stop.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# NN-critical carrier seed: backward natural first-stop corridor\n\nStatus: **{cert['status']}**.\n\nStart from one atom of the already-certified NN-critical heat law.  Its lower-scale smooth carrier `w=Q_Au` has terminal coefficient `z(t)` with\n\n`A|z(t)|^2 >= 8 pi^2/(25 c^2)`.\n\nInspect the **backward natural interval**\n\n`[max(0,t-c A^-2), t]`.\n\nUse elapsed backward time `ell=t-s`.  No common scalar clock is introduced.  Four native continuous/absolutely-continuous observables are monitored independently:\n\n- renewed strain action `K_A[s,t]`, with first contact at `1/30`;\n- magnitude of the nonaffine role-interface coefficient obstruction, threshold `|z(t)|/4`;\n- magnitude of the HH coefficient obstruction, threshold `|z(t)|/2`;\n- minus the minimum distance of the two retained NN heat-edge endpoint witnesses to the old-material boundary, threshold `0`.\n\nTheir first backward debut is measurable by the existing smooth first-hit theorem; exact ties remain one unsplit first-stop set.  Until the strain face is hit, the dual-transported smooth role stays above the low--low output because\n\n`(3/5)e^(-1/30)A > A/2`.\n\nThe exact outer moving-role equation and adjoint interaction picture give\n\n`z(t)=z(s)+I_HH[s,t]+I_interface[s,t]`.\n\nTherefore there are only three possibilities before the requested natural endpoint:\n\n1. a first obstruction occurs: strain, interface coefficient obstruction, HH coefficient obstruction, or material-boundary contact;\n2. the interval reaches `t=0`, which is absorbing;\n3. no stop occurs through a full `cA^-2` interval, in which case the exact triangle inequality forces\n\n`|z(s)| >= |z(t)|/4`.\n\nThe surviving carrier is still quantitatively critical:\n\n`A|z(s)|^2 >= (1/16) A|z(t)|^2 >= pi^2/(50 c^2)`.\n\nThe material statement remains exact but deliberately narrow.  The two NN heat-edge endpoint witnesses stay outside the old pool until their continuous boundary distance reaches zero.  A no-hit survivor therefore retains an NN **witness**, not a claim that the whole Fourier carrier energy is new material.\n\nA large interface coefficient impulse is not promoted to physical work.  It only locates reentry of the same smooth carrier into the physical-energy gate.  If that gate selects actual gauge-quotiented native interface work, the quadratic-carrier theorem routes it to conservative relink or existing strain.  A large HH coefficient impulse likewise only locates physical-energy reentry; actual HH generation is named there.\n\nStress: `{out.samples}` backward natural corridors\n- minimum clean survivor critical-mass margin: `{out.minimum_clean_critical_mass_margin:.3e}`\n- minimum survivor 1/4-coefficient margin: `{out.minimum_survivor_coefficient_margin:.3e}`\n- minimum safe low--low gap: `{out.minimum_low_low_gap:.6e}`\n- worst exact Duhamel residual: `{out.worst_duhamel_residual:.3e}`\n- monitor-order failures: `{out.order_invariance_failures}`\n- independent-unit failures: `{out.unit_rescaling_failures}`\n- maximum exact joint first-stop count: `{out.maximum_joint_first_stop_count}`\n- outcomes: `{out.branch_counts}`\n\nThis closes **temporal coefficient survival-or-first-stop for the high-strain carrier seed**.  The remaining high-strain seam is material/efficiency attachment and global physical-owner assembly.  Universal source/relink renewal remains open.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
