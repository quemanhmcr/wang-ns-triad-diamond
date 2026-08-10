from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Sequence

import numpy as np

from src.common_slice_coefficient_registration import (
    HH_COEFFICIENT_OBSTRUCTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
    registration_first_stop,
)


class ThresholdTopology(str, Enum):
    """Topology of a physical superlevel stop.

    CLOSED is used when the downstream certificate is valid at equality.
    STRICT is used for a state corridor whose safe set contains its boundary;
    the returned time is then the debut (infimum) of the open superlevel set.
    """

    CLOSED = "closed_superlevel"
    STRICT = "strict_superlevel_debut"


@dataclass(frozen=True)
class PhysicalPathMonitor:
    """One scalar physical observable in its native units.

    ``values`` are only a piecewise-linear regression representation.  The
    theorem stated by :func:`theorem_certificate` is a path-space statement for
    continuous/absolutely-continuous continuum observables and does not depend on
    this discretization.
    """

    cause: str
    threshold: float
    values: tuple[float, ...]
    topology: ThresholdTopology = ThresholdTopology.CLOSED

    @property
    def label(self) -> str:
        """Canonical first-stop label; ``cause`` is retained for API compatibility."""
        return self.cause


@dataclass(frozen=True)
class JointFirstExit:
    first_time: float | None
    joint_causes: tuple[str, ...]
    individual_debuts: dict[str, float | None]

    @property
    def joint_first_stops(self) -> tuple[str, ...]:
        """Canonical first-stop set; ``joint_causes`` is the compatibility field."""
        return self.joint_causes


def _validated_path(times: Sequence[float], values: Sequence[float], threshold: float) -> tuple[np.ndarray, np.ndarray]:
    t = np.asarray(times, float)
    v = np.asarray(values, float)
    if t.ndim != 1 or v.ndim != 1 or len(t) != len(v) or len(t) < 2:
        raise ValueError("matching one-dimensional path samples of length >=2 required")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(v)) or not math.isfinite(threshold):
        raise ValueError("finite path and threshold required")
    if np.any(np.diff(t) <= 0):
        raise ValueError("strictly increasing physical times required")
    return t, v


def superlevel_debut_piecewise_linear(
    times: Sequence[float],
    values: Sequence[float],
    threshold: float,
    topology: ThresholdTopology = ThresholdTopology.CLOSED,
) -> float | None:
    """Exact debut for the piecewise-linear interpolant of one physical path.

    No normalization against other observables occurs.  For a strict stop the
    returned time is ``inf{t:f(t)>threshold}``; at a genuine crossing continuity
    normally gives ``f(tau)=threshold``.  This boundary timestamp is exactly what
    is needed to order native first-stop contacts, while any
    dt-absolutely-continuous work law gives the singleton zero mass.
    """
    t, v = _validated_path(times, values, threshold)
    strict = topology is ThresholdTopology.STRICT
    if (v[0] > threshold) if strict else (v[0] >= threshold):
        return float(t[0])

    for i in range(1, len(t)):
        left = float(v[i - 1])
        right = float(v[i])
        fires = right > threshold if strict else right >= threshold
        if not fires:
            continue
        # The first endpoint in the target superlevel determines the first
        # crossing on this segment.  A strict crossing from equality has debut
        # at the left endpoint although the endpoint itself is still safe.
        if left >= threshold:
            return float(t[i - 1])
        frac = (threshold - left) / (right - left)
        return float(t[i - 1] + frac * (t[i] - t[i - 1]))
    return None


def first_physical_corridor_exit(
    times: Sequence[float],
    monitors: Sequence[PhysicalPathMonitor],
    *,
    tie_tolerance: float | None = None,
) -> JointFirstExit:
    """First boundary contact of a finite family of native-unit observables.

    Exact continuum ties are set-valued.  ``tie_tolerance`` exists only because
    floating regression paths may represent the same algebraic crossing with a
    few ulps of disagreement; it has no role in the analytic theorem.
    """
    t = np.asarray(times, float)
    if t.ndim != 1 or len(t) < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("strictly increasing times required")
    if not monitors:
        return JointFirstExit(None, (), {})
    if len({m.label for m in monitors}) != len(monitors):
        raise ValueError("one monitor per already-quotiented first-stop label required")

    debuts: dict[str, float | None] = {}
    for m in monitors:
        debuts[m.label] = superlevel_debut_piecewise_linear(t, m.values, m.threshold, m.topology)
    finite = [x for x in debuts.values() if x is not None]
    if not finite:
        return JointFirstExit(None, (), debuts)
    first = min(finite)
    if tie_tolerance is None:
        span = max(1.0, abs(float(t[-1] - t[0])), abs(first))
        tie_tolerance = 128.0 * math.ulp(span)
    if tie_tolerance < 0 or not math.isfinite(tie_tolerance):
        raise ValueError("finite nonnegative tie tolerance required")
    labels = tuple(sorted(c for c, x in debuts.items() if x is not None and abs(x - first) <= tie_tolerance))
    if not labels:
        raise AssertionError("finite first exit lost its first-stop set")
    return JointFirstExit(float(first), labels, debuts)


def rescale_monitor_units(monitor: PhysicalPathMonitor, factor: float) -> PhysicalPathMonitor:
    """Change units of one observable and its threshold together.

    Positive unit changes leave the physical stopping set exactly unchanged.
    This is the regression counterpart of the path-space invariance under any
    strictly increasing continuous reparameterization of a single observable.
    """
    if factor <= 0 or not math.isfinite(factor):
        raise ValueError("positive finite unit factor required")
    return PhysicalPathMonitor(
        cause=monitor.label,
        threshold=factor * monitor.threshold,
        values=tuple(factor * x for x in monitor.values),
        topology=monitor.topology,
    )


def moyal_cell_energy_rate_identity(
    phase_state: np.ndarray,
    phase_state_rate: np.ndarray,
    cell_mask: np.ndarray,
) -> dict[str, float]:
    """Finite Moyal regression of d||1_C F||_2^2/dt.

    The continuum identity is the same Hilbert-space calculation with
    ``F=V_{g(t)}u(t)`` and a fixed material cell ``C`` in intrinsic zeta.
    """
    F = np.asarray(phase_state, complex)
    dF = np.asarray(phase_state_rate, complex)
    mask = np.asarray(cell_mask, bool)
    if F.ndim != 1 or dF.shape != F.shape or mask.shape != F.shape:
        raise ValueError("matching one-dimensional phase state/rate/mask required")
    Fc = F[mask]
    dFc = dF[mask]
    energy = float(np.vdot(Fc, Fc).real)
    derivative = 2.0 * float(np.real(np.vdot(Fc, dFc)))
    local_rate_norm = float(np.linalg.norm(dFc))
    bound = 2.0 * math.sqrt(max(0.0, energy)) * local_rate_norm
    return {
        "cell_energy": energy,
        "cell_energy_derivative": derivative,
        "local_derivative_bound": bound,
        "margin": bound - abs(derivative),
    }


def transported_moyal_cell_rate_upper(
    *,
    cell_energy: float,
    u_rate_l2: float,
    u_l2: float,
    window_rate_l2: float,
    window_l2: float = 1.0,
) -> float:
    """Physical AC bound for an anchored material Moyal cell.

    Polarized Moyal and ``d V_g u = V_g u_t + V_{g_t}u`` give

      |E_C'| <= 2 sqrt(E_C) (||g||_2 ||u_t||_2 + ||g_t||_2 ||u||_2).

    No best-cell selector or cell-name derivative enters this estimate.
    """
    vals = (cell_energy, u_rate_l2, u_l2, window_rate_l2, window_l2)
    if any((not math.isfinite(x) or x < 0) for x in vals):
        raise ValueError("finite nonnegative Moyal rate data required")
    return 2.0 * math.sqrt(cell_energy) * (window_l2 * u_rate_l2 + window_rate_l2 * u_l2)


def circle_geodesic_holonomy(unit_holonomy: complex) -> float:
    """Branch-free distance on S^1 from the identity, in [0,pi].

    Unlike a chosen principal ``arg``, this state observable is continuous
    through the negative-real-axis branch cut.
    """
    z = complex(unit_holonomy)
    r = abs(z)
    if not math.isfinite(r) or r == 0:
        raise ValueError("nonzero finite holonomy required")
    x = max(-1.0, min(1.0, z.real / r))
    return math.acos(x)


def registration_no_hit_exhaustion(
    z_event: complex,
    z_slice: complex,
    i_hh: complex,
    i_r: complex,
    *,
    material_relink: bool = False,
) -> dict[str, object]:
    """Expose the existing exact common-slice theorem in survivor/stop language."""
    out = registration_first_stop(z_event, z_slice, i_hh, i_r, material_relink=material_relink)
    if bool(out["continuing"]):
        return {
            "classification": "registered_generated_survivor",
            "first_stops": (),
            "stop_causes": (),
            "registered_amplitude_lower": float(out["registered_amplitude_lower"]),
            "event_amplitude": float(out["event_amplitude"]),
            "requires_physical_energy_reentry": False,
            "coefficient_impulses_used_as_work": False,
        }
    stops = tuple(str(x) for x in out.get("first_stops", (str(out["branch"]),)))
    needs_energy_reentry = any(
        label in {ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION, HH_COEFFICIENT_OBSTRUCTION}
        for label in stops
    )
    return {
        "classification": (
            "coefficient_obstruction_energy_reentry"
            if needs_energy_reentry
            else "named_backward_physical_stop"
        ),
        "first_stops": stops,
        # Backward-compatible alias; coefficient labels are obstructions, not causes.
        "stop_causes": stops,
        "registered_amplitude_lower": 0.0,
        "event_amplitude": float(out["event_amplitude"]),
        "requires_physical_energy_reentry": needs_energy_reentry,
        "coefficient_impulses_used_as_work": False,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_SMOOTH_SGS_MEASURABLE_FIRST_EXIT_AND_LOCAL_NO_HIT_EXHAUSTION__RECURSIVE_REENTRY_REMAINS",
        "causal_filtration": (
            "event/support selection (Xi and fixed transfer cost); smooth slab observables; physical HH-work measure; "
            "backward common-slice registration; only then ancestry/reuse information"
        ),
        "closed_hit_measurability": (
            "for f in C([a,b]), tau_ge=inf{t:f(t)>=theta} is Borel because "
            "{tau_ge<=q}={sup_[a,q] f>=theta} and the sup functional is continuous in the uniform norm"
        ),
        "strict_exit_measurability": (
            "tau_>=inf{t:f(t)>theta} is Borel by rational-time sup sets; if f starts safe and tau_> is finite, "
            "continuity gives f(tau_>)=theta and the open-superlevel crossing germ is the physical state exit"
        ),
        "joint_ties": (
            "for finitely many native observables, tau=min_r tau_r and J={r:tau_r=tau} are measurable; "
            "J is retained as an unsplit first-stop set with no theorem-name priority or common-unit weights; "
            "coefficient members require physical-energy routing before they become causes"
        ),
        "unit_invariance": (
            "each stopping set is unchanged under its own strictly increasing continuous change of physical units; "
            "heterogeneous observables are never added or normalized against one another"
        ),
        "smooth_sgs_regularities": (
            "on every pre-singular smooth interval, strict-lowpass resolved fields are smooth; coherent X,L ODEs are C1; "
            "strain/deformation/source/work integrals are AC; singular-value/aspect/radius states are continuous"
        ),
        "material_cell": (
            "anchor a measurable coherent cell in intrinsic zeta=(L^-1 X/2,L^T k), which common affine transport fixes; "
            "E_C(t)=||1_C V_{g(t)}u(t)||_2^2 is AC with |E_C'|<=2 sqrt(E_C)(||g|| ||u_t||+||g_t|| ||u||)"
        ),
        "cell_boundaries": (
            "Moyal energy/work have phase-space densities, so dyadic cell boundaries are null; representation-cell boundary motion carries no physical relink mass"
        ),
        "phase": (
            "use the branch-free S1 geodesic holonomy acos(Re h) rather than a chosen arg branch, eliminating artificial phase jumps"
        ),
        "generated_no_hit": (
            "backward adjoint impulses are AC in the backward endpoint; if no role-interface coefficient, HH coefficient, material-relink, or t=0 boundary obstruction hits, "
            "the exact triangle identity registers the generated coefficient with the existing 1/4 lower factor; a coefficient hit only locates Q^2 physical-energy reentry"
        ),
        "flat_no_hit": (
            "on the retained low-transfer block, absence of the already-certified strain/source/deformation/aspect/radius/phase service exits leaves exactly the existing coherent Kelvin-flat alternative"
        ),
        "scope": (
            "this proves local measurable first-hit extraction for any recursively selected smooth-SGS block once that block is supplied; "
            "it does not yet prove that every routed source/dissipation/relink/HH-generation re-entry again satisfies the selector hypotheses"
        ),
    }


@dataclass(frozen=True)
class FirstHitStress:
    samples: int
    order_invariance_failures: int
    unit_rescaling_failures: int
    maximum_joint_first_stop_count: int
    minimum_moyal_rate_margin: float
    worst_phase_branch_cut_gap: float
    minimum_registered_fraction_margin: float


def stress(samples: int = 50_000, seed: int = 20260809) -> FirstHitStress:
    rng = np.random.default_rng(seed)
    order_fail = unit_fail = 0
    max_joint = 0
    min_moyal = float("inf")
    worst_phase = 0.0
    min_reg = float("inf")
    times = np.array([0.0, 0.2, 0.5, 0.7, 1.0])

    for _ in range(samples):
        # Two physically unrelated observables are engineered to touch their own
        # native thresholds at the same physical time.  No relative magnitude is
        # ever used to decide ownership.
        tau = 0.5
        th1 = float(math.exp(rng.uniform(-5.0, 4.0)))
        th2 = float(math.exp(rng.uniform(-5.0, 4.0)))
        s1 = float(math.exp(rng.uniform(-3.0, 3.0)))
        s2 = float(math.exp(rng.uniform(-3.0, 3.0)))
        v1 = tuple(th1 + s1 * (times - tau))
        v2 = tuple(th2 + s2 * (times - tau))
        late_th = float(math.exp(rng.uniform(-5.0, 4.0)))
        late_s = float(math.exp(rng.uniform(-3.0, 3.0)))
        v3 = tuple(late_th + late_s * (times - 0.8))
        mons = (
            PhysicalPathMonitor("strain_action", th1, v1, ThresholdTopology.CLOSED),
            PhysicalPathMonitor("material_state_exit", th2, v2, ThresholdTopology.STRICT),
            PhysicalPathMonitor("source_action", late_th, v3, ThresholdTopology.CLOSED),
        )
        out = first_physical_corridor_exit(times, mons, tie_tolerance=2e-12)
        if out.first_time is None or abs(out.first_time - tau) > 2e-12:
            raise AssertionError("piecewise-linear first physical exit moved away from the exact crossing")
        if set(out.joint_first_stops) != {"strain_action", "material_state_exit"}:
            raise AssertionError("exact joint first-stop set was not retained")
        max_joint = max(max_joint, len(out.joint_first_stops))

        perm = rng.permutation(len(mons))
        reordered = tuple(mons[int(i)] for i in perm)
        out2 = first_physical_corridor_exit(times, reordered, tie_tolerance=2e-12)
        if out2.first_time != out.first_time or out2.joint_first_stops != out.joint_first_stops:
            order_fail += 1
            raise AssertionError("first exit depended on monitor enumeration")

        factors = [float(math.exp(rng.uniform(-12.0, 12.0))) for _ in mons]
        rescaled = tuple(rescale_monitor_units(m, a) for m, a in zip(mons, factors))
        out3 = first_physical_corridor_exit(times, rescaled, tie_tolerance=2e-10)
        if out3.first_time is None or abs(out3.first_time - tau) > 2e-10 or set(out3.joint_first_stops) != set(out.joint_first_stops):
            unit_fail += 1
            raise AssertionError("first exit depended on independent physical unit choices")

        n = int(rng.integers(2, 24))
        F = rng.normal(size=n) + 1j * rng.normal(size=n)
        dF = rng.normal(size=n) + 1j * rng.normal(size=n)
        mask = rng.random(n) < 0.5
        if not bool(mask.any()):
            mask[int(rng.integers(0, n))] = True
        mo = moyal_cell_energy_rate_identity(F, dF, mask)
        min_moyal = min(min_moyal, float(mo["margin"]))
        if float(mo["margin"]) < -3e-12 * max(1.0, float(mo["local_derivative_bound"])):
            raise AssertionError("material Moyal cell energy rate violated Cauchy-Schwarz")

        eps = float(rng.uniform(1e-9, 0.2))
        left = circle_geodesic_holonomy(complex(math.cos(math.pi - eps), math.sin(math.pi - eps)))
        right = circle_geodesic_holonomy(complex(math.cos(-math.pi + eps), math.sin(-math.pi + eps)))
        gap = abs(left - right)
        worst_phase = max(worst_phase, gap)
        if gap > 3e-12:
            raise AssertionError("branch-free phase distance jumped at the principal-arg cut")

        A = float(math.exp(rng.uniform(-4.0, 4.0)))
        theta = float(rng.uniform(-math.pi, math.pi))
        ze = A * complex(math.cos(theta), math.sin(theta))
        rp = float(rng.uniform(-math.pi, math.pi))
        hp = float(rng.uniform(-math.pi, math.pi))
        ir = float(rng.uniform(0.0, 0.24)) * A * complex(math.cos(rp), math.sin(rp))
        ih = float(rng.uniform(0.0, 0.49)) * A * complex(math.cos(hp), math.sin(hp))
        zs = ze - ir - ih
        reg = registration_no_hit_exhaustion(ze, zs, ih, ir)
        if reg["classification"] != "registered_generated_survivor":
            raise AssertionError("no-hit exact adjoint event failed to register as generated survivor")
        margin = float(reg["registered_amplitude_lower"]) - 0.25 * float(reg["event_amplitude"])
        min_reg = min(min_reg, margin)
        if margin < -3e-12 * max(1.0, A):
            raise AssertionError("registered survivor lost the exact quarter coefficient")

    return FirstHitStress(
        samples=samples,
        order_invariance_failures=order_fail,
        unit_rescaling_failures=unit_fail,
        maximum_joint_first_stop_count=max_joint,
        minimum_moyal_rate_margin=min_moyal,
        worst_phase_branch_cut_gap=worst_phase,
        minimum_registered_fraction_margin=min_reg,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-smooth-sgs-first-hit-extraction"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "smooth_sgs_first_hit_extraction.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Smooth-SGS first-hit extraction from physical observables\n\nStatus: **{cert['status']}**.\n\nThe extraction is a **causal filtration**, not one artificial vector clock.  Event/support facts are resolved first; smooth resolved/coherent dynamics then generate continuous or absolutely-continuous physical observables; actual positive HH work creates the child transfer measure; common-slice registration is a backward adjoint problem; Shannon/Renyi information is only formed after the parent law exists.\n\nFor a finite continuous family `f_r(t)` in their own physical units, define each certified superlevel hit/debut separately and set `tau=min_r tau_r`.  Closed-threshold hitting times are Borel because `{{tau_r<=q}}={{sup_[a,q] f_r>=theta_r}}`; strict state exits are Borel by the analogous rational-time open-superlevel formula.  The simultaneous first-stop set `J={{r:tau_r=tau}}` is therefore measurable and finite.  No theorem-name order and no normalization of strain, source work, Moyal relink energy, or phase is needed.  Independent monotone changes of units leave every stopping set unchanged.\n\nMaterial coherence is anchored rather than re-selected.  In intrinsic `zeta=(L^-1 X/2,L^T k)`, common affine transport fixes a cell exactly.  With `F(t)=V_{{g(t)}}u(t)`,\n\n`E_C(t)=||1_C F(t)||_2^2`,\n\n`|E_C'(t)| <= 2 sqrt(E_C) [ ||g||_2 ||u_t||_2 + ||g_t||_2 ||u||_2 ]`.\n\nThus the physical Moyal content of the transported cell is AC.  Dyadic cell boundaries have zero Moyal energy/work because these measures have phase-space densities, so a representation boundary cannot create relink mass.  A genuine new material address/relink is an event of the physical cell measure, not an `argmax` label chatter.\n\nHelical phase is treated similarly: use the branch-free circle distance `acos(Re h)` for unit holonomy `h`, rather than a chosen principal-angle branch.\n\nFor generated events, the backward impulses in `z(t)=z(s)+I_HH+I_R` are AC in the endpoint.  If no role-interface coefficient, HH coefficient, material-relink, or initial-boundary obstruction hits, the existing exact triangle gate gives `|z(s)|>=|z(t)|/4`; hence the no-hit event is precisely a registered generated survivor.  At block level, absence of the already-certified service exits leaves the existing coherent Kelvin-flat alternative.\n\nStress: `{out.samples}` physical-path/Moyal/phase/registration regressions\n- monitor-order failures: `{out.order_invariance_failures}`\n- independent-unit-rescaling failures: `{out.unit_rescaling_failures}`\n- largest exact joint first-stop set sampled: `{out.maximum_joint_first_stop_count}`\n- minimum Moyal rate inequality margin: `{out.minimum_moyal_rate_margin:.3e}`\n- worst phase branch-cut gap: `{out.worst_phase_branch_cut_gap:.3e}`\n- minimum registered quarter-factor margin: `{out.minimum_registered_fraction_margin:.3e}`\n\nThis closes the **local measurable first-hit extraction once a smooth-SGS block has been recursively selected**.  It does not yet prove universal recursive re-entry: after a routed source, critical dissipation, material/new-ancestry, or HH-generation stop, one must still show that the next continuum selector again supplies a block satisfying the hypotheses.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
