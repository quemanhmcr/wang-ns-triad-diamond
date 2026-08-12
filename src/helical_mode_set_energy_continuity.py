from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from src.continuum_helical_edge_measure_registration import HelicalModeIdentity
from src.cyclic_helical_triad_donor_kernel import (
    CyclicTriadMeasureKernel,
    cyclic_triad_measure_kernel,
    register_closed_helical_triad,
    signed_good_integer_triad,
)
from src.helical import stable_norm3

STATUS = (
    "EXACT_HELICAL_MODE_SET_ENERGY_CONTINUITY__CYCLIC_DONOR_BOUNDARY_FLUX__"
    "INTERNAL_FLOW_CANCELLATION__MODEWISE_VISCOUS_STOCK_BALANCE__"
    "NO_FIFO_LIFO_OR_GROSS_TRANSFER_BUDGET"
)


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


def _native_residual(actual: float, expected: float, scale: float) -> float:
    a = float(actual)
    b = float(expected)
    s = float(scale)
    if not all(math.isfinite(v) for v in (a, b, s)) or s < 0.0:
        raise ValueError("finite nonnegative native comparison scale required")
    gap = abs(a - b)
    if s == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / s


@dataclass(frozen=True)
class HelicalModeFlowAtom:
    donor_mode: HelicalModeIdentity
    recipient_mode: HelicalModeIdentity
    physical_work_mass: float

    def __post_init__(self) -> None:
        if self.donor_mode == self.recipient_mode:
            raise ValueError("cyclic energy flow must connect distinct physical helical modes")
        _finite_nonnegative(self.physical_work_mass, "helical mode-flow work mass")
        if self.physical_work_mass <= 0.0:
            raise ValueError("helical mode-flow atom must carry positive physical work")


def flow_atoms_from_cyclic_kernel(kernel: CyclicTriadMeasureKernel) -> tuple[HelicalModeFlowAtom, ...]:
    if not kernel.numerically_resolved_transport:
        raise ValueError("helical mode flow refuses numerically unresolved cyclic Hahn signs")
    atoms = tuple(
        HelicalModeFlowAtom(
            donor_mode=atom.donor_child_mode,
            recipient_mode=atom.recipient_child_mode,
            physical_work_mass=atom.physical_work_mass,
        )
        for atom in kernel.atoms
    )
    if not atoms:
        raise ValueError("resolved cyclic kernel carries no physical donor flow")
    return atoms


@dataclass(frozen=True)
class HelicalModeSetBoundaryBalance:
    modes: frozenset[HelicalModeIdentity]
    internal_flow: float
    inward_boundary_flow: float
    outward_boundary_flow: float
    recipient_positive_work: float
    donor_negative_work: float
    signed_nonlinear_work: float
    recipient_decomposition_native_residual: float
    donor_decomposition_native_residual: float
    boundary_divergence_native_residual: float
    native_work_scale: float
    internal_flow_is_dissipation: bool = False
    internal_flow_creates_event_depth: bool = False
    internal_flow_creates_scale_progress: bool = False

    def __post_init__(self) -> None:
        if not self.modes:
            raise ValueError("nonempty physical helical mode set required")
        for name, value in (
            ("internal flow", self.internal_flow),
            ("inward boundary flow", self.inward_boundary_flow),
            ("outward boundary flow", self.outward_boundary_flow),
            ("recipient positive work", self.recipient_positive_work),
            ("donor negative work", self.donor_negative_work),
            ("native work scale", self.native_work_scale),
            ("recipient decomposition residual", self.recipient_decomposition_native_residual),
            ("donor decomposition residual", self.donor_decomposition_native_residual),
            ("boundary divergence residual", self.boundary_divergence_native_residual),
        ):
            _finite_nonnegative(value, name)
        if not math.isfinite(float(self.signed_nonlinear_work)):
            raise ValueError("finite signed nonlinear work required")
        if self.native_work_scale <= 0.0:
            raise ValueError("positive native physical work scale required")
        if max(
            self.recipient_decomposition_native_residual,
            self.donor_decomposition_native_residual,
            self.boundary_divergence_native_residual,
        ) > 5.0e-10:
            raise AssertionError("helical mode-set boundary law left the native physical work scale")
        if self.internal_flow_is_dissipation or self.internal_flow_creates_event_depth or self.internal_flow_creates_scale_progress:
            raise ValueError("internal same-time energy redistribution may not become dissipation, event depth, or scale progress")


def mode_set_boundary_balance(
    atoms: Sequence[HelicalModeFlowAtom],
    modes: Iterable[HelicalModeIdentity],
) -> HelicalModeSetBoundaryBalance:
    flow = tuple(atoms)
    selected = frozenset(modes)
    if not flow or not selected:
        raise ValueError("nonempty physical flow and mode set required")
    native = math.fsum(atom.physical_work_mass for atom in flow)
    if native <= 0.0:
        raise ValueError("positive physical flow mass required")
    internal = math.fsum(
        atom.physical_work_mass
        for atom in flow
        if atom.donor_mode in selected and atom.recipient_mode in selected
    )
    inward = math.fsum(
        atom.physical_work_mass
        for atom in flow
        if atom.donor_mode not in selected and atom.recipient_mode in selected
    )
    outward = math.fsum(
        atom.physical_work_mass
        for atom in flow
        if atom.donor_mode in selected and atom.recipient_mode not in selected
    )
    positive = math.fsum(atom.physical_work_mass for atom in flow if atom.recipient_mode in selected)
    negative = math.fsum(atom.physical_work_mass for atom in flow if atom.donor_mode in selected)
    signed = positive - negative
    return HelicalModeSetBoundaryBalance(
        modes=selected,
        internal_flow=internal,
        inward_boundary_flow=inward,
        outward_boundary_flow=outward,
        recipient_positive_work=positive,
        donor_negative_work=negative,
        signed_nonlinear_work=signed,
        recipient_decomposition_native_residual=_native_residual(positive, internal + inward, native),
        donor_decomposition_native_residual=_native_residual(negative, internal + outward, native),
        boundary_divergence_native_residual=_native_residual(signed, inward - outward, native),
        native_work_scale=native,
    )


@dataclass(frozen=True)
class HelicalModeEnergyStock:
    mode: HelicalModeIdentity
    energy: float
    wave_number_squared: float

    def __post_init__(self) -> None:
        _finite_nonnegative(self.energy, "helical modal energy")
        _finite_nonnegative(self.wave_number_squared, "helical modal wave-number squared")
        expected = math.fsum(float(v) * float(v) for v in self.mode.wavevector)
        scale = max(1.0, expected, self.wave_number_squared)
        if abs(self.wave_number_squared - expected) > 5.0e-13 * scale:
            raise AssertionError("helical modal stock wave-number changed from its physical mode")


@dataclass(frozen=True)
class HelicalModeSetIntervalContinuity:
    modes: frozenset[HelicalModeIdentity]
    initial_energy: float
    final_energy: float
    integrated_inward_flow: float
    integrated_outward_flow: float
    viscous_dissipation: float
    balance_native_residual: float
    native_energy_throughput_scale: float
    fifo_matching_used: bool = False
    lifo_matching_used: bool = False
    gross_transfer_declared_finite_resource: bool = False
    hard_interaction_cell_used_as_persistent_inventory: bool = False

    def __post_init__(self) -> None:
        if not self.modes:
            raise ValueError("nonempty physical helical mode set required")
        for name, value in (
            ("initial energy", self.initial_energy),
            ("final energy", self.final_energy),
            ("integrated inward flow", self.integrated_inward_flow),
            ("integrated outward flow", self.integrated_outward_flow),
            ("viscous dissipation", self.viscous_dissipation),
            ("balance native residual", self.balance_native_residual),
            ("native energy-throughput scale", self.native_energy_throughput_scale),
        ):
            _finite_nonnegative(value, name)
        if self.native_energy_throughput_scale <= 0.0:
            raise ValueError("positive native energy-throughput scale required")
        # The mathematical continuity law is exact.  This object records a
        # finite-precision / finite-time-quadrature residual but does not impose
        # a universal numerical tolerance: that tolerance belongs to the
        # concrete PDE discretization/audit which produced the observation.
        if (
            self.fifo_matching_used
            or self.lifo_matching_used
            or self.gross_transfer_declared_finite_resource
            or self.hard_interaction_cell_used_as_persistent_inventory
        ):
            raise ValueError("between-time continuity may not invent deposit matching, a gross-transfer budget, or hard-cell inventory")

    @property
    def initial_plus_inward(self) -> float:
        return self.initial_energy + self.integrated_inward_flow

    @property
    def final_plus_outward_plus_viscosity(self) -> float:
        return self.final_energy + self.integrated_outward_flow + self.viscous_dissipation


def interval_continuity_certificate(
    *,
    modes: Iterable[HelicalModeIdentity],
    initial_energy: float,
    final_energy: float,
    integrated_inward_flow: float,
    integrated_outward_flow: float,
    viscous_dissipation: float,
    native_energy_throughput_scale: float | None = None,
) -> HelicalModeSetIntervalContinuity:
    selected = frozenset(modes)
    lhs = float(final_energy) + float(integrated_outward_flow) + float(viscous_dissipation)
    rhs = float(initial_energy) + float(integrated_inward_flow)
    scale = (
        max(abs(lhs), abs(rhs), 1.0e-300)
        if native_energy_throughput_scale is None
        else float(native_energy_throughput_scale)
    )
    return HelicalModeSetIntervalContinuity(
        modes=selected,
        initial_energy=float(initial_energy),
        final_energy=float(final_energy),
        integrated_inward_flow=float(integrated_inward_flow),
        integrated_outward_flow=float(integrated_outward_flow),
        viscous_dissipation=float(viscous_dissipation),
        balance_native_residual=_native_residual(lhs, rhs, scale),
        native_energy_throughput_scale=scale,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "same_time_flow": "cyclic donor measure is pushed to physical helical-mode nodes; for mode set A, P_A=I_A+In_A and N_A=I_A+Out_A, hence P_A-N_A=In_A-Out_A",
        "between_time_continuity": "E_A(t1)+D_A[s,t1]+Phi_out(A;[s,t1]) = E_A(s)+Phi_in(A;[s,t1])",
        "singleton": "for one helical mode m, E_m(t1)+integral W_m^- + 2 nu |k|^2 integral E_m = E_m(s)+integral W_m^+, using gross canonical edge Hahn work, not Hahn of modal net work",
        "full_set": "for the full retained mode set, inward=outward=0 and the identity reduces to the ordinary viscous Navier-Stokes energy balance",
        "internal_flow": "F(AxA) is real same-time nonlinear redistribution and cancels exactly from the set divergence; it is neither dissipation nor event depth",
        "gross_transfer_anti_theorem": "a closed triad wholly inside A can carry arbitrarily rescaled positive internal donor/recipient work while Phi_in=Phi_out=0; the energy stock law is not a finite gross-transfer budget",
        "temporal_provenance": "the identity gives aggregate stock/flow balance only; it does not canonically pair an earlier positive deposit with a later withdrawal and uses no FIFO/LIFO rule",
        "state_ontology": "persistent state energy lives on physical helical modes; hard product cells remain interaction labels and are not wallets across time",
        "capacity_is_causal_law": False,
        "later_hahn_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class HelicalModeSetContinuityStress:
    samples: int
    resolved_cases: int
    unresolved_cases: int
    mode_sets_checked: int
    proper_boundary_cases: int
    full_closed_set_cases: int
    worst_recipient_decomposition_native_residual: float
    worst_donor_decomposition_native_residual: float
    worst_boundary_divergence_native_residual: float
    maximum_internal_to_boundary_ratio_when_boundary_positive: float
    closed_triad_internal_flow_base: float
    closed_triad_internal_flow_scaled: float
    closed_triad_boundary_flux_base: float
    closed_triad_boundary_flux_scaled: float
    gross_transfer_budget_counterexample_passed: bool


def amplitude_scaled_closed_triad(triad, factor: float):
    """Scale one actual closed-triad Fourier state, not its base measure."""
    lam = float(factor)
    if not math.isfinite(lam) or lam <= 0.0:
        raise ValueError("positive finite physical amplitude scaling required")
    return register_closed_helical_triad(
        wavevectors=tuple(mode.wavevector for mode in triad.modes),
        helicities=tuple(mode.helicity for mode in triad.modes),
        amplitudes=tuple(lam * a for a in triad.amplitudes),
    )


def _random_closed_triad(rng: np.random.Generator):
    while True:
        k0 = rng.normal(size=3)
        k1 = rng.normal(size=3)
        k2 = -(k0 + k1)
        if min(stable_norm3(k) for k in (k0, k1, k2)) > 0.08:
            break
    helicities = tuple(int(v) for v in rng.choice((-1, 1), size=3))
    amplitudes = tuple(complex(v) for v in (rng.normal(size=3) + 1j * rng.normal(size=3)))
    return register_closed_helical_triad(
        wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=amplitudes
    )


def stress(samples: int = 75_000, seed: int = 2026081205) -> HelicalModeSetContinuityStress:
    if samples <= 0:
        raise ValueError("positive mode-set continuity stress sample count required")
    rng = np.random.default_rng(seed)
    resolved = unresolved = sets = proper = full = 0
    wr = wd = wb = max_internal_boundary = 0.0
    for _ in range(int(samples)):
        triad = _random_closed_triad(rng)
        kernel = cyclic_triad_measure_kernel(
            triad, quotient_measure_mass=math.exp(float(rng.uniform(-8.0, 8.0)))
        )
        if not kernel.numerically_resolved_transport:
            unresolved += 1
            continue
        atoms = flow_atoms_from_cyclic_kernel(kernel)
        nodes = tuple(sorted({atom.donor_mode for atom in atoms} | {atom.recipient_mode for atom in atoms}))
        resolved += 1
        # Every physical triad is checked on one proper random nonempty subset and on
        # its full node set.  No synthetic graph is inserted.
        count = int(rng.integers(1, len(nodes)))
        chosen = frozenset(rng.choice(len(nodes), size=count, replace=False).tolist())
        proper_modes = frozenset(nodes[int(i)] for i in chosen)
        for selected, is_full in ((proper_modes, False), (frozenset(nodes), True)):
            balance = mode_set_boundary_balance(atoms, selected)
            sets += 1
            full += int(is_full)
            proper += int(not is_full and (balance.inward_boundary_flow + balance.outward_boundary_flow) > 0.0)
            wr = max(wr, balance.recipient_decomposition_native_residual)
            wd = max(wd, balance.donor_decomposition_native_residual)
            wb = max(wb, balance.boundary_divergence_native_residual)
            boundary = balance.inward_boundary_flow + balance.outward_boundary_flow
            if boundary > 0.0:
                max_internal_boundary = max(max_internal_boundary, balance.internal_flow / boundary)
            if is_full and (balance.inward_boundary_flow != 0.0 or balance.outward_boundary_flow != 0.0):
                raise AssertionError("closed physical triad leaked across its full helical-mode set")

    triad, _ = signed_good_integer_triad()
    base_kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    scaled_triad = amplitude_scaled_closed_triad(triad, 10.0)
    scaled_kernel = cyclic_triad_measure_kernel(scaled_triad, quotient_measure_mass=1.0)
    base_atoms = flow_atoms_from_cyclic_kernel(base_kernel)
    scaled_atoms = flow_atoms_from_cyclic_kernel(scaled_kernel)
    all_modes = frozenset({a.donor_mode for a in base_atoms} | {a.recipient_mode for a in base_atoms})
    base = mode_set_boundary_balance(base_atoms, all_modes)
    scaled = mode_set_boundary_balance(scaled_atoms, all_modes)
    anti = (
        base.internal_flow > 0.0
        and scaled.internal_flow > 999.999999 * base.internal_flow
        and base.inward_boundary_flow == 0.0
        and base.outward_boundary_flow == 0.0
        and scaled.inward_boundary_flow == 0.0
        and scaled.outward_boundary_flow == 0.0
        and abs(base.signed_nonlinear_work) <= 5.0e-10 * base.native_work_scale
        and abs(scaled.signed_nonlinear_work) <= 5.0e-10 * scaled.native_work_scale
    )
    if not anti:
        raise AssertionError("physical closed-triad internal circulation did not falsify a gross-transfer budget")

    return HelicalModeSetContinuityStress(
        samples=int(samples),
        resolved_cases=resolved,
        unresolved_cases=unresolved,
        mode_sets_checked=sets,
        proper_boundary_cases=proper,
        full_closed_set_cases=full,
        worst_recipient_decomposition_native_residual=wr,
        worst_donor_decomposition_native_residual=wd,
        worst_boundary_divergence_native_residual=wb,
        maximum_internal_to_boundary_ratio_when_boundary_positive=max_internal_boundary,
        closed_triad_internal_flow_base=base.internal_flow,
        closed_triad_internal_flow_scaled=scaled.internal_flow,
        closed_triad_boundary_flux_base=base.inward_boundary_flow + base.outward_boundary_flow,
        closed_triad_boundary_flux_scaled=scaled.inward_boundary_flow + scaled.outward_boundary_flow,
        gross_transfer_budget_counterexample_passed=anti,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=75_000)
    parser.add_argument("--seed", type=int, default=2026081205)
    parser.add_argument("--outdir", type=Path, default=Path("results-helical-mode-set-energy-continuity"))
    args = parser.parse_args()
    out = stress(samples=args.samples, seed=args.seed)
    cert = theorem_certificate()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "certificate.json").write_text(json.dumps({"theorem": cert, "stress": asdict(out)}, indent=2, sort_keys=True) + "\n")
    summary = f"""# Helical mode-set energy continuity

Status: **{STATUS}**.

The same-time cyclic donor law is read as a positive flow on physical helical-mode nodes.  For every mode set `A`, its internal flow cancels exactly from nonlinear divergence: `P_A=I_A+In_A`, `N_A=I_A+Out_A`, hence `P_A-N_A=In_A-Out_A`.  Combining this with the native modal Navier--Stokes energy equation gives the between-time continuity law `E_A(t1)+D_A+Phi_out=E_A(t0)+Phi_in`.

Stress: `{out.samples}` physical closed triads
- resolved / numerically unresolved: `{out.resolved_cases}` / `{out.unresolved_cases}`
- physical mode sets checked: `{out.mode_sets_checked}`
- proper boundary-flow cases: `{out.proper_boundary_cases}`
- full closed-triad sets: `{out.full_closed_set_cases}`
- worst recipient decomposition native residual: `{out.worst_recipient_decomposition_native_residual:.3e}`
- worst donor decomposition native residual: `{out.worst_donor_decomposition_native_residual:.3e}`
- worst boundary-divergence native residual: `{out.worst_boundary_divergence_native_residual:.3e}`
- maximum sampled internal/boundary flow ratio when boundary was nonzero: `{out.maximum_internal_to_boundary_ratio_when_boundary_positive:.12g}`
- physical closed-triad internal flow, base/scaled: `{out.closed_triad_internal_flow_base:.12g}` / `{out.closed_triad_internal_flow_scaled:.12g}`
- corresponding boundary flux, base/scaled: `{out.closed_triad_boundary_flux_base:.12g}` / `{out.closed_triad_boundary_flux_scaled:.12g}`
- gross-transfer-budget anti-theorem: `{out.gross_transfer_budget_counterexample_passed}`

Internal nonlinear flow is real physical redistribution, not dissipation, event depth, or scale progress.  The interval identity gives aggregate stock/flow conservation only; it creates no FIFO/LIFO matching of prior deposits to later withdrawals and no finite gross-transfer budget.  Persistent inventory lives on physical modes, not hard interaction cells.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
