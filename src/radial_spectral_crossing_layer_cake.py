from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.continuum_helical_edge_measure_registration import HelicalModeIdentity
from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    register_closed_helical_triad,
    signed_good_integer_triad,
)
from src.helical import coupling_g, stable_norm3
from src.helical_mode_set_energy_continuity import (
    HelicalModeFlowAtom,
    flow_atoms_from_cyclic_kernel,
    mode_set_boundary_balance,
)

STATUS = (
    "EXACT_RADIAL_SPECTRAL_CROSSING_LAYER_CAKE__PHYSICAL_HELICAL_DONOR_RECIPIENT_RADII__"
    "TAIL_INWARD_OUTWARD_FLUX__TRUNCATED_LOG_RADIUS_LAYER_CAKE__"
    "NO_REHAHN_DYADIC_PROGRESS_OR_GROSS_CROSSING_BUDGET"
)


def _finite_positive(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"positive finite {name} required")
    return out


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


def _native_residual(actual: float, expected: float, scale: float) -> float:
    a, b, s = float(actual), float(expected), float(scale)
    if not all(math.isfinite(v) for v in (a, b, s)) or s < 0.0:
        raise ValueError("finite nonnegative native comparison scale required")
    gap = abs(a - b)
    if s == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / s


def mode_radius(mode: HelicalModeIdentity) -> float:
    return stable_norm3(np.asarray(mode.wavevector, dtype=float))


def clipped_log_radius_potential(radius: float, lower: float, upper: float) -> float:
    rho = _finite_positive(radius, "mode radius")
    lo = _finite_positive(lower, "lower radial boundary")
    hi = _finite_positive(upper, "upper radial boundary")
    if not lo < hi:
        raise ValueError("truncated radial layer cake requires lower < upper")
    clipped = min(max(rho, lo), hi)
    return math.log(clipped / lo)


@dataclass(frozen=True)
class RadialSpectralCrossingBalance:
    radius: float
    low_internal_flow: float
    high_internal_flow: float
    upward_crossing_flow: float
    downward_crossing_flow: float
    tail_positive_work: float
    tail_negative_work: float
    tail_signed_work: float
    tail_mode_set_signed_work: float
    partition_native_residual: float
    tail_positive_native_residual: float
    tail_negative_native_residual: float
    tail_divergence_native_residual: float
    native_work_scale: float
    later_hahn_used: bool = False
    crossing_creates_event_depth: bool = False
    minimum_scale_progress_claimed: bool = False
    gross_crossing_declared_finite_resource: bool = False
    hard_shell_reweighting_used: bool = False

    def __post_init__(self) -> None:
        _finite_positive(self.radius, "radial boundary")
        for name, value in (
            ("low internal flow", self.low_internal_flow),
            ("high internal flow", self.high_internal_flow),
            ("upward crossing flow", self.upward_crossing_flow),
            ("downward crossing flow", self.downward_crossing_flow),
            ("tail positive work", self.tail_positive_work),
            ("tail negative work", self.tail_negative_work),
            ("native work scale", self.native_work_scale),
            ("partition residual", self.partition_native_residual),
            ("tail positive residual", self.tail_positive_native_residual),
            ("tail negative residual", self.tail_negative_native_residual),
            ("tail divergence residual", self.tail_divergence_native_residual),
        ):
            _finite_nonnegative(value, name)
        if not all(math.isfinite(float(v)) for v in (self.tail_signed_work, self.tail_mode_set_signed_work)):
            raise ValueError("finite signed tail work required")
        if self.native_work_scale <= 0.0:
            raise ValueError("positive native radial work scale required")
        if max(
            self.partition_native_residual,
            self.tail_positive_native_residual,
            self.tail_negative_native_residual,
            self.tail_divergence_native_residual,
        ) > 5.0e-10:
            raise AssertionError("radial spectral crossing left the native physical work scale")
        if (
            self.later_hahn_used
            or self.crossing_creates_event_depth
            or self.minimum_scale_progress_claimed
            or self.gross_crossing_declared_finite_resource
            or self.hard_shell_reweighting_used
        ):
            raise ValueError(
                "radial crossing may not create a later Hahn law, event depth, minimum progress, finite traffic budget, or shell reweighting"
            )


def radial_exterior_balance(
    atoms: Sequence[HelicalModeFlowAtom], *, radius: float
) -> RadialSpectralCrossingBalance:
    flow = tuple(atoms)
    R = _finite_positive(radius, "radial boundary")
    if not flow:
        raise ValueError("nonempty physical donor flow required")
    native = math.fsum(atom.physical_work_mass for atom in flow)
    if native <= 0.0:
        raise ValueError("positive physical donor flow mass required")

    def high(mode: HelicalModeIdentity) -> bool:
        return mode_radius(mode) > R

    low_internal = math.fsum(
        a.physical_work_mass for a in flow if not high(a.donor_mode) and not high(a.recipient_mode)
    )
    high_internal = math.fsum(
        a.physical_work_mass for a in flow if high(a.donor_mode) and high(a.recipient_mode)
    )
    upward = math.fsum(
        a.physical_work_mass for a in flow if not high(a.donor_mode) and high(a.recipient_mode)
    )
    downward = math.fsum(
        a.physical_work_mass for a in flow if high(a.donor_mode) and not high(a.recipient_mode)
    )
    tail_positive = high_internal + upward
    tail_negative = high_internal + downward
    tail_signed = upward - downward
    high_modes = frozenset(
        mode
        for atom in flow
        for mode in (atom.donor_mode, atom.recipient_mode)
        if high(mode)
    )
    if high_modes:
        mode_balance = mode_set_boundary_balance(flow, high_modes)
        set_signed = mode_balance.signed_nonlinear_work
    else:
        set_signed = 0.0
    partition = low_internal + high_internal + upward + downward
    return RadialSpectralCrossingBalance(
        radius=R,
        low_internal_flow=low_internal,
        high_internal_flow=high_internal,
        upward_crossing_flow=upward,
        downward_crossing_flow=downward,
        tail_positive_work=tail_positive,
        tail_negative_work=tail_negative,
        tail_signed_work=tail_signed,
        tail_mode_set_signed_work=set_signed,
        partition_native_residual=_native_residual(partition, native, native),
        tail_positive_native_residual=_native_residual(
            tail_positive,
            math.fsum(a.physical_work_mass for a in flow if high(a.recipient_mode)),
            native,
        ),
        tail_negative_native_residual=_native_residual(
            tail_negative,
            math.fsum(a.physical_work_mass for a in flow if high(a.donor_mode)),
            native,
        ),
        tail_divergence_native_residual=_native_residual(tail_signed, set_signed, native),
        native_work_scale=native,
    )


@dataclass(frozen=True)
class TruncatedRadialLayerCake:
    lower_radius: float
    upper_radius: float
    upward_log_action: float
    downward_log_action: float
    signed_log_action: float
    recipient_log_moment: float
    donor_log_moment: float
    marginal_log_difference: float
    upward_atomwise_identity_native_residual: float
    downward_atomwise_identity_native_residual: float
    signed_marginal_identity_native_residual: float
    native_log_action_scale: float
    dlog_radius_measure_used: bool = True
    existing_cyclic_flow_restricted: bool = True
    later_hahn_used: bool = False
    identified_with_single_edge_young_progress: bool = False
    gross_log_action_declared_finite_resource: bool = False

    def __post_init__(self) -> None:
        lo = _finite_positive(self.lower_radius, "lower radial boundary")
        hi = _finite_positive(self.upper_radius, "upper radial boundary")
        if not lo < hi:
            raise ValueError("truncated radial layer cake requires lower < upper")
        for name, value in (
            ("upward log action", self.upward_log_action),
            ("downward log action", self.downward_log_action),
            ("recipient log moment", self.recipient_log_moment),
            ("donor log moment", self.donor_log_moment),
            ("native log-action scale", self.native_log_action_scale),
            ("upward identity residual", self.upward_atomwise_identity_native_residual),
            ("downward identity residual", self.downward_atomwise_identity_native_residual),
            ("signed marginal identity residual", self.signed_marginal_identity_native_residual),
        ):
            _finite_nonnegative(value, name)
        if not all(math.isfinite(float(v)) for v in (self.signed_log_action, self.marginal_log_difference)):
            raise ValueError("finite signed radial log action required")
        if self.native_log_action_scale <= 0.0:
            raise ValueError("positive native radial log-action scale required")
        if max(
            self.upward_atomwise_identity_native_residual,
            self.downward_atomwise_identity_native_residual,
            self.signed_marginal_identity_native_residual,
        ) > 5.0e-10:
            raise AssertionError("truncated radial layer cake left its native physical log-action scale")
        if not self.dlog_radius_measure_used or not self.existing_cyclic_flow_restricted:
            raise ValueError("radial layer cake must use dR/R on the already-existing cyclic flow")
        if self.later_hahn_used or self.identified_with_single_edge_young_progress or self.gross_log_action_declared_finite_resource:
            raise ValueError("radial log action may not become a later Hahn law, Young progress, or finite traffic budget")


def truncated_radial_layer_cake(
    atoms: Sequence[HelicalModeFlowAtom], *, lower_radius: float, upper_radius: float
) -> TruncatedRadialLayerCake:
    flow = tuple(atoms)
    lo = _finite_positive(lower_radius, "lower radial boundary")
    hi = _finite_positive(upper_radius, "upper radial boundary")
    if not lo < hi:
        raise ValueError("truncated radial layer cake requires lower < upper")
    if not flow:
        raise ValueError("nonempty physical donor flow required")
    width = math.log(hi / lo)
    flow_mass = math.fsum(a.physical_work_mass for a in flow)
    if flow_mass <= 0.0:
        raise ValueError("positive physical flow mass required")

    upward_terms: list[float] = []
    downward_terms: list[float] = []
    recipient_terms: list[float] = []
    donor_terms: list[float] = []
    for atom in flow:
        pd = clipped_log_radius_potential(mode_radius(atom.donor_mode), lo, hi)
        pr = clipped_log_radius_potential(mode_radius(atom.recipient_mode), lo, hi)
        delta = pr - pd
        mass = atom.physical_work_mass
        upward_terms.append(mass * max(0.0, delta))
        downward_terms.append(mass * max(0.0, -delta))
        recipient_terms.append(mass * pr)
        donor_terms.append(mass * pd)
    upward = math.fsum(upward_terms)
    downward = math.fsum(downward_terms)
    recipient = math.fsum(recipient_terms)
    donor = math.fsum(donor_terms)
    signed = upward - downward
    marginal = recipient - donor
    native = max(flow_mass * width, upward + downward, abs(signed), 1.0e-300)

    # Atomwise indicator integration over dR/R has the same clipped-log length.
    # The explicit expressions above are the analytic integral, so these residuals
    # defend the stored identities rather than numerically quadraturing R.
    upward_expected = math.fsum(
        atom.physical_work_mass
        * max(
            0.0,
            clipped_log_radius_potential(mode_radius(atom.recipient_mode), lo, hi)
            - clipped_log_radius_potential(mode_radius(atom.donor_mode), lo, hi),
        )
        for atom in flow
    )
    downward_expected = math.fsum(
        atom.physical_work_mass
        * max(
            0.0,
            clipped_log_radius_potential(mode_radius(atom.donor_mode), lo, hi)
            - clipped_log_radius_potential(mode_radius(atom.recipient_mode), lo, hi),
        )
        for atom in flow
    )
    return TruncatedRadialLayerCake(
        lower_radius=lo,
        upper_radius=hi,
        upward_log_action=upward,
        downward_log_action=downward,
        signed_log_action=signed,
        recipient_log_moment=recipient,
        donor_log_moment=donor,
        marginal_log_difference=marginal,
        upward_atomwise_identity_native_residual=_native_residual(upward, upward_expected, native),
        downward_atomwise_identity_native_residual=_native_residual(downward, downward_expected, native),
        signed_marginal_identity_native_residual=_native_residual(signed, marginal, native),
        native_log_action_scale=native,
    )


@dataclass(frozen=True)
class FullFiniteRadialLogAction:
    upward_log_action: float
    downward_log_action: float
    signed_log_action: float
    recipient_log_moment: float
    donor_log_moment: float
    signed_marginal_identity_native_residual: float
    native_log_action_scale: float
    finite_atom_log_moment: bool = True
    continuum_extension_requires_log_moment: bool = True

    def __post_init__(self) -> None:
        for name, value in (
            ("upward full log action", self.upward_log_action),
            ("downward full log action", self.downward_log_action),
            ("native full log-action scale", self.native_log_action_scale),
            ("signed full identity residual", self.signed_marginal_identity_native_residual),
        ):
            _finite_nonnegative(value, name)
        if not all(math.isfinite(float(v)) for v in (
            self.signed_log_action, self.recipient_log_moment, self.donor_log_moment
        )):
            raise ValueError("finite signed full radial log action and log moments required")
        if self.native_log_action_scale <= 0.0 or not self.finite_atom_log_moment or not self.continuum_extension_requires_log_moment:
            raise ValueError("full radial action must remain finite-atom or explicitly log-moment conditioned")


def finite_radial_log_action(atoms: Sequence[HelicalModeFlowAtom]) -> FullFiniteRadialLogAction:
    flow = tuple(atoms)
    if not flow:
        raise ValueError("nonempty finite physical donor flow required")
    upward = downward = recipient = donor = 0.0
    action_scale = 0.0
    for atom in flow:
        rd = mode_radius(atom.donor_mode)
        rr = mode_radius(atom.recipient_mode)
        delta = math.log(rr / rd)
        m = atom.physical_work_mass
        upward += m * max(0.0, delta)
        downward += m * max(0.0, -delta)
        recipient += m * math.log(rr)
        donor += m * math.log(rd)
        action_scale += m * abs(delta)
    signed = upward - downward
    marginal = recipient - donor
    native = max(action_scale, abs(signed), abs(marginal), 1.0e-300)
    return FullFiniteRadialLogAction(
        upward_log_action=upward,
        downward_log_action=downward,
        signed_log_action=signed,
        recipient_log_moment=recipient,
        donor_log_moment=donor,
        signed_marginal_identity_native_residual=_native_residual(signed, marginal, native),
        native_log_action_scale=native,
    )


def equiradial_physical_transfer_triad():
    """A physical closed triad with nonzero energy transfer but zero radial displacement."""
    root3 = math.sqrt(3.0)
    k0 = np.asarray((1.0, 0.0, 0.0))
    k1 = np.asarray((-0.5, 0.5 * root3, 0.0))
    k2 = -(k0 + k1)
    helicities = (1, 1, -1)
    g = coupling_g(k1, k2, -k0, helicities[1], helicities[2], helicities[0])
    if abs(g) == 0.0:
        raise AssertionError("equiradial physical triad unexpectedly lost Waleffe coupling")
    amplitudes = (1.0, g / abs(g), 1.0)
    triad = register_closed_helical_triad(
        wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=amplitudes
    )
    if not triad.donor_kernel.numerically_resolved_transport:
        raise AssertionError("equiradial physical transfer fell below numerical sign resolution")
    return triad


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "radial_exterior": "H_R={(k,s):|k|>R}; existing cyclic donor flow restricted low->high is Phi_up(R), high->low is Phi_down(R), and W_HR^+-W_HR^-=Phi_up-Phi_down with high->high traffic cancelling exactly",
        "tail_continuity": "E_>R(t1)+D_>R+integral Phi_down = E_>R(t0)+integral Phi_up, inherited from the certified mode-set NS continuity law",
        "truncated_layer_cake": "for 0<R0<R1, integral_R0^R1 Phi_up(R)dR/R equals the existing donor-flow mass weighted by the positive clipped log-radius displacement; downward is the negative displacement, and signed action is recipient clipped-log moment minus donor clipped-log moment",
        "finite_full_layer_cake": "for finite physical flow atoms, integral_0^infty (Phi_up-Phi_down)dR/R = sum M(d->r) log(|k_r|/|k_d|); a continuum infinite-range version requires an explicit finite log moment",
        "crossing_semantics": "crossing direction is real Fourier-radius geometry; it adds no event depth and does not imply a minimum dyadic or positive logarithmic step",
        "equiradial_anti_theorem": "a regular closed helical triad can carry nonzero donor-recipient energy flow between distinct equal-radius modes, so physical transfer can have exactly zero radial displacement",
        "young_progress_distinction": "donor-recipient log-radius displacement is not the one-edge quantity log_+(|child|/max(|parent1|,|parent2|)) and is never identified with J_e or Young/Christ saturation",
        "causal_semantics": "Phi_up and Phi_down are positive restrictions of the already-certified cyclic donor disintegration of canonical dW-/dW+; no tail-net Hahn law is minted",
        "gross_crossing_budget": False,
        "later_hahn_used": False,
        "hard_shell_reweighting_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class RadialSpectralCrossingStress:
    samples: int
    resolved_cases: int
    unresolved_cases: int
    upward_crossing_cases: int
    downward_crossing_cases: int
    both_direction_cases_across_sampled_radii: int
    worst_partition_native_residual: float
    worst_tail_divergence_native_residual: float
    worst_truncated_layer_cake_native_residual: float
    worst_full_log_marginal_native_residual: float
    worst_uniform_dilation_normalized_action_residual: float
    equiradial_internal_transfer: float
    equiradial_full_radial_log_action: float
    equiradial_zero_progress_counterexample_passed: bool
    radial_vs_edge_progress_mismatch_cases: int


def _random_closed_triad(rng: np.random.Generator):
    while True:
        k0 = rng.normal(size=3)
        k1 = rng.normal(size=3)
        k2 = -(k0 + k1)
        if min(stable_norm3(k) for k in (k0, k1, k2)) > 0.08:
            break
    helicities = tuple(int(v) for v in rng.choice((-1, 1), size=3))
    amps = tuple(complex(v) for v in (rng.normal(size=3) + 1j * rng.normal(size=3)))
    return register_closed_helical_triad(wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=amps)


def _dilated_triad(triad, factor: float):
    lam = _finite_positive(factor, "wavevector dilation")
    return register_closed_helical_triad(
        wavevectors=tuple(lam * np.asarray(m.wavevector, float) for m in triad.modes),
        helicities=tuple(m.helicity for m in triad.modes),
        amplitudes=triad.amplitudes,
    )


def stress(samples: int = 75_000, seed: int = 2026081206) -> RadialSpectralCrossingStress:
    if samples <= 0:
        raise ValueError("positive radial crossing stress sample count required")
    rng = np.random.default_rng(seed)
    resolved = unresolved = up_cases = down_cases = both_cases = mismatch = 0
    worst_partition = worst_div = worst_layer = worst_full = worst_dilation = 0.0
    for _ in range(samples):
        triad = _random_closed_triad(rng)
        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
        if not kernel.numerically_resolved_transport:
            unresolved += 1
            continue
        resolved += 1
        atoms = flow_atoms_from_cyclic_kernel(kernel)
        radii = sorted({mode_radius(a.donor_mode) for a in atoms} | {mode_radius(a.recipient_mode) for a in atoms})
        if len(radii) < 2:
            continue
        probes = [0.5 * (radii[i] + radii[i + 1]) for i in range(len(radii) - 1)]
        saw_up = saw_down = False
        for R in probes:
            bal = radial_exterior_balance(atoms, radius=R)
            saw_up = saw_up or bal.upward_crossing_flow > 0.0
            saw_down = saw_down or bal.downward_crossing_flow > 0.0
            worst_partition = max(worst_partition, bal.partition_native_residual)
            worst_div = max(worst_div, bal.tail_divergence_native_residual)
        up_cases += int(saw_up)
        down_cases += int(saw_down)
        both_cases += int(saw_up and saw_down)
        lo = 0.8 * min(radii)
        hi = 1.2 * max(radii)
        layer = truncated_radial_layer_cake(atoms, lower_radius=lo, upper_radius=hi)
        full = finite_radial_log_action(atoms)
        worst_layer = max(
            worst_layer,
            layer.upward_atomwise_identity_native_residual,
            layer.downward_atomwise_identity_native_residual,
            layer.signed_marginal_identity_native_residual,
        )
        worst_full = max(worst_full, full.signed_marginal_identity_native_residual)

        lam = 1.7
        dilated = _dilated_triad(triad, lam)
        dk = cyclic_triad_measure_kernel(dilated, quotient_measure_mass=1.0)
        if dk.numerically_resolved_transport:
            dfull = finite_radial_log_action(flow_atoms_from_cyclic_kernel(dk))
            base_up = full.upward_log_action / kernel.total_mass
            base_down = full.downward_log_action / kernel.total_mass
            dil_up = dfull.upward_log_action / dk.total_mass
            dil_down = dfull.downward_log_action / dk.total_mass
            native_per_mass = max(
                (full.upward_log_action + full.downward_log_action) / kernel.total_mass,
                (dfull.upward_log_action + dfull.downward_log_action) / dk.total_mass,
                1.0e-300,
            )
            worst_dilation = max(
                worst_dilation,
                abs(base_up - dil_up) / native_per_mass,
                abs(base_down - dil_down) / native_per_mass,
            )

        for atom in atoms:
            rr = mode_radius(atom.recipient_mode)
            rd = mode_radius(atom.donor_mode)
            radial = math.log(rr / rd)
            # Recipient root edge progress uses both interaction parents, not just
            # this particular energy donor.  Generic mismatch is expected and
            # protects the two observables from accidental identification.
            slot = next(
                s for s in triad.slots if s.edge_identity.child == atom.recipient_mode
            )
            if abs(radial - slot.edge_registration.scale_progress) > 1.0e-8:
                mismatch += 1

    eq = equiradial_physical_transfer_triad()
    eq_kernel = cyclic_triad_measure_kernel(eq, quotient_measure_mass=1.0)
    eq_atoms = flow_atoms_from_cyclic_kernel(eq_kernel)
    eq_full = finite_radial_log_action(eq_atoms)
    eq_internal = math.fsum(a.physical_work_mass for a in eq_atoms)
    anti = eq_internal > 0.0 and abs(eq_full.upward_log_action) < 1e-13 * eq_internal and abs(eq_full.downward_log_action) < 1e-13 * eq_internal
    if not anti:
        raise AssertionError("physical equiradial transfer did not falsify minimum radial progress")
    if mismatch == 0:
        raise AssertionError("radial donor displacement was accidentally identified with recipient edge progress")
    return RadialSpectralCrossingStress(
        samples=int(samples),
        resolved_cases=resolved,
        unresolved_cases=unresolved,
        upward_crossing_cases=up_cases,
        downward_crossing_cases=down_cases,
        both_direction_cases_across_sampled_radii=both_cases,
        worst_partition_native_residual=worst_partition,
        worst_tail_divergence_native_residual=worst_div,
        worst_truncated_layer_cake_native_residual=worst_layer,
        worst_full_log_marginal_native_residual=worst_full,
        worst_uniform_dilation_normalized_action_residual=worst_dilation,
        equiradial_internal_transfer=eq_internal,
        equiradial_full_radial_log_action=eq_full.upward_log_action + eq_full.downward_log_action,
        equiradial_zero_progress_counterexample_passed=anti,
        radial_vs_edge_progress_mismatch_cases=mismatch,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=STATUS)
    parser.add_argument("--samples", type=int, default=75_000)
    parser.add_argument("--seed", type=int, default=2026081206)
    parser.add_argument("--outdir", type=Path, default=Path("results-radial-spectral-crossing"))
    args = parser.parse_args()
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    result = stress(args.samples, args.seed)
    cert = theorem_certificate()
    payload = {"certificate": cert, "stress": asdict(result)}
    (out / "result.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    summary = [
        "# Radial spectral crossing layer cake",
        "",
        f"Status: `{STATUS}`",
        f"- samples: `{result.samples}`",
        f"- resolved/unresolved: `{result.resolved_cases}` / `{result.unresolved_cases}`",
        f"- upward/downward crossing cases: `{result.upward_crossing_cases}` / `{result.downward_crossing_cases}`",
        f"- both directions across sampled radii: `{result.both_direction_cases_across_sampled_radii}`",
        f"- worst radial partition native residual: `{result.worst_partition_native_residual:.3e}`",
        f"- worst tail divergence native residual: `{result.worst_tail_divergence_native_residual:.3e}`",
        f"- worst truncated layer-cake native residual: `{result.worst_truncated_layer_cake_native_residual:.3e}`",
        f"- worst full log-marginal native residual: `{result.worst_full_log_marginal_native_residual:.3e}`",
        f"- worst dilation normalized-action residual: `{result.worst_uniform_dilation_normalized_action_residual:.3e}`",
        f"- equiradial physical transfer / radial action: `{result.equiradial_internal_transfer:.12g}` / `{result.equiradial_full_radial_log_action:.12g}`",
        f"- radial-vs-edge progress mismatch atoms: `{result.radial_vs_edge_progress_mismatch_cases}`",
    ]
    (out / "summary.md").write_text("\n".join(summary) + "\n")


if __name__ == "__main__":
    main()
