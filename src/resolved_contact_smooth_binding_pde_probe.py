from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.continuum_helical_edge_measure_pde_probe import (
    _divergence_norm,
    _index,
    _leray_dealias,
    _nonlinear_term,
    _spectral_average_inner,
    _trapezoid,
)
from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel
from src.cyclic_helical_triad_donor_kernel_pde_probe import _helical_amplitude
from src.hard_tail_true_upward_supply import hard_tail_upward_supply_split
from src.hard_tail_true_upward_supply_pde_probe import (
    _closed_triad_subspace_mask,
    _reregister_fixture_triad,
    _small_spectral_geometry,
    _triad_subspace_rk4_step,
)
from src.helical import helical_basis
from src.helical_physical_edge_registration import (
    direct_child_source_coefficient,
    register_helical_physical_edge,
)
from src.resolved_contact_smooth_binding import (
    SignedResolvedKSAtom,
    bind_canonical_mixed_submeasure_to_ks,
    canonical_positive_resolved_cutoff,
    deep_contact_smooth_repartition,
    strict_deep_resolved_mixed_fixture,
)
from src.hard_tail_true_upward_supply import deep_upward_resolved_contact_fixture


STATUS = (
    "ORTHOGONAL_FOURIER_GALERKIN_NS_RESOLVED_CONTACT_SMOOTH_BINDING__"
    "ACTUAL_U_EQUALS_V_PLUS_H_BILINEAR_REPARTITION__BORDERLINE_M4_TRANSITION_HH__"
    "STRICT_DEEP_M8_MIXED_VH__SAME_EDGE_KS_ADJOINT_REFEREE__CANONICAL_DWPLUS_NONCLONING"
)


def _native_relative(actual: complex | float, expected: complex | float, scale: float) -> float:
    s = float(scale)
    if not math.isfinite(s) or s <= 0.0:
        raise ValueError("positive finite native comparison scale required")
    return abs(complex(actual) - complex(expected)) / s


def _cutoff_grid(radius: np.ndarray, shell_scale: float) -> np.ndarray:
    M = float(shell_scale)
    if not math.isfinite(M) or M <= 0.0:
        raise ValueError("positive finite shell scale required")
    r = np.asarray(radius, dtype=float)
    q = np.zeros_like(r)
    core = M / 8.0
    outer = M / 4.0
    q[r <= core] = 1.0
    mid = (r > core) & (r < outer)
    if np.any(mid):
        t = (r[mid] - core) / (outer - core)
        a = np.exp(-1.0 / t)
        b = np.exp(-1.0 / (1.0 - t))
        q[mid] = b / (a + b)
    return q


def _mixed_nonlinear_term(
    V_hat: np.ndarray,
    h_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Direct Leray projection of (V.grad)h+(h.grad)V, not a subtraction identity."""
    V = np.fft.ifftn(V_hat, axes=(1, 2, 3)).real
    h = np.fft.ifftn(h_hat, axes=(1, 2, 3)).real
    cross = np.zeros_like(V)
    for component in range(3):
        for direction in range(3):
            dh = np.fft.ifftn(
                1j * k[direction] * h_hat[component], axes=(0, 1, 2)
            ).real
            dV = np.fft.ifftn(
                1j * k[direction] * V_hat[component], axes=(0, 1, 2)
            ).real
            cross[component] += V[direction] * dh + h[direction] * dV
    return _leray_dealias(
        np.fft.fftn(cross, axes=(1, 2, 3)), k, k2, dealias
    )


def _embed_registered_triad(*, resolution: int, cutoff: int, amplitude: float, base):
    n = int(resolution)
    k, k2, dealias, actual_cutoff = _small_spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("resolved-contact referee changed requested Galerkin cutoff")
    coeff = np.zeros((3, n, n, n), dtype=complex)
    for mode, a in zip(base.modes, base.amplitudes):
        kv = tuple(int(round(v)) for v in mode.wavevector)
        if max(abs(v) for v in kv) > actual_cutoff:
            raise ValueError("fixture mode lies outside requested dealiased Galerkin cutoff")
        h = helical_basis(np.asarray(kv, dtype=float), mode.helicity)
        value = complex(a) * h
        coeff[(slice(None),) + _index(kv, n)] += value
        neg = tuple(-v for v in kv)
        coeff[(slice(None),) + _index(neg, n)] += np.conjugate(value)
    energy = float(np.vdot(coeff, coeff).real)
    if not math.isfinite(energy) or energy <= 0.0:
        raise AssertionError("resolved-contact fixture lost positive Fourier energy")
    coeff *= float(amplitude) / math.sqrt(energy)
    state = _leray_dealias(coeff * float(n**3), k, k2, dealias)
    return state, k, k2, dealias, base


def _scaled_recipient_edge_work(triad, atom, donor_factor: float, other_factor: float):
    r = atom.recipient_closed_mode_index
    d = atom.donor_closed_mode_index
    slot = triad.slot_for_closed_mode_index(r)
    p0, p1 = slot.parent_closed_indices
    if d not in (p0, p1):
        raise AssertionError("cyclic radial donor is not a recipient interaction parent")
    factors = {d: float(donor_factor), (p1 if p0 == d else p0): float(other_factor)}
    xmode = triad.modes[p0]
    ymode = triad.modes[p1]
    ax = factors[p0] * triad.amplitudes[p0]
    ay = factors[p1] * triad.amplitudes[p1]
    az = np.conjugate(triad.amplitudes[r])
    return register_helical_physical_edge(
        x=np.asarray(xmode.wavevector, dtype=float),
        y=np.asarray(ymode.wavevector, dtype=float),
        z=np.asarray(slot.edge_identity.child.wavevector, dtype=float),
        sx=xmode.helicity,
        sy=ymode.helicity,
        sz=triad.modes[r].helicity,
        ax=ax,
        ay=ay,
        az=az,
    )


def _same_pair_ks_work(triad, atom, donor_cutoff_value: float) -> tuple[float, float, float, float, float]:
    """Actual two-high-mode block of the resolved linearized NS RHS operator."""
    r = atom.recipient_closed_mode_index
    d = atom.donor_closed_mode_index
    slot = triad.slot_for_closed_mode_index(r)
    p0, p1 = slot.parent_closed_indices
    o = p1 if p0 == d else p0
    kd = np.asarray(triad.modes[d].wavevector, dtype=float)
    ko = np.asarray(triad.modes[o].wavevector, dtype=float)
    z = np.asarray(slot.edge_identity.child.wavevector, dtype=float)
    sd = triad.modes[d].helicity
    so = triad.modes[o].helicity
    sz = triad.modes[r].helicity
    ad_v = float(donor_cutoff_value) * triad.amplitudes[d]
    a_o = triad.amplitudes[o]
    a_z = np.conjugate(triad.amplitudes[r])

    # L_{z,o}: one fixed resolved mode kd acts on the high mode ko.
    L_zo = direct_child_source_coefficient(
        kd, ko, z, sd, so, sz, ad_v, 1.0 + 0.0j
    )
    # The adjoint partner L_{o,z} uses the reality partner -kd of the same V.
    L_oz = direct_child_source_coefficient(
        -kd, z, ko, sd, sz, so, np.conjugate(ad_v), 1.0 + 0.0j
    )
    K_zo = 0.5 * (L_zo - np.conjugate(L_oz))
    S_zo = 0.5 * (L_zo + np.conjugate(L_oz))
    K_oz = 0.5 * (L_oz - np.conjugate(L_zo))
    S_oz = 0.5 * (L_oz + np.conjugate(L_zo))

    mixed = 2.0 * float(np.real(np.conjugate(a_z) * L_zo * a_o))
    skew = 2.0 * float(np.real(np.conjugate(a_z) * K_zo * a_o))
    strain = 2.0 * float(np.real(np.conjugate(a_z) * S_zo * a_o))
    reverse_skew = 2.0 * float(np.real(np.conjugate(a_o) * K_oz * a_z))
    reverse_strain = 2.0 * float(np.real(np.conjugate(a_o) * S_oz * a_z))
    return mixed, skew, strain, reverse_skew, reverse_strain


@dataclass(frozen=True)
class ResolvedContactGalerkinRun:
    label: str
    resolution: int
    cutoff: int
    steps: int
    duration: float
    viscosity: float
    radial_boundary: float
    recipient_shell_scale: float
    target_snapshots: int
    total_snapshots: int
    initial_canonical_target_mass: float
    final_canonical_target_mass: float
    minimum_mixed_fraction: float
    maximum_mixed_fraction: float
    minimum_transition_hh_fraction: float
    maximum_transition_hh_fraction: float
    worst_whole_pde_bilinear_repartition_relative_residual: float
    worst_high_shell_low_low_relative_mass: float
    worst_same_edge_signed_repartition_native_residual: float
    worst_same_edge_cutoff_scaling_native_residual: float
    worst_ks_to_mixed_edge_native_residual: float
    worst_ks_signed_identity_native_residual: float
    worst_ks_skew_pair_residual: float
    worst_ks_strain_pair_residual: float
    worst_canonical_ks_mass_residual: float
    worst_canonical_ks_domination_excess: float
    global_energy_balance_relative_residual: float
    maximum_divergence_relative_to_initial_l2: float


@dataclass(frozen=True)
class ResolvedContactSmoothBindingPDEProbe:
    status: str
    borderline_runs: tuple[ResolvedContactGalerkinRun, ...]
    strict_deep_runs: tuple[ResolvedContactGalerkinRun, ...]
    maximum_representation_spread: float


def _run_one(
    *,
    label: str,
    base,
    expected_shell_index: int,
    resolution: int,
    cutoff: int,
    steps: int,
    duration: float,
    viscosity: float,
    amplitude: float,
    radial_boundary: float,
) -> ResolvedContactGalerkinRun:
    n = int(resolution)
    count = int(steps)
    horizon = float(duration)
    nu = float(viscosity)
    N = float(radial_boundary)
    if count < 8 or not all(math.isfinite(v) and v > 0.0 for v in (horizon, nu, float(amplitude), N)):
        raise ValueError("positive finite Galerkin referee parameters and at least eight steps required")
    state, k, k2, dealias, base = _embed_registered_triad(
        resolution=n, cutoff=int(cutoff), amplitude=float(amplitude), base=base
    )
    subspace = _closed_triad_subspace_mask(n, base)
    state = state * subspace[None, ...]
    qmass = 1.0 / unitary_fourier_convolution_factor()
    M = (2.0 ** int(expected_shell_index)) * N
    radius_grid = np.sqrt(k2)
    qgrid = _cutoff_grid(radius_grid, M)
    high_shell = radius_grid > 0.5 * M
    dt = horizon / count
    times = tuple(j * dt for j in range(count + 1))

    energy: list[float] = []
    gradient: list[float] = []
    divergence: list[float] = []
    target_mass: list[float] = []
    mixed_fracs: list[float] = []
    hh_fracs: list[float] = []
    whole_res: list[float] = []
    lowlow_high: list[float] = []
    edge_split: list[float] = []
    edge_scale: list[float] = []
    ks_edge: list[float] = []
    ks_identity: list[float] = []
    ks_skew: list[float] = []
    ks_strain: list[float] = []
    ks_mass: list[float] = []
    ks_dom: list[float] = []

    for step in range(count + 1):
        weighted = np.sqrt(k2)[None, ...] * state
        energy.append(_spectral_average_inner(state, state, n))
        gradient.append(_spectral_average_inner(weighted, weighted, n))
        divergence.append(_divergence_norm(state, k, n))

        # Actual full-state Navier--Stokes bilinear repartition u=V+h.
        V = qgrid[None, ...] * state
        h = state - V
        F_u = -_nonlinear_term(state, k, k2, dealias)
        F_v = -_nonlinear_term(V, k, k2, dealias)
        F_h = -_nonlinear_term(h, k, k2, dealias)
        F_mixed = -_mixed_nonlinear_term(V, h, k, k2, dealias)
        reconstructed = F_v + F_mixed + F_h
        pde_scale = max(float(np.linalg.norm(F_u)), float(np.linalg.norm(F_v)), float(np.linalg.norm(F_h)), 1.0e-300)
        whole_res.append(float(np.linalg.norm(F_u - reconstructed)) / pde_scale)
        high_ll = float(np.linalg.norm(F_v[:, high_shell]))
        lowlow_high.append(high_ll / pde_scale)

        triad = _reregister_fixture_triad(state, base)
        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=qmass)
        if not kernel.numerically_resolved_transport:
            raise AssertionError("evolved resolved-contact triad fell below cyclic sign resolution")
        split = hard_tail_upward_supply_split(triad, kernel, boundary=N)
        targets = [
            atom
            for atom in split.atoms
            if atom.recipient_shell_index == int(expected_shell_index)
            and atom.donor_radius <= N + 5.0e-12 * N
        ]
        if not targets:
            if step == 0:
                raise AssertionError("engineered Galerkin fixture has no target canonical upward atom initially")
        else:
            atom = max(targets, key=lambda a: a.physical_work_mass)
            target_mass.append(atom.physical_work_mass)
            binding = deep_contact_smooth_repartition(atom)
            mixed_fracs.append(binding.mixed_vh_bound_mass / binding.canonical_positive_mass)
            hh_fracs.append(binding.transition_hh_bound_mass / binding.canonical_positive_mass)

            full_reg = triad.slot_for_closed_mode_index(atom.recipient_closed_mode_index).edge_registration
            mixed_reg = _scaled_recipient_edge_work(triad, atom, binding.donor_cutoff_value, 1.0)
            hh_reg = _scaled_recipient_edge_work(triad, atom, 1.0 - binding.donor_cutoff_value, 1.0)
            native = max(float(full_reg.native_modal_capacity), 1.0e-300)
            edge_split.append(
                _native_relative(
                    full_reg.signed_child_energy_work,
                    mixed_reg.signed_child_energy_work + hh_reg.signed_child_energy_work,
                    native,
                )
            )
            edge_scale.append(
                max(
                    _native_relative(
                        mixed_reg.signed_child_energy_work,
                        binding.donor_cutoff_value * full_reg.signed_child_energy_work,
                        native,
                    ),
                    _native_relative(
                        hh_reg.signed_child_energy_work,
                        (1.0 - binding.donor_cutoff_value) * full_reg.signed_child_energy_work,
                        native,
                    ),
                )
            )

            if binding.mixed_vh_bound_mass > 0.0:
                I, K, S, Krev, Srev = _same_pair_ks_work(triad, atom, binding.donor_cutoff_value)
                work_scale = max(native, abs(I), abs(K), abs(S), 1.0e-300)
                ks_edge.append(abs(I - mixed_reg.signed_child_energy_work) / work_scale)
                ks_identity.append(abs(I - K - S) / work_scale)
                ks_skew.append(abs(K + Krev) / work_scale)
                ks_strain.append(abs(S - Srev) / work_scale)
                same = SignedResolvedKSAtom(I, K, S)
                owner = bind_canonical_mixed_submeasure_to_ks(
                    binding.mixed_vh_bound_mass, same, common_unit_scale=N
                )
                ks_mass.append(owner.canonical_mass_residual)
                ks_dom.append(owner.maximum_domination_excess / work_scale)

        if step < count:
            state = _triad_subspace_rk4_step(state, dt, nu, k, k2, dealias, subspace)

    if not target_mass:
        raise AssertionError("target canonical upward branch vanished from every Galerkin snapshot")
    e0 = energy[0]
    balance = abs(energy[-1] - e0 + 2.0 * nu * _trapezoid(gradient, times)) / e0
    max_div = max(divergence) / math.sqrt(e0)
    target_count = len(target_mass)
    if target_count != count + 1:
        raise AssertionError("target canonical upward branch changed sign during the short actual-NS referee")
    if label == "borderline_M4_transition":
        if max(mixed_fracs) > 5.0e-12 or min(hh_fracs) < 1.0 - 5.0e-12:
            raise AssertionError("M=4N boundary-contact actual PDE fixture was incorrectly promoted to mixed/interface work")
    elif label == "strict_deep_M8_mixed":
        if min(mixed_fracs) < 1.0 - 5.0e-12 or max(hh_fracs) > 5.0e-12:
            raise AssertionError("M>=8N actual PDE fixture failed pure mixed V-h binding")
        if not ks_identity:
            raise AssertionError("strict-deep mixed branch did not reach its actual K/S referee")
    else:
        raise ValueError("unknown resolved-contact Galerkin referee label")
    if max(whole_res) > 2.0e-12:
        raise AssertionError("actual u=V+h nonlinear repartition lost the Navier--Stokes bilinear identity")
    if max(lowlow_high) > 2.0e-11:
        raise AssertionError("resolved low-low source leaked into the target high recipient shell")
    if max(edge_split, default=0.0) > 5.0e-10 or max(edge_scale, default=0.0) > 5.0e-10:
        raise AssertionError("same-edge smooth V/h work repartition left native helical work scale")
    if max(ks_edge, default=0.0) > 5.0e-10:
        raise AssertionError("actual resolved operator K/S block did not reconstruct the same mixed physical edge")
    if max(ks_identity, default=0.0) > 5.0e-10 or max(ks_skew, default=0.0) > 5.0e-10 or max(ks_strain, default=0.0) > 5.0e-10:
        raise AssertionError("actual resolved mixed operator lost K/S adjoint structure")
    if max(ks_mass, default=0.0) > 5.0e-10 or max(ks_dom, default=0.0) > 5.0e-10:
        raise AssertionError("canonical mixed dW+ failed non-cloning K/S positive binding")
    if balance > 8.0e-5 or max_div > 5.0e-11:
        raise AssertionError("orthogonal Galerkin referee lost native NS energy/divergence invariants")

    return ResolvedContactGalerkinRun(
        label=label,
        resolution=n,
        cutoff=int(cutoff),
        steps=count,
        duration=horizon,
        viscosity=nu,
        radial_boundary=N,
        recipient_shell_scale=M,
        target_snapshots=target_count,
        total_snapshots=count + 1,
        initial_canonical_target_mass=target_mass[0],
        final_canonical_target_mass=target_mass[-1],
        minimum_mixed_fraction=min(mixed_fracs),
        maximum_mixed_fraction=max(mixed_fracs),
        minimum_transition_hh_fraction=min(hh_fracs),
        maximum_transition_hh_fraction=max(hh_fracs),
        worst_whole_pde_bilinear_repartition_relative_residual=max(whole_res),
        worst_high_shell_low_low_relative_mass=max(lowlow_high),
        worst_same_edge_signed_repartition_native_residual=max(edge_split, default=0.0),
        worst_same_edge_cutoff_scaling_native_residual=max(edge_scale, default=0.0),
        worst_ks_to_mixed_edge_native_residual=max(ks_edge, default=0.0),
        worst_ks_signed_identity_native_residual=max(ks_identity, default=0.0),
        worst_ks_skew_pair_residual=max(ks_skew, default=0.0),
        worst_ks_strain_pair_residual=max(ks_strain, default=0.0),
        worst_canonical_ks_mass_residual=max(ks_mass, default=0.0),
        worst_canonical_ks_domination_excess=max(ks_dom, default=0.0),
        global_energy_balance_relative_residual=balance,
        maximum_divergence_relative_to_initial_l2=max_div,
    )


def _relative_spread(values: Sequence[float]) -> float:
    vals = tuple(float(v) for v in values)
    if not vals:
        return 0.0
    scale = max(max(abs(v) for v in vals), 1.0e-300)
    return (max(vals) - min(vals)) / scale


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    borderline_cutoff: int = 2,
    strict_cutoff: int = 6,
    steps: int = 24,
    duration: float = 2.0e-4,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
) -> ResolvedContactSmoothBindingPDEProbe:
    borderline_base, _, _ = deep_upward_resolved_contact_fixture()
    strict_base, _, _ = strict_deep_resolved_mixed_fixture()
    border = tuple(
        _run_one(
            label="borderline_M4_transition",
            base=borderline_base,
            expected_shell_index=2,
            resolution=int(n),
            cutoff=int(borderline_cutoff),
            steps=int(steps),
            duration=float(duration),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            radial_boundary=1.0,
        )
        for n in resolutions
    )
    strict = tuple(
        _run_one(
            label="strict_deep_M8_mixed",
            base=strict_base,
            expected_shell_index=3,
            resolution=int(n),
            cutoff=int(strict_cutoff),
            steps=int(steps),
            duration=float(duration),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            radial_boundary=1.0,
        )
        for n in resolutions
    )
    spread = max(
        _relative_spread([r.initial_canonical_target_mass for r in border]),
        _relative_spread([r.initial_canonical_target_mass for r in strict]),
        _relative_spread([r.final_canonical_target_mass for r in border]),
        _relative_spread([r.final_canonical_target_mass for r in strict]),
    )
    if spread > 5.0e-10:
        raise AssertionError("same six-mode Galerkin physics changed across FFT representations")
    return ResolvedContactSmoothBindingPDEProbe(
        status=STATUS,
        borderline_runs=border,
        strict_deep_runs=strict,
        maximum_representation_spread=spread,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", type=int, nargs="+", default=[24, 28])
    ap.add_argument("--borderline-cutoff", type=int, default=2)
    ap.add_argument("--strict-cutoff", type=int, default=6)
    ap.add_argument("--steps", type=int, default=24)
    ap.add_argument("--duration", type=float, default=2.0e-4)
    ap.add_argument("--viscosity", type=float, default=0.03)
    ap.add_argument("--amplitude", type=float, default=1.0)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-contact-smooth-binding-pde"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = run_probe(
        resolutions=args.resolutions,
        borderline_cutoff=args.borderline_cutoff,
        strict_cutoff=args.strict_cutoff,
        steps=args.steps,
        duration=args.duration,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
    )
    (args.outdir / "resolved_contact_smooth_binding_pde.json").write_text(
        json.dumps(asdict(out), indent=2, sort_keys=True) + "\n"
    )
    all_runs = (*out.borderline_runs, *out.strict_deep_runs)
    md = f"""# Actual Fourier--Galerkin Navier--Stokes referee: resolved-contact smooth binding\n\nStatus: **{STATUS}**.\n\nThe same real divergence-free six-mode data are evolved on FFT representations `{tuple(args.resolutions)}`. At every snapshot the referee reads the actual cyclic donor law, the actual smooth `u=V+h` decomposition, the same-edge helical work, and—on the strict-deep mixed branch—the adjoint `K/S` pair of that same resolved linearized NS operator.\n\n- maximum representation spread: `{out.maximum_representation_spread:.3e}`\n- worst full PDE bilinear repartition residual: `{max(r.worst_whole_pde_bilinear_repartition_relative_residual for r in all_runs):.3e}`\n- worst high-shell low-low leakage: `{max(r.worst_high_shell_low_low_relative_mass for r in all_runs):.3e}`\n- worst same-edge signed V/h repartition residual: `{max(r.worst_same_edge_signed_repartition_native_residual for r in all_runs):.3e}`\n- worst same-edge cutoff scaling residual: `{max(r.worst_same_edge_cutoff_scaling_native_residual for r in all_runs):.3e}`\n- worst K/S signed identity residual: `{max(r.worst_ks_signed_identity_native_residual for r in all_runs):.3e}`\n- worst K skew-pair residual: `{max(r.worst_ks_skew_pair_residual for r in all_runs):.3e}`\n- worst S symmetric-pair residual: `{max(r.worst_ks_strain_pair_residual for r in all_runs):.3e}`\n- worst canonical K/S mass residual: `{max(r.worst_canonical_ks_mass_residual for r in all_runs):.3e}`\n- worst Galerkin energy-balance residual: `{max(r.global_energy_balance_relative_residual for r in all_runs):.3e}`\n\nThe `M=4N` boundary-contact trajectory remains transition-HH rather than being renamed interface work. The `M=8N` trajectory is fully mixed `V-h` and its donor-restricted canonical `dW+` is bound without cloning into existing skew-relink/strain positive laws only after the signed K/S identity is checked.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(json.dumps(asdict(out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
