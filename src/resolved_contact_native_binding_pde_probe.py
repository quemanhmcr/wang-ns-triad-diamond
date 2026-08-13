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
from src.hard_tail_true_upward_supply import (
    deep_upward_resolved_contact_fixture,
    hard_tail_upward_supply_split,
)
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
from src.resolved_contact_native_binding import (
    SignedResolvedKSAtom,
    coalesce_recipient_mixed_cause,
    cover_canonical_mixed_submeasure_by_ks,
    interior_contact_fixture,
    resolved_contact_smooth_binding,
)


STATUS = (
    "ORTHOGONAL_FOURIER_GALERKIN_NS_RESOLVED_CONTACT_NATIVE_BINDING__"
    "SAME_PHYSICAL_DWPLUS_UNDER_TWO_ADMISSIBLE_CUTOFF_REPARTITIONS__"
    "ACTUAL_U_EQUALS_V_PLUS_H_AND_SAME_EDGE_VH_HH__MIXED_KS_ADJOINT_REFEREE"
)

TRANSITION_CORE_FRACTION = 1.0 / 8.0
PLATEAU_CORE_FRACTION = 1.0 / 5.0
SUPPORT_FRACTION = 1.0 / 4.0


def _native_relative(actual: complex | float, expected: complex | float, scale: float) -> float:
    s = float(scale)
    if not math.isfinite(s) or s <= 0.0:
        raise ValueError("positive finite native comparison scale required")
    return abs(complex(actual) - complex(expected)) / s


def _smooth_cutoff_value(radius: float, shell_scale: float, core_fraction: float) -> float:
    r = float(radius)
    M = float(shell_scale)
    a = float(core_fraction)
    if not all(math.isfinite(v) for v in (r, M, a)) or r < 0.0 or M <= 0.0 or not (0.0 <= a < SUPPORT_FRACTION):
        raise ValueError("finite cutoff geometry with 0<=core<1/4 required")
    core = a * M
    outer = SUPPORT_FRACTION * M
    if r <= core:
        return 1.0
    if r >= outer:
        return 0.0
    t = (r - core) / (outer - core)
    left = math.exp(-1.0 / t)
    right = math.exp(-1.0 / (1.0 - t))
    return right / (left + right)


def _smooth_cutoff_grid(radius: np.ndarray, shell_scale: float, core_fraction: float) -> np.ndarray:
    r = np.asarray(radius, dtype=float)
    M = float(shell_scale)
    a = float(core_fraction)
    if M <= 0.0 or not math.isfinite(M) or not (0.0 <= a < SUPPORT_FRACTION):
        raise ValueError("positive shell scale and 0<=core<1/4 required")
    core = a * M
    outer = SUPPORT_FRACTION * M
    q = np.zeros_like(r)
    q[r <= core] = 1.0
    mid = (r > core) & (r < outer)
    if np.any(mid):
        t = (r[mid] - core) / (outer - core)
        left = np.exp(-1.0 / t)
        right = np.exp(-1.0 / (1.0 - t))
        q[mid] = right / (left + right)
    return q


def _mixed_nonlinear_term(
    V_hat: np.ndarray,
    h_hat: np.ndarray,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
) -> np.ndarray:
    """Direct Leray projection of (V.grad)h+(h.grad)V."""
    V = np.fft.ifftn(V_hat, axes=(1, 2, 3)).real
    h = np.fft.ifftn(h_hat, axes=(1, 2, 3)).real
    cross = np.zeros_like(V)
    for component in range(3):
        for direction in range(3):
            dh = np.fft.ifftn(1j * k[direction] * h_hat[component], axes=(0, 1, 2)).real
            dV = np.fft.ifftn(1j * k[direction] * V_hat[component], axes=(0, 1, 2)).real
            cross[component] += V[direction] * dh + h[direction] * dV
    return _leray_dealias(np.fft.fftn(cross, axes=(1, 2, 3)), k, k2, dealias)


def _embed_registered_triad(*, resolution: int, cutoff: int, amplitude: float, base):
    n = int(resolution)
    k, k2, dealias, actual_cutoff = _small_spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("resolved-contact referee changed requested Galerkin cutoff")
    coeff = np.zeros((3, n, n, n), dtype=complex)
    for mode, modal_amplitude in zip(base.modes, base.amplitudes):
        kv = tuple(int(round(v)) for v in mode.wavevector)
        if max(abs(v) for v in kv) > actual_cutoff:
            raise ValueError("fixture mode lies outside requested Galerkin cutoff")
        fiber = helical_basis(np.asarray(kv, dtype=float), mode.helicity)
        value = complex(modal_amplitude) * fiber
        coeff[(slice(None),) + _index(kv, n)] += value
        neg = tuple(-v for v in kv)
        coeff[(slice(None),) + _index(neg, n)] += np.conjugate(value)
    energy = float(np.vdot(coeff, coeff).real)
    if not math.isfinite(energy) or energy <= 0.0:
        raise AssertionError("resolved-contact fixture lost positive Fourier energy")
    coeff *= float(amplitude) / math.sqrt(energy)
    state = _leray_dealias(coeff * float(n**3), k, k2, dealias)
    return state, k, k2, dealias, base


def _resolved_and_uv_parent_indices(triad, atom) -> tuple[int, int]:
    slot = triad.slot_for_closed_mode_index(atom.recipient_closed_mode_index)
    p0, p1 = slot.parent_closed_indices
    r0 = math.sqrt(sum(float(v) ** 2 for v in triad.modes[p0].wavevector))
    r1 = math.sqrt(sum(float(v) ** 2 for v in triad.modes[p1].wavevector))
    return (p0, p1) if r0 <= r1 else (p1, p0)


def _scaled_recipient_edge_work(triad, atom, resolved_factor: float, uv_factor: float):
    r = atom.recipient_closed_mode_index
    low, high = _resolved_and_uv_parent_indices(triad, atom)
    slot = triad.slot_for_closed_mode_index(r)
    factors = {low: float(resolved_factor), high: float(uv_factor)}
    p0, p1 = slot.parent_closed_indices
    xmode, ymode = triad.modes[p0], triad.modes[p1]
    return register_helical_physical_edge(
        x=np.asarray(xmode.wavevector, dtype=float),
        y=np.asarray(ymode.wavevector, dtype=float),
        z=np.asarray(slot.edge_identity.child.wavevector, dtype=float),
        sx=xmode.helicity,
        sy=ymode.helicity,
        sz=triad.modes[r].helicity,
        ax=factors[p0] * triad.amplitudes[p0],
        ay=factors[p1] * triad.amplitudes[p1],
        az=np.conjugate(triad.amplitudes[r]),
    )


def _same_pair_ks_work(triad, atom, resolved_cutoff_value: float) -> tuple[float, float, float, float, float]:
    """Two-high-mode block of the same physical resolved linearized NS operator."""
    r = atom.recipient_closed_mode_index
    low, high = _resolved_and_uv_parent_indices(triad, atom)
    slot = triad.slot_for_closed_mode_index(r)
    kl = np.asarray(triad.modes[low].wavevector, dtype=float)
    kh = np.asarray(triad.modes[high].wavevector, dtype=float)
    z = np.asarray(slot.edge_identity.child.wavevector, dtype=float)
    sl = triad.modes[low].helicity
    sh = triad.modes[high].helicity
    sz = triad.modes[r].helicity
    a_v = float(resolved_cutoff_value) * triad.amplitudes[low]
    a_h = triad.amplitudes[high]
    a_z = np.conjugate(triad.amplitudes[r])

    L_zh = direct_child_source_coefficient(kl, kh, z, sl, sh, sz, a_v, 1.0 + 0.0j)
    L_hz = direct_child_source_coefficient(-kl, z, kh, sl, sz, sh, np.conjugate(a_v), 1.0 + 0.0j)
    K_zh = 0.5 * (L_zh - np.conjugate(L_hz))
    S_zh = 0.5 * (L_zh + np.conjugate(L_hz))
    K_hz = 0.5 * (L_hz - np.conjugate(L_zh))
    S_hz = 0.5 * (L_hz + np.conjugate(L_zh))

    mixed = 2.0 * float(np.real(np.conjugate(a_z) * L_zh * a_h))
    skew = 2.0 * float(np.real(np.conjugate(a_z) * K_zh * a_h))
    strain = 2.0 * float(np.real(np.conjugate(a_z) * S_zh * a_h))
    reverse_skew = 2.0 * float(np.real(np.conjugate(a_h) * K_hz * a_z))
    reverse_strain = 2.0 * float(np.real(np.conjugate(a_h) * S_hz * a_z))
    return mixed, skew, strain, reverse_skew, reverse_strain


def _select_contact_atom(split, expected_low_ratio: float):
    candidates = [atom for atom in split.atoms if atom.resolved_scale_parent_contact]
    if not candidates:
        raise AssertionError("evolved fixture carries no resolved-contact upward atom")
    return min(
        candidates,
        key=lambda atom: abs(min(atom.interaction_parent_radii) / atom.recipient_shell_scale - float(expected_low_ratio)),
    )


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
    minimum_transition_profile_q: float
    maximum_transition_profile_q: float
    minimum_plateau_profile_q: float
    maximum_plateau_profile_q: float
    worst_whole_pde_bilinear_repartition_relative_residual: float
    worst_cutoff_repartition_gauge_residual: float
    worst_high_shell_low_low_relative_mass: float
    worst_same_edge_signed_repartition_native_residual: float
    worst_same_edge_cutoff_scaling_native_residual: float
    worst_ks_to_mixed_edge_native_residual: float
    worst_ks_signed_identity_native_residual: float
    worst_ks_skew_pair_residual: float
    worst_ks_strain_pair_residual: float
    worst_canonical_ks_positive_cover_defect: float
    global_energy_balance_relative_residual: float
    maximum_divergence_relative_to_initial_l2: float


@dataclass(frozen=True)
class ResolvedContactNativeBindingPDEProbe:
    status: str
    boundary_contact_runs: tuple[ResolvedContactGalerkinRun, ...]
    interior_contact_runs: tuple[ResolvedContactGalerkinRun, ...]
    maximum_representation_spread: float


def _run_one(
    *,
    label: str,
    base,
    expected_shell_index: int,
    expected_low_ratio: float,
    resolution: int,
    cutoff: int,
    steps: int,
    duration: float,
    viscosity: float,
    amplitude: float,
    radial_boundary: float,
) -> ResolvedContactGalerkinRun:
    n, count = int(resolution), int(steps)
    horizon, nu, N = float(duration), float(viscosity), float(radial_boundary)
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
    profiles = (TRANSITION_CORE_FRACTION, PLATEAU_CORE_FRACTION)
    qgrids = tuple(_smooth_cutoff_grid(radius_grid, M, a) for a in profiles)
    high_shell = radius_grid > 0.5 * M
    dt = horizon / count
    times = tuple(j * dt for j in range(count + 1))

    energy: list[float] = []
    gradient: list[float] = []
    divergence: list[float] = []
    target_mass: list[float] = []
    q_values = [[], []]
    whole_res: list[float] = []
    gauge_res: list[float] = []
    lowlow_high: list[float] = []
    edge_split: list[float] = []
    edge_scale: list[float] = []
    ks_edge: list[float] = []
    ks_identity: list[float] = []
    ks_skew: list[float] = []
    ks_strain: list[float] = []
    ks_cover_defect: list[float] = []

    for step in range(count + 1):
        weighted = np.sqrt(k2)[None, ...] * state
        energy.append(_spectral_average_inner(state, state, n))
        gradient.append(_spectral_average_inner(weighted, weighted, n))
        divergence.append(_divergence_norm(state, k, n))
        F_u = -_nonlinear_term(state, k, k2, dealias)
        pde_scale = max(float(np.linalg.norm(F_u)), 1.0e-300)
        reconstructed_profiles: list[np.ndarray] = []
        for qgrid in qgrids:
            V = qgrid[None, ...] * state
            h = state - V
            F_v = -_nonlinear_term(V, k, k2, dealias)
            F_h = -_nonlinear_term(h, k, k2, dealias)
            F_mixed = -_mixed_nonlinear_term(V, h, k, k2, dealias)
            reconstructed = F_v + F_mixed + F_h
            reconstructed_profiles.append(reconstructed)
            whole_res.append(float(np.linalg.norm(F_u - reconstructed)) / pde_scale)
            lowlow_high.append(float(np.linalg.norm(F_v[:, high_shell])) / pde_scale)
        gauge_res.append(float(np.linalg.norm(reconstructed_profiles[0] - reconstructed_profiles[1])) / pde_scale)

        triad = _reregister_fixture_triad(state, base)
        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=qmass)
        if not kernel.numerically_resolved_transport:
            raise AssertionError("evolved resolved-contact triad fell below cyclic sign resolution")
        split = hard_tail_upward_supply_split(triad, kernel, boundary=N)
        atom = _select_contact_atom(split, expected_low_ratio)
        if atom.recipient_shell_index != int(expected_shell_index):
            raise AssertionError("target resolved-contact atom changed recipient shell")
        recipient_atoms = tuple(
            row for row in split.atoms
            if row.recipient_closed_mode_index == atom.recipient_closed_mode_index
        )
        if not recipient_atoms:
            raise AssertionError("target recipient lost all canonical donor provenance")
        target_mass.append(math.fsum(row.physical_work_mass for row in recipient_atoms))
        low_radius = min(atom.interaction_parent_radii)
        full_reg = triad.slot_for_closed_mode_index(atom.recipient_closed_mode_index).edge_registration
        native = max(float(full_reg.native_modal_capacity), 1.0e-300)

        for profile_index, core_fraction in enumerate(profiles):
            q = _smooth_cutoff_value(low_radius, M, core_fraction)
            q_values[profile_index].append(q)
            bindings = tuple(
                resolved_contact_smooth_binding(row, resolved_parent_cutoff_value=q)
                for row in recipient_atoms
            )
            binding = next(
                row for row in bindings
                if row.donor_closed_mode_index == atom.donor_closed_mode_index
            )
            recipient_cause = coalesce_recipient_mixed_cause(bindings)
            mixed_reg = _scaled_recipient_edge_work(triad, atom, q, 1.0)
            hh_reg = _scaled_recipient_edge_work(triad, atom, 1.0 - q, 1.0)
            edge_split.append(
                _native_relative(
                    full_reg.signed_child_energy_work,
                    mixed_reg.signed_child_energy_work + hh_reg.signed_child_energy_work,
                    native,
                )
            )
            edge_scale.append(
                max(
                    _native_relative(mixed_reg.signed_child_energy_work, q * full_reg.signed_child_energy_work, native),
                    _native_relative(hh_reg.signed_child_energy_work, (1.0 - q) * full_reg.signed_child_energy_work, native),
                )
            )
            if binding.mixed_vh_submeasure_mass > 0.0:
                I, K, S, Krev, Srev = _same_pair_ks_work(triad, atom, q)
                work_scale = max(native, abs(I), abs(K), abs(S), 1.0e-300)
                ks_edge.append(abs(I - mixed_reg.signed_child_energy_work) / work_scale)
                ks_identity.append(abs(I - K - S) / work_scale)
                ks_skew.append(abs(K + Krev) / work_scale)
                ks_strain.append(abs(S - Srev) / work_scale)
                cover = cover_canonical_mixed_submeasure_by_ks(
                    recipient_cause,
                    SignedResolvedKSAtom(I, K, S),
                )
                ks_cover_defect.append(max(0.0, -cover.positive_cover_margin) / work_scale)

        if step < count:
            state = _triad_subspace_rk4_step(state, dt, nu, k, k2, dealias, subspace)

    if len(target_mass) != count + 1:
        raise AssertionError("target canonical upward branch changed sign during the short PDE referee")
    e0 = energy[0]
    balance = abs(energy[-1] - e0 + 2.0 * nu * _trapezoid(gradient, times)) / e0
    max_div = max(divergence) / math.sqrt(e0)

    if label == "M4_boundary_contact":
        if max(q_values[0]) > 5.0e-12 or max(q_values[1]) > 5.0e-12:
            raise AssertionError("M/4 boundary contact acquired false resolved mixed weight")
    elif label == "M4_interior_contact":
        if not (min(q_values[0]) > 1.0e-4 and max(q_values[0]) < 1.0 - 1.0e-4):
            raise AssertionError("transition profile failed to realize genuine 0<q<1 on interior contact")
        if min(q_values[1]) < 1.0 - 5.0e-12:
            raise AssertionError("plateau profile failed to resolve the same interior parent fully")
        if not ks_identity:
            raise AssertionError("interior resolved-contact mixed work did not reach K/S referee")
    else:
        raise ValueError("unknown resolved-contact PDE referee label")

    if max(whole_res) > 2.0e-12 or max(gauge_res) > 3.0e-12:
        raise AssertionError("actual u=V+h repartition lost cutoff-independent Navier--Stokes source")
    if max(lowlow_high) > 2.0e-11:
        raise AssertionError("resolved low-low source leaked into the recipient high shell")
    if max(edge_split) > 5.0e-10 or max(edge_scale) > 5.0e-10:
        raise AssertionError("same-edge V/h signed work repartition left native helical work scale")
    if max(ks_edge, default=0.0) > 5.0e-10:
        raise AssertionError("same resolved K/S operator failed to reconstruct mixed edge work")
    if max(ks_identity, default=0.0) > 5.0e-10 or max(ks_skew, default=0.0) > 5.0e-10 or max(ks_strain, default=0.0) > 5.0e-10:
        raise AssertionError("actual resolved mixed operator lost K/S adjoint structure")
    if max(ks_cover_defect, default=0.0) > 5.0e-10:
        raise AssertionError("canonical mixed dW+ escaped the existing K+/S+ positive cover")
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
        target_snapshots=len(target_mass),
        total_snapshots=count + 1,
        initial_canonical_target_mass=target_mass[0],
        final_canonical_target_mass=target_mass[-1],
        minimum_transition_profile_q=min(q_values[0]),
        maximum_transition_profile_q=max(q_values[0]),
        minimum_plateau_profile_q=min(q_values[1]),
        maximum_plateau_profile_q=max(q_values[1]),
        worst_whole_pde_bilinear_repartition_relative_residual=max(whole_res),
        worst_cutoff_repartition_gauge_residual=max(gauge_res),
        worst_high_shell_low_low_relative_mass=max(lowlow_high),
        worst_same_edge_signed_repartition_native_residual=max(edge_split),
        worst_same_edge_cutoff_scaling_native_residual=max(edge_scale),
        worst_ks_to_mixed_edge_native_residual=max(ks_edge, default=0.0),
        worst_ks_signed_identity_native_residual=max(ks_identity, default=0.0),
        worst_ks_skew_pair_residual=max(ks_skew, default=0.0),
        worst_ks_strain_pair_residual=max(ks_strain, default=0.0),
        worst_canonical_ks_positive_cover_defect=max(ks_cover_defect, default=0.0),
        global_energy_balance_relative_residual=balance,
        maximum_divergence_relative_to_initial_l2=max_div,
    )


def _relative_spread(values: Sequence[float]) -> float:
    vals = tuple(float(v) for v in values)
    scale = max(max(abs(v) for v in vals), 1.0e-300)
    return (max(vals) - min(vals)) / scale


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    boundary_cutoff: int = 2,
    interior_cutoff: int = 3,
    steps: int = 24,
    duration: float = 2.0e-4,
    viscosity: float = 0.05,
    amplitude: float = 1.0,
) -> ResolvedContactNativeBindingPDEProbe:
    boundary_base, _, boundary_split = deep_upward_resolved_contact_fixture()
    boundary_atom = min(
        (a for a in boundary_split.atoms if a.resolved_scale_parent_contact),
        key=lambda a: abs(min(a.interaction_parent_radii) / a.recipient_shell_scale - 0.25),
    )
    boundary_ratio = min(boundary_atom.interaction_parent_radii) / boundary_atom.recipient_shell_scale

    interior_base, _, interior_split = interior_contact_fixture()
    interior_atom = next(
        a
        for a in interior_split.atoms
        if a.resolved_scale_parent_contact
        and a.recipient_shell_index == 2
        and 0.125 < min(a.interaction_parent_radii) / a.recipient_shell_scale < 0.2
    )
    interior_ratio = min(interior_atom.interaction_parent_radii) / interior_atom.recipient_shell_scale

    boundary_runs = tuple(
        _run_one(
            label="M4_boundary_contact",
            base=boundary_base,
            expected_shell_index=boundary_atom.recipient_shell_index,
            expected_low_ratio=boundary_ratio,
            resolution=int(n),
            cutoff=int(boundary_cutoff),
            steps=int(steps),
            duration=float(duration),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            radial_boundary=boundary_atom.boundary,
        )
        for n in resolutions
    )
    interior_runs = tuple(
        _run_one(
            label="M4_interior_contact",
            base=interior_base,
            expected_shell_index=interior_atom.recipient_shell_index,
            expected_low_ratio=interior_ratio,
            resolution=int(n),
            cutoff=int(interior_cutoff),
            steps=int(steps),
            duration=float(duration),
            viscosity=float(viscosity),
            amplitude=float(amplitude),
            radial_boundary=interior_atom.boundary,
        )
        for n in resolutions
    )
    spread = max(
        _relative_spread([r.initial_canonical_target_mass for r in boundary_runs]),
        _relative_spread([r.final_canonical_target_mass for r in boundary_runs]),
        _relative_spread([r.initial_canonical_target_mass for r in interior_runs]),
        _relative_spread([r.final_canonical_target_mass for r in interior_runs]),
    )
    if spread > 5.0e-10:
        raise AssertionError("same six-mode Galerkin physics changed across FFT representations")
    return ResolvedContactNativeBindingPDEProbe(
        status=STATUS,
        boundary_contact_runs=boundary_runs,
        interior_contact_runs=interior_runs,
        maximum_representation_spread=spread,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--resolutions", type=int, nargs="+", default=[24, 28])
    ap.add_argument("--boundary-cutoff", type=int, default=2)
    ap.add_argument("--interior-cutoff", type=int, default=3)
    ap.add_argument("--steps", type=int, default=24)
    ap.add_argument("--duration", type=float, default=2.0e-4)
    ap.add_argument("--viscosity", type=float, default=0.05)
    ap.add_argument("--amplitude", type=float, default=1.0)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-contact-native-binding-pde"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = run_probe(
        resolutions=args.resolutions,
        boundary_cutoff=args.boundary_cutoff,
        interior_cutoff=args.interior_cutoff,
        steps=args.steps,
        duration=args.duration,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
    )
    payload = asdict(out)
    (args.outdir / "resolved_contact_native_binding_pde.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    all_runs = (*out.boundary_contact_runs, *out.interior_contact_runs)
    md = f"""# Actual Fourier--Galerkin Navier--Stokes referee: resolved-contact native binding\n\nStatus: **{STATUS}**.\n\nThe referee evolves the same real divergence-free six-mode Navier--Stokes data independently of the analysis cutoff.  At every snapshot it reconstructs canonical cyclic `dW+`, then reads two different admissible smooth `u=V+h` decompositions of that same physical interaction.  Both decompositions must reconstruct the identical full NS nonlinear source.  On interior contact one profile leaves genuine mixed+HH transition while the other resolves the same low parent completely; canonical cause and donor provenance do not change.\n\n- maximum representation spread: `{out.maximum_representation_spread:.3e}`\n- worst full PDE bilinear repartition residual: `{max(r.worst_whole_pde_bilinear_repartition_relative_residual for r in all_runs):.3e}`\n- worst cutoff-repartition gauge residual: `{max(r.worst_cutoff_repartition_gauge_residual for r in all_runs):.3e}`\n- worst high-shell low-low leakage: `{max(r.worst_high_shell_low_low_relative_mass for r in all_runs):.3e}`\n- worst same-edge signed V/h repartition residual: `{max(r.worst_same_edge_signed_repartition_native_residual for r in all_runs):.3e}`\n- worst K/S signed identity residual: `{max(r.worst_ks_signed_identity_native_residual for r in all_runs):.3e}`\n- worst K skew-pair residual: `{max(r.worst_ks_skew_pair_residual for r in all_runs):.3e}`\n- worst S symmetric-pair residual: `{max(r.worst_ks_strain_pair_residual for r in all_runs):.3e}`\n- worst canonical K/S positive-cover defect: `{max(r.worst_canonical_ks_positive_cover_defect for r in all_runs):.3e}`\n- worst Galerkin energy-balance residual: `{max(r.global_energy_balance_relative_residual for r in all_runs):.3e}`\n\nThis is a referee for the physical identities and type barriers, not a substitute for the continuum proof and not a Navier--Stokes regularity claim.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
