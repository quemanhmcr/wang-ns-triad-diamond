import dataclasses
import math

import numpy as np

from src.affine_coherent_moyal import periodic_discrete_stft
from src.canonical_positive_edge_work_routing_pde_probe import _ledger_from_actual_state
from src.coherent_localization_operators import trilinear_tensor_value
from src.continuum_helical_edge_measure_pde_probe import (
    _deterministic_smooth_initial_state,
    _nonlinear_term,
    _snapshot,
    _snapshot_with_ledger,
    _spectral_geometry,
)
from src.continuum_helical_edge_measure_registration import (
    _register_continuum_triad_fiber_with_source,
    leray_project,
    register_continuum_triad_fiber,
    unordered_parent_curl_source_vector,
)
from src.extremal_helicity_symplectic import (
    _transfer_relevant_strain_observability_validated,
    _transfer_relevant_strain_observability_validated_batch,
)
from src.full_strain_observability import tracefree_2x2
from src.helical import helical_basis
from src.helical_phase_holonomy import (
    _diamond_incidence_spin_holonomy_from_edges,
    _diamond_phase_residuals_from_edges,
    diamond_edge_data,
    diamond_incidence_spin_holonomy,
    diamond_phase_residuals,
)
from src.recursive_coherent_witness_extraction import bilinear_tensor_apply


def _reference_helical_basis(k: np.ndarray, s: int) -> np.ndarray:
    q = np.asarray(k, dtype=float)
    norm = float(math.hypot(float(q[0]), float(q[1]), float(q[2])))
    if norm == 0.0:
        raise ValueError("zero wavevector")
    sign = 1
    for value in q:
        if float(value) != 0.0:
            sign = 1 if float(value) > 0.0 else -1
            break
    kp = sign * q
    khat = kp / norm
    axis = np.eye(3)[int(np.argmin(np.abs(khat)))]
    e1 = axis - np.dot(axis, khat) * khat
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(khat, e1)
    h = (e1 + 1j * s * e2) / math.sqrt(2.0)
    return np.conjugate(h) if sign < 0 else h


def test_fixed_three_vector_helical_basis_matches_readable_reference_bitwise():
    rng = np.random.default_rng(2026081201)
    for _ in range(128):
        k = rng.normal(size=3) * math.exp(float(rng.uniform(-120.0, 120.0)))
        s = int(rng.choice((-1, 1)))
        assert np.array_equal(helical_basis(k, s), _reference_helical_basis(k, s))


def test_continuum_registration_source_sidecar_matches_public_reference_bitwise():
    rng = np.random.default_rng(2026081202)
    for _ in range(24):
        while True:
            x = rng.normal(size=3)
            y = rng.normal(size=3)
            z = x + y
            if min(np.linalg.norm(x), np.linalg.norm(y), np.linalg.norm(z)) > 0.2:
                break
        ux = leray_project(x, rng.normal(size=3) + 1j * rng.normal(size=3))
        uy = leray_project(y, rng.normal(size=3) + 1j * rng.normal(size=3))
        uz = leray_project(z, rng.normal(size=3) + 1j * rng.normal(size=3))
        qmass = float(math.exp(float(rng.uniform(-3.0, 3.0))))
        fast, source = _register_continuum_triad_fiber_with_source(
            x=x,
            y=y,
            z=z,
            ux=ux,
            uy=uy,
            uz=uz,
            quotient_measure_mass=qmass,
        )
        reference = register_continuum_triad_fiber(
            x=x,
            y=y,
            z=z,
            ux=ux,
            uy=uy,
            uz=uz,
            quotient_measure_mass=qmass,
        )
        assert dataclasses.asdict(fast) == dataclasses.asdict(reference)
        assert np.array_equal(source, unordered_parent_curl_source_vector(x, y, z, ux, uy))


def test_batched_periodic_stft_matches_reference_loop_bitwise():
    rng = np.random.default_rng(2026081203)
    for n in (8, 13, 21):
        f = rng.normal(size=n) + 1j * rng.normal(size=n)
        g = rng.normal(size=n) + 1j * rng.normal(size=n)
        reference = np.empty((n, n), complex)
        for m in range(n):
            reference[m] = np.fft.fft(f * np.conj(np.roll(g, m)))
        assert np.array_equal(periodic_discrete_stft(f, g), reference)


def test_resolved_bilinear_einsum_path_matches_numpy_optimizer_bitwise():
    rng = np.random.default_rng(2026081204)
    for n in (2, 4, 7):
        B = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
        f = rng.normal(size=n) + 1j * rng.normal(size=n)
        g = rng.normal(size=n) + 1j * rng.normal(size=n)
        reference = np.einsum("kij,i,j->k", B, f, g, optimize=True)
        assert np.array_equal(bilinear_tensor_apply(B, f, g), reference)


def test_resolved_trilinear_einsum_path_matches_numpy_optimizer_bitwise():
    rng = np.random.default_rng(2026081205)
    for n in (2, 4, 7):
        T = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
        f = rng.normal(size=n) + 1j * rng.normal(size=n)
        g = rng.normal(size=n) + 1j * rng.normal(size=n)
        h = rng.normal(size=n) + 1j * rng.normal(size=n)
        reference = np.einsum("ijk,i,j,k->", T, f, g, h, optimize=True)
        assert trilinear_tensor_value(T, f, g, h) == reference


def test_batched_transfer_observability_matches_scalar_reference_bitwise():
    rng = np.random.default_rng(2026081206)
    strains = []
    for _ in range(64):
        X = rng.normal(size=(3, 3))
        S = 0.5 * (X + X.T)
        S -= np.trace(S) / 3.0 * np.eye(3)
        strains.append(S)
    batch = np.stack(strains)
    q_fast, n_fast = _transfer_relevant_strain_observability_validated_batch(batch)
    scalar = [_transfer_relevant_strain_observability_validated(S) for S in strains]
    q_ref = np.array([x[0] for x in scalar])
    n_ref = np.array([x[1] for x in scalar])
    assert np.array_equal(q_fast, q_ref)
    assert np.array_equal(n_fast, n_ref)


def test_tracefree_2x2_matches_numpy_reference_bitwise():
    rng = np.random.default_rng(2026081207)
    for _ in range(128):
        M = rng.normal(size=(2, 2))
        reference = M - 0.5 * np.trace(M) * np.eye(2)
        assert np.array_equal(tracefree_2x2(M), reference)


def test_phase_holonomy_reuses_same_registered_edges_without_changing_result():
    a = np.array([0.9, 0.4, -0.2])
    b = np.array([-0.1, 0.8, 0.3])
    c = np.array([0.2, -0.3, 1.0])
    signs = (1, -1, 1, -1, 1, -1)
    phases = {name: value for name, value in zip(
        ("a", "b", "c", "m", "n", "d"),
        (0.1, -0.7, 1.1, -1.3, 0.4, 2.0),
    )}
    edges = diamond_edge_data(a, b, c, signs)
    assert _diamond_incidence_spin_holonomy_from_edges(a, b, c, signs, edges) == diamond_incidence_spin_holonomy(a, b, c, signs)
    assert _diamond_phase_residuals_from_edges(edges, phases) == diamond_phase_residuals(a, b, c, signs, phases)


def test_actual_ns_snapshot_ledger_reuse_matches_independent_rebuild_exactly():
    n = 20
    k, k2, dealias, cutoff = _spectral_geometry(n, None)
    state = _deterministic_smooth_initial_state(n, k, k2, dealias, 4.0)
    nonlinear = _nonlinear_term(state, k, k2, dealias)
    row, ledger = _snapshot_with_ledger(state, k, k2, dealias, cutoff, nonlinear_hat=nonlinear)
    rebuilt = _ledger_from_actual_state(state, cutoff=cutoff)
    assert dataclasses.asdict(ledger) == dataclasses.asdict(rebuilt)
    assert _snapshot(state, k, k2, dealias, cutoff, nonlinear_hat=nonlinear) == row
