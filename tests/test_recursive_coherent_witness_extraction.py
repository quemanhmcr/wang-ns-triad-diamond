import math

import numpy as np

from src.coherent_localization_operators import partition_operators, random_parseval_frame
from src.recursive_coherent_witness_extraction import (
    bilinear_partition_residual,
    bilinear_tensor_apply,
    binary_work_ledger,
    coherent_binary_work_atoms,
    excise_positive_xi,
    exact_work_reconstruction_residual,
    normalized_positive_binary_events,
    retained_generation_lower,
)


def _simple_partition(n: int):
    rng = np.random.default_rng(7)
    frame = random_parseval_frame(rng, n, n + 3)
    return partition_operators(frame, [[0, 1], list(range(2, n + 3))])


def test_coherent_parent_partitions_reconstruct_quadratic_source_exactly():
    n = 3
    rng = np.random.default_rng(11)
    B = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
    f = rng.normal(size=n) + 1j * rng.normal(size=n)
    g = rng.normal(size=n) + 1j * rng.normal(size=n)
    ops = _simple_partition(n)
    r = bilinear_partition_residual(B, f, g, ops, ops)
    assert np.linalg.norm(r) < 1e-12 * max(1.0, np.linalg.norm(bilinear_tensor_apply(B, f, g)))


def test_binary_coherent_work_atoms_reconstruct_actual_child_work():
    n = 3
    rng = np.random.default_rng(13)
    B = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
    f = rng.normal(size=n) + 1j * rng.normal(size=n)
    g = rng.normal(size=n) + 1j * rng.normal(size=n)
    h = rng.normal(size=n) + 1j * rng.normal(size=n)
    ops = _simple_partition(n)
    atoms = coherent_binary_work_atoms(B, f, g, h, ops, ops, ops)
    assert abs(exact_work_reconstruction_residual(B, f, g, h, atoms)) < 1e-11


def test_atomic_positive_work_dominates_positive_aggregate_and_defines_binary_law():
    atoms = np.array([[[3.0, -5.0], [4.0, -1.0]]])
    ledger = binary_work_ledger(atoms)
    assert math.isclose(ledger.signed_work, 1.0)
    assert math.isclose(ledger.positive_transfer_mass, 7.0)
    assert math.isclose(ledger.negative_backscatter_mass, 6.0)
    assert ledger.positive_transfer_mass >= ledger.aggregate_positive_work
    events = normalized_positive_binary_events(atoms)
    assert len(events) == 2
    assert math.isclose(sum(e.probability for e in events), 1.0)
    assert {(e.parent1, e.parent2, e.child) for e in events} == {(0, 0, 0), (0, 1, 0)}


def test_xi_excision_is_single_charge_partition_of_positive_work():
    atoms = np.array([[[3.0, -5.0], [4.0, 2.0]]])
    mask = np.array([[[True, False], [False, True]]])
    out = excise_positive_xi(atoms, mask)
    assert math.isclose(out.total_positive_mass, 9.0)
    assert math.isclose(out.xi_positive_mass, 5.0)
    assert math.isclose(out.retained_positive_mass, 4.0)


def test_energy_gate_and_relative_xi_leave_quantitative_binary_generation():
    assert math.isclose(retained_generation_lower(15.0, 0.25), 6.0)
