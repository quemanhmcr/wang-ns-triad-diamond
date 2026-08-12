from dataclasses import replace

import numpy as np
import pytest

from src.cyclic_hard_cell_single_charge_quotient import (
    FATE_BAD,
    FATE_GOOD,
    STATUS,
    aggregate_hard_cell_single_charge_quotients,
    fine_hard_role_map,
    hard_cell_single_charge_quotient,
    pushforward_restricted_hard_cell_donor_work,
    required_hard_modes,
    single_hard_role_map,
    theorem_certificate,
)
from src.cyclic_helical_triad_donor_kernel import (
    generic_two_donor_counterexample,
    register_closed_helical_triad,
    signed_good_integer_triad,
)


def _generic_data():
    k0 = np.array([-1.0, -2.0, 0.5])
    k1 = np.array([0.2, 1.3, -0.1])
    k2 = -(k0 + k1)
    s = (1, -1, 1)
    a = (0.7 + 0.4j, -0.3 + 1.1j, 0.9 - 0.2j)
    return (k0, k1, k2), s, a


def _atom_signature(quotient):
    return sorted(
        (
            atom.donor_cell,
            atom.recipient_cell,
            atom.recipient_fate,
            atom.physical_work_mass,
        )
        for atom in quotient.atoms
    )


def test_certificate_preserves_canonical_cause_and_forbids_synthetic_payment_or_time_matching():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "pi_# dW-" in cert["row_marginal"]
    assert "pi_# dW_G+" in cert["column_fate_marginal"]
    assert "one already-canonical recipient charge" in cert["single_charge"]
    assert "does not infer" in cert["failure_anti_shortcut"]
    assert "no FIFO/LIFO" in cert["time_anti_shortcut"]
    assert not cert["capacity_is_causal_law"]
    assert not cert["later_hahn_used"]
    assert not cert["claims_global_regularity"]


def test_signed_good_donor_splits_into_existing_good_and_bad_positive_recipient_causes():
    triad, side = signed_good_integer_triad()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=fine_hard_role_map(triad)
    )
    assert triad.donor_kernel.donor_count == 1
    assert triad.donor_kernel.recipient_count == 2
    assert len(quotient.donor_charges) == 1
    assert len(quotient.recipient_charges) == 2
    assert quotient.good_recipient_mass > 0.0
    assert quotient.bad_recipient_mass > 0.0
    assert quotient.total_negative_work_mass == pytest.approx(quotient.total_positive_work_mass)
    assert quotient.good_recipient_mass + quotient.bad_recipient_mass == pytest.approx(
        quotient.total_positive_work_mass
    )
    fates = {charge.fate: charge for charge in quotient.recipient_charges}
    assert set(fates) == {FATE_GOOD, FATE_BAD}
    assert fates[FATE_BAD].downstream_semantics == "existing_TRANSFER_WORK_LOSS_first_time_None"
    assert fates[FATE_GOOD].downstream_semantics == "existing_geometry_good_Young_eligible_route"
    assert fates[FATE_BAD].canonical_positive_work_mass / fates[FATE_GOOD].canonical_positive_work_mass == pytest.approx(
        side.side_positive_work / side.recipient_work
    )


def test_maximal_hard_coarsening_keeps_same_time_physical_self_loop_but_creates_no_event_or_progress():
    triad, _ = signed_good_integer_triad()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=single_hard_role_map(triad)
    )
    assert len(quotient.donor_charges) == 1
    assert len({charge.cell for charge in quotient.recipient_charges}) == 1
    assert quotient.self_loop_atom_count == len(quotient.atoms)
    assert quotient.self_loop_mass == pytest.approx(quotient.total_positive_work_mass)
    assert quotient.self_loops_zero_recursion_depth
    assert not quotient.self_loops_create_scale_progress
    assert not quotient.donor_work_is_new_owner
    assert not quotient.negative_failure_payment_inferred
    assert not quotient.between_time_inventory_matching_inferred
    assert all(atom.coarse_self_loop for atom in quotient.atoms)
    assert all(not atom.creates_new_event and not atom.creates_scale_progress for atom in quotient.atoms)


def test_generic_two_donor_case_recombines_to_one_canonical_recipient_charge_not_two_causes():
    triad = generic_two_donor_counterexample()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=fine_hard_role_map(triad)
    )
    assert triad.donor_kernel.donor_count == 2
    assert triad.donor_kernel.recipient_count == 1
    assert len(quotient.donor_charges) == 2
    assert len(quotient.recipient_charges) == 1
    recipient = quotient.recipient_charges[0]
    assert recipient.incoming_donor_cell_count == 2
    assert quotient.overlapping_recipient_charge_count == 1
    assert recipient.incoming_donor_mass == pytest.approx(recipient.canonical_positive_work_mass)
    assert recipient.charged_once_downstream
    assert not recipient.creates_new_causal_law


def test_restricting_one_of_two_donor_cells_gives_a_true_submeasure_of_the_same_recipient_charge():
    triad = generic_two_donor_counterexample()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=fine_hard_role_map(triad)
    )
    selected = quotient.donor_charges[0].cell
    restricted = pushforward_restricted_hard_cell_donor_work(quotient, donor_cells=(selected,))
    assert restricted.recipient_total_mass == pytest.approx(restricted.selected_negative_work_mass)
    assert restricted.mass_conservation_native_residual < 5e-10
    assert restricted.dominated_by_full_recipient_charge
    assert len(restricted.recipient_subcharges) == 1
    sub = restricted.recipient_subcharges[0]
    assert 0.0 < sub.submeasure_mass < sub.full_canonical_positive_work_mass
    assert sub.incoming_selected_donor_cell_count == 1
    assert not restricted.creates_new_event
    assert not restricted.creates_new_owner
    assert not restricted.later_hahn_used
    assert not restricted.between_time_matching_used


def test_disjoint_donor_cell_partition_recombines_exactly_to_full_recipient_fate_charges():
    triad = generic_two_donor_counterexample()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=fine_hard_role_map(triad)
    )
    pieces = [
        pushforward_restricted_hard_cell_donor_work(quotient, donor_cells=(charge.cell,))
        for charge in quotient.donor_charges
    ]
    full = {(charge.cell, charge.fate): charge.canonical_positive_work_mass for charge in quotient.recipient_charges}
    recombined = {key: 0.0 for key in full}
    for piece in pieces:
        for sub in piece.recipient_subcharges:
            recombined[(sub.cell, sub.fate)] += sub.submeasure_mass
    assert set(recombined) == set(full)
    for key in full:
        assert recombined[key] == pytest.approx(full[key])
    assert sum(piece.selected_negative_work_mass for piece in pieces) == pytest.approx(
        quotient.total_negative_work_mass
    )


def test_role_map_must_bind_every_and_only_physical_mode_identity():
    triad, _ = signed_good_integer_triad()
    roles = fine_hard_role_map(triad)
    missing = dict(roles)
    missing.pop(next(iter(missing)))
    with pytest.raises(ValueError, match="every and only physical mode"):
        hard_cell_single_charge_quotient(triad, quotient_measure_mass=1.0, mode_roles=missing)

    extra = dict(roles)
    k, s, a = _generic_data()
    other = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    alien = next(mode for mode in required_hard_modes(other) if mode not in roles)
    extra[alien] = "alien"
    with pytest.raises(ValueError, match="every and only physical mode"):
        hard_cell_single_charge_quotient(triad, quotient_measure_mass=1.0, mode_roles=extra)


def test_exact_zero_native_work_is_rejected_rather_than_minting_floating_donor_provenance():
    k, s, _ = _generic_data()
    triad = register_closed_helical_triad(
        wavevectors=k,
        helicities=s,
        amplitudes=(0.0 + 0.0j, 1.0 + 0.3j, -0.2 + 0.9j),
    )
    assert triad.donor_kernel.native_work_scale == 0.0
    with pytest.raises(ValueError, match="unresolved near-zero"):
        hard_cell_single_charge_quotient(
            triad, quotient_measure_mass=1.0, mode_roles=fine_hard_role_map(triad)
        )


def test_parent_permutation_changes_storage_only_not_physical_hard_cell_single_charge_table():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    roles = fine_hard_role_map(base)
    q0 = hard_cell_single_charge_quotient(base, quotient_measure_mass=0.7, mode_roles=roles)
    perm = (2, 0, 1)
    moved = register_closed_helical_triad(
        wavevectors=tuple(k[i] for i in perm),
        helicities=tuple(s[i] for i in perm),
        amplitudes=tuple(a[i] for i in perm),
    )
    q1 = hard_cell_single_charge_quotient(moved, quotient_measure_mass=0.7, mode_roles=roles)
    sig0 = _atom_signature(q0)
    sig1 = _atom_signature(q1)
    assert [(a, b, f) for a, b, f, _ in sig0] == [(a, b, f) for a, b, f, _ in sig1]
    assert [m for *_, m in sig0] == pytest.approx([m for *_, m in sig1], rel=3e-11, abs=3e-12)


def test_uniform_amplitude_scaling_preserves_cells_and_cubic_physical_work_charge():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    roles = fine_hard_role_map(base)
    q0 = hard_cell_single_charge_quotient(base, quotient_measure_mass=1.0, mode_roles=roles)
    lam = 1.7
    scaled = register_closed_helical_triad(
        wavevectors=k,
        helicities=s,
        amplitudes=tuple(lam * z for z in a),
    )
    q1 = hard_cell_single_charge_quotient(scaled, quotient_measure_mass=1.0, mode_roles=roles)
    assert [(a, b, f) for a, b, f, _ in _atom_signature(q0)] == [
        (a, b, f) for a, b, f, _ in _atom_signature(q1)
    ]
    assert q1.total_positive_work_mass == pytest.approx(lam**3 * q0.total_positive_work_mass, rel=3e-11)
    assert q1.total_negative_work_mass == pytest.approx(lam**3 * q0.total_negative_work_mass, rel=3e-11)


def test_negative_work_is_not_converted_into_failed_good_payment_even_when_donor_and_good_recipient_coarsen_together():
    triad, _ = signed_good_integer_triad()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=single_hard_role_map(triad)
    )
    assert quotient.good_recipient_mass > 0.0
    assert quotient.total_negative_work_mass > 0.0
    assert quotient.self_loop_mass > 0.0
    assert not quotient.negative_failure_payment_inferred
    assert all(charge.charged_once_downstream for charge in quotient.recipient_charges)


def test_multiple_closed_triad_measure_pieces_landing_on_same_hard_support_aggregate_to_one_charge_per_cell_fate():
    signed_good, _ = signed_good_integer_triad()
    two_donor = generic_two_donor_counterexample()
    q1 = hard_cell_single_charge_quotient(
        signed_good, quotient_measure_mass=0.4, mode_roles=single_hard_role_map(signed_good)
    )
    q2 = hard_cell_single_charge_quotient(
        two_donor, quotient_measure_mass=0.7, mode_roles=single_hard_role_map(two_donor)
    )
    combined = aggregate_hard_cell_single_charge_quotients((q1, q2))
    keys = [(charge.cell, charge.fate) for charge in combined.recipient_charges]
    assert len(keys) == len(set(keys))
    assert combined.total_positive_work_mass == pytest.approx(
        q1.total_positive_work_mass + q2.total_positive_work_mass
    )
    assert combined.total_negative_work_mass == pytest.approx(
        q1.total_negative_work_mass + q2.total_negative_work_mass
    )
    assert combined.self_loop_atom_count == len(combined.atoms)
    assert combined.self_loop_mass == pytest.approx(combined.total_positive_work_mass)
    assert all(charge.charged_once_downstream for charge in combined.recipient_charges)


def test_semantic_guards_reject_new_owner_progress_failure_payment_and_between_time_matching():
    triad, _ = signed_good_integer_triad()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=single_hard_role_map(triad)
    )
    with pytest.raises(ValueError, match="progress, ownership, failure payment, or time matching"):
        replace(quotient, self_loops_create_scale_progress=True)
    with pytest.raises(ValueError, match="progress, ownership, failure payment, or time matching"):
        replace(quotient, donor_work_is_new_owner=True)
    with pytest.raises(ValueError, match="progress, ownership, failure payment, or time matching"):
        replace(quotient, negative_failure_payment_inferred=True)
    with pytest.raises(ValueError, match="progress, ownership, failure payment, or time matching"):
        replace(quotient, between_time_inventory_matching_inferred=True)
    atom = quotient.atoms[0]
    with pytest.raises(ValueError, match="may not create cause, event depth, or scale progress"):
        replace(atom, creates_new_event=True)
    with pytest.raises(ValueError, match="may not create cause, event depth, or scale progress"):
        replace(atom, replaces_recipient_causal_law=True)


def test_restricted_submeasure_domination_is_certified_on_native_work_scale_not_tiny_realized_recipient_mass():
    triad = generic_two_donor_counterexample()
    quotient = hard_cell_single_charge_quotient(
        triad, quotient_measure_mass=1.0, mode_roles=fine_hard_role_map(triad)
    )
    restricted = pushforward_restricted_hard_cell_donor_work(
        quotient, donor_cells=(quotient.donor_charges[0].cell,)
    )
    sub = restricted.recipient_subcharges[0]
    # Native-scale roundoff may put a floating submeasure infinitesimally above
    # its independently reconstructed full recipient mass.  This is numerical
    # certification only; the exact theorem remains submeasure domination.
    replace(
        sub,
        submeasure_mass=sub.full_canonical_positive_work_mass + 1.0e-12 * sub.native_work_mass_scale,
    )
    with pytest.raises(AssertionError, match="native work scale"):
        replace(
            sub,
            submeasure_mass=sub.full_canonical_positive_work_mass + 1.0e-8 * sub.native_work_mass_scale,
        )
