from dataclasses import replace

import pytest

from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    signed_good_integer_triad,
)
from src.helical_mode_set_energy_continuity import (
    STATUS,
    amplitude_scaled_closed_triad,
    flow_atoms_from_cyclic_kernel,
    interval_continuity_certificate,
    mode_set_boundary_balance,
    theorem_certificate,
)


def _all_modes(atoms):
    return frozenset({a.donor_mode for a in atoms} | {a.recipient_mode for a in atoms})


def test_certificate_keeps_same_time_flow_and_between_time_stock_as_native_distinct_ledgers():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "P_A=I_A+In_A" in cert["same_time_flow"]
    assert "E_A(t1)+D_A" in cert["between_time_continuity"]
    assert "gross canonical edge Hahn work" in cert["singleton"]
    assert "inward=outward=0" in cert["full_set"]
    assert "neither dissipation nor event depth" in cert["internal_flow"]
    assert "not a finite gross-transfer budget" in cert["gross_transfer_anti_theorem"]
    assert "no FIFO/LIFO" in cert["temporal_provenance"]
    assert "physical helical modes" in cert["state_ontology"]
    assert not cert["capacity_is_causal_law"]
    assert not cert["later_hahn_used"]
    assert not cert["claims_global_regularity"]


def test_full_closed_triad_mode_set_has_only_internal_physical_flow_and_zero_boundary_divergence():
    triad, _ = signed_good_integer_triad()
    atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    balance = mode_set_boundary_balance(atoms, _all_modes(atoms))
    assert balance.internal_flow > 0.0
    assert balance.inward_boundary_flow == 0.0
    assert balance.outward_boundary_flow == 0.0
    assert balance.recipient_positive_work == pytest.approx(balance.internal_flow)
    assert balance.donor_negative_work == pytest.approx(balance.internal_flow)
    assert balance.signed_nonlinear_work == pytest.approx(0.0, abs=5e-12)
    assert not balance.internal_flow_is_dissipation
    assert not balance.internal_flow_creates_event_depth
    assert not balance.internal_flow_creates_scale_progress


def test_one_mode_recipient_set_reads_positive_work_as_inward_boundary_and_no_internal_flow():
    triad, _ = signed_good_integer_triad()
    atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    recipient = atoms[0].recipient_mode
    # Pick any actual recipient node.  In a regular closed triad donor/recipient
    # roots have distinct wavevectors, so a singleton has no internal cyclic flow.
    balance = mode_set_boundary_balance(atoms, (recipient,))
    assert balance.internal_flow == 0.0
    assert balance.recipient_positive_work == pytest.approx(balance.inward_boundary_flow)
    assert balance.donor_negative_work == pytest.approx(balance.outward_boundary_flow)
    assert balance.signed_nonlinear_work == pytest.approx(
        balance.inward_boundary_flow - balance.outward_boundary_flow
    )


def test_generic_two_donor_recipient_is_one_mode_node_with_two_physical_inflows():
    triad = generic_two_donor_counterexample()
    atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    recipients = {a.recipient_mode for a in atoms}
    assert len(recipients) == 1
    recipient = next(iter(recipients))
    incoming_atoms = [a for a in atoms if a.recipient_mode == recipient]
    assert len({a.donor_mode for a in incoming_atoms}) == 2
    balance = mode_set_boundary_balance(atoms, (recipient,))
    assert balance.inward_boundary_flow == pytest.approx(sum(a.physical_work_mass for a in incoming_atoms))
    assert balance.internal_flow == 0.0


def test_physical_closed_triad_rescaling_falsifies_any_finite_gross_transfer_budget_from_set_stock():
    triad, _ = signed_good_integer_triad()
    base_atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    scaled_triad = amplitude_scaled_closed_triad(triad, 10.0)
    scaled_atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(scaled_triad, quotient_measure_mass=1.0))
    modes = _all_modes(base_atoms)
    base = mode_set_boundary_balance(base_atoms, modes)
    scaled = mode_set_boundary_balance(scaled_atoms, modes)
    assert base.inward_boundary_flow == base.outward_boundary_flow == 0.0
    assert scaled.inward_boundary_flow == scaled.outward_boundary_flow == 0.0
    assert scaled.internal_flow == pytest.approx(1000.0 * base.internal_flow)
    assert scaled.signed_nonlinear_work == pytest.approx(0.0, abs=5e-9)


def test_interval_continuity_is_stock_plus_boundary_flow_plus_viscosity_not_deposit_matching():
    triad, _ = signed_good_integer_triad()
    atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    mode = atoms[0].recipient_mode
    cert = interval_continuity_certificate(
        modes=(mode,),
        initial_energy=2.0,
        final_energy=2.5,
        integrated_inward_flow=1.0,
        integrated_outward_flow=0.25,
        viscous_dissipation=0.25,
        native_energy_throughput_scale=3.0,
    )
    assert cert.initial_plus_inward == pytest.approx(3.0)
    assert cert.final_plus_outward_plus_viscosity == pytest.approx(3.0)
    assert cert.balance_native_residual == 0.0
    assert not cert.fifo_matching_used
    assert not cert.lifo_matching_used
    assert not cert.gross_transfer_declared_finite_resource
    assert not cert.hard_interaction_cell_used_as_persistent_inventory


def test_interval_semantic_guards_reject_fifo_lifo_gross_budget_and_hard_cell_wallet():
    triad, _ = signed_good_integer_triad()
    atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    mode = atoms[0].recipient_mode
    cert = interval_continuity_certificate(
        modes=(mode,), initial_energy=1.0, final_energy=1.0,
        integrated_inward_flow=0.2, integrated_outward_flow=0.1,
        viscous_dissipation=0.1, native_energy_throughput_scale=1.2,
    )
    for field in (
        "fifo_matching_used",
        "lifo_matching_used",
        "gross_transfer_declared_finite_resource",
        "hard_interaction_cell_used_as_persistent_inventory",
    ):
        with pytest.raises(ValueError, match="may not invent"):
            replace(cert, **{field: True})


def test_internal_flow_semantic_guards_reject_dissipation_event_depth_and_scale_progress():
    triad, _ = signed_good_integer_triad()
    atoms = flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))
    balance = mode_set_boundary_balance(atoms, _all_modes(atoms))
    for field in (
        "internal_flow_is_dissipation",
        "internal_flow_creates_event_depth",
        "internal_flow_creates_scale_progress",
    ):
        with pytest.raises(ValueError, match="may not become"):
            replace(balance, **{field: True})
