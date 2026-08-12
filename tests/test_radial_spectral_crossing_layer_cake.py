from __future__ import annotations

from dataclasses import replace
import math

import numpy as np
import pytest

from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    register_closed_helical_triad,
    signed_good_integer_triad,
)
from src.helical_mode_set_energy_continuity import flow_atoms_from_cyclic_kernel
from src.radial_spectral_crossing_layer_cake import (
    STATUS,
    clipped_log_radius_potential,
    equiradial_physical_transfer_triad,
    finite_radial_log_action,
    mode_radius,
    radial_exterior_balance,
    theorem_certificate,
    truncated_radial_layer_cake,
)


def _atoms(triad):
    return flow_atoms_from_cyclic_kernel(cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0))


def _sign_reversed(triad):
    return register_closed_helical_triad(
        wavevectors=tuple(np.asarray(m.wavevector, float) for m in triad.modes),
        helicities=tuple(m.helicity for m in triad.modes),
        amplitudes=tuple(-a for a in triad.amplitudes),
    )


def test_certificate_keeps_radial_crossing_as_existing_physical_flow_not_new_hahn_or_budget():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert cert["later_hahn_used"] is False
    assert cert["gross_crossing_budget"] is False
    assert cert["hard_shell_reweighting_used"] is False
    assert "not the one-edge" in cert["young_progress_distinction"]
    assert cert["claims_global_regularity"] is False


def test_signed_good_actual_triad_crosses_upward_at_radius_eight():
    triad, _ = signed_good_integer_triad()
    balance = radial_exterior_balance(_atoms(triad), radius=8.0)
    assert balance.upward_crossing_flow > 0.0
    assert balance.downward_crossing_flow == pytest.approx(0.0, abs=1e-13)
    assert balance.tail_signed_work == pytest.approx(balance.upward_crossing_flow)
    assert balance.tail_divergence_native_residual < 1e-12
    assert balance.minimum_scale_progress_claimed is False


def test_global_sign_reversal_turns_same_physical_triad_into_downward_crossing():
    triad, _ = signed_good_integer_triad()
    reversed_triad = _sign_reversed(triad)
    balance = radial_exterior_balance(_atoms(reversed_triad), radius=8.0)
    assert balance.downward_crossing_flow > 0.0
    assert balance.upward_crossing_flow == pytest.approx(0.0, abs=1e-13)
    assert balance.tail_signed_work == pytest.approx(-balance.downward_crossing_flow)


def test_truncated_layer_cake_is_exact_clipped_log_radius_marginal_identity():
    triad = generic_two_donor_counterexample()
    atoms = _atoms(triad)
    radii = [mode_radius(m) for a in atoms for m in (a.donor_mode, a.recipient_mode)]
    lo = 0.7 * min(radii)
    hi = 1.3 * max(radii)
    layer = truncated_radial_layer_cake(atoms, lower_radius=lo, upper_radius=hi)
    expected = sum(
        a.physical_work_mass
        * (
            clipped_log_radius_potential(mode_radius(a.recipient_mode), lo, hi)
            - clipped_log_radius_potential(mode_radius(a.donor_mode), lo, hi)
        )
        for a in atoms
    )
    assert layer.signed_log_action == pytest.approx(expected, rel=1e-12, abs=1e-14)
    assert layer.signed_log_action == pytest.approx(layer.marginal_log_difference, rel=1e-12, abs=1e-14)
    assert layer.signed_marginal_identity_native_residual < 1e-12


def test_full_finite_layer_cake_equals_donor_recipient_log_radius_displacement():
    triad, _ = signed_good_integer_triad()
    atoms = _atoms(triad)
    full = finite_radial_log_action(atoms)
    expected = sum(
        a.physical_work_mass * math.log(mode_radius(a.recipient_mode) / mode_radius(a.donor_mode))
        for a in atoms
    )
    assert full.signed_log_action == pytest.approx(expected, rel=1e-12, abs=1e-14)
    assert full.signed_marginal_identity_native_residual < 1e-12
    assert full.continuum_extension_requires_log_moment is True


def test_equiradial_actual_helical_transfer_has_positive_work_but_exactly_zero_radial_progress():
    triad = equiradial_physical_transfer_triad()
    atoms = _atoms(triad)
    assert sum(a.physical_work_mass for a in atoms) > 0.0
    radii = [mode_radius(m) for a in atoms for m in (a.donor_mode, a.recipient_mode)]
    assert max(radii) - min(radii) < 5e-13
    full = finite_radial_log_action(atoms)
    assert full.upward_log_action == pytest.approx(0.0, abs=1e-13)
    assert full.downward_log_action == pytest.approx(0.0, abs=1e-13)


def test_global_reality_partner_has_same_radial_crossing_and_log_action():
    triad, _ = signed_good_integer_triad()
    reality = register_closed_helical_triad(
        wavevectors=tuple(-np.asarray(m.wavevector, float) for m in triad.modes),
        helicities=tuple(m.helicity for m in triad.modes),
        amplitudes=tuple(np.conjugate(a) for a in triad.amplitudes),
    )
    base_atoms = _atoms(triad)
    reality_atoms = _atoms(reality)
    base = radial_exterior_balance(base_atoms, radius=8.0)
    partner = radial_exterior_balance(reality_atoms, radius=8.0)
    assert partner.upward_crossing_flow == pytest.approx(base.upward_crossing_flow, rel=1e-11, abs=1e-13)
    assert partner.downward_crossing_flow == pytest.approx(base.downward_crossing_flow, rel=1e-11, abs=1e-13)
    bfull = finite_radial_log_action(base_atoms)
    pfull = finite_radial_log_action(reality_atoms)
    assert pfull.upward_log_action == pytest.approx(bfull.upward_log_action, rel=1e-11, abs=1e-13)
    assert pfull.downward_log_action == pytest.approx(bfull.downward_log_action, rel=1e-11, abs=1e-13)


def test_radial_donor_displacement_is_not_recipient_edge_progress_generically():
    triad = generic_two_donor_counterexample()
    atoms = _atoms(triad)
    mismatches = 0
    for atom in atoms:
        slot = next(s for s in triad.slots if s.edge_identity.child == atom.recipient_mode)
        radial = math.log(mode_radius(atom.recipient_mode) / mode_radius(atom.donor_mode))
        if abs(radial - slot.edge_registration.scale_progress) > 1e-8:
            mismatches += 1
    assert mismatches > 0


def test_radial_balance_refuses_semantic_promotions():
    triad, _ = signed_good_integer_triad()
    balance = radial_exterior_balance(_atoms(triad), radius=8.0)
    with pytest.raises(ValueError):
        replace(balance, later_hahn_used=True)
    with pytest.raises(ValueError):
        replace(balance, crossing_creates_event_depth=True)
    with pytest.raises(ValueError):
        replace(balance, minimum_scale_progress_claimed=True)
    with pytest.raises(ValueError):
        replace(balance, gross_crossing_declared_finite_resource=True)
    with pytest.raises(ValueError):
        replace(balance, hard_shell_reweighting_used=True)


def test_layer_cake_refuses_young_identification_and_later_hahn():
    triad, _ = signed_good_integer_triad()
    layer = truncated_radial_layer_cake(_atoms(triad), lower_radius=5.0, upper_radius=12.0)
    with pytest.raises(ValueError):
        replace(layer, identified_with_single_edge_young_progress=True)
    with pytest.raises(ValueError):
        replace(layer, later_hahn_used=True)
    with pytest.raises(ValueError):
        replace(layer, gross_log_action_declared_finite_resource=True)


def test_boundary_outside_all_modes_has_zero_crossing_without_minting_progress():
    triad, _ = signed_good_integer_triad()
    atoms = _atoms(triad)
    largest = max(mode_radius(m) for a in atoms for m in (a.donor_mode, a.recipient_mode))
    balance = radial_exterior_balance(atoms, radius=2.0 * largest)
    assert balance.upward_crossing_flow == 0.0
    assert balance.downward_crossing_flow == 0.0
    assert balance.tail_signed_work == 0.0
    assert balance.minimum_scale_progress_claimed is False
