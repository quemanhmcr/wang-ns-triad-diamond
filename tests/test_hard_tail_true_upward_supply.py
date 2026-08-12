from dataclasses import replace

import pytest

from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel, signed_good_integer_triad
from src.hard_tail_true_upward_supply import (
    deep_upward_resolved_contact_fixture,
    hard_tail_upward_supply_split,
    tail_stock_upward_supply_certificate,
    theorem_certificate,
    upward_owner_support_alternative,
)
from src.helical_mode_set_energy_continuity import flow_atoms_from_cyclic_kernel
from src.radial_spectral_crossing_layer_cake import equiradial_physical_transfer_triad, radial_exterior_balance


def test_certificate_keeps_true_upward_supply_distinct_from_circulation_and_interface_overclaim():
    cert = theorem_certificate()
    assert "N E_>(t1)" in cert["tail_identity"]
    assert "low-to-high" in cert["upward_measure"]
    assert "M=2N" in cert["pure_uv_rigidity"]
    assert "not called an interface owner" in cert["interface_scope"]
    assert cert["later_hahn_used"] is False
    assert cert["fifo_lifo_used"] is False
    assert cert["claims_global_regularity"] is False


def test_signed_good_actual_upward_supply_is_pure_uv_first_shell_and_automatically_comparable():
    triad, _ = signed_good_integer_triad()
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    split = hard_tail_upward_supply_split(triad, kernel, boundary=8.0)
    pure = [a for a in split.atoms if a.pure_uv_hh_by_support]
    assert pure
    assert all(a.first_dyadic_shell for a in pure)
    assert all(a.recipient_shell_scale == 16.0 for a in pure)
    assert all(a.donor_is_interaction_parent for a in pure)
    assert all(0.25 < a.donor_radius / a.recipient_shell_scale <= 0.5 for a in pure)
    assert all(a.comparable_parent_upper_ratio <= 1.5 + 5e-12 for a in pure)
    assert split.pure_uv_nonfirst_shell_atoms == 0
    assert split.deep_pure_uv_atoms == 0
    assert split.radial_upward_binding_native_residual <= 1e-15


def test_deep_actual_upward_supply_has_resolved_scale_parent_contact_but_is_not_promoted_to_interface_owner():
    _, _, split = deep_upward_resolved_contact_fixture()
    deep = [a for a in split.atoms if a.deep_upward_shell]
    assert deep
    assert all(a.recipient_shell_index >= 2 for a in deep)
    assert all(a.resolved_scale_parent_contact for a in deep)
    assert all(not a.pure_uv_hh_by_support for a in deep)
    assert all(a.donor_radius <= 0.25 * a.recipient_shell_scale + 5e-12*a.recipient_shell_scale for a in deep)
    assert split.resolved_contact_is_interface_owner is False
    with pytest.raises(ValueError):
        replace(split, resolved_contact_is_interface_owner=True)


def test_exact_tail_identity_selects_inherited_or_true_upward_owner_in_common_N_unit():
    # N E0 + N Phi_up = 1.55 + 0.55 = 2.10, while 2 nu D=2.
    inherited = tail_stock_upward_supply_certificate(
        boundary=2.0,
        viscosity=0.5,
        initial_tail_energy=1.55/2.0,
        final_tail_energy=0.10/2.0,
        integrated_upward_work=0.55/2.0,
        integrated_downward_work=0.0,
        normalized_tail_dissipation=2.0,
    )
    assert inherited.inherited_owner
    assert not inherited.true_upward_owner
    assert inherited.owner_threshold == pytest.approx(1.0)
    assert inherited.continuity_native_residual <= 1e-15

    upward = tail_stock_upward_supply_certificate(
        boundary=2.0,
        viscosity=0.5,
        initial_tail_energy=0.55/2.0,
        final_tail_energy=0.10/2.0,
        integrated_upward_work=1.55/2.0,
        integrated_downward_work=0.0,
        normalized_tail_dissipation=2.0,
    )
    assert upward.true_upward_owner
    assert not upward.inherited_owner
    assert upward.upward_common_work == pytest.approx(1.55)
    assert upward.continuity_native_residual <= 1e-15

    joint = tail_stock_upward_supply_certificate(
        boundary=2.0,
        viscosity=0.5,
        initial_tail_energy=1.1/2.0,
        final_tail_energy=0.2/2.0,
        integrated_upward_work=1.1/2.0,
        integrated_downward_work=0.0,
        normalized_tail_dissipation=2.0,
    )
    assert joint.inherited_owner and joint.true_upward_owner


def test_true_upward_owner_splits_only_into_pure_uv_support_or_resolved_contact_without_new_causal_units():
    triad, _ = signed_good_integer_triad()
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    split = hard_tail_upward_supply_split(triad, kernel, boundary=8.0)
    alt = upward_owner_support_alternative(split, owner_threshold=0.9*split.upward_common_unit_work)
    assert alt.pure_uv_owner or alt.resolved_contact_owner
    assert alt.support_partition_native_residual <= 1e-15
    assert split.upward_common_unit_work == pytest.approx(8.0*split.upward_physical_work)
    with pytest.raises(ValueError):
        replace(split, later_hahn_used=True)
    with pytest.raises(ValueError):
        replace(split, recipient_shell_reweighting_used=True)


def test_internal_high_tail_circulation_is_real_work_but_exactly_not_upward_supply():
    triad = equiradial_physical_transfer_triad()
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    radii = [sum(v*v for v in mode.wavevector)**0.5 for mode in triad.modes]
    balance = radial_exterior_balance(flow_atoms_from_cyclic_kernel(kernel), radius=0.5*min(radii))
    assert balance.high_internal_flow > 0.0
    assert balance.upward_crossing_flow == 0.0
    assert balance.downward_crossing_flow == 0.0
    with pytest.raises(ValueError):
        # No true upward supply exists at a boundary lying below every mode.
        hard_tail_upward_supply_split(triad, kernel, boundary=0.5*min(radii))


def test_tail_supply_certificate_refuses_gross_positive_work_temporal_matching_and_own_shell_causal_units():
    cert = tail_stock_upward_supply_certificate(
        boundary=2.0,
        viscosity=0.5,
        initial_tail_energy=0.55/2.0,
        final_tail_energy=0.10/2.0,
        integrated_upward_work=1.55/2.0,
        integrated_downward_work=0.0,
        normalized_tail_dissipation=2.0,
    )
    with pytest.raises(ValueError):
        replace(cert, positive_tail_work_used_instead_of_upward_crossing=True)
    with pytest.raises(ValueError):
        replace(cert, internal_high_high_counted_as_supply=True)
    with pytest.raises(ValueError):
        replace(cert, fifo_matching_used=True)
    with pytest.raises(ValueError):
        replace(cert, own_shell_causal_reweighting_used=True)
