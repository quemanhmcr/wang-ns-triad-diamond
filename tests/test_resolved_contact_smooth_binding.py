from dataclasses import replace

import pytest

from src.hard_tail_true_upward_supply import deep_upward_resolved_contact_fixture
from src.resolved_contact_smooth_binding import (
    SignedResolvedKSAtom,
    bind_canonical_mixed_submeasure_to_ks,
    boundary_contact_counterexample,
    canonical_positive_resolved_cutoff,
    coarse_hahn_cancellation_counterexample,
    deep_contact_smooth_repartition,
    strict_deep_resolved_mixed_fixture,
    theorem_certificate,
)


def test_certificate_keeps_signed_first_canonical_cause_and_common_N_unit():
    cert = theorem_certificate()
    assert "before any new Hahn" in cert["signed_first"]
    assert "M>=8N" in cert["strict_deep"]
    assert "M=4N" in cert["borderline"]
    assert "before downstream aggregation" in cert["coarse_hahn_forbidden"]
    assert cert["later_hahn_used"] is False
    assert cert["cause_cloning_used"] is False
    assert cert["recipient_shell_reweighting_used"] is False
    assert cert["claims_global_regularity"] is False


def test_canonical_cutoff_is_positive_core_plateau_transition_and_zero_at_support_boundary():
    M = 8.0
    assert canonical_positive_resolved_cutoff(0.0, M) == 1.0
    assert canonical_positive_resolved_cutoff(M / 8.0, M) == 1.0
    q = canonical_positive_resolved_cutoff(3.0 * M / 16.0, M)
    assert 0.0 < q < 1.0
    assert canonical_positive_resolved_cutoff(M / 4.0, M) == 0.0
    assert canonical_positive_resolved_cutoff(M, M) == 0.0


def test_certified_boundary_contact_fixture_is_actual_HH_not_interface_under_smooth_cutoff():
    counter = boundary_contact_counterexample()
    assert counter["donor_to_shell_ratio"] == pytest.approx(0.25)
    assert counter["cutoff_value"] == 0.0
    assert counter["mixed_fraction"] == 0.0
    assert counter["hh_fraction"] == pytest.approx(1.0)
    assert counter["support_contact_is_not_interface"] is True


def test_deep_fixture_positive_cause_repartitions_without_low_low_clone_or_M_reweight():
    _, _, split = deep_upward_resolved_contact_fixture()
    deep = [a for a in split.atoms if a.deep_upward_shell]
    assert deep
    for atom in deep:
        bind = deep_contact_smooth_repartition(atom)
        assert bind.low_low_bound_mass == 0.0
        assert bind.mixed_vh_bound_mass + bind.transition_hh_bound_mass == pytest.approx(atom.physical_work_mass)
        assert bind.common_unit_scale == atom.boundary
        assert bind.later_hahn_used is False
        assert bind.coarse_hahn_push_used is False
        assert bind.cloned_owner_mass is False
        with pytest.raises(ValueError):
            replace(bind, recipient_shell_reweighting_used=True)


def test_actual_strict_deep_triad_geometry_is_fully_mixed_before_KS_handoff():
    _, _, split = strict_deep_resolved_mixed_fixture()
    strict = [a for a in split.atoms if a.recipient_shell_index >= 3 and a.donor_radius <= 1.0 + 5e-12]
    assert strict
    for atom in strict:
        bind = deep_contact_smooth_repartition(atom)
        assert atom.recipient_shell_scale >= 8.0 * atom.boundary
        assert bind.donor_cutoff_value == 1.0
        assert bind.other_parent_cutoff_value == 0.0
        assert bind.mixed_vh_bound_mass == pytest.approx(atom.physical_work_mass)
        assert bind.transition_hh_bound_mass == 0.0


def test_same_atom_signed_KS_identity_binds_canonical_submeasure_without_cloning():
    signed = SignedResolvedKSAtom(signed_mixed_work=3.0, signed_skew_work=4.5, signed_strain_work=-1.5)
    bind = bind_canonical_mixed_submeasure_to_ks(1.2, signed, common_unit_scale=7.0)
    assert bind.skew_bound_mass + bind.strain_bound_mass == pytest.approx(1.2)
    assert bind.skew_bound_mass <= bind.available_positive_skew_work
    assert bind.strain_bound_mass <= bind.available_positive_strain_work
    assert bind.strain_bound_mass == 0.0
    assert bind.common_unit_scale == 7.0
    assert bind.owner_mass_cloned is False


def test_KS_binding_rejects_wrong_signed_identity_unquotiented_gauge_and_overlarge_cause():
    with pytest.raises(ValueError):
        SignedResolvedKSAtom(1.0, 0.8, 0.3)
    with pytest.raises(ValueError):
        SignedResolvedKSAtom(1.0, 0.6, 0.4, observer_gauge_quotiented_or_fixed_event=False)
    signed = SignedResolvedKSAtom(1.0, 0.6, 0.4)
    with pytest.raises(ValueError):
        bind_canonical_mixed_submeasure_to_ks(1.1, signed, common_unit_scale=2.0)


def test_downstream_coarse_Hahn_cannot_receive_atomic_canonical_positive_cause():
    counter = coarse_hahn_cancellation_counterexample()
    assert counter["atomic_cause_exceeds_coarse_hahn"] is True
    assert counter["canonical_positive_first_atom"] == pytest.approx(1.0)
    assert counter["coarse_positive_hahn_mass"] == pytest.approx(0.1)
