from dataclasses import replace
import math

import pytest

from src.high_tail_natural_window_reentry import comparable_natural_window_common_work_upper
from src.resolved_contact_native_binding import (
    CONTACT_PARENT_UPPER_RATIO,
    HH_WINDOW,
    K_RELAY,
    S_STRAIN,
    SignedResolvedKSAtom,
    SingleChargedRecipientMixedCause,
    coalesce_recipient_mixed_cause,
    boundary_contact_counterexample,
    canonical_contact_hh_natural_window_capacity_upper,
    canonical_hh_edge_total_variation_young_bridge,
    coarse_hahn_cancellation_counterexample,
    contact_hh_direct_natural_window_reentry,
    cover_canonical_mixed_submeasure_by_ks,
    hard_tail_resolved_contact_route,
    interior_contact_fixture,
    resolved_contact_component_route,
    resolved_contact_smooth_binding,
    theorem_certificate,
)


def _interior_atom():
    _, _, split = interior_contact_fixture()
    return next(
        atom
        for atom in split.atoms
        if atom.resolved_scale_parent_contact
        and atom.recipient_shell_index == 2
        and 0.125 < min(atom.interaction_parent_radii) / atom.recipient_shell_scale < 0.2
    )


def _single_recipient_cause(mass: float) -> SingleChargedRecipientMixedCause:
    return SingleChargedRecipientMixedCause(
        recipient_closed_mode_index=0,
        donor_closed_mode_indices=(1,),
        canonical_mixed_submeasure_mass=mass,
        common_unit_scale=4.0,
    )


def test_certificate_is_cutoff_uniform_signed_first_and_keeps_common_N_unit():
    cert = theorem_certificate()
    assert "before any new Hahn" in cert["signed_first"]
    assert "no plateau" in cert["cutoff_uniformity"]
    assert "5M/4" in cert["contact_geometry"]
    assert "N dW" in cert["causal_unit"]
    assert cert["claims_global_regularity"] is False


def test_resolved_contact_has_unique_low_parent_and_upper_comparable_other_parent():
    atom = _interior_atom()
    out = resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=0.37)
    M = out.recipient_shell_scale
    assert out.resolved_parent_radius <= M / 4 + 5e-12 * M
    assert out.uv_parent_radius > M / 4 - 5e-12 * M
    assert out.maximum_parent_to_shell_ratio <= CONTACT_PARENT_UPPER_RATIO + 5e-12
    assert out.mixed_vh_submeasure_mass + out.hh_complement_submeasure_mass == pytest.approx(
        out.canonical_positive_mass
    )
    assert out.common_unit_scale == atom.boundary


def test_same_canonical_atom_can_repartition_between_mixed_and_hh_without_changing_cause():
    atom = _interior_atom()
    transition = resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=0.31)
    fully_resolved = resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=1.0)
    assert transition.canonical_positive_mass == pytest.approx(fully_resolved.canonical_positive_mass)
    assert 0.0 < transition.mixed_vh_submeasure_mass < transition.canonical_positive_mass
    assert 0.0 < transition.hh_complement_submeasure_mass < transition.canonical_positive_mass
    assert fully_resolved.mixed_vh_submeasure_mass == pytest.approx(fully_resolved.canonical_positive_mass)
    assert fully_resolved.hh_complement_submeasure_mass == pytest.approx(0.0)
    assert transition.same_time_donor_provenance_preserved
    assert fully_resolved.same_time_donor_provenance_preserved


def test_boundary_contact_is_actual_hh_counterexample_to_contact_equals_interface():
    out = boundary_contact_counterexample()
    assert out["resolved_parent_to_shell_ratio"] == pytest.approx(0.25)
    assert out["mixed_fraction"] == pytest.approx(0.0)
    assert out["hh_fraction"] == pytest.approx(1.0)
    assert out["contact_is_not_interface_owner"] is True

    from src.hard_tail_true_upward_supply import deep_upward_resolved_contact_fixture
    _, _, split = deep_upward_resolved_contact_fixture()
    atom = min(
        (a for a in split.atoms if a.resolved_scale_parent_contact),
        key=lambda a: abs(min(a.interaction_parent_radii) / a.recipient_shell_scale - 0.25),
    )
    with pytest.raises(ValueError):
        resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=0.2)


def test_binding_rejects_nonpositive_cutoff_semantics_and_causal_mutations():
    atom = _interior_atom()
    with pytest.raises(ValueError):
        resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=-0.1)
    good = resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=0.5)
    with pytest.raises(ValueError):
        replace(good, later_hahn_used=True)
    with pytest.raises(ValueError):
        replace(good, owner_mass_cloned=True)
    with pytest.raises(ValueError):
        replace(good, recipient_shell_reweighting_used=True)


def test_multiple_donor_sidecars_coalesce_before_one_recipient_KS_charge():
    atom = _interior_atom()
    full = resolved_contact_smooth_binding(atom, resolved_parent_cutoff_value=0.4)
    other_donor = next(
        i for i in (0, 1, 2)
        if i not in (full.recipient_closed_mode_index, full.donor_closed_mode_index)
    )
    first = replace(
        full,
        canonical_positive_mass=0.5 * full.canonical_positive_mass,
        mixed_vh_submeasure_mass=0.5 * full.mixed_vh_submeasure_mass,
        hh_complement_submeasure_mass=0.5 * full.hh_complement_submeasure_mass,
    )
    second = replace(first, donor_closed_mode_index=other_donor)
    cause = coalesce_recipient_mixed_cause((first, second))
    assert set(cause.donor_closed_mode_indices) == {full.donor_closed_mode_index, other_donor}
    assert cause.canonical_mixed_submeasure_mass == pytest.approx(full.mixed_vh_submeasure_mass)
    with pytest.raises(ValueError):
        coalesce_recipient_mixed_cause((first, first))
    with pytest.raises(ValueError):
        replace(cause, recipient_submeasure_single_charged=False)


def test_signed_KS_identity_covers_canonical_mixed_cause_without_fraction_matching():
    skew = cover_canonical_mixed_submeasure_by_ks(
        _single_recipient_cause(1.5),
        SignedResolvedKSAtom(signed_mixed_work=2.0, signed_skew_work=3.0, signed_strain_work=-1.0),
    )
    assert skew.positive_cover_margin >= 0.0
    assert skew.joint_owner_witnesses == (K_RELAY,)
    assert skew.canonical_cause_unsplit is True
    assert skew.owner_fraction_matching_used is False

    strain = cover_canonical_mixed_submeasure_by_ks(
        _single_recipient_cause(1.5),
        SignedResolvedKSAtom(signed_mixed_work=2.0, signed_skew_work=-1.0, signed_strain_work=3.0),
    )
    assert strain.joint_owner_witnesses == (S_STRAIN,)


def test_KS_cover_fails_closed_before_signed_identity_or_if_cause_exceeds_mixed_positive_work():
    with pytest.raises(ValueError):
        SignedResolvedKSAtom(2.0, 0.3, 0.4)
    with pytest.raises(ValueError):
        cover_canonical_mixed_submeasure_by_ks(
            SingleChargedRecipientMixedCause(
                recipient_closed_mode_index=0,
                donor_closed_mode_indices=(1,),
                canonical_mixed_submeasure_mass=2.1,
                common_unit_scale=1.0,
            ),
            SignedResolvedKSAtom(2.0, 0.5, 1.5),
        )


def test_integrated_contact_owner_routes_only_to_hh_window_K_relay_or_existing_strain():
    out = resolved_contact_component_route(
        required_contact_common_work_lower=1.0,
        actual_contact_common_work=1.2,
        hh_complement_common_work=0.6,
        mixed_common_work=0.6,
        positive_skew_common_work=0.3,
        positive_strain_common_work=0.3,
    )
    assert set(out["joint_physical_continuations"]) == {HH_WINDOW, K_RELAY, S_STRAIN}
    assert out["hh_owner_threshold"] == pytest.approx(0.5)
    assert out["ks_owner_threshold"] == pytest.approx(0.25)
    assert out["later_hahn_used"] is False
    assert out["recipient_shell_reweighting_used"] is False


def test_hard_tail_contact_corollary_keeps_exact_clean_nuD_thresholds():
    out = hard_tail_resolved_contact_route(
        physical_tail_dissipation=2.0,
        viscosity=0.5,
        actual_contact_common_work=0.6,
        hh_complement_common_work=0.3,
        mixed_common_work=0.3,
        positive_skew_common_work=0.15,
        positive_strain_common_work=0.15,
    )
    assert out["clean_contact_owner_lower"] == pytest.approx(0.5)
    assert out["clean_contact_HH_lower"] == pytest.approx(0.25)
    assert out["clean_K_or_S_lower"] == pytest.approx(0.125)
    assert set(out["joint_physical_continuations"]) == {HH_WINDOW, K_RELAY, S_STRAIN}


def test_canonical_HH_edge_positive_law_binds_through_total_variation_before_clean_Young():
    bridge = canonical_hh_edge_total_variation_young_bridge()
    assert bridge.canonical_positive_submeasure_dominated_by_edge_variation is True
    assert bridge.edge_variation_dominated_by_capacity_measure is True
    assert bridge.aggregate_hahn_used is False
    assert bridge.edge_variation_to_clean_young_ratio < 1.0
    assert bridge.exact_edge_variation_prefactor_over_A3 == pytest.approx(
        bridge.clean_pair_young_prefactor_over_A3
        * bridge.edge_variation_to_clean_young_ratio
    )
    assert bridge.edge_variation_to_clean_young_ratio == pytest.approx(
        math.sqrt(2.0) * bridge.unitary_fourier_factor
    )


def test_positive_cutoff_complement_contraction_removes_old_factor_four_from_HH_capacity():
    args = dict(
        window_peak_child_mass=1.7,
        parent_frequency=2.0,
        global_energy=3.0,
        scaled_lifetime=0.8,
        locality_radius=5.0 / 4.0,
    )
    new = canonical_contact_hh_natural_window_capacity_upper(**args)
    old = comparable_natural_window_common_work_upper(**args)
    assert 4.0 * new == pytest.approx(old)


def test_contact_hh_goes_directly_to_sliding_natural_window_with_R_5_over_4():
    N = 1.0
    c = 1.0
    j = 1
    M = 2.0
    T = c / (M * M)
    out = contact_hh_direct_natural_window_reentry(
        positive_hh_shell_common_work={1: 2.0, 2: 1.0},
        required_total_hh_common_work_lower=1.5,
        parent_frequency=N,
        global_energy=1.0,
        scaled_lifetime=c,
        viscosity=0.1,
        maximum_window_common_work=0.5,
        window_length=T,
        window_peak_child_mass=1.0,
    )
    assert out["selected_shell_level"] == j
    assert out["parent_upper_comparability_ratio"] == pytest.approx(5.0 / 4.0)
    assert out["output_scale_locality_theorem_used"] is False
    assert out["later_hahn_used"] is False
    assert out["scale_time_tradeoff_margin"] >= -1e-12
    assert math.isfinite(out["full_survivor_own_scale_service_lower"])


def test_downstream_coarse_Hahn_cannot_replace_atomic_canonical_positive_cause():
    out = coarse_hahn_cancellation_counterexample()
    assert out["canonical_positive_first_atom"] == pytest.approx(1.0)
    assert out["coarse_positive_hahn_mass"] == pytest.approx(0.1)
    assert out["atomic_cause_exceeds_coarse_hahn"] is True
