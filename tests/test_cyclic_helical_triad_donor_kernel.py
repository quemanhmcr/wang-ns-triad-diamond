import itertools
import math

import numpy as np
import pytest

from src.cyclic_helical_triad_donor_kernel import (
    CHILD_TO_DONOR_HI,
    CHILD_TO_DONOR_LO,
    CLOSED_TRIAD_ROOT_CHART_DENSITY,
    EDGE_ROOT_CHART_DENSITY,
    SAME_HELICITY_MULTIPLIER_UPPER,
    SIDE_TO_CHILD_HI,
    SIDE_TO_CHILD_LO,
    SIDE_TO_DONOR_HI,
    SIDE_TO_DONOR_LO,
    STATUS,
    closed_triad_radon_certificate,
    cyclic_sum_relative_reroot,
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    global_reality_negation,
    pushforward_restricted_negative_work,
    register_closed_helical_triad,
    signed_good_integer_triad,
    theorem_certificate,
    translate_closed_amplitudes,
)
from src.helical import coupling_g


def _generic_data():
    k0 = np.array([-1.0, -2.0, 0.5])
    k1 = np.array([0.2, 1.3, -0.1])
    k2 = -(k0 + k1)
    s = (1, -1, 1)
    a = (0.7 + 0.4j, -0.3 + 1.1j, 0.9 - 0.2j)
    return (k0, k1, k2), s, a


def _work_map(triad):
    return {
        (slot.edge_identity.child.wavevector, slot.edge_identity.child.helicity): slot.signed_work
        for slot in triad.slots
    }


def test_full_s3_closed_triad_quotient_recovers_existing_unordered_edge_density():
    cert = closed_triad_radon_certificate()
    assert cert["ordered_closed_triad_S3_quotient"] == "1/6"
    assert cert["root_marks"] == 3
    assert cert["rooted_parent_swap_factor"] == "1/2"
    assert cert["closed_triad_root_chart_density"] == "1/48"
    assert cert["canonical_edge_density_after_three_roots"] == "1/16"
    assert 3 * CLOSED_TRIAD_ROOT_CHART_DENSITY == EDGE_ROOT_CHART_DENSITY
    assert cert["cyclic_reroot_abs_jacobian"] == "1"
    assert not cert["reality_negation_quotiented"]


def test_cyclic_sum_relative_reroot_is_the_same_closed_triad_with_unit_jacobian_formula():
    x = np.array([0.7, 0.2, -0.1])
    y = np.array([0.4, -0.3, 0.5])
    z = x + y
    r = x - y
    zp, rp = cyclic_sum_relative_reroot(z, r)
    assert np.allclose(zp, -x)
    assert np.allclose(rp, y + z)
    # Recover the new parents from (z',r'): y and -z, up to the unordered swap.
    p = 0.5 * (zp + rp)
    q = 0.5 * (zp - rp)
    assert any(
        np.allclose(p, a) and np.allclose(q, b)
        for a, b in ((y, -z), (-z, y))
    )


def test_three_cyclic_physical_edges_share_one_waleffe_factor_and_conserve_signed_energy_before_hahn():
    k, s, a = _generic_data()
    triad = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    assert triad.signed_energy_conservation_native_residual < 5e-10
    assert triad.cyclic_coupling_native_residual < 5e-12
    assert sum(slot.signed_work for slot in triad.slots) == pytest.approx(0.0, abs=2e-11)
    for slot in triad.slots:
        assert slot.signed_work == pytest.approx(slot.expected_signed_work, rel=3e-10, abs=3e-12)
        assert slot.edge_identity.child.wavevector == tuple(-x for x in slot.closed_mode.wavevector)
        assert slot.edge_identity.child.helicity == slot.closed_mode.helicity


def test_full_parent_permutation_is_storage_gauge_not_physical_work():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    expected = _work_map(base)
    for perm in itertools.permutations(range(3)):
        out = register_closed_helical_triad(
            wavevectors=tuple(k[i] for i in perm),
            helicities=tuple(s[i] for i in perm),
            amplitudes=tuple(a[i] for i in perm),
        )
        assert _work_map(out) == pytest.approx(expected)
        assert out.modes == base.modes
        assert out.amplitudes == pytest.approx(base.amplitudes)


def test_spatial_translation_changes_modal_phases_but_not_closed_triad_work_or_kernel():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    moved_a = translate_closed_amplitudes(k, a, (0.37, -1.1, 0.6))
    moved = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=moved_a)
    assert _work_map(moved) == pytest.approx(_work_map(base), rel=3e-11, abs=3e-12)
    assert moved.donor_kernel.total_positive_work == pytest.approx(base.donor_kernel.total_positive_work)


def test_global_reality_negation_is_covariant_but_not_quotiented():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    nk, na = global_reality_negation(k, a)
    neg = register_closed_helical_triad(wavevectors=nk, helicities=s, amplitudes=na)
    assert sorted(slot.signed_work for slot in neg.slots) == pytest.approx(
        sorted(slot.signed_work for slot in base.slots), rel=4e-11, abs=4e-12
    )
    assert not base.reality_negation_quotiented
    assert not neg.reality_negation_quotiented


def test_exact_zero_work_has_no_fake_tie_flow():
    k, s, _a = _generic_data()
    # One exact zero modal amplitude kills the physical cubic work on all three
    # roots.  Do not use an approximately pi/2 phase to manufacture a float tie.
    out = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=(0.0, 1.0, 1.0))
    assert out.common_phase_work_factor == 0.0
    assert all(slot.signed_work == 0.0 for slot in out.slots)
    assert out.donor_kernel.total_positive_work == 0.0
    assert out.donor_kernel.flows == ()
    assert out.donor_kernel.donor_count == 0
    assert out.donor_kernel.recipient_count == 0


def test_generic_positive_recipient_can_have_two_energy_donors():
    out = generic_two_donor_counterexample()
    assert out.donor_kernel.donor_count == 2
    assert out.donor_kernel.recipient_count == 1
    assert len(out.donor_kernel.flows) == 2
    assert out.donor_kernel.transport_unique
    assert not out.donor_kernel.creates_new_event
    assert not out.donor_kernel.canonical_positive_law_replaced


def test_measure_kernel_binds_theorem_canonical_marginals_on_native_physical_scale():
    k, s, a = _generic_data()
    triad = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    out = cyclic_triad_measure_kernel(triad, quotient_measure_mass=2.75)
    assert out.numerically_resolved_transport
    assert out.balance_native_residual < 5e-10
    assert out.donor_marginal_native_residual < 5e-10
    assert out.recipient_marginal_native_residual < 5e-10
    assert out.theorem_recipient_is_submeasure_of_canonical_dW_plus
    assert out.theorem_donor_is_restriction_of_canonical_dW_minus
    assert not out.canonical_dW_plus_replaced
    assert not out.creates_new_event


def test_any_nonempty_negative_root_restriction_pushes_to_a_submeasure_of_same_triad_canonical_dW_plus():
    two = generic_two_donor_counterexample()
    kernel = cyclic_triad_measure_kernel(two, quotient_measure_mass=1.3)
    donors = tuple(i for i, mass in enumerate(kernel.donor_edge_negative_masses) if mass > 0.0)
    assert len(donors) == 2
    for selected in ((donors[0],), (donors[1],), donors):
        out = pushforward_restricted_negative_work(kernel, donor_closed_mode_indices=selected)
        assert out.mass_conservation_native_residual < 5e-10
        assert all(out.recipient_dominated_by_full_canonical_positive_mass)
        for got, full in zip(out.recipient_masses, kernel.recipient_edge_positive_masses):
            assert got <= full + 5e-10 * kernel.native_work_mass_scale
        assert not out.creates_new_event
        assert not out.capacity_used_as_causal_law
        assert not out.later_hahn_used


def test_near_zero_phase_cancellation_does_not_mint_floating_donor_provenance():
    k, s, _a = _generic_data()
    g = coupling_g(k[1], k[2], k[0], s[1], s[2], s[0])
    a0 = 1j * g / abs(g)
    triad = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=(a0, 1.0, 1.0))
    assert not triad.donor_kernel.numerically_resolved_transport
    assert triad.donor_kernel.flows == ()
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    assert not kernel.numerically_resolved_transport
    assert kernel.atoms == ()
    donors = tuple(i for i,mass in enumerate(kernel.donor_edge_negative_masses) if mass>0.0)
    if donors:
        with pytest.raises(ValueError, match="numerically resolved"):
            pushforward_restricted_negative_work(kernel, donor_closed_mode_indices=(donors[0],))


def test_signed_good_forward_recipient_has_unique_energy_donor_and_positive_nonforward_side_recipient():
    triad, out = signed_good_integer_triad()
    recipient = triad.slot_for_closed_mode_index(out.recipient_closed_mode_index)
    donor = triad.slot_for_closed_mode_index(out.energy_donor_closed_mode_index)
    side = triad.slot_for_closed_mode_index(out.side_recipient_closed_mode_index)
    assert recipient.signed_efficiency > 1.0 - 1.0e-4
    assert out.parent_helicities_opposite
    assert out.unique_energy_donor
    assert out.interaction_parents_remain_two
    assert donor.signed_work < 0.0
    assert side.signed_work > 0.0
    assert side.edge_registration.scale_progress == 0.0
    assert side.edge_registration.geometric_multiplier_J == 0.0
    assert out.side_is_positive_nonforward
    assert out.side_terminal_transfer_loss_is_existing_router_consequence
    assert not out.creates_new_event
    assert float(SIDE_TO_CHILD_LO) < out.side_to_recipient_ratio < float(SIDE_TO_CHILD_HI)
    assert float(CHILD_TO_DONOR_LO) < out.recipient_to_donor_ratio < float(CHILD_TO_DONOR_HI)
    assert float(SIDE_TO_DONOR_LO) < out.side_to_donor_ratio < float(SIDE_TO_DONOR_HI)
    assert out.donor_negative_work == pytest.approx(out.recipient_work + out.side_positive_work, rel=4e-11)


def test_signed_good_same_helicity_parent_shortcut_is_impossible_by_exact_sign_factor_gap():
    assert SAME_HELICITY_MULTIPLIER_UPPER == pytest.approx(9 / 1600)
    assert float(SAME_HELICITY_MULTIPLIER_UPPER) < 1.0 - 1.0e-4
    _triad, out = signed_good_integer_triad()
    assert out.parent_helicities_opposite


def test_positive_rescaling_preserves_kernel_fractions_and_cubic_work_homogeneity():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    lam = 7.0
    scaled = register_closed_helical_triad(
        wavevectors=k, helicities=s, amplitudes=tuple(lam * z for z in a)
    )
    assert scaled.donor_kernel.total_positive_work == pytest.approx(
        lam**3 * base.donor_kernel.total_positive_work, rel=5e-11
    )
    assert scaled.donor_kernel.total_negative_work == pytest.approx(
        lam**3 * base.donor_kernel.total_negative_work, rel=5e-11
    )
    for i in range(3):
        b = base.slot_for_closed_mode_index(i)
        q = scaled.slot_for_closed_mode_index(i)
        assert q.signed_work == pytest.approx(lam**3 * b.signed_work, rel=5e-11, abs=5e-12)


def test_uniform_wavevector_dilation_preserves_kernel_fractions_and_scales_work_linearly():
    k, s, a = _generic_data()
    base = register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=a)
    lam = 5.25
    scaled = register_closed_helical_triad(
        wavevectors=tuple(lam*x for x in k), helicities=s, amplitudes=a
    )
    assert scaled.donor_kernel.total_positive_work == pytest.approx(
        lam * base.donor_kernel.total_positive_work, rel=6e-10
    )
    assert scaled.donor_kernel.total_negative_work == pytest.approx(
        lam * base.donor_kernel.total_negative_work, rel=6e-10
    )


def test_closed_triad_registration_rejects_nonclosure_duplicate_modes_and_nonfinite_amplitudes():
    k, s, a = _generic_data()
    bad = list(k)
    bad[2] = bad[2] + np.array([0.1, 0.0, 0.0])
    with pytest.raises(ValueError, match="sum to zero"):
        register_closed_helical_triad(wavevectors=bad, helicities=s, amplitudes=a)
    with pytest.raises(ValueError, match="distinct wavevectors"):
        register_closed_helical_triad(
            wavevectors=(np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0]), np.array([-2.0, 0.0, 0.0])),
            helicities=(1, -1, 1),
            amplitudes=(1.0, 1.0, 1.0),
        )
    with pytest.raises(ValueError, match="finite amplitude"):
        register_closed_helical_triad(wavevectors=k, helicities=s, amplitudes=(complex(float("nan"), 0.0), a[1], a[2]))


def test_restricted_negative_pushforward_rejects_zero_or_non_donor_selection():
    triad, _side = signed_good_integer_triad()
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    recipient_indices = tuple(i for i,mass in enumerate(kernel.recipient_edge_positive_masses) if mass>0.0)
    assert recipient_indices
    with pytest.raises(ValueError, match="no canonical negative work"):
        pushforward_restricted_negative_work(kernel, donor_closed_mode_indices=(recipient_indices[0],))
    with pytest.raises(ValueError, match="nonempty donor restriction"):
        pushforward_restricted_negative_work(kernel, donor_closed_mode_indices=())


def test_theorem_certificate_keeps_physics_and_scope_typed():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "unique positive transport" in cert["donor_kernel"]
    assert "submeasure" in cert["negative_restriction"]
    assert "two energy donors" in cert["generic_anti_theorem"]
    assert "interaction parents remain two" in cert["interaction_ontology"]
    assert "not a general coherent POVM kernel" in cert["coherent_scope"]
    assert not cert["reality_scope"].startswith("quotient")
    assert not cert["claims_global_regularity"]
