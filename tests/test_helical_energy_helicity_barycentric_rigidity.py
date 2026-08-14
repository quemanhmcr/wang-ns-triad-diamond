import math

import numpy as np
import pytest

from src.cyclic_helical_triad_donor_kernel import (
    generic_two_donor_counterexample,
    register_closed_helical_triad,
    signed_good_integer_triad,
)
from src.helical import coupling_g
from src.helical_energy_helicity_barycentric_rigidity import (
    certify_helical_energy_helicity_rigidity,
    theorem_certificate,
)


def test_signed_good_is_barycentric_spread_with_positive_enstrophy_and_uv_frontier():
    triad, _ = signed_good_integer_triad()
    out = certify_helical_energy_helicity_rigidity(triad)
    assert out.transfer_orientation == "mean_preserving_spread"
    assert out.median_is_unique_singleton_side
    assert out.helicity_conservation_native_residual < 8e-10
    assert out.quadratic_moment_production > 0.0
    assert out.strict_uv_frontier_positive_slots
    assert out.strict_uv_frontier_slots_are_spreads


def test_generic_two_donor_example_is_the_reverse_contraction_and_reduces_enstrophy():
    out = certify_helical_energy_helicity_rigidity(generic_two_donor_counterexample())
    assert out.transfer_orientation == "mean_preserving_contraction"
    assert out.median_is_unique_singleton_side
    assert out.quadratic_moment_production < 0.0
    assert not out.strict_uv_frontier_positive_slots


def test_phase_reversal_turns_same_homochiral_geometry_into_spread_with_frontier_recipient():
    base = generic_two_donor_counterexample()
    amps = list(base.amplitudes)
    amps[0] = -amps[0]
    triad = register_closed_helical_triad(
        wavevectors=tuple(m.wavevector for m in base.modes),
        helicities=tuple(m.helicity for m in base.modes),
        amplitudes=tuple(amps),
    )
    out = certify_helical_energy_helicity_rigidity(triad)
    assert out.transfer_orientation == "mean_preserving_spread"
    assert out.quadratic_moment_production > 0.0
    assert out.strict_uv_frontier_positive_slots


def test_equiradial_opposite_helicity_exchange_has_zero_radial_convex_moments():
    r3 = math.sqrt(3.0) / 2.0
    k = (
        np.array([1.0, 0.0, 0.0]),
        np.array([-0.5, r3, 0.0]),
        np.array([-0.5, -r3, 0.0]),
    )
    s = (1, 1, -1)
    g = coupling_g(k[1], k[2], k[0], s[1], s[2], s[0])
    assert abs(g) > 0.0
    a0 = g / abs(g)
    triad = register_closed_helical_triad(
        wavevectors=k,
        helicities=s,
        amplitudes=(a0, 1.0, 1.0),
    )
    works = tuple(slot.signed_work for slot in triad.slots)
    assert max(abs(w) for w in works) > 0.0
    for p in (1.0, 2.0, 4.0):
        radial_moment_rate = sum(
            (np.linalg.norm(slot.closed_mode.wavevector) ** p) * slot.signed_work
            for slot in triad.slots
        )
        assert radial_moment_rate == pytest.approx(0.0, abs=2e-11)


def test_certificate_keeps_new_law_as_structure_of_existing_work_not_new_causality():
    cert = theorem_certificate()
    assert "lambda=s|k|" in cert["curl_coordinate"]
    assert "energy and helicity" in cert["two_invariants"]
    assert "does not replace canonical dW+" in cert["causal_scope"]
    assert not cert["claims_global_regularity"]

from src.helical_energy_helicity_barycentric_rigidity import critical_helicity_pair_balance


def test_critical_pair_balance_on_signed_good_spread():
    triad, _ = signed_good_integer_triad()
    out = critical_helicity_pair_balance(triad)
    assert not out.homochiral
    assert out.absolute_critical_rate > 0.0
    assert out.positive_helicity_critical_rate == pytest.approx(out.negative_helicity_critical_rate, rel=2e-10, abs=2e-12)
    assert out.absolute_critical_rate == pytest.approx(2.0*out.singleton_helicity_weighted_work, rel=2e-10, abs=2e-12)
    assert out.critical_growth_has_nonforward_side
    assert out.side_weighted_positive_work + 2e-10 >= 0.5*out.absolute_critical_rate


def test_critical_pair_balance_on_homochiral_contraction_is_zero():
    triad = generic_two_donor_counterexample()
    out = critical_helicity_pair_balance(triad)
    assert out.homochiral
    assert out.absolute_critical_rate == pytest.approx(0.0, abs=2e-11)
    assert out.singleton_helicity_closed_mode_index is None

from src.helical_energy_helicity_barycentric_rigidity import good_edge_critical_source_bridge


def test_good_edge_is_uniformly_equivalent_to_native_critical_source():
    triad, side = signed_good_integer_triad()
    out = good_edge_critical_source_bridge(
        triad, recipient_closed_mode_index=side.recipient_closed_mode_index
    )
    assert 18.0/49.0 < out.critical_source_to_child_scale_work_ratio < 20.0/49.0
    assert out.critical_source_rate > 0.0
    assert not out.canonical_causality_reweighted
    assert not out.critical_source_declared_causal_probability
