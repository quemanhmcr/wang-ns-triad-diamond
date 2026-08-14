import math

import pytest

from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.curl_spectral_curvature_balance import (
    aggregate_weighted_curl_spectral_production,
    critical_majority_helicity_multiplier,
    critical_singleton_child_multiplier,
    critical_tanaka_triad_balance,
    helicity_pair_production,
    sharp_critical_frontier_geometry,
    theorem_certificate,
    weighted_curl_spectral_production,
)
from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    signed_good_integer_triad,
)
from src.helical_energy_helicity_barycentric_rigidity import (
    certify_helical_energy_helicity_rigidity,
    curl_eigenvalue,
)


def _kernel(triad, q=1.0):
    return cyclic_triad_measure_kernel(triad, quotient_measure_mass=q)


def test_all_weighted_readings_use_same_cyclic_work_law():
    triad, _ = signed_good_integer_triad()
    kernel = _kernel(triad, 1.7)
    for phi in (
        lambda x: 1.0,
        lambda x: x,
        lambda x: x*x,
        abs,
        lambda x: x*x*x,
    ):
        out = weighted_curl_spectral_production(triad, kernel, phi)
        assert out.rooted_flow_native_residual < 1.5e-9
        assert out.affine_energy_residual < 1.5e-9
        assert out.affine_helicity_residual < 1.5e-9


def test_energy_and_helicity_are_the_affine_kernel_but_enstrophy_is_curvature():
    triad, _ = signed_good_integer_triad()
    kernel = _kernel(triad)
    energy = weighted_curl_spectral_production(triad, kernel, lambda _: 1.0)
    helicity = weighted_curl_spectral_production(triad, kernel, lambda x: x)
    enstrophy = weighted_curl_spectral_production(triad, kernel, lambda x: x*x)
    rigidity = certify_helical_energy_helicity_rigidity(triad)
    assert energy.rooted_signed_production == pytest.approx(0.0, abs=2e-11)
    assert helicity.rooted_signed_production == pytest.approx(0.0, abs=2e-10)
    assert enstrophy.rooted_signed_production == pytest.approx(unitary_fourier_convolution_factor()*rigidity.quadratic_moment_production, rel=2e-10, abs=2e-11)
    assert enstrophy.rooted_signed_production > 0.0


def test_critical_Hhalf_is_exact_tanaka_defect_and_homochiral_is_zero():
    spread, _ = signed_good_integer_triad()
    out = critical_tanaka_triad_balance(spread, _kernel(spread, 2.3))
    assert out.heterochiral
    assert out.spread
    assert out.tanaka_defect > 0.0
    assert out.critical_production_mass > 0.0
    assert out.critical_production_mass == pytest.approx(out.donor_flow_critical_mass, rel=2e-10, abs=2e-11)

    homo = generic_two_donor_counterexample()
    zero = critical_tanaka_triad_balance(homo, _kernel(homo, 0.7))
    assert zero.homochiral_zero
    assert not zero.heterochiral
    assert zero.tanaka_defect == pytest.approx(0.0, abs=2e-12)
    assert zero.critical_production_mass == pytest.approx(0.0, abs=2e-11)
    assert certify_helical_energy_helicity_rigidity(homo).quadratic_moment_production != 0.0


def test_global_sum_is_signed_not_a_fake_convex_monotone_budget():
    spread, _ = signed_good_integer_triad()
    contraction = generic_two_donor_counterexample()
    rows = ((spread, _kernel(spread, 1.0)), (contraction, _kernel(contraction, 2000.0)))
    total = aggregate_weighted_curl_spectral_production(rows, lambda x: x*x)
    expected = unitary_fourier_convolution_factor() * (
        certify_helical_energy_helicity_rigidity(spread).quadratic_moment_production
        + 2000.0 * certify_helical_energy_helicity_rigidity(contraction).quadratic_moment_production
    )
    assert total == pytest.approx(expected, rel=3e-10, abs=3e-11)
    # Both spread and contraction are legal; no theorem declares this sum monotone.
    assert expected < 0.0


def test_sharp_critical_geometry_is_not_the_existing_log_progress_extremizer():
    out = sharp_critical_frontier_geometry()
    assert out.same_helicity_parent_ratio == pytest.approx(0.4539303256551502, rel=2e-13)
    assert out.opposite_helicity_parent_ratio == pytest.approx(0.8242109621975628, rel=2e-13)
    assert out.critical_multiplier == pytest.approx(0.09902212964575472, rel=2e-13)
    assert out.critical_multiplier > out.boundary_upper
    assert out.log_progress_efficiency_ratio == pytest.approx(0.27152694, rel=2e-7)
    assert out.log_progress_efficiency_ratio < 0.3
    assert out.rational_interior_witness_squared > out.rational_boundary_upper_squared


def test_two_heterochiral_frontier_families_are_distinct():
    D = 0.4539303256551502
    S = 0.8242109621975628
    majority = critical_majority_helicity_multiplier(D, S)
    singleton = critical_singleton_child_multiplier(D, S)
    assert majority > singleton
    assert majority > 0.099


def test_certificate_states_one_curvature_law_not_new_causality():
    cert = theorem_certificate()
    assert "three-point" in cert["weak_nonlinear_law"]
    assert "lambda=0" in cert["tanaka"]
    assert "0.2715" in cert["anti_J_primitive"]
    assert "canonical edge dW+" in cert["causal_scope"]
    assert cert["claims_global_regularity"] is False


def test_full_family_critical_dynamics_has_one_common_helicity_pair_source():
    spread, _ = signed_good_integer_triad()
    contraction = generic_two_donor_counterexample()
    rows = ((spread, _kernel(spread, 1.3)), (contraction, _kernel(contraction, 0.8)))
    out = helicity_pair_production(rows)
    assert out.positive_helicity_critical_production == pytest.approx(out.negative_helicity_critical_production, rel=3e-10, abs=3e-11)
    assert out.absolute_critical_production == pytest.approx(2.0*out.common_pair_source, rel=3e-10, abs=3e-11)
    assert out.signed_helicity_production == pytest.approx(0.0, abs=3e-10)
    assert out.native_residual < 1.5e-9


def test_homochiral_traffic_cannot_feed_either_critical_pair_reservoir():
    homo = generic_two_donor_counterexample()
    out = helicity_pair_production(((homo, _kernel(homo, 3.0)),))
    assert out.common_pair_source == pytest.approx(0.0, abs=3e-11)
    assert out.absolute_critical_production == pytest.approx(0.0, abs=3e-11)
