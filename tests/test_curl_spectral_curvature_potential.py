import math

import pytest

from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    generic_two_donor_counterexample,
    signed_good_integer_triad,
)
from src.curl_spectral_curvature_potential import (
    CurlSpectralCurvatureFamily,
    critical_tanaka_scale_cocycle,
    theorem_certificate,
    triad_curl_spectral_potential,
)
from src.helical_energy_helicity_barycentric_rigidity import certify_helical_energy_helicity_rigidity


def _registered_potential(triad, qmass=1.0):
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=qmass)
    return triad_curl_spectral_potential(triad, kernel)


def test_signed_good_spread_is_exact_second_derivative_tent_with_affine_kernel():
    triad, _ = signed_good_integer_triad()
    pot = _registered_potential(triad, 1.7)
    assert pot.signed_tent_mass > 0.0
    assert pot.source_reconstruction_native_residual < 1.5e-9
    assert pot.weighted_nonlinear_production(lambda _: 1.0) == pytest.approx(0.0, abs=2e-11)
    assert pot.weighted_nonlinear_production(lambda x: x) == pytest.approx(0.0, abs=2e-11)


def test_tanaka_value_enstrophy_area_and_radial_layer_cake_are_same_tent_readings():
    triad, _ = signed_good_integer_triad()
    pot = _registered_potential(triad, 0.75)
    pcrit = pot.weighted_nonlinear_production(abs)
    penst = pot.weighted_nonlinear_production(lambda x: x*x)
    assert pot.critical_production == pytest.approx(pcrit, rel=3e-10, abs=2e-11)
    assert pot.enstrophy_production == pytest.approx(penst, rel=3e-10, abs=2e-11)
    assert pot.integrated_radial_first_moment_production == pytest.approx(pcrit, rel=3e-10, abs=2e-11)
    assert pcrit > 0.0
    assert penst > 0.0


def test_contraction_is_the_same_signed_potential_with_reversed_curvature():
    triad = generic_two_donor_counterexample()
    pot = _registered_potential(triad, 1.0)
    assert pot.signed_tent_mass < 0.0
    assert pot.critical_production == pytest.approx(pot.weighted_nonlinear_production(abs), rel=3e-10, abs=2e-11)
    assert pot.enstrophy_production == pytest.approx(pot.weighted_nonlinear_production(lambda x: x*x), rel=3e-10, abs=2e-11)
    with pytest.raises(ValueError):
        critical_tanaka_scale_cocycle(pot)


def test_critical_tanaka_scale_cocycle_uses_actual_donor_work_and_no_J():
    triad, _ = signed_good_integer_triad()
    pot = _registered_potential(triad, 2.0)
    out = critical_tanaka_scale_cocycle(pot)
    assert out.normalized_critical_production > 0.0
    assert out.normalized_critical_production <= out.radial_gap_fraction + 2e-12
    assert out.radial_gap_fraction <= out.log_scale_displacement + 2e-12
    assert out.log_parabolic_lifetime_expansion == pytest.approx(2.0*out.log_scale_displacement, rel=2e-13)
    assert not out.uses_log_progress_J
    assert not out.uses_capacity_as_causality
    assert not out.creates_budget


def test_family_superposes_signed_potentials_without_claiming_convex_monotonicity():
    spread, _ = signed_good_integer_triad()
    contraction = generic_two_donor_counterexample()
    p1 = _registered_potential(spread, 0.01)
    p2 = _registered_potential(contraction, 300.0)
    fam = CurlSpectralCurvatureFamily((p1, p2))
    expected = p1.weighted_nonlinear_production(lambda x: x*x) + p2.weighted_nonlinear_production(lambda x: x*x)
    assert fam.enstrophy_production == pytest.approx(expected, rel=3e-10, abs=2e-11)
    assert fam.enstrophy_production < 0.0


def test_physical_normalization_is_exactly_the_existing_CF_closed_triad_measure():
    triad, _ = signed_good_integer_triad()
    qmass = 3.25
    pot = _registered_potential(triad, qmass)
    assert pot.physical_fourier_factor == pytest.approx(unitary_fourier_convolution_factor()*qmass, rel=2e-14)
    rigidity = certify_helical_energy_helicity_rigidity(triad)
    expected = unitary_fourier_convolution_factor()*qmass*rigidity.quadratic_moment_production
    assert pot.enstrophy_production == pytest.approx(expected, rel=3e-10, abs=2e-11)


def test_certificate_keeps_constitutive_speed_as_the_real_frontier():
    cert = theorem_certificate()
    assert "partial_t rho" in cert["full_weak_form"]
    assert "eta_crit" in cert["cocycle"]
    assert "Waleffe" in cert["constitutive_frontier"]
    assert cert["claims_global_regularity"] is False
