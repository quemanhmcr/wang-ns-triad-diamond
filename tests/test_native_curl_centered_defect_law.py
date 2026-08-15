import math

import pytest

from src.native_curl_centered_defect_law import (
    curl_centered_state,
    defect_balance_rate,
    energy_dissipation_split,
    second_divided_difference,
    theorem_certificate,
    triad_curvature_current,
    triad_observable_response,
    viscous_defect_rate,
)


def test_curl_centered_defect_is_exact_signed_frequency_variance():
    a = (-2.0, -0.75, 1.25, 3.0)
    e = (0.7, 1.4, 2.1, 0.9)
    out = curl_centered_state(a, e)
    lam = sum(ai * ei for ai, ei in zip(a, e)) / sum(e)
    variance = sum(ei * (ai - lam) ** 2 for ai, ei in zip(a, e)) / sum(e)
    assert out.mean_signed_frequency == pytest.approx(lam)
    assert out.defect_variance == pytest.approx(variance)
    assert out.defect_energy == pytest.approx(out.energy * variance)
    assert out.defect_energy == pytest.approx(out.enstrophy - out.helicity**2 / out.energy)


def test_energy_dissipation_splits_into_beltrami_baseline_plus_defect():
    out = curl_centered_state((-2.0, 1.0, 3.0), (1.0, 2.0, 0.5))
    split = energy_dissipation_split(out, 0.17)
    assert split["minus_half_energy_rate"] == pytest.approx(0.17 * out.enstrophy)
    assert split["defect_dissipation"] == pytest.approx(0.17 * out.defect_energy)
    assert split["beltrami_baseline"] + split["defect_dissipation"] == pytest.approx(
        split["minus_half_energy_rate"]
    )


def test_zero_defect_is_exact_single_curl_eigenvalue_state():
    out = curl_centered_state((2.5, 2.5, 2.5), (0.5, 3.0, 1.2))
    assert out.mean_signed_frequency == pytest.approx(2.5)
    assert out.defect_energy == pytest.approx(0.0, abs=1.0e-13)
    assert out.curl_defect_energy == pytest.approx(0.0, abs=1.0e-13)
    assert viscous_defect_rate(out, 0.2) == pytest.approx(0.0, abs=1.0e-13)


def test_defect_balance_is_stretching_minus_defect_diffusion():
    assert defect_balance_rate(3.5, 4.0, 0.25) == pytest.approx(5.0)


def _closed_triad_data():
    # Exact scalar-current form T0=(a1-a2)R, etc.; no temporal ancestry is used.
    a = (0.7, -0.9, 1.2)
    R = 1.3
    T = ((a[1] - a[2]) * R, (a[2] - a[0]) * R, (a[0] - a[1]) * R)
    return a, T


def test_triad_curvature_current_annihilates_affine_observables():
    a, T = _closed_triad_data()
    triad_curvature_current(a, T)
    constant = math.fsum(T)
    affine = math.fsum((2.0 - 3.0 * ai) * ti for ai, ti in zip(a, T))
    assert constant == pytest.approx(0.0, abs=1.0e-13)
    assert affine == pytest.approx(0.0, abs=1.0e-13)


def test_second_divided_difference_represents_every_tested_observable_response():
    a, T = _closed_triad_data()
    for phi in (lambda x: x * x, lambda x: x**4, math.exp):
        out = triad_observable_response(a, T, phi)
        assert out["direct_response"] == pytest.approx(out["represented_response"])
        assert out["residual"] == pytest.approx(0.0, abs=1.0e-12)
    assert second_divided_difference(a, lambda x: x * x) == pytest.approx(1.0)


def test_quadratic_curvature_current_is_the_enstrophy_response():
    a, T = _closed_triad_data()
    Q = triad_curvature_current(a, T)
    response = math.fsum(ai * ai * ti for ai, ti in zip(a, T))
    assert Q == pytest.approx(response)


def test_certificate_keeps_the_claim_native_and_draft():
    cert = theorem_certificate()
    assert "actual viscous energy dissipation" in cert["energy"]
    assert "Beltrami" in cert["nonlinear_force"]
    assert "variance production" in cert["interpretation"]
    assert "not a master scalar" in cert["defect_role_guard"]
    assert "raw-power extremizer" in cert["transport_geometry_guard"]
    assert "phase/current geometry" in cert["open_question"]
    assert "without assuming the exact r_* value" in cert["open_question"]
    assert cert["temporal_matching_used"] is False
    assert cert["owner_bookkeeping_used"] is False
    assert cert["global_regularity_claimed"] is False
