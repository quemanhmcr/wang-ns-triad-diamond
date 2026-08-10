import math

import numpy as np
import pytest

from src.high_tail_natural_window_reentry import (
    STATUS,
    comparable_hh_temporal_shell_reentry,
    comparable_natural_window_common_work_upper,
    natural_window_geometry,
    sliding_window_piecewise_constant,
    temporal_concentration_statistics,
    theorem_certificate,
)
from src.high_tail_ultraviolet_locality import high_tail_hh_locality_tradeoff


def _comparable_locality_route(D: float = 2.0, nu: float = 1.0):
    # H=nu D/2=1; selected p=0.6; tiny child peak forces comparable owner.
    return high_tail_hh_locality_tradeoff(
        D,
        nu,
        {1: 0.6, 2: 0.4},
        {1: 1e-10, 2: 0.2},
        2.0,
    )


def test_hard_tail_geometry_has_real_scale_progress_and_time_shortening():
    g = natural_window_geometry(3.0, 1, 0.7)
    assert g["selected_shell_frequency"] == pytest.approx(6.0)
    assert g["forward_scale_ratio"] == pytest.approx(2.0)
    assert g["natural_time_ratio"] == pytest.approx(0.25)
    g3 = natural_window_geometry(3.0, 3, 0.7)
    assert g3["forward_scale_ratio"] == pytest.approx(8.0)
    assert g3["natural_time_ratio"] == pytest.approx(1.0 / 64.0)


def test_uniform_density_sliding_fraction_is_exact_natural_time_ratio():
    g = natural_window_geometry(2.0, 2, 1.0)
    L = g["parent_natural_duration"]
    T = g["selected_natural_window"]
    out = sliding_window_piecewise_constant([0.0], [L], [3.0], 0.0, L, T)
    assert out["p_time"] == pytest.approx(T / L)
    assert out["p_time"] == pytest.approx(1.0 / 16.0)


def test_sliding_measure_is_invariant_under_origin_units_and_refinement():
    a = np.array([0.0, 0.2, 0.55])
    b = np.array([0.2, 0.55, 1.0])
    d = np.array([1.0, 3.0, 0.5])
    base = sliding_window_piecewise_constant(a, b, d, 0.0, 1.0, 0.25)
    moved = sliding_window_piecewise_constant(a + 7.0, b + 7.0, d, 7.0, 8.0, 0.25)
    lam = 13.0
    scaled = sliding_window_piecewise_constant(lam * a, lam * b, d / lam, 0.0, lam, lam * 0.25)
    ar = np.array([0.0, 0.1, 0.2, 0.375, 0.55, 0.775])
    br = np.array([0.1, 0.2, 0.375, 0.55, 0.775, 1.0])
    dr = np.array([1.0, 1.0, 3.0, 3.0, 0.5, 0.5])
    refined = sliding_window_piecewise_constant(ar, br, dr, 0.0, 1.0, 0.25)
    assert moved["p_time"] == pytest.approx(base["p_time"], abs=1e-14)
    assert scaled["p_time"] == pytest.approx(base["p_time"], abs=1e-14)
    assert refined["p_time"] == pytest.approx(base["p_time"], abs=1e-14)


def test_temporal_statistics_are_measure_native_not_bin_native():
    out = temporal_concentration_statistics(4.0, 1.0, 0.25)
    assert out["p_time"] == pytest.approx(0.25)
    assert out["H_inf_time"] == pytest.approx(math.log(4.0))


def test_common_unit_window_capacity_has_required_parent_N_factor():
    mu = 0.04
    N = 5.0
    E = 2.0
    c = 0.5
    R = 2.0
    out = comparable_natural_window_common_work_upper(mu, N, E, c, R)
    expect = 24.0 * c * math.sqrt(math.pi) * N * E * math.sqrt(mu)
    assert out == pytest.approx(expect)


def test_composition_yields_actual_high_tail_shell_event_and_weighted_lower():
    D, nu = 2.0, 1.0
    route = _comparable_locality_route(D, nu)
    p_s = route["p_max"]
    W_lower = route["comparable_parent_common_work_lower"]
    W = 1.1 * W_lower
    p_t = 0.4
    Ww = p_t * W
    N, c = 2.0, 0.8
    R = route["locality_radius"]
    mu_block = route["child_peak_critical_mass"]
    # Pick fixture energy from the certified block peak so the window capacity
    # is feasible before testing the theorem composition itself.
    target_mu = 0.1 * mu_block
    E_needed = Ww / (12.0 * c * math.sqrt(math.pi) * R * N * math.sqrt(target_mu))
    E = 2.0 * E_needed
    C = 12.0 * c * math.sqrt(math.pi) * R * N * E
    mu_required = (Ww / C) ** 2
    mu = 1.1 * mu_required
    assert mu < mu_block
    T = natural_window_geometry(N, int(route["selected_shell_level"]), c)["selected_natural_window"]
    out = comparable_hh_temporal_shell_reentry(route, D, nu, N, E, c, W, Ww, T, mu)
    clean = nu * D / (48.0 * c * math.sqrt(math.pi) * R * N * E)
    assert out["entropy_weighted_sqrt_child_mass"] >= clean - 1e-14
    assert out["p_scale"] == pytest.approx(p_s)
    assert out["p_time"] == pytest.approx(p_t)
    assert out["forward_scale_ratio"] >= 2.0
    assert out["natural_time_ratio"] <= 0.25
    assert out["next_owner"] == "generic_critical_shell_first_stop"
    assert out["full_survivor_service_is_conditional"] is True
    assert out["time_partition_used"] is False
    assert out["packet_persistence_used"] is False


def test_composition_rejects_a_window_from_the_wrong_scale():
    D, nu = 2.0, 1.0
    route = _comparable_locality_route(D, nu)
    W = 1.1 * route["comparable_parent_common_work_lower"]
    Ww = 0.4 * W
    N, c = 2.0, 0.8
    R = route["locality_radius"]
    mu_block = route["child_peak_critical_mass"]
    target_mu = 0.1 * mu_block
    E_needed = Ww / (12.0 * c * math.sqrt(math.pi) * R * N * math.sqrt(target_mu))
    E = 2.0 * E_needed
    C = 12.0 * c * math.sqrt(math.pi) * R * N * E
    mu = 1.1 * (Ww / C) ** 2
    T = natural_window_geometry(N, int(route["selected_shell_level"]), c)["selected_natural_window"]
    with pytest.raises(ValueError, match="selected shell natural window"):
        comparable_hh_temporal_shell_reentry(route, D, nu, N, E, c, W, Ww, 2.0 * T, mu)


def test_certificate_keeps_cutoff_scope_and_no_time_bins():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "no time bins" in cert["sliding_measure"]
    assert "density/lambda" in cert["time_gauge"]
    assert "||(I-S)u||_2<=2||u||_2" in cert["cutoff_scope"]
    assert "12 c sqrt(pi) R N E_global" in cert["window_capacity"]
    assert "forward scale ratio>=2" in cert["hard_tail_progress"]
    assert "no packet persistence" in cert["scope"]
