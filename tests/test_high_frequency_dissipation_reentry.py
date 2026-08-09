import math

import pytest

from src.high_frequency_dissipation_reentry import (
    STATUS,
    classify_high_tail_energy_owners,
    classify_regeneration_work_owners,
    direct_high_enstrophy_shell_counterexample,
    high_tail_clean_reentry_thresholds,
    high_tail_scaled_gradient_bounds,
    inherited_branch_clean_shell_mass_lower,
    lp_high_clean_reentry_thresholds,
    physical_tail_dissipation_lower_from_lp,
    inherited_tail_shell_witness,
    integrated_high_lp_currency,
    positive_shell_work_disintegration,
    scaled_tail_energy_from_shell_masses,
    theorem_certificate,
)


def test_high_enstrophy_alone_has_no_uniform_shell_floor():
    c10 = direct_high_enstrophy_shell_counterexample(10, 1.0)
    c50 = direct_high_enstrophy_shell_counterexample(50, 1.0)
    assert c10["high_currency"] == pytest.approx(1.0)
    assert c50["high_currency"] == pytest.approx(1.0)
    assert c50["critical_shell_mass"] < 1e-10 * c10["critical_shell_mass"]


def test_hard_annulus_gradient_bridge_constants():
    U = [0.7, 0.2, 0.05]
    D = integrated_high_lp_currency(U)
    lo, hi = high_tail_scaled_gradient_bounds(U)
    assert lo == pytest.approx(D / 4.0)
    assert hi == pytest.approx(D)


def test_energy_gate_keeps_exact_tie_joint():
    D_tail = 0.8
    nu = 1.3
    threshold = nu * D_tail
    out = classify_high_tail_energy_owners(D_tail, nu, threshold, threshold)
    assert set(out["joint_owners"]) == {
        "inherited_tail_energy",
        "positive_nonlinear_regeneration",
    }
    assert out["resolved_DV_relabel"] is False


def test_inherited_tail_energy_exposes_real_high_shell():
    masses = [0.3, 0.8, 0.4, 0.1]
    tail = scaled_tail_energy_from_shell_masses(masses)
    out = inherited_tail_shell_witness(masses)
    assert out["critical_shell_mass"] >= tail
    assert out["shell_to_block_frequency_ratio"] >= 2.0


def test_clean_physical_tail_threshold_is_nu_Dtail():
    D_tail = 0.6
    nu = 0.9
    assert inherited_branch_clean_shell_mass_lower(D_tail, nu) == pytest.approx(nu * D_tail)
    th = high_tail_clean_reentry_thresholds(D_tail, nu)
    assert th["inherited_critical_shell_mass"] == pytest.approx(nu * D_tail)
    assert th["HH_or_interface_work_if_regeneration"] == pytest.approx(nu * D_tail)
    assert th["master_semantics"] == "RECURSE_CRITICAL"


def test_lp_supplier_constant_is_explicit_and_hard_annulus_gives_quarter():
    D_lp = 0.6
    c_lp = 0.25
    nu = 0.9
    D_tail = physical_tail_dissipation_lower_from_lp(D_lp, c_lp)
    assert D_tail == pytest.approx(D_lp / 4.0)
    th = lp_high_clean_reentry_thresholds(D_lp, c_lp, nu)
    assert th["physical_tail_dissipation_lower"] == pytest.approx(D_lp / 4.0)
    assert th["inherited_critical_shell_mass"] == pytest.approx(nu * D_lp / 4.0)


def test_positive_tail_work_becomes_own_scale_shell_work():
    W = 0.25
    # N W_j^+ sums exactly to the positive tail-work mass.
    shell = [0.1, 0.08, 0.07]
    out = positive_shell_work_disintegration(W, shell)
    assert out["scaled_shell_positive_work_sum"] == pytest.approx(W)
    assert out["own_scale_positive_shell_work"] >= 2.0 * W


def test_low_low_free_regeneration_split_keeps_joint_tie():
    S = 0.8
    out = classify_regeneration_work_owners(S, 0.4, 0.4)
    assert set(out["joint_owners"]) == {
        "positive_HH_regeneration",
        "positive_resolved_cross_interface",
    }
    assert out["HH_is_productivity_generated_branch"] is False
    assert out["interface_is_free"] is False


def test_certificate_forbids_two_false_shortcuts():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "neither smooth-LP d_high nor resolved low-pass D_V" in cert["native_currency"]
    assert "D_tail>=c_LP D_high" in cert["lp_supplier"]
    assert "no critical-shell floor" in cert["anti_relabel"]
    assert "not automatically" in cert["no_false_productivity"]
