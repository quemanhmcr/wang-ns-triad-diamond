import math

import numpy as np
import pytest

from src.high_tail_ultraviolet_locality import (
    STATUS,
    comparable_work_binary_reentry,
    hh_output_shell_law,
    high_tail_hh_locality_tradeoff,
    highpass_l32_gradient_constant,
    theorem_certificate,
    ultraviolet_common_work_upper,
    ultraviolet_hh_work_constant,
)


def _atoms(net: float) -> np.ndarray:
    return np.array([[[net + 3.0, -1.0], [-1.0, -1.0]]], dtype=float)


def test_exact_uv_geometric_constant_collapses_to_three_sqrt_pi():
    assert highpass_l32_gradient_constant() == pytest.approx((4.0 * math.pi / 3.0) ** (1.0 / 6.0))
    assert ultraviolet_hh_work_constant() == pytest.approx(3.0 * math.sqrt(math.pi), rel=2e-15)


def test_output_shell_law_is_actual_positive_HH_work_not_atomic_mass():
    out = hh_output_shell_law({1: 2.0, 2: 1.0})
    assert out["selected_shell_level"] == 1
    assert out["p_max"] == pytest.approx(2.0 / 3.0)
    assert out["total_positive_HH_common_work"] == pytest.approx(3.0)


def test_uv_bound_uses_child_mass_and_physical_Dtail_only():
    mu = 0.04
    D = 3.0
    out = ultraviolet_common_work_upper(mu, D, selected_shell_level=1, locality_radius=2.0)
    expect = 3.0 * math.sqrt(math.pi / 2.0) * math.sqrt(mu) * D
    assert out == pytest.approx(expect)


def test_continuous_R2_tradeoff_and_exact_balanced_mass_constant():
    D = 2.0
    nu = 1.0
    # HH owner lower is nu D/2=1.0. Selected p=1/2.
    work = {1: 0.5, 2: 0.5}
    p = 0.5
    clean_weighted_mass = nu * nu / (72.0 * math.pi)
    # Put weighted mass exactly at the balanced threshold.
    mu = clean_weighted_mass * p * p
    out = high_tail_hh_locality_tradeoff(D, nu, work, {1: mu, 2: 0.1}, 2.0)
    assert out["continuous_tradeoff_lhs"] >= 0.5 - 1e-13
    assert out["clean_entropy_weighted_child_mass_lower"] == pytest.approx(clean_weighted_mass)
    assert "critical_child_shell" in out["joint_clean_owners"]
    # Exact threshold also leaves a comparable quarter because sqrt(pi)<2.
    assert "comparable_parent_HH_work" in out["joint_clean_owners"]


def test_small_child_mass_forces_comparable_HH_work():
    D = 4.0
    nu = 0.8
    work = {1: 1.0, 2: 0.6}  # total=1.6=nu D/2, p=5/8
    out = high_tail_hh_locality_tradeoff(D, nu, work, {1: 1e-10, 2: 0.2}, 2.0)
    assert out["selected_shell_level"] == 1
    assert "comparable_parent_HH_work" in out["joint_clean_owners"]
    assert out["entropy_weighted_comparable_work_lower"] >= nu * D / 4.0 - 1e-11
    assert out["coherent_refinement_used_for_locality"] is False
    assert out["atomic_Hahn_mass_used_for_output_scale_law"] is False


def test_large_child_mass_routes_directly_to_generic_critical_shell():
    D = 1.0
    nu = 1.0
    work = {1: 0.3, 2: 0.2}
    out = high_tail_hh_locality_tradeoff(D, nu, work, {1: 0.5, 2: 0.5}, 2.0)
    assert "critical_child_shell" in out["joint_clean_owners"]
    assert out["next_owner_if_critical"] == "generic_critical_shell_first_stop"


def test_localized_comparable_source_is_atomized_only_after_fourier_locality():
    out = comparable_work_binary_reentry(2, 3.0, [_atoms(2.0), _atoms(1.0)])
    assert out["actual_comparable_positive_common_work"] == pytest.approx(3.0)
    assert out["binary_positive_common_work"] >= 3.0
    assert out["locality_was_read_before_coherent_refinement"] is True
    assert out["productivity_gate_supplied"] is False


def test_certificate_keeps_native_continuous_tradeoff_and_scope():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "3 sqrt(pi)" in cert["ultraviolet_bound"]
    assert "W_comp e^(Hinf_out)" in cert["native_tradeoff"]
    assert "before coherent Hahn refinement" in cert["output_scale_law"]
    assert "atomic Hahn positive mass is forbidden" in cert["representation_rule"]
    assert "no temporal natural-window concentration" in cert["scope"]
