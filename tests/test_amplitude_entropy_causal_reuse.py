import math

import numpy as np

from src.amplitude_entropy_causal_reuse import (
    anchor_coefficient_energy_fraction,
    closed_form_reuse_lower,
    entropy_energy_amplitude_upper,
    l2_normalized_dual_l32_scaled_upper,
    one_layer_log_product_lower,
    registered_coefficient_productivity_lower,
    reuse_information_lower_without_mass_floor,
    root_expected_log_lower,
)
from src.weighted_causal_reuse import entropy


def test_dual_probe_l32_scaling_is_scale_free_after_sqrt_N():
    c = l2_normalized_dual_l32_scaled_upper(0.0, 2.0 / 3.0)
    expected = (4.0 / 3.0) * math.pi ** 0.25 * math.sqrt(3.0 / 2.0)
    assert math.isclose(c, expected)


def test_registered_productivity_is_positive_without_any_root_mass_floor():
    assert registered_coefficient_productivity_lower(1.0, 1.0) > 0.0


def test_log_product_recursion_handles_extreme_parent_imbalance():
    lam = 0.01
    child = 3.0
    a = math.sqrt(lam * child) * 1e12
    b = math.sqrt(lam * child) * 1e-12
    slot_log = 0.5 * (math.log(a) + math.log(b))
    assert math.isclose(slot_log, one_layer_log_product_lower(math.log(child), lam), rel_tol=1e-12, abs_tol=1e-12)


def test_depth_solution_is_exact_for_equal_product_splits():
    lam = 0.03
    terminal = 7.0
    ell = math.log(terminal)
    for _ in range(9):
        ell = one_layer_log_product_lower(ell, lam)
    assert math.isclose(ell, root_expected_log_lower(9, lam, terminal), rel_tol=1e-13, abs_tol=1e-13)


def test_entropy_energy_amplitude_logsum_has_correct_direction():
    p = np.array([0.65, 0.25, 0.10])
    a = np.array([0.2, 4.0, 13.0])
    lhs = entropy(p) + 2.0 * float(np.dot(p, np.log(a)))
    rhs = math.log(float(np.dot(a, a)))
    assert lhs <= rhs + 1e-14
    assert math.isclose(entropy_energy_amplitude_upper(a, p), rhs - 2.0 * float(np.dot(p, np.log(a))))


def test_closed_form_matches_entropy_route_before_zero_truncation():
    L = 8
    N0 = 2.0
    E = 10.0
    P = 3.0
    lam = 0.02
    terminal = 0.7
    beta = anchor_coefficient_energy_fraction()
    out = reuse_information_lower_without_mass_floor(L, N0, E, P, lam, terminal, beta)
    closed = closed_form_reuse_lower(L, N0, E, P, lam, terminal, beta)
    if out["feasible"] and float(out["root_entropy_upper"]) > 0:
        assert math.isclose(float(out["shannon_reuse_lower"]), closed, rel_tol=1e-13, abs_tol=1e-13)
