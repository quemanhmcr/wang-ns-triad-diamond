import math
import numpy as np
import pytest

from src.affine_aspect_sticky import (
    EXTENDED_ASPECT,
    MILD_ASPECT,
    arb_extended_certificate,
    condition_growth_upper,
    extended_h1_no_escape_constants,
    extended_qpol_lower,
    fresh_critical_mass_from_condition_lower,
    log_condition_rate_exact,
    predecessor_condition_lower,
    shell_radius_from_condition_lower,
)


def test_extended_threshold_exact_relation():
    assert EXTENDED_ASPECT == (27 * MILD_ASPECT / 25)


def test_extended_qpol_clean():
    assert extended_qpol_lower(float(EXTENDED_ASPECT)) > 1.0 / 4000.0


def test_condition_growth_inverse():
    k0, K = 1.3, 0.02
    k1 = condition_growth_upper(k0, K)
    assert abs(predecessor_condition_lower(k1, K) - k0) < 1e-14


def test_viscosity_does_not_increase_condition_in_isotropic_zero_strain():
    A = np.zeros((3, 3))
    Sigma = np.diag([1.0, 2.0, 4.0])
    assert log_condition_rate_exact(A, Sigma, 1.0) < 0


def test_shell_and_fresh_mass_condition_laws():
    k = 8.0
    assert abs(shell_radius_from_condition_lower(k) - 4.0 / 3.0) < 1e-14
    assert abs(fresh_critical_mass_from_condition_lower(k) - 0.4) < 1e-14


def test_extended_constants():
    c = extended_h1_no_escape_constants()
    assert c["h1_pair_or_deficit"] == 1.0 / 28_800_000.0


def test_arb_optional():
    pytest.importorskip("flint")
    c = arb_extended_certificate()
    assert c["extended_qpol_lower"] == "1/4000"
