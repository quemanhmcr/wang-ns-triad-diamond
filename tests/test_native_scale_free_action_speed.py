import math

import pytest

from src.native_scale_free_action_speed import (
    HIGH_STRAIN_THRESHOLD,
    action_speed_lock,
    high_strain_global_rate_upper,
    maximum_action_vertices,
    objective_global_rate_upper,
    objective_low_coeffs,
    theorem_certificate,
)


def _data():
    return dict(
        global_energy=4.0,
        viscosity=0.2,
        scaled_lifetime=0.03,
        filter_kernel_l1=1.4,
        stress_reproducing_kernel_l32=0.8,
        smooth_objective_rate_coefficient=7.0,
        grad_u_linf=5.0,
    )


def test_high_strain_speed_is_scale_free():
    assert high_strain_global_rate_upper(1.4, 5.0) == pytest.approx(7.0)


def test_objective_low_high_bounds_lock_to_finite_scale_free_speed():
    data = _data()
    a, b = objective_low_coeffs(
        data["global_energy"],
        data["viscosity"],
        data["scaled_lifetime"],
        data["filter_kernel_l1"],
        data["stress_reproducing_kernel_l32"],
    )
    d = data["scaled_lifetime"] * data["smooth_objective_rate_coefficient"]
    expected = (2.0 * a) ** (2.0 / 5.0) * d ** (3.0 / 5.0)
    expected += (2.0 * b) ** (4.0 / 9.0) * d ** (5.0 / 9.0)
    got = objective_global_rate_upper(
        data["global_energy"],
        data["viscosity"],
        data["scaled_lifetime"],
        data["filter_kernel_l1"],
        data["stress_reproducing_kernel_l32"],
        data["smooth_objective_rate_coefficient"],
    )
    assert got == pytest.approx(expected)
    assert math.isfinite(got) and got > 0.0


def test_native_action_faces_have_positive_common_time_floor_and_finite_count():
    data = _data()
    out = action_speed_lock(0.05, **data)
    assert out.high_strain_time_floor == pytest.approx(HIGH_STRAIN_THRESHOLD / 7.0)
    assert out.objective_time_floor > 0.0
    assert out.common_time_floor == min(out.high_strain_time_floor, out.objective_time_floor)
    count = maximum_action_vertices(2.0, 0.05, **data)
    assert count >= 1
    assert count < math.inf


def test_certificate_states_native_no_zeno_not_scale_partition():
    cert = theorem_certificate()
    assert "uniformly in N" in cert["high_strain"]
    assert "sup_N" in cert["scale_lock"]
    assert "not an analyst scale partition" in cert["semantics"]
    assert cert["global_regularity_claimed"] is False
