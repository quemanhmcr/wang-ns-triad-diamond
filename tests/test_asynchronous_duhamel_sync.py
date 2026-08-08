import math
from fractions import Fraction

import pytest

from src.asynchronous_duhamel_sync import (
    BACKWARD_FRACTION,
    COMMON_SLICE_MARGIN,
    INITIAL_HALF_SPAN,
    LIFETIME_GROWTH_MIN,
    SYNC_CONE,
    SYNC_FIXED_POINT,
    choose_heavy_half,
    common_reference_slice,
    initial_parent_span_ratio,
    initial_root_count_upper,
    minimum_backward_displacement,
    next_span_ratio,
    registration_xi_upper,
)


def test_exact_parabolic_constants():
    assert initial_parent_span_ratio() == Fraction(25, 128)
    assert INITIAL_HALF_SPAN < SYNC_CONE
    assert COMMON_SLICE_MARGIN == Fraction(9, 40)
    assert SYNC_FIXED_POINT == Fraction(10, 39)
    assert Fraction(25, 64) * (SYNC_CONE + BACKWARD_FRACTION) == Fraction(155, 512)
    assert Fraction(155, 512) < SYNC_CONE


def test_sync_fixed_point():
    assert math.isclose(next_span_ratio(float(SYNC_FIXED_POINT)), float(SYNC_FIXED_POINT), rel_tol=0, abs_tol=1e-15)


def test_common_reference_slice_is_inside_every_natural_window():
    Tmin = 7.0
    a = 4.0
    b = a + float(SYNC_CONE) * Tmin
    s = common_reference_slice(a, b, Tmin)
    assert b - s < Tmin
    assert math.isclose(Tmin - (b - s), float(COMMON_SLICE_MARGIN) * Tmin, rel_tol=0, abs_tol=1e-14)


def test_backward_displacement_closed_form():
    T0 = 0.3
    L = 6
    brute = sum(float(BACKWARD_FRACTION) * T0 * float(LIFETIME_GROWTH_MIN) ** j for j in range(L))
    assert math.isclose(minimum_backward_displacement(T0, L), brute, rel_tol=0, abs_tol=1e-13)


def test_heavy_half_keeps_half_mass():
    row = choose_heavy_half([0.1, 0.3, 0.8], [1.0, 2.0, 4.0], 0.0, 1.0)
    assert row["mass"] >= row["total"] / 2


def test_initial_boundary_sobolev_root_count():
    # m=1 gives the expected inverse-frequency root count.
    assert math.isclose(initial_root_count_upper(10.0, 3.0, 0.2, 1.0), 1.5)


def test_registration_charge_contains_each_interface_currency_once():
    xi = registration_xi_upper(0.2, 3.0, 0.1)
    assert xi >= 0.2
    assert xi >= 0.1 * 3.0 / math.sqrt(2.0)


def test_common_slice_rejects_overwide_layer():
    with pytest.raises(ValueError):
        common_reference_slice(0.0, 0.5, 1.0)
