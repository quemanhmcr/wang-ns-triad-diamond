import math

from src.smooth_symbol_freezing import (
    freezing_certificate,
    infinite_quadratic_tail_upper,
    sharp_young_constant_3d,
    symbol_freezing_error,
)


def test_sharp_young_constant_matches_repository_value():
    assert math.isclose(sharp_young_constant_3d(), (math.sqrt(3) / 2) ** 3)


def test_symbol_freezing_error_is_linear_in_cell_size():
    a = symbol_freezing_error(2.0, 0.01)
    b = symbol_freezing_error(2.0, 0.02)
    assert math.isclose(b, 2 * a)


def test_quadratic_schedule_is_summable_by_integral_test():
    assert math.isclose(infinite_quadratic_tail_upper(3), 1 / 9 + 1 / 3)
    cert = freezing_certificate(1.0, 3)
    assert cert.infinite_normalized_transfer_error_upper < 0.3
