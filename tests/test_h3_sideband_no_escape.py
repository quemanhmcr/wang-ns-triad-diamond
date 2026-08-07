import math
from src.h3_sideband_no_escape import (
    NO_ESCAPE_COEFF, actual_daughter_lower, exact_constant_certificate,
    first_impulse_lower, no_escape_quadratic_cost,
)


def test_clean_constants():
    I=.2
    assert math.isclose(first_impulse_lower(I)**2,3*I*I/32,rel_tol=1e-14)
    assert math.isclose(actual_daughter_lower(I)**2,3*I*I/128,rel_tol=1e-14)
    assert math.isclose(no_escape_quadratic_cost(I),3*I*I/4096,rel_tol=1e-14)
    assert exact_constant_certificate()['pair_rescue_or_net_deficit']=='3 I^2/4096'
