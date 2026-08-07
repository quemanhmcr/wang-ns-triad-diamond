import math
from src.h1_swirl_no_escape import (
    exact_constant_certificate,full_mild_aspect_quadratic_cost,
    h1_first_impulse_lower,h1_quadratic_cost,
)


def test_h1_constants():
    I=.2
    assert math.isclose(h1_first_impulse_lower(I)**2,I*I/200,rel_tol=1e-14)
    assert math.isclose(h1_quadratic_cost(I),I*I/25600,rel_tol=1e-14)
    assert math.isclose(full_mild_aspect_quadratic_cost(I),I*I/102400,rel_tol=1e-14)
    assert exact_constant_certificate()['H1_pair_or_deficit']=='I_1^2/25600'
