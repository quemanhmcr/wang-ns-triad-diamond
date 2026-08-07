import math
from src.h1_swirl_no_escape import (
    arb_physical_transport_certificate,
    clean_first_impulse_square_lower,
    exact_constant_certificate,
    full_mild_aspect_quadratic_cost,
    h1_physical_quadratic_cost,
)


def test_physical_h1_constants():
    I = 0.2
    assert math.isclose(clean_first_impulse_square_lower(I), I * I / 480, rel_tol=1e-14)
    assert math.isclose(h1_physical_quadratic_cost(I), I * I / 184320, rel_tol=1e-14)
    assert math.isclose(full_mild_aspect_quadratic_cost(I), I * I / 737280, rel_tol=1e-14)
    assert exact_constant_certificate()["physical_H1_pair_or_deficit"] == "I_1^2/184320"


def test_transport_certificate_optional():
    import pytest
    pytest.importorskip("flint")
    assert arb_physical_transport_certificate()["status"] == "CERTIFIED_PHYSICAL_LOW_STRAIN_CONDITIONING"
