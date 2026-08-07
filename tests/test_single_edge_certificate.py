import itertools
import math

import numpy as np

from src.helical import coupling_magnitude_closed
from src.single_edge_certificate import (
    A_CUSP,
    B_TANGENT,
    GLOBAL_GAP,
    U0,
    V0,
    float_envelope,
    float_jstar,
    hodge_residual_energy,
    local_mixed_rhs,
    residual_coordinates,
)


def direct_J(x: float, y: float, signs: tuple[int, int, int]) -> float:
    sx, sy, sz = signs
    area_term = (x + y + 1) * (-x + y + 1) * (x - y + 1) * (x + y - 1)
    if area_term <= 0:
        return 0.0
    g = coupling_magnitude_closed(x, y, 1.0, sx, sy, sz)
    return math.log(1.0 / y) * abs(sx * x - sy * y) * g


def test_exact_sign_reduced_envelope_matches_exhaustion():
    rng = np.random.default_rng(20260807)
    for _ in range(1000):
        y = rng.uniform(0.5001, 0.999)
        x = rng.uniform(1.0 - y + 1e-8, y)
        brute = max(direct_J(x, y, s) for s in itertools.product((-1, 1), repeat=3))
        env = float_envelope(x, y)
        assert math.isclose(brute, env, rel_tol=2e-12, abs_tol=2e-12)


def test_reference_extremizer_and_residual_coordinates():
    jstar = float_jstar()
    assert abs(jstar - 0.10011017585619) < 5e-12
    # At the symmetric extremizer both structural residuals vanish.
    from src.single_edge_certificate import float_rstar
    r = float_rstar()
    u, v = residual_coordinates(r, r)
    assert abs(u) < 1e-14
    assert abs(v) < 1e-14


def test_mixed_bound_implies_half_hodge_on_local_box_randomly():
    rng = np.random.default_rng(19)
    from src.single_edge_certificate import float_rstar
    r = float_rstar()
    gamma = math.log(1.0 / r)
    jstar = float_jstar()
    for _ in range(5000):
        u = rng.uniform(0.0, float(U0))
        v = rng.uniform(-float(V0), float(V0))
        R = r * math.exp(-v)
        x = R * math.exp(-u / 2)
        y = R * math.exp(u / 2)
        deficit = 1.0 - float_envelope(x, y) / jstar
        assert deficit + 1e-12 >= local_mixed_rhs(x, y)
        assert deficit + 1e-12 >= 0.5 * hodge_residual_energy(x, y)
        uu, vv = residual_coordinates(x, y)
        assert abs(uu - u) < 2e-14
        assert abs(vv - v) < 2e-14


def test_constants_are_the_documented_rationals():
    assert float(U0) == 0.08
    assert float(V0) == 0.08
    assert float(A_CUSP) == 0.02
    assert float(B_TANGENT) == 1.0
    assert float(GLOBAL_GAP) == 0.01
