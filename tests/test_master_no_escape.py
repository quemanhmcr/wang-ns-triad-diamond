import math
import random

from src.master_no_escape import (
    MasterStep,
    balanced_entropy_lower_bound,
    cross_penalty,
    master_lower_bound,
    spherical_cap_mass_lower_bound,
    verify_trace,
)


def test_spherical_concentration_three_quarters():
    alpha = 0.2
    eta = (1.0 - math.cos(alpha)) / 4.0
    bnorm = 1.0 - eta
    assert spherical_cap_mass_lower_bound(bnorm, alpha) >= 0.75 - 1e-12


def test_balanced_entropy_positive():
    assert balanced_entropy_lower_bound(0.1) > 0


def test_cross_penalty_exact_algebra():
    c = 0.3
    eta = 0.01
    lhs = -math.log(math.exp(-c) + eta)
    rhs = c - cross_penalty(c, eta)
    assert abs(lhs - rhs) < 1e-14


def test_master_deterministic_trace():
    c0 = 0.08
    kappa = 0.2
    p0 = 0.45
    # At most two exact free steps are admissible.
    steps = [
        MasterStep(math.exp(-0.09), 0.09, 0.0, False),
        MasterStep(1.0, 0.0, 0.0, True),
        MasterStep(math.exp(-0.12), 0.12, 0.0, False),
        MasterStep(1.0, 0.0, 0.0, True),
        MasterStep(math.exp(-0.10), 0.10, 0.0, False),
    ]
    out = verify_trace(steps, c0=c0, kappa0=kappa, potential0=p0)
    assert out["trace_admissible"]
    assert out["master_margin"] >= -1e-12


def test_random_admissible_traces():
    rng = random.Random(8)
    for _ in range(1000):
        c0 = 0.05
        kappa = 0.17
        p0 = rng.uniform(0.0, 1.5)
        p = p0
        steps = []
        for _j in range(rng.randint(1, 30)):
            make_free = p >= kappa and rng.random() < 0.35
            if make_free:
                zeta = rng.uniform(0.0, 0.02)
                # Only accept it when the untruncated potential remains nonnegative.
                if p - kappa + zeta >= -1e-12:
                    p = p - kappa + zeta
                    steps.append(MasterStep(1.0, 0.0, 0.0, True, zeta))
                    continue
            cost = rng.uniform(c0, 0.4)
            eta = rng.uniform(0.0, 0.003)
            eff = min(1.0, math.exp(-cost) + eta)
            steps.append(MasterStep(eff, cost, eta, False))
        out = verify_trace(steps, c0=c0, kappa0=kappa, potential0=p0)
        assert out["trace_admissible"]
        assert out["master_margin"] >= -2e-12


def test_master_formula_only_uses_depth_and_budgets():
    steps = [MasterStep(0.9, 0.11, 0.0, False) for _ in range(10)]
    out = master_lower_bound(steps, c0=0.1, kappa0=0.2, potential0=0.4)
    assert abs(out["free_block_bound"] - 2.0) < 1e-12
    assert abs(out["log_efficiency_lower_bound"] - 0.8) < 1e-12
