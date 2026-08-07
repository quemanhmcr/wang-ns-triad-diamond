import math
import random

from src.master_no_escape import (
    MasterStep,
    cap_or_entropy_constants,
    cross_penalty,
    master_episode_bound,
    master_lower_bound,
    spherical_cap_mass_lower_bound,
    verify_episode_trace,
)


def test_cap_or_entropy_constants():
    alpha = 1.0
    eps = 0.1
    c = cap_or_entropy_constants(alpha, eps)
    bnorm = 1.0 - c["eta"]
    assert spherical_cap_mass_lower_bound(bnorm, alpha) >= 1.0 - eps - 1e-12
    assert c["entropy_floor"] > 0
    assert abs(c["potential_reset"] + math.log(math.cos(alpha))) < 1e-14


def test_cross_penalty_exact_algebra():
    c = 0.3
    eta = 0.01
    lhs = -math.log(math.exp(-c) + eta)
    rhs = c - cross_penalty(c, eta)
    assert abs(lhs - rhs) < 1e-14


def test_episode_count_identity():
    # 3 costly separators create at most 4 flat episodes.
    margin = master_episode_bound(
        depth=11, costly_steps=3, potential_reset=0.7,
        kappa0=0.3, potential_error=0.5
    )
    # free=8: (3+1)*.7+.5-8*.3=.9
    assert abs(margin - 0.9) < 1e-12


def test_master_rate_with_resets():
    steps = [MasterStep(math.exp(-0.1), 0.1, 0.0, False) for _ in range(10)]
    out = master_lower_bound(steps, c0=0.08, kappa0=0.2, potential_reset=0.6)
    expected_nk = (0.2 * 10 - 0.6) / 0.8
    assert abs(out["costly_step_lower_bound"] - expected_nk) < 1e-12
    assert out["effective_rate"] > 0


def test_adversarial_recharge_is_counted_by_episodes():
    # This is exactly the logical failure of the old P0-only theorem: every
    # costly separator may reset the potential. The corrected theorem counts
    # each reset as a new episode.
    low = MasterStep(1.0, 0.0, 0.0, True, 0.0)
    episodes = [[low, low], [low], [low, low], [low]]
    out = verify_episode_trace(episodes, c0=0.1, kappa0=0.2, potential_reset=0.5)
    assert out["admissible"]
    assert out["costly_steps"] == 3


def test_random_admissible_episode_trees():
    rng = random.Random(18)
    for _ in range(2000):
        kappa = rng.uniform(0.1, 0.4)
        reset = rng.uniform(kappa, 1.5)
        episodes = []
        for _e in range(rng.randint(1, 12)):
            p = reset
            ep = []
            for _ in range(30):
                z = rng.uniform(0.0, 0.02)
                if p - kappa + z < 0:
                    break
                p = p - kappa + z
                ep.append(MasterStep(1.0, 0.0, 0.0, True, z))
            episodes.append(ep)
        out = verify_episode_trace(episodes, c0=0.05, kappa0=kappa, potential_reset=reset)
        assert out["admissible"]
