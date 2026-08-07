from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class MasterStep:
    efficiency: float
    cost: float
    cross_error: float
    low_cost_flat: bool
    potential_error: float = 0.0



def barycentric_reset_constants(beta: float) -> dict[str, float]:
    if not (0.0 < beta < 1.0):
        raise ValueError("beta must lie in (0,1)")
    return {
        "entropy_floor": math.log(2.0 / (1.0 + beta)),
        "potential_reset": -math.log(beta),
    }

def cap_or_entropy_constants(alpha: float, epsilon_cap: float) -> dict[str, float]:
    if not (0.0 < alpha < math.pi / 2.0):
        raise ValueError("alpha must lie in (0,pi/2)")
    if not (0.0 < epsilon_cap < 1.0):
        raise ValueError("epsilon_cap must lie in (0,1)")
    eta = epsilon_cap * (1.0 - math.cos(alpha))
    return {
        "eta": eta,
        "entropy_floor": math.log(2.0 / (2.0 - eta)),
        "potential_reset": -math.log(math.cos(alpha)),
    }


def spherical_cap_mass_lower_bound(bnorm: float, alpha: float) -> float:
    if not (0.0 <= bnorm <= 1.0):
        raise ValueError("bnorm must be in [0,1]")
    if not (0.0 < alpha < math.pi):
        raise ValueError("alpha must lie in (0,pi)")
    eta = 1.0 - bnorm
    denom = 1.0 - math.cos(alpha)
    return max(0.0, 1.0 - eta / denom)


def cross_penalty(cost: float, eta: float) -> float:
    if cost < 0 or eta < 0:
        raise ValueError("nonnegative cost and error required")
    return math.log1p(eta * math.exp(cost))


def master_episode_bound(
    *, depth: int, costly_steps: int, potential_reset: float,
    kappa0: float, potential_error: float
) -> float:
    if depth < 0 or costly_steps < 0 or costly_steps > depth:
        raise ValueError("invalid counts")
    if potential_reset < 0 or kappa0 <= 0 or potential_error < 0:
        raise ValueError("invalid parameters")
    free = depth - costly_steps
    return (costly_steps + 1) * potential_reset + potential_error - free * kappa0


def master_lower_bound(
    steps: Iterable[MasterStep], *, c0: float, kappa0: float,
    potential_reset: float
) -> dict[str, float]:
    steps = list(steps)
    if c0 <= 0 or kappa0 <= 0 or potential_reset < 0:
        raise ValueError("invalid master constants")
    costly = sum(not s.low_cost_flat for s in steps)
    zeta = sum(s.potential_error for s in steps if s.low_cost_flat)
    xi = sum(cross_penalty(s.cost, s.cross_error) for s in steps)
    depth = len(steps)
    nk_lower = max(0.0, (kappa0 * depth - potential_reset - zeta) / (kappa0 + potential_reset))
    cost_lower = c0 * nk_lower
    lower = cost_lower - xi
    return {
        "depth": float(depth),
        "costly_steps": float(costly),
        "potential_error": zeta,
        "cross_penalty": xi,
        "costly_step_lower_bound": nk_lower,
        "cost_lower_bound": cost_lower,
        "log_efficiency_lower_bound": lower,
        "efficiency_upper_bound": math.exp(-lower),
        "effective_rate": c0 * kappa0 / (kappa0 + potential_reset),
    }


def verify_episode_trace(
    episodes: list[list[MasterStep]], *, c0: float, kappa0: float,
    potential_reset: float
) -> dict[str, float | bool]:
    # Each episode contains only low-cost flat steps. Costly separator blocks are
    # supplied separately as one block between consecutive episodes.
    zeta = 0.0
    free = 0
    admissible = True
    for ep in episodes:
        p = potential_reset
        for s in ep:
            if not s.low_cost_flat:
                raise ValueError("episodes may contain only flat steps")
            p = p - kappa0 + s.potential_error
            zeta += s.potential_error
            free += 1
            if p < -1e-10:
                admissible = False
    costly = max(0, len(episodes) - 1)
    depth = free + costly
    ledger = (costly + 1) * potential_reset + zeta - free * kappa0
    return {
        "episodes": float(len(episodes)),
        "free_steps": float(free),
        "costly_steps": float(costly),
        "depth": float(depth),
        "episode_ledger_margin": ledger,
        "admissible": admissible and ledger >= -1e-10,
    }
