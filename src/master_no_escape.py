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


def spherical_cap_mass_lower_bound(bnorm: float, alpha: float) -> float:
    if not (0.0 <= bnorm <= 1.0):
        raise ValueError("bnorm must be in [0,1]")
    if not (0.0 < alpha < math.pi):
        raise ValueError("alpha must lie in (0,pi)")
    eta = 1.0 - bnorm
    denom = 1.0 - math.cos(alpha)
    return max(0.0, 1.0 - eta / denom)


def balanced_entropy_lower_bound(eta: float) -> float:
    if not (0.0 <= eta < 1.0):
        raise ValueError("eta must be in [0,1)")
    return math.log(2.0 / (2.0 - eta))


def cross_penalty(cost: float, eta: float) -> float:
    if cost < 0 or eta < 0:
        raise ValueError("nonnegative cost and error required")
    return math.log1p(eta * math.exp(cost))


def master_lower_bound(
    steps: Iterable[MasterStep], *, c0: float, kappa0: float, potential0: float
) -> dict[str, float]:
    steps = list(steps)
    if c0 <= 0 or kappa0 <= 0 or potential0 < 0:
        raise ValueError("invalid master constants")
    zeta = sum(s.potential_error for s in steps if s.low_cost_flat)
    xi = sum(cross_penalty(s.cost, s.cross_error) for s in steps)
    free_bound = (potential0 + zeta) / kappa0
    lower = c0 * (len(steps) - free_bound) - xi
    return {
        "depth": float(len(steps)),
        "potential_error": zeta,
        "cross_penalty": xi,
        "free_block_bound": free_bound,
        "log_efficiency_lower_bound": lower,
        "efficiency_upper_bound": math.exp(-lower),
    }


def verify_trace(
    steps: Iterable[MasterStep], *, c0: float, kappa0: float, potential0: float
) -> dict[str, float | bool]:
    steps = list(steps)
    p = potential0
    actual_log_loss = 0.0
    costly = 0
    free = 0
    for s in steps:
        if not (0.0 < s.efficiency <= 1.0 + s.cross_error + 1e-12):
            raise ValueError("invalid efficiency")
        actual_log_loss += -math.log(min(1.0, s.efficiency))
        if s.low_cost_flat:
            free += 1
            p = max(0.0, p - kappa0 + s.potential_error)
        else:
            costly += 1
            if s.cost + 1e-12 < c0:
                raise ValueError("costly step below c0")
        ideal_upper = math.exp(-s.cost) + s.cross_error
        if s.efficiency > ideal_upper + 1e-10:
            raise ValueError("step violates its ideal block bound")
    bound = master_lower_bound(steps, c0=c0, kappa0=kappa0, potential0=potential0)
    # The abstract theorem only allows as many low-cost steps as the *untruncated*
    # potential ledger permits.
    zeta = sum(s.potential_error for s in steps if s.low_cost_flat)
    admissible = free * kappa0 <= potential0 + zeta + 1e-12
    return {
        **bound,
        "actual_log_loss": actual_log_loss,
        "costly_steps": float(costly),
        "free_steps": float(free),
        "trace_admissible": admissible,
        "master_margin": actual_log_loss - float(bound["log_efficiency_lower_bound"]),
    }
