from __future__ import annotations

import math

# The certified Arb theorem uses the corresponding exact rational endpoints.
RSTAR_LO = 0.61090410158
RSTAR_HI = 0.61090410160


def critical_equation(r: float) -> float:
    """Critical equation f(r)=-log(r)-4r^2+1 for the symmetric edge."""
    return -math.log(r) - 4.0 * r * r + 1.0


def symmetric_rstar() -> float:
    """Double-precision r* obtained from the exact monotone critical equation.

    The theorem-level containment of the root in [RSTAR_LO,RSTAR_HI] is proved
    by the Arb certificate; this helper merely avoids stale optimizer decimals
    in numerical/regression modules.
    """
    lo, hi = RSTAR_LO, RSTAR_HI
    if not (critical_equation(lo) > 0.0 and critical_equation(hi) < 0.0):
        raise RuntimeError("certified r* bracket no longer straddles the root")
    for _ in range(64):
        mid = 0.5 * (lo + hi)
        if critical_equation(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def symmetric_jstar(r: float | None = None) -> float:
    if r is None:
        r = symmetric_rstar()
    return math.sqrt(4.0 * r * r - 1.0) * math.log(1.0 / r) / (4.0 * math.sqrt(2.0) * r)


def symmetric_gamma(r: float | None = None) -> float:
    if r is None:
        r = symmetric_rstar()
    return math.log(1.0 / r)
