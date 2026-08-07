from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from .triad_extremizer import symmetric_gamma, symmetric_rstar


@dataclass(frozen=True)
class MellinFluxSegments:
    lower: float
    upper: float
    total: float


def sharp_cutoff_triad_flux(
    K: float,
    k: float,
    p: float,
    q: float,
    dE_k: float,
    dE_p: float,
    dE_q: float,
) -> float:
    """Outward nonlinear energy flux of one ordered conservative triad.

    The sharp low-pass energy is the sum of modal energies with magnitude <=K,
    and outward flux is minus its nonlinear time derivative. Boundary values
    are immaterial for the logarithmic K integral.
    """
    if not (0.0 < k <= p <= q):
        raise ValueError("require 0 < k <= p <= q")
    scale = max(1.0, abs(dE_k), abs(dE_p), abs(dE_q))
    if abs(dE_k + dE_p + dE_q) > 1e-12 * scale:
        raise ValueError("triad energy rates must sum to zero")
    if K < k:
        return 0.0
    if K < p:
        return -dE_k
    if K < q:
        return dE_q  # = -(dE_k+dE_p) by energy conservation
    return 0.0


def mellin_flux_segments(
    k: float,
    p: float,
    q: float,
    dE_k: float,
    dE_p: float,
    dE_q: float,
) -> MellinFluxSegments:
    """Exact integral of the sharp-cutoff triad flux against dK/K."""
    if not (0.0 < k <= p <= q):
        raise ValueError("require 0 < k <= p <= q")
    scale = max(1.0, abs(dE_k), abs(dE_p), abs(dE_q))
    if abs(dE_k + dE_p + dE_q) > 1e-12 * scale:
        raise ValueError("triad energy rates must sum to zero")
    lower = -dE_k * math.log(p / k)
    upper = dE_q * math.log(q / p)
    return MellinFluxSegments(lower, upper, lower + upper)


def maximizing_orbit_rate_coefficients(x: float, y: float, child_sign: int = -1) -> tuple[float, float, float]:
    """Energy-rate coefficients for parents (+,-) and child sign ±.

    Child magnitude is one, 0<x<=y<1.  The common phase/amplitude factor R is
    omitted: dE_j = coefficient_j * R.  For R>0 the child gains energy.
    """
    if child_sign not in (-1, 1):
        raise ValueError("child_sign must be ±1")
    if not (0.0 < x <= y < 1.0 and x + y > 1.0):
        raise ValueError("require an ordered forward triangle")
    sx, sy, sq = 1, -1, child_sign
    return sy * y - sq, sq - sx * x, sx * x - sy * y


def maximizing_orbit_mellin_coefficients(x: float, y: float, child_sign: int = -1) -> MellinFluxSegments:
    """Mellin flux coefficients after factoring out the common triad R."""
    rates = maximizing_orbit_rate_coefficients(x, y, child_sign)
    return mellin_flux_segments(x, y, 1.0, *rates)


def adverse_lower_to_upper_ratio(x: float, y: float) -> float:
    """Magnitude of the lower backscatter segment / upper forward segment.

    This is the child_sign=-1 maximizing orbit.  It is the adverse case; for
    child_sign=+1 the lower segment is forward and adds to the Mellin flux.
    """
    row = maximizing_orbit_mellin_coefficients(x, y, -1)
    if row.upper <= 0.0:
        raise ValueError("upper forward segment must be positive")
    return -row.lower / row.upper


def local_retention_stress(samples: int = 10000) -> dict[str, float]:
    """Deterministic grid stress only; rigorous 0.9 retention is Arb-certified."""
    r = symmetric_rstar()
    gamma = symmetric_gamma(r)
    worst = 1.0
    worst_ratio = 0.0
    n = max(2, int(math.sqrt(samples)))
    for i in range(n):
        u = 0.08 * i / (n - 1)
        for j in range(n):
            v = -0.08 + 0.16 * j / (n - 1)
            R = r * math.exp(-v)
            x = R * math.exp(-u / 2.0)
            y = R * math.exp(u / 2.0)
            ratio = adverse_lower_to_upper_ratio(x, y)
            worst_ratio = max(worst_ratio, ratio)
            worst = min(worst, 1.0 - ratio)
    return {"grid_points": n * n, "worst_adverse_ratio": worst_ratio, "minimum_retention": worst}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("results-log-scale-flux"))
    ap.add_argument("--samples", type=int, default=10000)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    r = symmetric_rstar()
    rates = maximizing_orbit_rate_coefficients(r, r, -1)
    symmetric = mellin_flux_segments(r, r, 1.0, *rates)
    stress = local_retention_stress(args.samples)
    result = {"rstar": r, "gamma": symmetric_gamma(r), "symmetric_mellin": asdict(symmetric), "stress": stress}
    (args.outdir / "log_scale_flux.json").write_text(json.dumps(result, indent=2))
    md = f'''# Sharp-cutoff log-scale flux bridge

For one ordered conservative triad `k<=p<=q`,

`Integral Pi_K dK/K = -dE_k log(p/k) + dE_q log(q/p)`.

At equal parent scales `k=p=r*`, the lower segment vanishes exactly and the
upper segment is the logarithmic progress factor used by the single-edge
functional.

- r*: `{r:.15f}`
- gamma*: `{symmetric_gamma(r):.15f}`
- symmetric lower segment: `{symmetric.lower:.3e}` times the common triad factor
- symmetric upper segment: `{symmetric.upper:.12f}` times the common triad factor
- deterministic local grid points: `{stress['grid_points']}`
- worst adverse lower/upper ratio on grid: `{stress['worst_adverse_ratio']:.9f}`
- minimum full/upper retention on grid: `{stress['minimum_retention']:.9f}`

The 90% retention theorem is certified separately by Arb in the single-edge
certificate. The grid values above are regression evidence only.
'''
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
