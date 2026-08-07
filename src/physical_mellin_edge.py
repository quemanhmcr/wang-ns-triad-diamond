from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution

from .helical import coupling_magnitude_closed
from .triad_extremizer import symmetric_jstar, symmetric_rstar


def forward_mellin_coefficient(x: float, y: float, sx: int, sy: int, sq: int) -> float:
    """Full signed Mellin flux coefficient for q=1, with phase chosen so child gains energy.

    For modal energy rates dE_j = a_j R and |R| carrying the common amplitude/phase
    factor including |g|, choose sign(R) so dE_q>=0.  The returned coefficient
    multiplies the positive common amplitude factor before |g| normalization.
    """
    if not (0.0 < x <= y < 1.0 and x + y > 1.0):
        return -math.inf
    ak = sy * y - sq
    aq = sx * x - sy * y
    if abs(aq) < 1e-15:
        return -math.inf
    phase_sign = 1.0 if aq > 0.0 else -1.0
    g = coupling_magnitude_closed(x, y, 1.0, sx, sy, sq)
    moment = -ak * math.log(y / x) + aq * math.log(1.0 / y)
    return g * phase_sign * moment


def upper_progress_coefficient(x: float, y: float, sx: int, sy: int, sq: int) -> float:
    if not (0.0 < x <= y < 1.0 and x + y > 1.0):
        return -math.inf
    aq = sx * x - sy * y
    if abs(aq) < 1e-15:
        return -math.inf
    g = coupling_magnitude_closed(x, y, 1.0, sx, sy, sq)
    return g * abs(aq) * math.log(1.0 / y)


def global_search(seed: int = 20260807) -> dict[str, object]:
    best = None
    rows = []
    for idx, signs in enumerate(itertools.product((-1, 1), repeat=3)):
        sx, sy, sq = signs

        def unpack(z: np.ndarray) -> tuple[float, float]:
            y = float(z[0])
            lam = float(z[1])
            x = (1.0 - y) + lam * (2.0 * y - 1.0)
            return x, y

        def objective(z: np.ndarray) -> float:
            x, y = unpack(z)
            val = forward_mellin_coefficient(x, y, sx, sy, sq)
            if not math.isfinite(val):
                return 10.0
            return -val

        res = differential_evolution(
            objective,
            [(0.500000001, 0.999999999), (1e-9, 1.0 - 1e-9)],
            seed=seed + idx,
            maxiter=500,
            popsize=24,
            tol=1e-11,
            polish=True,
            workers=1,
        )
        x, y = unpack(res.x)
        val = forward_mellin_coefficient(x, y, sx, sy, sq)
        row = {
            "signs": list(signs),
            "x": x,
            "y": y,
            "mellin": val,
            "upper": upper_progress_coefficient(x, y, sx, sy, sq),
            "lower_correction": val - upper_progress_coefficient(x, y, sx, sy, sq),
        }
        rows.append(row)
        if best is None or val > float(best["mellin"]):
            best = row
    assert best is not None
    r = symmetric_rstar()
    j = symmetric_jstar(r)
    return {
        "symmetric_rstar": r,
        "symmetric_Jstar": j,
        "best": best,
        "best_over_Jstar": float(best["mellin"]) / j,
        "all_signs": rows,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=Path("results-physical-mellin"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = global_search()
    (args.outdir / "physical_mellin_edge.json").write_text(json.dumps(result, indent=2))
    b = result["best"]
    md = f'''# Full physical Mellin single-edge search

This is a numerical adversarial search, not a theorem certificate.

- old symmetric r*: `{result['symmetric_rstar']:.15f}`
- old J*: `{result['symmetric_Jstar']:.15f}`
- best full Mellin coefficient: `{b['mellin']:.15f}`
- ratio best/J*: `{result['best_over_Jstar']:.9f}`
- best parent ratios: `x={b['x']:.12f}`, `y={b['y']:.12f}`
- best helicities: `{b['signs']}`
- upper progress part there: `{b['upper']:.15f}`
- lower-segment correction there: `{b['lower_correction']:.15f}`

If the ratio exceeds one materially, the old progress functional is not the
extremizer of the full signed Mellin flux and must not be promoted unchanged to
the PDE ledger.
'''
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
