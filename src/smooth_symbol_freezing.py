from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path


def sharp_young_constant_3d() -> float:
    return (math.sqrt(3.0) / 2.0) ** 3


def symbol_freezing_error(
    lipschitz_constant: float,
    cell_diameter: float,
    norm_f: float = 1.0,
    norm_g: float = 1.0,
    norm_h: float = 1.0,
) -> float:
    """Sharp-Young bound for freezing a bounded trilinear multiplier on cells.

    If |m-m_h| <= L h pointwise, then
      |T_m-T_mh| <= A_3 L h ||f||_{3/2}||g||_{3/2}||h||_{3/2}.
    """
    if min(lipschitz_constant, cell_diameter, norm_f, norm_g, norm_h) < 0:
        raise ValueError("all parameters must be nonnegative")
    return sharp_young_constant_3d() * lipschitz_constant * cell_diameter * norm_f * norm_g * norm_h


def quadratic_cell_schedule(depth: int, offset: int = 3) -> list[float]:
    if depth < 0 or offset <= 0:
        raise ValueError("invalid schedule")
    return [1.0 / (j + offset) ** 2 for j in range(depth)]


def infinite_quadratic_tail_upper(offset: int = 3) -> float:
    """Elementary integral-test upper bound sum_{j>=0}(j+offset)^-2.

    1/offset^2 + integral_offset^infty x^-2 dx = 1/o^2+1/o.
    """
    if offset <= 0:
        raise ValueError("offset must be positive")
    return 1.0 / offset**2 + 1.0 / offset


@dataclass(frozen=True)
class FreezingCertificate:
    lipschitz_constant: float
    schedule_offset: int
    infinite_cell_diameter_sum_upper: float
    infinite_normalized_transfer_error_upper: float


def freezing_certificate(lipschitz_constant: float, offset: int = 3) -> FreezingCertificate:
    tail = infinite_quadratic_tail_upper(offset)
    return FreezingCertificate(
        lipschitz_constant,
        offset,
        tail,
        sharp_young_constant_3d() * lipschitz_constant * tail,
    )


def stress(samples: int = 50_000, seed: int = 20260807) -> dict[str, float]:
    """Finite-dimensional sanity check of the sup-symbol freezing inequality.

    We use a scale-covariant toy symbol m(x,y)=tanh(a.x+b.y) on a compact
    dimensionless block. Its Euclidean Lipschitz constant is bounded by
    sqrt(|a|^2+|b|^2). Random cells verify the pointwise freezing estimate; the
    theorem itself is the analytic sup bound plus sharp Young, not this probe.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    worst_ratio = 0.0
    worst_transfer_bound = 0.0
    for _ in range(samples):
        a = rng.normal(size=3)
        b = rng.normal(size=3)
        L = float(np.sqrt(np.dot(a, a) + np.dot(b, b)))
        center_x = rng.uniform(-1.0, 1.0, size=3)
        center_y = rng.uniform(-1.0, 1.0, size=3)
        h = float(rng.uniform(1e-5, 0.2))
        # Sample one point within a 6D ball-box whose Euclidean displacement is
        # rescaled to be <= h.
        d = rng.normal(size=6)
        d /= max(np.linalg.norm(d), 1e-30)
        d *= rng.uniform(0.0, h)
        x = center_x + d[:3]
        y = center_y + d[3:]
        m0 = math.tanh(float(np.dot(a, center_x) + np.dot(b, center_y)))
        m1 = math.tanh(float(np.dot(a, x) + np.dot(b, y)))
        denom = L * h
        ratio = abs(m1 - m0) / denom if denom > 0 else 0.0
        if ratio > 1.0 + 1e-12:
            raise AssertionError("Lipschitz freezing estimate failed")
        worst_ratio = max(worst_ratio, ratio)
        worst_transfer_bound = max(worst_transfer_bound, symbol_freezing_error(L, h))
    cert = freezing_certificate(1.0, 3)
    return {
        "samples": samples,
        "worst_pointwise_lipschitz_ratio": worst_ratio,
        "largest_random_normalized_transfer_bound": worst_transfer_bound,
        "young_constant_A3": sharp_young_constant_3d(),
        "quadratic_schedule_sum_upper": cert.infinite_cell_diameter_sum_upper,
        "unit_lipschitz_infinite_transfer_error_upper": cert.infinite_normalized_transfer_error_upper,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-symbol-freezing"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = stress(args.samples)
    (args.outdir / "smooth_symbol_freezing.json").write_text(json.dumps(result, indent=2))
    md = f"""# Smooth SGS symbol freezing\n\n- random compact-cell checks: `{result['samples']}`\n- worst pointwise Lipschitz ratio: `{result['worst_pointwise_lipschitz_ratio']:.9f}`\n- sharp scalar Young constant A3: `{result['young_constant_A3']:.12f}`\n- quadratic relative-cell schedule sum upper: `{result['quadratic_schedule_sum_upper']:.9f}`\n- unit-Lipschitz infinite normalized freezing-error upper: `{result['unit_lipschitz_infinite_transfer_error_upper']:.9f}`\n\nThe theorem is the deterministic inequality `|T_m-T_mh| <= A3 L h` on one\nblock.  The random checks only stress the cellwise Lipschitz bookkeeping.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
