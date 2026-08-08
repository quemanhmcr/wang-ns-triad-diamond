from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.ancestor_reservoir_sync import CLEAN_GENERATION_PROGRESS, CLEAN_RESERVOIR_GROWTH


def pair_energy_service_ratio_upper() -> Fraction:
    """Per-unit sqrt(E_a E_b) pressure-third service ratio for one reused pair.

    A pair (M_a,M_b) contributes, before the final 3/2 power,
      M_max^3 sqrt(M_a M_b) sqrt(E_a E_b) / N^3.
    Each reservoir frequency grows by <=21/20 on a low-strain step while the
    lineage scale grows by >=8/5.
    """
    return CLEAN_RESERVOIR_GROWTH**4 / CLEAN_GENERATION_PROGRESS**3


def pressure_pair_service_coefficient(
    frequency_a: float,
    frequency_b: float,
    block_frequency: float,
) -> float:
    if min(frequency_a, frequency_b, block_frequency) <= 0:
        raise ValueError("positive frequencies required")
    mmax = max(frequency_a, frequency_b)
    return mmax**3 * math.sqrt(frequency_a * frequency_b) / block_frequency**3


def pressure_pair_service_capacity(
    frequency_a: float,
    frequency_b: float,
    block_frequency: float,
    energy_a: float,
    energy_b: float,
) -> float:
    if min(energy_a, energy_b) < 0:
        raise ValueError("nonnegative energies required")
    return pressure_pair_service_coefficient(frequency_a, frequency_b, block_frequency) * math.sqrt(energy_a * energy_b)


def amortized_pressure_pair_capacity_upper(
    generation: int,
    frequency_a0: float,
    frequency_b0: float,
    block_frequency0: float,
    global_energy_cap: float,
) -> float:
    if generation < 0 or global_energy_cap < 0:
        raise ValueError("invalid pressure-pair data")
    base = pressure_pair_service_coefficient(frequency_a0, frequency_b0, block_frequency0) * global_energy_cap
    return base * float(pair_energy_service_ratio_upper()) ** generation


def total_pressure_pair_capacity_upper(
    frequency_a0: float,
    frequency_b0: float,
    block_frequency0: float,
    global_energy_cap: float,
) -> float:
    base = pressure_pair_service_coefficient(frequency_a0, frequency_b0, block_frequency0) * global_energy_cap
    r = float(pair_energy_service_ratio_upper())
    return base / (1.0 - r)


def exact_pressure_pair_certificate() -> dict[str, str]:
    r = pair_energy_service_ratio_upper()
    if not (r < Fraction(1, 3)):
        raise AssertionError("pressure reservoir pair did not have one-third-life")
    if not (Fraction(1, 1) / (Fraction(1, 1) - r) < Fraction(3, 2)):
        raise AssertionError("total pressure pair capacity exceeded 3/2 base capacities")
    return {
        "pair_service_ratio": f"{r.numerator}/{r.denominator}",
        "clean_pair_service_ratio": "<1/3",
        "total_future_pair_capacity": "<3/2 times generation-0 pair energy-capacity coefficient",
        "status": "EXACT_GIVEN_SIGNED_GOOD_PROGRESS_AND_LOW_STRAIN_KELVIN_GROWTH",
    }


@dataclass(frozen=True)
class PressurePairStress:
    samples: int
    maximum_pair_ratio: float
    minimum_one_third_margin: float
    minimum_amortized_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> PressurePairStress:
    rng = np.random.default_rng(seed)
    mr = 0.0
    mm = ma = float("inf")
    g = float(CLEAN_RESERVOIR_GROWTH)
    r = float(CLEAN_GENERATION_PROGRESS)
    exact = float(pair_energy_service_ratio_upper())
    for _ in range(samples):
        ga = float(rng.uniform(0.2, g * (1.0 - 1e-10)))
        gb = float(rng.uniform(0.2, g * (1.0 - 1e-10)))
        gn = float(rng.uniform(r * (1.0 + 1e-10), 3.0))
        # Worst max-frequency growth is max(ga,gb)^3 and sqrt-product adds sqrt(ga gb).
        actual = max(ga, gb) ** 3 * math.sqrt(ga * gb) / gn**3
        mr = max(mr, actual)
        mm = min(mm, 1.0 / 3.0 - actual)
        if actual >= 1.0 / 3.0 + 1e-12:
            raise AssertionError("pressure pair one-third-life failed")

        ma0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        mb0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        n0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        E = float(math.exp(rng.uniform(-3.0, 2.0)))
        q = int(rng.integers(0, 16))
        cap = amortized_pressure_pair_capacity_upper(q, ma0, mb0, n0, E)
        third = pressure_pair_service_coefficient(ma0, mb0, n0) * E * (1.0 / 3.0) ** q
        ma = min(ma, third - cap)
        if cap > third + 2e-12 * max(1.0, third):
            raise AssertionError("pressure pair clean one-third envelope failed")
    return PressurePairStress(samples, mr, mm, ma)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-pressure-reservoir-sync"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = exact_pressure_pair_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "pressure_reservoir_sync.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Pressure reservoir synchronization: a low-low pair has one-third-life

Status: **{cert['status']}**.

For a fixed materially reused low-frequency pair `(a,b)`, the filtered pressure-third source generated by `V_a V_b` has normalized pre-`3/2` service coefficient

`M_max^3 sqrt(M_a M_b) / N^3`

per unit `sqrt(E_a E_b)`.  On a signed-good low-strain lineage each reservoir frequency grows by less than `21/20`, while the block scale advances by more than `8/5`.  Hence the coefficient contracts by

`(21/20)^4 (5/8)^3 = 194481/655360 < 1/3`.

Even allowing both reservoirs to carry the whole global energy cap at every service time, one fixed low-low pair has total future pressure service less than `3/2` times its generation-0 pair energy-capacity coefficient.

Persistent pressure-third service must therefore relink to new reservoir pairs, fragment over many pairs (an atomic/component entropy problem), or leave the low-strain material-reuse branch.  Together with band-limited source sampling, this connects the pressure near field to the same ancestry architecture rather than leaving a global low-pass `L^3` reservoir as a black box.

Stress: `{out.samples}`
- maximum sampled pair ratio: `{out.maximum_pair_ratio:.9f}`
- minimum margin below `1/3`: `{out.minimum_one_third_margin:.3e}`
- minimum clean one-third envelope margin: `{out.minimum_amortized_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
