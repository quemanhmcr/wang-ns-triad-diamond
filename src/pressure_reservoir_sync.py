from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.ancestor_reservoir_sync import CLEAN_GENERATION_PROGRESS, CLEAN_RESERVOIR_GROWTH


def pressure_hessian_pair_energy_service_ratio_upper() -> Fraction:
    """Per-unit sqrt(E_a E_b) objective pressure-Hessian pair ratio.

    The resolved objective source is `N^-4 ||grad^2 P||_inf`.  For one
    low--low pair, order-two differentiation plus L^(3/2)->Linf Bernstein gives
    the scale coefficient

      M_max^4 sqrt(M_a M_b) / N^4

    per unit `sqrt(E_a E_b)` (fixed Riesz/Bernstein constants cancel between
    generations).  Reservoir frequencies grow by <=21/20 and the signed-good
    block scale by >=8/5, hence the clean ratio below.
    """
    return CLEAN_RESERVOIR_GROWTH**5 / CLEAN_GENERATION_PROGRESS**4


def pressure_hessian_pair_service_coefficient(
    frequency_a: float,
    frequency_b: float,
    block_frequency: float,
) -> float:
    if min(frequency_a, frequency_b, block_frequency) <= 0:
        raise ValueError("positive frequencies required")
    mmax = max(frequency_a, frequency_b)
    return mmax**4 * math.sqrt(frequency_a * frequency_b) / block_frequency**4


def pressure_hessian_pair_service_capacity(
    frequency_a: float,
    frequency_b: float,
    block_frequency: float,
    energy_a: float,
    energy_b: float,
) -> float:
    if min(energy_a, energy_b) < 0:
        raise ValueError("nonnegative energies required")
    return pressure_hessian_pair_service_coefficient(
        frequency_a, frequency_b, block_frequency
    ) * math.sqrt(energy_a * energy_b)


def amortized_pressure_hessian_pair_capacity_upper(
    generation: int,
    frequency_a0: float,
    frequency_b0: float,
    block_frequency0: float,
    global_energy_cap: float,
) -> float:
    if generation < 0 or global_energy_cap < 0:
        raise ValueError("invalid pressure-Hessian pair data")
    base = (
        pressure_hessian_pair_service_coefficient(
            frequency_a0, frequency_b0, block_frequency0
        )
        * global_energy_cap
    )
    return base * float(pressure_hessian_pair_energy_service_ratio_upper()) ** generation


def total_pressure_hessian_pair_capacity_upper(
    frequency_a0: float,
    frequency_b0: float,
    block_frequency0: float,
    global_energy_cap: float,
) -> float:
    base = (
        pressure_hessian_pair_service_coefficient(
            frequency_a0, frequency_b0, block_frequency0
        )
        * global_energy_cap
    )
    r = float(pressure_hessian_pair_energy_service_ratio_upper())
    return base / (1.0 - r)


def pair_energy_service_ratio_upper() -> Fraction:
    """Per-unit sqrt(E_a E_b) **pressure-third** service ratio for one reused pair.

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
    rh = pressure_hessian_pair_energy_service_ratio_upper()
    r3 = pair_energy_service_ratio_upper()
    if not (rh < Fraction(1, 5)):
        raise AssertionError("objective pressure-Hessian pair did not have one-fifth-life")
    if not (Fraction(1, 1) / (Fraction(1, 1) - rh) < Fraction(5, 4)):
        raise AssertionError("objective pressure-Hessian total pair capacity exceeded 5/4 base")
    if not (r3 < Fraction(1, 3)):
        raise AssertionError("pressure-third reservoir pair did not have one-third-life")
    if not (Fraction(1, 1) / (Fraction(1, 1) - r3) < Fraction(3, 2)):
        raise AssertionError("pressure-third total pair capacity exceeded 3/2 base")
    return {
        "hessian_pair_service_ratio": f"{rh.numerator}/{rh.denominator}",
        "hessian_pair_service_clean": "<1/5",
        "hessian_total_future_pair_capacity": "<5/4 times generation-0 pair energy-capacity coefficient",
        "pressure_third_pair_service_ratio": f"{r3.numerator}/{r3.denominator}",
        "pressure_third_pair_service_clean": "<1/3",
        "pressure_third_total_future_pair_capacity": "<3/2 times generation-0 pair energy-capacity coefficient",
        "status": "EXACT_OBJECTIVE_PRESSURE_HESSIAN_AND_PRESSURE_THIRD_PAIR_EROSION_GIVEN_SIGNED_GOOD_LOW_STRAIN_LINEAGE",
    }


@dataclass(frozen=True)
class PressurePairStress:
    samples: int
    maximum_hessian_pair_ratio: float
    minimum_one_fifth_margin: float
    minimum_hessian_amortized_margin: float
    maximum_pressure_third_pair_ratio: float
    minimum_one_third_margin: float
    minimum_pressure_third_amortized_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> PressurePairStress:
    rng = np.random.default_rng(seed)
    mrh = mr3 = 0.0
    mh = mah = m3 = ma3 = float("inf")
    g = float(CLEAN_RESERVOIR_GROWTH)
    r = float(CLEAN_GENERATION_PROGRESS)
    exact_h = float(pressure_hessian_pair_energy_service_ratio_upper())
    exact_3 = float(pair_energy_service_ratio_upper())
    for _ in range(samples):
        ga = float(rng.uniform(0.2, g * (1.0 - 1e-10)))
        gb = float(rng.uniform(0.2, g * (1.0 - 1e-10)))
        gn = float(rng.uniform(r * (1.0 + 1e-10), 3.0))

        # Objective pressure Hessian: Mmax^4 sqrt(Ma Mb) / N^4.
        actual_h = max(ga, gb) ** 4 * math.sqrt(ga * gb) / gn**4
        mrh = max(mrh, actual_h)
        mh = min(mh, 1.0 / 5.0 - actual_h)
        if actual_h >= 1.0 / 5.0 + 1e-12:
            raise AssertionError("objective pressure-Hessian pair one-fifth-life failed")

        # H1 pressure-third: Mmax^3 sqrt(Ma Mb) / N^3.
        actual_3 = max(ga, gb) ** 3 * math.sqrt(ga * gb) / gn**3
        mr3 = max(mr3, actual_3)
        m3 = min(m3, 1.0 / 3.0 - actual_3)
        if actual_3 >= 1.0 / 3.0 + 1e-12:
            raise AssertionError("pressure-third pair one-third-life failed")

        ma0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        mb0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        n0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        E = float(math.exp(rng.uniform(-3.0, 2.0)))
        q = int(rng.integers(0, 16))

        cap_h = amortized_pressure_hessian_pair_capacity_upper(q, ma0, mb0, n0, E)
        fifth = pressure_hessian_pair_service_coefficient(ma0, mb0, n0) * E * (1.0 / 5.0) ** q
        mah = min(mah, fifth - cap_h)
        if cap_h > fifth + 2e-12 * max(1.0, fifth):
            raise AssertionError("objective pressure-Hessian clean one-fifth envelope failed")

        cap_3 = amortized_pressure_pair_capacity_upper(q, ma0, mb0, n0, E)
        third = pressure_pair_service_coefficient(ma0, mb0, n0) * E * (1.0 / 3.0) ** q
        ma3 = min(ma3, third - cap_3)
        if cap_3 > third + 2e-12 * max(1.0, third):
            raise AssertionError("pressure-third clean one-third envelope failed")

    # Exact clean constants must dominate all adversarial sampled ratios.
    if exact_h >= 1.0 / 5.0 or exact_3 >= 1.0 / 3.0:
        raise AssertionError("exact pressure pair clean ratio lost")
    return PressurePairStress(samples, mrh, mh, mah, mr3, m3, ma3)


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
    md = f"""# Pressure reservoir synchronization: derivative order sets the pair lifetime

Status: **{cert['status']}**.

There are two distinct pressure source objects and they must not be conflated.

For the **resolved objective pressure Hessian** `N^-4||grad^2 P||_inf`, one low--low pair `(a,b)` has, up to the fixed Riesz/Bernstein constant, per-unit `sqrt(E_a E_b)` coefficient

`M_max^4 sqrt(M_a M_b) / N^4`.

On a signed-good low-strain lineage, reservoir frequencies grow by less than `21/20` while the block scale grows by more than `8/5`.  Hence the objective-Hessian pair coefficient contracts by

`(21/20)^5 (5/8)^4 = {pressure_hessian_pair_energy_service_ratio_upper().numerator}/{pressure_hessian_pair_energy_service_ratio_upper().denominator} < 1/5`.

Even allowing both reservoirs the entire global energy cap at every future generation, one fixed objective-Hessian pair has total future capacity less than `5/4` times its generation-zero capacity.

For the separate **H1 pressure-third** source, the coefficient is

`M_max^3 sqrt(M_a M_b) / N^3`,

and its previously certified contraction remains

`(21/20)^4 (5/8)^3 = {pair_energy_service_ratio_upper().numerator}/{pair_energy_service_ratio_upper().denominator} < 1/3`,

with total future fixed-pair capacity `<3/2` of generation zero.

Thus the derivative order is physical provenance: objective Hessian and pressure-third share the same material-pair erosion mechanism but have different exact lifetimes.  Persistent objective pressure service must relink pairs, fragment into pair/component entropy, leave the supplied signed-good low-strain lineage, or use its SGS-stress alternative.

Stress: `{out.samples}`
- maximum sampled objective-Hessian pair ratio: `{out.maximum_hessian_pair_ratio:.9f}`
- minimum margin below `1/5`: `{out.minimum_one_fifth_margin:.3e}`
- minimum objective-Hessian one-fifth envelope margin: `{out.minimum_hessian_amortized_margin:.3e}`
- maximum sampled pressure-third pair ratio: `{out.maximum_pressure_third_pair_ratio:.9f}`
- minimum margin below `1/3`: `{out.minimum_one_third_margin:.3e}`
- minimum pressure-third one-third envelope margin: `{out.minimum_pressure_third_amortized_margin:.3e}`
"""

    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
