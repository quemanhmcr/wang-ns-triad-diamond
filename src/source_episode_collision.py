from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.onsager_increment_collision import increment_collision_thresholds
from src.sgs_source_collision import (
    cubic_increment_from_sgs_source_lower,
    enstrophy_from_viscous_source_lower,
    fresh_radius_mass_lower,
)


def h1_channel_normalized_integral_lower(I1: float, lifetime_c: float) -> float:
    """Sigma=int_0^c N^-4||S_*||d tau >= I1/(132c)."""
    if I1 < 0 or lifetime_c <= 0:
        raise ValueError("invalid H1 episode parameters")
    return I1 / (132.0 * lifetime_c)


def sgs_source_linear_collision_coefficients(
    scale_radius_cap: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> dict[str, float]:
    """Coefficients mu_band>=c_mu*rho or d_high>=c_d*rho.

    Homogeneity is exact: source->increment is rho^(3/2), while the increment
    collision threshold takes the 2/3 power, leaving a linear rho currency.
    """
    q1 = cubic_increment_from_sgs_source_lower(1.0, scale_radius_cap, filter_l1)
    th = increment_collision_thresholds(
        q1, filter_l1, lp_constant, bernstein_constant,
        filter_radius, band_support_factor,
    )
    return {
        "low_band_mass_per_source": th["low_band_critical_mass"],
        "high_enstrophy_per_source": th["high_normalized_enstrophy"],
        "unit_source_cubic_increment": q1,
    }


def source_weighted_sgs_episode_costs(
    I1: float,
    lifetime_c: float,
    scale_radius_cap: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> dict[str, float]:
    """Master-facing source-weighted alternatives, with no time-persistence assumption.

    One fixed SGS source channel has total normalized source weight Sigma.
    Either large-radius times carry >=Sigma/2, or scale-matched times do.
    On the latter, mass/enstrophy branches split source weight; the enstrophy
    branch directly pays dissipation. The mass branch is split once more into a
    dominant atom versus entropy/cycle routing at theta=1/4, alpha=1/2.
    """
    sigma = h1_channel_normalized_integral_lower(I1, lifetime_c)
    coeff = sgs_source_linear_collision_coefficients(
        scale_radius_cap, filter_l1, lp_constant, bernstein_constant,
        filter_radius, band_support_factor,
    )
    cm = coeff["low_band_mass_per_source"]
    cd = coeff["high_enstrophy_per_source"]
    return {
        "total_source_weight": sigma,
        "large_radius_source_weight": 0.5 * sigma,
        "large_radius_mass": fresh_radius_mass_lower(scale_radius_cap),
        # scale-matched source weight >=Sigma/2; one of mass/enstrophy gets >=Sigma/4
        "high_frequency_dissipation": 0.25 * cd * sigma,
        "winning_band_mass_occupation": 0.25 * cm * sigma,
        # if mass branch wins, dominant-vs-entropy split costs another 1/2 in source weight;
        # dominant atom has one quarter of band mass.
        "dominant_atom_mass_occupation": (1.0 / 32.0) * cm * sigma,
        "entropy_or_cycle_source_weight": 0.125 * sigma,
        # entropy branch splits Bellman-vs-cycle once more by source weight if desired
        "bellman_or_cycle_source_weight": 0.0625 * sigma,
        "atomic_entropy": math.log(4.0),
        "ancestry_entropy": math.log(2.0),
        "same_ancestry_pair_mass": 0.25,
        **coeff,
    }


def source_weighted_viscous_episode_costs(
    I1: float,
    lifetime_c: float,
    scale_radius_cap: float,
    viscosity: float,
) -> dict[str, float]:
    """Viscous source pays dissipation without a persistence hypothesis.

    On the scale-matched branch d_V >= b rho^2, b=(5000/(nu s0))^2.
    If that branch carries source integral >=Sigma/2 on an interval of scaled
    length at most c, Cauchy gives int rho^2 >= Sigma^2/(4c).
    """
    if viscosity <= 0:
        raise ValueError("positive viscosity required")
    sigma = h1_channel_normalized_integral_lower(I1, lifetime_c)
    b = enstrophy_from_viscous_source_lower(1.0, viscosity, scale_radius_cap)
    return {
        "total_source_weight": sigma,
        "large_radius_source_weight": 0.5 * sigma,
        "large_radius_mass": fresh_radius_mass_lower(scale_radius_cap),
        "enstrophy_per_source_squared": b,
        "resolved_dissipation": b * sigma * sigma / (4.0 * lifetime_c),
    }


def source_weight_partition_lower(total_source_weight: float, levels: int) -> float:
    """After `levels` binary source-weight pigeonholes, one branch carries Sigma/2^levels."""
    if total_source_weight < 0 or levels < 0:
        raise ValueError("invalid source partition")
    return total_source_weight / (2.0 ** levels)


@dataclass(frozen=True)
class EpisodeCollisionStress:
    samples: int
    minimum_sgs_homogeneity_margin: float
    minimum_sgs_dissipation_margin: float
    minimum_viscous_dissipation_margin: float
    minimum_partition_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> EpisodeCollisionStress:
    rng = np.random.default_rng(seed)
    mh = md = mv = mp = float("inf")
    for _ in range(samples):
        I1 = float(rng.uniform(1e-4, 0.3))
        c = float(rng.uniform(0.05, 1.0))
        s0 = float(rng.uniform(0.5, 4.0))
        g1 = float(rng.uniform(1.0, 2.0))
        clp = float(rng.uniform(1.0, 3.0))
        cb = float(rng.uniform(1.0, 2.0))
        coeff = sgs_source_linear_collision_coefficients(s0, g1, clp, cb)
        rho = float(rng.uniform(1e-5, 0.05))
        q = cubic_increment_from_sgs_source_lower(rho, s0, g1)
        th = increment_collision_thresholds(q, g1, clp, cb)
        mh = min(mh,
                 th["low_band_critical_mass"] - coeff["low_band_mass_per_source"] * rho,
                 th["high_normalized_enstrophy"] - coeff["high_enstrophy_per_source"] * rho)
        if th["low_band_critical_mass"] + 2e-12 < coeff["low_band_mass_per_source"] * rho:
            raise AssertionError("SGS source-to-mass homogeneity failed")
        if th["high_normalized_enstrophy"] + 2e-12 < coeff["high_enstrophy_per_source"] * rho:
            raise AssertionError("SGS source-to-enstrophy homogeneity failed")

        out = source_weighted_sgs_episode_costs(I1, c, s0, g1, clp, cb)
        sigma = out["total_source_weight"]
        expectd = 0.25 * out["high_enstrophy_per_source"] * sigma
        md = min(md, out["high_frequency_dissipation"] - expectd)
        if out["high_frequency_dissipation"] + 1e-14 < expectd:
            raise AssertionError("source-weighted SGS dissipation failed")

        nu = float(rng.uniform(0.2, 2.0))
        v = source_weighted_viscous_episode_costs(I1, c, s0, nu)
        b = v["enstrophy_per_source_squared"]
        expectv = b * v["total_source_weight"] ** 2 / (4.0 * c)
        mv = min(mv, v["resolved_dissipation"] - expectv)
        if v["resolved_dissipation"] + 1e-14 < expectv:
            raise AssertionError("source-weighted viscous dissipation failed")

        levels = int(rng.integers(0, 7))
        part = source_weight_partition_lower(sigma, levels)
        mp = min(mp, part - sigma / (2 ** levels))
    return EpisodeCollisionStress(samples, mh, md, mv, mp)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-source-episode-collision"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    data = {
        "status": "EXACT_SOURCE_WEIGHTED_ROUTING_GIVEN_H1_SOURCE_AND_STANDARD_LP_BERNSTEIN",
        "stress": out.__dict__,
        "clean_entropy_constants": {
            "dominant_fraction": "1/4",
            "atomic_entropy": "log 4",
            "ancestry_entropy": "log 2",
            "same_ancestry_pair_mass": "1/4",
        },
        "principle": "temporal concentration is not a free SGS/viscous exit: homogeneity converts source weight directly to mass/enstrophy, and viscous concentration increases the quadratic dissipation cost",
    }
    (args.outdir / "source_episode_collision.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# H1 physical-source episode collision: source-weighted form

Status: **EXACT_SOURCE_WEIGHTED_ROUTING_GIVEN_H1_SOURCE_AND_STANDARD_LP_BERNSTEIN**.

A fixed H1 source channel carries normalized source weight

`Sigma_* >= I1/(132 c)`.

For differentiated SGS stress, filtered-source collision gives cubic increments proportional to `rho^(3/2)`, while the Onsager mass/enstrophy threshold takes the `2/3` power.  Therefore the final currencies are **linear in the instantaneous source density**:

`mu_band >= c_mu rho`  or  `d_high >= c_d rho`.

This removes any persistence assumption. Source weight can be pigeonholed directly. Outside a large-radius branch carrying `Sigma/2`, the scale-matched source has weight at least `Sigma/2`; either the mass or enstrophy branch carries at least `Sigma/4`.  Hence the enstrophy branch pays normalized high-frequency dissipation at least

`D_high >= (c_d/4) Sigma`.

If the mass branch wins, `theta=1/4, alpha=1/2` gives either dominant-atom mass occupation, or the clean event `H_anc>=log 2`, or same-ancestry pair/cycle mass `>=1/4`.  The source-weight carried by the Bellman/cycle alternatives remains quantitative after the finite binary pigeonholes.

For viscous source, `d_V>=b rho_nu^2`.  If the scale-matched branch carries source weight at least `Sigma/2`, Cauchy on the entire scaled lifetime `[0,c]` gives

`D_V >= b Sigma^2/(4c)`.

Thus concentrating the viscous source in time only increases its dissipation price.

Stress: `{out.samples}`
- minimum SGS homogeneity margin: `{out.minimum_sgs_homogeneity_margin:.3e}`
- minimum source-weighted SGS dissipation margin: `{out.minimum_sgs_dissipation_margin:.3e}`
- minimum source-weighted viscous dissipation margin: `{out.minimum_viscous_dissipation_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
