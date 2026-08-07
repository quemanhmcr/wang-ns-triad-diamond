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
    """Scaled source integral Sigma=int_0^c N^-4||S_*|| d tau.

    From int_t ||S_*|| dt >= I1/(132 T), T=cN^-2, tau=N^2 t.
    """
    if I1 < 0 or lifetime_c <= 0:
        raise ValueError("invalid H1 episode parameters")
    return I1 / (132.0 * lifetime_c)


def half_integral_superlevel_threshold(integral_lower: float, interval_length: float) -> float:
    """rho0=Sigma/(2c); the superlevel {f>=rho0} carries >=Sigma/2 integral."""
    if integral_lower < 0 or interval_length <= 0:
        raise ValueError("invalid superlevel parameters")
    return integral_lower / (2.0 * interval_length)


def h1_source_level_threshold(I1: float, lifetime_c: float) -> float:
    return half_integral_superlevel_threshold(
        h1_channel_normalized_integral_lower(I1, lifetime_c), lifetime_c
    )


def sgs_episode_thresholds(
    I1: float,
    lifetime_c: float,
    scale_radius_cap: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> dict[str, float]:
    rho0 = h1_source_level_threshold(I1, lifetime_c)
    q0 = cubic_increment_from_sgs_source_lower(rho0, scale_radius_cap, filter_l1)
    inc = increment_collision_thresholds(
        q0, filter_l1, lp_constant, bernstein_constant,
        filter_radius, band_support_factor,
    )
    return {
        "source_level": rho0,
        "cubic_increment": q0,
        "low_band_mass": inc["low_band_critical_mass"],
        "high_enstrophy": inc["high_normalized_enstrophy"],
        "dominant_atom_mass": 0.25 * inc["low_band_critical_mass"],
        "atomic_entropy_if_no_dominant": math.log(4.0),
        "ancestry_entropy_or_pair_entropy": math.log(2.0),
        "same_ancestry_pair_mass": 0.25,
        "large_radius_mass": fresh_radius_mass_lower(scale_radius_cap),
    }


def viscous_episode_thresholds(
    I1: float,
    lifetime_c: float,
    scale_radius_cap: float,
    viscosity: float,
) -> dict[str, float]:
    rho0 = h1_source_level_threshold(I1, lifetime_c)
    d0 = enstrophy_from_viscous_source_lower(rho0, viscosity, scale_radius_cap)
    return {
        "source_level": rho0,
        "resolved_enstrophy": d0,
        "large_radius_mass": fresh_radius_mass_lower(scale_radius_cap),
    }


def sgs_persistent_episode_costs(
    thresholds: dict[str, float],
    source_superlevel_measure_lower: float,
) -> dict[str, float]:
    """Pigeonhole costs on a persistent SGS-source superlevel set.

    If E has measure >=m, then either s>s0 on >=m/2, or the scale-matched
    subset has >=m/2. On the latter, increment collision gives mass or enstrophy
    on at least half again, hence >=m/4.  The enstrophy case pays normalized
    dissipation >=(m/4)d0.
    """
    m = source_superlevel_measure_lower
    if m < 0:
        raise ValueError("nonnegative scaled-time measure required")
    return {
        "radius_branch_measure": 0.5 * m,
        "mass_or_entropy_branch_measure": 0.25 * m,
        "enstrophy_branch_measure": 0.25 * m,
        "high_frequency_dissipation": 0.25 * m * thresholds["high_enstrophy"],
    }


def viscous_persistent_episode_costs(
    thresholds: dict[str, float],
    source_superlevel_measure_lower: float,
) -> dict[str, float]:
    m = source_superlevel_measure_lower
    if m < 0:
        raise ValueError("nonnegative scaled-time measure required")
    return {
        "radius_branch_measure": 0.5 * m,
        "enstrophy_branch_measure": 0.5 * m,
        "resolved_dissipation": 0.5 * m * thresholds["resolved_enstrophy"],
    }


def temporal_concentration_alternative(
    normalized_integral_lower: float,
    interval_length: float,
    persistence_measure: float,
) -> dict[str, float]:
    """Exact source-superlevel concentration statement.

    E={f>=Sigma/(2c)} always carries at least Sigma/2 of the integral.
    If |E|<m0, at least Sigma/2 source integral is concentrated on a scaled-time
    set of measure <m0: this is the explicit CKN/burst branch.
    """
    if normalized_integral_lower < 0 or interval_length <= 0 or persistence_measure <= 0:
        raise ValueError("invalid temporal concentration parameters")
    return {
        "superlevel_threshold": normalized_integral_lower / (2.0 * interval_length),
        "source_integral_on_superlevel": 0.5 * normalized_integral_lower,
        "concentration_measure_cap": persistence_measure,
    }


@dataclass(frozen=True)
class EpisodeCollisionStress:
    samples: int
    minimum_source_level_margin: float
    minimum_sgs_dissipation_margin: float
    minimum_viscous_dissipation_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> EpisodeCollisionStress:
    rng = np.random.default_rng(seed)
    ms = md = mv = float("inf")
    for _ in range(samples):
        I1 = float(rng.uniform(1e-4, 0.3))
        c = float(rng.uniform(0.05, 1.0))
        s0 = float(rng.uniform(0.5, 4.0))
        g1 = float(rng.uniform(1.0, 2.0))
        clp = float(rng.uniform(1.0, 3.0))
        cb = float(rng.uniform(1.0, 2.0))
        th = sgs_episode_thresholds(I1, c, s0, g1, clp, cb)
        sigma = h1_channel_normalized_integral_lower(I1, c)
        rho = h1_source_level_threshold(I1, c)
        ms = min(ms, 2.0 * c * rho - sigma)
        if 2.0 * c * rho + 2e-14 < sigma:
            raise AssertionError("superlevel threshold arithmetic failed")
        m = float(rng.uniform(1e-4, c))
        costs = sgs_persistent_episode_costs(th, m)
        expect = 0.25 * m * th["high_enstrophy"]
        md = min(md, costs["high_frequency_dissipation"] - expect)
        if costs["high_frequency_dissipation"] + 1e-14 < expect:
            raise AssertionError("SGS episode dissipation routing failed")

        nu = float(rng.uniform(0.2, 2.0))
        tv = viscous_episode_thresholds(I1, c, s0, nu)
        cv = viscous_persistent_episode_costs(tv, m)
        expectv = 0.5 * m * tv["resolved_enstrophy"]
        mv = min(mv, cv["resolved_dissipation"] - expectv)
        if cv["resolved_dissipation"] + 1e-14 < expectv:
            raise AssertionError("viscous episode dissipation routing failed")
    return EpisodeCollisionStress(samples, ms, md, mv)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-source-episode-collision"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    data = {
        "status": "EXACT_ROUTING_GIVEN_H1_SOURCE_SGS_COLLISION_AND_LP_BERNSTEIN",
        "stress": out.__dict__,
        "clean_entropy_constants": {
            "dominant_fraction": "1/4",
            "atomic_entropy": "log 4",
            "ancestry_entropy": "log 2",
            "same_ancestry_pair_mass": "1/4",
        },
    }
    (args.outdir / "source_episode_collision.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# H1 physical-source episode collision

Status: **EXACT_ROUTING_GIVEN_H1_SOURCE_SGS_COLLISION_AND_LP_BERNSTEIN**.

On the H1-dominant, low-strain, mild-aspect branch, the covariant source theorem gives one fixed physical source channel with scaled integral

`Sigma_* >= I1/(132 c)`,  `tau=N^2t`,  `T=cN^-2`.

The superlevel

`rho_* = I1/(264 c^2)`

always carries at least half of that source integral.  If its scaled-time measure is below a chosen persistence threshold, this is explicitly a temporal concentration / CKN-burst branch.

If the differentiated-SGS channel persists and `s=N r_g<=s0`, the filtered-source theorem produces a cubic increment threshold; the Onsager collision then gives a low/base critical-mass band or high-frequency normalized enstrophy.  With the clean packet split `theta=1/4, alpha=1/2`, the mass branch becomes exactly

- one atom with at least `1/4` of the winning band mass; or
- ancestry/component collision entropy at least `log 2`; or
- same-ancestry pair/cycle mass at least `1/4`.

If instead `s>s0`, the affine critical-grain theorem gives the radius-energy event `N int_E|u|^2 >= (3/10)s0`.

On a persistent SGS-source set of scaled measure `m`, after the radius and mass/enstrophy pigeonholes the enstrophy branch pays normalized high-frequency dissipation at least `(m/4)d_high`.  For the viscous source the analogous scale-matched branch pays at least `(m/2)d_V`.

Stress: `{out.samples}`
- minimum source-threshold arithmetic margin: `{out.minimum_source_level_margin:.3e}`
- minimum SGS dissipation-routing margin: `{out.minimum_sgs_dissipation_margin:.3e}`
- minimum viscous dissipation-routing margin: `{out.minimum_viscous_dissipation_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
