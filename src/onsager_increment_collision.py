from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.atomic_component_entropy import collision_chain


def coarse_increment_weight_upper(
    low_mass_max: float,
    high_normalized_enstrophy: float,
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> float:
    """Upper bound for the LP square-function mass entering cubic increments.

    For dyadic q_j=2^j and |r|<=R_G/N,
      sum_j min(4,(beta R_G q_j)^2) mu_j
      <= (4/3)(beta R_G)^2 max_{j<=0} mu_j + 2 sum_{j>=1}q_j mu_j.
    """
    if min(low_mass_max, high_normalized_enstrophy) < 0:
        raise ValueError("nonnegative inputs required")
    if filter_radius <= 0 or band_support_factor <= 0:
        raise ValueError("positive geometry constants required")
    a = (4.0 / 3.0) * (band_support_factor * filter_radius) ** 2
    return a * low_mass_max + 2.0 * high_normalized_enstrophy


def cubic_increment_upper(
    low_mass_max: float,
    high_normalized_enstrophy: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> float:
    """Physical cubic increment upper bound given standard LP/Bernstein inputs."""
    if filter_l1 < 0 or lp_constant <= 0 or bernstein_constant <= 0:
        raise ValueError("invalid analytic constants")
    w = coarse_increment_weight_upper(
        low_mass_max, high_normalized_enstrophy, filter_radius, band_support_factor
    )
    c = lp_constant * bernstein_constant
    return filter_l1 * c ** 3 * w ** 1.5


def increment_collision_thresholds(
    cubic_charge: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> dict[str, float]:
    """Mass/enstrophy thresholds forced by a cubic increment event.

    If Q >= cubic_charge, then either
      mu_low_max >= X/(2a), or d_high >= X/4,
    where X=(Q/[g1(C_LP C_B)^3])^(2/3) and a=4(beta R_G)^2/3.
    """
    if cubic_charge <= 0 or filter_l1 <= 0 or lp_constant <= 0 or bernstein_constant <= 0:
        raise ValueError("positive collision parameters required")
    geom = band_support_factor * filter_radius
    if geom <= 0:
        raise ValueError("positive increment geometry required")
    c = lp_constant * bernstein_constant
    X = (cubic_charge / (filter_l1 * c ** 3)) ** (2.0 / 3.0)
    a = (4.0 / 3.0) * geom ** 2
    return {
        "square_mass_threshold": X,
        "low_band_critical_mass": X / (2.0 * a),
        "high_normalized_enstrophy": X / 4.0,
    }


def exact_dyadic_square_weight(
    masses: dict[int, float],
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> float:
    if filter_radius <= 0 or band_support_factor <= 0:
        raise ValueError("positive geometry required")
    total = 0.0
    for j, mu in masses.items():
        if mu < 0:
            raise ValueError("critical masses must be nonnegative")
        q = 2.0 ** int(j)
        m2 = min(4.0, (band_support_factor * filter_radius * q) ** 2)
        total += m2 * mu
    return total


def coarse_dyadic_bound_from_masses(
    masses: dict[int, float],
    filter_radius: float = 1.0,
    band_support_factor: float = 1.0,
) -> float:
    low = max([mu for j, mu in masses.items() if j <= 0] or [0.0])
    high_d = sum((2.0 ** j) * mu for j, mu in masses.items() if j >= 1)
    return coarse_increment_weight_upper(low, high_d, filter_radius, band_support_factor)


def packet_mass_entropy_route(
    packet_masses: Sequence[float],
    ancestry_labels: Sequence[object] | None = None,
    dominant_fraction: float = 0.25,
    ancestry_alpha: float = 0.5,
) -> dict[str, float | str]:
    """Route a large aggregate band mass into one atom or entropy/cycle mass.

    If no atom has fraction >=theta, Q_at<=theta and H_at>=-log(theta).
    With ancestry labels, either H_anc>=alpha h or hidden same-ancestry pair mass
    is at least theta^alpha-theta.
    """
    w = np.asarray(packet_masses, float)
    if np.any(w < 0) or float(w.sum()) <= 0:
        raise ValueError("positive aggregate packet mass required")
    if not (0 < dominant_fraction < 1) or not (0 < ancestry_alpha < 1):
        raise ValueError("invalid entropy thresholds")
    total = float(w.sum())
    probs = w / total
    imax = int(np.argmax(probs))
    pmax = float(probs[imax])
    h0 = -math.log(dominant_fraction)
    if pmax >= dominant_fraction:
        return {
            "branch": "dominant_packet",
            "aggregate_mass": total,
            "dominant_mass": float(w[imax]),
            "dominant_fraction": pmax,
            "entropy_threshold": h0,
        }
    q_at = float(np.dot(probs, probs))
    h_at = -math.log(q_at)
    if h_at + 1e-14 < h0:
        raise AssertionError("no-dominant-atom collision entropy bound failed")
    if ancestry_labels is None:
        return {
            "branch": "atomic_collision_entropy",
            "aggregate_mass": total,
            "H_atomic": h_at,
            "entropy_threshold": h0,
        }
    if len(ancestry_labels) != len(w):
        raise ValueError("ancestry label length mismatch")
    chain = collision_chain(probs, ancestry_labels)
    if chain["h_ancestry"] >= ancestry_alpha * h0 - 1e-14:
        return {
            "branch": "ancestry_Bellman_entropy",
            "aggregate_mass": total,
            "H_atomic": h_at,
            "H_ancestry": chain["h_ancestry"],
            "ancestry_entropy_lower": ancestry_alpha * h0,
        }
    pair_lower = dominant_fraction ** ancestry_alpha - dominant_fraction
    if chain["hidden_pair_mass"] + 2e-14 < pair_lower:
        raise AssertionError("same-ancestry pair lower bound failed")
    return {
        "branch": "same_ancestry_pair_cycle",
        "aggregate_mass": total,
        "H_atomic": h_at,
        "H_ancestry": chain["h_ancestry"],
        "hidden_pair_mass": chain["hidden_pair_mass"],
        "hidden_pair_lower": pair_lower,
    }


def persistent_dissipation_lower(enstrophy_threshold: float, scaled_time_measure: float) -> float:
    """Normalized dissipation D=N int||grad P_>N u||^2 dt on a persistent event set.

    Scaled time is tau=N^2 t, so D=int d_>(tau) d tau.
    """
    if enstrophy_threshold < 0 or scaled_time_measure < 0:
        raise ValueError("nonnegative inputs required")
    return enstrophy_threshold * scaled_time_measure


@dataclass(frozen=True)
class IncrementCollisionStress:
    samples: int
    worst_exact_over_coarse_weight: float
    minimum_collision_margin: float
    minimum_entropy_margin: float
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260808) -> IncrementCollisionStress:
    rng = np.random.default_rng(seed)
    wr = 0.0
    mc = me = float("inf")
    branches: dict[str, int] = {}
    for _ in range(samples):
        masses = {j: float(rng.lognormal(mean=-2.0, sigma=1.0)) for j in range(-8, 9)}
        Rg = float(rng.uniform(0.4, 1.8))
        beta = float(rng.uniform(1.0, 1.5))
        exact = exact_dyadic_square_weight(masses, Rg, beta)
        coarse = coarse_dyadic_bound_from_masses(masses, Rg, beta)
        if exact > coarse + 2e-12 * max(1.0, coarse):
            raise AssertionError("dyadic low/high coarse bound failed")
        if coarse > 1e-16:
            wr = max(wr, exact / coarse)

        g1 = float(rng.uniform(1.0, 2.0))
        clp = float(rng.uniform(1.0, 3.0))
        cb = float(rng.uniform(1.0, 2.0))
        Q = 0.97 * g1 * (clp * cb) ** 3 * exact ** 1.5
        if Q <= 1e-18:
            continue
        th = increment_collision_thresholds(Q, g1, clp, cb, Rg, beta)
        low = max(masses[j] for j in masses if j <= 0)
        dh = sum((2.0 ** j) * masses[j] for j in masses if j >= 1)
        margin = max(low - th["low_band_critical_mass"], dh - th["high_normalized_enstrophy"])
        mc = min(mc, margin)
        if margin < -2e-12:
            raise AssertionError("increment mass/enstrophy dichotomy failed")

        n = int(rng.integers(4, 30))
        pm = rng.dirichlet(np.full(n, 1.2))
        labels = rng.integers(0, max(2, n // 4), size=n).tolist()
        route = packet_mass_entropy_route(pm, labels, 0.25, 0.5)
        branches[route["branch"]] = branches.get(route["branch"], 0) + 1
        if route["branch"] == "ancestry_Bellman_entropy":
            me = min(me, float(route["H_ancestry"]) - float(route["ancestry_entropy_lower"]))
        elif route["branch"] == "same_ancestry_pair_cycle":
            me = min(me, float(route["hidden_pair_mass"]) - float(route["hidden_pair_lower"]))
        else:
            me = min(me, 0.0)
    return IncrementCollisionStress(samples, wr, mc, me, branches)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-onsager-increment-collision"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    data = {
        "status": "EXACT_SEQUENCE_ALGEBRA_GIVEN_STANDARD_LP_BERNSTEIN",
        "stress": out.__dict__,
        "theorem": {
            "increment_upper": "Q_N <= g1(C_LP C_B)^3[(4/3)(beta R_G)^2 mu_low_max+2 d_high]^(3/2)",
            "mass_threshold": "mu_low_max >= 3X/[8(beta R_G)^2]",
            "enstrophy_threshold": "d_high >= X/4",
            "X": "[Q_N/(g1(C_LP C_B)^3)]^(2/3)",
            "atom_entropy": "no atom >=theta aggregate mass => H_atomic>=-log theta",
            "ancestry_pair": "if H_anc<alpha(-log theta), hidden pair mass >=theta^alpha-theta",
        },
    }
    (args.outdir / "onsager_increment_collision.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Onsager cubic-increment to grain/enstrophy collision

Status: **EXACT_SEQUENCE_ALGEBRA_GIVEN_STANDARD_LP_BERNSTEIN**.

For a dyadic Littlewood--Paley decomposition at `N_j=2^j N`, standard `L^3` square-function and Bernstein estimates give

`Q_N <= g1(C_LP C_B)^3[(4/3)(beta R_G)^2 mu_low_max + 2 d_high]^(3/2)`,

where `mu_j=N_j||u_j||_2^2` and `d_high=sum_(j>=1)2^j mu_j = N^-1||grad P_>N u||_2^2` up to the fixed LP partition constants.

Thus if `X=[Q_N/(g1(C_LP C_B)^3)]^(2/3)`, every cubic increment event has the exact alternative

`mu_low_max >= 3X/[8(beta R_G)^2]`

or

`d_high >= X/4`.

A large aggregate band mass is not allowed to hide in many tiny packets.  For packet fractions `w_a`, either one atom has `w_a>=theta`, or `H_atomic>=-log(theta)`.  With ancestry labels the latter routes, by the exact existing collision chain rule, into component Bellman entropy or same-ancestry pair/cycle mass `>=theta^alpha-theta`.

If the enstrophy branch persists on a scaled-time set of measure `m`, normalized dissipation pays at least `m d_high`; failure of persistence is explicitly a temporal-concentration/CKN branch rather than a free escape.

Stress: `{out.samples}`
- worst exact dyadic square weight / coarse bound: `{out.worst_exact_over_coarse_weight:.9f}`
- minimum mass/enstrophy routing margin: `{out.minimum_collision_margin:.3e}`
- entropy branch counts: `{out.branch_counts}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
