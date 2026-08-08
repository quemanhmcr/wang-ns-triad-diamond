from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.ancestor_reservoir_sync import physical_energy_service_ratio_upper
from src.nested_grains import TriadEdge, incidence_components


def old_pool_service_capacity_upper(
    generation: int,
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    frame_energy_bound: float,
    global_energy: float,
    beta_filter_radius: float = 1.0,
) -> float:
    """Maximum increment-service capacity of an old materially transported pool.

    Assumptions at every generation q:
      * every old atom began with M_a,0 <= alpha N_0;
      * low-strain Kelvin transport gives M_a,q <= (21/20)^q M_a,0;
      * signed-good lineage gives N_q >= (8/5)^q N_0;
      * packet/Bessel budget sum_a E_a,q <= P E_global.

    Since each atom contributes beta^2 M^3 E/N^2, the whole old pool has the
    same geometric erosion ratio as one reservoir, independent of the number of atoms.
    """
    if generation < 0:
        raise ValueError("generation must be nonnegative")
    if min(initial_low_cut_ratio, initial_block_frequency, frame_energy_bound, beta_filter_radius) <= 0:
        raise ValueError("positive pool parameters required")
    if global_energy < 0:
        raise ValueError("nonnegative global energy required")
    base = (
        beta_filter_radius**2
        * initial_low_cut_ratio**3
        * initial_block_frequency
        * frame_energy_bound
        * global_energy
    )
    return base * float(physical_energy_service_ratio_upper()) ** generation


def total_old_pool_service_upper(
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    frame_energy_bound: float,
    global_energy: float,
    beta_filter_radius: float = 1.0,
) -> float:
    base = old_pool_service_capacity_upper(
        0,
        initial_low_cut_ratio,
        initial_block_frequency,
        frame_energy_bound,
        global_energy,
        beta_filter_radius,
    )
    r = float(physical_energy_service_ratio_upper())
    return base / (1.0 - r)


def first_forced_relink_generation(
    uniform_service_threshold: float,
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    frame_energy_bound: float,
    global_energy: float,
    beta_filter_radius: float = 1.0,
    new_fraction: float = 0.5,
) -> int:
    """First q where old pool cannot supply (1-new_fraction) of required service."""
    if uniform_service_threshold <= 0 or not (0.0 < new_fraction < 1.0):
        raise ValueError("invalid relink threshold")
    target = (1.0 - new_fraction) * uniform_service_threshold
    base = old_pool_service_capacity_upper(
        0,
        initial_low_cut_ratio,
        initial_block_frequency,
        frame_energy_bound,
        global_energy,
        beta_filter_radius,
    )
    if base < target:
        return 0
    r = float(physical_energy_service_ratio_upper())
    # Need the first integer q with base*r^q < target.
    q = int(math.floor(math.log(target / base) / math.log(r))) + 1
    return max(0, q)


def forced_new_service_lower(required_service: float, old_pool_capacity: float) -> float:
    if min(required_service, old_pool_capacity) < 0:
        raise ValueError("nonnegative service data required")
    return max(0.0, required_service - old_pool_capacity)


def relinking_incidence_route(edges: list[TriadEdge]) -> dict[str, float | int | str]:
    """Apply the existing exact 3-uniform fresh-or-cycle Euler identity.

    This is a finite-atomic corollary: relinking events must already have been
    represented by active quadratic triad edges after the existing tiny-edge pruning.
    """
    rows = incidence_components(edges)
    if not rows:
        return {"triads": 0, "fresh_units": 0, "cycle_rank": 0, "regime": "none"}
    if len(rows) != 1:
        raise ValueError("relinking route expects one connected incidence component")
    row = rows[0]
    m = int(row["triads"])
    fresh = int(row["fresh_units"])
    cyc = int(row["cycle_rank"])
    if max(fresh, cyc) < m:
        raise AssertionError("fresh-or-cycle relinking dichotomy failed")
    return row


def exact_pool_certificate() -> dict[str, str]:
    r = physical_energy_service_ratio_upper()
    if not (r < 1 / 2):
        raise AssertionError("old-pool erosion ratio is not below one half")
    return {
        "old_pool_ratio": f"{r.numerator}/{r.denominator}",
        "clean_old_pool_ratio": "<1/2",
        "frame_input": "sum packet energies <= P ||u||_2^2",
        "relink_graph_identity": "(n-1)+beta=2m; hence fresh_units>=m or cycle_rank>=m",
        "status": "EXACT_GIVEN_BESSEL_FRAME_BUDGET_AND_EXISTING_INCIDENCE_THEOREM",
    }


@dataclass(frozen=True)
class PoolErosionStress:
    samples: int
    maximum_generation_ratio: float
    minimum_half_life_margin: float
    minimum_total_budget_margin: float
    minimum_forced_new_service_margin: float
    incidence_checks: int


def stress(samples: int = 50_000, seed: int = 20260808) -> PoolErosionStress:
    rng = np.random.default_rng(seed)
    mr = 0.0
    mh = mt = mn = float("inf")
    r = float(physical_energy_service_ratio_upper())
    incidence_checks = 0
    for _ in range(samples):
        alpha = float(rng.uniform(0.05, 1.0))
        N0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        P = float(rng.uniform(1.0, 5.0))
        E = float(math.exp(rng.uniform(-4.0, 2.0)))
        beta = float(rng.uniform(0.2, 2.0))
        q = int(rng.integers(0, 20))
        c0 = old_pool_service_capacity_upper(0, alpha, N0, P, E, beta)
        cq = old_pool_service_capacity_upper(q, alpha, N0, P, E, beta)
        ratio = cq / c0 if c0 else 0.0
        expected = r**q
        mr = max(mr, ratio ** (1.0 / q) if q > 0 and ratio > 0 else 0.0)
        mh = min(mh, (0.5**q) * c0 - cq)
        if cq > c0 * 0.5**q + 2e-12 * max(1.0, c0):
            raise AssertionError("old pool clean half-life failed")
        total = total_old_pool_service_upper(alpha, N0, P, E, beta)
        mt = min(mt, 2.0 * c0 - total)
        if total >= 2.0 * c0 + 2e-12 * max(1.0, c0):
            raise AssertionError("old pool geometric budget failed")

        req = float(math.exp(rng.uniform(-5.0, 1.0)))
        old = float(rng.uniform(0.0, req))
        new = forced_new_service_lower(req, old)
        mn = min(mn, new - (req - old))
        if abs(new - (req - old)) > 1e-12 * max(1.0, req):
            raise AssertionError("forced new-service identity failed")

        if incidence_checks < 5000:
            # Connected chain: exact existing Euler accounting must route it.
            m = int(rng.integers(1, 12))
            edges: list[TriadEdge] = []
            prev = "v0"
            next_id = 1
            for j in range(m):
                a = prev
                b = f"v{next_id}"; next_id += 1
                c = f"v{next_id}"; next_id += 1
                edges.append(TriadEdge((a, b, c), 0.0, 1.0, 1.0))
                prev = c
            row = relinking_incidence_route(edges)
            if max(int(row["fresh_units"]), int(row["cycle_rank"])) < m:
                raise AssertionError("incidence relinking route failed")
            incidence_checks += 1
    return PoolErosionStress(samples, mr, mh, mt, mn, incidence_checks)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-reservoir-pool-erosion"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = exact_pool_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "reservoir_pool_erosion.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Reservoir pool erosion: old capacity forces relinking

Status: **{cert['status']}**.

Suppose a transfer-adapted old reservoir pool at generation zero has band frequencies `M_a,0<=alpha N_0` and a Bessel/frame energy budget

`sum_a E_a,q <= P E_global`

at every later service time.  On a signed-good low-strain lineage, every materially transported old atom grows by at most `21/20` per generation while the block scale grows by at least `8/5`.  Summing the low-band increment service `M^3 E/N^2` over the entire old pool gives

`C_old(q) <= beta^2 alpha^3 N_0 P E_global (231525/512000)^q < 2^-q C_old(0)`.

Thus the **whole old pool**, not just one chosen atom, has a geometric service half-life and total future capacity `<2 C_old(0)`.  If every efficient generation requires a uniform service threshold, after finitely many generations a quantitative fraction must come from newly relinked/fresh spectral atoms.

In the finite quadratic atomic model, once those relinking events are represented by active triad edges, the existing exact incidence identity

`(n-1)+beta=2m`

forces either at least `m` fresh units or cycle rank at least `m` in every connected component with `m` relinking triads.  Relabeling old capacity is therefore not a neutral third regime.

Stress: `{out.samples}`
- maximum effective one-step pool ratio: `{out.maximum_generation_ratio:.9f}`
- minimum half-life envelope margin: `{out.minimum_half_life_margin:.3e}`
- minimum total-budget margin: `{out.minimum_total_budget_margin:.3e}`
- minimum forced-new-service margin: `{out.minimum_forced_new_service_margin:.3e}`
- incidence checks: `{out.incidence_checks}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
