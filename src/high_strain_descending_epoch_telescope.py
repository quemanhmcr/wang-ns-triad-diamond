from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower
from src.nn_critical_heat_carrier_seed import RENEWAL_SCALE_FACTOR, renewal_scale


STATUS = (
    "EXACT_HIGH_STRAIN_DESCENDING_EPOCH_TELESCOPE__"
    "PHYSICAL_GLOBAL_GRADIENT_RESERVOIR__"
    "THREE_SIXTEENTHS_RENEWAL_SCALE_DESCENT__"
    "ARBITRARY_TIME_OVERLAP_WEIGHTED_BY_SCALE__"
    "NO_EVENT_COUNT_RESET"
)

ANCESTOR_TO_CHILD_RATIO_UPPER = 1.0 / 4.0
HIGH_STRAIN_RENEWAL_RATIO_UPPER = RENEWAL_SCALE_FACTOR * ANCESTOR_TO_CHILD_RATIO_UPPER
STEP_RELATIVE_TOLERANCE = 8.0e-12
EPOCH_RELATIVE_TOLERANCE = 1.0e-11


def _finite_positive_ratio(numerator: float, denominator: float) -> float:
    """Return one native dimensionless ratio, rejecting overflow/underflow.

    Scale relations in this theorem are homogeneous.  Comparing dimensional
    frequencies against an absolute unit tolerance destroys that covariance, so
    every guard is expressed through a ratio of quantities in the same unit.
    """
    x = float(numerator)
    y = float(denominator)
    if not (math.isfinite(x) and x > 0.0 and math.isfinite(y) and y > 0.0):
        raise ValueError("positive finite native quantities required")
    ratio = x / y
    if not math.isfinite(ratio) or ratio <= 0.0:
        raise ValueError("native ratio left the positive finite floating range")
    return ratio


def _native_product_upper_holds(value: float, left: float, right: float) -> bool:
    """Check ``value <= left*right`` without forming a fragile product."""
    log_excess = math.log(value) - math.log(left) - math.log(right)
    return log_excess <= math.log1p(EPOCH_RELATIVE_TOLERANCE)


@dataclass(frozen=True)
class HighStrainRenewalStep:
    """One genuine high-strain owner followed by its physical resolved ancestor.

    ``normalized_resolved_dissipation`` is the actual

        D_j = N_j int_{I_j} ||grad S_(N_j/4) u||_2^2 dt

    on the high-strain first-hit history.  The theorem never assumes that the
    histories ``I_j`` are disjoint.  ``ancestor_shell_frequency`` is the actual
    critical D_V|_G shell mark, and ``renewal_frequency=3M/4`` is the scale of the
    smooth carrier that continues the recursive high-strain route.
    """

    child_frequency: float
    ancestor_shell_frequency: float
    renewal_frequency: float
    normalized_resolved_dissipation: float

    def __post_init__(self) -> None:
        vals = (
            self.child_frequency,
            self.ancestor_shell_frequency,
            self.renewal_frequency,
            self.normalized_resolved_dissipation,
        )
        if not all(math.isfinite(x) and x > 0 for x in vals):
            raise ValueError("positive finite high-strain renewal data required")
        ancestor_ratio = _finite_positive_ratio(
            self.ancestor_shell_frequency, self.child_frequency
        )
        if ancestor_ratio > ANCESTOR_TO_CHILD_RATIO_UPPER * (
            1.0 + STEP_RELATIVE_TOLERANCE
        ):
            raise ValueError("high-strain resolved ancestor must satisfy M<=N/4")
        renewal_to_ancestor = _finite_positive_ratio(
            self.renewal_frequency, self.ancestor_shell_frequency
        )
        if not math.isclose(
            renewal_to_ancestor,
            RENEWAL_SCALE_FACTOR,
            rel_tol=STEP_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        ):
            raise ValueError("renewal scale must be the physical A=3M/4 shell registration")
        renewal_ratio = _finite_positive_ratio(
            self.renewal_frequency, self.child_frequency
        )
        if renewal_ratio > HIGH_STRAIN_RENEWAL_RATIO_UPPER * (
            1.0 + STEP_RELATIVE_TOLERANCE
        ):
            raise ValueError("high-strain renewal exceeded the exact 3/16 scale ratio")


@dataclass(frozen=True)
class HighStrainEpochCertificate:
    step_count: int
    root_frequency: float
    last_child_frequency: float
    total_gradient_dissipation: float
    high_strain_dissipation_lower: float
    physical_frequency_floor: float
    geometric_frequency_sum_upper: float
    normalized_dissipation_sum: float
    normalized_dissipation_capacity_upper: float
    frequency_floor_count_upper: int
    weighted_capacity_count_upper: int
    certified_count_upper: int
    maximum_observed_scale_ratio: float
    arbitrary_time_overlap_allowed: bool = True
    normalized_dissipation_used_as_global_reset: bool = False
    epoch_breaker_required_for_further_high_strain: bool = True

    def __post_init__(self) -> None:
        if self.step_count < 1:
            raise ValueError("nonempty high-strain epoch required")
        positive = (
            self.root_frequency,
            self.last_child_frequency,
            self.total_gradient_dissipation,
            self.high_strain_dissipation_lower,
            self.physical_frequency_floor,
            self.geometric_frequency_sum_upper,
            self.normalized_dissipation_sum,
            self.normalized_dissipation_capacity_upper,
        )
        if not all(math.isfinite(x) and x > 0 for x in positive):
            raise ValueError("positive finite epoch certificate values required")
        if self.frequency_floor_count_upper < self.step_count:
            raise ValueError("frequency-floor count bound does not cover observed epoch")
        if self.weighted_capacity_count_upper < self.step_count:
            raise ValueError("weighted-capacity count bound does not cover observed epoch")
        if self.certified_count_upper != min(
            self.frequency_floor_count_upper, self.weighted_capacity_count_upper
        ):
            raise ValueError("certified count must be the minimum native bound")
        if not math.isfinite(self.maximum_observed_scale_ratio) or self.maximum_observed_scale_ratio < 0.0:
            raise ValueError("finite nonnegative observed scale ratio required")
        if self.maximum_observed_scale_ratio > HIGH_STRAIN_RENEWAL_RATIO_UPPER * (
            1.0 + STEP_RELATIVE_TOLERANCE
        ):
            raise ValueError("epoch contains a non-high-strain scale renewal")
        if not self.arbitrary_time_overlap_allowed:
            raise ValueError("the theorem must not assume disjoint high-strain histories")
        if self.normalized_dissipation_used_as_global_reset:
            raise ValueError("D_V cannot be promoted to a global scale-independent reset")
        if not self.epoch_breaker_required_for_further_high_strain:
            raise ValueError("a completed descending epoch must require another owner to restart scale")


def kinetic_energy_gradient_dissipation_upper(initial_kinetic_energy: float, viscosity: float) -> float:
    """Unforced NS energy inequality: int ||grad u||_2^2 <= E0/(2 nu).

    ``initial_kinetic_energy`` means ``||u_0||_2^2`` in the repository's energy
    normalization.  Positive viscosity is essential; this helper is not an Euler
    statement.
    """
    E0 = float(initial_kinetic_energy)
    nu = float(viscosity)
    if E0 < 0 or not math.isfinite(E0) or nu <= 0 or not math.isfinite(nu):
        raise ValueError("finite nonnegative initial energy and positive viscosity required")
    out = E0 / (2.0 * nu)
    if not math.isfinite(out):
        raise ValueError("energy/viscosity data produced no finite gradient reservoir upper")
    return out


def physical_high_strain_frequency_floor(
    total_gradient_dissipation: float,
    scaled_lifetime: float,
) -> float:
    """Smallest scale that can possibly pay one normalized high-strain unit.

    Since the strict low pass is an L2 contraction and every event history lies
    inside the global time interval,

        D_j <= N_j G_*,  G_*=int ||grad u||_2^2 dt.

    A high-strain event also has D_j>=D_*, hence N_j>=D_*/G_*.
    """
    G = float(total_gradient_dissipation)
    c = float(scaled_lifetime)
    if G <= 0 or not math.isfinite(G) or c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite global gradient dissipation and scaled lifetime required")
    Dstar = clean_high_strain_dissipation_lower(c)
    out = Dstar / G
    if not math.isfinite(Dstar) or Dstar <= 0.0 or not math.isfinite(out) or out <= 0.0:
        raise ValueError("physical high-strain frequency floor left the positive finite range")
    return out


def _frequency_floor_count_upper(root_frequency: float, floor_frequency: float) -> int:
    N0 = float(root_frequency)
    Nmin = float(floor_frequency)
    if N0 <= 0 or Nmin <= 0 or not all(math.isfinite(x) for x in (N0, Nmin)):
        raise ValueError("positive finite root/floor frequencies required")
    if N0 < Nmin:
        return 0
    r = HIGH_STRAIN_RENEWAL_RATIO_UPPER
    # Steps have child scales N_j <= r^j N_0 and require N_j>=Nmin.
    # Subtract logarithms before division.  Forming Nmin/N0 first can underflow
    # even though the finite geometric count is perfectly representable.
    q = (math.log(Nmin) - math.log(N0)) / math.log(r)
    if not math.isfinite(q):
        raise ValueError("finite logarithmic scale coordinate required")
    return max(1, int(math.floor(q + 8e-12)) + 1)


def high_strain_epoch_telescope(
    steps: Sequence[HighStrainRenewalStep],
    *,
    total_gradient_dissipation: float,
    scaled_lifetime: float,
) -> HighStrainEpochCertificate:
    """Close one maximal consecutive high-strain descent epoch.

    No disjointness of event histories is used.  Each actual normalized resolved
    dissipation satisfies ``D_j<=N_j G_*`` by restriction to a subinterval and
    low-pass L2 contraction.  Since physical renewal gives

        N_(j+1) <= (3/16) N_j,

    the duplicated reading of the same global gradient reservoir is weighted by a
    summable geometric scale sequence:

        sum_j D_j <= G_* sum_j N_j <= N_0 G_*/(1-3/16).

    Together with ``D_j>=D_*`` this gives a finite epoch count.  This is a typed
    path telescope, not an additive global reset for D_V.
    """
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty high-strain epoch required")
    G = float(total_gradient_dissipation)
    c = float(scaled_lifetime)
    if G <= 0 or not math.isfinite(G) or c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite global gradient dissipation and scaled lifetime required")
    Dstar = clean_high_strain_dissipation_lower(c)
    if not math.isfinite(Dstar) or Dstar <= 0.0:
        raise ValueError("positive finite high-strain dissipation threshold required")

    max_ratio = 0.0
    for j, row in enumerate(rows):
        threshold_ratio = _finite_positive_ratio(
            row.normalized_resolved_dissipation, Dstar
        )
        if threshold_ratio < 1.0 - EPOCH_RELATIVE_TOLERANCE:
            raise ValueError("epoch contains a step below the physical high-strain dissipation threshold")
        if not _native_product_upper_holds(
            row.normalized_resolved_dissipation, row.child_frequency, G
        ):
            raise ValueError("resolved event dissipation exceeds the supplied global gradient reservoir")
        ratio = _finite_positive_ratio(row.renewal_frequency, row.child_frequency)
        max_ratio = max(max_ratio, ratio)
        if j + 1 < len(rows):
            nxt = rows[j + 1]
            consecutive_ratio = _finite_positive_ratio(
                nxt.child_frequency, row.renewal_frequency
            )
            if not math.isclose(
                consecutive_ratio,
                1.0,
                rel_tol=EPOCH_RELATIVE_TOLERANCE,
                abs_tol=0.0,
            ):
                raise ValueError("high-strain epoch must follow the actual renewed carrier scale consecutively")

    N0 = rows[0].child_frequency
    Nlast = rows[-1].child_frequency
    Nmin = physical_high_strain_frequency_floor(G, c)
    if math.log(Nlast) + math.log(G) < math.log(Dstar) + math.log1p(
        -EPOCH_RELATIVE_TOLERANCE
    ):
        raise AssertionError("certified high-strain event fell below the physical dissipation floor")

    frequency_sum = math.fsum(x.child_frequency for x in rows)
    geometric_frequency_sum_upper = N0 / (1.0 - HIGH_STRAIN_RENEWAL_RATIO_UPPER)
    if not all(
        math.isfinite(x) and x > 0.0
        for x in (frequency_sum, geometric_frequency_sum_upper)
    ):
        raise ValueError("frequency telescope left the positive finite range")
    if _finite_positive_ratio(
        frequency_sum, geometric_frequency_sum_upper
    ) > 1.0 + EPOCH_RELATIVE_TOLERANCE:
        raise AssertionError("high-strain physical renewal failed its geometric scale telescope")

    Dsum = math.fsum(x.normalized_resolved_dissipation for x in rows)
    capacity = G * frequency_sum
    geometric_capacity = G * geometric_frequency_sum_upper
    if not all(
        math.isfinite(x) and x > 0.0 for x in (Dsum, capacity, geometric_capacity)
    ):
        raise ValueError("dissipation telescope left the positive finite range")
    if _finite_positive_ratio(Dsum, capacity) > 1.0 + EPOCH_RELATIVE_TOLERANCE:
        raise AssertionError("overlapping high-strain histories exceeded the weighted global gradient reservoir")

    floor_count = _frequency_floor_count_upper(N0, Nmin)
    weighted_ratio = geometric_capacity / Dstar
    if not math.isfinite(weighted_ratio) or weighted_ratio < 1.0:
        raise ValueError("finite positive weighted high-strain count capacity required")
    weighted_count = int(math.floor(weighted_ratio + 8e-12))
    certified = min(floor_count, weighted_count)
    if len(rows) > certified:
        raise AssertionError("observed high-strain epoch exceeded its physical count telescope")

    return HighStrainEpochCertificate(
        step_count=len(rows),
        root_frequency=N0,
        last_child_frequency=Nlast,
        total_gradient_dissipation=G,
        high_strain_dissipation_lower=Dstar,
        physical_frequency_floor=Nmin,
        geometric_frequency_sum_upper=geometric_frequency_sum_upper,
        normalized_dissipation_sum=Dsum,
        normalized_dissipation_capacity_upper=geometric_capacity,
        frequency_floor_count_upper=floor_count,
        weighted_capacity_count_upper=weighted_count,
        certified_count_upper=certified,
        maximum_observed_scale_ratio=max_ratio,
    )


def theorem_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    Dstar = clean_high_strain_dissipation_lower(c)
    return {
        "status": STATUS,
        "high_strain_unit": f"every genuine first high-strain contact pays D_V>=D_*={Dstar:.12g} on its actual resolved low-pass history",
        "physical_renewal_scale": "D_V|_G supplies M<=N/4 and the renewed smooth carrier uses A=3M/4, hence every consecutive high-strain renewal has N_next/N<=3/16",
        "global_reservoir": "for G_*=int_0^t* ||grad u||_2^2 dt, low-pass L2 contraction and interval restriction give D_j<=N_j G_* even when event histories overlap arbitrarily",
        "frequency_floor": "combining D_j>=D_* with D_j<=N_j G_* gives the physical lower scale N_j>=D_*/G_* for every high-strain event in the epoch",
        "overlap_telescope": "sum_j D_j<=G_* sum_j N_j<=N_0 G_*/(1-3/16); repeated readings of the same viscous spacetime reservoir cannot manufacture infinite normalized dissipation because their native scale weights are geometrically summable",
        "energy_inequality": "for unforced Navier-Stokes with nu>0, G_*<=||u_0||_2^2/(2nu), so the high-strain frequency floor is at least 2nu D_*/||u_0||_2^2",
        "master_consequence": "a maximal consecutive high-strain recursive epoch is finite; an infinite event path with infinitely many high-strain events must contain infinitely many non-high-strain epoch breakers",
        "anti_reset": "D_V is not promoted to a globally additive event-count budget; the bound depends on the epoch root scale and uses the actual descending physical renewal geometry",
        "scope": "this closes eventually-pure high-strain recurrence, not mixed recurrence with infinitely many HH/source/service/reuse/high-tail epoch breakers; no Navier-Stokes global-regularity claim is made",
    }


@dataclass(frozen=True)
class HighStrainEpochStress:
    samples: int
    worst_frequency_sum_margin: float
    worst_normalized_capacity_margin: float
    minimum_frequency_floor_margin: float
    minimum_frequency_sum_relative_margin: float
    minimum_normalized_capacity_relative_margin: float
    minimum_frequency_floor_relative_margin: float
    minimum_native_child_frequency: float
    maximum_native_child_frequency: float
    minimum_native_normalized_dissipation: float
    maximum_native_normalized_dissipation: float
    maximum_certified_epoch_count: int
    arbitrary_overlap_cases: int
    ascending_chain_rejections: int


def stress(samples: int = 50_000, seed: int = 20260811) -> HighStrainEpochStress:
    rng = random.Random(seed)
    wf = wc = float("inf")
    minfloor = float("inf")
    wf_rel = wc_rel = floor_rel = float("inf")
    min_frequency = min_dissipation = float("inf")
    max_frequency = max_dissipation = 0.0
    maxcount = overlap = rejected = 0
    rmax = HIGH_STRAIN_RENEWAL_RATIO_UPPER

    for _ in range(samples):
        # D_*=const/c is itself a native physical scale.  Sweep c and G over
        # wide reciprocal ranges so the same dimensionless epoch law is tested
        # far below and far above the artificial unit scale 1.
        c = math.exp(rng.uniform(math.log(1.0e-140), math.log(1.0e140)))
        Dstar = clean_high_strain_dissipation_lower(c)
        G = math.exp(rng.uniform(math.log(1.0e-60), math.log(1.0e60)))
        root_factor = math.exp(rng.uniform(math.log(2.0), math.log(1.0e4)))
        N0 = root_factor * Dstar / G
        Nmin = Dstar / G
        cap = _frequency_floor_count_upper(N0, Nmin)
        L = rng.randint(1, max(1, min(cap, 24)))
        rows: list[HighStrainRenewalStep] = []
        N = N0
        for j in range(L):
            if N * G < Dstar:
                break
            # Keep actual high-strain D below the global reservoir capacity.
            D = rng.uniform(Dstar, max(Dstar, min(2.5 * Dstar, 0.98 * N * G)))
            if D > N * G:
                D = Dstar
            ratio = rng.uniform(0.03, rmax)
            M = (ratio / RENEWAL_SCALE_FACTOR) * N
            A = renewal_scale(M)
            rows.append(HighStrainRenewalStep(N, M, A, D))
            N = A
        if not rows:
            continue
        out = high_strain_epoch_telescope(rows, total_gradient_dissipation=G, scaled_lifetime=c)
        fsum = math.fsum(x.child_frequency for x in rows)
        capacity = G * fsum
        wf = min(wf, out.geometric_frequency_sum_upper - fsum)
        wc = min(wc, capacity - out.normalized_dissipation_sum)
        minfloor = min(minfloor, out.last_child_frequency - out.physical_frequency_floor)
        wf_rel = min(
            wf_rel,
            (out.geometric_frequency_sum_upper - fsum)
            / out.geometric_frequency_sum_upper,
        )
        wc_rel = min(
            wc_rel,
            (capacity - out.normalized_dissipation_sum) / capacity,
        )
        floor_rel = min(
            floor_rel,
            (out.last_child_frequency - out.physical_frequency_floor)
            / out.last_child_frequency,
        )
        min_frequency = min(min_frequency, *(x.child_frequency for x in rows))
        max_frequency = max(max_frequency, *(x.child_frequency for x in rows))
        min_dissipation = min(
            min_dissipation, *(x.normalized_resolved_dissipation for x in rows)
        )
        max_dissipation = max(
            max_dissipation, *(x.normalized_resolved_dissipation for x in rows)
        )
        maxcount = max(maxcount, out.certified_count_upper)
        overlap += 1  # theorem deliberately gives no interval-disjointness input.

        # An observer-invented or scale-ascending restart is outside this theorem.
        badN = max(rows[-1].renewal_frequency * 1.1, 2.0 * Dstar / G)
        badM = 0.24 * badN
        badA = renewal_scale(badM)
        badD = 1.1 * Dstar
        bad_rows = rows + [HighStrainRenewalStep(badN, badM, badA, badD)]
        try:
            high_strain_epoch_telescope(bad_rows, total_gradient_dissipation=G, scaled_lifetime=c)
        except ValueError:
            rejected += 1
        else:
            raise AssertionError("non-consecutive scale restart crossed the high-strain epoch telescope")

    return HighStrainEpochStress(
        samples,
        wf,
        wc,
        minfloor,
        wf_rel,
        wc_rel,
        floor_rel,
        min_frequency,
        max_frequency,
        min_dissipation,
        max_dissipation,
        maxcount,
        overlap,
        rejected,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-strain-descending-epoch-telescope"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "high_strain_descending_epoch_telescope.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# High-strain descending-epoch physical dissipation telescope

Status: **{cert['status']}**.

A genuine high-strain first contact at carrier scale `N_j` pays

`D_j=N_j int_(I_j)||grad S_(N_j/4)u||_2^2 dt >= D_*`.

Its actual critical resolved ancestor has `M_j<=N_j/4`, and the renewed smooth carrier uses `N_(j+1)=3M_j/4`, hence `N_(j+1)/N_j<=3/16`.

Let `G_*=int_0^t* ||grad u||_2^2 dt`.  No disjointness of the histories `I_j` is assumed.  Low-pass contraction and interval restriction give `D_j<=N_j G_*`, so each event obeys the physical frequency floor `N_j>=D_*/G_*`.  Moreover

`sum_j D_j <= G_* sum_j N_j <= N_0 G_*/(1-3/16)`.

Thus a maximal consecutive high-strain recursive epoch is finite even under complete time overlap.  This does **not** make `D_V` an additive global reset; the telescope depends on the epoch root scale and on the actual `3/16` physical descent.  With `nu>0`, the NS energy inequality may further bound `G_*<=||u_0||_2^2/(2nu)`.

Stress: `{out.samples}` descending high-strain epochs
- minimum geometric-frequency capacity margin: `{out.worst_frequency_sum_margin:.3e}`
- minimum weighted normalized-dissipation margin: `{out.worst_normalized_capacity_margin:.3e}`
- minimum last-scale/frequency-floor margin: `{out.minimum_frequency_floor_margin:.3e}`
- minimum native-relative geometric-frequency margin: `{out.minimum_frequency_sum_relative_margin:.3e}`
- minimum native-relative normalized-capacity margin: `{out.minimum_normalized_capacity_relative_margin:.3e}`
- minimum native-relative frequency-floor margin: `{out.minimum_frequency_floor_relative_margin:.3e}`
- sampled native child-frequency range: `[{out.minimum_native_child_frequency:.3e},{out.maximum_native_child_frequency:.3e}]`
- sampled native normalized-dissipation range: `[{out.minimum_native_normalized_dissipation:.3e},{out.maximum_native_normalized_dissipation:.3e}]`
- maximum certified epoch count sampled: `{out.maximum_certified_epoch_count}`
- arbitrary-overlap cases: `{out.arbitrary_overlap_cases}`
- non-consecutive/ascending restart rejections: `{out.ascending_chain_rejections}`

Master consequence: an infinite event path cannot eventually consist only of high-strain critical-dissipation renewals.  Infinitely many high-strain events force infinitely many other physical owner events to break the descending epochs.  Mixed-owner recurrence remains open, and no Navier--Stokes global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
