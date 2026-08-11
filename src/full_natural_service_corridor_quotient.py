from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from src.critical_annular_carrier_service_reentry import transported_annular_support_ratios
from src.heat_edge_material_ownership import partition_positive_edge_measure


STATUS = (
    "EXACT_FULL_NATURAL_SERVICE_CORRIDOR_QUOTIENT__"
    "OWN_SCALE_SERVICE_IS_SAME_INTERVAL_WITNESS_NOT_NEW_EVENT__"
    "MATERIAL_DISINTEGRATION_ZERO_RECURSION_DEPTH__"
    "ENDPOINT_SMOOTH_CARRIER_HAS_COMPARABLE_HARD_SHELL_WITNESS_SET"
)

FULL_NATURAL_SERVICE_WITNESS = "full_natural_own_scale_service"
RENEWAL_TO_PARENT_SHELL_RATIO = 3.0 / 4.0
TWO_HARD_SHELL_COVER_FACTOR = 2.0 / 3.0
RELATIVE_CERTIFICATE_TOLERANCE = 6.0e-12


def _relative_certificate_close(left: float, right: float, *, tolerance: float = RELATIVE_CERTIFICATE_TOLERANCE) -> bool:
    """Compare like-dimensional positive data without an absolute unit floor.

    The corridor is used precisely on UV tails where ``c A^-2`` can be
    arbitrarily small.  A tolerance proportional to ``max(1, ...)`` would turn
    one unit of the observer's time/mass convention into certificate slack and
    would therefore destroy parabolic/unit covariance.
    """
    a = float(left)
    b = float(right)
    eps = float(tolerance)
    if not all(math.isfinite(x) for x in (a, b, eps)) or eps < 0:
        return False
    scale = max(abs(a), abs(b))
    if scale == 0.0:
        return a == b
    return abs(a - b) <= eps * scale


@dataclass(frozen=True)
class FullNaturalServiceCorridor:
    """One completed physical natural corridor with its attached service law.

    The service is physically real, but it lives on the interval which has already
    been traversed in physical time.  Reading or disintegrating that same law does
    not create another event time or another recursion edge.
    """

    terminal_time: float
    endpoint_time: float
    renewal_frequency: float
    scaled_lifetime: float
    uniform_service_lower: float
    integrated_service_lower: float
    endpoint_carrier_critical_mass_lower: float

    def __post_init__(self) -> None:
        vals = (
            self.terminal_time,
            self.endpoint_time,
            self.renewal_frequency,
            self.scaled_lifetime,
            self.uniform_service_lower,
            self.integrated_service_lower,
            self.endpoint_carrier_critical_mass_lower,
        )
        if not all(math.isfinite(x) for x in vals):
            raise ValueError("finite corridor/service data required")
        if self.terminal_time <= 0 or self.endpoint_time <= 0 or self.endpoint_time >= self.terminal_time:
            raise ValueError("full-natural service requires one nontrivial backward physical interval")
        if self.renewal_frequency <= 0 or self.scaled_lifetime <= 0:
            raise ValueError("positive renewal frequency and lifetime required")
        if min(
            self.uniform_service_lower,
            self.integrated_service_lower,
            self.endpoint_carrier_critical_mass_lower,
        ) <= 0:
            raise ValueError("positive own-scale service and endpoint carrier lower required")
        expected = self.scaled_lifetime / (self.renewal_frequency * self.renewal_frequency)
        actual = self.terminal_time - self.endpoint_time
        if not _relative_certificate_close(actual, expected):
            raise ValueError("service witness interval is not the completed natural corridor")
        expected_integrated = self.scaled_lifetime * self.uniform_service_lower
        if not _relative_certificate_close(self.integrated_service_lower, expected_integrated):
            raise ValueError("integrated service lower is not the normalized service of this corridor")

    @property
    def physical_time_drop(self) -> float:
        return self.terminal_time - self.endpoint_time

    @property
    def parent_shell_frequency(self) -> float:
        return self.renewal_frequency / RENEWAL_TO_PARENT_SHELL_RATIO


def quotient_full_natural_service_outcome(
    outcome: dict[str, object],
    *,
    event_time: float,
    renewal_frequency: float,
    scaled_lifetime: float,
) -> FullNaturalServiceCorridor:
    """Re-type the generic shell full-survivor output as a same-corridor witness.

    No physical service mass is discarded.  The quotient removes only a false
    *additional event depth*: the service lower was proved on the very corridor
    whose physical-time drop has already occurred.
    """
    if str(outcome.get("classification", "")) != FULL_NATURAL_SERVICE_WITNESS:
        raise ValueError("only a completed full-natural own-scale service outcome can be quotiented")
    stops = tuple(outcome.get("joint_first_stops", outcome.get("joint_causes", ())))
    if stops:
        raise ValueError("a full-natural service witness cannot simultaneously contain a first stop")
    if bool(outcome.get("requires_physical_energy_reentry", False)):
        raise ValueError("unresolved coefficient obstruction cannot masquerade as full-natural service")
    if bool(outcome.get("coefficient_impulses_used_as_work", False)):
        raise ValueError("coefficient impulse cannot be used as service work")
    if not bool(outcome.get("service_same_corridor_witness", False)):
        raise ValueError("service outcome has not certified same-corridor witness semantics")
    if bool(outcome.get("service_adds_recursion_depth", True)):
        raise ValueError("service outcome still claims an additional recursion edge")

    t = float(event_time)
    A = float(renewal_frequency)
    c = float(scaled_lifetime)
    if t <= 0 or A <= 0 or c <= 0 or not all(math.isfinite(x) for x in (t, A, c)):
        raise ValueError("positive finite event time, renewal frequency and lifetime required")
    required = c / (A * A)
    if required >= t:
        raise ValueError("a corridor reaching t=0 is an absorbing boundary, not a full-natural service witness")
    reported = float(outcome.get("required_elapsed", -1.0))
    if reported <= 0 or not _relative_certificate_close(reported, required):
        raise ValueError("service outcome does not report the same completed natural interval")
    horizon = float(outcome.get("observed_elapsed_end", -1.0))
    horizon_tol = RELATIVE_CERTIFICATE_TOLERANCE * max(abs(horizon), required)
    if horizon <= 0 or horizon + horizon_tol < required:
        raise ValueError("service outcome was inferred from an incomplete monitor horizon")

    provenance = {
        "corridor_terminal_time": t,
        "corridor_endpoint_time": t - required,
        "physical_time_drop": required,
        "renewal_frequency": A,
        "scaled_lifetime": c,
        "parent_shell_frequency": A / RENEWAL_TO_PARENT_SHELL_RATIO,
    }
    for field, expected in provenance.items():
        actual = float(outcome.get(field, math.nan))
        if not _relative_certificate_close(actual, expected):
            raise ValueError(f"service outcome {field} does not match the physical corridor provenance")

    return FullNaturalServiceCorridor(
        terminal_time=t,
        endpoint_time=t - required,
        renewal_frequency=A,
        scaled_lifetime=c,
        uniform_service_lower=float(outcome["uniform_square_service_lower"]),
        integrated_service_lower=float(outcome["integrated_bounded_heat_service_lower"]),
        endpoint_carrier_critical_mass_lower=float(outcome["endpoint_carrier_critical_mass_lower"]),
    )


def endpoint_hard_shell_cover_from_full_natural_outcome(
    outcome: dict[str, object],
    *,
    parent_shell_frequency: float,
) -> dict[str, object]:
    """Read the endpoint hard-shell cover directly from a certified survivor output."""
    if str(outcome.get("classification", "")) != FULL_NATURAL_SERVICE_WITNESS:
        raise ValueError("endpoint hard-shell cover requires a full-natural survivor output")
    if not bool(outcome.get("service_same_corridor_witness", False)) or bool(outcome.get("service_adds_recursion_depth", True)):
        raise ValueError("full-natural outcome has not certified same-corridor service semantics")
    A = float(outcome.get("renewal_frequency", -1.0))
    c = float(outcome.get("scaled_lifetime", -1.0))
    reported_parent = float(outcome.get("parent_shell_frequency", -1.0))
    if min(A, c, reported_parent) <= 0 or not all(math.isfinite(x) for x in (A, c, reported_parent)):
        raise ValueError("full-natural outcome supplied no physical scale provenance")
    derived_parent = A / RENEWAL_TO_PARENT_SHELL_RATIO
    requested_parent = float(parent_shell_frequency)
    if not _relative_certificate_close(reported_parent, derived_parent) or not _relative_certificate_close(
        requested_parent, derived_parent
    ):
        raise ValueError("parent shell frequency does not match the certified carrier scale provenance")
    mu = float(outcome.get("endpoint_carrier_critical_mass_lower", -1.0))
    if mu <= 0 or not math.isfinite(mu):
        raise ValueError("full-natural outcome supplied no positive endpoint carrier mass lower")
    return endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=derived_parent,
        endpoint_carrier_critical_mass=mu,
    )


def material_partition_is_same_corridor_measure(
    corridor: FullNaturalServiceCorridor,
    edge_weights: Sequence[float],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
) -> dict[str, object]:
    """OO/ON/NN disintegrates the attached positive service law, not physical time.

    The output is deliberately a witness record.  It creates no second causal
    charge, no second service mass, and no recursion edge.  A downstream theorem
    may use one of these positive submeasures to certify a new physical state, but
    that theorem must supply the actual state/time it creates.
    """
    part = partition_positive_edge_measure(edge_weights, old_here, old_neighbor)
    total = float(part["total"])
    lower = corridor.integrated_service_lower
    tol = RELATIVE_CERTIFICATE_TOLERANCE * max(total, lower)
    if total <= 0 or total + tol < lower:
        raise ValueError("edge measure does not realize the positive integrated service lower of the same corridor")
    return {
        **part,
        "corridor_terminal_time": corridor.terminal_time,
        "corridor_endpoint_time": corridor.endpoint_time,
        "physical_time_drop_already_counted": corridor.physical_time_drop,
        "recursion_edges_added": 0,
        "causal_charge_created": False,
        "service_mass_duplicated": False,
        "same_positive_measure": True,
        "integrated_service_lower": lower,
        "service_lower_margin": total - lower,
    }


def endpoint_comparable_hard_shell_cover(
    *,
    parent_shell_frequency: float,
    endpoint_carrier_critical_mass: float,
) -> dict[str, object]:
    """The full-survivor smooth carrier already contains an actual endpoint shell.

    The generic shell uses A=3M/4.  Its transported smooth carrier remains inside
    (A/2,2A), while |Q_A|<=1.  Splitting that annulus into the exact hard shells
    (A/2,A] and (A,2A] gives

        A ||Q_A u||_2^2 <= mu_A + mu_(2A)/2
                           <= (3/2) max(mu_A,mu_(2A)).

    Hence at the *same physical endpoint* at least one hard shell has critical
    mass >=(2/3) A||Q_Au||^2.  The theorem returns the witness set rather than
    imposing an observer-chosen tie break.
    """
    M = float(parent_shell_frequency)
    mu = float(endpoint_carrier_critical_mass)
    if M <= 0 or mu <= 0 or not all(math.isfinite(x) for x in (M, mu)):
        raise ValueError("positive finite parent shell frequency and carrier mass required")
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    lo, hi = transported_annular_support_ratios()
    if lo <= 0.5 or hi >= 2.0:
        raise AssertionError("transported smooth carrier escaped the two-hard-shell endpoint cover")
    return {
        "renewal_frequency": A,
        "parent_shell_frequency": M,
        "hard_shell_candidates": (A, 2.0 * A),
        "next_corridor_renewal_candidates": (
            RENEWAL_TO_PARENT_SHELL_RATIO * A,
            RENEWAL_TO_PARENT_SHELL_RATIO * 2.0 * A,
        ),
        "candidate_ratios_to_parent": (RENEWAL_TO_PARENT_SHELL_RATIO, 2.0 * RENEWAL_TO_PARENT_SHELL_RATIO),
        "guaranteed_max_hard_shell_critical_mass_lower": TWO_HARD_SHELL_COVER_FACTOR * mu,
        "same_physical_endpoint": True,
        "tie_rule": "retain_all_exact_maximizers_as_joint_witnesses; no causal priority is created",
        "new_causal_charge_created": False,
        "new_physical_time_edge_created": False,
    }


def realized_endpoint_hard_shell_witnesses(
    cover: dict[str, object],
    candidate_critical_masses: Sequence[float],
    *,
    tie_tolerance: float = 2e-12,
) -> dict[str, object]:
    """Read all maximizing endpoint shell witnesses without lexicographic priority."""
    freqs = tuple(float(x) for x in cover["hard_shell_candidates"])
    masses = tuple(float(x) for x in candidate_critical_masses)
    if len(freqs) != 2 or len(masses) != 2 or any((not math.isfinite(x) or x < 0) for x in masses):
        raise ValueError("two finite nonnegative endpoint hard-shell masses required")
    lower = float(cover["guaranteed_max_hard_shell_critical_mass_lower"])
    mx = max(masses)
    cert_tol = RELATIVE_CERTIFICATE_TOLERANCE * max(mx, lower)
    if mx + cert_tol < lower:
        raise ValueError("actual endpoint shell masses do not realize the certified smooth-carrier cover")
    tie_eps = float(tie_tolerance)
    if not math.isfinite(tie_eps) or tie_eps < 0:
        raise ValueError("finite nonnegative tie tolerance required")
    tie_tol = tie_eps * mx
    ids = tuple(i for i, x in enumerate(masses) if mx - x <= tie_tol)
    return {
        "joint_witness_frequencies": tuple(freqs[i] for i in ids),
        "joint_witness_critical_masses": tuple(masses[i] for i in ids),
        "joint_next_corridor_renewal_frequencies": tuple(
            RENEWAL_TO_PARENT_SHELL_RATIO * freqs[i] for i in ids
        ),
        "maximum_critical_mass": mx,
        "recursion_edges_added": 0,
        "physical_time_drop_added": 0.0,
        "causal_primary_selected": False,
    }


def theorem_certificate() -> dict[str, object]:
    lo, hi = transported_annular_support_ratios()
    return {
        "status": STATUS,
        "time_ontology": "full_natural_own_scale_service is a positive law carried by the already-completed physical interval [t-cA^-2,t]; reading it does not create a second event time or recursion edge",
        "scale_provenance": "the producer carries parent shell M, renewal carrier A=3M/4, scaled lifetime c, and both physical endpoints; downstream adapters may verify but never rebind those values",
        "unit_covariance": "all corridor, service, and hard-shell certificate comparisons use relative native-unit slack with no max(1,...) absolute floor",
        "material_ontology": "OO/ON/NN is a positive disintegration of that same service measure and adds zero causal charge and zero recursion depth",
        "material_mass_gate": "the complete integrated edge measure must realize at least the corridor's positive integrated service lower before it can be marked as the same physical law",
        "endpoint_carrier": "the surviving Q_A carrier is present at the corridor endpoint itself; service is a sidecar/witness, not a prerequisite for inventing endpoint persistence",
        "hard_shell_cover": f"transported support lies in ({lo:.12g}A,{hi:.12g}A) subset (A/2,2A); therefore max(mu_A,mu_2A)>=(2/3) A||Q_Au||^2 at the same endpoint",
        "scale_geometry": "with A=3M/4, the endpoint hard-shell witness ratios are 3/4 or 3/2 relative to M; the next corridor registrations are 3/4 of those witnessed shell scales. Parent shell, carrier, witness shell, and next carrier remain separately typed",
        "ties": "equal hard-shell witnesses remain joint; no theorem-name or frequency-order priority is introduced",
        "master_quotient": "a chain which only alternates full-natural corridors with their own-scale service/material witness layers is just a chain of full-natural physical corridors; bounded-scale such chains hit t=0 by physical time",
        "scope": "this removes service-theorem depth from the master and closes the endpoint-service attachment seam; it does not terminate UV-unbounded full-survivor chains or genuine first-hit/work/source/reuse owner recurrence",
    }


@dataclass(frozen=True)
class ServiceCorridorStress:
    samples: int
    worst_time_identity_relative_residual: float
    worst_material_partition_relative_residual: float
    minimum_hard_shell_cover_relative_margin: float
    joint_tie_witnesses: int


def stress(samples: int = 50_000, seed: int = 20260811) -> ServiceCorridorStress:
    rng = random.Random(seed)
    wt = wp = 0.0
    mh = math.inf
    ties = 0
    for i in range(samples):
        A = 10.0 ** rng.uniform(-40.0, 40.0)
        c = 10.0 ** rng.uniform(-12.0, 12.0)
        T = c / (A * A)
        t = T * (1.0 + math.exp(rng.uniform(-2.0, 2.0)))
        y = 10.0 ** rng.uniform(-120.0, 80.0)
        outcome = {
            "classification": FULL_NATURAL_SERVICE_WITNESS,
            "joint_first_stops": (),
            "required_elapsed": T,
            "observed_elapsed_end": T * (1.0 + rng.uniform(0.0, 0.25)),
            "uniform_square_service_lower": y,
            "integrated_bounded_heat_service_lower": c * y,
            "endpoint_carrier_critical_mass_lower": 2.0 * y,
            "corridor_terminal_time": t,
            "corridor_endpoint_time": t - T,
            "physical_time_drop": T,
            "renewal_frequency": A,
            "scaled_lifetime": c,
            "parent_shell_frequency": A / RENEWAL_TO_PARENT_SHELL_RATIO,
            "service_same_corridor_witness": True,
            "service_adds_recursion_depth": False,
            "requires_physical_energy_reentry": False,
            "coefficient_impulses_used_as_work": False,
        }
        corridor = quotient_full_natural_service_outcome(
            outcome,
            event_time=t,
            renewal_frequency=A,
            scaled_lifetime=c,
        )
        time_relative_residual = abs(corridor.physical_time_drop - T) / T
        wt = max(wt, time_relative_residual)
        if time_relative_residual > 8e-12:
            raise AssertionError("full-natural service acquired a second physical-time edge")

        n = rng.randint(2, 40)
        weights = [math.exp(rng.uniform(-6.0, 1.0)) for _ in range(n)]
        weight_total = sum(weights)
        if weight_total < corridor.integrated_service_lower:
            factor = 1.1 * corridor.integrated_service_lower / weight_total
            weights = [factor * value for value in weights]
        old0 = [bool(rng.getrandbits(1)) for _ in range(n)]
        old1 = [bool(rng.getrandbits(1)) for _ in range(n)]
        part = material_partition_is_same_corridor_measure(corridor, weights, old0, old1)
        partition_scale = max(float(part["total"]), corridor.integrated_service_lower)
        wp = max(wp, abs(float(part["partition_residual"])) / partition_scale)
        if int(part["recursion_edges_added"]) != 0 or bool(part["causal_charge_created"]):
            raise AssertionError("material rereading manufactured recursion depth or causal charge")

        parent_M = A / RENEWAL_TO_PARENT_SHELL_RATIO
        carrier_mu = 10.0 ** rng.uniform(-220.0, 80.0)
        cover = endpoint_comparable_hard_shell_cover(
            parent_shell_frequency=parent_M,
            endpoint_carrier_critical_mass=carrier_mu,
        )
        lower = float(cover["guaranteed_max_hard_shell_critical_mass_lower"])
        # Construct exact hard masses whose cover dominates carrier_mu.
        x = rng.random()
        mu0 = x * carrier_mu
        mu1 = 2.0 * (1.0 - x) * carrier_mu
        if not _relative_certificate_close(mu0 + 0.5 * mu1, carrier_mu):
            raise AssertionError("stress fixture lost exact two-shell cover")
        actual = max(mu0, mu1)
        mh = min(mh, (actual - lower) / lower)
        wit = realized_endpoint_hard_shell_witnesses(cover, (mu0, mu1))
        if len(tuple(wit["joint_witness_frequencies"])) > 1:
            ties += 1
        if actual + 3e-12 * max(actual, lower) < lower:
            raise AssertionError("endpoint hard-shell witness lower failed")
        if float(wit["physical_time_drop_added"]) != 0.0:
            raise AssertionError("endpoint shell witness created artificial time")

        if i == 0:
            tie_cover = endpoint_comparable_hard_shell_cover(
                parent_shell_frequency=parent_M,
                endpoint_carrier_critical_mass=carrier_mu,
            )
            tie_mass = max(lower, carrier_mu)
            tie = realized_endpoint_hard_shell_witnesses(tie_cover, (tie_mass, tie_mass))
            if len(tuple(tie["joint_witness_frequencies"])) != 2:
                raise AssertionError("exact endpoint shell tie was broken by priority")
            ties += 1

    return ServiceCorridorStress(samples, wt, wp, mh, ties)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-full-natural-service-corridor-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "full_natural_service_corridor_quotient.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    md = f"""# Full-natural service corridor quotient

Status: **{cert['status']}**.

A generic critical shell which survives its complete backward natural interval has already traversed one physical corridor

`[t-c A^-2,t]`.

The positive bounded heat/increment service proved on that interval is physically real, but it is **not a second event**.  Reading the service, applying Moyal, and partitioning the same positive law into OO/ON/NN adds zero physical-time drop, zero recursion edge and zero second causal charge.

The endpoint smooth carrier is also already present at `t-cA^-2`.  Its transported support lies inside `(A/2,2A)`.  Therefore the two exact hard shells `(A/2,A]` and `(A,2A]` satisfy

`max(mu_A,mu_2A) >= (2/3) A||Q_Au||_2^2`.

With the generic-shell registration `A=3M/4`, the endpoint witness shell set lies at ratios `3/4` or `3/2` relative to the incoming shell.  Equal witnesses remain joint; the theorem introduces no frequency-order priority.  This is comparable endpoint geometry, not a monotone-scale theorem.

The scale types are never collapsed.  A witnessed hard shell `H` registers the next smooth corridor at `3H/4`; its natural duration is therefore computed from that next carrier frequency, not by silently reusing either the old parent shell or the old carrier scale.  The producer carries `M`, `A`, `c`, and the two physical endpoints, and downstream readers are forbidden to rebind them.

All numerical certificate slack is relative to the native quantity being checked.  There is no `max(1,...)` floor: changing physical time, frequency, or service units cannot turn an incomplete UV corridor or zero hard-shell mass into a valid witness.  Likewise OO/ON/NN may be marked as the same corridor law only when the complete positive edge measure realizes the corridor's integrated service lower.

Hence a master path which only inserts

`full natural corridor -> own-scale service -> Moyal/material rereading -> endpoint survivor`

has acquired only **one** physical recursion edge: the natural corridor itself.  All intermediate service/material theorem layers are same-corridor witness maps.  If such a corridor chain remains bounded in frequency, the existing physical-time telescope forces `t=0`; if it avoids `t=0`, its unresolved possibility is genuinely UV-unbounded rather than service-theorem recursion.

Stress: `{out.samples}` corridor/service/material/endpoint-shell states
- worst natural-time identity relative residual: `{out.worst_time_identity_relative_residual:.3e}`
- worst OO/ON/NN same-measure partition relative residual: `{out.worst_material_partition_relative_residual:.3e}`
- minimum two-hard-shell cover relative margin: `{out.minimum_hard_shell_cover_relative_margin:.3e}`
- exact joint shell-witness ties retained: `{out.joint_tie_witnesses}`

This theorem removes service-theorem depth from the continuum master.  It does not prove UV no-escape and does not terminate genuine first-hit, work, source or material-reuse owner recurrence.  No Navier--Stokes global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
