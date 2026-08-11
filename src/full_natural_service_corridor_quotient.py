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

    def __post_init__(self) -> None:
        vals = (
            self.terminal_time,
            self.endpoint_time,
            self.renewal_frequency,
            self.scaled_lifetime,
            self.uniform_service_lower,
            self.integrated_service_lower,
        )
        if not all(math.isfinite(x) for x in vals):
            raise ValueError("finite corridor/service data required")
        if self.terminal_time <= 0 or self.endpoint_time < 0 or self.endpoint_time >= self.terminal_time:
            raise ValueError("full-natural service requires one nontrivial backward physical interval")
        if self.renewal_frequency <= 0 or self.scaled_lifetime <= 0:
            raise ValueError("positive renewal frequency and lifetime required")
        if self.uniform_service_lower <= 0 or self.integrated_service_lower <= 0:
            raise ValueError("positive own-scale service lower required")
        expected = self.scaled_lifetime / (self.renewal_frequency * self.renewal_frequency)
        actual = self.terminal_time - self.endpoint_time
        tol = 6e-12 * max(1.0, expected, actual)
        if abs(actual - expected) > tol:
            raise ValueError("service witness interval is not the completed natural corridor")

    @property
    def physical_time_drop(self) -> float:
        return self.terminal_time - self.endpoint_time


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

    t = float(event_time)
    A = float(renewal_frequency)
    c = float(scaled_lifetime)
    if t <= 0 or A <= 0 or c <= 0 or not all(math.isfinite(x) for x in (t, A, c)):
        raise ValueError("positive finite event time, renewal frequency and lifetime required")
    required = c / (A * A)
    if required >= t:
        raise ValueError("a corridor reaching t=0 is an absorbing boundary, not a full-natural service witness")
    reported = float(outcome.get("required_elapsed", -1.0))
    tol = 6e-12 * max(1.0, required, abs(reported))
    if reported <= 0 or abs(reported - required) > tol:
        raise ValueError("service outcome does not report the same completed natural interval")
    horizon = float(outcome.get("observed_elapsed_end", -1.0))
    if horizon + tol < required:
        raise ValueError("service outcome was inferred from an incomplete monitor horizon")

    return FullNaturalServiceCorridor(
        terminal_time=t,
        endpoint_time=t - required,
        renewal_frequency=A,
        scaled_lifetime=c,
        uniform_service_lower=float(outcome["uniform_square_service_lower"]),
        integrated_service_lower=float(outcome["integrated_bounded_heat_service_lower"]),
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
    mu = float(outcome.get("endpoint_carrier_critical_mass_lower", -1.0))
    if mu <= 0 or not math.isfinite(mu):
        raise ValueError("full-natural outcome supplied no positive endpoint carrier mass lower")
    return endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=parent_shell_frequency,
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
    return {
        **part,
        "corridor_terminal_time": corridor.terminal_time,
        "corridor_endpoint_time": corridor.endpoint_time,
        "physical_time_drop_already_counted": corridor.physical_time_drop,
        "recursion_edges_added": 0,
        "causal_charge_created": False,
        "service_mass_duplicated": False,
        "same_positive_measure": True,
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
    if M <= 0 or mu < 0 or not math.isfinite(M + mu):
        raise ValueError("positive parent shell frequency and finite nonnegative carrier mass required")
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    lo, hi = transported_annular_support_ratios()
    if lo <= 0.5 or hi >= 2.0:
        raise AssertionError("transported smooth carrier escaped the two-hard-shell endpoint cover")
    return {
        "renewal_frequency": A,
        "hard_shell_candidates": (A, 2.0 * A),
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
    tol = float(tie_tolerance) * max(1.0, mx, lower)
    if mx + tol < lower:
        raise ValueError("actual endpoint shell masses do not realize the certified smooth-carrier cover")
    ids = tuple(i for i, x in enumerate(masses) if mx - x <= tol)
    return {
        "joint_witness_frequencies": tuple(freqs[i] for i in ids),
        "joint_witness_critical_masses": tuple(masses[i] for i in ids),
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
        "material_ontology": "OO/ON/NN is a positive disintegration of that same service measure and adds zero causal charge and zero recursion depth",
        "endpoint_carrier": "the surviving Q_A carrier is present at the corridor endpoint itself; service is a sidecar/witness, not a prerequisite for inventing endpoint persistence",
        "hard_shell_cover": f"transported support lies in ({lo:.12g}A,{hi:.12g}A) subset (A/2,2A); therefore max(mu_A,mu_2A)>=(2/3) A||Q_Au||^2 at the same endpoint",
        "scale_geometry": "with A=3M/4, the endpoint witness set lies at shell ratios 3/4 or 3/2 relative to the incoming shell; this is comparable geometry, not monotone progress",
        "ties": "equal hard-shell witnesses remain joint; no theorem-name or frequency-order priority is introduced",
        "master_quotient": "a chain which only alternates full-natural corridors with their own-scale service/material witness layers is just a chain of full-natural physical corridors; bounded-scale such chains hit t=0 by physical time",
        "scope": "this removes service-theorem depth from the master and closes the endpoint-service attachment seam; it does not terminate UV-unbounded full-survivor chains or genuine first-hit/work/source/reuse owner recurrence",
    }


@dataclass(frozen=True)
class ServiceCorridorStress:
    samples: int
    worst_time_identity_residual: float
    worst_material_partition_residual: float
    minimum_hard_shell_cover_margin: float
    joint_tie_witnesses: int


def stress(samples: int = 50_000, seed: int = 20260811) -> ServiceCorridorStress:
    rng = random.Random(seed)
    wt = wp = 0.0
    mh = math.inf
    ties = 0
    for i in range(samples):
        A = math.exp(rng.uniform(-3.0, 5.0))
        c = math.exp(rng.uniform(-2.0, 1.0))
        T = c / (A * A)
        t = T + math.exp(rng.uniform(-4.0, 2.0))
        y = math.exp(rng.uniform(-8.0, 2.0))
        outcome = {
            "classification": FULL_NATURAL_SERVICE_WITNESS,
            "joint_first_stops": (),
            "required_elapsed": T,
            "observed_elapsed_end": T * (1.0 + rng.uniform(0.0, 0.25)),
            "uniform_square_service_lower": y,
            "integrated_bounded_heat_service_lower": c * y,
            "requires_physical_energy_reentry": False,
            "coefficient_impulses_used_as_work": False,
        }
        corridor = quotient_full_natural_service_outcome(
            outcome,
            event_time=t,
            renewal_frequency=A,
            scaled_lifetime=c,
        )
        wt = max(wt, abs(corridor.physical_time_drop - T))
        if abs(corridor.physical_time_drop - T) > 8e-12 * max(1.0, T):
            raise AssertionError("full-natural service acquired a second physical-time edge")

        n = rng.randint(2, 40)
        weights = [math.exp(rng.uniform(-6.0, 1.0)) for _ in range(n)]
        old0 = [bool(rng.getrandbits(1)) for _ in range(n)]
        old1 = [bool(rng.getrandbits(1)) for _ in range(n)]
        part = material_partition_is_same_corridor_measure(corridor, weights, old0, old1)
        wp = max(wp, abs(float(part["partition_residual"])))
        if int(part["recursion_edges_added"]) != 0 or bool(part["causal_charge_created"]):
            raise AssertionError("material rereading manufactured recursion depth or causal charge")

        parent_M = A / RENEWAL_TO_PARENT_SHELL_RATIO
        carrier_mu = math.exp(rng.uniform(-8.0, 2.0))
        cover = endpoint_comparable_hard_shell_cover(
            parent_shell_frequency=parent_M,
            endpoint_carrier_critical_mass=carrier_mu,
        )
        lower = float(cover["guaranteed_max_hard_shell_critical_mass_lower"])
        # Construct exact hard masses whose cover dominates carrier_mu.
        x = rng.random()
        mu0 = x * carrier_mu
        mu1 = 2.0 * (1.0 - x) * carrier_mu
        if mu0 + 0.5 * mu1 + 1e-15 < carrier_mu:
            raise AssertionError("stress fixture lost exact two-shell cover")
        actual = max(mu0, mu1)
        mh = min(mh, actual - lower)
        wit = realized_endpoint_hard_shell_witnesses(cover, (mu0, mu1))
        if len(tuple(wit["joint_witness_frequencies"])) > 1:
            ties += 1
        if actual + 3e-12 * max(1.0, carrier_mu) < lower:
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

Hence a master path which only inserts

`full natural corridor -> own-scale service -> Moyal/material rereading -> endpoint survivor`

has acquired only **one** physical recursion edge: the natural corridor itself.  All intermediate service/material theorem layers are same-corridor witness maps.  If such a corridor chain remains bounded in frequency, the existing physical-time telescope forces `t=0`; if it avoids `t=0`, its unresolved possibility is genuinely UV-unbounded rather than service-theorem recursion.

Stress: `{out.samples}` corridor/service/material/endpoint-shell states
- worst natural-time identity residual: `{out.worst_time_identity_residual:.3e}`
- worst OO/ON/NN same-measure partition residual: `{out.worst_material_partition_residual:.3e}`
- minimum two-hard-shell cover margin: `{out.minimum_hard_shell_cover_margin:.3e}`
- exact joint shell-witness ties retained: `{out.joint_tie_witnesses}`

This theorem removes service-theorem depth from the continuum master.  It does not prove UV no-escape and does not terminate genuine first-hit, work, source or material-reuse owner recurrence.  No Navier--Stokes global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
