from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path

import numpy as np

from src.physical_branch_compiler import CAUSE_TO_CURRENCY, CauseHit, PhysicalCause, PhysicalCurrency
from src.physical_pair_weighted_productivity import conditioned_physical_log_productivity_constant


class RecursiveInternalRoute(str, Enum):
    HH_REGENERATION = "earlier_high_high_regeneration"
    CONTINUE_PRODUCTIVITY = "registered_physical_productivity_continuation"


@dataclass(frozen=True)
class RegenerationHit:
    """Earlier genuine HH-generation stop, measured against the same dT law.

    ``weight`` is required only on an exact first-time tie.  It must already be
    a positive Radon--Nikodym stopping weight in the same physical-transfer units
    as any tied ``CauseHit`` weights.  This module never manufactures a
    dimensionless comparison between heterogeneous observables.
    """

    time: float
    weight: float
    witness: str = ""


@dataclass(frozen=True)
class GeneratedPairEvent:
    """One atom/quadrature element of retained positive physical HH work.

    ``marking_good`` means the event lies on the retained complex-Young/phase
    branch.  If false, a physical first-stop cause (normally transfer loss) must
    already be supplied.  ``registration_good`` means both selected parent marks
    reach the common slice; if false, a physical cause or HH-regeneration hit
    must already be supplied.  No theorem-name fallback is invented here.
    """

    mass: float
    pair_cell: int
    marking_good: bool
    registration_good: bool
    physical_hits: tuple[CauseHit, ...] = ()
    regeneration_hits: tuple[RegenerationHit, ...] = ()


@dataclass(frozen=True)
class GeneratedMeasurePartition:
    total_mass: float
    xi_mass: float
    retained_generated_mass: float
    currency_mass: dict[str, float]
    regeneration_mass: float
    continuation_mass: float
    continuation_fraction: float
    pair_cells: int
    conditioned_productivity: float | None
    majority_continues: bool
    exact_tie_events: int


def _validate_event(event: GeneratedPairEvent) -> None:
    if not math.isfinite(event.mass) or event.mass <= 0:
        raise ValueError("generated physical event mass must be positive and finite")
    if event.pair_cell < 0:
        raise ValueError("pair-cell index must be nonnegative")
    for hit in event.physical_hits:
        if not math.isfinite(hit.time) or not math.isfinite(hit.weight) or hit.weight <= 0:
            raise ValueError("physical cause hits require finite time and positive RN weight")
    for hit in event.regeneration_hits:
        if not math.isfinite(hit.time) or not math.isfinite(hit.weight) or hit.weight <= 0:
            raise ValueError("regeneration hits require finite time and positive RN weight")
    if not event.marking_good:
        if not any(h.cause is PhysicalCause.TRANSFER_WORK_LOSS for h in event.physical_hits):
            raise ValueError(
                "a complex-Young/phase-bad event must carry its certified physical transfer-loss first stop"
            )
    if not event.registration_good and not event.physical_hits and not event.regeneration_hits:
        raise ValueError(
            "a failed common-slice registration must carry its physical cause or earlier HH-regeneration stop"
        )


def _first_event_split(event: GeneratedPairEvent) -> tuple[dict[str, float], float, bool]:
    """Split one event at its first physical/internal causal time.

    Returns route fractions, first time and whether an exact independent tie was
    split.  Duplicate physical manifestations of one cause are quotiented before
    normalization.  The initial boundary at t=0 remains absorbing.
    """
    _validate_event(event)
    all_times = [h.time for h in event.physical_hits] + [h.time for h in event.regeneration_hits]
    if not all_times:
        if not (event.marking_good and event.registration_good):
            raise AssertionError("unregistered event escaped without a causal stop")
        return {RecursiveInternalRoute.CONTINUE_PRODUCTIVITY.value: 1.0}, math.inf, False

    first = min(all_times)
    tol = 64.0 * math.ulp(max(1.0, abs(first)))
    ph = [h for h in event.physical_hits if abs(h.time-first) <= tol]
    rh = [h for h in event.regeneration_hits if abs(h.time-first) <= tol]

    if first <= tol and any(h.cause is PhysicalCause.INITIAL_BOUNDARY for h in ph):
        return {PhysicalCurrency.INITIAL_BOUNDARY.value: 1.0}, first, len(ph)+len(rh) > 1

    cause_weights: dict[PhysicalCause, float] = {}
    for h in ph:
        cause_weights[h.cause] = cause_weights.get(h.cause, 0.0) + h.weight
    regen_weight = sum(h.weight for h in rh)
    total = sum(cause_weights.values()) + regen_weight
    if total <= 0:
        raise AssertionError("first-stop RN quotient has zero mass")

    routes: dict[str, float] = {}
    for cause, weight in cause_weights.items():
        key = CAUSE_TO_CURRENCY[cause].value
        routes[key] = routes.get(key, 0.0) + weight/total
    if regen_weight:
        key = RecursiveInternalRoute.HH_REGENERATION.value
        routes[key] = routes.get(key, 0.0) + regen_weight/total
    if not math.isclose(sum(routes.values()), 1.0, rel_tol=2e-14, abs_tol=2e-14):
        raise AssertionError("first-stop route fractions do not conserve event mass")
    return routes, first, len(cause_weights) + int(regen_weight > 0) > 1


def compile_generated_pair_measure(
    *,
    events: tuple[GeneratedPairEvent, ...],
    xi_mass: float = 0.0,
    pair_cells_upper: int | None = None,
    scaled_lifetime: float = 1.0,
) -> GeneratedMeasurePartition:
    """Exact single-charge partition of a generated positive HH work measure.

    The input events already represent the retained actual positive HH
    child-work law **after** support-level Xi excision.  ``xi_mass`` is the
    disjoint physical/interface measure removed before these events.  No event
    is rescaled: total physical mass is exactly Xi plus the retained event mass.
    """
    if not events:
        raise ValueError("at least one generated physical-work event is required")
    for event in events:
        _validate_event(event)
    retained = sum(e.mass for e in events)
    if not math.isfinite(xi_mass) or xi_mass < 0:
        raise ValueError("Xi mass must be finite and nonnegative")
    raw_total = retained + xi_mass

    currency: dict[str, float] = {}
    if xi_mass:
        currency[PhysicalCurrency.XI.value] = xi_mass
    regen = 0.0
    cont = 0.0
    ties = 0

    for event in events:
        routes, _, tied = _first_event_split(event)
        mass = event.mass
        ties += int(tied)
        for key, frac in routes.items():
            amount = mass*frac
            if key == RecursiveInternalRoute.HH_REGENERATION.value:
                regen += amount
            elif key == RecursiveInternalRoute.CONTINUE_PRODUCTIVITY.value:
                cont += amount
            else:
                currency[key] = currency.get(key, 0.0) + amount

    accounted = xi_mass + regen + cont + sum(v for k,v in currency.items() if k != PhysicalCurrency.XI.value)
    if not math.isclose(accounted, raw_total, rel_tol=3e-13, abs_tol=3e-13*raw_total):
        raise AssertionError("generated first-stop constructor lost or double-charged physical work")

    # Survival is measured inside the retained generated HH law, because that is
    # the law used by the physical KL productivity theorem.
    q = cont/retained if retained > 0 else 0.0
    cells_seen = 1 + max(e.pair_cell for e in events)
    M = cells_seen if pair_cells_upper is None else int(pair_cells_upper)
    if M < cells_seen:
        raise ValueError("pair_cells_upper cannot be smaller than observed event cells")
    lam = None
    if q > 0:
        lam = conditioned_physical_log_productivity_constant(q, M, scaled_lifetime)

    return GeneratedMeasurePartition(
        total_mass=raw_total,
        xi_mass=xi_mass,
        retained_generated_mass=retained,
        currency_mass=currency,
        regeneration_mass=regen,
        continuation_mass=cont,
        continuation_fraction=q,
        pair_cells=M,
        conditioned_productivity=lam,
        majority_continues=q >= 0.5,
        exact_tie_events=ties,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_GENERATED_PHYSICAL_MEASURE_FIRST_STOP_AND_SURVIVAL_PRODUCTIVITY__CAUSE_RN_EXTRACTION_REMAINS_PDE_BRIDGE",
        "domain": "the constructor acts only on retained actual positive HH child-energy work after the physical-energy generation gate; it does not synthesize a new probability law",
        "first_stop": "each event is routed at its first causal time to an existing physical currency, earlier HH regeneration, or registered continuation; duplicate manifestations of one physical cause are quotiented",
        "ties": "exact independent ties are split only from positive RN stopping weights already expressed against the same physical-transfer law; heterogeneous theorem thresholds are never normalized by fiat",
        "boundary": "t=0 is absorbing and wins any zero-time tie; it is never relabeled as fresh interior ancestry",
        "young_failure": "a complex-Young/phase-bad event cannot disappear: it must already carry its certified physical transfer-loss stop",
        "registration_failure": "a failed common-slice parent mark must carry classified source/relink or earlier HH-regeneration provenance",
        "survival": "if continuation carries fraction q of retained generated physical work, its exact physical weighted productivity is Lambda_survivor=q Lambda_full",
        "half_gate": "q<1/2 means a majority of current physical HH work already left free continuation through named first stops or earlier generation; q>=1/2 costs at most log2 in the productivity offset",
        "single_charge": "Xi is excised once; currency mass + HH-regeneration mass + continuation mass reconstruct the original generated physical work exactly",
        "continuum_status": "remaining local PDE bridge is to extract measurable eventwise cause RN stopping densities on the actual dT space, especially on exact tie sets; the constructor refuses lexicographic or unit-incompatible tie weights",
    }


@dataclass(frozen=True)
class ConstructorStress:
    samples: int
    worst_mass_residual: float
    minimum_positive_conditioned_productivity: float
    maximum_tie_fraction: float
    majority_continue_samples: int
    majority_stopped_or_regenerated_samples: int
    zero_survivor_samples: int


def _random_event(rng: np.random.Generator, cell: int) -> GeneratedPairEvent:
    mass = float(math.exp(rng.uniform(-4.0, 3.0)))
    mode = int(rng.integers(0, 7))
    t = float(rng.uniform(0.01, 1.0))
    if mode == 0:
        return GeneratedPairEvent(mass, cell, True, True)
    if mode == 1:
        return GeneratedPairEvent(
            mass, cell, False, True,
            physical_hits=(CauseHit(t, PhysicalCause.TRANSFER_WORK_LOSS, float(rng.uniform(0.1,2.0)), "Young/phase deficit"),),
        )
    if mode == 2:
        return GeneratedPairEvent(
            mass, cell, True, False,
            physical_hits=(CauseHit(t, PhysicalCause.RESOLVED_SOURCE, float(rng.uniform(0.1,2.0)), "classified residual"),),
        )
    if mode == 3:
        return GeneratedPairEvent(
            mass, cell, True, False,
            regeneration_hits=(RegenerationHit(t, float(rng.uniform(0.1,2.0)), "earlier HH generation"),),
        )
    if mode == 4:
        # Exact independent tie between source and regeneration.
        return GeneratedPairEvent(
            mass, cell, True, False,
            physical_hits=(CauseHit(t, PhysicalCause.RESOLVED_SOURCE, float(rng.uniform(0.1,2.0)), "source tie"),),
            regeneration_hits=(RegenerationHit(t, float(rng.uniform(0.1,2.0)), "HH tie"),),
        )
    if mode == 5:
        return GeneratedPairEvent(
            mass, cell, True, False,
            physical_hits=(CauseHit(t, PhysicalCause.MATERIAL_RELINK, float(rng.uniform(0.1,2.0)), "material relink"),),
        )
    return GeneratedPairEvent(
        mass, cell, True, False,
        physical_hits=(CauseHit(0.0, PhysicalCause.INITIAL_BOUNDARY, 1.0, "initial boundary"),),
        regeneration_hits=(RegenerationHit(0.0, 1.0, "zero-time synthetic tie"),),
    )


def stress(samples: int = 50_000, seed: int = 20260809) -> ConstructorStress:
    rng = np.random.default_rng(seed)
    worst = 0.0
    minlam = float("inf")
    maxtie = 0.0
    majc = majs = zero = 0
    for _ in range(samples):
        M = int(rng.integers(1, 12))
        n = int(rng.integers(1, 24))
        events = tuple(_random_event(rng, int(rng.integers(0,M))) for _ in range(n))
        total = sum(e.mass for e in events)
        xi = float(rng.uniform(0.0, 0.15))*total
        out = compile_generated_pair_measure(events=events, xi_mass=xi, pair_cells_upper=M)
        currency_nonxi = sum(v for k,v in out.currency_mass.items() if k != PhysicalCurrency.XI.value)
        recon = out.xi_mass + currency_nonxi + out.regeneration_mass + out.continuation_mass
        resid = abs(recon-out.total_mass)/max(1.0,out.total_mass)
        worst=max(worst,resid)
        if resid>3e-12:
            raise AssertionError("generated first-stop mass reconstruction failed")
        if out.conditioned_productivity is not None:
            minlam=min(minlam,out.conditioned_productivity)
            if out.conditioned_productivity<=0:
                raise AssertionError("surviving physical law lost positive productivity")
        if out.majority_continues:
            majc += 1
        else:
            majs += 1
        if out.continuation_fraction==0:
            zero += 1
        maxtie=max(maxtie,out.exact_tie_events/max(1,n))
        # Zero-time initial boundary is absorbing even against regeneration.
        b = GeneratedPairEvent(
            1.0,0,True,False,
            physical_hits=(CauseHit(0.0,PhysicalCause.INITIAL_BOUNDARY,0.2,"boundary"),),
            regeneration_hits=(RegenerationHit(0.0,5.0,"regen"),),
        )
        bo=compile_generated_pair_measure(events=(b,),pair_cells_upper=1)
        if bo.currency_mass.get(PhysicalCurrency.INITIAL_BOUNDARY.value,0.0)!=1.0 or bo.regeneration_mass!=0:
            raise AssertionError("initial boundary did not absorb zero-time tie")
    if minlam==float("inf"):
        minlam=0.0
    return ConstructorStress(samples,worst,minlam,maxtie,majc,majs,zero)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--outdir",type=Path,default=Path("results-recursive-physical-witness-constructor"))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); out=stress(args.samples)
    (args.outdir/"recursive_physical_witness_constructor.json").write_text(
        json.dumps({"certificate":cert,"stress":asdict(out)},indent=2),encoding="utf-8"
    )
    md=f"""# Recursive generated physical-work first-stop constructor\n\nStatus: **{cert['status']}**.\n\nThe constructor acts on the retained **actual positive high--high child-energy work measure** after the physical-energy gate.  It never invents a packet probability.  Each physical event is routed at its first causal time to exactly one of: an existing physical currency (with exact RN splitting on independent ties), an earlier HH-generation recursion, or a registered Young-good continuation.\n\n`Xi` is excised once before causal routing.  On the remaining generated measure, duplicate theorem manifestations of one physical cause are quotiented.  Exact ties are split only when all supplied stopping weights are Radon--Nikodym weights against the same physical `dT`; the constructor refuses to normalize heterogeneous threshold excesses.  At `t=0` the initial boundary is absorbing.\n\nLet `q` be the fraction of retained generated physical work which survives every earlier causal stop and both complex-Young/phase and common-slice registration.  Restricting the physical KL productivity theorem to that survivor law gives exactly\n\n`Lambda_survivor = q Lambda_full`.\n\nThus if `q>=1/2`, the entire survival conditioning costs at most `log 2` in the logarithmic amplitude offset.  If `q<1/2`, a majority of the current physical HH work has already left free continuation through an existing cause or an earlier HH-generation recursion.  There is no fourth free branch.\n\nStress: `{out.samples}` random physical-work measures\n- worst total-mass reconstruction residual: `{out.worst_mass_residual:.3e}`\n- minimum positive survivor productivity: `{out.minimum_positive_conditioned_productivity:.3e}`\n- maximum exact-tie event fraction: `{out.maximum_tie_fraction:.3f}`\n- majority-continuation samples: `{out.majority_continue_samples}`\n- majority-stopped/regenerated samples: `{out.majority_stopped_or_regenerated_samples}`\n- zero-survivor samples: `{out.zero_survivor_samples}`\n\nThis closes the **measure-level survival/first-stop algebra** of the generated branch.  The remaining local continuum bridge is now sharply stated: construct measurable eventwise cause stopping densities on the actual `dT` space, in common transfer units, for exact tie sets.  Source/relink/transfer causes already possess physical work or Moyal measures in their own theorems, but the universal eventwise RN assembly has not yet been proved.  No lexicographic substitute is used, and no Navier--Stokes global-regularity conclusion is claimed.\n"""
    (args.outdir/"summary.md").write_text(md,encoding="utf-8"); print(md)


if __name__=="__main__":
    main()
