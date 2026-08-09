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
from src.joint_causal_stop_projection import InternalHit, joint_stop_master_projection
from src.physical_branch_compiler import MasterDisposition, UniformResourceCertificate


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


@dataclass(frozen=True)
class GeneratedMasterPartition:
    total_mass: float
    xi_mass: float
    retained_generated_mass: float
    master_mass: dict[str, float]
    continuation_mass: float
    continuation_fraction: float
    pair_cells: int
    conditioned_productivity: float | None
    majority_continues: bool
    joint_stop_events: int
    maximum_joint_cause_count: int


def _validate_event(event: GeneratedPairEvent, *, require_fine_weights: bool = True) -> None:
    if not math.isfinite(event.mass) or event.mass <= 0:
        raise ValueError("generated physical event mass must be positive and finite")
    if event.pair_cell < 0:
        raise ValueError("pair-cell index must be nonnegative")
    for hit in event.physical_hits:
        if not math.isfinite(hit.time):
            raise ValueError("physical cause hits require finite time")
        if require_fine_weights and (not math.isfinite(hit.weight) or hit.weight <= 0):
            raise ValueError("fine RN physical cause hits require positive finite weights")
    for hit in event.regeneration_hits:
        if not math.isfinite(hit.time):
            raise ValueError("regeneration hits require finite time")
        if require_fine_weights and (not math.isfinite(hit.weight) or hit.weight <= 0):
            raise ValueError("fine RN regeneration hits require positive finite weights")
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


def compile_generated_pair_master_measure(
    *,
    events: tuple[GeneratedPairEvent, ...],
    xi_mass: float = 0.0,
    pair_cells_upper: int | None = None,
    scaled_lifetime: float = 1.0,
    uniform_certificates: dict[PhysicalCause, UniformResourceCertificate] | None = None,
) -> GeneratedMasterPartition:
    """Preferred unsplit generated measure compiler for the no-escape master.

    Exact simultaneous causes are kept as one joint stop.  CauseHit/RegenerationHit
    weights are ignored by the master path; they exist only for the optional fine
    RN subledger.  Thus changing tie weights cannot change any master mass.
    """
    if not events:
        raise ValueError("at least one generated physical-work event is required")
    for event in events:
        _validate_event(event, require_fine_weights=False)
    retained = sum(e.mass for e in events)
    if not math.isfinite(xi_mass) or xi_mass < 0:
        raise ValueError("Xi mass must be finite and nonnegative")
    total = retained + xi_mass
    master: dict[str, float] = {}
    if xi_mass:
        master[MasterDisposition.XI.value] = xi_mass
    cont = 0.0
    joint_events = 0
    max_joint = 0

    for event in events:
        if not event.physical_hits and not event.regeneration_hits:
            if not (event.marking_good and event.registration_good):
                raise AssertionError("unregistered event escaped without a joint stop")
            cont += event.mass
            continue
        internal = tuple(InternalHit(h.time, witness=h.witness) for h in event.regeneration_hits)
        out = joint_stop_master_projection(
            physical_hits=event.physical_hits,
            internal_hits=internal,
            uniform_certificates=uniform_certificates,
        )
        key = out.master_disposition
        master[key] = master.get(key, 0.0) + event.mass
        joint_events += 1
        max_joint = max(max_joint, len(out.joint_physical_causes) + len(out.joint_internal_causes))

    accounted = cont + sum(master.values())
    if not math.isclose(accounted, total, rel_tol=3e-13, abs_tol=3e-13*max(1.0,total)):
        raise AssertionError("unsplit joint-stop master lost or duplicated generated physical work")
    q = cont/retained if retained > 0 else 0.0
    cells_seen = 1 + max(e.pair_cell for e in events)
    M = cells_seen if pair_cells_upper is None else int(pair_cells_upper)
    if M < cells_seen:
        raise ValueError("pair_cells_upper cannot be smaller than observed event cells")
    lam = conditioned_physical_log_productivity_constant(q,M,scaled_lifetime) if q>0 else None
    return GeneratedMasterPartition(
        total_mass=total, xi_mass=xi_mass, retained_generated_mass=retained,
        master_mass=master, continuation_mass=cont, continuation_fraction=q,
        pair_cells=M, conditioned_productivity=lam, majority_continues=q>=0.5,
        joint_stop_events=joint_events, maximum_joint_cause_count=max_joint,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_GENERATED_PHYSICAL_FIRST_STOP_SURVIVAL__JOINT_MASTER_UNSPLIT__FIRST_HIT_EXTRACTION_REMAINS",
        "domain": "the constructor acts only on retained actual positive HH child-energy work after the physical-energy generation gate",
        "preferred_ties": "exact simultaneous causes remain one unsplit joint physical stop; master projection is invariant under CauseHit weights and requires no common-unit RN density",
        "fine_subledger": "the older RN split remains available only as optional fine-currency bookkeeping when physical common-unit weights happen to exist",
        "boundary": "t=0 is absorbing",
        "young_failure": "a complex-Young/phase-bad event must already carry its certified physical transfer-loss first stop",
        "registration_failure": "a failed common-slice parent mark must carry classified source/relink or earlier HH-regeneration provenance",
        "survival": "if continuation carries fraction q of retained generated physical work, its exact physical weighted productivity is Lambda_survivor=q Lambda_full",
        "half_gate": "q<1/2 means a majority of current physical HH work has already stopped; q>=1/2 loses at most log2 in productivity offset",
        "single_charge": "Xi plus one unsplit master disposition per stopped event plus survivor mass reconstruct the physical generated work exactly",
        "continuum_status": "remaining PDE bridge is measurable first-hit cause-set extraction from actual smooth-SGS observables; pair weights and exact-tie weights are no longer required",
    }


@dataclass(frozen=True)
class ConstructorStress:
    samples: int
    worst_mass_residual: float
    worst_unsplit_master_mass_residual: float
    minimum_positive_conditioned_productivity: float
    maximum_tie_fraction: float
    majority_continue_samples: int
    majority_stopped_or_regenerated_samples: int
    zero_survivor_samples: int
    unsplit_weight_invariance_failures: int


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
    worst_master = 0.0
    weight_fail = 0
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

        # Preferred master path keeps ties unsplit and must ignore all dummy weights.
        mout = compile_generated_pair_master_measure(events=events, xi_mass=xi, pair_cells_upper=M)
        mrecon = mout.continuation_mass + sum(mout.master_mass.values())
        mresid = abs(mrecon-mout.total_mass)/max(1.0,mout.total_mass)
        worst_master=max(worst_master,mresid)
        if mresid>3e-12:
            raise AssertionError("unsplit generated master mass reconstruction failed")
        mutated=[]
        for e in events:
            ph=tuple(CauseHit(h.time,h.cause,float(math.exp(rng.uniform(-20,20))),h.witness) for h in e.physical_hits)
            rh=tuple(RegenerationHit(h.time,float(math.exp(rng.uniform(-20,20))),h.witness) for h in e.regeneration_hits)
            mutated.append(GeneratedPairEvent(e.mass,e.pair_cell,e.marking_good,e.registration_good,ph,rh))
        mout2=compile_generated_pair_master_measure(events=tuple(mutated), xi_mass=xi, pair_cells_upper=M)
        if mout.master_mass != mout2.master_mass or not math.isclose(mout.continuation_mass,mout2.continuation_mass):
            weight_fail += 1
            raise AssertionError("preferred generated master depended on arbitrary tie weights")
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
    return ConstructorStress(samples,worst,worst_master,minlam,maxtie,majc,majs,zero,weight_fail)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--outdir",type=Path,default=Path("results-recursive-physical-witness-constructor"))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); out=stress(args.samples)
    (args.outdir/"recursive_physical_witness_constructor.json").write_text(
        json.dumps({"certificate":cert,"stress":asdict(out)},indent=2),encoding="utf-8"
    )
    md=f"""# Recursive generated physical-work first-stop constructor\n\nStatus: **{cert['status']}**.\n\nThe constructor acts on the retained **actual positive high--high child-energy work measure** after the physical-energy gate.  It never invents a packet probability.  The preferred master path keeps every exact simultaneous first hit as one unsplit joint physical stop; legacy RN splitting is optional fine bookkeeping only.\n\n`Xi` is excised once before causal routing.  The joint-stop theorem projects the entire stopped event mass to one coarse master fate by terminal semantics; source/strain/relink/HH-regeneration ties simply remain recursive.  Arbitrary dummy tie weights cannot alter the preferred master.  At `t=0` the initial boundary is absorbing.\n\nLet `q` be the fraction of retained generated physical work which survives every earlier causal stop and both complex-Young/phase and common-slice registration.  Restricting the physical KL productivity theorem to that survivor law gives exactly\n\n`Lambda_survivor = q Lambda_full`.\n\nThus if `q>=1/2`, the entire survival conditioning costs at most `log 2` in the logarithmic amplitude offset.  If `q<1/2`, a majority of the current physical HH work has already left free continuation through an existing cause or an earlier HH-generation recursion.  There is no fourth free branch.\n\nStress: `{out.samples}` random physical-work measures\n- worst total-mass reconstruction residual: `{out.worst_mass_residual:.3e}`\n- minimum positive survivor productivity: `{out.minimum_positive_conditioned_productivity:.3e}`\n- maximum exact-tie event fraction: `{out.maximum_tie_fraction:.3f}`\n- majority-continuation samples: `{out.majority_continue_samples}`\n- majority-stopped/regenerated samples: `{out.majority_stopped_or_regenerated_samples}`\n- zero-survivor samples: `{out.zero_survivor_samples}`\n\nThis closes the **measure-level survival/first-stop algebra** of the generated branch without requiring exact-tie fractions.  The remaining local continuum bridge is now only to construct the measurable first-hit **cause set** from actual smooth-SGS observables and verify that its no-hit set is the registered survivor set.  No common-unit RN assembly is required, and no Navier--Stokes global-regularity conclusion is claimed.\n"""
    (args.outdir/"summary.md").write_text(md,encoding="utf-8"); print(md)


if __name__=="__main__":
    main()
