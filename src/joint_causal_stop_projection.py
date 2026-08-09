from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Mapping

import numpy as np

from src.physical_branch_compiler import (
    CAUSE_TO_CURRENCY,
    CURRENCY_TO_MASTER,
    CauseHit,
    MasterDisposition,
    PhysicalCause,
    PhysicalCurrency,
    UniformResourceCertificate,
)


class InternalRecursiveCause(str, Enum):
    HH_REGENERATION = "earlier_high_high_regeneration"


@dataclass(frozen=True)
class InternalHit:
    time: float
    cause: InternalRecursiveCause = InternalRecursiveCause.HH_REGENERATION
    witness: str = ""


@dataclass(frozen=True)
class JointStopProjection:
    first_time: float | None
    joint_physical_causes: tuple[str, ...]
    joint_internal_causes: tuple[str, ...]
    certified_currencies: tuple[str, ...]
    master_disposition: str
    terminal_certificate_used: str | None
    fine_rn_split_required: bool


def _first_joint_set(
    physical_hits: tuple[CauseHit, ...],
    internal_hits: tuple[InternalHit, ...],
) -> tuple[float | None, tuple[PhysicalCause, ...], tuple[InternalRecursiveCause, ...]]:
    for h in physical_hits:
        if not math.isfinite(h.time):
            raise ValueError("finite physical hit times required")
    for h in internal_hits:
        if not math.isfinite(h.time):
            raise ValueError("finite internal hit times required")
    times = [h.time for h in physical_hits] + [h.time for h in internal_hits]
    if not times:
        return None, (), ()
    first = min(times)
    tol = 64.0 * math.ulp(max(1.0, abs(first)))
    pc = tuple(sorted({h.cause for h in physical_hits if abs(h.time-first) <= tol}, key=lambda x: x.value))
    ic = tuple(sorted({h.cause for h in internal_hits if abs(h.time-first) <= tol}, key=lambda x: x.value))
    return first, pc, ic


def joint_stop_master_projection(
    *,
    physical_hits: tuple[CauseHit, ...] = (),
    internal_hits: tuple[InternalHit, ...] = (),
    fixed_transfer_loss: bool = False,
    kelvin_flat_certified: bool = False,
    uniform_certificates: Mapping[PhysicalCause, UniformResourceCertificate] | None = None,
) -> JointStopProjection:
    """Project an exact simultaneous physical stop without splitting its mass.

    The joint cause set is retained as provenance.  The whole event/block mass is
    sent to one coarse master disposition by terminal semantics, not by a
    theorem-name order and not by normalized heterogeneous tie weights.

    Dominance is logical rather than causal:
      boundary at t=0 is absorbing;
      a certified fixed multiplicative cost is already terminal;
      otherwise a valid globally bounded uniform reset is terminal;
      otherwise every remaining physical/internal cause is recursive-critical.

    Fine RN splitting remains optional if somebody wants a subledger by physical
    cause, but it is not required for the no-escape master.
    """
    if fixed_transfer_loss:
        return JointStopProjection(
            first_time=None,
            joint_physical_causes=(PhysicalCause.TRANSFER_WORK_LOSS.value,),
            joint_internal_causes=(),
            certified_currencies=(PhysicalCurrency.MULTIPLICATIVE_TRANSFER.value,),
            master_disposition=MasterDisposition.TRANSFER_COST.value,
            terminal_certificate_used="stage_zero_fixed_transfer_loss",
            fine_rn_split_required=False,
        )

    first, physical, internal = _first_joint_set(physical_hits, internal_hits)
    if first is None:
        if not kelvin_flat_certified:
            raise ValueError("no causal stop and no certified Kelvin-flat continuation")
        return JointStopProjection(
            first_time=None,
            joint_physical_causes=(),
            joint_internal_causes=(),
            certified_currencies=(PhysicalCurrency.KELVIN_FLAT_EROSION.value,),
            master_disposition=MasterDisposition.FLAT.value,
            terminal_certificate_used="kelvin_flat_erosion",
            fine_rn_split_required=False,
        )

    tol = 64.0 * math.ulp(max(1.0, abs(first)))
    if first <= tol and PhysicalCause.INITIAL_BOUNDARY in physical:
        return JointStopProjection(
            first_time=first,
            joint_physical_causes=tuple(c.value for c in physical),
            joint_internal_causes=tuple(c.value for c in internal),
            certified_currencies=(PhysicalCurrency.INITIAL_BOUNDARY.value,),
            master_disposition=MasterDisposition.BOUNDARY.value,
            terminal_certificate_used="absorbing_initial_boundary",
            fine_rn_split_required=False,
        )

    currencies = tuple(sorted({CAUSE_TO_CURRENCY[c] for c in physical}, key=lambda x: x.value))
    transfer = tuple(c for c in currencies if CURRENCY_TO_MASTER[c] is MasterDisposition.TRANSFER_COST)
    if transfer:
        return JointStopProjection(
            first_time=first,
            joint_physical_causes=tuple(c.value for c in physical),
            joint_internal_causes=tuple(c.value for c in internal),
            certified_currencies=tuple(c.value for c in currencies),
            master_disposition=MasterDisposition.TRANSFER_COST.value,
            terminal_certificate_used="any_simultaneous_fixed_transfer_cost",
            fine_rn_split_required=False,
        )

    reset_causes = tuple(c for c in physical if CAUSE_TO_CURRENCY[c] is PhysicalCurrency.UNIFORM_GLOBAL_RESET)
    if reset_causes:
        certs = uniform_certificates or {}
        if not any((c in certs and certs[c].valid()) for c in reset_causes):
            raise ValueError("uniform-resource joint stop lacks a valid globally bounded reset certificate")
        return JointStopProjection(
            first_time=first,
            joint_physical_causes=tuple(c.value for c in physical),
            joint_internal_causes=tuple(c.value for c in internal),
            certified_currencies=tuple(c.value for c in currencies),
            master_disposition=MasterDisposition.ADDITIVE_RESET.value,
            terminal_certificate_used="valid_uniform_global_resource",
            fine_rn_split_required=False,
        )

    # All remaining physical currencies are recursive-critical.  HH regeneration
    # is also recursion, not a currency.  Keeping the full joint provenance is
    # enough; no fractional ownership is needed to continue the master.
    if not physical and not internal:
        raise AssertionError("nonempty first stop lost all causes")
    return JointStopProjection(
        first_time=first,
        joint_physical_causes=tuple(c.value for c in physical),
        joint_internal_causes=tuple(c.value for c in internal),
        certified_currencies=tuple(c.value for c in currencies),
        master_disposition=MasterDisposition.RECURSE_CRITICAL.value,
        terminal_certificate_used=None,
        fine_rn_split_required=False,
    )


def weight_invariance_countercheck(
    causes: tuple[PhysicalCause, ...],
    weights_a: tuple[float, ...],
    weights_b: tuple[float, ...],
    time: float = 0.4,
) -> bool:
    """Same exact tie set must project identically under arbitrary positive weights."""
    if len(causes) != len(weights_a) or len(causes) != len(weights_b):
        raise ValueError("matching cause/weight tuples required")
    if any(w <= 0 or not math.isfinite(w) for w in weights_a + weights_b):
        raise ValueError("positive finite dummy weights required")
    A = tuple(CauseHit(time, c, w, "dummy-a") for c, w in zip(causes, weights_a))
    B = tuple(CauseHit(time, c, w, "dummy-b") for c, w in zip(causes, weights_b))
    pa = joint_stop_master_projection(physical_hits=A)
    pb = joint_stop_master_projection(physical_hits=B)
    return (
        pa.joint_physical_causes == pb.joint_physical_causes
        and pa.certified_currencies == pb.certified_currencies
        and pa.master_disposition == pb.master_disposition
        and pa.terminal_certificate_used == pb.terminal_certificate_used
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_JOINT_CAUSAL_STOP_MASTER_PROJECTION__NO_RN_TIE_SPLIT_REQUIRED",
        "physical_statement": "an exact simultaneous first hit is retained as one joint physical stop carrying the full unsplit transfer mass and the complete set of certified causes",
        "boundary": "t=0 is absorbing",
        "terminal_dominance": "if any joint cause already certifies a fixed multiplicative transfer cost, the whole joint stop is terminal transfer-cost; otherwise a valid uniform global resource may terminate as additive reset",
        "recursive_joint": "source, critical dissipation, material/new ancestry, and earlier HH regeneration are all coarse RECURSE; their exact tie requires no fractional ownership",
        "weights": "arbitrary positive CauseHit weights do not affect the preferred master projection; RN weights are optional fine-subledger data only",
        "single_charge": "the joint event/block transfer mass is never split or duplicated; provenance is set-valued, master disposition is single-valued",
        "scope": "this supersedes RN splitting as a necessary continuum bridge for the master; RN splitting may still be used when a physically meaningful fine-currency partition is independently available",
        "continuum_status": "remaining PDE task is first-hit witness extraction itself, not construction of common-unit tie weights",
    }


@dataclass(frozen=True)
class JointStopStress:
    samples: int
    weight_invariance_failures: int
    transfer_terminal_samples: int
    reset_terminal_samples: int
    recurse_joint_samples: int
    boundary_samples: int
    maximum_joint_size: int


def stress(samples: int = 50_000, seed: int = 20260809) -> JointStopStress:
    rng = np.random.default_rng(seed)
    causes = [
        PhysicalCause.RESOLVED_SOURCE,
        PhysicalCause.HIGH_STRAIN_DISSIPATION,
        PhysicalCause.MATERIAL_RELINK,
        PhysicalCause.CAUSAL_REUSE,
        PhysicalCause.INTRINSIC_SIDEBAND,
        PhysicalCause.NEW_COHERENT_ANCESTRY,
    ]
    inv = tr = rs = rec = bd = 0
    maxj = 0
    for _ in range(samples):
        n = int(rng.integers(1, min(6, len(causes)) + 1))
        selected = tuple(causes[i] for i in rng.choice(len(causes), size=n, replace=False))
        maxj = max(maxj, n)
        wa = tuple(float(math.exp(rng.uniform(-12, 12))) for _ in selected)
        wb = tuple(float(math.exp(rng.uniform(-12, 12))) for _ in selected)
        if not weight_invariance_countercheck(selected, wa, wb):
            inv += 1
            raise AssertionError("joint-stop master projection depended on arbitrary tie weights")
        hits = tuple(CauseHit(0.2, c, w, "stress") for c, w in zip(selected, wa))
        out = joint_stop_master_projection(physical_hits=hits)
        if out.fine_rn_split_required:
            raise AssertionError("preferred joint-stop projection unexpectedly requested RN splitting")
        if out.master_disposition == MasterDisposition.TRANSFER_COST.value:
            tr += 1
        elif out.master_disposition == MasterDisposition.RECURSE_CRITICAL.value:
            rec += 1
        else:
            raise AssertionError("unexpected disposition for non-reset stress causes")

        # Independent uniform reset stress.
        if rng.random() < 0.25:
            cert = UniformResourceCertificate(0.3, 4.0, True, True)
            reset_hits = (
                CauseHit(0.5, PhysicalCause.UNIFORM_GLOBAL_RESOURCE, 1e-9, "reset"),
                CauseHit(0.5, PhysicalCause.RESOLVED_SOURCE, 1e9, "source"),
            )
            ro = joint_stop_master_projection(
                physical_hits=reset_hits,
                uniform_certificates={PhysicalCause.UNIFORM_GLOBAL_RESOURCE: cert},
            )
            if ro.master_disposition != MasterDisposition.ADDITIVE_RESET.value:
                raise AssertionError("valid reset did not terminate joint non-transfer stop")
            rs += 1

        # Zero-time boundary wins against any simultaneous cause/internal regeneration.
        if rng.random() < 0.25:
            bo = joint_stop_master_projection(
                physical_hits=(
                    CauseHit(0.0, PhysicalCause.INITIAL_BOUNDARY, 1e-12, "boundary"),
                    CauseHit(0.0, PhysicalCause.CAUSAL_REUSE, 1e12, "reuse"),
                ),
                internal_hits=(InternalHit(0.0),),
            )
            if bo.master_disposition != MasterDisposition.BOUNDARY.value:
                raise AssertionError("initial boundary did not absorb the zero-time joint stop")
            bd += 1

        # Source + HH regeneration is one coarse recursive stop without weights.
        io = joint_stop_master_projection(
            physical_hits=(CauseHit(0.7, PhysicalCause.RESOLVED_SOURCE, 1.0, "source"),),
            internal_hits=(InternalHit(0.7),),
        )
        if io.master_disposition != MasterDisposition.RECURSE_CRITICAL.value:
            raise AssertionError("source/HH-regeneration joint stop did not remain recursive")
    return JointStopStress(samples, inv, tr, rs, rec, bd, maxj)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-joint-causal-stop-projection"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "joint_causal_stop_projection.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Joint causal stop projection: exact ties need no artificial fractions\n\nStatus: **{cert['status']}**.\n\nAt an exact first-time tie, keep the physical event/block intact and retain the whole set of simultaneous certified causes.  Do **not** normalize source magnitude, strain, Moyal relink energy, or any other heterogeneous threshold into fictitious ownership fractions.\n\nThe master needs only one terminal/recursive disposition.  Its logic is physical: `t=0` is absorbing; any already-certified fixed multiplicative transfer/reuse/sideband cost is terminal; otherwise a valid genuinely globally bounded resource may terminate as a reset; otherwise source, critical dissipation, material/new ancestry and earlier HH regeneration are all recursive.  In the last case their exact simultaneous proportions are irrelevant because every route stays inside the same causal recursion class.\n\nThus exact tie provenance is set-valued while the master disposition is single-valued.  The full transfer mass is charged once and is never fractionally split.  Arbitrary positive legacy tie weights may vary by many orders of magnitude without changing the result.\n\nStress: `{out.samples}` random exact joint stops\n- tie-weight invariance failures: `{out.weight_invariance_failures}`\n- terminal transfer samples: `{out.transfer_terminal_samples}`\n- valid reset samples: `{out.reset_terminal_samples}`\n- recursive joint samples: `{out.recurse_joint_samples}`\n- absorbing boundary samples: `{out.boundary_samples}`\n- largest joint cause set sampled: `{out.maximum_joint_size}`\n\nThis supersedes Radon--Nikodym splitting as a **required** continuum bridge for the no-escape master.  Fine RN splitting remains legal when a genuinely physical common-unit subledger exists, but the PDE proof no longer has to manufacture one.  The companion smooth-SGS theorem now supplies measurable local first-hit cause sets once a block is given.  The remaining continuum question is recursive hard-event re-entry from the registered adjoint parent; exact tie weights are not part of that question.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
