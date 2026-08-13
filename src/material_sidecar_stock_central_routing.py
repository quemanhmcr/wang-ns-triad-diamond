from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Mapping

from src.joint_causal_stop_projection import InternalHit, JointStopProjection, joint_stop_master_projection
from src.physical_branch_compiler import CauseHit, MasterDisposition, PhysicalCause, UniformResourceCertificate
from src.material_label_carrier_quotient import MATERIAL_MEMBERSHIP_EVENT, SELECTED_FAMILY_EVENT
from src.material_sidecar_stock_owner_decomposition import (
    MEMBERSHIP_PROVENANCE_CURRENCY,
    SELECTED_FAMILY_MOYAL_CURRENCY,
    MaterialSidecarStockDecomposition,
    material_sidecar_stock_decomposition,
)
from src.same_carrier_inherited_energy_relay import (
    SAME_CARRIER_INHERITED_STOCK_RELAY,
    SameCarrierInheritedEnergyRelayCertificate,
)


STATUS = (
    "EXACT_MATERIAL_SIDECAR_STOCK_CENTRAL_RELAY__"
    "ONE_INHERITED_STOCK_CHARGE__NON_EVENT_MOYAL_BOUNDARY_PRESERVED__"
    "GENUINE_MATERIAL_SERVICE_REMAINS_SEPARATE_PHYSICAL_CAUSE"
)


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


def _relative_match(a: float, b: float, factor: float = 5.0e-12) -> bool:
    aa = float(a)
    bb = float(b)
    return abs(aa - bb) <= factor * max(abs(aa), abs(bb), 1.0e-300)


@dataclass(frozen=True)
class MaterialSidecarStockCentralRelayCertificate:
    """Non-event material sidecars carried next to one inherited stock charge.

    This object is deliberately not a physical CauseHit.  It preserves the
    observation/selected-service boundary data that survive the stock quotient,
    while certifying that those data alone mint neither work nor event depth.
    """

    carrier_id: str
    inherited_stock_mass: float
    stock_relay_label: str
    sidecar_events: tuple[str, ...]
    sidecar_currencies: tuple[str, ...]
    selected_family_boundary_energy: float
    decomposition: MaterialSidecarStockDecomposition
    stock_charge_count: int = 1
    recursive_generation_created: bool = False
    physical_work_created: bool = False
    sidecar_promoted_to_physical_hit: bool = False
    smooth_k_phys_identification_used: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        stock = _finite_nonnegative(self.inherited_stock_mass, "inherited carrier stock")
        boundary = _finite_nonnegative(self.selected_family_boundary_energy, "selected-family boundary energy")
        if not self.carrier_id or stock <= 0.0:
            raise ValueError("positive inherited stock on a named physical carrier required")
        if self.stock_relay_label != SAME_CARRIER_INHERITED_STOCK_RELAY or self.stock_charge_count != 1:
            raise ValueError("central material sidecar relay must preserve exactly one inherited stock charge")
        if tuple(sorted(set(self.sidecar_events))) != self.sidecar_events or not self.sidecar_events:
            raise ValueError("central sidecar event set must be a nonempty sorted quotient")
        if tuple(sorted(set(self.sidecar_currencies))) != self.sidecar_currencies:
            raise ValueError("central sidecar currency set must be sorted and quotiented")
        dec = self.decomposition
        if dec.carrier_id != self.carrier_id or not _relative_match(dec.inherited_stock_mass, stock):
            raise ValueError("central sidecar relay changed carrier identity or inherited stock mass")
        if dec.stock_relay_label != self.stock_relay_label or dec.stock_charge_count != 1:
            raise ValueError("central sidecar relay changed the one-stock-charge decomposition")
        events = tuple(c.event for c in dec.charges)
        currencies = tuple(sorted(c.currency for c in dec.charges))
        if events != self.sidecar_events or currencies != self.sidecar_currencies:
            raise ValueError("central sidecar relay changed the typed material events/currencies")
        if not _relative_match(dec.selected_family_switch_energy, boundary):
            raise ValueError("central sidecar relay changed the selected-family Moyal boundary energy")
        if SELECTED_FAMILY_EVENT in self.sidecar_events:
            if not dec.family_switch_moyal_certificate_bound:
                raise ValueError("selected-family central relay lacks the exact Phase-A Moyal binding")
        elif boundary != 0.0 or dec.family_switch_moyal_certificate_bound:
            raise ValueError("selected-family boundary energy exists without selected-family sidecar")
        membership = tuple(c for c in dec.charges if c.event == MATERIAL_MEMBERSHIP_EVENT)
        if any(c.currency != MEMBERSHIP_PROVENANCE_CURRENCY or c.charge != 0.0 for c in membership):
            raise ValueError("membership rereading changed from zero-charge provenance")
        family = tuple(c for c in dec.charges if c.event == SELECTED_FAMILY_EVENT)
        if any(c.currency != SELECTED_FAMILY_MOYAL_CURRENCY or not _relative_match(c.charge, boundary) for c in family):
            raise ValueError("selected-family sidecar changed from exact Moyal boundary currency")
        if (
            self.recursive_generation_created
            or self.physical_work_created
            or self.sidecar_promoted_to_physical_hit
            or self.smooth_k_phys_identification_used
            or self.later_hahn_used
        ):
            raise ValueError("non-event material boundary sidecar was promoted into physical generation/work/relink")


def material_sidecar_stock_central_relay(
    inherited_certificate: SameCarrierInheritedEnergyRelayCertificate,
    decomposition: MaterialSidecarStockDecomposition,
) -> MaterialSidecarStockCentralRelayCertificate:
    """Bind a Phase-A sidecar decomposition to the exact inherited-stock carrier."""
    if not isinstance(inherited_certificate, SameCarrierInheritedEnergyRelayCertificate):
        raise TypeError("typed same-carrier inherited-stock certificate required")
    if not isinstance(decomposition, MaterialSidecarStockDecomposition):
        raise TypeError("typed Phase-A material-sidecar stock decomposition required")
    if not inherited_certificate.material_sidecars:
        raise TypeError("material-sidecar central relay is only for sidecar-bearing inherited stock")
    if decomposition.carrier_id != inherited_certificate.carrier_id:
        raise TypeError("material-sidecar decomposition was transplanted to a different carrier")
    if not _relative_match(decomposition.inherited_stock_mass, inherited_certificate.initial_energy):
        raise TypeError("material-sidecar decomposition changed the inherited E0 stock mass")
    events = tuple(c.event for c in decomposition.charges)
    if events != inherited_certificate.material_sidecars:
        raise TypeError("material-sidecar decomposition does not belong to this inherited certificate")
    if not _relative_match(
        decomposition.selected_family_switch_energy,
        inherited_certificate.selected_family_switch_energy,
    ):
        raise TypeError("selected-family Moyal boundary charge was transplanted or changed")
    if decomposition.stock_cloned_per_sidecar or decomposition.sidecar_charge_added_to_stock:
        raise TypeError("material sidecars may not clone or augment inherited carrier stock")
    if decomposition.sidecar_charge_added_to_physical_work or decomposition.k_phys_identification_used:
        raise TypeError("material boundary sidecar may not be promoted to dW or smooth K_phys")
    return MaterialSidecarStockCentralRelayCertificate(
        carrier_id=inherited_certificate.carrier_id,
        inherited_stock_mass=inherited_certificate.initial_energy,
        stock_relay_label=SAME_CARRIER_INHERITED_STOCK_RELAY,
        sidecar_events=events,
        sidecar_currencies=tuple(sorted(c.currency for c in decomposition.charges)),
        selected_family_boundary_energy=decomposition.selected_family_switch_energy,
        decomposition=decomposition,
    )


@dataclass(frozen=True)
class MaterialSidecarJointStopProjection:
    """Certified physical joint stop with non-event material boundary data beside it."""

    physical_projection: JointStopProjection
    relay_certificate: MaterialSidecarStockCentralRelayCertificate
    non_event_material_sidecar_events: tuple[str, ...]
    non_event_material_sidecar_currencies: tuple[str, ...]
    selected_family_boundary_energy: float

    def __post_init__(self) -> None:
        relay = self.relay_certificate
        if self.non_event_material_sidecar_events != relay.sidecar_events:
            raise ValueError("joint wrapper changed material sidecar provenance")
        if self.non_event_material_sidecar_currencies != relay.sidecar_currencies:
            raise ValueError("joint wrapper changed material sidecar currency")
        if not _relative_match(self.selected_family_boundary_energy, relay.selected_family_boundary_energy):
            raise ValueError("joint wrapper changed selected-family Moyal boundary energy")
        if relay.recursive_generation_created or relay.sidecar_promoted_to_physical_hit:
            raise ValueError("joint wrapper received a material boundary already promoted into a physical event")

    @property
    def first_time(self) -> float | None:
        return self.physical_projection.first_time

    @property
    def joint_physical_causes(self) -> tuple[str, ...]:
        return self.physical_projection.joint_physical_causes

    @property
    def joint_internal_causes(self) -> tuple[str, ...]:
        return self.physical_projection.joint_internal_causes

    @property
    def certified_currencies(self) -> tuple[str, ...]:
        return self.physical_projection.certified_currencies

    @property
    def master_disposition(self) -> str:
        return self.physical_projection.master_disposition

    @property
    def terminal_certificate_used(self) -> str | None:
        return self.physical_projection.terminal_certificate_used

    @property
    def fine_rn_split_required(self) -> bool:
        return self.physical_projection.fine_rn_split_required


def material_sidecar_joint_stop_projection(
    relay_certificate: MaterialSidecarStockCentralRelayCertificate,
    *,
    physical_hits: tuple[CauseHit, ...] = (),
    internal_hits: tuple[InternalHit, ...] = (),
    fixed_transfer_loss: bool = False,
    kelvin_flat_certified: bool = False,
    uniform_certificates: Mapping[PhysicalCause, UniformResourceCertificate] | None = None,
) -> MaterialSidecarJointStopProjection:
    """Run the certified joint-stop law without ever inserting sidecars as hits.

    The decisive point is structural: ``relay_certificate`` is not converted to a
    ``CauseHit`` and is not passed to the first-time minimizer.  The physical joint
    stop is computed first from independently witnessed PDE-facing causes; only
    then is the already-typed material boundary data carried beside the result.
    """
    if not isinstance(relay_certificate, MaterialSidecarStockCentralRelayCertificate):
        raise TypeError("typed non-event material-sidecar relay required")
    physical = joint_stop_master_projection(
        physical_hits=physical_hits,
        internal_hits=internal_hits,
        fixed_transfer_loss=fixed_transfer_loss,
        kelvin_flat_certified=kelvin_flat_certified,
        uniform_certificates=uniform_certificates,
    )
    return MaterialSidecarJointStopProjection(
        physical_projection=physical,
        relay_certificate=relay_certificate,
        non_event_material_sidecar_events=relay_certificate.sidecar_events,
        non_event_material_sidecar_currencies=relay_certificate.sidecar_currencies,
        selected_family_boundary_energy=relay_certificate.selected_family_boundary_energy,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "stock": "the physical inherited E0 remains exactly one between-time carrier-stock charge",
        "membership": "intrinsic membership rereading is retained beside stock as zero-charge provenance and never becomes a physical CauseHit",
        "selected_family": "the exact Phase-A R_switch is retained beside stock as non-event Moyal selected-service boundary currency; positive R_switch alone cannot create a first stop or recursive child",
        "genuine_material_service": "an independently certified physical material/source service event is a different object and remains eligible for the existing recursive physical-cause route",
        "joint_stop": "the wrapper never converts a material sidecar to CauseHit: the certified core joint-stop law runs first on independently witnessed physical/internal hits, then the sidecar is preserved beside that result; therefore it cannot alter first time, tie ownership, or disposition",
        "currency_separation": "no stock cloning, no dW addition, no smooth K_phys identification, no later Hahn split",
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class MaterialSidecarCentralStress:
    samples: int
    membership_cases: int
    selected_family_cases: int
    mixed_cases: int
    positive_boundary_cases: int
    transplant_rejections: int
    recursive_generation_events_created: int
    physical_hits_promoted: int
    stock_clone_violations: int
    joint_sidecar_projections: int
    sidecar_only_stop_rejections: int
    genuine_material_service_preserved: int
    maximum_boundary_energy: float


def stress(samples: int = 50_000, seed: int = 2026081306) -> MaterialSidecarCentralStress:
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = random.Random(int(seed))
    membership = family = mixed = positive = rejected = recursive = promoted = clone = 0
    joint_count = sidecar_only_rejected = genuine_material_preserved = 0
    max_boundary = 0.0
    for j in range(count):
        mode = j % 3
        membership_change = mode in (0, 2)
        family_change = mode in (1, 2)
        sidecars = tuple(sorted(
            ([MATERIAL_MEMBERSHIP_EVENT] if membership_change else [])
            + ([SELECTED_FAMILY_EVENT] if family_change else [])
        ))
        e0 = 10.0 ** rng.uniform(-8.0, 4.0)
        e1 = e0 * rng.uniform(0.55, 4.8)  # always e0 >= e1/5 with room
        threshold = 0.2 * e1
        residual = threshold * rng.uniform(0.0, 0.70)
        boundary = (10.0 ** rng.uniform(-18.0, 3.0)) if family_change else 0.0
        inherited = SameCarrierInheritedEnergyRelayCertificate(
            carrier_id=f"carrier-{j}",
            initial_time=1.0,
            terminal_time=1.25,
            initial_energy=e0,
            terminal_energy=e1,
            inherited_fraction=e0 / e1,
            residual_positive_work=residual,
            residual_owner_threshold=threshold,
            observed_elapsed=0.25,
            analysis_segments=1,
            inserted_checkpoint_boundaries=0,
            material_sidecars=sidecars,
            selected_family_switch_energy=boundary,
        )
        switch = None
        if family_change:
            jump = boundary * rng.uniform(-0.95, 0.95)
            switch = {
                "selected_family_changed": True,
                "symmetric_difference_energy": boundary,
                "selection_energy_jump": jump,
                "jump_bound_margin": boundary - abs(jump),
            }
        dec = material_sidecar_stock_decomposition(
            inherited,
            selected_family_switch_certificate=switch,
        )
        out = material_sidecar_stock_central_relay(inherited, dec)
        membership += int(membership_change)
        family += int(family_change)
        mixed += int(membership_change and family_change)
        positive += int(out.selected_family_boundary_energy > 0.0)
        max_boundary = max(max_boundary, out.selected_family_boundary_energy)
        recursive += int(out.recursive_generation_created)
        promoted += int(out.sidecar_promoted_to_physical_hit)
        clone += int(out.stock_charge_count != 1 or out.decomposition.stock_cloned_per_sidecar)

        # The sidecar is carried beside a real stop, never injected into its cause set.
        cause = (
            PhysicalCause.MATERIAL_RELINK
            if j % 5 == 0
            else (PhysicalCause.HIGH_STRAIN_DISSIPATION if j % 5 == 1 else PhysicalCause.RESOLVED_SOURCE)
        )
        physical_hits = (CauseHit(0.4, cause, 10.0 ** rng.uniform(-9.0, 9.0), "independent physical witness"),)
        joint = material_sidecar_joint_stop_projection(out, physical_hits=physical_hits)
        joint_count += 1
        if joint.joint_physical_causes != (cause.value,) or joint.master_disposition != MasterDisposition.RECURSE_CRITICAL.value:
            raise AssertionError("material boundary sidecar changed the independently witnessed physical joint stop")
        if cause is PhysicalCause.MATERIAL_RELINK:
            genuine_material_preserved += int(joint.joint_physical_causes == (PhysicalCause.MATERIAL_RELINK.value,))

        # A positive/zero boundary sidecar by itself cannot manufacture a stop.
        if j % 17 == 0:
            try:
                material_sidecar_joint_stop_projection(out)
            except ValueError:
                sidecar_only_rejected += 1
            else:
                raise AssertionError("non-event material sidecar manufactured a causal first stop")

        if j % 11 == 0:
            foreign = replace(dec, carrier_id=f"foreign-{j}")
            try:
                material_sidecar_stock_central_relay(inherited, foreign)
            except TypeError:
                rejected += 1
            else:
                raise AssertionError("transplanted material-sidecar decomposition was accepted")

    if recursive or promoted or clone:
        raise AssertionError("central material-sidecar relay manufactured physical event depth or stock multiplicity")
    return MaterialSidecarCentralStress(
        samples=count,
        membership_cases=membership,
        selected_family_cases=family,
        mixed_cases=mixed,
        positive_boundary_cases=positive,
        transplant_rejections=rejected,
        recursive_generation_events_created=recursive,
        physical_hits_promoted=promoted,
        stock_clone_violations=clone,
        joint_sidecar_projections=joint_count,
        sidecar_only_stop_rejections=sidecar_only_rejected,
        genuine_material_service_preserved=genuine_material_preserved,
        maximum_boundary_energy=max_boundary,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081306)
    ap.add_argument("--outdir", type=Path, default=Path("results-material-sidecar-stock-central-routing"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples, args.seed)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "material_sidecar_stock_central_routing.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    md = f"""# Material-sidecar inherited-stock central routing\n\nStatus: **{STATUS}**.\n\nOne inherited carrier stock charge survives. Membership rereading is zero-charge provenance. Selected-family `R_switch` survives as non-event Moyal boundary currency and is never promoted into `dW`, smooth `K_phys`, or a physical first-stop hit. Independently evidenced material/source service remains a separate physical cause.\n\nStress: `{out.samples}` typed sidecar-bearing stock relays\n- membership / selected-family / mixed: `{out.membership_cases}` / `{out.selected_family_cases}` / `{out.mixed_cases}`\n- positive Moyal-boundary cases: `{out.positive_boundary_cases}`\n- transplanted-decomposition rejections: `{out.transplant_rejections}`\n- recursive-generation events created: `{out.recursive_generation_events_created}`\n- sidecars promoted to physical hits: `{out.physical_hits_promoted}`\n- stock-clone violations: `{out.stock_clone_violations}`\n- joint stops carrying sidecars after physical classification: `{out.joint_sidecar_projections}`\n- sidecar-only no-stop rejections: `{out.sidecar_only_stop_rejections}`\n- genuine material-service hits preserved: `{out.genuine_material_service_preserved}`\n- maximum sampled boundary energy: `{out.maximum_boundary_energy:.12e}`\n\nNo global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
