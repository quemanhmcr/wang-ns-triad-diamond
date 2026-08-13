from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from src.material_label_carrier_quotient import (
    SELECTED_FAMILY_EVENT,
    carrier_registration_with_material_sidecars,
)
from src.physical_energy_causal_bridge import INHERIT_ENERGY_FRACTION, route_physical_energy_causality
from src.same_carrier_checkpoint_segmentation_quotient import (
    SameCarrierMonitorSegment,
    join_same_carrier_segments,
    partition_same_carrier_path,
    same_carrier_first_exit,
)

STATUS = (
    "EXACT_SAME_CARRIER_INHERITED_ENERGY_RELAY__"
    "PHYSICAL_STOCK_CONTINUATION_NOT_GENERATION__"
    "NO_EVENT_CHECKPOINT__MATERIAL_SIDECARS_PRESERVED_SEPARATELY__"
    "ROLE_OR_PROBE_CHANGE_FAILS_CLOSED__NO_TEMPORAL_DEPOSIT_MATCHING"
)

SAME_CARRIER_INHERITED_STOCK_RELAY = "same_carrier_inherited_energy_stock_relay"


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


def _tol(*values: float) -> float:
    return 2.0e-11 * max(1.0, *(abs(float(x)) for x in values))


@dataclass(frozen=True)
class SameCarrierInheritedEnergyRelayCertificate:
    carrier_id: str
    initial_time: float
    terminal_time: float
    initial_energy: float
    terminal_energy: float
    inherited_fraction: float
    residual_positive_work: float
    residual_owner_threshold: float
    observed_elapsed: float
    analysis_segments: int
    inserted_checkpoint_boundaries: int
    material_sidecars: tuple[str, ...] = ()
    selected_family_switch_energy: float = 0.0
    initial_endpoint_is_non_event_carrier_slice: bool = True
    relay_label: str = SAME_CARRIER_INHERITED_STOCK_RELAY
    recursive_generation_created: bool = False
    new_causal_charge_created: bool = False
    monitor_reset_used: bool = False
    later_hahn_used: bool = False
    between_time_deposit_matching_used: bool = False
    hard_cell_inventory_used: bool = False

    def __post_init__(self) -> None:
        if not self.carrier_id:
            raise ValueError("same-carrier inherited stock requires a carrier id")
        s = _finite_nonnegative(self.initial_time, "initial physical time")
        t = _finite_nonnegative(self.terminal_time, "terminal physical time")
        if not t > s:
            raise ValueError("inherited stock relay requires positive physical elapsed time")
        e0 = _finite_nonnegative(self.initial_energy, "initial carrier energy")
        e1 = _finite_nonnegative(self.terminal_energy, "terminal carrier energy")
        if e0 <= 0.0 or e1 <= 0.0:
            raise ValueError("positive carrier stock required at both relay endpoints")
        frac = _finite_nonnegative(self.inherited_fraction, "inherited energy fraction")
        if abs(frac - e0 / e1) > _tol(frac, e0 / e1):
            raise AssertionError("stored inherited fraction changed from physical carrier energies")
        inherit_threshold = float(INHERIT_ENERGY_FRACTION) * e1
        if e0 + _tol(e0, e1) < inherit_threshold:
            raise AssertionError("same-carrier relay does not satisfy the physical inherited-energy gate")
        wr = _finite_nonnegative(self.residual_positive_work, "classified residual positive work")
        wr_face = _finite_nonnegative(self.residual_owner_threshold, "classified residual owner threshold")
        if abs(wr_face - inherit_threshold) > _tol(wr_face, inherit_threshold):
            raise AssertionError("residual owner threshold changed from E1/5")
        # The stock relay is a one-component quotient, not a lexicographic erasure
        # of simultaneous physical work.  If classified residual work reaches its
        # own E1/5 face, preserve that event-facing owner instead of declaring the
        # whole block to be only inherited stock.  Near ties fail closed.
        if wr + _tol(wr, wr_face, e1) >= wr_face:
            raise TypeError("classified residual physical work reaches its owner face; stock-only quotient is inapplicable")
        elapsed = _finite_nonnegative(self.observed_elapsed, "observed elapsed time")
        if abs(elapsed - (t - s)) > _tol(elapsed, t - s):
            raise AssertionError("same-carrier monitor interval changed physical elapsed time")
        if self.analysis_segments < 1 or self.inserted_checkpoint_boundaries != self.analysis_segments - 1:
            raise ValueError("checkpoint count changed from one fixed-carrier path segmentation")
        if tuple(sorted(set(self.material_sidecars))) != self.material_sidecars:
            raise ValueError("material sidecars must be a sorted quotiented set")
        switch_energy = _finite_nonnegative(self.selected_family_switch_energy, "selected-family Moyal switch energy")
        has_family_switch = SELECTED_FAMILY_EVENT in self.material_sidecars
        if not has_family_switch and switch_energy > 5.0e-15 * max(e0, e1, 1.0e-300):
            raise ValueError("selected-family switch energy requires the selected-family sidecar label")
        if not self.initial_endpoint_is_non_event_carrier_slice:
            raise TypeError("a genuine physical event at the earlier endpoint cannot be erased into inherited stock")
        if self.relay_label != SAME_CARRIER_INHERITED_STOCK_RELAY:
            raise ValueError("same-carrier inherited relay label changed")
        if (
            self.recursive_generation_created
            or self.new_causal_charge_created
            or self.monitor_reset_used
            or self.later_hahn_used
            or self.between_time_deposit_matching_used
            or self.hard_cell_inventory_used
        ):
            raise ValueError(
                "inherited stock continuity may not create event depth/charge, reset monitors, re-Hahn, temporally match deposits, or use hard cells as wallets"
            )


def same_carrier_inherited_energy_relay(
    segments: Sequence[SameCarrierMonitorSegment],
    *,
    initial_time: float,
    terminal_time: float,
    initial_energy: float,
    terminal_energy: float,
    residual_positive_work: float,
    strain_action: float,
    material_registration: Mapping[str, object] | None = None,
    initial_endpoint_is_non_event_carrier_slice: bool,
) -> SameCarrierInheritedEnergyRelayCertificate:
    """Quotient a same-carrier inheritance branch into between-time stock continuity.

    This theorem is intentionally narrow.  The fixed carrier must have no named
    first stop on the observed interval, cumulative monitors must be the same
    event-anchored paths across all inserted checkpoints, and the *actual physical
    energy gate* must select ``material_energy_inheritance``.  Optional material
    membership/family changes are allowed only when the material-label quotient
    certifies that the same smooth role and analysis probe continue.

    The conclusion is not that inherited energy is free.  It is persistent
    physical stock at the earlier endpoint of the same carrier.  The stock
    component adds zero recursive generation depth and performs no FIFO/LIFO
    matching of earlier deposits to later withdrawals.  Any material sidecars
    remain explicit and require their existing ancestry/relink routing.
    """
    if not bool(initial_endpoint_is_non_event_carrier_slice):
        raise TypeError("same-carrier stock relay requires a certified non-event earlier endpoint")
    path = join_same_carrier_segments(segments)
    exit_record = same_carrier_first_exit(segments)
    if exit_record["classification"] != "same_carrier_no_hit_continuation":
        raise TypeError("same-carrier inherited stock quotient requires a no-physical-stop interval")
    if tuple(exit_record["joint_first_stops"]):
        raise AssertionError("no-hit inheritance interval unexpectedly carries a first-stop label")
    if bool(exit_record["requires_physical_energy_reentry"]):
        raise AssertionError("coefficient obstruction cannot be quotiented as inherited stock")
    if int(exit_record["carrier_restarts"]) or int(exit_record["monitor_resets"]):
        raise TypeError("same-carrier inherited stock cannot cross a carrier or monitor reset")

    K = _finite_nonnegative(strain_action, "strain action")
    path_K = float(tuple(path["strain_action"])[-1])
    if abs(K - path_K) > _tol(K, path_K):
        raise TypeError("energy-gate strain action must be the same cumulative carrier monitor")

    e1_for_face = _finite_nonnegative(terminal_energy, "terminal carrier energy")
    wr = _finite_nonnegative(residual_positive_work, "classified residual positive work")
    residual_face = float(INHERIT_ENERGY_FRACTION) * e1_for_face
    if wr + _tol(wr, residual_face, e1_for_face) >= residual_face:
        raise TypeError("classified residual physical work is a simultaneous owner; do not quotient the block to stock only")
    gate = route_physical_energy_causality(
        terminal_energy=e1_for_face,
        initial_energy=float(initial_energy),
        residual_positive_work=wr,
        strain_action=K,
    )
    if gate.get("branch") != "material_energy_inheritance":
        raise TypeError("physical energy gate did not select inherited carrier stock")
    if abs(float(gate["value"]) - float(initial_energy)) > _tol(float(gate["value"]), float(initial_energy)):
        raise AssertionError("inheritance gate value changed from initial carrier energy")

    sidecars: tuple[str, ...] = ()
    switch_energy = 0.0
    if material_registration is not None:
        if not bool(material_registration.get("quotient_applicable", False)):
            raise TypeError("role/probe-changing material registration cannot be quotiented as same-carrier stock")
        if not bool(material_registration.get("carrier_continuation_certified", False)):
            raise TypeError("material registration did not certify same-carrier continuation")
        if tuple(material_registration.get("carrier_first_stops", ())):
            raise TypeError("material registration contains a physical carrier stop")
        if not bool(material_registration.get("same_carrier_reusable_after_sidecar", False)):
            raise TypeError("material sidecar did not preserve the same PDE carrier")
        sidecars = tuple(sorted(set(str(x) for x in material_registration.get("sidecar_events", ()) if str(x))))
        switch_energy = _finite_nonnegative(
            float(material_registration.get("selected_family_switch_energy", 0.0)),
            "registered selected-family switch energy",
        )

    s = _finite_nonnegative(initial_time, "initial physical time")
    t = _finite_nonnegative(terminal_time, "terminal physical time")
    observed = float(exit_record["observed_elapsed_end"])
    if abs(observed - (t - s)) > _tol(observed, t - s):
        raise TypeError("same-carrier monitor path does not cover the claimed physical interval")

    e0 = _finite_nonnegative(initial_energy, "initial carrier energy")
    e1 = _finite_nonnegative(terminal_energy, "terminal carrier energy")
    return SameCarrierInheritedEnergyRelayCertificate(
        carrier_id=str(exit_record["carrier_id"]),
        initial_time=s,
        terminal_time=t,
        initial_energy=e0,
        terminal_energy=e1,
        inherited_fraction=e0 / e1,
        residual_positive_work=wr,
        residual_owner_threshold=float(INHERIT_ENERGY_FRACTION) * e1,
        observed_elapsed=observed,
        analysis_segments=int(exit_record["analysis_segments"]),
        inserted_checkpoint_boundaries=int(exit_record["inserted_checkpoint_boundaries"]),
        material_sidecars=sidecars,
        selected_family_switch_energy=switch_energy,
        initial_endpoint_is_non_event_carrier_slice=True,
    )


@dataclass(frozen=True)
class SameCarrierInheritanceMasterProjection:
    physical_measure: str
    mass: float
    certificate: SameCarrierInheritedEnergyRelayCertificate
    sidecar_events: tuple[str, ...]
    between_time_stock_relays: tuple[str, ...] = (SAME_CARRIER_INHERITED_STOCK_RELAY,)
    stock_owner_bundle_created: bool = False
    stock_recursive_event_created: bool = False
    sidecars_quotiented_as_stock: bool = False

    def __post_init__(self) -> None:
        if not self.physical_measure or not math.isfinite(float(self.mass)) or float(self.mass) < 0.0:
            raise ValueError("valid physical stock measure and mass required")
        if self.between_time_stock_relays != (SAME_CARRIER_INHERITED_STOCK_RELAY,):
            raise ValueError("same-carrier inheritance master projection lost its stock-relay type")
        if self.sidecar_events != self.certificate.material_sidecars:
            raise ValueError("master projection must preserve every material sidecar separately from stock")
        if self.stock_owner_bundle_created or self.stock_recursive_event_created:
            raise ValueError("between-time inherited stock may not be projected as a generation owner")
        if self.sidecars_quotiented_as_stock:
            raise ValueError("material ancestry/relink sidecars may not be erased by the stock quotient")


def same_carrier_inheritance_master_projection(
    physical_measure: str,
    mass: float,
    reentry: Mapping[str, object],
    certificate: SameCarrierInheritedEnergyRelayCertificate,
) -> SameCarrierInheritanceMasterProjection:
    """Typed master-facing projection, before wiring into the central master.

    The reentry must be the exact physical energy-gate inheritance law from which
    ``certificate`` was constructed.  This function intentionally lives outside
    the central master until its theorem is independently certified; that keeps
    theorem proof and master wiring as separate reviewable commits.
    """
    if str(reentry.get("branch", "")) != "material_energy_inheritance":
        raise TypeError("same-carrier stock projection requires the material-energy inheritance gate")
    if bool(reentry.get("coefficient_impulse_used_as_physical_work", False)):
        raise TypeError("coefficient impulse magnitude may not be used as inherited physical stock")
    if bool(reentry.get("observer_partition_motion_charged_as_physics", False)):
        raise TypeError("observer partition motion may not be charged as inherited physical stock")
    value = float(reentry.get("value", math.nan))
    threshold = float(reentry.get("threshold", math.nan))
    residual = float(reentry.get("classified_residual_positive_work", math.nan))
    if not all(math.isfinite(v) for v in (value, threshold, residual)):
        raise TypeError("same-carrier stock projection requires finite gate value, threshold and actual classified residual work")
    tol = _tol(value, threshold, residual, certificate.initial_energy, certificate.terminal_energy)
    if abs(value - certificate.initial_energy) > tol:
        raise TypeError("relay certificate does not belong to this physical inheritance-gate value")
    expected = float(INHERIT_ENERGY_FRACTION) * certificate.terminal_energy
    if abs(threshold - expected) > tol:
        raise TypeError("relay certificate does not belong to this physical inheritance-gate threshold")
    if abs(residual - certificate.residual_positive_work) > tol:
        raise TypeError("relay certificate does not belong to this classified residual physical-work law")
    return SameCarrierInheritanceMasterProjection(
        physical_measure=str(physical_measure),
        mass=float(mass),
        certificate=certificate,
        sidecar_events=certificate.material_sidecars,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "physical_gate": "on one fixed low-strain carrier, E1<=exp(2K)(E0+W_HH+W_R); the branch E0>=E1/5 is earlier physical carrier stock, not a new positive work law",
        "same_carrier_requirement": "Q/psi identity and cumulative carrier-first-hit monitors must survive all analysis checkpoints; any role/probe change or named carrier first stop fails closed",
        "material_sidecars": "intrinsic membership and selected-family sidecars may be reread only through the certified same-Q/same-psi quotient; they do not create a second coefficient impulse or carrier, but their ancestry/relink currency is preserved for separate routing",
        "master_ontology": "the same-carrier inherited-stock component is a between-time relay with zero generation depth only when classified residual work stays strictly below its E1/5 owner face; material sidecars remain separately routed, and legacy/untyped inheritance remains event-facing until this certificate is supplied",
        "temporal_provenance": "the theorem uses endpoint stock continuity only and performs no FIFO/LIFO/proportional matching of earlier deposits to later withdrawals",
        "event_boundary": "the earlier endpoint must be certified as a non-event slice of the continuing carrier; if it is a new physical interaction, or the carrier role/probe changes, this quotient is inapplicable and the existing physical-event route must be used",
        "later_hahn_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class SameCarrierInheritanceStress:
    samples: int
    relays_checked: int
    sidecar_relays: int
    maximum_checkpoint_count: int
    minimum_inheritance_margin: float
    maximum_elapsed_residual: float
    stock_recursive_events_created: int


def stress(samples: int = 50_000, seed: int = 2026081304) -> SameCarrierInheritanceStress:
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = random.Random(int(seed))
    checked = sidecars = max_checkpoints = recursive = 0
    min_margin = math.inf
    max_elapsed = 0.0
    for j in range(count):
        horizon = 10.0 ** rng.uniform(-5.0, -1.0)
        n = rng.randint(5, 13)
        elapsed = tuple(horizon * i / (n - 1) for i in range(n))
        K_end = rng.uniform(0.0, 0.92 / 30.0)
        strain_power = rng.uniform(0.7, 1.8)
        strain = tuple(K_end * (i / (n - 1)) ** strain_power for i in range(n))
        amp = 10.0 ** rng.uniform(-2.0, 2.0)
        residual = tuple((0.20 * amp) * (i / (n - 1)) * rng.uniform(0.0, 0.95) for i in range(n))
        hh = tuple((0.42 * amp) * (i / (n - 1)) * rng.uniform(0.0, 0.95) for i in range(n))
        # Impulse magnitudes need not be monotone physically.  Keep this stress on
        # the no-hit side while allowing some cancellation after an interior peak.
        if n >= 7 and j % 3 == 0:
            residual = list(residual)
            residual[-2] = min(0.22 * amp, residual[-2] + 0.04 * amp)
            residual[-1] = min(residual[-1], residual[-2] * 0.8)
            residual = tuple(residual)
        interior = list(range(1, n - 1))
        rng.shuffle(interior)
        cuts = tuple(sorted(interior[: rng.randint(0, min(3, len(interior)))]))
        segments = partition_same_carrier_path(
            carrier_id=f"carrier-{j}",
            terminal_amplitude=amp,
            elapsed_times=elapsed,
            strain_action=strain,
            residual_impulse_abs=residual,
            hh_impulse_abs=hh,
            checkpoint_indices=cuts,
        )
        e1 = 10.0 ** rng.uniform(-4.0, 2.0)
        ratio = rng.uniform(0.205, 1.8)
        e0 = ratio * e1
        wR = rng.uniform(0.0, 0.19) * e1
        use_sidecar = (j % 2 == 0)
        material = None
        if use_sidecar:
            i_hh = complex(rng.uniform(-0.10, 0.10) * amp, rng.uniform(-0.10, 0.10) * amp)
            i_if = complex(rng.uniform(-0.06, 0.06) * amp, rng.uniform(-0.06, 0.06) * amp)
            z_event = complex(amp, 0.0)
            z_slice = z_event - i_hh - i_if
            material = carrier_registration_with_material_sidecars(
                z_event,
                z_slice,
                i_hh,
                i_if,
                intrinsic_material_membership_change=True,
                selected_family_change=(j % 4 == 0),
                selected_family_switch_energy=(0.01 if j % 4 == 0 else 0.0),
                same_smooth_role=True,
                same_analysis_probe=True,
            )
        cert = same_carrier_inherited_energy_relay(
            segments,
            initial_time=3.0,
            terminal_time=3.0 + horizon,
            initial_energy=e0,
            terminal_energy=e1,
            residual_positive_work=wR,
            strain_action=K_end,
            material_registration=material,
            initial_endpoint_is_non_event_carrier_slice=True,
        )
        checked += 1
        sidecars += int(bool(cert.material_sidecars))
        max_checkpoints = max(max_checkpoints, cert.inserted_checkpoint_boundaries)
        min_margin = min(min_margin, e0 - float(INHERIT_ENERGY_FRACTION) * e1)
        max_elapsed = max(max_elapsed, abs(cert.observed_elapsed - horizon))
        recursive += int(cert.recursive_generation_created)
    if checked != count or recursive:
        raise AssertionError("same-carrier inheritance stress changed stock relay into recursive generation")
    return SameCarrierInheritanceStress(
        samples=count,
        relays_checked=checked,
        sidecar_relays=sidecars,
        maximum_checkpoint_count=max_checkpoints,
        minimum_inheritance_margin=min_margin,
        maximum_elapsed_residual=max_elapsed,
        stock_recursive_events_created=recursive,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081304)
    ap.add_argument("--outdir", type=Path, default=Path("results-same-carrier-inherited-energy-relay"))
    args = ap.parse_args()
    out = stress(args.samples, args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)
    payload = {"certificate": theorem_certificate(), "stress": asdict(out)}
    (args.outdir / "same_carrier_inherited_energy_relay.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = f"""# Same-carrier inherited-energy relay

Status: **{STATUS}**.

On a fixed no-hit carrier interval, if `E0>=E1/5` while classified residual physical work stays strictly below its own `E1/5` owner face, the earlier endpoint is carrier stock rather than a new generation charge.  Checkpoint cuts leave cumulative monitors unchanged.  Same-Q/same-psi material sidecars may be reread without replacing the carrier, but remain separately routed ancestry/relink currency.  Only the stock component has zero recursive generation depth; no temporal deposit matching is performed.

Stress: `{out.samples}` exact same-carrier paths
- relays checked: `{out.relays_checked}`
- stock relays carrying material sidecars preserved for separate routing: `{out.sidecar_relays}`
- maximum inserted checkpoints: `{out.maximum_checkpoint_count}`
- minimum inheritance-gate margin: `{out.minimum_inheritance_margin:.3e}`
- maximum physical elapsed-time residual: `{out.maximum_elapsed_residual:.3e}`
- recursive events created by the stock component: `{out.stock_recursive_events_created}`

Untyped inheritance, role/probe change, or any named physical first stop remains on the existing event route.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")


if __name__ == "__main__":
    main()
