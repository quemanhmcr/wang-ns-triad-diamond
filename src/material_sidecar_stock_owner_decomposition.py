from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

from src.coherent_transfer_cells import coherent_ledger, service_no_escape
from src.material_label_carrier_quotient import (
    MATERIAL_MEMBERSHIP_EVENT,
    SELECTED_FAMILY_EVENT,
    selected_family_switch_sidecar,
)
from src.same_carrier_inherited_energy_relay import (
    SAME_CARRIER_INHERITED_STOCK_RELAY,
    SameCarrierInheritedEnergyRelayCertificate,
)

STATUS = (
    "EXACT_MATERIAL_SIDECAR_CHARGE_PRESERVATION_ON_INHERITED_STOCK__"
    "MEMBERSHIP_ZERO_CHARGE_PROVENANCE__SELECTED_FAMILY_EXACT_MOYAL_ENERGY__"
    "STOCK_NOT_CLONED__NO_DW_OR_KPHYS_IDENTIFICATION__SEPARATE_ANCESTRY_SERVICE"
)

MEMBERSHIP_PROVENANCE_CURRENCY = "zero_charge_material_membership_provenance"
SELECTED_FAMILY_MOYAL_CURRENCY = "selected_family_moyal_symmetric_difference_energy"


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
class MaterialSidecarCharge:
    event: str
    currency: str
    charge: float
    recursion_classification: str
    physical_work: bool = False
    carrier_stock: bool = False
    smooth_k_phys_relink: bool = False
    creates_second_coefficient_impulse: bool = False

    def __post_init__(self) -> None:
        q = _finite_nonnegative(self.charge, "material sidecar charge")
        if self.event == MATERIAL_MEMBERSHIP_EVENT:
            if self.currency != MEMBERSHIP_PROVENANCE_CURRENCY or q != 0.0:
                raise ValueError("material membership rereading is zero-charge provenance only")
            if self.recursion_classification != "zero_charge_provenance":
                raise ValueError("membership sidecar changed from zero-charge provenance")
        elif self.event == SELECTED_FAMILY_EVENT:
            if self.currency != SELECTED_FAMILY_MOYAL_CURRENCY:
                raise ValueError("selected-family switch lost its Moyal symmetric-difference currency")
            if self.recursion_classification != "zero_generation_depth_selected_service_boundary":
                raise ValueError("selected-family switch changed from zero-generation-depth service-boundary currency")
        else:
            raise ValueError("unknown material sidecar event")
        if self.physical_work or self.carrier_stock or self.smooth_k_phys_relink or self.creates_second_coefficient_impulse:
            raise ValueError(
                "material sidecar charge may not be identified with dW, inherited stock, smooth K_phys, or a second carrier impulse"
            )


@dataclass(frozen=True)
class MaterialSidecarStockDecomposition:
    carrier_id: str
    inherited_stock_mass: float
    stock_relay_label: str
    charges: tuple[MaterialSidecarCharge, ...]
    selected_family_switch_energy: float
    family_switch_moyal_certificate_bound: bool
    stock_charge_count: int = 1
    stock_cloned_per_sidecar: bool = False
    sidecar_charge_added_to_stock: bool = False
    sidecar_charge_added_to_physical_work: bool = False
    k_phys_identification_used: bool = False
    later_hahn_used: bool = False

    def __post_init__(self) -> None:
        stock = _finite_nonnegative(self.inherited_stock_mass, "inherited stock mass")
        family = _finite_nonnegative(self.selected_family_switch_energy, "selected-family switch energy")
        if stock <= 0.0 or not self.carrier_id:
            raise ValueError("positive inherited physical stock on a named carrier required")
        if self.stock_relay_label != SAME_CARRIER_INHERITED_STOCK_RELAY or self.stock_charge_count != 1:
            raise ValueError("inherited stock must remain one same-carrier physical charge")
        if tuple(sorted(c.event for c in self.charges)) != tuple(c.event for c in self.charges):
            raise ValueError("material sidecar charges must use sorted event order")
        if len({c.event for c in self.charges}) != len(self.charges):
            raise ValueError("material sidecar event was charged more than once")
        selected = tuple(c for c in self.charges if c.event == SELECTED_FAMILY_EVENT)
        if selected:
            if len(selected) != 1 or not _relative_match(selected[0].charge, family):
                raise AssertionError("stored selected-family charge changed from exact Moyal switch energy")
            if not self.family_switch_moyal_certificate_bound:
                raise AssertionError("selected-family sidecar lacks exact Moyal symmetric-difference certificate")
        elif family != 0.0 or self.family_switch_moyal_certificate_bound:
            raise AssertionError("family-switch energy/certificate exists without selected-family sidecar")
        if any(c.event == MATERIAL_MEMBERSHIP_EVENT and c.charge != 0.0 for c in self.charges):
            raise AssertionError("membership rereading created sidecar energy")
        if (
            self.stock_cloned_per_sidecar
            or self.sidecar_charge_added_to_stock
            or self.sidecar_charge_added_to_physical_work
            or self.k_phys_identification_used
            or self.later_hahn_used
        ):
            raise ValueError("material sidecar decomposition changed physical currency or causal multiplicity")

    @property
    def membership_only_or_no_sidecar(self) -> bool:
        return all(c.event == MATERIAL_MEMBERSHIP_EVENT for c in self.charges)

    @property
    def requires_selected_family_ancestry_routing(self) -> bool:
        return any(c.event == SELECTED_FAMILY_EVENT for c in self.charges)


@dataclass(frozen=True)
class SameStateSelectedFamilySwitchAntiTheorem:
    symmetric_difference_energy: float
    coherent_ledger_relink_energy: float
    positive_increment_work: float
    negative_increment_work: float
    identical_state_energy_residual: float
    generation_event_inferred: bool = False
    physical_work_inferred: bool = False

    def __post_init__(self) -> None:
        R = _finite_nonnegative(self.symmetric_difference_energy, "Moyal symmetric-difference energy")
        LR = _finite_nonnegative(self.coherent_ledger_relink_energy, "coherent-ledger relink energy")
        pp = _finite_nonnegative(self.positive_increment_work, "positive coherent increment")
        pm = _finite_nonnegative(self.negative_increment_work, "negative coherent increment")
        residual = _finite_nonnegative(self.identical_state_energy_residual, "identical-state energy residual")
        if not _relative_match(R, LR, 2.0e-13):
            raise AssertionError("material sidecar and coherent ledger disagree on selected-family Moyal charge")
        if pp > 2.0e-13 * max(R, 1.0e-300) or pm > 2.0e-13 * max(R, 1.0e-300) or residual > 2.0e-13 * max(R, 1.0e-300):
            raise AssertionError("same-state selected-family switch unexpectedly created coherent evolution/work")
        if self.generation_event_inferred or self.physical_work_inferred:
            raise ValueError("selected-family boundary charge alone cannot certify physical work or a generation event")


def same_state_selected_family_switch_anti_theorem(
    cell_energies: Sequence[float],
    old_selected: Sequence[int],
    new_selected: Sequence[int],
) -> SameStateSelectedFamilySwitchAntiTheorem:
    energies = tuple(_finite_nonnegative(x, "cell energy") for x in cell_energies)
    switch = selected_family_switch_sidecar(energies, old_selected, new_selected)
    # Use real amplitudes sqrt(E_C) at two identical observations.  The coherent
    # state is literally unchanged; only the selected family is reread.
    state = tuple(complex(math.sqrt(E), 0.0) for E in energies)
    ledger = coherent_ledger((state, state), (tuple(old_selected), tuple(new_selected)))
    initial = math.fsum(abs(z) ** 2 for z in state)
    final = math.fsum(abs(z) ** 2 for z in state)
    return SameStateSelectedFamilySwitchAntiTheorem(
        symmetric_difference_energy=float(switch["symmetric_difference_energy"]),
        coherent_ledger_relink_energy=float(ledger.relink_energy),
        positive_increment_work=float(ledger.positive_work),
        negative_increment_work=float(ledger.negative_work),
        identical_state_energy_residual=abs(final - initial),
    )


@dataclass(frozen=True)
class SelectedFamilyServiceNoEscapeBinding:
    selected_family_charge: float
    positive_selected_service: float
    negative_selected_service: float
    final_selected_energy: float
    branch: str
    branch_value: float
    threshold: float
    margin: float
    identified_with_physical_work: bool = False
    identified_with_carrier_stock: bool = False
    identified_with_smooth_k_phys: bool = False

    def __post_init__(self) -> None:
        q = _finite_nonnegative(self.selected_family_charge, "selected-family charge")
        for name, value in (
            ("positive selected service", self.positive_selected_service),
            ("negative selected service", self.negative_selected_service),
            ("final selected energy", self.final_selected_energy),
            ("branch value", self.branch_value),
            ("threshold", self.threshold),
            ("margin", self.margin),
        ):
            _finite_nonnegative(value, name)
        if self.branch not in {
            "zero_service",
            "terminal_coherent_energy",
            "backflow_or_cancellation",
            "relink_symmetric_difference",
        }:
            raise ValueError("unknown selected-cell service no-escape branch")
        if self.branch == "relink_symmetric_difference" and not _relative_match(self.branch_value, q, 2.0e-13):
            raise AssertionError("selected-family service branch changed from exact symmetric-difference charge")
        if self.identified_with_physical_work or self.identified_with_carrier_stock or self.identified_with_smooth_k_phys:
            raise ValueError("Moyal selected-family boundary energy is not dW, carrier stock, or smooth K_phys role flux")


def material_sidecar_stock_decomposition(
    certificate: SameCarrierInheritedEnergyRelayCertificate,
    *,
    selected_family_switch_certificate: Mapping[str, object] | None = None,
) -> MaterialSidecarStockDecomposition:
    events = tuple(certificate.material_sidecars)
    charges: list[MaterialSidecarCharge] = []
    family_bound = False
    family_energy = _finite_nonnegative(
        certificate.selected_family_switch_energy,
        "stored selected-family switch energy",
    )
    for event in events:
        if event == MATERIAL_MEMBERSHIP_EVENT:
            charges.append(
                MaterialSidecarCharge(
                    event=event,
                    currency=MEMBERSHIP_PROVENANCE_CURRENCY,
                    charge=0.0,
                    recursion_classification="zero_charge_provenance",
                )
            )
        elif event == SELECTED_FAMILY_EVENT:
            if selected_family_switch_certificate is None:
                raise TypeError("selected-family stock sidecar requires exact Moyal switch certificate")
            changed = bool(selected_family_switch_certificate.get("selected_family_changed", False))
            exact = _finite_nonnegative(
                float(selected_family_switch_certificate.get("symmetric_difference_energy", math.nan)),
                "exact Moyal symmetric-difference energy",
            )
            margin = float(selected_family_switch_certificate.get("jump_bound_margin", math.nan))
            if not changed or not math.isfinite(margin) or margin < -5.0e-13 * max(exact, 1.0e-300):
                raise TypeError("selected-family Moyal certificate is not an admissible changed-family charge")
            if not _relative_match(exact, family_energy):
                raise TypeError("inherited-stock certificate lost the exact selected-family Moyal charge")
            family_bound = True
            charges.append(
                MaterialSidecarCharge(
                    event=event,
                    currency=SELECTED_FAMILY_MOYAL_CURRENCY,
                    charge=exact,
                    recursion_classification="zero_generation_depth_selected_service_boundary",
                )
            )
        else:
            raise TypeError("same-carrier stock certificate contains an unrecognized material sidecar")
    return MaterialSidecarStockDecomposition(
        carrier_id=certificate.carrier_id,
        inherited_stock_mass=certificate.initial_energy,
        stock_relay_label=SAME_CARRIER_INHERITED_STOCK_RELAY,
        charges=tuple(sorted(charges, key=lambda c: c.event)),
        selected_family_switch_energy=family_energy,
        family_switch_moyal_certificate_bound=family_bound,
    )


def selected_family_service_no_escape_binding(
    decomposition: MaterialSidecarStockDecomposition,
    *,
    positive_selected_service: float,
    negative_selected_service: float,
    final_selected_energy: float,
) -> SelectedFamilyServiceNoEscapeBinding:
    selected = tuple(c for c in decomposition.charges if c.event == SELECTED_FAMILY_EVENT)
    if len(selected) != 1:
        raise TypeError("selected-family service binding requires exactly one Moyal sidecar charge")
    route = service_no_escape(
        float(positive_selected_service),
        float(negative_selected_service),
        float(final_selected_energy),
        selected[0].charge,
    )
    return SelectedFamilyServiceNoEscapeBinding(
        selected_family_charge=selected[0].charge,
        positive_selected_service=float(positive_selected_service),
        negative_selected_service=float(negative_selected_service),
        final_selected_energy=float(final_selected_energy),
        branch=str(route["branch"]),
        branch_value=float(route.get("branch_value", 0.0)),
        threshold=float(route["threshold"]),
        margin=float(route["margin"]),
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "membership": "intrinsic material membership rereading preserves the existing positive service total and contributes exactly zero sidecar energy; its label is retained as zero-charge provenance",
        "selected_family": "a selected-family change must bind the stock certificate's stored switch energy to the exact Moyal symmetric-difference energy before downstream use; this R_switch is separate ancestry/service currency",
        "service_binding": "the exact selected-cell law P_plus <= E_final + P_minus + R_switch reads the same R_switch in its relink_symmetric_difference candidate.  A same coherent state with only S_old->S_new changed has the identical R_switch but zero increments, so this is selected-family boundary energy, not smooth K_phys or generation work",
        "currency_separation": "E0 remains one inherited carrier-stock charge; membership provenance and selected-family Moyal charge are never added to dW, stock E0, or K/S work and never clone stock once per sidecar",
        "recurrence_scope": "membership is zero-charge provenance.  The selected-family Moyal charge itself also has zero generation depth: the exact same-state anti-theorem produces R_switch>0 with P_plus=P_minus=0 and no state change.  Genuine material/source recursion therefore requires separate physical service/source evidence, not the selected-family boundary charge alone",
        "later_hahn_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class MaterialSidecarStockStress:
    samples: int
    membership_cases: int
    selected_family_cases: int
    mixed_sidecar_cases: int
    zero_charge_membership_violations: int
    maximum_moyal_charge_native_residual: float
    maximum_stock_mass_native_residual: float
    selected_family_relink_branch_cases: int
    same_state_positive_switch_cases: int
    maximum_same_state_increment_native_residual: float
    smooth_k_phys_identifications: int


def stress(samples: int = 50_000, seed: int = 2026081305) -> MaterialSidecarStockStress:
    # This stress is intentionally representation-native: it attacks the exact
    # Moyal sidecar charge and stock/currency decomposition.  Actual Navier--Stokes
    # carrier dynamics are refereed separately by the certified inherited-stock PDE probe.
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = random.Random(int(seed))
    from src.material_label_carrier_quotient import carrier_registration_with_material_sidecars
    from src.same_carrier_checkpoint_segmentation_quotient import partition_same_carrier_path
    from src.same_carrier_inherited_energy_relay import same_carrier_inherited_energy_relay

    membership = family = mixed = membership_bad = relink_branch = same_state_positive = kbad = 0
    worst_charge = worst_stock = worst_same_state = 0.0
    for j in range(count):
        mode = j % 3
        membership_change = mode in (0, 2)
        family_change = mode in (1, 2)
        switch_cert = None
        switch_energy = 0.0
        if family_change:
            ncell = rng.randint(4, 12)
            energies = [10.0 ** rng.uniform(-5.0, 1.0) for _ in range(ncell)]
            old = {i for i in range(ncell) if rng.random() < 0.45}
            new = {i for i in range(ncell) if rng.random() < 0.45}
            if old == new:
                new = set(new)
                k = rng.randrange(ncell)
                if k in new:
                    new.remove(k)
                else:
                    new.add(k)
            switch_cert = selected_family_switch_sidecar(energies, sorted(old), sorted(new))
            switch_energy = float(switch_cert["symmetric_difference_energy"])
            anti = same_state_selected_family_switch_anti_theorem(energies, sorted(old), sorted(new))
            same_state_positive += int(anti.symmetric_difference_energy > 0.0)
            same_scale = max(anti.symmetric_difference_energy, 1.0e-300)
            worst_same_state = max(
                worst_same_state,
                anti.positive_increment_work / same_scale,
                anti.negative_increment_work / same_scale,
                anti.identical_state_energy_residual / same_scale,
            )

        amp = 1.0
        ih = complex(rng.uniform(-0.12, 0.12), rng.uniform(-0.12, 0.12))
        ii = complex(rng.uniform(-0.08, 0.08), rng.uniform(-0.08, 0.08))
        z_event = complex(amp, 0.0)
        material = carrier_registration_with_material_sidecars(
            z_event,
            z_event - ih - ii,
            ih,
            ii,
            intrinsic_material_membership_change=membership_change,
            selected_family_change=family_change,
            selected_family_switch_energy=switch_energy,
            same_smooth_role=True,
            same_analysis_probe=True,
        )
        if not material["carrier_continuation_certified"]:
            raise AssertionError("sidecar stress unexpectedly crossed a carrier coefficient face")

        horizon = 10.0 ** rng.uniform(-5.0, -2.0)
        segments = partition_same_carrier_path(
            carrier_id=f"material-sidecar-{j}",
            terminal_amplitude=1.0,
            elapsed_times=(0.0, 0.25 * horizon, 0.5 * horizon, 0.75 * horizon, horizon),
            strain_action=(0.0, 0.004, 0.008, 0.012, 0.016),
            residual_impulse_abs=(0.0, 0.03, 0.05, 0.04, 0.06),
            hh_impulse_abs=(0.0, 0.08, 0.12, 0.10, 0.16),
            checkpoint_indices=(2,),
        )
        e1 = 10.0 ** rng.uniform(-3.0, 1.0)
        e0 = rng.uniform(0.25, 1.5) * e1
        wr = rng.uniform(0.0, 0.15) * e1
        stock = same_carrier_inherited_energy_relay(
            segments,
            initial_time=2.0,
            terminal_time=2.0 + horizon,
            initial_energy=e0,
            terminal_energy=e1,
            residual_positive_work=wr,
            strain_action=0.016,
            material_registration=material,
            initial_endpoint_is_non_event_carrier_slice=True,
        )
        out = material_sidecar_stock_decomposition(
            stock,
            selected_family_switch_certificate=switch_cert,
        )
        membership += int(membership_change)
        family += int(family_change)
        mixed += int(membership_change and family_change)
        membership_bad += sum(int(c.event == MATERIAL_MEMBERSHIP_EVENT and c.charge != 0.0) for c in out.charges)
        worst_stock = max(worst_stock, abs(out.inherited_stock_mass - e0) / max(abs(e0), 1.0e-300))
        if family_change:
            exact = float(switch_cert["symmetric_difference_energy"])
            worst_charge = max(worst_charge, abs(out.selected_family_switch_energy - exact) / max(exact, 1.0e-300))
            # Force an actual relink_symmetric_difference service branch in half
            # the family cases, while the other half exercises competing terminal/backflow branches.
            if j % 2 == 0:
                pplus = max(1.0e-300, 2.0 * exact)
                pminus = 0.1 * pplus
                efinal = 0.1 * pplus
                if pplus > pminus + efinal + exact:
                    pplus = pminus + efinal + exact
            else:
                pplus = max(1.0e-300, exact)
                pminus = 1.2 * pplus
                efinal = 0.2 * pplus
            binding = selected_family_service_no_escape_binding(
                out,
                positive_selected_service=pplus,
                negative_selected_service=pminus,
                final_selected_energy=efinal,
            )
            relink_branch += int(binding.branch == "relink_symmetric_difference")
            kbad += int(binding.identified_with_smooth_k_phys)
    if membership_bad or kbad or worst_charge > 5.0e-12 or worst_stock > 5.0e-12 or worst_same_state > 5.0e-12:
        raise AssertionError("material sidecar stock stress changed physical currency/provenance")
    if family and (relink_branch == 0 or same_state_positive == 0):
        raise AssertionError("stress did not exercise both Moyal service branch and positive same-state switch anti-theorem")
    return MaterialSidecarStockStress(
        samples=count,
        membership_cases=membership,
        selected_family_cases=family,
        mixed_sidecar_cases=mixed,
        zero_charge_membership_violations=membership_bad,
        maximum_moyal_charge_native_residual=worst_charge,
        maximum_stock_mass_native_residual=worst_stock,
        selected_family_relink_branch_cases=relink_branch,
        same_state_positive_switch_cases=same_state_positive,
        maximum_same_state_increment_native_residual=worst_same_state,
        smooth_k_phys_identifications=kbad,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--seed", type=int, default=2026081305)
    ap.add_argument("--outdir", type=Path, default=Path("results-material-sidecar-stock-owner-decomposition"))
    args = ap.parse_args()
    out = stress(args.samples, args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)
    payload = {"certificate": theorem_certificate(), "stress": asdict(out)}
    (args.outdir / "material_sidecar_stock_owner_decomposition.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = f"""# Material-sidecar / inherited-stock charge decomposition

Status: **{STATUS}**.

Inherited carrier stock remains one charge.  Intrinsic membership rereading is zero-charge provenance.  A selected-family change is accepted downstream only after its stored switch energy is bound to the exact Moyal symmetric-difference certificate; that charge remains separate ancestry/service currency.  It is not `dW`, carrier stock, a second coefficient impulse, or smooth `K_phys` relink.

Stress: `{out.samples}` exact sidecar/stock registrations
- membership / selected-family / mixed cases: `{out.membership_cases}` / `{out.selected_family_cases}` / `{out.mixed_sidecar_cases}`
- zero-charge membership violations: `{out.zero_charge_membership_violations}`
- maximum Moyal charge native residual: `{out.maximum_moyal_charge_native_residual:.3e}`
- maximum inherited-stock mass native residual: `{out.maximum_stock_mass_native_residual:.3e}`
- selected-family `relink_symmetric_difference` service branches: `{out.selected_family_relink_branch_cases}`
- positive-charge same-state selected-family switches: `{out.same_state_positive_switch_cases}`
- maximum same-state coherent increment residual: `{out.maximum_same_state_increment_native_residual:.3e}`
- smooth `K_phys` identifications: `{out.smooth_k_phys_identifications}`

The selected-family **boundary charge itself** is zero generation depth.  Genuine material/source recurrence still requires independent physical service/source evidence.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(summary, encoding="utf-8")


if __name__ == "__main__":
    main()
