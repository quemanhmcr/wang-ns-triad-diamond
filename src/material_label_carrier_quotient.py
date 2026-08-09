from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.coherent_transfer_cells import selection_jump, symmetric_difference_energy
from src.common_slice_coefficient_registration import (
    GENERATED_FRACTION,
    INHERIT_FRACTION,
    RESIDUAL_FRACTION,
    exact_adjoint_residual,
    registration_first_stop,
)
from src.heat_edge_material_ownership import partition_positive_edge_measure


MATERIAL_MEMBERSHIP_EVENT = "intrinsic_material_membership_update"
SELECTED_FAMILY_EVENT = "selected_family_moyal_switch"
ROLE_DELEGATE_EVENT = "role_or_probe_change_requires_interface_registration"
INTERFACE_STOP = "classified_role_interface_impulse"
HH_STOP = "hh_regeneration_impulse"


def reclassify_positive_service_sidecar(
    edge_weights: Sequence[float],
    old_here_before: Sequence[bool],
    old_neighbor_before: Sequence[bool],
    old_here_after: Sequence[bool],
    old_neighbor_after: Sequence[bool],
) -> dict[str, object]:
    """Re-read OO/ON/NN without creating or destroying the positive service law.

    The edge weights are held fixed.  Only the two intrinsic endpoint membership
    indicators are updated.  Hence before/after OO+ON+NN totals are exactly the
    same positive Moyal/heat service.  Category movement is bookkeeping of actual
    material ownership, not a new carrier source or coefficient impulse.
    """
    w = np.asarray(edge_weights, float)
    a0 = np.asarray(old_here_before, bool)
    b0 = np.asarray(old_neighbor_before, bool)
    a1 = np.asarray(old_here_after, bool)
    b1 = np.asarray(old_neighbor_after, bool)
    if w.ndim != 1 or any(x.shape != w.shape for x in (a0, b0, a1, b1)):
        raise ValueError("matching one-dimensional service/ownership arrays required")
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("finite nonnegative service weights required")
    before = partition_positive_edge_measure(w, a0, b0)
    after = partition_positive_edge_measure(w, a1, b1)
    total_residual = float(after["total"] - before["total"])
    category_delta = {
        key: float(after[key] - before[key])
        for key in ("old_old", "old_new_interface", "new_new")
    }
    changed = np.array(
        [
            (bool(x0), bool(y0)) != (bool(x1), bool(y1))
            for x0, y0, x1, y1 in zip(a0, b0, a1, b1)
        ],
        dtype=bool,
    )
    return {
        "before": before,
        "after": after,
        "total_service_residual": total_residual,
        "category_delta": category_delta,
        "category_delta_sum": float(sum(category_delta.values())),
        "changed_endpoint_pair_service": float(w[changed].sum()),
        "service_created_by_relabel": 0.0,
    }


def selected_family_switch_sidecar(
    cell_energies: Sequence[float],
    old_selected: Sequence[int],
    new_selected: Sequence[int],
) -> dict[str, float | bool]:
    """Keep the exact Moyal symmetric-difference charge when the selected family changes.

    This theorem does *not* make selected-family switching free.  It only says
    that, if the smooth PDE role Q and analysis probe psi are unchanged, this
    observation-set jump is sidecar ancestry/service bookkeeping rather than a
    second term in the carrier coefficient Duhamel identity.
    """
    jump = float(selection_jump(cell_energies, old_selected, new_selected))
    charge = float(symmetric_difference_energy(cell_energies, old_selected, new_selected))
    margin = charge - abs(jump)
    tol = 3e-13 * max(1.0, charge, abs(jump))
    if margin < -tol:
        raise AssertionError("selected-family energy jump exceeded symmetric-difference Moyal charge")
    return {
        "selection_energy_jump": jump,
        "symmetric_difference_energy": charge,
        "jump_bound_margin": margin,
        "selected_family_changed": set(map(int, old_selected)) != set(map(int, new_selected)),
    }


def carrier_registration_with_material_sidecars(
    z_event: complex,
    z_slice: complex,
    i_hh: complex,
    i_interface: complex,
    *,
    intrinsic_material_membership_change: bool = False,
    selected_family_change: bool = False,
    selected_family_switch_energy: float = 0.0,
    same_smooth_role: bool = True,
    same_analysis_probe: bool = True,
) -> dict[str, object]:
    """Carrier-level quotient of pure material sidecars from true role changes.

    When Q and psi are the same objects on both sides, the exact registered
    coefficient identity is still

        z_event = z_slice + I_HH + I_interface.

    It contains no old-pool indicator and no selected-family characteristic
    function.  Therefore material membership updates and selected-family Moyal
    switches add *no second coefficient impulse*.  They remain explicit sidecar
    physical events/currencies while carrier continuation is decided only by the
    existing interface and HH amplitude faces.

    If Q or psi really changes, this quotient is inapplicable.  The result is
    delegated to event-role / nonaffine-interface registration rather than being
    declared transparent or assigned an artificial jump here.
    """
    switch_energy = float(selected_family_switch_energy)
    if switch_energy < 0 or not math.isfinite(switch_energy):
        raise ValueError("finite nonnegative selected-family switch energy required")
    if not selected_family_change and switch_energy > 1e-15:
        raise ValueError("nonzero switch energy requires an explicit selected-family change")

    sidecars: list[str] = []
    if intrinsic_material_membership_change:
        sidecars.append(MATERIAL_MEMBERSHIP_EVENT)
    if selected_family_change:
        sidecars.append(SELECTED_FAMILY_EVENT)

    if not same_smooth_role or not same_analysis_probe:
        return {
            "branch": "delegate_role_or_probe_change",
            "quotient_applicable": False,
            "carrier_continuation_certified": False,
            "carrier_stop_causes": (ROLE_DELEGATE_EVENT,),
            "sidecar_events": tuple(sidecars),
            "joint_physical_events": tuple(sidecars + [ROLE_DELEGATE_EVENT]),
            "selected_family_switch_energy": switch_energy,
            "sidecar_requires_accounting": bool(sidecars),
            "primary_selected": False,
        }

    ze = complex(z_event)
    zs = complex(z_slice)
    ih = complex(i_hh)
    ir = complex(i_interface)
    amp = abs(ze)
    if amp <= 0 or not math.isfinite(amp):
        raise ValueError("nonzero finite event coefficient required")
    res = abs(exact_adjoint_residual(ze, zs, ih, ir))
    tol = 4e-12 * max(1.0, amp, abs(zs), abs(ih), abs(ir))
    if res > tol:
        raise ValueError("same-role carrier Duhamel decomposition is not exact")

    stops: list[str] = []
    if abs(ir) >= RESIDUAL_FRACTION * amp - tol:
        stops.append(INTERFACE_STOP)
    if abs(ih) >= GENERATED_FRACTION * amp - tol:
        stops.append(HH_STOP)

    joint = tuple(sidecars + stops)
    if stops:
        return {
            "branch": "carrier_first_stop_with_material_sidecars" if sidecars else "carrier_first_stop",
            "quotient_applicable": True,
            "carrier_continuation_certified": False,
            "carrier_stop_causes": tuple(stops),
            "sidecar_events": tuple(sidecars),
            "joint_physical_events": joint,
            "selected_family_switch_energy": switch_energy,
            "duhamel_residual": res,
            "sidecar_requires_accounting": bool(sidecars),
            "primary_selected": False,
        }

    inherited = abs(zs)
    clean = INHERIT_FRACTION * amp
    if inherited < clean - tol:
        raise AssertionError("same-role no-hit carrier lost the one-quarter coefficient")
    return {
        "branch": "carrier_continues_with_material_sidecars" if sidecars else "carrier_continues",
        "quotient_applicable": True,
        "carrier_continuation_certified": True,
        "carrier_stop_causes": (),
        "sidecar_events": tuple(sidecars),
        "joint_physical_events": joint,
        "selected_family_switch_energy": switch_energy,
        "event_amplitude": amp,
        "slice_amplitude": inherited,
        "clean_slice_amplitude_lower": clean,
        "duhamel_residual": res,
        "sidecar_requires_accounting": bool(sidecars),
        "same_carrier_reusable_after_sidecar": True,
        "primary_selected": False,
    }


def legacy_relink_refinement_certificate() -> dict[str, object]:
    """Exhibit the precise subtype refined relative to the older conservative stop.

    The old common-slice theorem remains valid: it stopped on any boolean
    ``material_relink``.  The new theorem refines only cases where Q and psi are
    unchanged and all dynamics are already in I_interface.  It does not rewrite
    the old theorem or declare genuine role changes transparent.
    """
    ze = 1.0 + 0j
    ih = 0.1 + 0j
    ir = 0.1j
    zs = ze - ih - ir
    legacy = registration_first_stop(ze, zs, ih, ir, material_relink=True)
    refined = carrier_registration_with_material_sidecars(
        ze,
        zs,
        ih,
        ir,
        intrinsic_material_membership_change=True,
        same_smooth_role=True,
        same_analysis_probe=True,
    )
    if bool(legacy["continuing"]):
        raise AssertionError("legacy common-slice material-relink stop unexpectedly continued")
    if not bool(refined["carrier_continuation_certified"]):
        raise AssertionError("pure-label refined quotient did not continue the same carrier")
    return {
        "legacy_branch": str(legacy["branch"]),
        "refined_branch": str(refined["branch"]),
        "old_theorem_status": "retained_as_conservative_superset",
        "refined_subtype": "same_Q_same_psi_material_sidecar_only",
    }


def theorem_certificate() -> dict[str, object]:
    legacy = legacy_relink_refinement_certificate()
    return {
        "status": "EXACT_MATERIAL_LABEL_SIDECAR_QUOTIENT__PURE_LABEL_AND_SELECTED_FAMILY_SWITCH_ADD_NO_SECOND_CARRIER_IMPULSE__ROLE_CHANGE_REMAINS_INTERFACE",
        "coefficient_identity": "for the same smooth role Q and same registered analysis probe psi, z_event=z_slice+I_HH+I_interface contains no old-pool membership indicator or selected-family characteristic function",
        "no_double_count": "material-address/membership motion may be physically caused by nonaffine dynamics, but that dynamics is already represented by I_interface; crossing a material bookkeeping boundary does not add a second independent coefficient impulse",
        "service_reclassification": "holding the positive renewed service weights fixed, rereading OO/ON/NN with updated intrinsic endpoint memberships preserves the total service exactly",
        "selected_family": "a selected coherent-family change keeps its exact symmetric-difference Moyal energy R_switch and |Delta E_selected|<=R_switch; this charge remains ancestry/service currency even when the same smooth carrier is reused",
        "carrier_quotient": "with same Q and psi, pure material sidecars do not join the carrier-stop set; carrier continuation is still decided only by the existing interface |I|>=|z|/4 and HH |I_HH|>=|z|/2 faces",
        "joint_events": "sidecar material events and carrier stops may occur simultaneously and are returned together without lexicographic primary selection, while sidecar events are separately marked so they do not kill the reusable carrier",
        "nonquotient": "if Q or psi changes, no transparency is claimed: delegate to event-role/nonaffine-interface/relink registration and preserve any Moyal symmetric-difference charge",
        "master": "the quotient removes only an unnecessary carrier reconstruction; it does not erase material ancestry, R_switch, entropy/cycle accounting, or promote relink to a free/reset resource",
        "legacy_refinement": legacy,
        "scope": "this justifies the three-monitor shortest critical-shell corridor and eventwise rereading of material ownership after service; source/pressure routing and final continuum master assembly remain open",
    }


@dataclass(frozen=True)
class MaterialLabelCarrierQuotientStress:
    samples: int
    minimum_continuing_quarter_margin: float
    worst_same_role_duhamel_residual: float
    worst_service_reclassification_residual: float
    minimum_switch_jump_margin: float
    label_transparency_failures: int
    switch_charge_losses: int
    role_delegate_failures: int
    maximum_carrier_stop_count: int
    maximum_joint_physical_event_count: int
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260809) -> MaterialLabelCarrierQuotientStress:
    rng = np.random.default_rng(seed)
    mq = float("inf")
    wr = ws = 0.0
    mj = float("inf")
    label_fail = switch_fail = role_fail = 0
    max_stops = max_joint = 0
    branches: dict[str, int] = {}

    for _ in range(samples):
        amp = float(math.exp(rng.uniform(-6.0, 6.0)))
        ph = float(rng.uniform(-math.pi, math.pi))
        ze = amp * complex(math.cos(ph), math.sin(ph))
        mode = int(rng.integers(0, 5))
        membership = bool(rng.integers(0, 2))
        family_change = bool(rng.integers(0, 2))

        # Exact selected-family Moyal charge sidecar.
        ne = int(rng.integers(3, 30))
        e = rng.lognormal(mean=-2.0, sigma=1.5, size=ne)
        ids = np.arange(ne)
        rng.shuffle(ids)
        cut0 = int(rng.integers(0, ne + 1))
        cut1 = int(rng.integers(0, ne + 1))
        old = ids[:cut0].tolist()
        new = ids[-cut1:].tolist() if cut1 else []
        sw = selected_family_switch_sidecar(e, old, new)
        mj = min(mj, float(sw["jump_bound_margin"]))
        actual_change = bool(sw["selected_family_changed"])
        if family_change and not actual_change:
            # Force an actual family change for the branch flag without inventing energy.
            if ne > 0:
                new = list(set(new).symmetric_difference({int(ids[0])}))
                sw = selected_family_switch_sidecar(e, old, new)
                actual_change = bool(sw["selected_family_changed"])
        family_change = actual_change if family_change else False
        switch_energy = float(sw["symmetric_difference_energy"]) if family_change else 0.0

        # Service law is fixed while endpoint old/new ownership is reread.
        nw = int(rng.integers(1, 80))
        weights = rng.lognormal(mean=-2.0, sigma=1.7, size=nw)
        a0 = rng.random(nw) < 0.5
        b0 = rng.random(nw) < 0.5
        a1 = rng.random(nw) < 0.5
        b1 = rng.random(nw) < 0.5
        rec = reclassify_positive_service_sidecar(weights, a0, b0, a1, b1)
        rr = max(abs(float(rec["total_service_residual"])), abs(float(rec["category_delta_sum"])))
        ws = max(ws, rr)
        if rr > 3e-12 * max(1.0, float(weights.sum())):
            raise AssertionError("material relabel created or destroyed positive service")

        if mode == 4:
            # A true role/probe change is never quotiented.
            out = carrier_registration_with_material_sidecars(
                ze,
                0j,
                0j,
                0j,
                intrinsic_material_membership_change=membership,
                selected_family_change=family_change,
                selected_family_switch_energy=switch_energy,
                same_smooth_role=bool(rng.integers(0, 2)),
                same_analysis_probe=False,
            )
            if bool(out["quotient_applicable"]) or bool(out["carrier_continuation_certified"]):
                role_fail += 1
                raise AssertionError("genuine role/probe change was incorrectly quotiented")
        else:
            rfrac = float(rng.uniform(0.0, 0.20))
            hfrac = float(rng.uniform(0.0, 0.40))
            if mode == 1:
                rfrac = float(rng.uniform(0.26, 0.45))
            elif mode == 2:
                hfrac = float(rng.uniform(0.51, 0.80))
            elif mode == 3:
                rfrac = RESIDUAL_FRACTION
                hfrac = GENERATED_FRACTION
            pr = float(rng.uniform(-math.pi, math.pi))
            phh = float(rng.uniform(-math.pi, math.pi))
            ir = rfrac * amp * complex(math.cos(pr), math.sin(pr))
            ih = hfrac * amp * complex(math.cos(phh), math.sin(phh))
            zs = ze - ir - ih
            out = carrier_registration_with_material_sidecars(
                ze,
                zs,
                ih,
                ir,
                intrinsic_material_membership_change=membership,
                selected_family_change=family_change,
                selected_family_switch_energy=switch_energy,
                same_smooth_role=True,
                same_analysis_probe=True,
            )
            wr = max(wr, float(out["duhamel_residual"]))
            max_stops = max(max_stops, len(tuple(out["carrier_stop_causes"])))
            max_joint = max(max_joint, len(tuple(out["joint_physical_events"])))

            # Sidecars may change the master ledger but not the same-role carrier decision.
            base = carrier_registration_with_material_sidecars(ze, zs, ih, ir)
            if (
                bool(out["carrier_continuation_certified"]) != bool(base["carrier_continuation_certified"])
                or tuple(out["carrier_stop_causes"]) != tuple(base["carrier_stop_causes"])
            ):
                label_fail += 1
                raise AssertionError("pure material sidecar changed the same-role carrier first-stop decision")
            if family_change and abs(float(out["selected_family_switch_energy"]) - switch_energy) > 2e-14 * max(1.0, switch_energy):
                switch_fail += 1
                raise AssertionError("carrier quotient lost selected-family Moyal switch charge")
            if bool(out["carrier_continuation_certified"]):
                margin = float(out["slice_amplitude"]) - INHERIT_FRACTION * amp
                mq = min(mq, margin)
                if margin < -3e-12 * max(1.0, amp):
                    raise AssertionError("transparent material sidecar lost quarter coefficient")

        b = str(out["branch"])
        branches[b] = branches.get(b, 0) + 1

    if not math.isfinite(mq):
        mq = 0.0
    return MaterialLabelCarrierQuotientStress(
        samples,
        mq,
        wr,
        ws,
        mj,
        label_fail,
        switch_fail,
        role_fail,
        max_stops,
        max_joint,
        branches,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-material-label-carrier-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "material_label_carrier_quotient.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Material-label sidecar quotient for smooth carrier continuation

Status: **{cert['status']}**.

For a fixed smooth PDE role `Q` and the same registered analysis probe `psi`, the exact coefficient identity is

`z_event = z_slice + I_HH + I_interface`.

There is no old-pool indicator and no selected coherent-family characteristic function in this identity.  A material address/membership crossing may be caused by genuine nonaffine motion, but the effect of that motion on the coefficient is already in `I_interface`; the bookkeeping crossing does not generate a second independent Duhamel impulse.

Therefore, with the same `Q` and `psi`, pure material sidecars do **not** enter the carrier-stop set.  The carrier still stops only when

`|I_interface| >= |z_event|/4`

or

`|I_HH| >= |z_event|/2`.

If neither face is hit, the exact triangle gives the same

`|z_slice| >= |z_event|/4`

whether or not intrinsic old/new membership changed and whether or not a selected coherent family was changed for the service/ancestry ledger.

Nothing is made free.  If the selected coherent family changes, its exact positive Moyal symmetric-difference charge is retained:

`|E(S_new)-E(S_old)| <= R_switch = E(S_new symmetric_difference S_old)`.

`R_switch` remains ancestry/service currency and remains available to the existing service/ancestry routing.  The theorem only says that master may reuse the same smooth carrier rather than reconstructing it merely because the observation/material sidecar changed.

Likewise, rereading OO/ON/NN after an intrinsic material-membership update holds the positive service weights fixed.  Before and after classification have exactly the same total service; only the OO/ON/NN allocation changes.  Material ownership is therefore eventwise exact without becoming a carrier source term.

If `Q` or `psi` truly changes, this quotient refuses to act.  Such a case is delegated to event-role/nonaffine-interface/relink registration, with every Moyal switch charge preserved.  Thus the theorem distinguishes **material-label/selected-family bookkeeping** from **physical role change** instead of calling both `material_relink` at the carrier level.

The older common-slice theorem remains a correct conservative superset: it stopped on any boolean material relink.  The new theorem refines only the subtype `same Q + same psi + material sidecar`; it does not invalidate the old route.

Stress: `{out.samples}` exact coefficient/service/switch states
- minimum transparent quarter-coefficient margin: `{out.minimum_continuing_quarter_margin:.3e}`
- worst same-role Duhamel residual: `{out.worst_same_role_duhamel_residual:.3e}`
- worst service-total/reclassification residual: `{out.worst_service_reclassification_residual:.3e}`
- minimum selected-family switch-jump margin: `{out.minimum_switch_jump_margin:.3e}`
- label-transparency failures: `{out.label_transparency_failures}`
- switch-charge losses: `{out.switch_charge_losses}`
- role-delegate failures: `{out.role_delegate_failures}`
- maximum carrier-stop count: `{out.maximum_carrier_stop_count}`
- maximum simultaneous physical-event count including sidecars: `{out.maximum_joint_physical_event_count}`
- outcomes: `{out.branch_counts}`

This removes a non-PDE obstruction from the shortest carrier architecture while preserving material ancestry, symmetric-difference energy, and true role-interface physics.  Source/pressure routing and the final continuum master assembly remain open.  No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
