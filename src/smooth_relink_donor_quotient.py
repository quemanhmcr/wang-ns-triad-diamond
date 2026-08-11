from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.resolved_interface_donor_quotient import skew_donor_closure, skew_flux_balance
from src.smooth_quadratic_carrier_interface import (
    RELINK_OWNER,
    GaugeQuotientedInterfaceWork,
    positive_smooth_interface_split,
)


STATUS = (
    "EXACT_SMOOTH_PHYSICAL_RELINK_DONOR_QUOTIENT__"
    "GAUGE_QUOTIENTED_KPHYS_PAIR_FLUX__"
    "FINITE_SAME_EVENT_NEGATIVE_DONOR_CLOSURE__"
    "SMOOTH_RELINK_ZERO_RECURSION_DEPTH"
)

SMOOTH_RELINK_SAME_EVENT_RELAY = "smooth_physical_relink_same_event_donor_relay"


@dataclass(frozen=True)
class SmoothRelinkDonorCertificate:
    """Master-facing certificate that smooth K_phys relink is same-event flux.

    The certificate is bound to the exact pair matrix stored in the already
    gauge-quotiented smooth-interface work object.  It does not identify the
    smooth measure with the hard event-role measure; it uses only the common
    finite antisymmetric-flux lemma after each measure has established its own
    pair law.
    """

    relink_owner: str
    recipient_roles: tuple[int, ...]
    terminal_negative_net_donor_roles: tuple[int, ...]
    maximum_shortest_donor_path_length: int
    role_count: int
    positive_relink_work: float
    recipient_positive_incoming_flux: float
    pair_antisymmetry_residual: float
    row_binding_residual: float
    total_relink_work_residual: float
    same_physical_event: bool = True
    same_physical_time: bool = True
    new_causal_charge_created: bool = False
    recursive_generation_created: bool = False
    scale_progress_created: bool = False

    def __post_init__(self) -> None:
        if self.relink_owner != RELINK_OWNER:
            raise ValueError("smooth relink donor certificate must bind the canonical relink owner")
        if self.role_count < 2:
            raise ValueError("at least two smooth roles required for relink donor closure")
        if not self.recipient_roles or not self.terminal_negative_net_donor_roles:
            raise ValueError("positive recipients and negative-net donors required")
        if any(x < 0 or x >= self.role_count for x in (*self.recipient_roles, *self.terminal_negative_net_donor_roles)):
            raise ValueError("smooth relink role index outside certified partition")
        if not (0 <= self.maximum_shortest_donor_path_length <= self.role_count - 1):
            raise ValueError("finite donor path exceeds smooth role count")
        vals = (
            self.positive_relink_work,
            self.recipient_positive_incoming_flux,
            self.pair_antisymmetry_residual,
            self.row_binding_residual,
            self.total_relink_work_residual,
        )
        if not all(math.isfinite(x) and x >= 0 for x in vals):
            raise ValueError("finite nonnegative smooth relink donor diagnostics required")
        if self.positive_relink_work <= 0:
            raise ValueError("positive smooth relink work required")
        if not self.same_physical_event or not self.same_physical_time:
            raise ValueError("smooth relink donor trace must remain at the same physical event/time")
        if self.new_causal_charge_created or self.recursive_generation_created or self.scale_progress_created:
            raise ValueError("smooth conservative donor trace cannot create charge, recursion, or scale progress")


def _bound_smooth_relink_pair_matrix(
    interface_work: GaugeQuotientedInterfaceWork,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    """Validate that the pair matrix is the exact K_phys law of this certificate."""
    if not isinstance(interface_work, GaugeQuotientedInterfaceWork):
        raise TypeError("gauge-quotiented smooth-interface work certificate required")
    relink = np.asarray(interface_work.signed_physical_relink_atoms, dtype=float)
    T = np.asarray(interface_work.signed_physical_relink_pair_matrix, dtype=float)
    if relink.ndim != 1 or len(relink) < 2:
        raise ValueError("at least two smooth relink role atoms required")
    if T.shape != (len(relink), len(relink)):
        raise ValueError("smooth relink donor quotient requires the bound square K_phys pair matrix")
    if np.any(~np.isfinite(T)) or np.any(~np.isfinite(relink)):
        raise ValueError("finite smooth relink atoms and pair matrix required")
    scale = max(1.0, float(np.max(np.abs(T))), float(np.max(np.abs(relink))))
    antisym = float(np.max(np.abs(T + T.T)))
    row = float(np.max(np.abs(T.sum(axis=1) - relink)))
    total = abs(float(relink.sum()))
    if max(antisym, row, total) > 8e-11 * scale:
        raise ValueError("smooth relink pair matrix is not the exact conservative K_phys row law")
    return T, relink, scale, antisym, row


def smooth_relink_donor_quotient(
    interface_work: GaugeQuotientedInterfaceWork,
) -> dict[str, object]:
    """Trace positive smooth relink to finite negative-net donors at one event.

    Let eta_a=A_a^2 with sum eta_a=I and let K_phys*=-K_phys.  The certified
    smooth interface already provides

        T_ab=-2 Re <eta_a u, K_phys eta_b u>,
        T_ab=-T_ba,
        R_a=sum_b T_ab.

    For F[b->a]=[T_ab]_+, R_a=inflow_a-outflow_a.  Starting from every role with
    R_a>0 and closing backward under positive inflow must encounter a role with
    R_a<0.  Otherwise the finite closure has positive total net work and no
    positive external inflow, contradicting the exact subset divergence law.

    The trace is therefore provenance inside one physical interaction, not a new
    Navier--Stokes generation.
    """
    split = positive_smooth_interface_split(interface_work)
    if RELINK_OWNER not in tuple(split["joint_physical_owners"]):
        raise ValueError("smooth relink donor quotient requires relink to be a realized positive interface owner")

    T, relink, scale, antisym, row = _bound_smooth_relink_pair_matrix(interface_work)
    tol = 32.0 * math.ulp(scale)
    recipients = tuple(int(i) for i, x in enumerate(relink) if float(x) > tol)
    if not recipients:
        raise AssertionError("positive smooth relink owner has no positive-net recipient role")

    closure = skew_donor_closure(T, recipients)
    flux = skew_flux_balance(T)
    positive_relink = float(np.maximum(relink, 0.0).sum())
    incoming = float(closure["recipient_positive_incoming_flux"])
    total_residual = abs(float(relink.sum()))
    if abs(float(closure["recipient_net_gain"]) - positive_relink) > 8e-11 * scale:
        raise AssertionError("all positive smooth relink rows were not included in donor closure")
    if incoming + 8e-11 * scale < positive_relink:
        raise AssertionError("smooth relink recipient gain exceeds actual same-event donor inflow")
    if abs(float(flux["total_net_skew_work"])) > 8e-11 * scale:
        raise AssertionError("smooth relink pair law created net energy")

    donors = tuple(int(x) for x in closure["terminal_negative_net_donor_roles"])
    cert = SmoothRelinkDonorCertificate(
        relink_owner=RELINK_OWNER,
        recipient_roles=recipients,
        terminal_negative_net_donor_roles=donors,
        maximum_shortest_donor_path_length=int(closure["maximum_shortest_donor_path_length"]),
        role_count=len(relink),
        positive_relink_work=positive_relink,
        recipient_positive_incoming_flux=incoming,
        pair_antisymmetry_residual=antisym,
        row_binding_residual=row,
        total_relink_work_residual=total_residual,
    )
    return {
        "relay_kind": SMOOTH_RELINK_SAME_EVENT_RELAY,
        "certificate": cert,
        "recipient_roles": recipients,
        "backward_donor_closure": tuple(int(x) for x in closure["backward_donor_closure"]),
        "terminal_negative_net_donor_roles": donors,
        "maximum_shortest_donor_path_length": int(closure["maximum_shortest_donor_path_length"]),
        "positive_relink_work": positive_relink,
        "recipient_positive_incoming_flux": incoming,
        "pair_antisymmetry_residual": antisym,
        "row_binding_residual": row,
        "total_relink_work_residual": total_residual,
        "same_physical_event": True,
        "same_physical_time": True,
        "new_causal_charge_created": False,
        "recursive_generation_created": False,
        "scale_progress_created": False,
        "hard_smooth_measure_identification_used": False,
        "primary_selected": False,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "input_type": "only a GaugeQuotientedInterfaceWork carrying its exact K_phys pair matrix is admissible; arbitrary relink atom arrays cannot enter the theorem",
        "measure_separation": "smooth eta_a=A_a^2 roles and hard orthogonal event roles remain distinct physical disintegrations; only the abstract antisymmetric finite-flux lemma is shared",
        "pair_law": "T_ab=-2 Re<eta_a u,K_phys eta_b u> is antisymmetric and its row sums are exactly the signed smooth physical-relink work atoms",
        "divergence": "with F[b->a]=[T_ab]_+, each relink row is incoming minus outgoing positive same-event flux and the total relink work is zero",
        "donor_exhaustion": "the backward positive-inflow closure of all positive-net smooth roles contains a negative-net donor in at most (#roles-1) edges after quotienting cycles",
        "event_semantics": "smooth relink donor tracing stays at one physical event/time and creates no causal charge, recursive generation, or scale progress",
        "strain_separation": "the simultaneous S branch remains existing resolved strain/deformation and is not removed by the skew donor quotient",
        "master_boundary": "conservative_smooth_role_relink is a same-event relay label, not a recursive owner label; a pure-relink interface reentry therefore creates no child RecursiveEventState",
        "scope": "this closes smooth K_phys relink as an independent recursion mechanism; it does not terminate genuine strain, HH, source, dissipation, service, or reuse owners and does not prove Navier-Stokes regularity",
    }


@dataclass(frozen=True)
class SmoothRelinkStress:
    samples: int
    worst_antisymmetry_residual: float
    worst_row_binding_residual: float
    worst_total_relink_residual: float
    minimum_incoming_minus_gain_margin: float
    donor_existence_failures: int
    maximum_shortest_donor_path_length: int
    binding_rejection_failures: int


def _random_certificate(rng: np.random.Generator, m: int) -> GaugeQuotientedInterfaceWork:
    A = rng.normal(size=(m, m))
    T = A - A.T
    # Avoid the measure-zero all-zero pair law and keep pure relink as the owner.
    if float(np.max(np.abs(T))) < 1e-8:
        T[0, 1] = 1.0
        T[1, 0] = -1.0
    relink = T.sum(axis=1)
    if float(np.max(np.abs(relink))) < 1e-8:
        T[:] = 0.0
        T[0, 1] = 1.0
        T[1, 0] = -1.0
        relink = T.sum(axis=1)
    native = relink.copy()
    strain = np.zeros(m, dtype=float)
    return GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=tuple(float(x) for x in native),
        signed_physical_relink_atoms=tuple(float(x) for x in relink),
        signed_existing_strain_atoms=tuple(float(x) for x in strain),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=tuple(tuple(float(x) for x in row) for row in T),
    )


def stress(samples: int = 50_000, seed: int = 20260811) -> SmoothRelinkStress:
    rng = np.random.default_rng(seed)
    wa = wr = wt = 0.0
    minmargin = float("inf")
    donor_fail = bind_fail = 0
    maxpath = 0

    for _ in range(samples):
        m = int(rng.integers(2, 18))
        work = _random_certificate(rng, m)
        try:
            out = smooth_relink_donor_quotient(work)
        except Exception:
            donor_fail += 1
            raise
        wa = max(wa, float(out["pair_antisymmetry_residual"]))
        wr = max(wr, float(out["row_binding_residual"]))
        wt = max(wt, float(out["total_relink_work_residual"]))
        minmargin = min(
            minmargin,
            float(out["recipient_positive_incoming_flux"]) - float(out["positive_relink_work"]),
        )
        maxpath = max(maxpath, int(out["maximum_shortest_donor_path_length"]))
        cert = out["certificate"]
        if not isinstance(cert, SmoothRelinkDonorCertificate):
            raise AssertionError("smooth relink quotient lost typed donor certificate")
        if bool(out["recursive_generation_created"]) or bool(out["new_causal_charge_created"]):
            raise AssertionError("smooth relink donor trace created recursive physics")

        # Pair/row binding is part of the physical certificate and must fail closed.
        rows = [list(row) for row in work.signed_physical_relink_pair_matrix]
        rows[0][1] += 0.25
        bad = GaugeQuotientedInterfaceWork(
            work.signed_native_interface_atoms,
            work.signed_physical_relink_atoms,
            work.signed_existing_strain_atoms,
            work.gauge_transport_operator_residual,
            work.skew_decomposition_residual,
            tuple(tuple(float(x) for x in row) for row in rows),
        )
        try:
            smooth_relink_donor_quotient(bad)
        except ValueError:
            pass
        else:
            bind_fail += 1
            raise AssertionError("unbound smooth relink pair matrix crossed donor quotient")

    return SmoothRelinkStress(samples, wa, wr, wt, minmargin, donor_fail, maxpath, bind_fail)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-smooth-relink-donor-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "smooth_relink_donor_quotient.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Smooth physical relink donor quotient

Status: **{cert['status']}**.

For the already gauge-quotiented smooth square partition `eta_a=A_a^2`, the exact residual-skew pair work is

`T_ab=-2 Re<eta_a u,K_phys eta_b u>`.

The bound work certificate verifies `T_ab=-T_ba` and `R_a=sum_b T_ab`, where `R_a` is the signed physical relink row.  Thus `F[b->a]=[T_ab]_+` is a genuine same-event conservative flux and `R_a=inflow_a-outflow_a`.

Starting from all positive-net relink recipients and closing backward under positive inflow must meet a negative-net donor in finitely many roles.  Internal cycles cancel from subset divergence and create no PDE generation.  Smooth and hard role measures are never identified; only the finite antisymmetric-flux lemma is shared.

Stress: `{out.samples}` bound smooth relink laws
- worst pair antisymmetry residual: `{out.worst_antisymmetry_residual:.3e}`
- worst row-binding residual: `{out.worst_row_binding_residual:.3e}`
- worst total relink residual: `{out.worst_total_relink_residual:.3e}`
- minimum incoming-minus-recipient-gain margin: `{out.minimum_incoming_minus_gain_margin:.3e}`
- donor-existence failures: `{out.donor_existence_failures}`
- maximum shortest donor path: `{out.maximum_shortest_donor_path_length}`
- pair-binding rejection failures: `{out.binding_rejection_failures}`

The master consequence is type-level: `conservative_smooth_role_relink` is same-event provenance, not a recursive generation owner.  A simultaneous strain branch remains genuine resolved strain/deformation.  This theorem does not terminate strain, HH, source, dissipation, service, or reuse recurrence and makes no Navier--Stokes regularity claim.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
