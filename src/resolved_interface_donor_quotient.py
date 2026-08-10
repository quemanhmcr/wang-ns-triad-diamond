from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.nonaffine_role_interface_work import adjoint_split, projector_residual


STATUS = (
    "EXACT_RESOLVED_INTERFACE_DONOR_QUOTIENT__"
    "POSITIVE_INTERFACE_TO_CONSERVATIVE_SKEW_DONOR_OR_EXISTING_STRAIN__"
    "FINITE_SAME_EVENT_DONOR_EXHAUSTION__CIRCULATION_ZERO_RECURSION_DEPTH"
)

SKEW_OWNER = "conservative_skew_role_redistribution"
SYMMETRIC_OWNER = "existing_symmetric_strain_deformation"


def resolved_role_work_decomposition(
    projectors: Sequence[np.ndarray],
    resolved_operator: np.ndarray,
    state: np.ndarray,
) -> dict[str, object]:
    """Exact K+S split of the full resolved low-high role energy work.

    This is the master-facing form needed by high-tail mixed work. Let L be
    the actual resolved linearized RHS operator on the selected high field,
    L=K+S with K*=-K and S*=S, and let {P_a} be a complete orthogonal event-role
    partition. For w_a=P_a u define

        R_ab = 2 Re <w_a, L w_b>.

    Its K pair matrix is antisymmetric and its S pair matrix is symmetric.
    Hence the K row sums are conservative role flux, while every S contribution
    has the provenance of the same physical resolved strain/deformation operator.

    The moving-projector commutator theorem is a companion specialization: its
    skew interface is the same conservative pair flux and its symmetric interface
    is the off-diagonal part of this strain work.
    """
    L = np.asarray(resolved_operator, dtype=complex)
    u = np.asarray(state, dtype=complex)
    Ps = [np.asarray(P, dtype=complex) for P in projectors]
    n = len(u)
    if not Ps or L.shape != (n, n) or any(P.shape != (n, n) for P in Ps):
        raise ValueError("complete matching resolved-role data required")
    total = sum(Ps, np.zeros((n, n), dtype=complex))
    if float(np.linalg.norm(total - np.eye(n))) > 3e-10:
        raise ValueError("event role projectors must resolve identity")
    for i, P in enumerate(Ps):
        if projector_residual(P) > 2e-10:
            raise ValueError("orthogonal event-role projectors required")
        for Q in Ps[i + 1 :]:
            if float(np.linalg.norm(P @ Q)) > 2e-10:
                raise ValueError("event roles must be mutually orthogonal")

    K, S = adjoint_split(L)
    ws = [P @ u for P in Ps]
    m = len(ws)
    TK = np.zeros((m, m), dtype=float)
    DS = np.zeros((m, m), dtype=float)
    for a in range(m):
        for b in range(m):
            TK[a, b] = 2.0 * float(np.real(np.vdot(ws[a], K @ ws[b])))
            DS[a, b] = 2.0 * float(np.real(np.vdot(ws[a], S @ ws[b])))

    total_work = np.array([2.0 * float(np.real(np.vdot(w, L @ u))) for w in ws])
    skew_work = TK.sum(axis=1)
    symmetric_work = DS.sum(axis=1)
    scale = max(
        1.0,
        float(np.max(np.abs(total_work))),
        float(np.max(np.abs(TK))),
        float(np.max(np.abs(DS))),
    )
    split_residual = float(np.max(np.abs(total_work - skew_work - symmetric_work)))
    skew_residual = float(np.max(np.abs(TK + TK.T)))
    symmetric_residual = float(np.max(np.abs(DS - DS.T)))
    if max(split_residual, skew_residual, symmetric_residual) > 8e-11 * scale:
        raise AssertionError("resolved role work lost its exact adjoint K+S structure")
    return {
        "signed_resolved_role_work": total_work,
        "signed_skew_role_work": skew_work,
        "signed_symmetric_strain_work": symmetric_work,
        "skew_pair_matrix": TK,
        "symmetric_pair_matrix": DS,
        "work_split_residual": split_residual,
        "skew_antisymmetry_residual": skew_residual,
        "symmetric_symmetry_residual": symmetric_residual,
        "total_skew_work": float(skew_work.sum()),
        "symmetric_is_existing_strain_provenance": True,
        "new_interface_source_created": False,
    }


def positive_interface_component_split(
    signed_interface_atoms: Sequence[float],
    signed_skew_atoms: Sequence[float],
    signed_symmetric_atoms: Sequence[float],
) -> dict[str, object]:
    """Split positive resolved cross/interface work without inventing a new currency.

    At each physical work atom the exact adjoint decomposition is

        I = I_K + I_S,

    with K*=-K and S*=S.  Therefore [I]_+ <= [I_K]_+ + [I_S]_+.
    Summing any event/time disintegration of the *same physical work law* gives
    the corresponding positive-measure cover.  If positive interface work is
    nonzero, skew redistribution or symmetric strain carries at least one half;
    exact ties remain joint.
    """
    total = np.asarray(tuple(float(x) for x in signed_interface_atoms), dtype=float)
    skew = np.asarray(tuple(float(x) for x in signed_skew_atoms), dtype=float)
    sym = np.asarray(tuple(float(x) for x in signed_symmetric_atoms), dtype=float)
    if total.ndim != 1 or len(total) == 0 or skew.shape != total.shape or sym.shape != total.shape:
        raise ValueError("matching nonempty one-dimensional interface work atoms required")
    if np.any(~np.isfinite(total)) or np.any(~np.isfinite(skew)) or np.any(~np.isfinite(sym)):
        raise ValueError("finite interface work atoms required")

    identity_residual = float(np.max(np.abs(total - skew - sym)))
    identity_scale = max(1.0, float(np.max(np.abs(total))), float(np.max(np.abs(skew))), float(np.max(np.abs(sym))))
    if identity_residual > 6e-12 * identity_scale:
        raise ValueError("signed interface atoms do not realize I=I_K+I_S")

    W = float(np.maximum(total, 0.0).sum())
    WK = float(np.maximum(skew, 0.0).sum())
    WS = float(np.maximum(sym, 0.0).sum())
    cover_margin = WK + WS - W
    tol = 8e-13 * max(1.0, W, WK, WS)
    if cover_margin < -tol:
        raise AssertionError("positive interface work escaped skew+strain positive cover")

    owners: list[str] = []
    threshold = 0.5 * W
    if W > 0:
        if WK + tol >= threshold:
            owners.append(SKEW_OWNER)
        if WS + tol >= threshold:
            owners.append(SYMMETRIC_OWNER)
        if not owners:
            raise AssertionError("positive interface half-pigeonhole lost both physical continuations")

    return {
        "positive_interface_work": W,
        "positive_skew_redistribution_work": WK,
        "positive_symmetric_strain_work": WS,
        "owner_threshold": threshold,
        "joint_owners": tuple(owners),
        "signed_identity_max_residual": identity_residual,
        "positive_cover_margin": cover_margin,
        "new_source_currency_created": False,
        "primary_selected": False,
    }


def _validated_skew_pair_matrix(skew_pair_matrix: np.ndarray) -> np.ndarray:
    T = np.asarray(skew_pair_matrix, dtype=float)
    if T.ndim != 2 or T.shape[0] != T.shape[1] or T.shape[0] < 2:
        raise ValueError("square skew pair matrix of size at least two required")
    if np.any(~np.isfinite(T)):
        raise ValueError("finite skew pair matrix required")
    scale = max(1.0, float(np.max(np.abs(T))))
    if float(np.max(np.abs(T + T.T))) > 8e-12 * scale:
        raise ValueError("resolved skew pair work must be antisymmetric")
    if float(np.max(np.abs(np.diag(T)))) > 8e-12 * scale:
        raise ValueError("skew role self-transfer must vanish")
    return T


def skew_flux_balance(skew_pair_matrix: np.ndarray) -> dict[str, object]:
    """Discrete divergence form of conservative resolved role transfer.

    For any physical signed-work convention, the skew pair matrix obeys

        T_ab = -T_ba,
        R_a^K = sum_b T_ab.

    Read F[b->a]=[T_ab]_+.  Then exactly

        R_a^K = incoming_a - outgoing_a,
        sum_a R_a^K = 0.

    Thus positive skew interface work is a gain from other roles, not generation.
    """
    T = _validated_skew_pair_matrix(skew_pair_matrix)
    F = np.maximum(T, 0.0)
    net = T.sum(axis=1)
    incoming = F.sum(axis=1)
    outgoing = F.sum(axis=0)
    residual = net - (incoming - outgoing)
    return {
        "positive_flux_matrix_recipient_by_donor": F,
        "net_skew_role_work": net,
        "incoming_positive_flux": incoming,
        "outgoing_positive_flux": outgoing,
        "worst_role_divergence_residual": float(np.max(np.abs(residual))),
        "total_net_skew_work": float(net.sum()),
        "total_positive_directed_flux": float(F.sum()),
    }


def skew_subset_balance(skew_pair_matrix: np.ndarray, roles: Sequence[int]) -> dict[str, float]:
    """Exact boundary-flux identity for any role set.

    Internal role circulation cancels pairwise.  Only flux through the boundary
    can change the total skew work of the selected role set.
    """
    T = _validated_skew_pair_matrix(skew_pair_matrix)
    m = T.shape[0]
    C = tuple(sorted({int(x) for x in roles}))
    if not C or any(x < 0 or x >= m for x in C):
        raise ValueError("nonempty valid role subset required")
    inside = np.zeros(m, dtype=bool)
    inside[list(C)] = True
    outside = ~inside
    net = T.sum(axis=1)
    subset_net = float(net[inside].sum())
    boundary_signed = float(T[np.ix_(inside, outside)].sum())
    boundary_in = float(np.maximum(T[np.ix_(inside, outside)], 0.0).sum())
    boundary_out = float(np.maximum(T[np.ix_(outside, inside)], 0.0).sum())
    return {
        "subset_net_skew_work": subset_net,
        "boundary_signed_flux": boundary_signed,
        "boundary_positive_inflow": boundary_in,
        "boundary_positive_outflow": boundary_out,
        "boundary_balance_residual": subset_net - (boundary_in - boundary_out),
        "internal_circulation_contribution": 0.0,
    }


def skew_donor_closure(skew_pair_matrix: np.ndarray, recipient_roles: Sequence[int]) -> dict[str, object]:
    """Trace any positive skew gain to actual donating roles at the same event.

    Starting from one or several roles with strictly positive net skew work, add
    every role that sends positive skew flux into the current set.  The resulting
    backward donor closure is finite.  It must contain a role with negative net
    skew work.  Otherwise the closure has strictly positive total net work, while
    by construction it has no positive boundary inflow, contradicting the exact
    subset divergence identity.

    No donor is selected as primary.  All reachable negative-net donors are
    returned.  Cycles may occur inside the closure, but they are circulation and
    contribute exactly zero to the subset balance.
    """
    T = _validated_skew_pair_matrix(skew_pair_matrix)
    m = T.shape[0]
    starts = tuple(sorted({int(x) for x in recipient_roles}))
    if not starts or any(x < 0 or x >= m for x in starts):
        raise ValueError("nonempty valid recipient role set required")
    net = T.sum(axis=1)
    scale = max(1.0, float(np.max(np.abs(T))), float(np.max(np.abs(net))))
    tol = 32.0 * math.ulp(scale)
    if any(float(net[a]) <= tol for a in starts):
        raise ValueError("every starting recipient must have strictly positive net skew work")

    closure = set(starts)
    distance = {a: 0 for a in starts}
    frontier = list(starts)
    while frontier:
        a = frontier.pop(0)
        donors = [b for b in range(m) if float(T[a, b]) > tol]
        for b in donors:
            if b not in closure:
                closure.add(b)
                distance[b] = distance[a] + 1
                frontier.append(b)

    C = tuple(sorted(closure))
    terminal = tuple(a for a in C if float(net[a]) < -tol)
    if not terminal:
        balance = skew_subset_balance(T, C)
        raise AssertionError(
            "positive skew recipient donor-closure contains no negative-net donor; "
            f"subset_net={balance['subset_net_skew_work']}, boundary_in={balance['boundary_positive_inflow']}"
        )

    balance = skew_subset_balance(T, C)
    b_in = float(balance["boundary_positive_inflow"])
    if b_in > 8e-12 * scale:
        raise AssertionError("backward donor closure still has positive skew inflow from outside")
    max_shortest = max(distance[a] for a in terminal)
    if max_shortest > m - 1:
        raise AssertionError("finite donor trace exceeded simple-path role bound")

    start_gain = float(net[list(starts)].sum())
    start_incoming = float(np.maximum(T[np.ix_(list(starts), range(m))], 0.0).sum())
    if start_incoming + 8e-12 * scale < start_gain:
        raise AssertionError("positive role gain exceeds its actual incoming donor flux")

    return {
        "recipient_roles": starts,
        "backward_donor_closure": C,
        "terminal_negative_net_donor_roles": terminal,
        "shortest_donor_path_lengths": tuple((a, int(distance[a])) for a in terminal),
        "maximum_shortest_donor_path_length": int(max_shortest),
        "simple_path_length_upper": m - 1,
        "recipient_net_gain": start_gain,
        "recipient_positive_incoming_flux": start_incoming,
        "closure_subset_net_skew_work": float(balance["subset_net_skew_work"]),
        "closure_boundary_positive_inflow": b_in,
        "closure_boundary_positive_outflow": float(balance["boundary_positive_outflow"]),
        "closure_balance_residual": float(balance["boundary_balance_residual"]),
        "same_physical_event": True,
        "same_physical_time": True,
        "new_causal_charge_created": False,
        "recursive_generation_created": False,
        "scale_progress_created": False,
        "primary_selected": False,
    }


def high_tail_interface_component_route(
    physical_tail_dissipation_lower: float,
    viscosity: float,
    positive_interface_common_work: float,
    positive_skew_common_work: float,
    positive_symmetric_common_work: float,
) -> dict[str, object]:
    """High-tail interface owner -> skew donor or existing strain, in common N dW.

    If high-tail regeneration names the interface owner, it carries at least
    nu D_tail/2 in the parent block's common physical work unit.  The exact
    K+S interface split then forces skew redistribution or symmetric strain to
    carry at least one half of the realized interface work and hence, cleanly,
    at least nu D_tail/4.  No M/N reweighting or new probability law appears.
    """
    D = float(physical_tail_dissipation_lower)
    nu = float(viscosity)
    W = float(positive_interface_common_work)
    WK = float(positive_skew_common_work)
    WS = float(positive_symmetric_common_work)
    if min(D, nu) <= 0 or min(W, WK, WS) < 0 or not all(math.isfinite(x) for x in (D, nu, W, WK, WS)):
        raise ValueError("positive finite tail data and nonnegative common interface work required")
    clean_interface = 0.5 * nu * D
    tol = 8e-13 * max(1.0, clean_interface, W, WK, WS)
    if W + tol < clean_interface:
        raise ValueError("supplied interface branch is not a clean high-tail interface owner")
    if WK + WS + tol < W:
        raise ValueError("skew+symmetric positive common work does not cover interface positive work")

    threshold = 0.5 * W
    clean_component = 0.25 * nu * D
    owners: list[str] = []
    if WK + tol >= threshold:
        owners.append(SKEW_OWNER)
    if WS + tol >= threshold:
        owners.append(SYMMETRIC_OWNER)
    if not owners:
        raise AssertionError("clean high-tail interface owner lost both native interface continuations")
    if threshold + tol < clean_component:
        raise AssertionError("interface component threshold lost nu D_tail/4 clean lower")

    return {
        "positive_interface_common_work": W,
        "positive_skew_common_work": WK,
        "positive_symmetric_common_work": WS,
        "owner_threshold": threshold,
        "clean_interface_owner_lower": clean_interface,
        "clean_component_owner_lower": clean_component,
        "joint_owners": tuple(owners),
        "common_work_unit": "parent_block_frequency_times_physical_energy_work",
        "shell_scale_reweighting_used": False,
        "new_causal_probability_created": False,
        "primary_selected": False,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "operator_bridge": "full high-tail resolved mixed work uses the actual linearized low-high operator L_V on the event role partition; the moving-projector commutator is a companion observable with the same K/S provenance, not an identified measure",
        "signed_split": "for the actual resolved linearized low-high operator L=K+S with K*=-K and S*=S, every signed work atom satisfies R=R_K+R_S, hence [R]_+<=[R_K]_+ + [R_S]_+",
        "positive_owner_split": "positive resolved cross/interface work is covered by conservative skew redistribution and existing symmetric strain/deformation work; one carries at least half, exact ties joint",
        "skew_divergence": "for the signed resolved-work convention T_ab=-T_ba and R_a^K=sum_b T_ab; with F[b->a]=[T_ab]_+, R_a^K=inflow_a-outflow_a and sum_a R_a^K=0; an overall PDE sign only reverses edge orientation",
        "subset_identity": "for every role set C, sum_(a in C) R_a^K equals positive boundary inflow minus positive boundary outflow; all internal circulation cancels exactly",
        "donor_exhaustion": "the backward closure of any positive-gain role under positive incoming skew flux contains a negative-net donor; otherwise it has positive net gain but no positive external inflow, contradicting the subset identity",
        "finite_same_event": "after removing cycles, a donor path has at most (#roles-1) edges and stays at the same physical event time; role-transfer depth is not recursive PDE depth",
        "symmetric_semantics": "the S-branch is the existing physical symmetric strain/deformation work; the moving-commutator interface is its off-diagonal specialization; delegate once to coherent-deformation/high-strain/objective-source/critical-D_V owners, never to a new source currency",
        "skew_semantics": "the K-branch is conservative donor tracing; it creates neither energy nor a second causal charge, and internal role circulation is quotiented before master recursion",
        "high_tail": "a clean high-tail interface owner >=nu D_tail/2 therefore gives conservative skew donor work or existing symmetric strain work >=nu D_tail/4 in the same common N dW unit",
        "anti_reweight": "no role index, theorem name, shell-dependent M/N factor, or normalized interface fraction is allowed to redefine the physical causal law",
        "scale": "same-event interface donor tracing supplies no scale progress; any later scale motion must come from the physical owner reached after the quotient",
        "scope": "this closes resolved interface as an independent recursive-generation loophole; it does not prove that the reached donor/strain owner terminates globally and it makes no Navier-Stokes regularity claim",
    }


@dataclass(frozen=True)
class InterfaceDonorStress:
    samples: int
    worst_signed_split_residual: float
    minimum_positive_cover_margin: float
    worst_role_divergence_residual: float
    worst_total_skew_residual: float
    worst_closure_balance_residual: float
    minimum_recipient_incoming_margin: float
    donor_existence_failures: int
    maximum_shortest_donor_path_length: int
    high_tail_component_failures: int


def stress(samples: int = 50_000, seed: int = 20260810) -> InterfaceDonorStress:
    rng = np.random.default_rng(seed)
    wsplit = wdiv = wtotal = wclosure = 0.0
    mcover = mincoming = float("inf")
    donor_fail = high_fail = 0
    maxpath = 0

    for _ in range(samples):
        # Signed K+S work atoms: cancellation may make positive component mass
        # strictly larger than the positive total interface work.
        n = int(rng.integers(1, 30))
        k = rng.normal(size=n)
        s = rng.normal(size=n)
        total = k + s
        split = positive_interface_component_split(total, k, s)
        wsplit = max(wsplit, float(split["signed_identity_max_residual"]))
        mcover = min(mcover, float(split["positive_cover_margin"]))

        # Any antisymmetric pair-work matrix obeys the same finite donor theorem
        # as the one produced by the exact nonaffine role partition identity.
        m = int(rng.integers(2, 14))
        A = rng.normal(size=(m, m))
        T = A - A.T
        bal = skew_flux_balance(T)
        scale = max(1.0, float(np.max(np.abs(T))))
        wdiv = max(wdiv, float(bal["worst_role_divergence_residual"]) / scale)
        wtotal = max(wtotal, abs(float(bal["total_net_skew_work"])) / scale)
        net = np.asarray(bal["net_skew_role_work"], dtype=float)
        positives = [a for a in range(m) if net[a] > 64.0 * math.ulp(max(1.0, scale))]
        if positives:
            a = int(positives[int(rng.integers(0, len(positives)))])
            try:
                donor = skew_donor_closure(T, (a,))
            except AssertionError:
                donor_fail += 1
                continue
            wclosure = max(wclosure, abs(float(donor["closure_balance_residual"])) / scale)
            mincoming = min(
                mincoming,
                float(donor["recipient_positive_incoming_flux"]) - float(donor["recipient_net_gain"]),
            )
            maxpath = max(maxpath, int(donor["maximum_shortest_donor_path_length"]))
            if not tuple(donor["terminal_negative_net_donor_roles"]):
                donor_fail += 1

        D = float(math.exp(rng.uniform(-4.0, 3.0)))
        nu = float(math.exp(rng.uniform(-3.0, 1.0)))
        clean = 0.5 * nu * D
        W = clean * float(rng.uniform(1.0, 3.0))
        frac = float(rng.uniform(0.0, 1.0))
        extra = W * float(rng.uniform(0.0, 1.0))
        WK = frac * (W + extra)
        WS = (1.0 - frac) * (W + extra)
        route = high_tail_interface_component_route(D, nu, W, WK, WS)
        if not tuple(route["joint_owners"]):
            high_fail += 1

    if mincoming == float("inf"):
        mincoming = 0.0
    return InterfaceDonorStress(
        samples,
        wsplit,
        mcover,
        wdiv,
        wtotal,
        wclosure,
        mincoming,
        donor_fail,
        maxpath,
        high_fail,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-resolved-interface-donor-quotient"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "resolved_interface_donor_quotient.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Resolved interface donor/circulation quotient

Status: **{cert['status']}**.

The resolved low--high cross/interface work is not an independent energy source.  Write its actual physical linearized operator as `L=K+S`, with `K*=-K` and `S*=S`.  Event by event,

`R = R_K + R_S`,

so positive work obeys

`[R]_+ <= [R_K]_+ + [R_S]_+`.

Therefore positive resolved cross/interface work has only two native continuations: conservative skew role redistribution, or the already existing symmetric strain/deformation work.  One carries at least half of the positive interface law; exact ties remain joint.  The symmetric branch is not charged again: it delegates to the existing coherent-deformation / high-strain / objective-source / critical-`D_V` owners.

For a complete orthogonal event-role partition, the skew pair work is

`T_ab=2 Re <w_a,K w_b>=-T_ba`,

and the net skew work of role `a` is `R_a^K=sum_b T_ab`.  With directed physical flux `F[b->a]=[T_ab]_+`,

`R_a^K = incoming_a - outgoing_a`,

and `sum_a R_a^K=0`.  More generally, for every role set `C`, all internal transfers cancel and

`sum_(a in C) R_a^K = boundary inflow(C) - boundary outflow(C)`.

This immediately gives a finite donor theorem.  Start from any role with positive skew gain and close the set backward under every positive incoming donor edge.  If every role in that closure had nonnegative net skew work, the closure would have strictly positive total gain.  But by construction no positive flux enters it from outside, so its exact boundary balance is nonpositive.  Contradiction.  Hence the closure contains an actual negative-net donor role.  Removing cycles gives a donor path of at most `#roles-1` edges.

The whole donor trace occurs at the **same physical event time**.  A role cycle is circulation, not a new Navier--Stokes generation.  It creates no extra causal charge and no master recursion depth.  No donor is selected by role index; all reachable negative-net donors remain a set-valued physical provenance mark.

For the high-tail interface owner this yields, in the unchanged common `N dW` unit,

`W_interface^+ >= nu D_tail/2`

implies

`W_skew^+ >= nu D_tail/4`

or

`W_symmetric^+ >= nu D_tail/4`,

with ties joint.  The first is same-event donor tracing; the second is existing strain/deformation provenance.  No `M/N` reweighting is introduced.

Stress: `{out.samples}` split/flux/donor/high-tail states
- worst signed `R=R_K+R_S` residual: `{out.worst_signed_split_residual:.3e}`
- minimum positive-cover margin: `{out.minimum_positive_cover_margin:.3e}`
- worst role-divergence residual: `{out.worst_role_divergence_residual:.3e}`
- worst total skew-work residual: `{out.worst_total_skew_residual:.3e}`
- worst donor-closure balance residual: `{out.worst_closure_balance_residual:.3e}`
- minimum recipient incoming-flux margin: `{out.minimum_recipient_incoming_margin:.3e}`
- donor-existence failures: `{out.donor_existence_failures}`
- maximum sampled shortest donor path: `{out.maximum_shortest_donor_path_length}`
- high-tail component failures: `{out.high_tail_component_failures}`

This theorem closes **resolved interface work as an independent recursive-generation loophole**.  It does not claim that the eventual donor or strain owner globally terminates, and it makes no claim of Navier--Stokes global regularity.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
