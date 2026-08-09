from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.material_coherent_labels import intrinsic_zeta
from src.physical_energy_causal_bridge import route_physical_energy_causality


def registered_coefficient_energy_lower(coefficient: complex, probe_l2: float = 1.0) -> float:
    """Cauchy lower bound for the energy of the same smooth carrier.

    If z=<w,psi>, w=Q(t,D)u and ||psi||_2=probe_l2, then
        ||w||_2^2 >= |z|^2 / ||psi||_2^2.
    Nothing is projected to a new hard role at the common slice.
    """
    p = float(probe_l2)
    if p <= 0 or not math.isfinite(p):
        raise ValueError("positive finite probe norm required")
    return abs(complex(coefficient)) ** 2 / (p * p)


def coefficient_energy_margin(carrier: np.ndarray, probe: np.ndarray) -> float:
    w = np.asarray(carrier, complex)
    psi = np.asarray(probe, complex)
    if w.ndim != 1 or psi.shape != w.shape:
        raise ValueError("matching one-dimensional carrier/probe required")
    p = float(np.linalg.norm(psi))
    if p <= 0:
        raise ValueError("nonzero probe required")
    z = complex(np.vdot(psi, w))
    return float(np.vdot(w, w).real) - registered_coefficient_energy_lower(z, p)


def two_step_affine_reanchor_zeta_residual(
    M1: np.ndarray,
    M2: np.ndarray,
    L: np.ndarray,
    X: np.ndarray,
    k: np.ndarray,
) -> float:
    """Re-anchoring a material carrier is composition of the same affine gauge.

    First apply M1 and then M2 to physical space, with the dual Kelvin update
    k -> M^{-T}k.  The intrinsic zeta after the two steps equals both the
    original label and the one-step update by M2 M1.
    """
    M1 = np.asarray(M1, float)
    M2 = np.asarray(M2, float)
    L = np.asarray(L, float)
    X = np.asarray(X, float)
    k = np.asarray(k, float)
    for M in (M1, M2, L):
        if M.shape != (3, 3):
            raise ValueError("3x3 affine matrices required")
        if abs(float(np.linalg.det(M))) < 1e-12:
            raise ValueError("invertible affine matrices required")
    if X.shape != (3,) or k.shape != (3,):
        raise ValueError("X,k must be three-vectors")

    z0 = intrinsic_zeta(L, X, k)
    L1 = M1 @ L
    X1 = M1 @ X
    k1 = np.linalg.solve(M1.T, k)
    z1 = intrinsic_zeta(L1, X1, k1)

    L2 = M2 @ L1
    X2 = M2 @ X1
    k2 = np.linalg.solve(M2.T, k1)
    z2 = intrinsic_zeta(L2, X2, k2)

    M = M2 @ M1
    zdirect = intrinsic_zeta(M @ L, M @ X, np.linalg.solve(M.T, k))
    return float(max(np.linalg.norm(z1 - z0), np.linalg.norm(z2 - z0), np.linalg.norm(z2 - zdirect)))


def hahn_weighted_generation_chain(signed_work_density: Sequence[float], carrier_symbol: Sequence[float]) -> dict[str, float]:
    """Positive smooth-carrier work is dominated by the same physical Hahn law.

    For the signed physical HH work density r and a scalar Fourier carrier
    0<=q<=1,

        [integral q^2 r]_+ <= integral q^2 [r]_+ <= integral [r]_+.

    The first inequality is Hahn/triangle algebra; the second is contraction.
    Hence Q need not define a new causal probability measure.  The fine physical
    positive work measure already dominates every positive generation of Q.
    """
    r = np.asarray(signed_work_density, float)
    q = np.asarray(carrier_symbol, float)
    if r.ndim != 1 or q.shape != r.shape or not np.all(np.isfinite(r)) or not np.all(np.isfinite(q)):
        raise ValueError("matching finite one-dimensional work/symbol data required")
    if np.any(q < 0) or np.any(q > 1):
        raise ValueError("carrier multiplier must lie in [0,1]")
    weighted_signed = float(np.dot(q * q, r))
    aggregate_positive = max(0.0, weighted_signed)
    weighted_hahn_positive = float(np.dot(q * q, np.maximum(r, 0.0)))
    physical_hahn_positive = float(np.maximum(r, 0.0).sum())
    return {
        "weighted_signed_work": weighted_signed,
        "carrier_positive_generation": aggregate_positive,
        "weighted_hahn_positive": weighted_hahn_positive,
        "physical_hahn_positive": physical_hahn_positive,
        "hahn_margin": weighted_hahn_positive - aggregate_positive,
        "contraction_margin": physical_hahn_positive - weighted_hahn_positive,
    }


def orthogonal_hard_event_resolution(
    carrier: np.ndarray,
    forcing_atoms: Sequence[np.ndarray],
    hard_projectors: Sequence[np.ndarray],
) -> dict[str, float]:
    """Exact hard resolution is needed only when actual work is read at an event.

    The projectors are an orthogonal resolution of identity.  For F=sum_alpha F_a,

      2 Re<w,F> = sum_{C,alpha} 2 Re<P_C w,P_C F_alpha>.

    This is event algebra; no P_C is differentiated or propagated between events.
    """
    w = np.asarray(carrier, complex)
    if w.ndim != 1:
        raise ValueError("one-dimensional carrier vector required")
    n = len(w)
    Fs = tuple(np.asarray(F, complex) for F in forcing_atoms)
    Ps = tuple(np.asarray(P, complex) for P in hard_projectors)
    if not Fs or not Ps or any(F.shape != (n,) for F in Fs) or any(P.shape != (n, n) for P in Ps):
        raise ValueError("nonempty matching forcing atoms/projectors required")
    I = np.eye(n, dtype=complex)
    if np.linalg.norm(sum(Ps, np.zeros((n, n), complex)) - I) > 1e-10:
        raise ValueError("hard projectors must resolve identity")
    for i, P in enumerate(Ps):
        if np.linalg.norm(P - P.conj().T) > 1e-10 or np.linalg.norm(P @ P - P) > 1e-10:
            raise ValueError("self-adjoint projections required")
        for j, R in enumerate(Ps):
            if i != j and np.linalg.norm(P @ R) > 1e-10:
                raise ValueError("orthogonal hard projections required")
    Ftot = sum(Fs, np.zeros(n, complex))
    exact = 2.0 * float(np.real(np.vdot(w, Ftot)))
    atom_sum = 0.0
    positive_atoms = 0.0
    for P in Ps:
        Pw = P @ w
        for F in Fs:
            a = 2.0 * float(np.real(np.vdot(Pw, P @ F)))
            atom_sum += a
            positive_atoms += max(0.0, a)
    return {
        "signed_work": exact,
        "hard_atom_sum": atom_sum,
        "reconstruction_residual": atom_sum - exact,
        "hard_positive_mass": positive_atoms,
        "aggregate_positive_work": max(0.0, exact),
        "positive_dominance_margin": positive_atoms - max(0.0, exact),
    }


def relay_energy_route(
    *,
    terminal_carrier_energy: float,
    initial_carrier_energy: float,
    residual_positive_work: float,
    strain_action: float,
) -> dict[str, float | str]:
    """Apply the existing energy gate to the same smooth carrier on the previous slab."""
    return route_physical_energy_causality(
        terminal_energy=terminal_carrier_energy,
        initial_energy=initial_carrier_energy,
        residual_positive_work=residual_positive_work,
        strain_action=strain_action,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_COMMON_SLICE_SMOOTH_MATERIAL_CARRIER_RELAY__HARDEN_ONLY_AT_PHYSICAL_GENERATION__UNIVERSAL_SERVICE_REENTRY_REMAINS",
        "coefficient_to_energy": "z=<Q u,psi> implies ||Q u||_2^2>=|z|^2/||psi||_2^2; no hard role is needed at the common slice",
        "reanchor": "re-anchoring the same carrier is composition of common affine/Kelvin transport; intrinsic zeta is exactly unchanged, so no relink or new Xi is created",
        "previous_slab": "the existing outer moving-role equation and physical energy gate apply to the same smooth Q carrier whenever the previous slab satisfies their service/scale hypotheses",
        "generation_hahn": "for signed physical HH density r and 0<=q<=1, [int q^2 r]_+ <= int q^2[r]_+ <= int[r]_+; smooth-carrier generation is dominated by actual positive physical HH work",
        "event_hardening": "only at an actual HH interaction is the physical work resolved into orthogonal hard Fourier/helical roles; exact hard work reconstruction needs no persistent hard projector",
        "causal_law": "Q is never promoted to a new transfer probability law; actual positive physical HH work remains the causal law",
        "representation": "hard-cell symbol freezing at the actual event reuses the existing summable physical-multiplier Xi; relay/re-anchoring itself costs zero Xi",
        "scope": "this removes intermediate hard-event reselection from recursion architecture but does not prove that every renewed previous slab satisfies low-strain/service/moat/natural-lifetime hypotheses",
    }


@dataclass(frozen=True)
class CarrierRelayStress:
    samples: int
    minimum_coefficient_energy_margin: float
    worst_two_step_zeta_residual: float
    minimum_hahn_margin: float
    minimum_hahn_contraction_margin: float
    worst_hard_resolution_residual: float
    minimum_hard_positive_dominance_margin: float
    generation_routes: int


def _coordinate_projectors(n: int, groups: Sequence[Sequence[int]]) -> tuple[np.ndarray, ...]:
    out = []
    seen: set[int] = set()
    for group in groups:
        P = np.zeros((n, n), complex)
        for j in group:
            j = int(j)
            if j < 0 or j >= n or j in seen:
                raise ValueError("groups must partition coordinates")
            seen.add(j)
            P[j, j] = 1.0
        out.append(P)
    if seen != set(range(n)):
        raise ValueError("groups must cover all coordinates")
    return tuple(out)


def stress(samples: int = 50_000, seed: int = 20260809) -> CarrierRelayStress:
    rng = np.random.default_rng(seed)
    mce = mh = mcon = mhard = float("inf")
    wz = wr = 0.0
    generated = 0
    for _ in range(samples):
        n = int(rng.integers(2, 12))
        w = rng.normal(size=n) + 1j * rng.normal(size=n)
        psi = rng.normal(size=n) + 1j * rng.normal(size=n)
        cm = coefficient_energy_margin(w, psi)
        mce = min(mce, cm)
        if cm < -3e-12 * max(1.0, float(np.vdot(w, w).real)):
            raise AssertionError("registered coefficient exceeded same-carrier energy")

        def invmat() -> np.ndarray:
            A = np.eye(3) + 0.12 * rng.normal(size=(3, 3))
            if abs(float(np.linalg.det(A))) < 0.2:
                A += np.eye(3)
            return A
        M1, M2, L = invmat(), invmat(), invmat()
        X, k = rng.normal(size=3), rng.normal(size=3)
        zr = two_step_affine_reanchor_zeta_residual(M1, M2, L, X, k)
        wz = max(wz, zr)
        if zr > 5e-11:
            raise AssertionError("common-affine reanchor changed intrinsic material label")

        m = int(rng.integers(2, 100))
        r = rng.normal(size=m) * rng.lognormal(mean=0.0, sigma=1.0, size=m)
        q = rng.random(m)
        hh = hahn_weighted_generation_chain(r, q)
        mh = min(mh, float(hh["hahn_margin"]))
        mcon = min(mcon, float(hh["contraction_margin"]))
        if float(hh["hahn_margin"]) < -3e-12 * max(1.0, float(hh["physical_hahn_positive"])):
            raise AssertionError("carrier aggregate positive work exceeded weighted Hahn positive mass")
        if float(hh["contraction_margin"]) < -3e-12 * max(1.0, float(hh["physical_hahn_positive"])):
            raise AssertionError("carrier contraction created physical positive work")

        # Hard roles appear only at this synthetic work event.
        nf = int(rng.integers(1, 6))
        Fs = [rng.normal(size=n) + 1j * rng.normal(size=n) for _ in range(nf)]
        perm = rng.permutation(n)
        cuts = sorted(set([0, n] + [int(x) for x in rng.integers(1, n, size=int(rng.integers(0, min(5, n))))]))
        groups = [perm[cuts[j]:cuts[j+1]].tolist() for j in range(len(cuts)-1) if cuts[j] < cuts[j+1]]
        Ps = _coordinate_projectors(n, groups)
        ho = orthogonal_hard_event_resolution(w, Fs, Ps)
        rr = abs(float(ho["reconstruction_residual"]))
        wr = max(wr, rr)
        mhard = min(mhard, float(ho["positive_dominance_margin"]))
        if rr > 5e-12 * max(1.0, abs(float(ho["signed_work"]))):
            raise AssertionError("event-only hard resolution lost signed work")
        if float(ho["positive_dominance_margin"]) < -5e-12 * max(1.0, float(ho["hard_positive_mass"])):
            raise AssertionError("hard positive event mass failed to dominate aggregate work")

        E1 = float(rng.lognormal(mean=0.0, sigma=1.0))
        route = relay_energy_route(
            terminal_carrier_energy=E1,
            initial_carrier_energy=float(rng.uniform(0.0, 0.19)) * E1,
            residual_positive_work=float(rng.uniform(0.0, 0.19)) * E1,
            strain_action=float(rng.uniform(0.0, 1.0 / 30.0)),
        )
        if route["branch"] == "physical_high_high_transfer_generation":
            generated += 1

    return CarrierRelayStress(samples, mce, wz, mh, mcon, wr, mhard, generated)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-smooth-material-carrier-relay"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "smooth_material_carrier_relay.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Smooth material-carrier relay: hard roles only at physical events\n\nStatus: **{cert['status']}**.\n\nThe common slice does not need a fresh hard packet.  If the registered parent is the smooth moving role `w=Q(t,D)u` and its adjoint coefficient is `z=<w,psi>`, then simply\n\n`||w||_2^2 >= |z|^2/||psi||_2^2`.\n\nThus the registered coefficient already supplies nonzero terminal energy for the **same physical material carrier**.  Re-anchor that carrier at the common slice, rather than choosing a new maximizing Fourier/coherent cell.  Common affine/Kelvin maps compose, and the intrinsic coordinate `zeta=(L^-1 X/2,L^T k)` is exactly invariant under the composition.  Re-anchoring therefore changes coordinates, not material identity: it carries zero relink mass and zero new `Xi`.\n\nOn a previous slab for which the existing smooth-role/service hypotheses hold, apply the same outer-role equation and physical-energy gate to `Q u`.  If inherited energy, classified residual work or critical strain fires, the existing cause owns the block.  Otherwise the carrier has definite positive HH generation.  This does **not** make `Q` a new physical transfer law.  If `r` is the signed physical HH work density and `0<=q<=1` is the scalar carrier multiplier,\n\n`[ integral q^2 r ]_+ <= integral q^2 [r]_+ <= integral [r]_+`.\n\nSo positive generation of the smooth carrier is dominated by actual positive physical HH work on the same event support.  Only **at that actual nonlinear event** do we resolve the work into orthogonal hard Fourier/helical roles.  Exact orthogonal work reconstruction then gives the parent/child event labels; no hard boundary is ever differentiated or assumed to persist between events.  Cellwise freezing of the physical multiplier at that event is the already existing summable symbol `Xi`, not a relay cost.\n\nThe natural architecture is therefore:\n\n`hard physical interaction event -> smooth material carrier -> common-slice relay -> smooth material carrier -> next actual HH event -> hard physical interaction event`.\n\nStress: `{out.samples}` carrier/affine/Hahn/hard-event states\n- minimum coefficient-to-carrier-energy margin: `{out.minimum_coefficient_energy_margin:.3e}`\n- worst two-step intrinsic-zeta residual: `{out.worst_two_step_zeta_residual:.3e}`\n- minimum aggregate-to-Hahn margin: `{out.minimum_hahn_margin:.3e}`\n- minimum carrier-to-physical positive-work contraction margin: `{out.minimum_hahn_contraction_margin:.3e}`\n- worst event hard-resolution residual: `{out.worst_hard_resolution_residual:.3e}`\n- minimum hard-positive dominance margin: `{out.minimum_hard_positive_dominance_margin:.3e}`\n- low-strain generation routes sampled: `{out.generation_routes}`\n\nThis removes **intermediate hard-event reselection** as a necessary recursive object.  It does not yet prove universal recursion.  The remaining continuum question is whether every renewed backward slab of the relayed carrier satisfies, or physically exits through, the already named scale/moat/service/natural-lifetime hypotheses.  In particular the relay may not be extended through an arbitrary interval by fiat.  No Navier--Stokes global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
