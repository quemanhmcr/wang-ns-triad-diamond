from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.material_coherent_labels import intrinsic_zeta


@dataclass(frozen=True)
class IntrinsicBox:
    lower: tuple[float, ...]
    upper: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.lower) != len(self.upper) or not self.lower:
            raise ValueError("matching nonempty box endpoints required")
        if any((not math.isfinite(a) or not math.isfinite(b) or not a < b) for a, b in zip(self.lower, self.upper)):
            raise ValueError("finite strict half-open box intervals required")

    def contains(self, zeta: Sequence[float]) -> bool:
        z = tuple(float(x) for x in zeta)
        if len(z) != len(self.lower):
            raise ValueError("zeta dimension does not match material box")
        return all(a <= x < b for x, a, b in zip(z, self.lower, self.upper))


def old_pool_membership(zeta: Sequence[float], boxes: Sequence[IntrinsicBox]) -> bool:
    """Membership in a fixed old material pool in intrinsic-zeta coordinates.

    The pool is a union of half-open measurable boxes.  Half-open conventions
    name null boundaries only; they do not create material service mass.
    """
    return any(box.contains(zeta) for box in boxes)


def heat_edge_intrinsic_endpoints(
    L: np.ndarray,
    X: np.ndarray,
    k: np.ndarray,
    displacement: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Intrinsic endpoints of the coherent increment edge.

    Translation covariance reads
      V_g(delta_r f)(X,k)=exp(-ik.r)V_g f(X-r,k)-V_g f(X,k),
    so the positive edge atom connects (X,k) and (X-r,k).
    """
    L = np.asarray(L, float)
    X = np.asarray(X, float)
    k = np.asarray(k, float)
    r = np.asarray(displacement, float)
    if X.shape != (3,) or k.shape != (3,) or r.shape != (3,):
        raise ValueError("X,k,r must be three-vectors")
    return intrinsic_zeta(L, X, k), intrinsic_zeta(L, X - r, k)


def affine_heat_edge_endpoint_residual(
    M: np.ndarray,
    L: np.ndarray,
    X: np.ndarray,
    k: np.ndarray,
    displacement: np.ndarray,
) -> float:
    """Common affine/Kelvin transport preserves both edge endpoint labels.

    Under L->ML, X->MX, r->Mr and k->M^{-T}k, both intrinsic endpoint
    coordinates are exactly unchanged.  This certifies ownership of an edge,
    not invariance of the isotropic heat-kernel weight under arbitrary M.
    """
    M = np.asarray(M, float)
    if M.shape != (3, 3) or abs(float(np.linalg.det(M))) < 1e-12:
        raise ValueError("invertible 3x3 affine map required")
    L = np.asarray(L, float)
    X = np.asarray(X, float)
    k = np.asarray(k, float)
    r = np.asarray(displacement, float)
    z0, z1 = heat_edge_intrinsic_endpoints(L, X, k, r)
    z0p, z1p = heat_edge_intrinsic_endpoints(M @ L, M @ X, np.linalg.solve(M.T, k), M @ r)
    return float(max(np.linalg.norm(z0p - z0), np.linalg.norm(z1p - z1)))


def ownership_class(old_here: bool, old_neighbor: bool) -> str:
    if old_here and old_neighbor:
        return "old_old"
    if old_here or old_neighbor:
        return "old_new_interface"
    return "new_new"


def partition_positive_edge_measure(
    edge_weights: Sequence[float],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
) -> dict[str, float]:
    """Exact OO/ON/NN partition of a positive edge measure.

    This acts *after* Moyal/heat disintegration.  It never writes f=f_old+f_new,
    hence there are no quadratic cross terms to estimate or discard.
    """
    w = np.asarray(edge_weights, float)
    a = np.asarray(old_here, bool)
    b = np.asarray(old_neighbor, bool)
    if w.ndim != 1 or a.shape != w.shape or b.shape != w.shape:
        raise ValueError("matching one-dimensional edge ownership data required")
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("finite nonnegative edge weights required")
    oo = a & b
    on = np.logical_xor(a, b)
    nn = (~a) & (~b)
    old_old = float(w[oo].sum())
    interface = float(w[on].sum())
    new_new = float(w[nn].sum())
    total = float(w.sum())
    return {
        "total": total,
        "old_old": old_old,
        "old_new_interface": interface,
        "new_new": new_new,
        "partition_residual": old_old + interface + new_new - total,
    }


def ownership_local_capacity(
    here_amplitudes: Sequence[complex],
    neighbor_amplitudes: Sequence[complex],
    phases: Sequence[complex],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
    positive_prefactors: Sequence[float] | None = None,
) -> dict[str, float]:
    """Ownership-wise endpoint-energy control for coherent increment edges.

    For |phase|=1,
      |phase*A_neighbor-A_here|^2 <= 2(|A_here|^2+|A_neighbor|^2).
    Multiplying by the OO/ON/NN indicators preserves this inequality class by
    class.  This is a local capacity statement only; it does not import the
    signed-good old-pool half-life through a high-strain slab.
    """
    A0 = np.asarray(here_amplitudes, complex)
    A1 = np.asarray(neighbor_amplitudes, complex)
    ph = np.asarray(phases, complex)
    oh = np.asarray(old_here, bool)
    on = np.asarray(old_neighbor, bool)
    if A0.ndim != 1 or A1.shape != A0.shape or ph.shape != A0.shape or oh.shape != A0.shape or on.shape != A0.shape:
        raise ValueError("matching one-dimensional endpoint arrays required")
    if np.any(~np.isfinite(A0)) or np.any(~np.isfinite(A1)) or np.any(~np.isfinite(ph)):
        raise ValueError("finite endpoint amplitudes/phases required")
    if np.any(np.abs(np.abs(ph) - 1.0) > 2e-12):
        raise ValueError("unit-modulus translation phases required")
    if positive_prefactors is None:
        p = np.ones(len(A0), float)
    else:
        p = np.asarray(positive_prefactors, float)
        if p.shape != A0.shape or np.any(~np.isfinite(p)) or np.any(p < 0):
            raise ValueError("finite nonnegative prefactors required")
    edge = p * np.abs(ph * A1 - A0) ** 2
    cap = 2.0 * p * (np.abs(A0) ** 2 + np.abs(A1) ** 2)
    masks = {
        "old_old": oh & on,
        "old_new_interface": np.logical_xor(oh, on),
        "new_new": (~oh) & (~on),
    }
    out: dict[str, float] = {}
    for name, mask in masks.items():
        service = float(edge[mask].sum())
        capacity = float(cap[mask].sum())
        out[f"{name}_service"] = service
        out[f"{name}_endpoint_capacity"] = capacity
        out[f"{name}_margin"] = capacity - service
    out["total_service"] = float(edge.sum())
    out["total_endpoint_capacity"] = float(cap.sum())
    return out


def ownership_change_requires_boundary_contact(
    start_zeta: Sequence[float],
    end_zeta: Sequence[float],
    box: IntrinsicBox,
) -> bool:
    """Finite-endpoint diagnostic for the continuum boundary-crossing fact.

    For a continuous intrinsic endpoint path, if membership in a half-open box
    changes, the path must meet the topological boundary.  This helper returns
    whether endpoint membership differs; the theorem certificate states the
    continuum intermediate-value conclusion.
    """
    return box.contains(start_zeta) != box.contains(end_zeta)


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_POSITIVE_HEAT_EDGE_MATERIAL_OWNERSHIP_PARTITION__AFFINE_ENDPOINT_LABELS_INVARIANT__CAPACITY_ROUTING_REMAINS",
        "edge": "a coherent heat-increment atom connects intrinsic endpoints zeta(X,k) and zeta(X-r,k)",
        "partition": "for any fixed old material set O, 1_O(z0)1_O(z1)+1_O(z0)1_Oc(z1)+1_Oc(z0)1_O(z1)+1_Oc(z0)1_Oc(z1)=1 pointwise on the positive edge density",
        "no_cross_term": "OO/ON/NN classification is performed after positive Moyal disintegration, not by decomposing the velocity field; therefore the three edge measures add exactly",
        "affine_invariance": "under L->ML, X->MX, r->Mr, k->M^-T k, both intrinsic edge endpoints are unchanged; common affine motion cannot change material ownership",
        "interface": "old-new is an actual material-interface edge mark; a time-dependent ownership change under continuous nonaffine evolution can occur only through the material-pool boundary",
        "local_capacity": "classwise |e^{-ik.r}A_neighbor-A_here|^2 <= 2(|A_here|^2+|A_neighbor|^2)",
        "boundary": "half-open dyadic conventions only assign null boundaries; they do not manufacture edge service or relink mass",
        "scope": "this closes exact ownership classification of the high-strain heat-edge law; it does not prove old-old lifetime capacity on a high-strain slab, nor that new-new service is already a selected transfer parent",
    }


@dataclass(frozen=True)
class OwnershipStress:
    samples: int
    worst_affine_endpoint_residual: float
    worst_partition_residual: float
    minimum_class_capacity_margin: float
    orientation_failures: int
    membership_invariance_failures: int


def stress(samples: int = 50_000, seed: int = 20260809) -> OwnershipStress:
    rng = np.random.default_rng(seed)
    wa = wp = 0.0
    mm = float("inf")
    orient_fail = membership_fail = 0
    for _ in range(samples):
        # Exact material endpoint covariance.
        Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        L = Q @ np.diag(np.exp(rng.uniform(-0.8, 0.8, size=3)))
        M = np.eye(3) + 0.15 * rng.normal(size=(3, 3))
        if abs(float(np.linalg.det(M))) < 0.2:
            M += np.eye(3)
        X = rng.normal(size=3)
        k = rng.normal(size=3)
        r = rng.normal(size=3)
        ar = affine_heat_edge_endpoint_residual(M, L, X, k, r)
        wa = max(wa, ar)
        if ar > 3e-11 * max(1.0, np.linalg.norm(X), np.linalg.norm(k), np.linalg.norm(r)):
            raise AssertionError("common affine transport changed a heat-edge material endpoint")

        z0, z1 = heat_edge_intrinsic_endpoints(L, X, k, r)
        center = rng.normal(size=6)
        width = np.exp(rng.uniform(-1.0, 0.8, size=6))
        box = IntrinsicBox(tuple(center - width), tuple(center + width))
        o0 = box.contains(z0)
        o1 = box.contains(z1)
        z0p, z1p = heat_edge_intrinsic_endpoints(M @ L, M @ X, np.linalg.solve(M.T, k), M @ r)
        if box.contains(z0p) != o0 or box.contains(z1p) != o1:
            membership_fail += 1
            raise AssertionError("affine gauge changed old-pool endpoint ownership")

        n = int(rng.integers(2, 100))
        weights = rng.lognormal(mean=-1.0, sigma=1.5, size=n)
        a = rng.random(n) < 0.5
        b = rng.random(n) < 0.5
        out = partition_positive_edge_measure(weights, a, b)
        wp = max(wp, abs(float(out["partition_residual"])))
        if abs(float(out["partition_residual"])) > 3e-12 * max(1.0, float(out["total"])):
            raise AssertionError("positive heat-edge ownership partition lost mass")
        swapped = partition_positive_edge_measure(weights, b, a)
        if any(abs(float(out[key]) - float(swapped[key])) > 2e-12 * max(1.0, float(out["total"])) for key in ("old_old", "old_new_interface", "new_new")):
            orient_fail += 1
            raise AssertionError("unoriented heat-edge ownership depended on endpoint order")

        A0 = rng.normal(size=n) + 1j * rng.normal(size=n)
        A1 = rng.normal(size=n) + 1j * rng.normal(size=n)
        theta = rng.uniform(-math.pi, math.pi, size=n)
        ph = np.exp(1j * theta)
        pref = rng.lognormal(mean=-1.0, sigma=1.0, size=n)
        cap = ownership_local_capacity(A0, A1, ph, a, b, pref)
        for name in ("old_old", "old_new_interface", "new_new"):
            margin = float(cap[f"{name}_margin"])
            mm = min(mm, margin)
            if margin < -4e-12 * max(1.0, float(cap[f"{name}_endpoint_capacity"])):
                raise AssertionError("ownership-specific edge service exceeded endpoint Moyal capacity")

    return OwnershipStress(samples, wa, wp, mm, orient_fail, membership_fail)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-heat-edge-material-ownership"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "heat_edge_material_ownership.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Material ownership of positive heat-increment edges\n\nStatus: **{cert['status']}**.\n\nThe high-strain heat theorem already supplies a positive coherent edge measure.  Material ownership should be assigned to that **measure**, not to a decomposition of the velocity field.\n\nFor one coherent heat edge, translation covariance identifies its endpoints as `(X,k)` and `(X-r,k)`.  In the intrinsic material coordinate\n\n`zeta=(L^-1 X/2,L^T k)`,\n\nlet `O` be the fixed transported old material pool and put `chi_i=1_O(zeta_i)`.  Pointwise on every nonnegative edge weight `s`,\n\n`s_OO=chi_0 chi_1 s`,\n\n`s_ON=[chi_0(1-chi_1)+(1-chi_0)chi_1]s`,\n\n`s_NN=(1-chi_0)(1-chi_1)s`,\n\nand exactly\n\n`s_OO+s_ON+s_NN=s`.\n\nThere is no quadratic cross term because no identity `V=V_old+V_new` is used.  The partition is performed **after** the Moyal/heat edge density has become positive.  It is also unoriented: swapping the two physical endpoints leaves OO, interface and NN weights unchanged.\n\nUnder common affine/Kelvin transport\n\n`L->ML,  X->MX,  r->Mr,  k->M^-T k`,\n\nboth endpoint labels are individually invariant:\n\n`zeta(ML,MX,M^-T k)=zeta(L,X,k)`,\n\n`zeta(ML,M(X-r),M^-T k)=zeta(L,X-r,k)`.\n\nThus common affine motion cannot convert an old--old edge into an interface or new--new edge.  A change of endpoint membership along a continuous nonaffine material evolution must meet the boundary of the old material set.  Half-open dyadic conventions merely assign that null boundary and cannot manufacture service mass.\n\nThe exact coherent increment identity also yields, class by class,\n\n`|e^(-ik.r)A_1-A_0|^2 <= 2(|A_0|^2+|A_1|^2)`.\n\nHence each ownership class is locally supported by the Moyal energy of its own two endpoints.  This is the correct starting point for old-pool capacity; the theorem does **not** import the signed-good low-strain half-life through a high-strain slab.\n\nStress: `{out.samples}` material-edge states\n- worst affine two-endpoint residual: `{out.worst_affine_endpoint_residual:.3e}`\n- worst positive-measure partition residual: `{out.worst_partition_residual:.3e}`\n- minimum ownership-specific endpoint-capacity margin: `{out.minimum_class_capacity_margin:.3e}`\n- orientation failures: `{out.orientation_failures}`\n- affine membership-invariance failures: `{out.membership_invariance_failures}`\n\nThis closes **material ownership classification** of the high-strain heat-edge seed.  The remaining high-strain problem is quantitative routing: bound repeated OO heat service using the physically valid history of the transported old pool; route ON service as actual material-interface/relink provenance; and show what NN service creates without postulating a packet mass floor.  Universal slab renewal remains open, and no global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
