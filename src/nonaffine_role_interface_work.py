from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


def adjoint_split(L: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return K=(L-L*)/2 and S=(L+L*)/2."""
    L = np.asarray(L, complex)
    if L.ndim != 2 or L.shape[0] != L.shape[1]:
        raise ValueError("square operator required")
    K = 0.5 * (L - L.conj().T)
    S = 0.5 * (L + L.conj().T)
    return K, S


def commutator(L: np.ndarray, Q: np.ndarray) -> np.ndarray:
    L = np.asarray(L, complex)
    Q = np.asarray(Q, complex)
    if L.shape != Q.shape or L.ndim != 2 or L.shape[0] != L.shape[1]:
        raise ValueError("matching square operators required")
    return L @ Q - Q @ L


def interface_work(Q: np.ndarray, L: np.ndarray, u: np.ndarray) -> float:
    """2 Re <Q u,[L,Q]u>, the moving-role interface energy work."""
    Q = np.asarray(Q, complex)
    L = np.asarray(L, complex)
    u = np.asarray(u, complex)
    if Q.shape != L.shape or Q.shape != (len(u), len(u)):
        raise ValueError("operator/vector dimension mismatch")
    w = Q @ u
    return 2.0 * float(np.real(np.vdot(w, commutator(L, Q) @ u)))


def projector_residual(Q: np.ndarray) -> float:
    Q = np.asarray(Q, complex)
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
        raise ValueError("square projector candidate required")
    return max(
        float(np.linalg.norm(Q - Q.conj().T)),
        float(np.linalg.norm(Q @ Q - Q)),
    )


def binary_interface_decomposition(Q: np.ndarray, L: np.ndarray, u: np.ndarray) -> dict[str, float]:
    """Exact two-role decomposition for Q and I-Q.

    If K*=-K and S*=S are the adjoint parts of L, z=Qu and y=(I-Q)u,

      I_Q^K = -2 Re <z,K y>,       I_(I-Q)^K = -I_Q^K,
      I_Q^S = -2 Re <z,S y>,       I_(I-Q)^S =  I_Q^S.

    The skew piece is conservative role-to-role transfer.  The symmetric piece
    is the off-diagonal physical strain/deformation work.
    """
    Q = np.asarray(Q, complex)
    L = np.asarray(L, complex)
    u = np.asarray(u, complex)
    n = len(u)
    if Q.shape != (n, n) or L.shape != (n, n):
        raise ValueError("matching operator/vector data required")
    pres = projector_residual(Q)
    if pres > 2e-10:
        raise ValueError("orthogonal projector required")
    P = np.eye(n, dtype=complex) - Q
    K, S = adjoint_split(L)
    z = Q @ u
    y = P @ u
    iqk = interface_work(Q, K, u)
    ipk = interface_work(P, K, u)
    iqs = interface_work(Q, S, u)
    ips = interface_work(P, S, u)
    pair_k = -2.0 * float(np.real(np.vdot(z, K @ y)))
    pair_s = -2.0 * float(np.real(np.vdot(z, S @ y)))
    return {
        "Q_skew_work": iqk,
        "complement_skew_work": ipk,
        "skew_pair_formula": pair_k,
        "Q_symmetric_work": iqs,
        "complement_symmetric_work": ips,
        "symmetric_pair_formula": pair_s,
        "skew_conservation_residual": iqk + ipk,
        "Q_skew_formula_residual": iqk - pair_k,
        "Q_symmetric_formula_residual": iqs - pair_s,
        "symmetric_pair_equality_residual": iqs - ips,
    }


def partition_interface_balance(
    projectors: Sequence[np.ndarray],
    L: np.ndarray,
    u: np.ndarray,
) -> dict[str, object]:
    """Exact balance for a complete orthogonal role partition.

    For K*=-K, pair transfers
      T_ab=-2 Re <w_a,K w_b>
    are antisymmetric and row sums are the K-interface works.

    For S*=S, pair cross works
      D_ab=-2 Re <w_a,S w_b>
    are symmetric and row sums are the S-interface works.  Their total obeys

      sum_a I_a^S
       = 2 sum_a Re<w_a,S w_a> - 2 Re<u,S u>.

    Thus the symmetric interface is exactly the off-diagonal part of the same
    physical strain work, not a new source.
    """
    L = np.asarray(L, complex)
    u = np.asarray(u, complex)
    Ps = [np.asarray(P, complex) for P in projectors]
    n = len(u)
    if not Ps or any(P.shape != (n, n) for P in Ps) or L.shape != (n, n):
        raise ValueError("complete matching role partition required")
    total = sum(Ps, np.zeros((n, n), complex))
    if np.linalg.norm(total - np.eye(n)) > 3e-10:
        raise ValueError("projectors must resolve identity")
    for i, P in enumerate(Ps):
        if projector_residual(P) > 2e-10:
            raise ValueError("orthogonal projectors required")
        for Q in Ps[i + 1 :]:
            if np.linalg.norm(P @ Q) > 2e-10:
                raise ValueError("projector ranges must be orthogonal")

    K, S = adjoint_split(L)
    ws = [P @ u for P in Ps]
    m = len(Ps)
    TK = np.zeros((m, m), float)
    DS = np.zeros((m, m), float)
    for a in range(m):
        for b in range(m):
            if a == b:
                continue
            TK[a, b] = -2.0 * float(np.real(np.vdot(ws[a], K @ ws[b])))
            DS[a, b] = -2.0 * float(np.real(np.vdot(ws[a], S @ ws[b])))

    ik = np.array([interface_work(P, K, u) for P in Ps])
    isym = np.array([interface_work(P, S, u) for P in Ps])
    diag_strain = 2.0 * sum(float(np.real(np.vdot(w, S @ w))) for w in ws)
    full_strain = 2.0 * float(np.real(np.vdot(u, S @ u)))
    return {
        "skew_pair_matrix": TK,
        "symmetric_pair_matrix": DS,
        "skew_interface_work": ik,
        "symmetric_interface_work": isym,
        "skew_antisymmetry_residual": float(np.linalg.norm(TK + TK.T)),
        "symmetric_symmetry_residual": float(np.linalg.norm(DS - DS.T)),
        "skew_row_sum_residual": float(np.linalg.norm(TK.sum(axis=1) - ik)),
        "symmetric_row_sum_residual": float(np.linalg.norm(DS.sum(axis=1) - isym)),
        "total_skew_interface": float(ik.sum()),
        "symmetric_global_balance_residual": float(isym.sum() - (diag_strain - full_strain)),
        "diagonal_strain_work": diag_strain,
        "full_strain_work": full_strain,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_NONAFFINE_ROLE_INTERFACE_SPLIT__SKEW_TRANSFER_CONSERVATIVE__SYMMETRIC_WORK_IS_STRAIN_PROVENANCE",
        "operator_split": "L_nonaff=K+S with K*=-K and S*=S; for incompressible resolved transport K is advection+rotation and S is the symmetric deformation operator",
        "skew": "T_ab=-2 Re<w_a,K w_b> is antisymmetric; its row sums are role-interface works and total skew interface work is exactly zero",
        "symmetric": "D_ab=-2 Re<w_a,S w_b> is symmetric; its row sums are the symmetric interface works",
        "global_strain_identity": "sum interface_S = 2 sum_a Re<w_a,S w_a> - 2 Re<u,S u>; interface_S is exactly off-diagonal strain work",
        "single_charge": "skew interface is conservative role relinking/redistribution, while symmetric interface delegates to the already existing strain/coherent-deformation cause; neither is a new source currency and neither is representation Xi",
        "continuum_status": "this is the exact hard-projector/event-role interface lemma; a non-idempotent smooth PDE envelope must instead use the quadratic Q^2 carrier-energy theorem before quantitative K_coh/strain/D_V or material-relink routing",
    }


@dataclass(frozen=True)
class InterfaceStress:
    samples: int
    worst_binary_residual: float
    worst_skew_antisymmetry_residual: float
    worst_symmetric_symmetry_residual: float
    worst_skew_row_sum_residual: float
    worst_symmetric_row_sum_residual: float
    worst_total_skew_interface: float
    worst_symmetric_global_balance_residual: float


def _random_partition(rng: np.random.Generator, n: int) -> list[np.ndarray]:
    Z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    U, _ = np.linalg.qr(Z)
    cuts = (1, max(2, n // 2))
    groups = [range(0, cuts[0]), range(cuts[0], cuts[1]), range(cuts[1], n)]
    out: list[np.ndarray] = []
    for g in groups:
        cols = U[:, list(g)]
        out.append(cols @ cols.conj().T)
    return out


def stress(samples: int = 50_000, seed: int = 20260809) -> InterfaceStress:
    rng = np.random.default_rng(seed)
    wb = wa = ws = wrk = wrs = wt = wg = 0.0
    n = 6
    for _ in range(samples):
        Ps = _random_partition(rng, n)
        Q = Ps[0] + Ps[1]
        L = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        u = rng.normal(size=n) + 1j * rng.normal(size=n)
        b = binary_interface_decomposition(Q, L, u)
        vals = [
            abs(float(b["skew_conservation_residual"])),
            abs(float(b["Q_skew_formula_residual"])),
            abs(float(b["Q_symmetric_formula_residual"])),
            abs(float(b["symmetric_pair_equality_residual"])),
        ]
        bscale = max(1.0, np.linalg.norm(L) * np.linalg.norm(u) ** 2)
        wb = max(wb, max(vals) / bscale)
        if max(vals) > 8e-11 * bscale:
            raise AssertionError("binary role-interface decomposition failed")

        p = partition_interface_balance(Ps, L, u)
        scale = max(1.0, np.linalg.norm(L) * np.linalg.norm(u) ** 2)
        wa = max(wa, float(p["skew_antisymmetry_residual"]) / scale)
        ws = max(ws, float(p["symmetric_symmetry_residual"]) / scale)
        wrk = max(wrk, float(p["skew_row_sum_residual"]) / scale)
        wrs = max(wrs, float(p["symmetric_row_sum_residual"]) / scale)
        wt = max(wt, abs(float(p["total_skew_interface"])) / scale)
        wg = max(wg, abs(float(p["symmetric_global_balance_residual"])) / scale)
        if max(wa, ws, wrk, wrs, wt, wg) > 2e-10:
            raise AssertionError("partition role-interface balance failed")
    return InterfaceStress(samples, wb, wa, ws, wrk, wrs, wt, wg)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-nonaffine-role-interface-work"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "nonaffine_role_interface_work.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Nonaffine hard-projector role-interface work\n\nStatus: **{cert['status']}**.\n\nFor the nonaffine resolved linearized operator write `L=K+S`, `K*=-K`, `S*=S`.  For any complete orthogonal **event-role** partition `u=sum_a w_a`, the Heisenberg interface work splits exactly.\n\nThe skew/advection+rotation piece has pair flux\n\n`T_ab=-2 Re <w_a,K w_b>`,\n\nwith `T_ab=-T_ba`.  Its row sums are the individual role-interface works, so the total skew interface work is **exactly zero**.  It is conservative physical role redistribution/relinking, never energy generation.\n\nThe symmetric piece has\n\n`D_ab=-2 Re <w_a,S w_b>=D_ba`.\n\nIts row sums are the symmetric interface works and\n\n`sum_a I_a^S = 2 sum_a Re<w_a,S w_a> - 2 Re<u,S u>`.\n\nThus the symmetric hard-role interface is precisely the off-diagonal part of the **same physical strain/deformation work** already present in the resolved transporter.  It delegates to coherent deformation / strain / critical `D_V`; it is not a new currency and not representation `Xi`.\n\nStress: `{out.samples}` random complex Hilbert-space partitions/operators\n- worst binary identity residual: `{out.worst_binary_residual:.3e}`\n- worst skew antisymmetry residual: `{out.worst_skew_antisymmetry_residual:.3e}`\n- worst symmetric symmetry residual: `{out.worst_symmetric_symmetry_residual:.3e}`\n- worst skew row-sum residual: `{out.worst_skew_row_sum_residual:.3e}`\n- worst symmetric row-sum residual: `{out.worst_symmetric_row_sum_residual:.3e}`\n- worst total skew-interface residual: `{out.worst_total_skew_interface:.3e}`\n- worst symmetric global strain-balance residual: `{out.worst_symmetric_global_balance_residual:.3e}`\n\nThis theorem is deliberately a hard-projector lemma.  It must not be applied directly to the non-idempotent smooth PDE envelope.  The companion smooth quadratic-carrier theorem reads the propagated energy at `Q^2`, recombines the outer commutator with diagonal resolved-role work, and only then obtains conservative relink plus existing strain.  Quantitative first stopping remains governed by actual energy work, never raw coefficient impulse.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
