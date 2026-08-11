from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.coherent_localization_operators import partition_operators, random_parseval_frame, random_partition
from src.physical_energy_causal_bridge import PHYSICAL_HH_WORK_FRACTION


def bilinear_tensor_apply(B: np.ndarray, f: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Finite-dimensional model of one exact quadratic Navier--Stokes source role.

    B[k,i,j] is the bilinear source tensor.  No symmetry is assumed: the actual
    designated parent-role source may already include the j<->k symmetrization.
    """
    B = np.asarray(B, complex)
    f = np.asarray(f, complex)
    g = np.asarray(g, complex)
    if B.ndim != 3 or B.shape[0] != B.shape[1] or B.shape[1] != B.shape[2]:
        raise ValueError("cubic bilinear tensor required")
    n = B.shape[0]
    if f.shape != (n,) or g.shape != (n,):
        raise ValueError("parent vectors have wrong dimension")
    return np.einsum("kij,i,j->k", B, f, g, optimize=True)


def _validate_localization_partition(ops: Sequence[np.ndarray], n: int, tol: float = 5e-11) -> tuple[np.ndarray, ...]:
    out = tuple(np.asarray(A, complex) for A in ops)
    if not out:
        raise ValueError("nonempty localization partition required")
    for A in out:
        if A.shape != (n, n):
            raise ValueError("localization operator dimension mismatch")
        if np.linalg.norm(A - A.conj().T) > tol:
            raise ValueError("coherent localization operators must be self-adjoint")
        ev = np.linalg.eigvalsh(0.5 * (A + A.conj().T))
        if float(ev.min(initial=0.0)) < -tol:
            raise ValueError("coherent localization operator is not positive")
    S = sum(out, np.zeros((n, n), complex))
    if np.linalg.norm(S - np.eye(n)) > tol:
        raise ValueError("localization operators must resolve the identity")
    return out


def bilinear_partition_residual(
    B: np.ndarray,
    f: np.ndarray,
    g: np.ndarray,
    parent1_ops: Sequence[np.ndarray],
    parent2_ops: Sequence[np.ndarray],
) -> np.ndarray:
    """Residual in N(f,g)=sum_CD N(A_C f,A_D g)."""
    B = np.asarray(B, complex)
    n = B.shape[0]
    p1 = _validate_localization_partition(parent1_ops, n)
    p2 = _validate_localization_partition(parent2_ops, n)
    exact = bilinear_tensor_apply(B, f, g)
    expanded = np.zeros(n, complex)
    for A in p1:
        af = A @ np.asarray(f, complex)
        for D in p2:
            expanded += bilinear_tensor_apply(B, af, D @ np.asarray(g, complex))
    return exact - expanded


def coherent_binary_work_atoms(
    B: np.ndarray,
    parent1: np.ndarray,
    parent2: np.ndarray,
    child: np.ndarray,
    parent1_ops: Sequence[np.ndarray],
    parent2_ops: Sequence[np.ndarray],
    child_ops: Sequence[np.ndarray],
) -> np.ndarray:
    """Exact signed child-energy work atoms indexed by (parent C,parent D,child E).

    W_CDE = 2 Re <A_E child, N(A_C parent1,A_D parent2)>.

    Since every localization family resolves the identity and A_E is self-adjoint,
    sum_CDE W_CDE = 2 Re <child,N(parent1,parent2)> exactly.  Each atom therefore
    has two physical parent labels and one child label before any positive-part
    selection.  The coherent pieces are analysis/synthesis devices; no assertion
    that an individual A_C f has compact Fourier support is made here.
    """
    B = np.asarray(B, complex)
    n = B.shape[0]
    p1 = _validate_localization_partition(parent1_ops, n)
    p2 = _validate_localization_partition(parent2_ops, n)
    pc = _validate_localization_partition(child_ops, n)
    f = np.asarray(parent1, complex)
    g = np.asarray(parent2, complex)
    h = np.asarray(child, complex)
    if f.shape != (n,) or g.shape != (n,) or h.shape != (n,):
        raise ValueError("role vector dimension mismatch")
    out = np.zeros((len(p1), len(p2), len(pc)), float)
    child_pieces = [E @ h for E in pc]
    for i, A in enumerate(p1):
        af = A @ f
        for j, D in enumerate(p2):
            source = bilinear_tensor_apply(B, af, D @ g)
            for k, eh in enumerate(child_pieces):
                out[i, j, k] = 2.0 * float(np.real(np.vdot(eh, source)))
    return out


@dataclass(frozen=True)
class BinaryWorkLedger:
    signed_work: float
    positive_transfer_mass: float
    negative_backscatter_mass: float
    aggregate_positive_work: float
    positive_dominance_margin: float
    cancellation_residual: float


def binary_work_ledger(atoms: np.ndarray) -> BinaryWorkLedger:
    """Positive/negative Hahn split of the exact coherent triple work atoms."""
    w = np.asarray(atoms, float)
    if w.ndim != 3 or not np.all(np.isfinite(w)):
        raise ValueError("finite three-index work atom array required")
    signed = float(w.sum())
    pos = float(np.maximum(w, 0.0).sum())
    neg = float(np.maximum(-w, 0.0).sum())
    aggregate_positive = max(0.0, signed)
    residual = (pos - neg) - signed
    margin = pos - aggregate_positive
    if abs(residual) > 3e-13 * max(1.0, pos + neg, abs(signed)):
        raise AssertionError("positive/negative coherent-work split lost signed work")
    if margin < -3e-13 * max(1.0, pos):
        raise AssertionError("atomic positive transfer failed to dominate aggregate positive work")
    return BinaryWorkLedger(signed, pos, neg, aggregate_positive, margin, residual)


@dataclass(frozen=True)
class PositiveBinaryEvent:
    parent1: int
    parent2: int
    child: int
    mass: float
    probability: float


def normalized_positive_binary_events(atoms: np.ndarray) -> tuple[PositiveBinaryEvent, ...]:
    """Normalize the positive Hahn atoms of this coherent signed-work representation.

    This is an exact representation-level probability law.  After canonical
    continuum edge registration it is *not* automatically the master causal law:
    that identification requires a positive mass-preserving kernel from the
    already-fixed physical edge law ``dW+``.
    """
    w = np.asarray(atoms, float)
    pos = np.maximum(w, 0.0)
    total = float(pos.sum())
    if total <= 0:
        raise ValueError("positive physical work mass required")
    events: list[PositiveBinaryEvent] = []
    for i, j, k in np.argwhere(pos > 0):
        m = float(pos[i, j, k])
        events.append(PositiveBinaryEvent(int(i), int(j), int(k), m, m / total))
    if not math.isclose(sum(e.probability for e in events), 1.0, rel_tol=2e-14, abs_tol=2e-14):
        raise AssertionError("binary causal law failed to normalize")
    return tuple(events)


@dataclass(frozen=True)
class XiExcision:
    total_positive_mass: float
    retained_positive_mass: float
    xi_positive_mass: float
    retained_fraction: float
    retained_events: int
    omitted_events: int


def excise_positive_xi(atoms: np.ndarray, omit_mask: np.ndarray) -> XiExcision:
    """Excise selected physical cross-cell positive work exactly once into Xi."""
    w = np.asarray(atoms, float)
    mask = np.asarray(omit_mask, bool)
    if w.shape != mask.shape or w.ndim != 3:
        raise ValueError("Xi mask must match work atoms")
    pos = np.maximum(w, 0.0)
    total = float(pos.sum())
    xi = float(pos[mask].sum())
    retained = total - xi
    frac = retained / total if total > 0 else 0.0
    if abs((retained + xi) - total) > 2e-13 * max(1.0, total):
        raise AssertionError("Xi excision duplicated or lost positive work")
    return XiExcision(
        total,
        retained,
        xi,
        frac,
        int(np.count_nonzero((~mask) & (pos > 0))),
        int(np.count_nonzero(mask & (pos > 0))),
    )


def retained_generation_lower(
    terminal_child_energy: float,
    xi_fraction_of_positive_hh: float,
    clean_hh_fraction: float = float(PHYSICAL_HH_WORK_FRACTION),
) -> float:
    """Clean retained physical generation after a relative positive-work Xi excision.

    The physical-energy causal gate supplies W_HH^+ >= clean_hh_fraction E_1.
    If the selected cross-cell moat removes at most rho of positive HH work, the
    retained binary physical-transfer mass is at least (1-rho)*clean*E_1.
    """
    if terminal_child_energy < 0 or not math.isfinite(terminal_child_energy):
        raise ValueError("finite nonnegative child energy required")
    if not (0.0 <= xi_fraction_of_positive_hh <= 1.0):
        raise ValueError("Xi fraction must lie in [0,1]")
    if not (0.0 <= clean_hh_fraction <= 1.0):
        raise ValueError("clean HH fraction must lie in [0,1]")
    return (1.0 - xi_fraction_of_positive_hh) * clean_hh_fraction * terminal_child_energy


def exact_work_reconstruction_residual(
    B: np.ndarray,
    parent1: np.ndarray,
    parent2: np.ndarray,
    child: np.ndarray,
    atoms: np.ndarray,
) -> float:
    exact = 2.0 * float(np.real(np.vdot(np.asarray(child, complex), bilinear_tensor_apply(B, parent1, parent2))))
    return float(np.asarray(atoms, float).sum()) - exact


@dataclass(frozen=True)
class RecursiveWitnessStress:
    samples: int
    worst_bilinear_partition_residual: float
    worst_work_reconstruction_residual: float
    minimum_positive_dominance_margin: float
    worst_probability_residual: float
    worst_xi_partition_residual: float
    minimum_retained_generation_margin: float


def stress(samples: int = 20_000, seed: int = 20260808) -> RecursiveWitnessStress:
    rng = np.random.default_rng(seed)
    wb = ww = wp = wx = 0.0
    mp = mg = float("inf")
    for _ in range(samples):
        n = int(rng.integers(2, 7))
        m1 = int(rng.integers(n, 3 * n + 3))
        m2 = int(rng.integers(n, 3 * n + 3))
        m3 = int(rng.integers(n, 3 * n + 3))
        F1 = random_parseval_frame(rng, n, m1)
        F2 = random_parseval_frame(rng, n, m2)
        F3 = random_parseval_frame(rng, n, m3)
        cells1 = random_partition(rng, m1, int(rng.integers(1, min(m1, 5) + 1)))
        cells2 = random_partition(rng, m2, int(rng.integers(1, min(m2, 5) + 1)))
        cells3 = random_partition(rng, m3, int(rng.integers(1, min(m3, 5) + 1)))
        A1 = partition_operators(F1, cells1)
        A2 = partition_operators(F2, cells2)
        A3 = partition_operators(F3, cells3)
        B = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
        f = rng.normal(size=n) + 1j * rng.normal(size=n)
        g = rng.normal(size=n) + 1j * rng.normal(size=n)
        h = rng.normal(size=n) + 1j * rng.normal(size=n)

        br = bilinear_partition_residual(B, f, g, A1, A2)
        bscale = max(1.0, np.linalg.norm(bilinear_tensor_apply(B, f, g)))
        bres = float(np.linalg.norm(br)) / bscale
        wb = max(wb, bres)
        if bres > 3e-10:
            raise AssertionError("coherent parent partition did not reconstruct quadratic source")

        atoms = coherent_binary_work_atoms(B, f, g, h, A1, A2, A3)
        wr = abs(exact_work_reconstruction_residual(B, f, g, h, atoms)) / max(
            1.0, abs(float(atoms.sum())), np.linalg.norm(h) * np.linalg.norm(bilinear_tensor_apply(B, f, g))
        )
        ww = max(ww, wr)
        if wr > 4e-10:
            raise AssertionError("coherent binary work atoms failed exact reconstruction")
        ledger = binary_work_ledger(atoms)
        mp = min(mp, ledger.positive_dominance_margin)

        if ledger.positive_transfer_mass > 1e-14:
            events = normalized_positive_binary_events(atoms)
            pres = abs(sum(e.probability for e in events) - 1.0)
            wp = max(wp, pres)
            if pres > 3e-13:
                raise AssertionError("positive binary event law did not normalize")

        mask = rng.random(size=atoms.shape) < float(rng.uniform(0.0, 0.5))
        ex = excise_positive_xi(atoms, mask)
        xres = abs(ex.retained_positive_mass + ex.xi_positive_mass - ex.total_positive_mass)
        wx = max(wx, xres)
        if xres > 3e-13 * max(1.0, ex.total_positive_mass):
            raise AssertionError("positive Xi split failed")

        E1 = float(rng.lognormal(mean=0.0, sigma=1.0))
        rho = float(rng.uniform(0.0, 0.4))
        lower = retained_generation_lower(E1, rho)
        exact_clean = (1.0 - rho) * float(PHYSICAL_HH_WORK_FRACTION) * E1
        mg = min(mg, lower - exact_clean)
        if abs(lower - exact_clean) > 2e-14 * max(1.0, exact_clean):
            raise AssertionError("retained physical generation lower bound changed")

    return RecursiveWitnessStress(samples, wb, ww, mp, wp, wx, mg)


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_BINARY_COHERENT_WORK_ATOMIZATION__EVENT_OUTER_ROLES_SUPPLIED",
        "parent_source": "N(w1,w2)=sum_CD N(A_C w1,A_D w2) exactly",
        "work_atoms": "W_CDE=2 Re <A_E w3,N(A_C w1,A_D w2)>",
        "work_identity": "sum_CDE W_CDE=2 Re <w3,N(w1,w2)>",
        "positive_measure": "P=sum[W_CDE]_+, N=sum[-W_CDE]_+, P-N=W and P>=[W]_+; this is the Hahn law of the coherent signed representation, not automatically the canonical causal law",
        "causal_event": "conditional only: a positive coherent atom becomes master-facing causal mass only after a positive mass-preserving pushforward/disintegration from canonical edge dW+ is proved",
        "xi_rule": "omitted positive cross-cell atoms are excised once; retained+Xi=total positive work",
        "energy_bridge": "W_HH^+>=8E1/15 implies atomic positive mass at least 8E1/15 before Xi; after relative Xi rho retain >=(1-rho)8E1/15",
        "label_rule": "use the canonical intrinsic-zeta dyadic addresses; no new packet label is introduced",
        "important_scope": "A_C pieces need not be Fourier compact; scale/helical representative errors remain in the existing symbol/covariance Xi ledgers",
        "continuum_status": "deterministic hard Fourier/helicity cells may inherit canonical dW+ by pushforward; a general coherent POVM remains signed diagnostic until a positive dW+ kernel theorem is supplied",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=20_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-recursive-coherent-witness-extraction"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "recursive_coherent_witness_extraction.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Recursive coherent witness extraction\n\nStatus: **{cert['status']}**.\n\nOnce the exact selected divergence-free parent/child roles exist, no packet synthesis is needed to create binary causal events.  For any coherent localization partitions with `sum_C A_C=I`, bilinearity gives\n\n`N(w1,w2)=sum_(C,D) N(A_C w1,A_D w2)`\n\nand the actual child-energy work decomposes exactly as\n\n`W_CDE=2 Re <A_E w3,N(A_C w1,A_D w2)>`,\n`sum_CDE W_CDE=2 Re <w3,N(w1,w2)>`.\n\nTaking positive and negative parts is a Hahn decomposition of the **coherent signed-work representation**: `P-N=W` and `P>=[W]_+`.  The signed identity remains exact.  However, after the canonical continuum edge law has been fixed, these positive coherent Hahn atoms are not automatically a second master causal law.  A general POVM localization must first be supplied with a positive mass-preserving kernel from canonical edge `dW+`; absent that theorem, the coherent positive atoms are representation diagnostics rather than newly minted cause.\n\nCombining with the physical-energy causal gate, the generated low-strain branch has atomic positive mass at least `8E1/15` before selected cross-cell excision.  If the physical defect moat removes a relative positive-work fraction `rho`, retained binary generation is at least `(1-rho)8E1/15`.  Omitted positive atoms enter `Xi` once; negative atoms remain physical backscatter/cancellation.\n\nThe labels are the existing intrinsic-zeta material coherent addresses.  No assertion is made that an individual `A_C w` has compact Fourier support; scale/helical representative information is still supplied by the outer selected role and the already summable symbol/covariance representation ledgers.\n\nStress: `{out.samples}` finite Parseval/POVM bilinear-work states\n- worst quadratic-source partition residual: `{out.worst_bilinear_partition_residual:.3e}`\n- worst work reconstruction residual: `{out.worst_work_reconstruction_residual:.3e}`\n- minimum atomic-positive dominance margin: `{out.minimum_positive_dominance_margin:.3e}`\n- worst binary probability residual: `{out.worst_probability_residual:.3e}`\n- worst Xi positive-mass partition residual: `{out.worst_xi_partition_residual:.3e}`\n- minimum retained-generation formula margin: `{out.minimum_retained_generation_margin:.3e}`\n\nThe deterministic hard Fourier/helicity event-role path is now the preferred master-facing bridge because it is a measurable label map on the canonical edge space.  General coherent POVM positivity remains a separate open kernel/disintegration seam; no re-Hahn of coherent cells may replace inherited canonical `dW+`.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
