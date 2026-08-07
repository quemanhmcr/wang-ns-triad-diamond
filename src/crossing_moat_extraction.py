from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

from .single_edge_certificate import RSTAR_LO, RSTAR_HI
from .triad_extremizer import symmetric_gamma, symmetric_rstar

ETA0 = Fraction(1, 10_000)
GAP_DEV = Fraction(1, 80)          # sqrt(eta)+25 eta at eta=1e-4
PARENT_IMBAL = Fraction(1, 200)    # u <= 50 eta
SMOOTH_DELTA = Fraction(1, 20)
SHELL_TARGET = Fraction(2, 25)
MOAT_TARGET = Fraction(1, 25)
BINS = 4
PHYSICAL_HODGE_BASE = Fraction(25, 106)
CROSSING_HODGE_COEFF = PHYSICAL_HODGE_BASE / BINS


@dataclass(frozen=True)
class CrossingEdge:
    lower_parent_log: float
    top_parent_log: float
    child_log: float
    child_transfer: float
    hodge_residual: float = 0.0

    @property
    def gap(self) -> float:
        return self.child_log - self.top_parent_log

    @property
    def imbalance(self) -> float:
        return self.top_parent_log - self.lower_parent_log

    @property
    def midgap(self) -> float:
        return 0.5 * (self.top_parent_log + self.child_log)


@dataclass(frozen=True)
class BinCore:
    index: int
    transfer: float
    hodge_numerator: float
    parent_center: float
    child_center: float
    max_parent_deviation: float
    max_child_deviation: float
    edge_count: int


@dataclass(frozen=True)
class CrossingExtraction:
    total_transfer: float
    total_hodge_numerator: float
    transfer_core: BinCore
    hodge_core: BinCore
    theoretical_shell_halfwidth: float
    theoretical_moat_margin: float


def theoretical_shell_halfwidth(gamma: float) -> float:
    a = float(GAP_DEV)
    u = float(PARENT_IMBAL)
    # Midgaps of crossing edges lie in an interval of length gamma+a. Four
    # equal bins have halfwidth (gamma+a)/8. Moving from the bin center to the
    # ideal parent/child shell adds a/2, and the lower parent adds u.
    return gamma / 8.0 + 5.0 * a / 8.0 + u


def theoretical_moat_margin(gamma: float) -> float:
    sigma = theoretical_shell_halfwidth(gamma)
    return gamma / 2.0 - 2.0 * sigma - float(SMOOTH_DELTA)


def _validate_edge(e: CrossingEdge, tau0: float, gamma: float) -> None:
    if not (e.child_transfer >= 0.0 and e.hodge_residual >= 0.0):
        raise ValueError("weights/residuals must be nonnegative")
    if not (e.lower_parent_log <= e.top_parent_log <= tau0 <= e.child_log):
        raise ValueError("edge must be ordered and cross the reference cut")
    if e.imbalance > float(PARENT_IMBAL) + 1e-12:
        raise ValueError("parent imbalance exceeds signed-good theorem bound")
    if abs(e.gap - gamma) > float(GAP_DEV) + 1e-12:
        raise ValueError("parent-child gap exceeds signed-good theorem bound")


def extract_four_bin_core(edges: Sequence[CrossingEdge], tau0: float, gamma: float | None = None) -> CrossingExtraction:
    if not edges:
        raise ValueError("need at least one crossing edge")
    gamma = symmetric_gamma(symmetric_rstar()) if gamma is None else float(gamma)
    for e in edges:
        _validate_edge(e, tau0, gamma)

    total_transfer = sum(e.child_transfer for e in edges)
    if total_transfer <= 0.0:
        raise ValueError("positive child-transfer mass is required")
    total_hodge = sum(e.child_transfer * e.hodge_residual for e in edges)

    a = float(GAP_DEV)
    lo = tau0 - 0.5 * (gamma + a)
    hi = tau0 + 0.5 * (gamma + a)
    width = (hi - lo) / BINS
    bins: list[list[CrossingEdge]] = [[] for _ in range(BINS)]
    for e in edges:
        if not (lo - 2e-12 <= e.midgap <= hi + 2e-12):
            raise AssertionError("crossing geometry failed to confine the midgap")
        idx = min(BINS - 1, max(0, int((e.midgap - lo) / width)))
        bins[idx].append(e)

    def core(idx: int) -> BinCore:
        bucket = bins[idx]
        c = lo + (idx + 0.5) * width
        pc = c - gamma / 2.0
        cc = c + gamma / 2.0
        pdev = 0.0
        cdev = 0.0
        for e in bucket:
            pdev = max(pdev, abs(e.top_parent_log - pc), abs(e.lower_parent_log - pc))
            cdev = max(cdev, abs(e.child_log - cc))
        return BinCore(
            index=idx,
            transfer=sum(e.child_transfer for e in bucket),
            hodge_numerator=sum(e.child_transfer * e.hodge_residual for e in bucket),
            parent_center=pc,
            child_center=cc,
            max_parent_deviation=pdev,
            max_child_deviation=cdev,
            edge_count=len(bucket),
        )

    cores = [core(i) for i in range(BINS)]
    transfer_core = max(cores, key=lambda x: x.transfer)
    hodge_core = max(cores, key=lambda x: x.hodge_numerator)
    sigma = theoretical_shell_halfwidth(gamma)
    moat = theoretical_moat_margin(gamma)

    if transfer_core.transfer + 1e-12 < total_transfer / BINS:
        raise AssertionError("four-bin transfer pigeonhole failed")
    if hodge_core.hodge_numerator + 1e-12 < total_hodge / BINS:
        raise AssertionError("four-bin Hodge pigeonhole failed")
    for c in cores:
        if c.edge_count and max(c.max_parent_deviation, c.max_child_deviation) > sigma + 2e-12:
            raise AssertionError("selected shell exceeded theorem halfwidth")

    return CrossingExtraction(total_transfer, total_hodge, transfer_core, hodge_core, sigma, moat)


def arb_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover - Actions certificate
        raise RuntimeError("python-flint is required for rigorous certification") from exc
    ctx.prec = 160
    lo = arb(f"{RSTAR_LO.numerator}/{RSTAR_LO.denominator}")
    hi = arb(f"{RSTAR_HI.numerator}/{RSTAR_HI.denominator}")
    r = lo.union(hi)
    gamma = -r.log()
    a = arb(1) / 80
    u = arb(1) / 200
    delta = arb(1) / 20
    sigma = gamma / 8 + 5 * a / 8 + u
    moat = gamma / 2 - 2 * sigma - delta
    if not (sigma < arb(2) / 25):
        raise AssertionError(f"crossing shell not inside certified shell: {sigma}")
    if not (moat > arb(1) / 25):
        raise AssertionError(f"crossing smooth moat too small: {moat}")
    return {
        "rstar_ball": str(r),
        "gamma_ball": str(gamma),
        "shell_halfwidth_ball": str(sigma),
        "shell_target": "2/25",
        "smooth_delta": "1/20",
        "moat_margin_ball": str(moat),
        "moat_lower_bound": "1/25",
        "transfer_fraction": "1/4",
        "hodge_numerator_fraction": "1/4",
        "crossing_physical_hodge_coefficient": "25/424",
        "status": "CERTIFIED",
    }


def stress(samples: int = 50_000, seed: int = 20260807) -> dict[str, float]:
    import numpy as np

    rng = np.random.default_rng(seed)
    gamma = symmetric_gamma(symmetric_rstar())
    worst_transfer = 1.0
    worst_hodge = 1.0
    worst_shell_slack = float("inf")
    for _ in range(samples):
        n = int(rng.integers(4, 40))
        tau0 = float(rng.normal())
        rows = []
        for _j in range(n):
            gap = gamma + rng.uniform(-float(GAP_DEV), float(GAP_DEV))
            # crossing p<=tau<=q means p-tau in [-gap,0]
            p = tau0 + rng.uniform(-gap, 0.0)
            u = rng.uniform(0.0, float(PARENT_IMBAL))
            q = p + gap
            rows.append(CrossingEdge(p - u, p, q, float(rng.exponential(1.0)), float(rng.exponential(1.0))))
        ex = extract_four_bin_core(rows, tau0, gamma)
        worst_transfer = min(worst_transfer, ex.transfer_core.transfer / ex.total_transfer)
        if ex.total_hodge_numerator > 0:
            worst_hodge = min(worst_hodge, ex.hodge_core.hodge_numerator / ex.total_hodge_numerator)
        worst_shell_slack = min(
            worst_shell_slack,
            ex.theoretical_shell_halfwidth - max(ex.transfer_core.max_parent_deviation, ex.transfer_core.max_child_deviation),
        )
    return {
        "samples": samples,
        "worst_transfer_fraction": worst_transfer,
        "worst_hodge_numerator_fraction": worst_hodge,
        "minimum_shell_slack": worst_shell_slack,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-crossing-moat"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_certificate()
    probe = stress(args.samples)
    out = {"certificate": cert, "stress": probe}
    (args.outdir / "crossing_moat_extraction.json").write_text(json.dumps(out, indent=2))
    md = f"""# Crossing-to-common-moat extraction\n\nStatus: **{cert['status']}**.\n\n- crossing transfer fraction: `>= {cert['transfer_fraction']}`\n- crossing Hodge numerator fraction: `>= {cert['hodge_numerator_fraction']}`\n- shell halfwidth enclosure: `{cert['shell_halfwidth_ball']}` < `{cert['shell_target']}`\n- smooth moat enclosure: `{cert['moat_margin_ball']}` > `{cert['moat_lower_bound']}`\n- conservative crossing physical-Hodge coefficient: `{cert['crossing_physical_hodge_coefficient']}`\n- stress samples: `{probe['samples']}`\n- worst transfer fraction seen: `{probe['worst_transfer_fraction']:.9f}`\n- worst Hodge fraction seen: `{probe['worst_hodge_numerator_fraction']:.9f}`\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
