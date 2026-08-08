from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.smooth_symbol_freezing import sharp_young_constant_3d


def complex_young_efficiency_lower(
    weighted_efficiency: float,
    normalized_symbol_freezing_error: float,
) -> float:
    """Reduce one selected frozen scalar/helical cell to ordinary complex Young.

    Normalize ||f||_(3/2)=||g||_(3/2)=||h||_(3/2)=1 and let A3 be the sharp
    scalar Young constant.  Suppose the physical weighted form T_m obeys

      |T_m| >= R m_* A3,
      |T_m-m0 T_1| <= xi m_* A3,
      |m0| <= m_*.

    Then m_* |T_1| >= |m0 T_1| >= (R-xi)m_* A3, hence
      |T_1|/A3 >= R-xi.

    No lower bound on |m0| is required; near weighted saturation itself prevents
    a tiny frozen symbol unless the freezing error is already large.
    """
    R = float(weighted_efficiency)
    xi = float(normalized_symbol_freezing_error)
    if not (0.0 <= R <= 1.0) or xi < 0 or not math.isfinite(xi):
        raise ValueError("R in [0,1] and finite nonnegative freezing error required")
    return max(0.0, R - xi)


def complex_young_deficit_upper(
    weighted_deficit: float,
    normalized_symbol_freezing_error: float,
) -> float:
    if weighted_deficit < 0 or normalized_symbol_freezing_error < 0:
        raise ValueError("nonnegative deficits required")
    return min(1.0, weighted_deficit + normalized_symbol_freezing_error)


def convolution_pair_efficiency_lower(
    weighted_efficiency: float,
    normalized_symbol_freezing_error: float,
) -> float:
    """The same lower bound holds for ||f*g||_3/A3 by Hölder against h."""
    return complex_young_efficiency_lower(weighted_efficiency, normalized_symbol_freezing_error)


def christ_complex_parent_mark_available(
    *,
    weighted_deficit: float,
    normalized_symbol_freezing_error: float,
    christ_modulus_for_target_distance: float,
) -> bool:
    """Logical interface to Christ's complex-valued Young stability theorem.

    Christ's modulus is external and is intentionally not numerically invented.
    If delta_weighted+xi <= delta_Christ(epsilon), the complex parent pair is
    epsilon-close in L^(3/2) to an extremizing complex Gaussian pair.
    """
    if christ_modulus_for_target_distance <= 0:
        raise ValueError("positive external Christ modulus required")
    return complex_young_deficit_upper(weighted_deficit, normalized_symbol_freezing_error) <= christ_modulus_for_target_distance


def complex_parent_marking_budget(
    *,
    weighted_deficit: float,
    symbol_lipschitz_constant: float,
    relative_cell_diameter: float,
    christ_modulus_for_target_distance: float,
) -> dict[str, float | bool]:
    """Insert the sharp symbol-freezing error A3 L h into the complex Young gate.

    After dividing by A3 and normalized role norms, xi=L h.
    """
    vals = (weighted_deficit, symbol_lipschitz_constant, relative_cell_diameter)
    if any(v < 0 or not math.isfinite(v) for v in vals):
        raise ValueError("finite nonnegative block data required")
    xi = symbol_lipschitz_constant * relative_cell_diameter
    deficit = complex_young_deficit_upper(weighted_deficit, xi)
    return {
        "normalized_symbol_freezing_error": xi,
        "complex_young_deficit_upper": deficit,
        "complex_parent_mark_available": deficit <= christ_modulus_for_target_distance,
    }


def weighted_cell_counterexample_to_magnitude_only_phase_claim() -> dict[str, float]:
    """Toy reminder: magnitude extremality alone does not identify a linear phase.

    f(x)=G(x)e^{i x^2} has |f|=G exactly but is not a Gaussian extremizer with
    affine phase.  The numerical values here are only labels for the exact logical
    counterexample: magnitude distance zero while a nonlinear phase is present.
    """
    return {
        "magnitude_profile_distance": 0.0,
        "nonlinear_phase_quadratic_coefficient": 1.0,
        "lesson": 1.0,
    }


@dataclass(frozen=True)
class ComplexYoungMarkingStress:
    samples: int
    minimum_efficiency_margin: float
    minimum_pair_efficiency_margin: float
    minimum_modulus_gate_margin: float
    magnitude_only_counterexample_phase: float


def stress(samples: int = 50_000, seed: int = 20260808) -> ComplexYoungMarkingStress:
    rng = np.random.default_rng(seed)
    me = mp = mm = float("inf")
    for _ in range(samples):
        R = float(rng.uniform(0.5, 1.0))
        xi = float(rng.uniform(0.0, min(0.2, R)))
        lower = complex_young_efficiency_lower(R, xi)
        exact_algebra = max(0.0, R - xi)
        me = min(me, lower - exact_algebra)
        if abs(lower - exact_algebra) > 2e-15:
            raise AssertionError("weighted-to-complex Young algebra changed")
        pair = convolution_pair_efficiency_lower(R, xi)
        mp = min(mp, pair - exact_algebra)
        if pair + 2e-15 < exact_algebra:
            raise AssertionError("convolution pair efficiency lost trilinear lower bound")

        target_modulus = float(rng.uniform(1e-5, 0.2))
        d = float(rng.uniform(0.0, target_modulus))
        x = float(rng.uniform(0.0, target_modulus - d))
        ok = christ_complex_parent_mark_available(
            weighted_deficit=d,
            normalized_symbol_freezing_error=x,
            christ_modulus_for_target_distance=target_modulus,
        )
        margin = target_modulus - (d + x)
        mm = min(mm, margin)
        if not ok or margin < -2e-15:
            raise AssertionError("Christ modulus interface rejected admissible deficit split")
    cm = weighted_cell_counterexample_to_magnitude_only_phase_claim()
    return ComplexYoungMarkingStress(samples, me, mp, mm, float(cm["nonlinear_phase_quadratic_coefficient"]))


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_WEIGHTED_TO_COMPLEX_YOUNG_PARENT_REDUCTION__CHRIST_MODULUS_EXTERNAL",
        "algebra": "|T_m|>=R m*A3 and |T_m-m0 T1|<=xi m*A3 with |m0|<=m imply |T1|/A3>=R-xi",
        "pair": "||f*g||_3/A3>=R-xi, so the actual complex parent pair is a Young near-extremizer",
        "external_input": "Michael Christ, Near Equality in the Young Inequality for Convolution, Theorem 1.1: complex-valued near-extremizing pairs are close in norm to extremizing pairs",
        "modulus_rule": "choose the physical block deficit plus one normalized symbol-freezing Xi below delta_Christ(epsilon); no numerical Christ modulus is invented",
        "phase_consequence": "for the parent roles, complex Gaussian modulation/phase comes from the complex Young theorem rather than a new phase-persistence assumption",
        "single_charge": "the same symbol-freezing discrepancy is Xi; it is not charged again as a phase defect",
        "countermodel": "magnitude Gaussian proximity alone cannot imply affine phase (e.g. G e^{ix^2}); complex reduction is essential",
        "continuum_status": "requires the recursive selector to expose a scalar/helical frozen-symbol parent cell triple carrying the stated physical weighted efficiency; transfer-cell/material-label alignment remains separate",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-complex-young-parent-marking"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "complex_young_parent_marking.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    A3 = sharp_young_constant_3d()
    md = f"""# Complex Young parent marking after physical symbol freezing\n\nStatus: **{cert['status']}**.\n\nOn one selected scalar/helical frozen-symbol parent cell triple, normalize the three `L^(3/2)` role norms to one.  If the physical weighted transfer obeys\n\n`|T_m| >= R m_* A3`\n\nand the already existing symbol-freezing ledger gives\n\n`|T_m-m0 T_1| <= xi m_* A3`, `|m0|<=m_*`,\n\nthen exactly\n\n`|T_1(f,g,h)|/A3 >= R-xi`,\n\nand hence\n\n`||f*g||_3/A3 >= R-xi`.\n\nHere `A3={A3:.12g}`.  Thus the **actual complex-valued parent pair** `(f,g)` is an ordinary Young near-extremizer once the physical block deficit plus normalized symbol-freezing error is small.\n\nMichael Christ's Theorem 1.1 is stated for complex-valued functions.  Therefore, for every target `epsilon>0`, the external theorem supplies a modulus `delta_Christ(epsilon)>0` such that\n\n`weighted_deficit + xi <= delta_Christ(epsilon)`\n\nimplies that the two complex parent roles are `epsilon`-close in `L^(3/2)` to an extremizing complex Gaussian pair.  No numerical value for that modulus is invented here.\n\nThis closes the parent **phase lift in principle** on frozen selected cells: modulation/phase is part of the complex extremizing pair.  Magnitude proximity alone would not suffice (`G e^(i x^2)` has exactly Gaussian magnitude but nonlinear phase), so the complex reduction is essential.\n\nThe symbol-freezing discrepancy remains the same summable `Xi`; it is not paid again as a phase defect.\n\nStress: `{out.samples}` algebra/modulus states\n- minimum weighted-to-complex efficiency margin: `{out.minimum_efficiency_margin:.3e}`\n- minimum convolution-pair margin: `{out.minimum_pair_efficiency_margin:.3e}`\n- minimum admissible Christ-modulus margin: `{out.minimum_modulus_gate_margin:.3e}`\n\nThe remaining PDE issue is not an abstract phase theorem.  It is to ensure the recursive physical selector exposes the scalar/helical frozen cell triple with the required weighted efficiency and then registers its complex Gaussian parent marks with the same material coherent labels used by causal transfer.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
