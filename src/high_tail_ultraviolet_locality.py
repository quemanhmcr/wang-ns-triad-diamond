from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from src.high_tail_binary_work_reentry import (
    STATUS as HIGH_TAIL_BINARY_STATUS,
    binary_hh_common_work_law,
)
from src.physical_pair_weighted_productivity import (
    physical_work_capacity_constant,
    shell_l32_energy_constant,
)


STATUS = (
    "EXACT_HIGH_TAIL_HH_OUTPUT_SCALE_LAW__ULTRAVIOLET_WORK_TO_PHYSICAL_DTAIL__"
    "CONTINUOUS_LOCALITY_RADIUS_MASS_TRADEOFF__COMPARABLE_REMAINDER_BINARY_READY__"
    "NO_PAIR_COUNT_OR_COHERENT_LOCALITY_SELECTION"
)

DEFAULT_LOCALITY_RADIUS = 2.0
HH_OWNER_FRACTION = 0.5
BALANCED_LOCAL_WORK_FRACTION = 0.25


def highpass_l32_gradient_constant() -> float:
    """C_hp in ||fhat||_(3/2)<=C_hp K^(-1/2)||grad f||_2 for supp fhat>{K}.

    Holder with exponents 4/3 and 4 gives
      int |fhat|^(3/2)
        <= (int |xi|^2|fhat|^2)^(3/4)
           (int_{|xi|>K}|xi|^-6)^(1/4),
    and int_{|xi|>K}|xi|^-6 dxi = 4pi/(3K^3).
    """
    return (4.0 * math.pi / 3.0) ** (1.0 / 6.0)


def ultraviolet_hh_work_constant() -> float:
    """Exact geometric constant 3 sqrt(pi) in the UV HH-work estimate.

    Physical work contributes C_Y=4 A_3 for child derivative <=M.  The child
    hard-shell L^(3/2) mass bound and the two high-pass parent estimates each
    contribute the same C_hp.  Therefore

      C_Y C_hp^3
      = (3 sqrt(3)/2) * sqrt(4pi/3)
      = 3 sqrt(pi).
    """
    C_y = physical_work_capacity_constant(1.0)
    C_hp = highpass_l32_gradient_constant()
    return C_y * C_hp**3


def hh_output_shell_law(positive_hh_shell_common_work: Mapping[int, float]) -> dict[str, object]:
    """Normalize the actual positive HH work pushed only to hard output shells.

    This law is read *before* coherent Hahn refinement.  Its mass at shell j is
      H_j = N int [r_HH,j(t)]_+ dt.
    Coherent atomic positive mass is deliberately not accepted here because it
    may exceed H_j through cancellation.
    """
    if not positive_hh_shell_common_work:
        raise ValueError("nonempty positive HH output-shell law required")
    items = sorted((int(j), float(w)) for j, w in positive_hh_shell_common_work.items())
    if any(j < 1 or w < 0 or not math.isfinite(w) for j, w in items):
        raise ValueError("high-tail shell levels j>=1 with finite nonnegative common work required")
    total = sum(w for _, w in items)
    if total <= 0 or not math.isfinite(total):
        raise ValueError("positive finite HH common work required")
    probs = [(j, w / total) for j, w in items if w > 0]
    j_star, p_star = max(probs, key=lambda jw: jw[1])
    return {
        "total_positive_HH_common_work": total,
        "selected_shell_level": int(j_star),
        "selected_shell_common_work": float(dict(items)[j_star]),
        "p_max": float(p_star),
        "H_inf_output_scale": -math.log(float(p_star)),
        "effective_output_scale_count": 1.0 / float(p_star),
    }


def ultraviolet_common_work_upper(
    child_peak_critical_mass: float,
    physical_tail_dissipation: float,
    selected_shell_level: int,
    locality_radius: float = DEFAULT_LOCALITY_RADIUS,
) -> float:
    """Common-unit positive/absolute HH work with one parent above R M.

    For output |xi|<=M, if one parent is above R M then the other is above
    (R-1)M.  Sharp Young, the hard-child mass bound and the two high-pass
    gradient estimates give

      N W_UV^abs <= [3 sqrt(pi)/sqrt(R(R-1))] sqrt(mu) D_tail.

    Because M/N=2^j and j>=1, requiring (R-1)2^j>=1 makes the lower parent
    high-pass a subset of P_>N, so its scaled gradient budget is <=D_tail.
    The strict resolved transporter S_(M/4)u vanishes identically on these UV
    parent frequencies, hence h=u there and no cutoff multiplier loss appears.
    """
    mu = float(child_peak_critical_mass)
    D = float(physical_tail_dissipation)
    j = int(selected_shell_level)
    R = float(locality_radius)
    if mu < 0 or D < 0 or not all(math.isfinite(x) for x in (mu, D, R)):
        raise ValueError("finite nonnegative mass/dissipation and finite locality radius required")
    if j < 1 or R <= 1.0:
        raise ValueError("high-tail shell j>=1 and locality radius R>1 required")
    if (R - 1.0) * (2.0**j) < 1.0 - 1e-14:
        raise ValueError("lower UV parent cutoff must lie above the parent block scale N")
    C = ultraviolet_hh_work_constant() / math.sqrt(R * (R - 1.0))
    return C * math.sqrt(mu) * D


def high_tail_hh_locality_tradeoff(
    physical_tail_dissipation: float,
    viscosity: float,
    positive_hh_shell_common_work: Mapping[int, float],
    child_peak_mass_by_shell: Mapping[int, float],
    locality_radius: float = DEFAULT_LOCALITY_RADIUS,
) -> dict[str, object]:
    """Continuous locality-radius / child-mass tradeoff for an HH regeneration owner.

    Let H=sum_j H_j be actual positive HH common work and assume the upstream HH
    owner lower H>=nu D_tail/2.  Select the maximal output-shell atom p=H_j/H.
    Split that shell's signed HH source by actual Fourier triads into
      local_R: max(parent frequencies)<=R M,
      UV_R:    max(parent frequencies)> R M.
    Then
      [r_HH,j]_+ <= [r_local_R]_+ + |r_UV_R|,
    and the UV estimate yields the native inequality

      (W_local_R e^Hinf)/D_tail
      + C_R sqrt(mu_peak) e^Hinf >= nu/2,

    C_R=3 sqrt(pi)/sqrt(R(R-1)).

    This is the theorem core.  The balanced quarter alternatives are only a
    readable corollary and do not create new master stop classes.
    """
    D = float(physical_tail_dissipation)
    nu = float(viscosity)
    R = float(locality_radius)
    if D <= 0 or nu <= 0 or not all(math.isfinite(x) for x in (D, nu, R)):
        raise ValueError("positive finite D_tail, viscosity and locality radius required")
    law = hh_output_shell_law(positive_hh_shell_common_work)
    H = float(law["total_positive_HH_common_work"])
    clean_H = HH_OWNER_FRACTION * nu * D
    tol = 5e-13 * max(1e-300, H, clean_H)
    if H + tol < clean_H:
        raise ValueError("positive HH law does not realize the upstream nu D_tail/2 owner lower")

    j = int(law["selected_shell_level"])
    p = float(law["p_max"])
    if j not in child_peak_mass_by_shell:
        raise ValueError("selected output shell is missing its physical peak critical mass")
    mu = float(child_peak_mass_by_shell[j])
    if mu < 0 or not math.isfinite(mu):
        raise ValueError("finite nonnegative selected child peak mass required")

    selected = float(law["selected_shell_common_work"])
    uv = ultraviolet_common_work_upper(mu, D, j, R)
    comparable_lower = max(0.0, selected - uv)
    C_R = ultraviolet_hh_work_constant() / math.sqrt(R * (R - 1.0))
    exp_h = 1.0 / p
    weighted_comparable = comparable_lower * exp_h
    weighted_sqrt_mass = math.sqrt(mu) * exp_h
    lhs = weighted_comparable / D + C_R * weighted_sqrt_mass
    tradeoff_margin = lhs - HH_OWNER_FRACTION * nu
    if tradeoff_margin < -8e-12 * max(1e-300, abs(nu), abs(lhs)):
        raise AssertionError("continuous UV locality/mass tradeoff failed")

    # Balanced readable corollary: one of the two nonnegative terms is >=nu/4.
    clean_weighted_comparable = BALANCED_LOCAL_WORK_FRACTION * nu * D
    clean_weighted_sqrt_mass = (BALANCED_LOCAL_WORK_FRACTION * nu) / C_R
    clean_weighted_mass = clean_weighted_sqrt_mass**2
    weighted_mass = mu * exp_h * exp_h
    owners: list[str] = []
    if weighted_mass + 8e-13 * max(1e-300, weighted_mass, clean_weighted_mass) >= clean_weighted_mass:
        owners.append("critical_child_shell")
    if weighted_comparable + 8e-13 * max(1e-300, weighted_comparable, clean_weighted_comparable) >= clean_weighted_comparable:
        owners.append("comparable_parent_HH_work")
    if not owners:
        raise AssertionError("balanced locality corollary lost both physical owners")

    dyadic_identity = None
    if abs(R - 2.0) <= 1e-14:
        dyadic_identity = {
            "C_R": C_R,
            "clean_entropy_weighted_mass_lower": clean_weighted_mass,
            "exact_clean_mass_formula": nu * nu / (72.0 * math.pi),
            "clean_entropy_weighted_comparable_work_lower": clean_weighted_comparable,
        }

    return {
        **law,
        "locality_radius": R,
        "child_peak_critical_mass": mu,
        "ultraviolet_common_work_upper": uv,
        "comparable_parent_common_work_lower": comparable_lower,
        "locality_constant": C_R,
        "entropy_weighted_comparable_work_lower": weighted_comparable,
        "entropy_weighted_sqrt_child_mass": weighted_sqrt_mass,
        "entropy_weighted_child_mass": weighted_mass,
        "continuous_tradeoff_lhs": lhs,
        "continuous_tradeoff_rhs": HH_OWNER_FRACTION * nu,
        "continuous_tradeoff_margin": tradeoff_margin,
        "clean_entropy_weighted_comparable_work_lower": clean_weighted_comparable,
        "clean_entropy_weighted_child_mass_lower": clean_weighted_mass,
        "joint_clean_owners": tuple(owners),
        "dyadic_R2_corollary": dyadic_identity,
        "comparable_parent_frequency_upper_over_child": R,
        "coherent_refinement_used_for_locality": False,
        "pair_count_used": False,
        "atomic_Hahn_mass_used_for_output_scale_law": False,
        "next_owner_if_critical": "generic_critical_shell_first_stop",
        "next_owner_if_comparable": "comparable_HH_physical_work",
        "master_semantics": "JOINT_RECURSIVE_LOCALITY_OR_CRITICAL_SHELL",
        "status": STATUS,
    }


def comparable_work_binary_reentry(
    selected_shell_level: int,
    actual_comparable_positive_common_work: float,
    comparable_event_atom_arrays: Sequence[np.ndarray],
) -> dict[str, object]:
    """Atomize the already-localized comparable HH source, after locality is proven.

    The Fourier locality split precedes coherent refinement.  Once the restricted
    source has actual positive work W_comp, exact coherent Hahn atomization applies
    to that same bilinear source and yields binary positive mass >=W_comp.
    """
    j = int(selected_shell_level)
    W = float(actual_comparable_positive_common_work)
    if j < 1 or W <= 0 or not math.isfinite(W) or not comparable_event_atom_arrays:
        raise ValueError("positive comparable work, high shell and event atoms required")
    out = binary_hh_common_work_law([j] * len(comparable_event_atom_arrays), comparable_event_atom_arrays)
    reconstructed = float(out["aggregate_positive_HH_common_work"])
    tol = 6e-13 * max(1e-300, W, reconstructed)
    if abs(reconstructed - W) > tol:
        raise ValueError("comparable coherent atoms do not realize the supplied comparable positive HH work")
    if float(out["binary_positive_common_work"]) + tol < W:
        raise AssertionError("comparable binary positive work lost localized physical work")
    return {
        "selected_shell_level": j,
        "actual_comparable_positive_common_work": W,
        "binary_law": out,
        "binary_positive_common_work": float(out["binary_positive_common_work"]),
        "locality_was_read_before_coherent_refinement": True,
        "productivity_gate_supplied": False,
        "next_owner": "binary_comparable_HH_physical_work_law",
    }


def theorem_certificate() -> dict[str, object]:
    C = ultraviolet_hh_work_constant()
    if abs(C - 3.0 * math.sqrt(math.pi)) > 2e-14:
        raise AssertionError("UV geometric constant did not collapse to 3 sqrt(pi)")
    R = 2.0
    C_R = C / math.sqrt(2.0)
    clean_mass_unit_nu = (1.0 / (4.0 * C_R)) ** 2
    if abs(clean_mass_unit_nu - 1.0 / (72.0 * math.pi)) > 2e-14:
        raise AssertionError("dyadic balanced mass constant identity failed")
    return {
        "status": STATUS,
        "upstream": HIGH_TAIL_BINARY_STATUS,
        "output_scale_law": "H_j=N int [r_HH,j(t)]_+dt on hard output shells before coherent Hahn refinement; p_j=H_j/sum H_j",
        "highpass_l32": "supp fhat in {|xi|>K} => ||fhat||_(3/2)<=(4pi/3)^(1/6)K^(-1/2)||grad f||_2",
        "triad_geometry": "for child |xi|<=M, max(parent)>R M forces the other parent >(R-1)M",
        "ultraviolet_bound": "N W_UV^abs <= [3 sqrt(pi)/sqrt(R(R-1))] sqrt(mu_peak) D_tail whenever the lower UV parent cutoff lies above N",
        "native_tradeoff": "W_comp e^(Hinf_out)/D_tail + [3 sqrt(pi)/sqrt(R(R-1))] sqrt(mu_peak)e^(Hinf_out) >= nu/2",
        "dyadic_R2": "R=2 gives balanced corollary mu_peak e^(2Hinf_out)>=nu^2/(72pi) OR W_comp e^(Hinf_out)>=nu D_tail/4; exact ties may satisfy both",
        "representation_rule": "locality is read on the signed PDE hard-shell source before coherent refinement; atomic Hahn positive mass is forbidden as a substitute for H_j",
        "peak_semantics": "on the smooth observed block, t->P_Mu is L2-continuous and the hard-shell critical mass attains its maximum on the compact slab; weaker essential-sup formulations may use epsilon slack",
        "binary_after_locality": "only after comparable physical work is isolated is exact coherent Hahn atomization applied to that restricted source",
        "scope": "no temporal natural-window concentration, signed-good parent progress, Young near-extremality or productivity gate is claimed here",
        "master_rule": "critical-shell and comparable-HH owners are recursive; the balanced split is a corollary of one continuous inequality, not a new stop taxonomy",
    }


@dataclass(frozen=True)
class UltravioletLocalityStress:
    samples: int
    worst_constant_identity_residual: float
    minimum_continuous_tradeoff_margin: float
    minimum_balanced_owner_margin: float
    minimum_binary_localized_dominance_margin: float
    maximum_joint_clean_owner_count: int
    atomic_scale_substitution_failures: int


def _event_atoms(rng: np.random.Generator, signed_positive_work: float) -> np.ndarray:
    h = float(signed_positive_work)
    neg = rng.lognormal(mean=-1.7, sigma=0.7, size=7)
    neg *= float(rng.uniform(0.02, 1.4)) * h / float(neg.sum())
    a = np.empty(8, float)
    a[0] = h + float(neg.sum())
    a[1:] = -neg
    return a.reshape(2, 2, 2)


def stress(samples: int = 50_000, seed: int = 20260810) -> UltravioletLocalityStress:
    rng = np.random.default_rng(seed)
    wc = 0.0
    mt = mb = ma = float("inf")
    max_joint = 0
    substitution_failures = 0
    C = ultraviolet_hh_work_constant()
    wc = abs(C - 3.0 * math.sqrt(math.pi))

    for _ in range(samples):
        D = float(math.exp(rng.uniform(-8.0, 3.0)))
        nu = float(rng.uniform(0.03, 3.0))
        n = int(rng.integers(1, 9))
        levels = np.arange(1, n + 1)
        probs = rng.dirichlet(np.ones(n))
        H = 0.5 * nu * D * float(rng.uniform(1.0, 2.5))
        work = {int(j): float(H * p) for j, p in zip(levels, probs)}
        law = hh_output_shell_law(work)
        j = int(law["selected_shell_level"])
        p = float(law["p_max"])
        R = float(rng.uniform(1.5, 4.0))
        C_R = C / math.sqrt(R * (R - 1.0))
        balance_mu = (p * nu / (4.0 * C_R)) ** 2
        mu = balance_mu * float(math.exp(rng.uniform(-4.0, 4.0)))
        masses = {int(k): float(math.exp(rng.uniform(-10.0, 1.0))) for k in levels}
        masses[j] = mu
        out = high_tail_hh_locality_tradeoff(D, nu, work, masses, R)
        mt = min(mt, float(out["continuous_tradeoff_margin"]))
        max_joint = max(max_joint, len(tuple(out["joint_clean_owners"])))

        weighted_mass = float(out["entropy_weighted_child_mass"])
        clean_mass = float(out["clean_entropy_weighted_child_mass_lower"])
        weighted_work = float(out["entropy_weighted_comparable_work_lower"])
        clean_work = float(out["clean_entropy_weighted_comparable_work_lower"])
        mass_relative_margin = weighted_mass / max(clean_mass, 1e-300) - 1.0
        work_relative_margin = weighted_work / max(clean_work, 1e-300) - 1.0
        owner_margin = max(mass_relative_margin, work_relative_margin)
        mb = min(mb, owner_margin)
        if owner_margin < -5e-11:
            raise AssertionError("balanced UV locality owner corollary failed")

        # On comparable-owner samples, construct the restricted source only after
        # Fourier locality has been read, then atomize its actual positive work.
        if "comparable_parent_HH_work" in tuple(out["joint_clean_owners"]):
            Wcomp = max(float(out["comparable_parent_common_work_lower"]), 1e-14)
            parts = Wcomp * rng.dirichlet(np.ones(int(rng.integers(1, 5))))
            atoms = tuple(_event_atoms(rng, float(x)) for x in parts)
            b = comparable_work_binary_reentry(j, Wcomp, atoms)
            ma = min(ma, float(b["binary_positive_common_work"]) - Wcomp)

        # Representation guard: a fictitious atomic-positive shell law generally
        # changes p_max; locality must not consume it.  Count such differences as
        # evidence that the two laws are genuinely distinct observables.
        fake = {k: v * float(rng.uniform(1.0, 5.0)) for k, v in work.items()}
        fake_law = hh_output_shell_law(fake)
        if abs(float(fake_law["p_max"]) - p) > 1e-10:
            substitution_failures += 1

    if not math.isfinite(ma):
        ma = 0.0
    return UltravioletLocalityStress(samples, wc, mt, mb, ma, max_joint, substitution_failures)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-tail-ultraviolet-locality"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "high_tail_ultraviolet_locality.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# High-tail HH ultraviolet locality from physical dissipation

Status: **{cert['status']}**.

Start only after the common-unit theorem has identified HH as a physical regeneration owner, so

`H=sum_j H_j >= nu D_tail/2`,

where `H_j=N int [r_HH,j(t)]_+ dt` is the actual positive HH child-work of hard output shell `M_j=2^jN`.  This output-scale law is read **before coherent Hahn refinement**.

Normalize `p_j=H_j/H`, choose the maximal shell `p_*`, and write `H_inf^out=-log p_*`.  For that selected child shell, let `mu_peak=max_t M||P_Mu(t)||_2^2`.

If one HH parent lies above `R M`, Fourier triad closure forces the other above `(R-1)M`.  The exact high-pass estimate

`||fhat||_(3/2) <= (4pi/3)^(1/6) K^(-1/2)||grad f||_2`

combines with sharp physical Young and the child hard-shell mass bound.  The constants collapse exactly:

`C_Y C_hp^3 = 3 sqrt(pi)`.

Therefore

`N W_UV^abs <= [3 sqrt(pi)/sqrt(R(R-1))] sqrt(mu_peak) D_tail`.

Since `[r_HH]_+ <= [r_comp]_+ + |r_UV|`, the native theorem is the continuous relation

`W_comp exp(H_inf^out)/D_tail + [3 sqrt(pi)/sqrt(R(R-1))] sqrt(mu_peak) exp(H_inf^out) >= nu/2`.

No locality threshold is fundamental here.  The balanced quarter split is only a readable corollary.  For `R=2`,

`mu_peak exp(2 H_inf^out) >= nu^2/(72 pi)`

or

`W_comp exp(H_inf^out) >= nu D_tail/4`,

with exact ties allowed to satisfy both.  The first owner enters the generic critical-shell theorem.  On the second owner, the comparable source has both parent frequencies at most `2M`; **only then** is coherent Hahn atomization applied to that restricted physical source, producing binary positive work at least as large as its aggregate comparable positive work.

The theorem uses no parent-pair count, no coherent-cell locality selector, no signed-good ratio, no Young near-extremality and no generated-energy productivity gate.  It also does not solve temporal natural-window concentration; that is the next, logically separate seam.

Stress: `{out.samples}` output-scale/locality/binary-ready states
- worst `3 sqrt(pi)` constant identity residual: `{out.worst_constant_identity_residual:.3e}`
- minimum continuous tradeoff margin: `{out.minimum_continuous_tradeoff_margin:.3e}`
- minimum balanced-owner margin: `{out.minimum_balanced_owner_margin:.3e}`
- minimum localized binary positive-dominance margin: `{out.minimum_binary_localized_dominance_margin:.3e}`
- maximum joint clean owner count: `{out.maximum_joint_clean_owner_count}`
- sampled cases where fake atomic reweighting changed output-scale p_max: `{out.atomic_scale_substitution_failures}`

No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
