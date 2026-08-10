from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.high_frequency_dissipation_reentry import (
    STATUS as HIGH_TAIL_STATUS,
    high_tail_energy_owner_threshold,
    positive_shell_work_disintegration,
)
from src.recursive_coherent_witness_extraction import binary_work_ledger
from src.resolved_interface_donor_quotient import STATUS as RESOLVED_INTERFACE_DONOR_STATUS


STATUS = (
    "EXACT_HIGH_TAIL_REGENERATION_COMMON_UNIT_CAUSAL_OWNERS__"
    "HH_TO_BINARY_PHYSICAL_WORK_LAW__OWN_SCALE_REWEIGHTING_NOT_CAUSAL__"
    "PRODUCTIVITY_AND_LOCALITY_REMAIN_CONDITIONAL"
)


def common_unit_regeneration_owners(
    physical_tail_dissipation_lower: float,
    viscosity: float,
    positive_scaled_tail_work: float,
    positive_scaled_shell_works: Sequence[float],
    positive_hh_common_work: float,
    positive_interface_common_work: float,
) -> dict[str, object]:
    """Route high-tail regeneration using only the common physical work unit N dW.

    The upstream regeneration owner supplies

        N W_>^+ >= nu D_tail.

    Orthogonal hard-shell disintegration gives

        sum_j N W_j^+ >= N W_>^+.

    At each shell low--low work is support-excluded and signed shell work is HH
    plus resolved mixed/interface work.  Taking positive parts and summing gives

        H_N^+ + I_N^+ >= sum_j N W_j^+.

    Therefore HH or interface common-unit work carries at least half of the
    realized shell positive work and, cleanly, at least nu D_tail/2.  Exact ties
    are retained jointly.  No factor M_j/N is used in this causal owner split.
    """
    D = float(physical_tail_dissipation_lower)
    nu = float(viscosity)
    W = float(positive_scaled_tail_work)
    H = float(positive_hh_common_work)
    I = float(positive_interface_common_work)
    if min(D, nu) <= 0 or min(W, H, I) < 0 or not all(math.isfinite(x) for x in (D, nu, W, H, I)):
        raise ValueError("finite positive tail data and nonnegative work owners required")
    clean_tail = high_tail_energy_owner_threshold(D, nu)
    tol_tail = 4e-13 * max(1.0, clean_tail, W)
    if W + tol_tail < clean_tail:
        raise ValueError("positive nonlinear regeneration owner does not carry nu D_tail")

    shell = positive_shell_work_disintegration(W, positive_scaled_shell_works)
    S = float(shell["scaled_shell_positive_work_sum"])
    tol = 5e-13 * max(1.0, S, H + I)
    if H + I + tol < S:
        raise ValueError("common-unit HH/interface work does not cover shell positive work")

    threshold = 0.5 * S
    clean_threshold = 0.5 * clean_tail
    owners: list[str] = []
    if H + tol >= threshold:
        owners.append("positive_HH_regeneration")
    if I + tol >= threshold:
        owners.append("positive_resolved_cross_interface")
    if not owners:
        raise AssertionError("common-unit HH/interface owner pigeonhole failed")
    if threshold + tol < clean_threshold:
        raise AssertionError("realized common-unit threshold lost the clean nu D_tail/2 lower")

    return {
        "scaled_tail_positive_work": W,
        "scaled_shell_positive_work_sum": S,
        "clean_tail_regeneration_lower": clean_tail,
        "owner_threshold": threshold,
        "clean_owner_threshold": clean_threshold,
        "joint_owners": tuple(owners),
        "positive_HH_common_work": H,
        "positive_interface_common_work": I,
        "own_scale_shell_work_diagnostic": float(shell["own_scale_positive_shell_work"]),
        "own_scale_diagnostic_clean_lower": float(shell["clean_own_scale_lower"]),
        "causal_probability_reweighted_by_shell_scale": False,
        "interface_is_free": False,
        "HH_is_productivity_generated_branch": False,
    }


def _validate_atom_events(
    event_shell_levels: Sequence[int],
    event_atom_arrays: Sequence[np.ndarray],
) -> tuple[tuple[int, ...], tuple[np.ndarray, ...]]:
    levels = tuple(int(j) for j in event_shell_levels)
    arrays = tuple(np.asarray(a, float) for a in event_atom_arrays)
    if not arrays or len(levels) != len(arrays):
        raise ValueError("matching nonempty shell-event atom family required")
    if any(j < 1 for j in levels):
        raise ValueError("high-tail shell levels must satisfy j>=1, M_j/N=2^j")
    for a in arrays:
        if a.ndim != 3 or not np.all(np.isfinite(a)):
            raise ValueError("every event requires a finite three-index coherent work array")
    return levels, arrays


def binary_hh_common_work_law(
    event_shell_levels: Sequence[int],
    event_atom_arrays: Sequence[np.ndarray],
) -> dict[str, object]:
    """Exact positive binary law from time/event-resolved signed HH work atoms.

    For each physical shell-time event the coherent atom theorem gives

        sum_CDE W_CDE = W_HH(event),
        sum_CDE [W_CDE]_+ >= [W_HH(event)]_+.

    Inputs here are already in the **same common work unit**, canonically N times
    physical child-energy work.  Multiplying every event by the common factor N
    leaves the normalized causal law unchanged.  Shell-dependent M_j factors are
    not inserted into the probabilities.
    """
    levels, arrays = _validate_atom_events(event_shell_levels, event_atom_arrays)
    aggregate_positive = 0.0
    atomic_positive = 0.0
    negative_backscatter = 0.0
    worst_cancel = 0.0
    shell_positive: dict[int, float] = {}
    raw_events: list[dict[str, object]] = []

    for event_index, (level, atoms) in enumerate(zip(levels, arrays)):
        led = binary_work_ledger(atoms)
        aggregate_positive += float(led.aggregate_positive_work)
        atomic_positive += float(led.positive_transfer_mass)
        negative_backscatter += float(led.negative_backscatter_mass)
        worst_cancel = max(worst_cancel, abs(float(led.cancellation_residual)))
        pos = np.maximum(atoms, 0.0)
        event_pos = float(pos.sum())
        shell_positive[level] = shell_positive.get(level, 0.0) + event_pos
        for i, j, k in np.argwhere(pos > 0):
            raw_events.append(
                {
                    "event": int(event_index),
                    "shell_level": int(level),
                    "parent1": int(i),
                    "parent2": int(j),
                    "child": int(k),
                    "mass": float(pos[i, j, k]),
                }
            )

    if atomic_positive <= 0:
        raise ValueError("positive binary physical work mass required")
    tol = 5e-13 * max(1.0, aggregate_positive, atomic_positive)
    if atomic_positive + tol < aggregate_positive:
        raise AssertionError("binary positive physical work lost aggregate positive HH work")

    for row in raw_events:
        row["probability"] = float(row["mass"]) / atomic_positive
    probability_residual = sum(float(row["probability"]) for row in raw_events) - 1.0
    if abs(probability_residual) > 3e-13:
        raise AssertionError("common-unit binary causal law failed to normalize")

    scale_total = sum(shell_positive.values())
    if abs(scale_total - atomic_positive) > tol:
        raise AssertionError("shell pushforward lost binary positive work")
    level_star, shell_star = max(sorted(shell_positive.items()), key=lambda kv: kv[1])
    pmax = shell_star / atomic_positive
    Hinf = -math.log(pmax)
    ownscale_selected = (2.0**level_star) * shell_star
    hinf_weighted_ownscale = ownscale_selected / pmax
    clean_two_atomic = 2.0 * atomic_positive
    if hinf_weighted_ownscale + tol < clean_two_atomic:
        raise AssertionError("selected-shell own-scale diagnostic lost the j>=1 factor two")

    return {
        "aggregate_positive_HH_common_work": aggregate_positive,
        "binary_positive_common_work": atomic_positive,
        "negative_binary_backscatter_common_work": negative_backscatter,
        "atomic_positive_dominance_margin": atomic_positive - aggregate_positive,
        "worst_event_cancellation_residual": worst_cancel,
        "events": tuple(raw_events),
        "probability_residual": probability_residual,
        "shell_positive_binary_common_work": dict(sorted(shell_positive.items())),
        "selected_shell_level": int(level_star),
        "selected_shell_probability": pmax,
        "H_inf_binary_scale": Hinf,
        "selected_shell_own_scale_binary_work": ownscale_selected,
        "H_inf_weighted_selected_own_scale_binary_work": hinf_weighted_ownscale,
        "clean_two_times_binary_common_work": clean_two_atomic,
        "causal_probability_uses_common_N_work": True,
        "causal_probability_uses_Mj_reweighting": False,
    }


def high_tail_hh_binary_reentry(
    physical_tail_dissipation_lower: float,
    viscosity: float,
    positive_scaled_tail_work: float,
    positive_scaled_shell_works: Sequence[float],
    positive_hh_common_work: float,
    positive_interface_common_work: float,
    event_shell_levels: Sequence[int],
    event_atom_arrays: Sequence[np.ndarray],
) -> dict[str, object]:
    """Compose high-tail regeneration ownership with exact HH binary atomization.

    If HH is a primary common-unit owner, its time/event-resolved positive work is
    exactly the aggregate positive work reconstructed by the atom family.  Atomic
    positive mass therefore carries at least nu D_tail/2.  The normalized atom law
    is a genuine binary physical child-work causal law.  No statement is made that
    its parents are comparable-scale, Young-good, or productivity-generating.
    """
    route = common_unit_regeneration_owners(
        physical_tail_dissipation_lower,
        viscosity,
        positive_scaled_tail_work,
        positive_scaled_shell_works,
        positive_hh_common_work,
        positive_interface_common_work,
    )
    H = float(positive_hh_common_work)
    binary: dict[str, object] | None = None
    if H > 0:
        binary = binary_hh_common_work_law(event_shell_levels, event_atom_arrays)
        reconstructed = float(binary["aggregate_positive_HH_common_work"])
        tol = 6e-13 * max(1.0, H, reconstructed)
        if abs(reconstructed - H) > tol:
            raise ValueError("time/event-resolved coherent atoms do not realize the supplied positive HH work")

    hh_primary = "positive_HH_regeneration" in tuple(route["joint_owners"])
    clean_binary = 0.5 * float(viscosity) * float(physical_tail_dissipation_lower)
    if hh_primary:
        if binary is None:
            raise AssertionError("HH primary owner requires a positive binary work law")
        P = float(binary["binary_positive_common_work"])
        diag = float(binary["H_inf_weighted_selected_own_scale_binary_work"])
        tol = 7e-13 * max(1.0, clean_binary, P, diag)
        if P + tol < clean_binary:
            raise AssertionError("HH primary owner lost the clean binary positive-work lower")
        if diag + tol < 2.0 * clean_binary:
            raise AssertionError("binary scale diagnostic lost the clean nu D_tail lower")

    return {
        "owner_route": route,
        "binary_HH_law": binary,
        "HH_primary": hh_primary,
        "clean_binary_positive_common_work_if_HH_owner": clean_binary,
        "clean_Hinf_weighted_selected_ownscale_work_if_HH_owner": 2.0 * clean_binary,
        "next_owner_if_interface": RESOLVED_INTERFACE_DONOR_STATUS,
        "next_owner_if_HH": "binary_physical_child_work_law",
        "productivity_energy_gate_supplied": False,
        "Young_near_extremality_supplied": False,
        "parent_child_scale_locality_supplied": False,
        "master_semantics": "JOINT_PHYSICAL_WORK_OWNERS__INTERFACE_QUOTIENT_BEFORE_RECURSION",
        "status": STATUS,
    }


def own_scale_reweighting_counterexample() -> dict[str, object]:
    """Equal physical causes become unequal if one incorrectly weights by M_j/N."""
    levels = (1, 10)
    common = np.array([1.0, 1.0])
    common_prob = common / common.sum()
    own = common * np.array([2.0**levels[0], 2.0**levels[1]])
    own_prob = own / own.sum()
    return {
        "shell_levels": levels,
        "common_unit_probabilities": tuple(float(x) for x in common_prob),
        "own_scale_reweighted_probabilities": tuple(float(x) for x in own_prob),
        "maximum_probability_distortion": float(np.max(np.abs(common_prob - own_prob))),
        "lesson": "M_j/N weighting changes causal probabilities across shells and is diagnostic only",
    }


def theorem_certificate() -> dict[str, object]:
    cex = own_scale_reweighting_counterexample()
    if float(cex["maximum_probability_distortion"]) < 0.49:
        raise AssertionError("scale-reweighting counterexample is not sufficiently separating")
    return {
        "status": STATUS,
        "upstream": HIGH_TAIL_STATUS,
        "common_unit": "all cross-shell causal weights are N times actual positive child-energy work; the common factor N cancels under normalization",
        "owner_cover": "regeneration owner N W_>^+>=nu D_tail; sum_j N W_j^+>=N W_>^+; low-low exclusion gives H_N^+ + I_N^+ >= sum_j N W_j^+",
        "clean_owner": "HH or resolved interface common-unit work carries at least nu D_tail/2, exact ties joint",
        "binary_atomization": "for each shell-time HH event, exact coherent work atoms satisfy P-N=W_HH and P>=[W_HH]_+; after integration binary positive common work dominates positive HH common work",
        "causal_law": "normalize positive coherent HH atoms in their common N dW unit; every event has exactly two parent coherent labels and one child label",
        "anti_reweight": "M_j/N weighting may be used only to read own-scale strength after the causal law exists; it must not redefine probabilities across shells",
        "scale_diagnostic": "for the common-unit binary law on high shells j>=1, selecting the largest shell atom gives ownscale_selected*exp(H_inf_scale)>=2*binary_positive_common_work; on an HH owner this is >=nu D_tail",
        "productivity_scope": "no W_HH>=8E1/15 child-energy generation gate is assumed; KL/log-productivity remains conditional on that separate gate",
        "locality_scope": "no parent/child scale-comparability or Young near-extremality is inferred from generic high-tail HH work; nonlocal K>>M geometry remains a separate physical seam",
        "interface_continuation": "the resolved-interface owner delegates to the resolved donor/circulation quotient: symmetric work is existing strain/deformation provenance, while skew work is same-event conservative donor tracing with no new recursive generation",
        "master_rule": "binary-HH remains a recursive physical work owner; resolved interface is first quotiented into existing strain/deformation or same-event conservative donor provenance; exact ties remain joint with no additive reset or lexicographic priority",
        "counterexample": cex,
    }


@dataclass(frozen=True)
class HighTailBinaryWorkStress:
    samples: int
    minimum_common_owner_margin: float
    minimum_clean_common_owner_margin: float
    minimum_atomic_positive_dominance_margin: float
    worst_binary_probability_residual: float
    minimum_HH_binary_clean_margin: float
    minimum_Hinf_ownscale_diagnostic_margin: float
    worst_common_scale_probability_invariance_residual: float
    maximum_joint_owner_count: int
    scale_reweighting_counterexample_distortion: float


def _event_atoms_with_positive_signed_work(rng: np.random.Generator, signed_positive_work: float) -> np.ndarray:
    h = float(signed_positive_work)
    if h <= 0:
        raise ValueError("positive event HH work required")
    neg = rng.lognormal(mean=-1.5, sigma=0.8, size=7)
    neg *= float(rng.uniform(0.05, 2.0)) * h / float(neg.sum())
    a = np.empty(8, float)
    a[0] = h + float(neg.sum())
    a[1:] = -neg
    return a.reshape(2, 2, 2)


def stress(samples: int = 50_000, seed: int = 20260810) -> HighTailBinaryWorkStress:
    rng = np.random.default_rng(seed)
    mo = mc = ma = mb = md = float("inf")
    wp = wi = 0.0
    max_joint = 0

    for _ in range(samples):
        D = float(math.exp(rng.uniform(-8.0, 3.0)))
        nu = float(rng.uniform(0.03, 3.0))
        clean = nu * D
        W = clean * float(rng.uniform(1.0, 2.5))
        n_shell = int(rng.integers(1, 9))
        shell_sum = W * float(rng.uniform(1.0, 1.8))
        shell_work = shell_sum * rng.dirichlet(np.ones(n_shell))

        total_owner_cover = shell_sum * float(rng.uniform(1.0, 1.8))
        frac = float(rng.uniform(0.02, 0.98))
        H = frac * total_owner_cover
        I = (1.0 - frac) * total_owner_cover
        route = common_unit_regeneration_owners(D, nu, W, shell_work, H, I)
        threshold = float(route["owner_threshold"])
        mo = min(mo, max(H, I) - threshold)
        mc = min(mc, threshold - 0.5 * clean)
        max_joint = max(max_joint, len(tuple(route["joint_owners"])))

        n_event = int(rng.integers(1, 7))
        hparts = H * rng.dirichlet(np.ones(n_event))
        levels = tuple(int(rng.integers(1, n_shell + 1)) for _ in range(n_event))
        atoms = tuple(_event_atoms_with_positive_signed_work(rng, float(h)) for h in hparts)
        binary = binary_hh_common_work_law(levels, atoms)
        ma = min(ma, float(binary["atomic_positive_dominance_margin"]))
        wp = max(wp, abs(float(binary["probability_residual"])))

        lam = float(math.exp(rng.uniform(-7.0, 7.0)))
        scaled = binary_hh_common_work_law(levels, tuple(lam * a for a in atoms))
        p0 = [float(e["probability"]) for e in tuple(binary["events"])]
        p1 = [float(e["probability"]) for e in tuple(scaled["events"])]
        if len(p0) != len(p1):
            raise AssertionError("common work-unit scaling changed binary event support")
        inv = max((abs(a - b) for a, b in zip(p0, p1)), default=0.0)
        wi = max(wi, inv)
        if inv > 3e-13:
            raise AssertionError("common work-unit scaling changed causal probabilities")

        out = high_tail_hh_binary_reentry(D, nu, W, shell_work, H, I, levels, atoms)
        if out["HH_primary"]:
            b = out["binary_HH_law"]
            assert isinstance(b, dict)
            clean_binary = float(out["clean_binary_positive_common_work_if_HH_owner"])
            mb = min(mb, float(b["binary_positive_common_work"]) - clean_binary)
            target = float(out["clean_Hinf_weighted_selected_ownscale_work_if_HH_owner"])
            md = min(md, float(b["H_inf_weighted_selected_own_scale_binary_work"]) - target)

    if not math.isfinite(mb):
        mb = 0.0
    if not math.isfinite(md):
        md = 0.0
    cex = own_scale_reweighting_counterexample()
    return HighTailBinaryWorkStress(
        samples,
        mo,
        mc,
        ma,
        wp,
        mb,
        md,
        wi,
        max_joint,
        float(cex["maximum_probability_distortion"]),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-tail-binary-work-reentry"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "high_tail_binary_work_reentry.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# High-tail regeneration -> common-unit binary physical work

Status: **{cert['status']}**.

The high-tail energy theorem already supplies a physical regeneration owner

`N W_>^+ >= nu D_tail`.

The causal law must stay in this common physical unit.  Orthogonal hard-shell work gives

`sum_j N W_j^+ >= N W_>^+`,

and low--low support exclusion at each shell gives

`H_N^+ + I_N^+ >= sum_j N W_j^+`.

Hence HH or resolved interface common-unit work carries at least `nu D_tail/2`; exact ties remain joint.  The stronger own-scale quantity `sum_j M_j W_j^+` remains available as a diagnostic, but it is **not** used to redefine causal probabilities across shells.

The resolved-interface continuation is now the donor/circulation quotient.  In the same common `N dW` unit its actual resolved operator splits as `K+S`, so an interface owner gives either conservative skew donor work or existing symmetric strain/deformation work at least `nu D_tail/4`.  Skew role traversal is same-event provenance, not a new recursive generation; no shell-scale reweighting is used.

On the HH owner, exact coherent atomization is already enough.  Event by event,

`sum_CDE W_CDE = W_HH(event)`,

and the Hahn split gives

`sum_CDE [W_CDE]_+ >= [W_HH(event)]_+`.

Therefore the integrated positive coherent atoms carry at least `nu D_tail/2` common-unit physical child work whenever HH is a primary owner.  Normalizing those atoms gives a genuine binary causal law with two parent coherent labels and one child label.  No generated-energy gate is needed merely to create this law.

The distinction is essential.  Two equal physical work atoms at shell levels `1` and `10` have common-unit probabilities `(1/2,1/2)`.  Multiplying them by `M_j/N=2^j` changes those probabilities to `{cert['counterexample']['own_scale_reweighted_probabilities']}`.  That changes the apparent cause by changing the observer's scale unit, so it is forbidden for causal weighting.

After the common-unit law exists, its shell pushforward may be used diagnostically.  Since every high shell has `M_j/N>=2`, the maximal binary-work shell satisfies

`ownscale_selected * exp(H_inf^binary-scale) >= 2 P_binary`,

and hence at least `nu D_tail` on an HH owner.  This is a strength-versus-scale-concentration relation, not a new probability law.

No claim is made that generic high-tail HH parents are comparable to the child scale.  In particular, nearly cancelling parents may live at `K>>M`; that nonlocal high-high-to-low geometry is the next physical seam.  Likewise no Young near-extremality or `W_HH>=8E1/15` child-energy productivity gate is assumed.  KL/log-productivity may be invoked only if that independent generated-energy hypothesis is later supplied.

Stress: `{out.samples}` common-unit owner / coherent Hahn / probability states
- minimum common HH/interface owner margin: `{out.minimum_common_owner_margin:.3e}`
- minimum clean `nu D_tail/2` owner margin: `{out.minimum_clean_common_owner_margin:.3e}`
- minimum atomic positive-dominance margin: `{out.minimum_atomic_positive_dominance_margin:.3e}`
- worst binary probability residual: `{out.worst_binary_probability_residual:.3e}`
- minimum HH-owner binary clean margin: `{out.minimum_HH_binary_clean_margin:.3e}`
- minimum `H_inf` own-scale diagnostic margin: `{out.minimum_Hinf_ownscale_diagnostic_margin:.3e}`
- worst common-unit probability-invariance residual: `{out.worst_common_scale_probability_invariance_residual:.3e}`
- maximum joint owner count: `{out.maximum_joint_owner_count}`
- forbidden scale-reweighting example probability distortion: `{out.scale_reweighting_counterexample_distortion:.6f}`

This theorem adds no packet, no reset currency, and no locality hypothesis.  It only keeps the physical work law in its natural unit and exposes the binary Navier--Stokes cause already present in the quadratic term.  No global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
