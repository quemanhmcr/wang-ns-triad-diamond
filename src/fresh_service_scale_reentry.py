from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from src.critical_shell_service_reentry import (
    critical_shell_bounded_service_lower,
    critical_shell_integrated_service_lower,
)

STATUS = (
    "EXACT_FRESH_MATERIAL_SERVICE_TO_REFINEMENT_INVARIANT_SCALE_LAW__"
    "CANONICAL_ANNULAR_LP_TO_TWO_HARD_SHELLS__NO_COHERENT_CELL_DOMINANCE_REQUIRED"
)

ANNULAR_SUPPORT_LOWER = 0.5
ANNULAR_SUPPORT_UPPER = 2.0
FRESH_SERVICE_FRACTION = 0.25
INCREMENT_TO_BAND_ENERGY = 4.0
TWO_HARD_SHELL_FACTOR = 2.0 / 3.0


def canonical_annular_frame_registration() -> dict[str, object]:
    """Register one global smooth square-normalized dyadic analysis frame.

    Choose a real nonnegative bump q supported in (1/2,2) and positive on
    [3/4,3/2].  For M_j=2^j N set q_j(xi)=q(|xi|/M_j) and

        phi_j=q_j/(sum_k q_k^2)^(1/2).

    The dyadic positive cores cover every nonzero frequency, so the denominator
    is nonzero there.  Normalization preserves support, gives sum phi_j^2=1 and
    |phi_j|<=1.  Thus u_j=phi_j(D)u and u=sum phi_j(D)u_j exactly in L2.
    """
    return {
        "support_lower_ratio": ANNULAR_SUPPORT_LOWER,
        "support_upper_ratio": ANNULAR_SUPPORT_UPPER,
        "square_partition": "sum_j phi_j(xi)^2=1 for xi!=0",
        "calderon_reconstruction": "u=sum_j phi_j(D)[phi_j(D)u]",
        "pointwise_multiplier_bound": "|phi_j|<=1",
        "construction": "square-normalize one smooth dyadic annular cover q_j with supp q_j subset (M_j/2,2M_j)",
    }


def pushforward_fresh_edges_to_bands(
    edge_weights: Sequence[float],
    band_indices: Sequence[int],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
) -> dict[int, float]:
    """Push the positive NN/fresh service measure to the canonical band index.

    Ownership is pointwise on the positive service law: an edge is fresh exactly
    when both intrinsic endpoints lie outside the transported old material set.
    Summing by the deterministic LP band index therefore quotients every coherent
    cell refinement.  Cell names never enter the output.
    """
    w = np.asarray(edge_weights, dtype=float)
    j = np.asarray(band_indices, dtype=int)
    a = np.asarray(old_here, dtype=bool)
    b = np.asarray(old_neighbor, dtype=bool)
    if w.ndim != 1 or j.shape != w.shape or a.shape != w.shape or b.shape != w.shape:
        raise ValueError("matching one-dimensional fresh-edge data required")
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("finite nonnegative edge weights required")
    if np.any(j > 0):
        raise ValueError("fresh low/base service theorem expects canonical bands j<=0")
    out: dict[int, float] = {}
    fresh = (~a) & (~b)
    for weight, band, keep in zip(w, j, fresh):
        if keep and weight > 0:
            key = int(band)
            out[key] = out.get(key, 0.0) + float(weight)
    return dict(sorted(out.items()))


def scale_law_statistics(fresh_band_weights: Mapping[int, float]) -> dict[str, float | int]:
    """Concentration data of the refinement-invariant fresh scale law."""
    if not fresh_band_weights:
        raise ValueError("positive fresh band law required")
    items = sorted((int(k), float(v)) for k, v in fresh_band_weights.items())
    if any(k > 0 or not math.isfinite(v) or v < 0 for k, v in items):
        raise ValueError("finite nonnegative low/base band weights required")
    total = sum(v for _, v in items)
    if total <= 0 or not math.isfinite(total):
        raise ValueError("positive finite fresh band service required")
    probs = [(k, v / total) for k, v in items if v > 0]
    if not probs:
        raise AssertionError("positive total fresh service lost all atoms")
    # A countable positive probability law has an attained maximum: otherwise
    # infinitely many atoms would exceed half a positive supremum.
    jmax, pmax = max(probs, key=lambda kv: kv[1])
    q2 = sum(p * p for _, p in probs)
    h_inf = -math.log(pmax)
    h2 = -math.log(q2)
    if pmax + 2e-15 < q2:
        raise AssertionError("scale-law max atom failed to dominate collision mass")
    return {
        "fresh_service": total,
        "selected_band": int(jmax),
        "p_max": pmax,
        "H_inf_scale": h_inf,
        "H2_scale": h2,
        "effective_scale_count_inf": 1.0 / pmax,
        "effective_scale_count_2": 1.0 / q2,
    }


def smooth_band_to_hard_shell_mass_lower(smooth_band_critical_mass: float) -> float:
    """One annular analysis band forces one of two exact hard-shell masses.

    With supp phi_j subset {M/2<|xi|<2M} and |phi_j|<=1, split into
      A_0={M/2<|xi|<=M}, A_1={M<|xi|<=2M}.
    If mu_0=M||P_A0 u||_2^2 and mu_1=2M||P_A1 u||_2^2, then

      M||phi_j(D)u||_2^2 <= mu_0 + mu_1/2
                              <= (3/2) max(mu_0,mu_1).

    Hence one actual hard shell has mass at least 2/3 of the smooth-band mass.
    """
    b = float(smooth_band_critical_mass)
    if b < 0 or not math.isfinite(b):
        raise ValueError("finite nonnegative smooth-band mass required")
    return TWO_HARD_SHELL_FACTOR * b


def selected_band_peak_smooth_mass_lower(fresh_band_service: float, scaled_lifetime: float) -> float:
    """Integrated fresh service -> some-time smooth-band critical mass.

    Pointwise Moyal plus ||delta_r f||_2<=2||f||_2 gives
      F_j <= 4 int M_j||u_j||_2^2 d tau.
    On a scaled interval of length c, some time therefore has
      M_j||u_j||_2^2 >= F_j/(4c).
    """
    f = float(fresh_band_service)
    c = float(scaled_lifetime)
    if f < 0 or c <= 0 or not all(math.isfinite(x) for x in (f, c)):
        raise ValueError("finite nonnegative service and positive lifetime required")
    return f / (INCREMENT_TO_BAND_ENERGY * c)


def selected_band_hard_shell_mass_lower(fresh_band_service: float, scaled_lifetime: float) -> float:
    return smooth_band_to_hard_shell_mass_lower(
        selected_band_peak_smooth_mass_lower(fresh_band_service, scaled_lifetime)
    )


def fresh_service_scale_route(
    integrated_square_service_threshold: float,
    scaled_lifetime: float,
    block_frequency: float,
    fresh_band_weights: Mapping[int, float],
    *,
    viscosity: float = 1.0,
) -> dict[str, object]:
    """Every sufficiently large fresh NN service law supplies a hard-shell seed.

    Existing aggregate coherent-service routing supplies integrated fresh NN service F>=Y/4 over the parent scaled interval of length c.
    Push it to the fixed LP band index.  If p_max is the largest scale mass, then
    that band carries at least p_max F.  Capacity, time averaging, and the two-hard
    shell cover give

      mu_hard >= p_max F/(6c) >= p_max Y/(24c).

    Equivalently mu_hard exp(H_inf^scale)>=Y/(24c).  Since H2>=H_inf, the same
    lower also holds after multiplying by exp(H2).  This is a scale-concentration
    relation, not a causal probability and not a new stop class.
    """
    Y = float(integrated_square_service_threshold)
    c = float(scaled_lifetime)
    N = float(block_frequency)
    nu = float(viscosity)
    if Y <= 0 or c <= 0 or N <= 0 or nu < 0 or not all(math.isfinite(x) for x in (Y, c, N, nu)):
        raise ValueError("positive finite service/lifetime/frequency and nonnegative viscosity required")
    stats = scale_law_statistics(fresh_band_weights)
    fresh = float(stats["fresh_service"])
    cover_tol = 4e-13 * max(1.0, Y, fresh)
    if fresh + cover_tol < FRESH_SERVICE_FRACTION * Y:
        raise ValueError("fresh band law does not realize the certified Y/4 service lower")
    j = int(stats["selected_band"])
    pmax = float(stats["p_max"])
    h_inf = float(stats["H_inf_scale"])
    h2 = float(stats["H2_scale"])
    selected_service = float(fresh_band_weights[j])
    actual_mu = selected_band_hard_shell_mass_lower(selected_service, c)
    clean_mu = pmax * Y / (24.0 * c)
    mass_tol = 4e-13 * max(1.0, actual_mu, clean_mu)
    if actual_mu + mass_tol < clean_mu:
        raise AssertionError("fresh scale law lost the clean hard-shell lower")
    M = N * (2.0 ** j)
    if M <= 0 or M > N + 2e-13 * max(1.0, N):
        raise AssertionError("selected low/base LP band escaped its canonical scale range")

    # Generic shell theorem is material-free at entrance.  Fresh NN provenance is
    # retained only as a sidecar and is reread from renewed service afterwards.
    y_shell = critical_shell_bounded_service_lower(actual_mu, c, nu)
    s_shell = critical_shell_integrated_service_lower(actual_mu, c, nu)
    # Use the effective counts directly rather than exponentiating large logs.
    exp_hinf = float(stats["effective_scale_count_inf"])
    exp_h2 = float(stats["effective_scale_count_2"])
    clean_unweighted_mass = Y / (24.0 * c)
    clean_weighted_y = critical_shell_bounded_service_lower(clean_unweighted_mass, c, nu)
    clean_weighted_s = critical_shell_integrated_service_lower(clean_unweighted_mass, c, nu)
    weighted_y_inf = y_shell * exp_hinf
    weighted_s_inf = s_shell * exp_hinf
    service_tol = 2e-12 * max(1.0, clean_weighted_y, clean_weighted_s, weighted_y_inf, weighted_s_inf)
    if actual_mu * exp_hinf + service_tol < clean_unweighted_mass:
        raise AssertionError("fresh H-infinity shell tradeoff failed")
    if actual_mu * exp_h2 + service_tol < clean_unweighted_mass:
        raise AssertionError("fresh H2 shell corollary failed")
    if weighted_y_inf + service_tol < clean_weighted_y or weighted_s_inf + service_tol < clean_weighted_s:
        raise AssertionError("fresh service-complexity conjugacy failed")

    return {
        "fresh_service": fresh,
        "selected_band": j,
        "selected_band_frequency": M,
        "selected_band_service": selected_service,
        "p_max": pmax,
        "H_inf_scale": h_inf,
        "H2_scale": h2,
        "smooth_band_peak_mass_lower": selected_band_peak_smooth_mass_lower(selected_service, c),
        "hard_shell_mass_lower": actual_mu,
        "clean_hard_shell_mass_lower": clean_mu,
        "H_inf_weighted_hard_shell_mass_lower": actual_mu * exp_hinf,
        "H2_weighted_hard_shell_mass_lower": actual_mu * exp_h2,
        "clean_weighted_hard_shell_mass_lower": clean_unweighted_mass,
        "hard_shell_candidates": (M, 2.0 * M),
        "maximum_candidate_shell_over_block_scale": 2.0 * M / N,
        "full_survivor_own_scale_service_lower": y_shell,
        "full_survivor_integrated_service_lower": s_shell,
        "H_inf_weighted_full_survivor_service_lower": weighted_y_inf,
        "H_inf_weighted_full_survivor_integrated_service_lower": weighted_s_inf,
        "clean_H_inf_weighted_full_survivor_service_lower": clean_weighted_y,
        "clean_H_inf_weighted_full_survivor_integrated_service_lower": clean_weighted_s,
        "material_semantics": "fresh NN service selects only scale provenance; the whole hard u-shell is not declared fresh material",
        "probability_semantics": "the normalized band law is a deterministic service-scale diagnostic, not a child-energy causal probability",
        "next_owner": "generic_critical_shell_first_stop",
        "master_semantics": "RECURSE_CRITICAL_VIA_GENERIC_SHELL",
        "status": STATUS,
    }


def theorem_certificate() -> dict[str, object]:
    frame = canonical_annular_frame_registration()
    return {
        "status": STATUS,
        "material_quotient": "fresh NN is the pointwise indicator 1_Oc(z0)1_Oc(z1) on the positive coherent edge measure; pushforward to band j is invariant under coherent-cell refinement",
        "canonical_frame": frame,
        "band_capacity": "F_j<=M_j int||delta_r u_j||_2^2 d tau<=4 int M_j||u_j||_2^2 d tau",
        "two_hard_shell_cover": "supp phi_j subset (M_j/2,2M_j), |phi_j|<=1 => M_j||u_j||_2^2<=mu_M+(1/2)mu_2M<=(3/2)max(mu_M,mu_2M)",
        "native_tradeoff": "fresh F>=Y/4 => mu_hard>=pmax Y/(24c), equivalently mu_hard exp(H_inf_scale)>=Y/(24c)",
        "renyi2_corollary": "because pmax>=sum p_j^2, mu_hard exp(H2_scale)>=Y/(24c)",
        "hard_shell_selection": "at the witness time choose the larger of the two exact hard-shell masses at M_j and 2M_j; this finite measurable selection is physical scale selection, not a coherent-cell argmax",
        "generic_shell": "every fresh scale law enters the existing material-free critical-shell first-stop theorem; full-survivor service inherits the same H_inf weighting by linearity",
        "no_cell_dominance": "coherent-cell quarter dominance/entropy/cycle may remain useful fine ancestry accounting, but it is not required for fresh-service renewal entrance",
        "material_scope": "fresh edge provenance is not promoted to whole-shell freshness; material OO/ON/NN is reread only from subsequent actual renewed service",
        "scale_scope": "for j<=0 the two candidate hard-shell frequencies are M_j and 2M_j<=2N; supplier-specific signed-good scale progress is not asserted",
        "master_rule": "H_inf/H2 here are deterministic scale-concentration coordinates, not causal Shannon/Renyi probabilities and not additive resets",
    }


@dataclass(frozen=True)
class FreshScaleStress:
    samples: int
    worst_refinement_pushforward_residual: float
    minimum_two_hard_shell_margin: float
    minimum_clean_hard_shell_margin: float
    minimum_Hinf_tradeoff_margin: float
    minimum_H2_tradeoff_margin: float
    minimum_full_survivor_service_conjugacy_margin: float
    maximum_candidate_shell_over_block_scale: float


def stress(samples: int = 50_000, seed: int = 20260810) -> FreshScaleStress:
    rng = np.random.default_rng(seed)
    wr = 0.0
    mh = mc = mi = m2 = ms = float("inf")
    max_ratio = 0.0

    for _ in range(samples):
        # Cell refinement cannot alter the fresh pushforward to the band index.
        n = int(rng.integers(1, 24))
        w = rng.lognormal(mean=-3.0, sigma=1.5, size=n)
        bands = rng.integers(-12, 1, size=n)
        oh = rng.integers(0, 2, size=n).astype(bool)
        on = rng.integers(0, 2, size=n).astype(bool)
        coarse = pushforward_fresh_edges_to_bands(w, bands, oh, on)
        # Split every physical atom into two representation records with identical
        # ownership and band marks.  Pushforward must be exactly unchanged.
        lam = rng.random(n)
        w2 = np.concatenate([lam * w, (1.0 - lam) * w])
        b2 = np.concatenate([bands, bands])
        oh2 = np.concatenate([oh, oh])
        on2 = np.concatenate([on, on])
        fine = pushforward_fresh_edges_to_bands(w2, b2, oh2, on2)
        keys = set(coarse) | set(fine)
        resid = max((abs(coarse.get(k, 0.0) - fine.get(k, 0.0)) for k in keys), default=0.0)
        wr = max(wr, resid)
        if resid > 2e-12 * max(1.0, float(w.sum())):
            raise AssertionError("fresh band pushforward changed under coherent-cell refinement")

        # Exact two-hard-shell cover algebra.
        B = float(rng.lognormal(mean=-2.0, sigma=1.5))
        x = float(rng.random())
        # Construct nonnegative hard masses with mu0+mu1/2 >= B.
        mu0 = x * B
        mu1 = 2.0 * (1.0 - x) * B * float(rng.uniform(1.0, 2.0))
        cover = mu0 + 0.5 * mu1
        if cover < B:
            mu1 += 2.0 * (B - cover)
        actual = max(mu0, mu1)
        lower = smooth_band_to_hard_shell_mass_lower(B)
        mh = min(mh, actual - lower)
        if actual + 2e-12 * max(1.0, B) < lower:
            raise AssertionError("two-hard-shell cover lower failed")

        # Fresh service route.  Build a positive low-band law with total >=Y/4.
        Y = float(rng.lognormal(mean=-1.5, sigma=1.4))
        c = float(rng.lognormal(mean=-0.2, sigma=1.0))
        N = float(rng.lognormal(mean=1.0, sigma=1.2))
        nu = float(rng.uniform(0.0, 2.0))
        k = int(rng.integers(1, 12))
        probs = rng.dirichlet(np.ones(k))
        fresh = FRESH_SERVICE_FRACTION * Y * float(rng.uniform(1.0, 3.0))
        labels = list(range(-k + 1, 1))
        law = {j: float(p * fresh) for j, p in zip(labels, probs)}
        route = fresh_service_scale_route(Y, c, N, law, viscosity=nu)
        pmax = float(route["p_max"])
        clean = pmax * Y / (24.0 * c)
        mu = float(route["hard_shell_mass_lower"])
        mc = min(mc, mu - clean)
        if mu + 3e-12 * max(1.0, clean) < clean:
            raise AssertionError("fresh clean hard-shell lower failed")
        base = Y / (24.0 * c)
        hinf_mass = float(route["H_inf_weighted_hard_shell_mass_lower"])
        h2_mass = float(route["H2_weighted_hard_shell_mass_lower"])
        mi = min(mi, hinf_mass - base)
        m2 = min(m2, h2_mass - base)
        if hinf_mass + 3e-12 * max(1.0, base) < base or h2_mass + 3e-12 * max(1.0, base) < base:
            raise AssertionError("fresh scale entropy tradeoff failed")
        weighted_s = float(route["H_inf_weighted_full_survivor_integrated_service_lower"])
        clean_s = float(route["clean_H_inf_weighted_full_survivor_integrated_service_lower"])
        ms = min(ms, weighted_s - clean_s)
        if weighted_s + 3e-12 * max(1.0, clean_s) < clean_s:
            raise AssertionError("fresh full-survivor service conjugacy failed")
        max_ratio = max(max_ratio, float(route["maximum_candidate_shell_over_block_scale"]))
        if route["maximum_candidate_shell_over_block_scale"] > 2.0 + 2e-12:
            raise AssertionError("fresh hard-shell candidate escaped 2N")

    return FreshScaleStress(samples, wr, mh, mc, mi, m2, ms, max_ratio)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-fresh-service-scale-reentry"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "fresh_service_scale_reentry.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Fresh coherent service -> refinement-invariant scale law -> hard critical shell

Status: **{cert['status']}**.

The fresh/new-new branch of coherent increment service should not depend on how an observer subdivides phase space into coherent cells.  Let `O` be the already transported measurable old material set.  On the positive Moyal increment-edge measure classify an edge as fresh by the pointwise indicator

`1_(O^c)(zeta_0) 1_(O^c)(zeta_1)`.

Push this positive fresh measure only to the fixed Littlewood--Paley band index `j`.  Its weights `F_j` are unchanged by every coherent-cell refinement: splitting one physical edge into many representation records changes no band total.

Fix the canonical smooth square-normalized annular frame once and for all from a dyadic bump cover with

`supp phi_j subset {{M_j/2<|xi|<2M_j}}`, `sum_j phi_j^2=1`, `|phi_j|<=1`,

and exact Calderon reconstruction `u=sum_j phi_j(D)[phi_j(D)u]`.  This merely fixes the standard annular representative of the smooth dyadic-cover construction; it introduces no packet.  The upstream theorems use the finite square-function/Bernstein constants `C_LP,C_B` of whichever canonical frame is fixed, so they are now understood to be the constants of this representative; no scale-dependent price is introduced.

For a fresh band,

`F_j <= M_j int ||delta_r u_j||_2^2 d tau <= 4 int M_j||u_j||_2^2 d tau`.

Over a scaled lifetime `c`, some time therefore has

`M_j||u_j||_2^2 >= F_j/(4c)`.

Because the smooth band touches only the two hard annuli

`A_0={{M_j/2<|xi|<=M_j}}`, `A_1={{M_j<|xi|<=2M_j}}`,

and `|phi_j|<=1`,

`M_j||u_j||_2^2 <= mu_(M_j) + (1/2)mu_(2M_j) <= (3/2) max(mu_(M_j),mu_(2M_j))`.

Hence one **actual hard shell of u** carries

`mu_hard >= F_j/(6c)`.

If total **integrated** fresh service on a parent scaled interval of length `c` satisfies the already-certified `F>=Y/4`, normalize `p_j=F_j/F`.  Let `p_max` be the largest scale atom and `H_inf^scale=-log p_max`.  Then

`mu_hard >= p_max Y/(24c)`,

or equivalently

`mu_hard exp(H_inf^scale) >= Y/(24c)`.

Since `p_max>=sum_j p_j^2`, the weaker collision-entropy corollary also holds:

`mu_hard exp(H_2^scale) >= Y/(24c)`.

There is no quarter coherent-cell dominance in this renewal entrance.  `H_inf^scale` and `H_2^scale` are deterministic concentration coordinates of the **canonical frequency pushforward**, not child-energy causal probabilities and not new stop classes.  The selected fresh NN edge law remains material provenance only; the whole hard `u` shell is not declared new material.  It enters the existing material-free critical-shell first-stop theorem, and material OO/ON/NN is reread only from subsequent actual renewed service.

On a full no-hit natural shell corridor, linearity of the generic shell service in `mu` preserves the same concentration tradeoff for both instantaneous own-scale service and integrated service.

Stress: `{out.samples}` refinement/scale/shell/service states
- worst refinement-pushforward residual: `{out.worst_refinement_pushforward_residual:.3e}`
- minimum two-hard-shell margin: `{out.minimum_two_hard_shell_margin:.3e}`
- minimum clean hard-shell margin: `{out.minimum_clean_hard_shell_margin:.3e}`
- minimum H_inf tradeoff margin: `{out.minimum_Hinf_tradeoff_margin:.3e}`
- minimum H2 tradeoff margin: `{out.minimum_H2_tradeoff_margin:.3e}`
- minimum full-survivor service-conjugacy margin: `{out.minimum_full_survivor_service_conjugacy_margin:.3e}`
- maximum candidate shell / block scale: `{out.maximum_candidate_shell_over_block_scale:.9f}`

This theorem is local renewal geometry.  It does not prove supplier-specific signed-good scale progress, does not turn scale concentration into a finite reset, and does not convert fresh service into near-extremal HH transfer.  No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
