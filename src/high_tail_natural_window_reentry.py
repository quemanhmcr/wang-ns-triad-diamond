from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.critical_shell_service_reentry import (
    critical_shell_bounded_service_lower,
    critical_shell_integrated_service_lower,
)
from src.high_tail_ultraviolet_locality import (
    STATUS as ULTRAVIOLET_LOCALITY_STATUS,
    high_tail_hh_locality_tradeoff,
    ultraviolet_hh_work_constant,
)


STATUS = (
    "EXACT_HIGH_TAIL_COMPARABLE_HH_TO_SLIDING_NATURAL_WINDOW__"
    "TIME_ORIGIN_AND_UNIT_INVARIANT__SCALE_TIME_CONCENTRATION_TO_CRITICAL_SHELL__"
    "NO_PACKET_PERSISTENCE_OR_TIME_BINNING"
)

COMPARABLE_OWNER = "comparable_parent_HH_work"


def natural_window_geometry(
    parent_frequency: float,
    selected_shell_level: int,
    scaled_lifetime: float,
) -> dict[str, float | int]:
    """Parent block and selected hard-shell natural-time geometry.

    The high-tail shell is M=2^j N with j>=1.  Therefore
      T_N=cN^-2,
      T_M=cM^-2=4^-j T_N <= T_N/4.
    This is exact support geometry, not signed-good Young progress.
    """
    N = float(parent_frequency)
    j = int(selected_shell_level)
    c = float(scaled_lifetime)
    if N <= 0 or c <= 0 or j < 1 or not all(math.isfinite(x) for x in (N, c)):
        raise ValueError("positive finite parent frequency/lifetime and high-tail level j>=1 required")
    M = N * (2.0**j)
    T_parent = c / (N * N)
    T_child = c / (M * M)
    ratio = T_child / T_parent
    return {
        "parent_frequency": N,
        "selected_shell_level": j,
        "selected_shell_frequency": M,
        "forward_scale_ratio": M / N,
        "parent_natural_duration": T_parent,
        "selected_natural_window": T_child,
        "natural_time_ratio": ratio,
    }


def _piecewise_integral(
    starts: np.ndarray,
    ends: np.ndarray,
    densities: np.ndarray,
    left: float,
    right: float,
) -> float:
    overlap = np.maximum(0.0, np.minimum(ends, right) - np.maximum(starts, left))
    return float(np.dot(densities, overlap))


def sliding_window_piecewise_constant(
    segment_starts: Sequence[float],
    segment_ends: Sequence[float],
    densities: Sequence[float],
    block_start: float,
    block_end: float,
    window_length: float,
) -> dict[str, object]:
    """Exact finite representation of a sliding positive-work measure.

    For a piecewise-constant nonnegative density rho, F(s)=int_s^(s+T)rho is
    piecewise linear.  Its breakpoints occur when s or s+T crosses a segment
    boundary, so the maximum is attained among block endpoints, segment
    boundaries, and those boundaries shifted by -T.  This is a stress/helper
    representation only; the continuum theorem uses the sliding measure itself.
    """
    a = np.asarray(segment_starts, float)
    b = np.asarray(segment_ends, float)
    d = np.asarray(densities, float)
    L = float(block_start)
    U = float(block_end)
    T = float(window_length)
    if a.ndim != 1 or b.shape != a.shape or d.shape != a.shape or len(a) == 0:
        raise ValueError("matching nonempty one-dimensional piecewise density data required")
    if not all(math.isfinite(x) for x in (L, U, T)) or U <= L or T <= 0 or T > U - L:
        raise ValueError("finite block and sliding window with 0<T<=block length required")
    if np.any(~np.isfinite(a)) or np.any(~np.isfinite(b)) or np.any(~np.isfinite(d)):
        raise ValueError("finite segment data required")
    if np.any(d < 0) or np.any(b <= a):
        raise ValueError("nonnegative densities on positive-length segments required")
    order = np.argsort(a)
    a, b, d = a[order], b[order], d[order]
    tol_t = 2e-13 * max(1e-300, abs(L), abs(U), abs(T), U - L)
    if a[0] < L - tol_t or b[-1] > U + tol_t:
        raise ValueError("piecewise density must lie inside the observed block")
    if len(a) > 1 and np.any(a[1:] < b[:-1] - tol_t):
        raise ValueError("piecewise density segments must be disjoint")

    total = float(np.dot(d, b - a))
    if total <= 0 or not math.isfinite(total):
        raise ValueError("positive finite comparable-work measure required")

    lo = L
    hi = U - T
    boundaries = np.concatenate((a, b))
    candidates = [lo, hi]
    for x in boundaries:
        for s in (float(x), float(x) - T):
            if lo - tol_t <= s <= hi + tol_t:
                candidates.append(min(hi, max(lo, s)))
    candidates = sorted(set(candidates))
    values = [_piecewise_integral(a, b, d, s, s + T) for s in candidates]
    imax = int(np.argmax(values))
    s_star = float(candidates[imax])
    w_star = float(values[imax])
    p_t = w_star / total
    if not (0.0 < p_t <= 1.0 + 3e-13):
        raise AssertionError("sliding positive-work concentration failed to define a probability fraction")
    p_t = min(1.0, p_t)
    return {
        "total_comparable_common_work": total,
        "maximum_window_common_work": w_star,
        "selected_window_start": s_star,
        "selected_window_end": s_star + T,
        "window_length": T,
        "p_time": p_t,
        "H_inf_time": -math.log(p_t),
        "candidate_count": len(candidates),
    }


def temporal_concentration_statistics(
    total_comparable_common_work: float,
    maximum_window_common_work: float,
    window_length: float,
) -> dict[str, float]:
    """Measure-native sliding concentration with its physical window attached."""
    W = float(total_comparable_common_work)
    Ww = float(maximum_window_common_work)
    T = float(window_length)
    if W <= 0 or Ww <= 0 or T <= 0 or not all(math.isfinite(x) for x in (W, Ww, T)):
        raise ValueError("positive finite total/window work and window length required")
    tol = 4e-13 * max(1e-300, W, Ww)
    if Ww > W + tol:
        raise ValueError("a sliding-window submeasure cannot exceed total comparable work")
    p = min(1.0, Ww / W)
    return {
        "p_time": p,
        "H_inf_time": -math.log(p),
        "effective_time_window_count": 1.0 / p,
        "window_length": T,
    }


def comparable_natural_window_common_work_upper(
    window_peak_child_mass: float,
    parent_frequency: float,
    global_energy: float,
    scaled_lifetime: float,
    locality_radius: float,
) -> float:
    """Common-unit positive comparable HH work capacity on one M-natural window.

    The canonical strict low-pass only supplies |S|<=1, hence
      ||h||_2=||(I-S)u||_2 <= 2||u||_2.
    For child shell M and both comparable parents <=R M,
      a_c <= C_hp sqrt(mu),
      a_1 a_2 <= 4 C_hp^2 R M E_global.
    Physical work rate therefore obeys
      r_comp <= 12 sqrt(pi) R M^2 E_global sqrt(mu).
    Integrating cM^-2 cancels M^2.  The high-tail causal law uses common unit
    N dW, so multiply by parent block N:
      N W_window <= 12 c sqrt(pi) R N E_global sqrt(mu_window).
    """
    mu = float(window_peak_child_mass)
    N = float(parent_frequency)
    E = float(global_energy)
    c = float(scaled_lifetime)
    R = float(locality_radius)
    if mu < 0 or min(N, E, c) <= 0 or R <= 1.0 or not all(math.isfinite(x) for x in (mu, N, E, c, R)):
        raise ValueError("finite nonnegative shell mass and positive frequency/energy/lifetime/R>1 required")
    # ultraviolet_hh_work_constant = 3 sqrt(pi).  The factor 4R is exactly the
    # two unresolved-parent L2 bounds (2 each) and their Fourier ball radius R.
    return 4.0 * R * ultraviolet_hh_work_constant() * c * N * E * math.sqrt(mu)


def comparable_hh_temporal_shell_reentry(
    locality_route: dict[str, object],
    physical_tail_dissipation: float,
    viscosity: float,
    parent_frequency: float,
    global_energy: float,
    scaled_lifetime: float,
    total_comparable_common_work: float,
    maximum_window_common_work: float,
    window_length: float,
    window_peak_child_mass: float,
) -> dict[str, object]:
    """Compose certified Fourier locality with sliding natural-time concentration.

    The locality route must already carry the clean comparable-HH owner.  Let
      p_s = selected output-shell fraction,
      p_t = maximal sliding M-natural-window fraction of actual comparable work.
    Since W_comp/p_s >= nu D_tail/4 and W_win=p_t W_comp, while the natural-window
    capacity is C_t sqrt(mu_win),

      sqrt(mu_win)/(p_s p_t)
        >= nu D_tail / [48 c sqrt(pi) R N E_global].

    No packet persistence or observer-chosen time partition occurs.
    """
    if locality_route.get("status") != ULTRAVIOLET_LOCALITY_STATUS:
        raise ValueError("certified ultraviolet-locality route required")
    owners = tuple(locality_route.get("joint_clean_owners", ()))
    if COMPARABLE_OWNER not in owners:
        raise ValueError("locality route does not carry the clean comparable-HH owner")

    D = float(physical_tail_dissipation)
    nu = float(viscosity)
    N = float(parent_frequency)
    E = float(global_energy)
    c = float(scaled_lifetime)
    W = float(total_comparable_common_work)
    Ww = float(maximum_window_common_work)
    Tw = float(window_length)
    mu = float(window_peak_child_mass)
    if min(D, nu, N, E, c, W, Ww, Tw) <= 0 or mu <= 0 or not all(math.isfinite(x) for x in (D, nu, N, E, c, W, Ww, Tw, mu)):
        raise ValueError("positive finite dissipation/viscosity/frequency/energy/lifetime/work/window/mass required")

    p_s = float(locality_route["p_max"])
    j = int(locality_route["selected_shell_level"])
    R = float(locality_route["locality_radius"])
    mu_block = float(locality_route["child_peak_critical_mass"])
    W_lower = float(locality_route["comparable_parent_common_work_lower"])
    tol_work = 6e-13 * max(1e-300, W, W_lower)
    if W + tol_work < W_lower:
        raise ValueError("actual comparable work is below the certified locality lower")
    if mu > mu_block + 6e-13 * max(1e-300, mu, mu_block):
        raise ValueError("window shell peak cannot exceed the certified block shell peak")

    geometry = natural_window_geometry(N, j, c)
    Tnatural = float(geometry["selected_natural_window"])
    tol_time = 8e-13 * max(1e-300, abs(Tw), abs(Tnatural))
    if abs(Tw - Tnatural) > tol_time:
        raise ValueError("temporal concentration window is not the selected shell natural window c M^-2")
    temporal = temporal_concentration_statistics(W, Ww, Tw)
    p_t = float(temporal["p_time"])
    capacity = comparable_natural_window_common_work_upper(mu, N, E, c, R)
    if Ww > capacity + 8e-13 * max(1e-300, Ww, capacity):
        raise ValueError("sliding natural-window comparable work exceeds the physical energy capacity")

    scale_weighted_actual = W / p_s
    clean_scale_weighted_work = 0.25 * nu * D
    if scale_weighted_actual + 8e-13 * max(1e-300, scale_weighted_actual, clean_scale_weighted_work) < clean_scale_weighted_work:
        raise AssertionError("actual comparable work lost the certified locality clean lower")

    weighted_sqrt_mass = math.sqrt(mu) / (p_s * p_t)
    clean_weighted_sqrt_mass = nu * D / (
        48.0 * c * math.sqrt(math.pi) * R * N * E
    )
    margin = weighted_sqrt_mass - clean_weighted_sqrt_mass
    if margin < -1e-11 * max(1e-300, weighted_sqrt_mass, clean_weighted_sqrt_mass):
        raise AssertionError("scale-time concentration to critical-shell mass failed")

    weighted_mass = mu / ((p_s * p_t) ** 2)
    clean_weighted_mass = clean_weighted_sqrt_mass**2
    Hs = -math.log(p_s)
    Ht = float(temporal["H_inf_time"])

    # The actual window peak is an actual hard-shell event, so the generic shell
    # theorem may consume mu directly.  Service is conditional on its full no-hit
    # natural corridor, exactly as in the generic theorem.
    full_survivor_service = critical_shell_bounded_service_lower(mu, c, nu)
    full_survivor_integrated_service = critical_shell_integrated_service_lower(mu, c, nu)

    return {
        "selected_shell_level": j,
        "selected_shell_frequency": float(geometry["selected_shell_frequency"]),
        "forward_scale_ratio": float(geometry["forward_scale_ratio"]),
        "natural_time_ratio": float(geometry["natural_time_ratio"]),
        "selected_natural_window": Tnatural,
        "p_scale": p_s,
        "H_inf_output_scale": Hs,
        "p_time": p_t,
        "H_inf_time": Ht,
        "scale_time_concentration_product": p_s * p_t,
        "total_comparable_common_work": W,
        "maximum_window_common_work": Ww,
        "window_peak_child_mass": mu,
        "natural_window_common_work_capacity": capacity,
        "entropy_weighted_sqrt_child_mass": weighted_sqrt_mass,
        "clean_entropy_weighted_sqrt_child_mass_lower": clean_weighted_sqrt_mass,
        "entropy_weighted_child_mass": weighted_mass,
        "clean_entropy_weighted_child_mass_lower": clean_weighted_mass,
        "scale_time_tradeoff_margin": margin,
        "scale_critical_global_energy": N * E,
        "next_owner": "generic_critical_shell_first_stop",
        "full_survivor_own_scale_service_lower": full_survivor_service,
        "full_survivor_integrated_service_lower": full_survivor_integrated_service,
        "full_survivor_service_is_conditional": True,
        "time_partition_used": False,
        "packet_persistence_used": False,
        "signed_good_progress_used": False,
        "master_semantics": "RECURSE_CRITICAL_VIA_ACTUAL_HIGH_TAIL_SHELL_EVENT",
        "status": STATUS,
    }


def theorem_certificate() -> dict[str, object]:
    C = 4.0 * 2.0 * ultraviolet_hh_work_constant()
    if abs(C - 24.0 * math.sqrt(math.pi)) > 3e-14:
        raise AssertionError("dyadic natural-window coefficient did not equal 24 sqrt(pi)")
    return {
        "status": STATUS,
        "upstream": ULTRAVIOLET_LOCALITY_STATUS,
        "sliding_measure": "p_t=sup_s mu_comp([s,s+cM^-2])/mu_comp(I); no time bins or chosen bin origin",
        "attainment": "for smooth NS the positive comparable-work density is continuous, hence the sliding integral is continuous on the compact admissible start interval and attains its maximum",
        "time_gauge": "translation of time origin and simultaneous rescaling t->lambda t, density->density/lambda, window->lambda window leave p_t unchanged",
        "refinement_rule": "subdividing a piecewise representation without changing the positive measure leaves the sliding maximum and p_t unchanged",
        "cutoff_scope": "only |S|<=1 is used, so ||(I-S)u||_2<=2||u||_2; no unregistered nonnegative-cutoff contraction is assumed",
        "window_capacity": "N W_win <= 12 c sqrt(pi) R N E_global sqrt(mu_win); for R=2 the coefficient is 24 c sqrt(pi)",
        "native_tradeoff": "sqrt(mu_win)/(p_scale p_time)>=nu D_tail/[48 c sqrt(pi) R N E_global]",
        "log_coordinate": "equivalently sqrt(mu_win) exp(H_inf^out+H_inf^time)>=nu D_tail/[48 c sqrt(pi) R mathcal_E_N], mathcal_E_N=N E_global",
        "hard_tail_progress": "selected hard shell M=2^jN with j>=1 gives forward scale ratio>=2 and natural-time ratio<=1/4 independently of signed-good Young geometry",
        "generic_shell": "the actual window peak is a hard-shell event and enters the existing material-free critical-shell first-stop theorem; own-scale service remains conditional on a full no-hit corridor",
        "scope": "no packet persistence, fixed time bins, signed-good parent ratio, Young near-extremality, productivity gate, or additive reset is asserted",
    }


@dataclass(frozen=True)
class NaturalWindowStress:
    samples: int
    worst_time_translation_invariance_residual: float
    worst_time_unit_invariance_residual: float
    worst_representation_refinement_residual: float
    minimum_scale_time_tradeoff_relative_margin: float
    minimum_forward_scale_margin: float
    minimum_natural_time_shortening_margin: float
    minimum_full_survivor_service_lower: float


def _random_partition_density(
    rng: np.random.Generator,
    block_start: float,
    block_end: float,
    pieces: int,
    total_work: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cuts = np.sort(rng.uniform(block_start, block_end, size=max(0, pieces - 1)))
    points = np.concatenate(([block_start], cuts, [block_end]))
    a = points[:-1]
    b = points[1:]
    raw = rng.lognormal(mean=0.0, sigma=1.2, size=pieces)
    integral = float(np.dot(raw, b - a))
    d = raw * (total_work / integral)
    return a, b, d


def _refine_segments(
    a: np.ndarray,
    b: np.ndarray,
    d: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aa: list[float] = []
    bb: list[float] = []
    dd: list[float] = []
    for left, right, rho in zip(a, b, d):
        mid = 0.5 * (float(left) + float(right))
        aa.extend((float(left), mid))
        bb.extend((mid, float(right)))
        dd.extend((float(rho), float(rho)))
    return np.asarray(aa), np.asarray(bb), np.asarray(dd)


def stress(samples: int = 50_000, seed: int = 20260810) -> NaturalWindowStress:
    rng = np.random.default_rng(seed)
    wt = wu = wr = 0.0
    mm = ms = mt = float("inf")
    min_service = float("inf")

    for _ in range(samples):
        D = float(math.exp(rng.uniform(-7.0, 2.0)))
        nu = float(rng.uniform(0.05, 2.5))
        N = float(math.exp(rng.uniform(-2.0, 4.0)))
        c = float(rng.uniform(0.08, 1.6))
        R = float(rng.uniform(1.5, 3.5))
        j_star = int(rng.integers(1, 5))
        n_shell = int(rng.integers(j_star, j_star + 5))
        levels = np.arange(1, n_shell + 1)

        # Make j_star the unique maximal output-scale atom without any coherent
        # representation.  Total HH realizes the upstream clean owner lower.
        base = rng.uniform(0.05, 0.8, size=n_shell)
        base[j_star - 1] = float(base.max() + rng.uniform(0.2, 1.0))
        probs = base / float(base.sum())
        H = 0.5 * nu * D * float(rng.uniform(1.0, 2.2))
        work = {int(k): float(H * p) for k, p in zip(levels, probs)}
        p_s = float(probs[j_star - 1])

        # Choose block peak far below the locality balanced mass threshold, so the
        # comparable-HH owner is certainly present.  The temporal window peak will
        # later be chosen below this same actual block peak.
        C_R = ultraviolet_hh_work_constant() / math.sqrt(R * (R - 1.0))
        clean_weighted_mass = (0.25 * nu / C_R) ** 2
        mu_block = 0.02 * clean_weighted_mass * p_s * p_s
        masses = {int(k): float(math.exp(rng.uniform(-10.0, 0.0))) for k in levels}
        masses[j_star] = mu_block
        locality = high_tail_hh_locality_tradeoff(D, nu, work, masses, R)
        if COMPARABLE_OWNER not in tuple(locality["joint_clean_owners"]):
            raise AssertionError("stress construction lost the comparable-HH locality owner")

        W_lower = float(locality["comparable_parent_common_work_lower"])
        W = W_lower * float(rng.uniform(1.0, 1.7))
        geom = natural_window_geometry(N, j_star, c)
        T_parent = float(geom["parent_natural_duration"])
        T_child = float(geom["selected_natural_window"])
        shift0 = float(rng.uniform(-3.0, 3.0)) * T_parent
        pieces = int(rng.integers(2, 9))
        a, b, d = _random_partition_density(rng, shift0, shift0 + T_parent, pieces, W)
        slide = sliding_window_piecewise_constant(a, b, d, shift0, shift0 + T_parent, T_child)
        Ww = float(slide["maximum_window_common_work"])

        # Choose global energy large enough that the physical window capacity can
        # be realized strictly below the already certified block peak.
        C0 = 4.0 * R * ultraviolet_hh_work_constant() * c * N
        target_mu = 0.35 * mu_block
        E_needed = Ww / (C0 * math.sqrt(target_mu))
        E = E_needed * float(rng.uniform(1.05, 2.5))
        capacity_coeff = C0 * E
        mu_required = (Ww / capacity_coeff) ** 2
        mu_window = min(mu_block, mu_required * float(rng.uniform(1.0, 1.8)))
        if mu_window + 1e-14 * max(mu_window, mu_required, 1e-300) < mu_required:
            raise AssertionError("stress construction lost natural-window capacity")

        out = comparable_hh_temporal_shell_reentry(
            locality,
            D,
            nu,
            N,
            E,
            c,
            W,
            Ww,
            T_child,
            mu_window,
        )
        clean = float(out["clean_entropy_weighted_sqrt_child_mass_lower"])
        got = float(out["entropy_weighted_sqrt_child_mass"])
        mm = min(mm, got / max(clean, 1e-300) - 1.0)
        ms = min(ms, float(out["forward_scale_ratio"]) - 2.0)
        mt = min(mt, 0.25 - float(out["natural_time_ratio"]))
        min_service = min(min_service, float(out["full_survivor_own_scale_service_lower"]))

        # Time-origin translation gauge.
        delta = float(rng.uniform(-20.0, 20.0)) * T_parent
        moved = sliding_window_piecewise_constant(a + delta, b + delta, d, shift0 + delta, shift0 + delta + T_parent, T_child)
        wt = max(wt, abs(float(moved["p_time"]) - float(slide["p_time"])))

        # Unit rescaling t'=lambda t, rho'=rho/lambda leaves the same measure.
        lam = float(math.exp(rng.uniform(-5.0, 5.0)))
        scaled = sliding_window_piecewise_constant(
            lam * a,
            lam * b,
            d / lam,
            lam * shift0,
            lam * (shift0 + T_parent),
            lam * T_child,
        )
        wu = max(wu, abs(float(scaled["p_time"]) - float(slide["p_time"])))

        # Representation refinement with unchanged positive measure.
        ar, br, dr = _refine_segments(a, b, d)
        refined = sliding_window_piecewise_constant(ar, br, dr, shift0, shift0 + T_parent, T_child)
        wr = max(wr, abs(float(refined["p_time"]) - float(slide["p_time"])))

    return NaturalWindowStress(samples, wt, wu, wr, mm, ms, mt, min_service)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-tail-natural-window-reentry"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "high_tail_natural_window_reentry.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# High-tail comparable HH work -> sliding natural-window shell reentry

Status: **{cert['status']}**.

The ultraviolet theorem first isolates actual positive comparable HH work on one hard output shell `M=2^jN`, `j>=1`.  Time is then read from **that same positive work measure**, not from a packet history and not from a chosen time grid.

Let `T_M=cM^-2` and define

`p_time = sup_s mu_comp([s,s+T_M]) / mu_comp(I)`.

For smooth Navier--Stokes the positive comparable-work density is continuous, so the sliding integral is continuous on the compact admissible start interval and its maximum is attained.  Translation of the time origin, rescaling the time unit, or refining a representation of the same measure leaves `p_time` unchanged.

The canonical strict low-pass currently gives only `|S|<=1`, hence `||h||_2<=2||u||_2`.  For comparable parents `<=R M`, sharp physical Young and the global physical energy cap therefore give the common-unit natural-window capacity

`N W_window <= 12 c sqrt(pi) R N E_global sqrt(mu_window)`.

Combining this with the certified locality owner

`W_comp exp(H_inf^out) >= nu D_tail/4`

and `W_window=p_time W_comp` yields the native scale--time relation

`sqrt(mu_window)/(p_scale p_time) >= nu D_tail/[48 c sqrt(pi) R N E_global]`.

Equivalently, with `mathcal_E_N=N E_global`,

`sqrt(mu_window) exp(H_inf^out+H_inf^time) >= nu D_tail/[48 c sqrt(pi) R mathcal_E_N]`.

For the dyadic locality corollary `R=2`, the denominator is `96 c sqrt(pi) mathcal_E_N`.

The peak `mu_window` is an actual hard-shell event, so it enters the existing generic critical-shell first-stop theorem.  Its full-no-hit own-scale service lower is recorded only conditionally.  In addition, hard-tail support itself gives genuine forward scale progress `M/N=2^j>=2` and natural-time shortening `T_M/T_N=4^-j<=1/4`; neither fact uses signed-good Young geometry.

Stress: `{out.samples}` sliding-measure / scale-time / shell-reentry states
- worst time-origin invariance residual: `{out.worst_time_translation_invariance_residual:.3e}`
- worst time-unit invariance residual: `{out.worst_time_unit_invariance_residual:.3e}`
- worst representation-refinement residual: `{out.worst_representation_refinement_residual:.3e}`
- minimum scale-time tradeoff relative margin: `{out.minimum_scale_time_tradeoff_relative_margin:.3e}`
- minimum forward scale-ratio margin above 2: `{out.minimum_forward_scale_margin:.3e}`
- minimum natural-time shortening margin below 1/4: `{out.minimum_natural_time_shortening_margin:.3e}`
- minimum conditional full-survivor own-scale service lower: `{out.minimum_full_survivor_service_lower:.3e}`

No packet persistence, fixed time bins, generated-energy productivity gate, or additive reset is introduced.  No 3D Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
