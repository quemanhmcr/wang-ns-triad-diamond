from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

PARENT_CHILD_LO = Fraction(3, 5)
PARENT_CHILD_HI = Fraction(5, 8)
LIFETIME_GROWTH_MIN = Fraction(64, 25)
LIFETIME_GROWTH_MAX = Fraction(25, 9)
TIME_CONTRACTION_MAX = Fraction(25, 64)
INITIAL_HALF_SPAN = Fraction(25, 128)
BACKWARD_FRACTION = Fraction(2, 5)
SYNC_CONE = Fraction(3, 8)  # legacy loose cone; useful as a coarse check
SYNC_FIXED_POINT = Fraction(10, 39)
SHARP_SYNC_BOUND = SYNC_FIXED_POINT
COMMON_SLICE_MARGIN = Fraction(67, 195)
MIN_REFERENCE_BACKSTEP = Fraction(1792, 4875)
BACKWARD_SUM_COEFF = Fraction(1792, 7605)


def lifetime_ratio_bounds() -> tuple[Fraction, Fraction]:
    return LIFETIME_GROWTH_MIN, LIFETIME_GROWTH_MAX


def initial_parent_span_ratio() -> Fraction:
    """Half a child slab, normalized by the shortest signed-good parent lifetime."""
    return Fraction(1, 2) * TIME_CONTRACTION_MAX


def next_span_ratio(alpha: float, backward_fraction: float = float(BACKWARD_FRACTION)) -> float:
    if alpha < 0 or backward_fraction < 0:
        raise ValueError("nonnegative span data required")
    return float(TIME_CONTRACTION_MAX) * (alpha + backward_fraction)


def synchronization_cone_margin(alpha: float = float(SHARP_SYNC_BOUND)) -> float:
    """Margin to the edge of a natural packet lifetime for the common slice."""
    return 1.0 - alpha - float(BACKWARD_FRACTION)


def common_reference_slice(interval_start: float, interval_end: float, min_lifetime: float) -> float:
    if min_lifetime <= 0 or interval_end < interval_start:
        raise ValueError("invalid interval/lifetime")
    alpha = (interval_end - interval_start) / min_lifetime
    if alpha > float(SHARP_SYNC_BOUND) + 1e-14:
        raise ValueError("event slab lies outside the sharp synchronization cone")
    s = interval_start - float(BACKWARD_FRACTION) * min_lifetime
    if interval_end - s >= min_lifetime - 1e-14 * max(1.0, min_lifetime):
        raise AssertionError("common reference slice left a natural backward window")
    return s


def recursive_span_sequence(alpha0: float, layers: int) -> list[float]:
    if layers < 0:
        raise ValueError("nonnegative layer count required")
    out = [float(alpha0)]
    for _ in range(layers):
        out.append(next_span_ratio(out[-1]))
    return out


def minimum_backward_displacement(initial_lifetime: float, layers: int) -> float:
    """Correct minimum decrease of the common reference slices after L steps.

    If H_j=[a_j,b_j] and s_j=a_j-(2/5)T_j, generated parent support only obeys
    H_(j+1) subset [s_j,b_j]; it need not start at s_j.  Since alpha_j<=10/39
    and T_(j+1)>=(64/25)T_j,

      s_j-s_(j+1) >= (1792/4875) T_j.

    Summing the geometric lifetime growth gives the coefficient 1792/7605.
    """
    if initial_lifetime <= 0 or layers < 0:
        raise ValueError("positive lifetime and nonnegative layer count required")
    if layers == 0:
        return 0.0
    g = float(LIFETIME_GROWTH_MIN)
    return float(BACKWARD_SUM_COEFF) * initial_lifetime * (g**layers - 1.0)


def interior_depth_upper(current_time: float, initial_lifetime: float) -> int:
    """Largest number of certified interior backward steps before t=0 must be met."""
    if current_time < 0 or initial_lifetime <= 0:
        raise ValueError("invalid time/lifetime")
    g = float(LIFETIME_GROWTH_MIN)
    x = 1.0 + current_time / (float(BACKWARD_SUM_COEFF) * initial_lifetime)
    return max(0, int(math.floor(math.log(x) / math.log(g) + 1e-14)))


def initial_root_count_upper(
    band_frequency: float,
    homogeneous_sobolev_norm_sq: float,
    critical_mass_threshold: float,
    sobolev_order: float,
) -> float:
    """Initial-boundary root count from an H^m Fourier/coherent band budget.

    ||P_M u_0||_2^2 <= M^{-2m} ||u_0||_{Hdot^m}^2 and every root uses
    E_root >= eta/M, hence #roots <= eta^{-1} M^{1-2m} ||u_0||^2.
    """
    if band_frequency <= 0 or homogeneous_sobolev_norm_sq < 0 or critical_mass_threshold <= 0:
        raise ValueError("invalid initial-boundary data")
    return (
        homogeneous_sobolev_norm_sq
        * band_frequency ** (1.0 - 2.0 * sobolev_order)
        / critical_mass_threshold
    )


def registration_xi_upper(
    switch_moyal_energy: float,
    total_energy: float,
    covariance_log_distance: float,
) -> float:
    """Single-boundary coherent registration charge.

    Exact common affine/Kelvin motion is free. A material cell-selection switch
    costs its symmetric-difference Moyal energy, while changing the normalized
    Gaussian covariance window costs d_log ||f||_2^2/sqrt(2).
    """
    if min(switch_moyal_energy, total_energy, covariance_log_distance) < 0:
        raise ValueError("nonnegative registration data required")
    return switch_moyal_energy + covariance_log_distance * total_energy / math.sqrt(2.0)


def choose_heavy_half(times: Sequence[float], weights: Sequence[float], slab_start: float, slab_end: float) -> dict[str, float | int]:
    """Pigeonhole positive generation mass into one half-child slab."""
    t = np.asarray(times, float)
    w = np.asarray(weights, float)
    if t.shape != w.shape or t.ndim != 1 or np.any(w < 0) or slab_end <= slab_start:
        raise ValueError("invalid positive generation law")
    if np.any(t < slab_start) or np.any(t > slab_end):
        raise ValueError("event outside child slab")
    mid = 0.5 * (slab_start + slab_end)
    left = float(np.sum(w[t <= mid]))
    right = float(np.sum(w[t > mid]))
    total = left + right
    if left >= right:
        return {"half": 0, "start": slab_start, "end": mid, "mass": left, "total": total}
    return {"half": 1, "start": mid, "end": slab_end, "mass": right, "total": total}


@dataclass(frozen=True)
class AsyncSyncStress:
    samples: int
    worst_cone_ratio: float
    minimum_cone_margin: float
    minimum_common_slice_margin: float
    minimum_half_mass_margin: float
    minimum_initial_boundary_margin: float
    minimum_registration_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> AsyncSyncStress:
    rng = np.random.default_rng(seed)
    worst = 0.0
    mcone = mslice = mhalf = minit = mreg = float("inf")
    cone = float(SHARP_SYNC_BOUND)
    for _ in range(samples):
        # Start from the certified first parent-event span and perturb downward.
        alpha = float(INITIAL_HALF_SPAN) * float(rng.uniform(0.0, 1.0))
        depth = int(rng.integers(1, 30))
        seq = recursive_span_sequence(alpha, depth)
        for a in seq:
            worst = max(worst, a / cone)
            mcone = min(mcone, cone - a)
            if a > cone + 2e-14:
                raise AssertionError("parabolic synchronization cone failed")
            mslice = min(mslice, synchronization_cone_margin(a))
            if synchronization_cone_margin(a) <= 0:
                raise AssertionError("common slice escaped natural lifetime")

        # Direct physical-time check for one common slice.
        Tmin = float(math.exp(rng.uniform(-5.0, 3.0)))
        a0 = float(rng.uniform(0.0, 10.0))
        width = float(rng.uniform(0.0, cone)) * Tmin
        b0 = a0 + width
        s = common_reference_slice(a0, b0, Tmin)
        mslice = min(mslice, Tmin - (b0 - s))

        # Correct asynchronous geometry: the next generated support may begin
        # anywhere in [s,b0], not necessarily at s.  Even in the worst case
        # a_next=b0, the next common reference slice moves left by a fixed
        # fraction of the current minimum lifetime.
        Tnext = float(LIFETIME_GROWTH_MIN) * Tmin * float(rng.uniform(1.0, 1.4))
        anext = float(rng.uniform(s, b0))
        snext = anext - float(BACKWARD_FRACTION) * Tnext
        backstep = s - snext
        required = float(MIN_REFERENCE_BACKSTEP) * Tmin
        minit = min(minit, backstep - required)
        if backstep + 2e-12 * max(1.0, Tmin) < required:
            raise AssertionError("asynchronous common slices did not move backward fast enough")

        # Positive aligned Duhamel mass loses at most one factor 1/2 initially.
        n = int(rng.integers(2, 80))
        slab0 = float(rng.uniform(-2.0, 2.0))
        T = float(math.exp(rng.uniform(-3.0, 2.0)))
        times = slab0 + rng.random(n) * T
        weights = rng.lognormal(mean=-1.0, sigma=1.0, size=n)
        hh = choose_heavy_half(times, weights, slab0, slab0 + T)
        margin = float(hh["mass"]) - 0.5 * float(hh["total"])
        mhalf = min(mhalf, margin)
        if margin < -2e-12 * max(1.0, float(hh["total"])):
            raise AssertionError("half-slab generation pigeonhole failed")

        # The geometric lower displacement really forces t=0 after finite depth.
        T0 = float(math.exp(rng.uniform(-8.0, 0.0)))
        current = float(math.exp(rng.uniform(-2.0, 5.0)))
        L = interior_depth_upper(current, T0)
        dL = minimum_backward_displacement(T0, L)
        dnext = minimum_backward_displacement(T0, L + 1)
        minit = min(minit, current - dL, dnext - current)
        if dL > current + 2e-11 * max(1.0, current):
            raise AssertionError("interior depth upper bound undershot t=0")
        if dnext <= current - 2e-11 * max(1.0, current):
            raise AssertionError("initial boundary was not forced at next layer")

        # Registration is a one-boundary positive charge, never below either term.
        Esw = float(rng.lognormal(mean=-3.0, sigma=1.0))
        E = float(rng.lognormal(mean=0.0, sigma=1.0))
        dlog = float(rng.uniform(0.0, 0.2))
        xi = registration_xi_upper(Esw, E, dlog)
        mreg = min(mreg, xi - Esw, xi - dlog * E / math.sqrt(2.0))
        if xi + 1e-14 < max(Esw, dlog * E / math.sqrt(2.0)):
            raise AssertionError("registration charge lost one interface term")

    return AsyncSyncStress(samples, worst, mcone, mslice, mhalf, minit, mreg)


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_PARABOLIC_ASYNC_SYNCHRONIZATION_GIVEN_ONE_STEP_ADJOINT_GATE",
        "scale_window": "3/5 < N_parent/N_child < 5/8",
        "lifetime_window": "64/25 < T_parent/T_child < 25/9",
        "initial_parent_span": "alpha_1 <= 25/128",
        "recurrence": "alpha_next <= (25/64)(alpha+2/5)",
        "loose_invariant_cone": "alpha <= 3/8",
        "sharp_post_half_slab_bound": "alpha <= 10/39",
        "common_slice_margin": "1-10/39-2/5 = 67/195",
        "fixed_point": "10/39",
        "one_step_reference_backshift": "s_j-s_(j+1) >= (1792/4875) T_j",
        "backward_displacement": "Delta s_L >= (1792/7605)T_0[(64/25)^L-1]",
        "registration": "common Kelvin transport free; Xi_boundary=E_switch+d_log E/sqrt(2)",
        "initial_boundary": "truncate exact adjoint Duhamel gate at t=0; do not count as interior fresh grain",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-asynchronous-duhamel-sync"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "asynchronous_duhamel_sync.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Asynchronous Duhamel synchronization cone\n\nStatus: **{cert['status']}**.\n\nThe signed-good scale window gives `64/25<T_parent/T_child<25/9`. After retaining a half-child slab with at least half of the positive aligned high--high Duhamel mass, parent events occupy normalized span at most `25/128` of the shortest parent lifetime.\n\nFor every later generated layer choose the common reference slice `s_j=a_j-(2/5)T_j^min`. If `alpha_j=W_j/T_j^min`, then\n\n`alpha_(j+1) <= (25/64)(alpha_j+2/5)`.\n\nThe loose cone `alpha_j<=3/8` is invariant, but after the first half-slab the sharper bound `alpha_j<=10/39` holds because `25/128<10/39` and `10/39` is the recurrence fixed point.  The common-window margin is therefore `67/195`. Therefore all **generated** nodes can be placed on common physical reference slices by exact adjoint Duhamel gates; inherited nodes stop as old material roots and classified residual nodes stop in the existing interface/source ledger. No frozen-packet persistence hypothesis is used.\n\nBecause generated parent support may start anywhere inside the previous common interval, the corrected one-step reference-slice shift is `s_j-s_(j+1)>=(1792/4875)T_j`.  Hence the minimum cumulative backward displacement after `L` interior synchronization steps is\n\n`(1792/7605)T_0[(64/25)^L-1]`.\n\nThus at finite physical time an ancestry either stops earlier or reaches `t=0` after finite depth. The gate is then truncated exactly at the initial surface, where the node is an initial-boundary root rather than a fresh interior grain. For `u_0 in Hdot^m`, a band-`M` initial root family with critical mass `M E>=eta` obeys `#roots <= eta^-1 M^(1-2m)||u_0||_Hdot^m^2`.\n\nMaterial coherent labels are back-transported to the common slice. Exact common Kelvin/affine motion is free. At a layer boundary, a cell-selection switch is charged once by its symmetric-difference Moyal energy and a covariance-window change by `d_log E/sqrt(2)`; large covariance changes remain in the existing strain/fresh/source branch.\n\nStress: `{out.samples}` randomized layer/slab/boundary checks\n- worst synchronization-cone ratio: `{out.worst_cone_ratio:.9f}`\n- minimum cone margin: `{out.minimum_cone_margin:.3e}`\n- minimum common-slice margin: `{out.minimum_common_slice_margin:.3e}`\n- minimum half-mass pigeonhole margin: `{out.minimum_half_mass_margin:.3e}`\n- minimum initial-boundary bracket margin: `{out.minimum_initial_boundary_margin:.3e}`\n- minimum registration margin: `{out.minimum_registration_margin:.3e}`\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
