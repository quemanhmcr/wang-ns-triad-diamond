from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

from src.ancestor_reservoir_sync import (
    CLEAN_GENERATION_PROGRESS,
    CLEAN_RESERVOIR_GROWTH,
    LOW_STRAIN_ACTION,
)
from src.high_strain_heat_increment_service import (
    high_strain_heat_service_lower,
    heat_increment_multiplier,
)

# The first-hit history gives the sharper physical Kelvin factor exp(1/30).
# The rational constants are only a clean certified envelope already present in
# ancestor_reservoir_sync.
PHYSICAL_HEAT_POOL_RATIO_UPPER = math.exp(2.0 * float(LOW_STRAIN_ACTION)) / float(CLEAN_GENERATION_PROGRESS)
CLEAN_HEAT_POOL_RATIO = CLEAN_RESERVOIR_GROWTH**2 / CLEAN_GENERATION_PROGRESS


def first_hit_kelvin_growth_upper(strain_action: float) -> float:
    """Kelvin frequency growth on a first-high-strain pre-hit history.

    The first high-strain contact has K=1/30 and K(t)<=1/30 at every earlier
    time.  A value above the canonical threshold is therefore not a first-hit
    history for this theorem.
    """
    K = float(strain_action)
    threshold = float(LOW_STRAIN_ACTION)
    if K < 0 or K > threshold + 2e-15 or not math.isfinite(K):
        raise ValueError("first-hit strain action must lie in [0,1/30]")
    return math.exp(K)


def heat_band_service_upper(
    *,
    child_frequency: float,
    band_frequency_upper: float,
    scaled_lifetime: float,
    energy_upper: float,
) -> float:
    """Normalized heat service capacity of a deterministic band |xi|<=M.

    For the NS heat probe at time 1/(2N^2),

      E_H ||delta_r f||_2^2 <= (M/N)^2 ||f||_2^2.

    On the canonical natural lifetime T=cN^-2 this yields

      N^3 int_I E_H||delta_r f||_2^2 dt <= c M^2 E/N.

    The bound is applied shellwise before material ownership is summed; it does
    not construct an ``old velocity field``.
    """
    N = float(child_frequency)
    M = float(band_frequency_upper)
    c = float(scaled_lifetime)
    E = float(energy_upper)
    if N <= 0 or M < 0 or c <= 0 or E < 0 or not all(math.isfinite(x) for x in (N, M, c, E)):
        raise ValueError("positive finite N,c and nonnegative finite M,E required")
    return c * M * M * E / N


def old_pool_heat_capacity_upper(
    *,
    generation: int,
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    frame_energy_bound: float,
    global_energy: float,
    scaled_lifetime: float = 1.0,
) -> float:
    """Clean whole-old-pool OO heat capacity on a supplied signed-good epoch.

    At material age q, first-hit Kelvin transport and signed-good scale progress
    give the sharper coefficient ratio

      (5/8) exp(1/15) < 441/640.

    We use the clean rational envelope 441/640.  The old pool is adversarially
    allowed the whole frame-energy budget P E_global at every service time.
    Since OO is selected only after the positive shell heat/Moyal law exists,
    its shellwise mass cannot exceed that shell's total positive heat service.
    Summing the deterministic orthogonal shell capacities therefore introduces
    no old/new field decomposition and no quadratic cross term.
    """
    q = int(generation)
    alpha = float(initial_low_cut_ratio)
    N0 = float(initial_block_frequency)
    P = float(frame_energy_bound)
    E = float(global_energy)
    c = float(scaled_lifetime)
    if q < 0:
        raise ValueError("nonnegative generation required")
    if min(alpha, N0, P, c) <= 0 or E < 0 or not all(math.isfinite(x) for x in (alpha, N0, P, E, c)):
        raise ValueError("positive finite pool geometry and nonnegative finite energy required")
    base = c * alpha * alpha * N0 * P * E
    return base * float(CLEAN_HEAT_POOL_RATIO) ** q


def total_future_old_pool_heat_capacity_upper(
    *,
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    frame_energy_bound: float,
    global_energy: float,
    scaled_lifetime: float = 1.0,
) -> float:
    base = old_pool_heat_capacity_upper(
        generation=0,
        initial_low_cut_ratio=initial_low_cut_ratio,
        initial_block_frequency=initial_block_frequency,
        frame_energy_bound=frame_energy_bound,
        global_energy=global_energy,
        scaled_lifetime=scaled_lifetime,
    )
    r = float(CLEAN_HEAT_POOL_RATIO)
    return base / (1.0 - r)


def first_forced_non_oo_generation(
    *,
    high_strain_heat_threshold: float,
    initial_old_capacity: float,
    forced_non_oo_fraction: float = 0.5,
) -> int:
    """First material age where OO cannot cover (1-f) of a high-strain event."""
    S = float(high_strain_heat_threshold)
    C0 = float(initial_old_capacity)
    f = float(forced_non_oo_fraction)
    if S <= 0 or C0 < 0 or not (0 < f < 1) or not all(math.isfinite(x) for x in (S, C0, f)):
        raise ValueError("positive threshold, nonnegative capacity and fraction in (0,1) required")
    target = (1.0 - f) * S
    if C0 <= target:
        return 0
    r = float(CLEAN_HEAT_POOL_RATIO)
    q = max(0, int(math.ceil(math.log(target / C0) / math.log(r))))
    # Remove floating endpoint ambiguity while preserving the exact monotone
    # definition of the first stopping generation.
    while C0 * r**q > target:
        q += 1
    while q > 0 and C0 * r ** (q - 1) <= target:
        q -= 1
    return q


def forced_non_oo_service_lower(*, total_heat_service: float, old_old_service: float) -> float:
    """Exact positive remainder S_ON+S_NN=S_heat-S_OO."""
    S = float(total_heat_service)
    OO = float(old_old_service)
    if S < 0 or OO < 0 or OO > S or not math.isfinite(S + OO):
        raise ValueError("finite 0<=OO<=total heat service required")
    return S - OO


def first_hit_epoch_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    clean = CLEAN_HEAT_POOL_RATIO
    if clean != Fraction(441, 640):
        raise AssertionError("clean heat-reservoir ratio algebra changed")
    if not (PHYSICAL_HEAT_POOL_RATIO_UPPER < float(clean) < 0.7):
        raise AssertionError("physical/clean heat-reservoir contraction ordering failed")
    total_factor = Fraction(1, 1) / (Fraction(1, 1) - clean)
    if total_factor != Fraction(640, 199) or not (total_factor < Fraction(13, 4)):
        raise AssertionError("clean heat-reservoir geometric budget changed")
    Sstar = high_strain_heat_service_lower(c)
    return {
        "status": "EXACT_FIRST_HIT_HEAT_RESERVOIR_EROSION__CLEAN_OO_RATIO_441_640__FINITE_NON_OO_FORCE_ON_SUPPLIED_SIGNED_GOOD_EPOCH",
        "first_hit": "a first high-strain contact occurs at K=1/30 with K(t)<=1/30 on its entire closed pre-hit history, so Kelvin growth is <=exp(1/30)",
        "band_heat_capacity": "for deterministic supp fhat subset B_M and T=cN^-2, N^3 int_I E_H||delta_r f||_2^2 dt <= c M^2 sup_I||f||_2^2/N",
        "positive_submeasure": "apply the band estimate shellwise to the already-existing positive heat/Moyal law; OO is a positive material submeasure, so no old/new velocity-field decomposition is introduced",
        "epoch_geometry": "on a supplied signed-good material epoch, material frequencies grow by <=exp(1/30) per first-hit block while N grows by >8/5",
        "physical_ratio": f"rho_phys <= (5/8)exp(1/15)={PHYSICAL_HEAT_POOL_RATIO_UPPER:.12g}",
        "clean_ratio": "rho_phys < 441/640 < 7/10",
        "total_old_budget": "sum_q C_OO(q) <= (640/199) C_OO(0) < (13/4) C_OO(0)",
        "high_strain_threshold": f"every first high-strain contact has S_heat>=S_*={Sstar:.12g} for the fixed natural-lifetime coefficient c={c:.12g}",
        "finite_force": "after the first q with C_OO(q)<=(1-f)S_*, exact ownership forces S_ON+S_NN>=f S_*",
        "scope": "this closes repeated-OO heat capacity only after the PDE has supplied a signed-good material recursive epoch with the canonical T(N)=cN^-2; it does not prove universal slab renewal or the final ON versus NN destination",
    }


@dataclass(frozen=True)
class HeatReservoirStress:
    samples: int
    maximum_spectral_capacity_ratio: float
    maximum_physical_epoch_ratio: float
    minimum_clean_ratio_margin: float
    minimum_first_hit_growth_margin: float
    minimum_forced_generation_margin: float
    minimum_previous_generation_margin: float
    minimum_non_oo_margin: float


def stress(samples: int = 50_000, seed: int = 20260809) -> HeatReservoirStress:
    rng = np.random.default_rng(seed)
    max_spec = max_epoch = 0.0
    min_ratio = min_growth = min_force = min_prev = min_nonoo = float("inf")
    clean_growth = float(CLEAN_RESERVOIR_GROWTH)
    clean_progress = float(CLEAN_GENERATION_PROGRESS)
    rclean = float(CLEAN_HEAT_POOL_RATIO)

    for _ in range(samples):
        # Direct spectral regression of
        # 2(1-e^-x)<=|xi|^2/N^2<=M^2/N^2.
        N = float(math.exp(rng.uniform(-2.0, 6.0)))
        M = float(rng.uniform(0.0, 0.25)) * N
        n = int(rng.integers(2, 100))
        k = rng.random(n) * M
        e = rng.lognormal(mean=-1.0, sigma=1.5, size=n)
        exact = float(np.dot(np.array([heat_increment_multiplier(x, N) / (N * N) for x in k]), e))
        E = float(e.sum())
        upper = (M / N) ** 2 * E
        if upper > 0:
            max_spec = max(max_spec, exact / upper)
        if exact > upper + 4e-12 * max(1.0, upper):
            raise AssertionError("band heat increment exceeded M^2/N^2 energy capacity")

        # First-hit boundary remains inside the closed low-strain transport history.
        K = float(rng.uniform(0.0, float(LOW_STRAIN_ACTION)))
        growth = first_hit_kelvin_growth_upper(K)
        min_growth = min(min_growth, clean_growth - growth)
        if growth >= clean_growth + 2e-14:
            raise AssertionError("first-hit Kelvin growth exceeded clean 21/20 envelope")

        # Use the actual first-hit exponential growth, not the rational envelope,
        # against an arbitrary signed-good scale step.
        mratio = math.exp(float(rng.uniform(0.0, float(LOW_STRAIN_ACTION))))
        nratio = float(rng.uniform(clean_progress * (1.0 + 1e-12), 4.0))
        eratio = mratio * mratio / nratio
        max_epoch = max(max_epoch, eratio)
        min_ratio = min(min_ratio, rclean - eratio)
        if eratio >= rclean + 2e-13:
            raise AssertionError("physical heat-reservoir epoch ratio exceeded clean 441/640 envelope")

        Sstar = float(math.exp(rng.uniform(-3.0, 3.0)))
        C0 = float(math.exp(rng.uniform(-4.0, 5.0)))
        f = float(rng.uniform(0.1, 0.9))
        q = first_forced_non_oo_generation(
            high_strain_heat_threshold=Sstar,
            initial_old_capacity=C0,
            forced_non_oo_fraction=f,
        )
        target = (1.0 - f) * Sstar
        Cq = C0 * rclean**q
        min_force = min(min_force, target - Cq)
        if Cq > target + 3e-13 * max(1.0, target):
            raise AssertionError("forced generation still permits too much OO heat capacity")
        if q > 0:
            Cp = C0 * rclean ** (q - 1)
            min_prev = min(min_prev, Cp - target)
            if Cp <= target - 3e-13 * max(1.0, target):
                raise AssertionError("forced non-OO generation is not minimal")
        else:
            min_prev = min(min_prev, 0.0)

        # S_heat>=Sstar and S_OO<=Cq imply the exact positive remainder is >=fSstar.
        total = float(rng.uniform(1.0, 3.0)) * Sstar
        oo = float(rng.uniform(0.0, 1.0)) * min(Cq, total)
        nonoo = forced_non_oo_service_lower(total_heat_service=total, old_old_service=oo)
        margin = nonoo - f * Sstar
        min_nonoo = min(min_nonoo, margin)
        if margin < -4e-13 * max(1.0, Sstar):
            raise AssertionError("ON+NN remainder fell below forced high-strain fraction")

    return HeatReservoirStress(
        samples,
        max_spec,
        max_epoch,
        min_ratio,
        min_growth,
        min_force,
        min_prev,
        min_nonoo,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-first-hit-heat-reservoir-erosion"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = first_hit_epoch_certificate()
    out = stress(args.samples)
    (args.outdir / "first_hit_heat_reservoir_erosion.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# First-hit heat reservoir erosion: repeated OO service cannot hide forever\n\nStatus: **{cert['status']}**.\n\nA high-strain event is read at its first boundary contact, not by retrospectively declaring its whole history high strain.  If `K(t)=int ||S_V||_op`, then at the first contact `K=1/30` and on the entire closed pre-hit history `K<=1/30`.  Kelvin transport therefore gives `M(t)/M(t0)<=exp(1/30)`.\n\nFor a deterministic band `|xi|<=M`, the NS heat probe at time `1/(2N^2)` obeys `E_H||delta_r f||_2^2 <= (M/N)^2||f||_2^2`.  On the canonical natural lifetime `T(N)=cN^-2`,\n\n`S_heat(f)=N^3 int_I E_H||delta_r f||_2^2dt <= c M^2 E/N`.\n\nApply this to the deterministic resolved shell law before material summation.  Moyal makes that shell service positive, and OO ownership is then merely a positive submeasure.  Hence OO cannot exceed the shell's total heat service.  No velocity decomposition `u=u_old+u_new` appears.\n\nOn any supplied signed-good material epoch the block scale advances by more than `8/5`, while a reused material frequency grows by at most `exp(1/30)` during each first-hit history.  Thus the physical heat-capacity coefficient contracts by\n\n`rho_phys <= (5/8)exp(1/15) = {PHYSICAL_HEAT_POOL_RATIO_UPPER:.12g}`.\n\nThe existing rational Kelvin envelope gives the clean bound\n\n`rho_phys < (21/20)^2/(8/5) = 441/640 < 7/10`.\n\nTherefore\n\n`C_OO(q) <= c alpha^2 N_0 P E_global (441/640)^q`,\n\nand the entire future old--old heat capacity satisfies\n\n`sum_(q>=0) C_OO(q) <= (640/199) C_OO(0) < (13/4) C_OO(0)`.\n\nEvery first high-strain contact simultaneously has the existing normalized heat lower `S_heat>=S_*>0`.  Once material age reaches the first `q_*` with `C_OO(q_*)<=(1-f)S_*`, exact positive ownership gives\n\n`S_ON+S_NN=S_heat-S_OO >= fS_*`.\n\nFor `f=1/2`, at least half of every sufficiently old high-strain heat event is therefore genuine interface-or-new service.  The statement consumes no new currency and assumes no packet mass floor.\n\nStress: `{out.samples}` spectral/first-hit/epoch states\n- maximum exact band heat service / analytic capacity ratio: `{out.maximum_spectral_capacity_ratio:.9f}`\n- maximum sampled physical one-step epoch ratio: `{out.maximum_physical_epoch_ratio:.9f}`\n- minimum clean `441/640` ratio margin: `{out.minimum_clean_ratio_margin:.3e}`\n- minimum first-hit Kelvin-growth margin to `21/20`: `{out.minimum_first_hit_growth_margin:.3e}`\n- minimum forced-generation capacity margin: `{out.minimum_forced_generation_margin:.3e}`\n- minimum previous-generation minimality margin: `{out.minimum_previous_generation_margin:.3e}`\n- minimum forced ON+NN service margin: `{out.minimum_non_oo_margin:.3e}`\n\nThis closes quantitative repeated-OO heat capacity **once the PDE has supplied the signed-good material epoch with the canonical natural lifetime**.  Universal slab renewal remains open.  The forced remainder is still only `ON+NN`; attaching ON to temporal material-interface/relink work and deriving the correct NN renewal law are the next continuum questions.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
