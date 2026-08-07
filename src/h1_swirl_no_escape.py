from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

PHYSICAL_H1_COST = Fraction(1, 184320)
FULL_MILD_COST = Fraction(1, 737280)


def physical_first_impulse_lower(hook_impulse: float, K: float = 1 / 30) -> float:
    """Physical direct-sum first Duhamel daughter on the low-strain branch.

    Three-role physical forcing norm is >= ||B_hook||/10. Pulling to the
    interaction picture costs exp(-K), and pushing the impulse back to physical
    roles costs a second exp(-K).
    """
    if hook_impulse < 0 or K < 0:
        raise ValueError("nonnegative data required")
    return math.exp(-2 * K) * hook_impulse / 20.0


def clean_first_impulse_square_lower(hook_impulse: float) -> float:
    if hook_impulse < 0:
        raise ValueError("nonnegative impulse required")
    return hook_impulse * hook_impulse / 480.0


def h1_physical_quadratic_cost(hook_impulse: float) -> float:
    if hook_impulse < 0:
        raise ValueError("nonnegative impulse required")
    return float(PHYSICAL_H1_COST) * hook_impulse * hook_impulse


def h1_dephasing_threshold(hook_impulse: float, lifetime: float) -> float:
    if hook_impulse < 0 or lifetime <= 0:
        raise ValueError("invalid data")
    return hook_impulse / (11.0 * lifetime)


def classify_h1_mild_no_escape(
    hook_impulse: float,
    lifetime: float,
    interaction_forcing_variation: float,
    physical_first_duhamel_norm: float,
    physical_nonlinear_feedback_norm: float,
    largest_role_critical_sigma: float,
    pair_rescue: float,
    net_transfer_deficit: float,
) -> dict[str, float | str]:
    """Physical H1 no-escape after nonunitary conditioning and role splitting.

    Assumptions: cond(L)<=21/20, each base role has action budget K<=1/30,
    and the certified physical three-role forcing energy is >=||B_hook||^2/100.
    """
    vals = [
        hook_impulse,
        interaction_forcing_variation,
        physical_first_duhamel_norm,
        physical_nonlinear_feedback_norm,
        largest_role_critical_sigma,
        pair_rescue,
        net_transfer_deficit,
    ]
    if lifetime <= 0 or any(v < 0 for v in vals):
        raise ValueError("invalid H1 no-escape data")
    I = hook_impulse
    varth = h1_dephasing_threshold(I, lifetime)
    req = math.sqrt(clean_first_impulse_square_lower(I))
    cost = h1_physical_quadratic_cost(I)
    if interaction_forcing_variation >= varth - 1e-12 * max(1.0, varth):
        branch = "H1_covariant_dephasing"
    else:
        if physical_first_duhamel_norm + 1e-11 * max(1.0, req) < req:
            raise ValueError("coherent physical H1 branch violates conditioned Duhamel lower")
        if physical_nonlinear_feedback_norm >= 0.5 * physical_first_duhamel_norm - 1e-12 * max(1.0, physical_first_duhamel_norm):
            branch = "nonlinear_sideband_feedback"
        else:
            total = max(0.0, physical_first_duhamel_norm - physical_nonlinear_feedback_norm)
            one_role = total / math.sqrt(3.0)
            if largest_role_critical_sigma + 1e-11 < one_role:
                raise ValueError("largest physical-role critical sigma too small")
            if largest_role_critical_sigma >= 1 / 80 - 1e-14:
                branch = "large_daughter_capacity"
            else:
                d0 = largest_role_critical_sigma**2 / 16.0
                if pair_rescue >= 0.5 * d0 - 1e-13 * max(1.0, d0):
                    branch = "pair_sideband_rescue"
                    if pair_rescue + 2e-12 < cost:
                        raise AssertionError("physical H1 pair cost failed")
                else:
                    branch = "transfer_deficit"
                    if net_transfer_deficit + 2e-12 < 0.5 * d0 or net_transfer_deficit + 2e-12 < cost:
                        raise AssertionError("physical H1 deficit cost failed")
    return {
        "branch": branch,
        "dephasing_threshold": varth,
        "physical_first_duhamel_required": req,
        "clean_physical_quadratic_cost": cost,
    }


def full_curvature_channel(I_full: float, I_h3: float, I_hook: float) -> str:
    if min(I_full, I_h3, I_hook) < 0:
        raise ValueError("nonnegative impulses required")
    if I_full > math.sqrt(6) * I_h3 + I_hook + 1e-11 * max(1.0, I_full):
        raise ValueError("impulses violate curvature split upper bound")
    if I_h3 >= I_full / (2 * math.sqrt(6)) - 1e-13:
        return "H3"
    if I_hook + 1e-12 < I_full / 2:
        raise AssertionError("full curvature channel dichotomy failed")
    return "H1_hook"


def full_mild_aspect_quadratic_cost(I_full: float) -> float:
    if I_full < 0:
        raise ValueError("nonnegative impulse required")
    return float(FULL_MILD_COST) * I_full * I_full


def arb_physical_transport_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    # K<=1/30. Interaction forcing L1 >= exp(-K) I/10 > I/11.
    if not ((-arb(1) / 30).exp() / 10 > arb(1) / 11):
        raise AssertionError("interaction variation conditioning failed")
    # Physical first Duhamel >= exp(-2K)I/20; square > I^2/480 iff exp(-4K)>5/6.
    if not ((-arb(2) / 15).exp() > arb(5) / 6):
        raise AssertionError("physical Duhamel conditioning failed")
    return {
        "base_role_action_budget": "K<=1/30",
        "interaction_variation_clean": "I_1/(11 T)",
        "physical_first_duhamel_square": "I_1^2/480",
        "post_feedback_total_square": "I_1^2/1920",
        "one_role_square": "I_1^2/5760",
        "single_role_deficit_before_rescue": "I_1^2/92160",
        "pair_or_deficit": "I_1^2/184320",
        "full_mild_pair_or_deficit": "I_B^2/737280",
        "status": "CERTIFIED_PHYSICAL_LOW_STRAIN_CONDITIONING",
    }


def exact_constant_certificate() -> dict[str, str]:
    # first physical Duhamel^2 1/480; half feedback -> 1/1920;
    # select one of 3 roles -> 1/5760; odd-role deficit /16 -> 1/92160;
    # split against pair rescue -> 1/184320.
    if Fraction(1, 480) * Fraction(1, 4) * Fraction(1, 3) * Fraction(1, 16) * Fraction(1, 2) != PHYSICAL_H1_COST:
        raise AssertionError("physical H1 constant mismatch")
    if PHYSICAL_H1_COST * Fraction(1, 4) != FULL_MILD_COST:
        raise AssertionError("full mild cost mismatch")
    if not Fraction(1, 32768) > FULL_MILD_COST:
        raise AssertionError("H3 branch should dominate full physical constant")
    return {
        "physical_H1_pair_or_deficit": "I_1^2/184320",
        "full_mild_pair_or_deficit": "I_B^2/737280",
        "status": "EXACT_AFTER_PHYSICAL_ROLE_SPLIT",
    }


@dataclass(frozen=True)
class H1NoEscapeStress:
    samples: int
    minimum_H1_pair_margin: float
    minimum_H1_deficit_margin: float
    minimum_full_channel_margin: float
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260807) -> H1NoEscapeStress:
    rng = np.random.default_rng(seed)
    mp = md = mf = float("inf")
    counts: dict[str, int] = {}
    for _ in range(samples):
        I = float(10 ** rng.uniform(-4, -0.5))
        T = float(10 ** rng.uniform(-3, 0))
        req = math.sqrt(clean_first_impulse_square_lower(I))
        cost = h1_physical_quadratic_cost(I)
        mode = int(rng.integers(0, 5))
        if mode == 0:
            J = (1 + rng.random()) * h1_dephasing_threshold(I, T)
            d1 = fb = sigma = rescue = deficit = 0.0
        else:
            J = float(rng.uniform(0, 0.99)) * h1_dephasing_threshold(I, T)
            d1 = (1 + rng.random()) * req
            if mode == 1:
                fb = float(rng.uniform(0.5, 1.2)) * d1
                sigma = rescue = deficit = 0.0
            else:
                fb = float(rng.uniform(0, 0.49)) * d1
                total = d1 - fb
                one = total / math.sqrt(3)
                if mode == 2:
                    sigma = max(1 / 80, one) * (1 + rng.random())
                    rescue = deficit = 0.0
                else:
                    upper = 0.99 / 80
                    if one >= upper:
                        I = 1e-4
                        req = math.sqrt(clean_first_impulse_square_lower(I))
                        d1 = 1.2 * req
                        fb = 0.1 * d1
                        total = d1 - fb
                        one = total / math.sqrt(3)
                        cost = h1_physical_quadratic_cost(I)
                        J = 0.1 * h1_dephasing_threshold(I, T)
                    sigma = float(rng.uniform(max(one, 1e-12), upper))
                    d0 = sigma * sigma / 16
                    if mode == 3:
                        rescue = float(rng.uniform(0.5, 1)) * d0
                        deficit = 0.0
                        mp = min(mp, rescue - cost)
                    else:
                        rescue = float(rng.uniform(0, 0.49)) * d0
                        deficit = max(0.5 * d0, cost) * (1 + rng.random())
                        md = min(md, deficit - cost)
        out = classify_h1_mild_no_escape(I, T, J, d1, fb, sigma, rescue, deficit)
        counts[out["branch"]] = counts.get(out["branch"], 0) + 1
        Ifull = float(10 ** rng.uniform(-4, 0))
        Ih3 = float(rng.uniform(0, Ifull / math.sqrt(6)))
        Ihook = max(0.0, Ifull - math.sqrt(6) * Ih3) + float(rng.uniform(0, 0.2)) * Ifull
        ch = full_curvature_channel(Ifull, Ih3, Ihook)
        margin = Ih3 - Ifull / (2 * math.sqrt(6)) if ch == "H3" else Ihook - Ifull / 2
        mf = min(mf, margin)
        if margin < -2e-12:
            raise AssertionError("full channel stress failed")
    return H1NoEscapeStress(samples, mp, md, mf, counts)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-h1-swirl-no-escape"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    arb = arb_physical_transport_certificate()
    exact = exact_constant_certificate()
    out = stress(args.samples)
    (args.outdir / "h1_swirl_no_escape.json").write_text(
        json.dumps({"arb_certificate": arb, "exact_certificate": exact, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Physical mild-aspect H1/swirl local no-escape theorem

Status: **{arb['status']}** + **{exact['status']}**.

Two corrections are built into this final physical theorem. First, `Q_pol/2` is an auxiliary relative-coordinate forcing norm; the three physical Young roles only satisfy `sum_i||F_i^H1||^2 >= Q_pol/4 >= ||B_hook||^2/100`. Second, the base polarization propagators are non-unitary. On the existing low-strain lifetime branch each physical role has action budget `K<=1/30`, so pullback and pushforward singular values are controlled by `exp(+-K)` rather than treated as isometries.

Let `I1=int||B_hook||dt`. The conditioned interaction-picture variation theorem gives

`J1 >= I1/(11 T)`

or a **physical three-role** first-Duhamel daughter with

`delta1^2 >= I1^2/480`.

If nonlinear physical feedback is less than half, total surviving daughter energy is at least `I1^2/1920`; one of the three roles therefore has energy at least `I1^2/5760`. Below critical sideband size `1/80`, odd-Hermite convexity and the pair-rescue split yield

`net transfer deficit >= I1^2/184320`

or

`pair-sideband rescue >= I1^2/184320`.

Combining with the H3 branch and `I_B<=sqrt(6)I3+I1` gives the clean physical mild-aspect full-curvature cost

`pair rescue or transfer deficit >= I_B^2/737280`

outside dephasing/source, nonlinear-feedback and large-daughter branches.

This supersedes the idealized relative-coordinate/isometric-pullback constants `1/25600` and `1/102400` from run `31195130386`. The pointwise mild-aspect bridge `Q_pol>=1/25||B_hook||^2` is unchanged. H1 covariant dephasing source calculus remains open, and high-aspect grains remain ancestry/reuse rather than an aspect defect.

Stress: `{out.samples}`
- branch counts: `{out.branch_counts}`
- minimum H1 pair-cost margin: `{out.minimum_H1_pair_margin:.3e}`
- minimum H1 deficit-cost margin: `{out.minimum_H1_deficit_margin:.3e}`
- minimum full-channel margin: `{out.minimum_full_channel_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
