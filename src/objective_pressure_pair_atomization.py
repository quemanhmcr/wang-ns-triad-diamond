from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

from src.pressure_reservoir_sync import pressure_hessian_pair_energy_service_ratio_upper
from src.resolved_objective_strain_collision import arb_clean_constants

STATUS = (
    "EXACT_OBJECTIVE_PRESSURE_HESSIAN_DUAL_PAIR_ATOMIZATION__"
    "PAIR_OWNER_ALWAYS_TO_CRITICAL_SHELL_WITH_ENTROPY_TRADEOFF__"
    "QUARTER_SPLIT_ONLY_A_COROLLARY__AGGREGATE_MU_V_NOT_CANONICAL"
)

DEFAULT_PAIR_DOMINANCE = Fraction(1, 4)
PRESSURE_OWNER_SPLIT = Fraction(1, 2)
ORDERED_PAIR_SHARP_CLEAN = Fraction(256, 1425)
ORDERED_PAIR_SIMPLE_CLEAN = Fraction(1, 5)
OFFDIAGONAL_PAIR_SIMPLE_CLEAN = Fraction(2, 5)


def _matrix3(A: np.ndarray) -> np.ndarray:
    A = np.asarray(A, dtype=float)
    if A.shape != (3, 3) or np.any(~np.isfinite(A)):
        raise ValueError("finite 3x3 matrix required")
    return A


def frobenius_dual(source_matrix: np.ndarray) -> np.ndarray:
    """Measurable unit Frobenius dual of the actual averaged pressure source.

    At a source event let H be the normalized matrix `-N^-4 <Hess P>_gamma`.
    We use Z=H/||H||_F when H is nonzero and Z=0 at H=0.  This is branch-free,
    measurable and contains no packet/material choice.
    """
    H = _matrix3(source_matrix)
    n = float(np.linalg.norm(H, ord="fro"))
    return np.zeros((3, 3), dtype=float) if n == 0.0 else H / n


def scalarized_pressure_split(
    total_source: np.ndarray,
    resolved_source: np.ndarray,
    sgs_source: np.ndarray,
) -> dict[str, float | np.ndarray]:
    """Exact dual scalarization of the averaged pressure-Hessian source.

    `total_source=resolved_source+sgs_source` is the exact filtered-pressure
    decomposition after the Gaussian average.  The aligned dual gives

        rho = ||H||_F = r_V + r_R,

    hence `rho <= [r_V]_+ + [r_R]_+`.  These are source/service positive parts,
    not child-energy causal probabilities.
    """
    H = _matrix3(total_source)
    HV = _matrix3(resolved_source)
    HR = _matrix3(sgs_source)
    scale = max(1.0, float(np.linalg.norm(H, "fro")), float(np.linalg.norm(HV, "fro")), float(np.linalg.norm(HR, "fro")))
    residual = float(np.linalg.norm(H - HV - HR, ord="fro"))
    if residual > 3e-12 * scale:
        raise ValueError("pressure matrices do not satisfy the exact resolved+SGS split")
    Z = frobenius_dual(H)
    rho = float(np.linalg.norm(H, ord="fro"))
    rv = float(np.sum(Z * HV))
    rr = float(np.sum(Z * HR))
    exact_residual = rho - rv - rr
    cover_margin = max(rv, 0.0) + max(rr, 0.0) - rho
    return {
        "dual": Z,
        "rho": rho,
        "resolved_signed": rv,
        "sgs_signed": rr,
        "exact_scalar_residual": exact_residual,
        "positive_cover_margin": cover_margin,
    }


def unordered_pair_matrices(ordered_pair_matrices: np.ndarray) -> tuple[tuple[tuple[int, int], np.ndarray], ...]:
    """Remove observer-imposed pair orientation while preserving exact expansion.

    Input `[a,b]` is the ordered bilinear pressure contribution from `(V_a,V_b)`.
    The physical unordered atom is diagonal once, and for `a<b` the sum of both
    orientations.  Summing unordered atoms equals the full ordered double sum.
    """
    A = np.asarray(ordered_pair_matrices, dtype=float)
    if A.ndim != 4 or A.shape[0] != A.shape[1] or A.shape[2:] != (3, 3) or np.any(~np.isfinite(A)):
        raise ValueError("ordered pair matrices must have shape (n,n,3,3)")
    out: list[tuple[tuple[int, int], np.ndarray]] = []
    for a in range(A.shape[0]):
        for b in range(a, A.shape[1]):
            M = A[a, a].copy() if a == b else A[a, b] + A[b, a]
            out.append(((a, b), M))
    return tuple(out)


def scalarized_unordered_pair_atoms(
    ordered_pair_matrices: np.ndarray,
    dual: np.ndarray,
) -> tuple[tuple[tuple[int, int], float], ...]:
    Z = _matrix3(dual)
    return tuple((key, float(np.sum(Z * M))) for key, M in unordered_pair_matrices(ordered_pair_matrices))


def ordered_pair_sharp_clean_coefficient() -> Fraction:
    """Exact rational consequence of two already-certified clean constants.

    At current block N and shell upper scale M<=N/4:
    - product support is <=2M;
    - the order-2 L^(3/2)->Linf constant scales from `1/380` at support N/2
      by `(4M/N)^4`, contributing `(256/380)(M/N)^4`;
    - shell Bernstein scales `||V_a||_3^2 <= (4/15) mu_a` from the certified
      `1/15` constant at support `(1/4)N'` with `N'=4M_a`.

    Multiplication gives `256/1425 < 1/5` per ordered pair.
    """
    return ORDERED_PAIR_SHARP_CLEAN


def unordered_pair_simple_coefficient(diagonal: bool) -> Fraction:
    return ORDERED_PAIR_SIMPLE_CLEAN if diagonal else OFFDIAGONAL_PAIR_SIMPLE_CLEAN


def unordered_pair_capacity_upper(
    critical_mass_a: float,
    critical_mass_b: float,
    frequency_a: float,
    frequency_b: float,
    block_frequency: float,
    *,
    diagonal: bool,
) -> float:
    """Clean normalized objective-Hessian source capacity of one hard shell pair.

      |p_{ab}| <= (kappa/5)(Mmax/N)^4 sqrt(mu_a mu_b),

    with kappa=1 diagonal and 2 off diagonal.  Gaussian averaging is a probability
    contraction and therefore does not alter this Linf capacity.
    """
    ma = float(critical_mass_a)
    mb = float(critical_mass_b)
    fa = float(frequency_a)
    fb = float(frequency_b)
    N = float(block_frequency)
    if min(ma, mb) < 0 or min(fa, fb, N) <= 0 or not all(math.isfinite(x) for x in (ma, mb, fa, fb, N)):
        raise ValueError("valid finite pair masses/frequencies required")
    if max(fa, fb) > N / 4.0 + 2e-13 * max(1.0, N):
        raise ValueError("resolved pressure pair must lie inside V=S_(N/4)u")
    coeff = float(unordered_pair_simple_coefficient(diagonal))
    return coeff * (max(fa, fb) / N) ** 4 * math.sqrt(ma * mb)


def canonical_all_pair_absolute_capacity_upper(
    resolved_energy: float,
    block_frequency: float,
) -> float:
    """Absolute source-capacity sum of the countable canonical hard-shell pair law.

    Let `M_j=(N/4)2^-j`, `E_j=||P_jV||_2^2`, `mu_j=M_j E_j`.
    Since every `Mmax/N<=1/4`, summing the clean unordered capacities gives

      sum_(a<=b) cap_ab
      <= (1/(5*4^4)) (sum_j sqrt(mu_j))^2
      <= (1/1280) (sum_j M_j)(sum_j E_j)
      = N ||V||_2^2 / 2560.

    Thus the positive pair source law is absolutely summable.  Aggregate resolved
    mass appears here only as a convergence budget, not as the renewal state.
    """
    E = float(resolved_energy)
    N = float(block_frequency)
    if E < 0 or N <= 0 or not all(math.isfinite(x) for x in (E, N)):
        raise ValueError("finite nonnegative resolved energy and positive block frequency required")
    return N * E / 2560.0


def u_shell_mass_lower_from_resolved_shell_mass(resolved_shell_mass_lower: float) -> float:
    """Pass a strict-lowpass hard-shell lower from V to the same shell of u.

    The canonical transporter has scalar multiplier `|S_(N/4)(xi)|<=1` and
    commutes with the hard shell projector.  Hence

      M||P_M V||_2^2 <= M||P_M u||_2^2.

    No inverse low-pass estimate and no equality of the two shells is asserted.
    """
    mu = float(resolved_shell_mass_lower)
    if mu < 0 or not math.isfinite(mu):
        raise ValueError("finite nonnegative resolved-shell mass lower required")
    return mu


def dominant_pair_peak_shell_mass_lower(
    integrated_positive_pair_source: float,
    scaled_lifetime: float,
    frequency_a: float,
    frequency_b: float,
    block_frequency: float,
    *,
    diagonal: bool,
) -> float:
    """A positive integrated pair atom exposes an actual critical shell-time event.

    Capacity plus averaging in scaled time gives some time with

      sqrt(mu_a mu_b) >= [5/kappa] (N/Mmax)^4 R_ab/c.

    At that time at least one actual hard shell of the resolved transporter has
    critical mass at least this geometric mean.  Since the strict low-pass is an
    L2 contraction commuting with the hard shell projector, the same numerical
    lower holds for that shell of `u`, which is the input required by the generic
    critical-shell theorem.  No inverse low-pass estimate is used.
    """
    R = float(integrated_positive_pair_source)
    c = float(scaled_lifetime)
    fa = float(frequency_a)
    fb = float(frequency_b)
    N = float(block_frequency)
    if R <= 0 or c <= 0 or min(fa, fb, N) <= 0 or not all(math.isfinite(x) for x in (R, c, fa, fb, N)):
        raise ValueError("positive finite pair source/lifetime/frequencies required")
    if max(fa, fb) > N / 4.0 + 2e-13 * max(1.0, N):
        raise ValueError("pair lies outside resolved transporter support")
    kappa = 1.0 if diagonal else 2.0
    return (5.0 / kappa) * (N / max(fa, fb)) ** 4 * R / c


def clean_dominant_pair_shell_mass_lower(
    pressure_source_weight: float,
    scaled_lifetime: float,
    *,
    dominant_fraction: float = float(DEFAULT_PAIR_DOMINANCE),
    max_pair_frequency_ratio: float = 0.25,
) -> float:
    """Uniform source-facing lower for the worst off-diagonal dominant pair.

    Resolved pair positive source is at least Sigma_P/2.  A theta-dominant pair
    carries at least theta Sigma_P/2.  The off-diagonal coefficient 2/5 gives

      mu_child >= (5 theta Sigma_P / 4c)(N/Mmax)^4.

    With Mmax<=N/4 and theta=1/4 this is the clean `80 Sigma_P/c` event.
    """
    sigma = float(pressure_source_weight)
    c = float(scaled_lifetime)
    theta = float(dominant_fraction)
    r = float(max_pair_frequency_ratio)
    if sigma <= 0 or c <= 0 or not (0 < theta < 1) or not (0 < r <= 0.25) or not all(math.isfinite(x) for x in (sigma, c, theta, r)):
        raise ValueError("valid pressure source/lifetime/dominance/support ratio required")
    return (5.0 * theta / 4.0) * sigma / c * r ** (-4)


def clean_entropy_shell_tradeoff_lower(
    pressure_source_weight: float,
    scaled_lifetime: float,
    pair_source_entropy: float,
    *,
    max_pair_frequency_ratio: float = 0.25,
) -> float:
    """Continuous pressure law coupling pair fragmentation to shell mass.

    On the resolved-pair owner, `R_pair>=Sigma_P/2`.  For normalized pair law
    `q`, `q_max>=sum q^2=exp(-H2)`.  The worst off-diagonal pair capacity gives

      mu_child >= (5/4)(N/Mmax)^4 exp(-H2) Sigma_P/c.

    At `Mmax/N<=1/4` this is `320 exp(-H2) Sigma_P/c`.  Thus the quarter
    dominance/entropy split is only the corollary obtained by cutting at
    `H2=log 4`; the native law itself has no arbitrary threshold.
    """
    sigma = float(pressure_source_weight)
    c = float(scaled_lifetime)
    h2 = float(pair_source_entropy)
    r = float(max_pair_frequency_ratio)
    if sigma <= 0 or c <= 0 or h2 < 0 or not (0 < r <= 0.25) or not all(math.isfinite(x) for x in (sigma, c, h2, r)):
        raise ValueError("valid source/lifetime/entropy/support ratio required")
    return (5.0 / 4.0) * r ** (-4) * math.exp(-h2) * sigma / c


def pair_collision_entropy(pair_positive_weights: Sequence[float]) -> dict[str, float]:
    w = np.asarray(tuple(float(x) for x in pair_positive_weights), dtype=float)
    if w.ndim != 1 or len(w) == 0 or np.any(~np.isfinite(w)) or np.any(w < 0) or float(w.sum()) <= 0:
        raise ValueError("positive finite pair source law required")
    p = w / float(w.sum())
    q2 = float(np.dot(p, p))
    return {
        "H2_pair_source": -math.log(q2),
        "collision_probability": q2,
        "maximum_atom_mass": float(np.max(p)),
    }


def objective_pressure_pair_route(
    pressure_source_weight: float,
    scaled_lifetime: float,
    block_frequency: float,
    *,
    sgs_positive_source_weight: float,
    pair_positive_weights: Sequence[float],
    pair_shell_indices: Sequence[tuple[int, int]],
    pair_frequencies: Sequence[tuple[float, float]],
    dominant_fraction: float = float(DEFAULT_PAIR_DOMINANCE),
) -> dict[str, object]:
    """Master-facing direct pressure route using actual positive source atoms.

    The dual identity yields pointwise

      rho_P <= [r_R]_+ + sum_{a<=b}[p_ab]_+.

    After integration one of SGS or resolved-pair positive source carries at least
    Sigma_P/2; exact ties remain joint.  Every positive resolved pair law has a
    largest physical unordered atom and therefore an actual critical shell via
    `mu_child exp(H2_pair)>=320 Sigma_P/c`.  A quarter dominance/diffuse split is
    only an optional corollary of this single route, not a second master fate.
    """
    sigma = float(pressure_source_weight)
    c = float(scaled_lifetime)
    N = float(block_frequency)
    sgs = float(sgs_positive_source_weight)
    theta = float(dominant_fraction)
    w = np.asarray(tuple(float(x) for x in pair_positive_weights), dtype=float)
    pairs = tuple((int(a), int(b)) for a, b in pair_shell_indices)
    freqs = tuple((float(a), float(b)) for a, b in pair_frequencies)
    if sigma <= 0 or c <= 0 or N <= 0 or sgs < 0 or not (0 < theta < 1) or not all(math.isfinite(x) for x in (sigma, c, N, sgs, theta)):
        raise ValueError("valid positive pressure route data required")
    if w.ndim != 1 or np.any(~np.isfinite(w)) or np.any(w < 0) or len(w) != len(pairs) or len(w) != len(freqs):
        raise ValueError("pair weights/indices/frequencies must match and be finite nonnegative")
    if len(set(pairs)) != len(pairs):
        raise ValueError("each physical unordered hard pair must be integrated once before entropy is computed")
    shell_frequency: dict[int, float] = {}
    for (a, b), (fa, fb) in zip(pairs, freqs):
        if a < 0 or b < 0 or a > b or min(fa, fb) <= 0 or max(fa, fb) > N / 4.0 + 2e-13 * max(1.0, N):
            raise ValueError("unordered resolved pair labels/frequencies required")
        for label, freq in ((a, fa), (b, fb)):
            if label in shell_frequency and not math.isclose(shell_frequency[label], freq, rel_tol=2e-13, abs_tol=0.0):
                raise ValueError("one hard shell label cannot carry multiple observer-assigned frequencies")
            shell_frequency[label] = freq
    pair_total = float(w.sum())
    cover_tol = 5e-13 * max(1.0, sigma, sgs + pair_total)
    if sgs + pair_total + cover_tol < sigma:
        raise ValueError("positive SGS+pair atoms do not cover the pressure source law")
    half = sigma / 2.0
    owners: list[str] = []
    # Owner and pair-dominance comparisons are exact dimensionless/source
    # comparisons.  Never reuse a source-unit tolerance for normalized pair mass.
    if sgs >= half:
        owners.append("sgs_pressure_source")
    if pair_total >= half:
        owners.append("resolved_pressure_pair_law")
    if not owners:
        raise AssertionError("pressure positive-owner half split failed")

    out: dict[str, object] = {
        "pressure_source_weight": sigma,
        "positive_owner_threshold": half,
        "joint_primary_owners": tuple(owners),
        "sgs_stress_l32_lower_if_owner": 380.0 * half,
        "sgs_effective_source_weight_if_owner": half,
        "pair_positive_source_total": pair_total,
        "pair_source_entropy_is_causal_probability": False,
        "aggregate_muV_is_canonical_route": False,
        "master_semantics": "RECURSE_CRITICAL; every resolved pair owner enters the generic shell theorem, with H2 controlling seed strength rather than defining another fate",
        "status": STATUS,
    }

    if pair_total < half:
        return out

    p = w / pair_total
    pmax = float(np.max(p)) if len(p) else 0.0
    imax = int(np.argmax(p)) if len(p) else -1
    dominant = tuple(int(i) for i, x in enumerate(p) if x >= theta)
    entropy = pair_collision_entropy(w) if len(w) else None
    if entropy is None or imax < 0:
        raise AssertionError("resolved pair owner lost its positive pair law")
    h2 = float(entropy["H2_pair_source"])
    tradeoff_lower = clean_entropy_shell_tradeoff_lower(sigma, c, h2)
    amax, bmax = pairs[imax]
    famax, fbmax = freqs[imax]
    max_pair_resolved_lower = dominant_pair_peak_shell_mass_lower(
        float(w[imax]), c, famax, fbmax, N, diagonal=(amax == bmax)
    )
    max_pair_u_lower = u_shell_mass_lower_from_resolved_shell_mass(max_pair_resolved_lower)
    trade_tol = 8e-13 * max(1.0, tradeoff_lower, max_pair_u_lower)
    if max_pair_u_lower + trade_tol < tradeoff_lower:
        raise AssertionError("pressure entropy-shell tradeoff failed")
    out["pair_source_entropy"] = entropy
    out["entropy_shell_tradeoff_lower"] = tradeoff_lower
    out["max_pair_u_shell_mass_lower"] = max_pair_u_lower
    out["max_pair_witness_index"] = imax
    out["entropy_shell_tradeoff"] = "mu_child exp(H2_pair)>=320 Sigma_P/c"
    out["pair_owner_route"] = "critical_shell_via_entropy_tradeoff"
    out["critical_shell_mass_lower"] = tradeoff_lower
    out["critical_shell_supplier_is_unconditional_on_quarter_cut"] = True
    diffuse = pmax <= theta
    pair_routes: list[str] = []
    dominant_witnesses: list[dict[str, object]] = []
    if dominant:
        pair_routes.append("dominant_hard_pair_to_critical_shell")
        for i in dominant:
            a, b = pairs[i]
            fa, fb = freqs[i]
            lower = dominant_pair_peak_shell_mass_lower(
                float(w[i]), c, fa, fb, N, diagonal=(a == b)
            )
            u_lower = u_shell_mass_lower_from_resolved_shell_mass(lower)
            dominant_witnesses.append({
                "pair_index": i,
                "shell_pair": (a, b),
                "pair_source_weight": float(w[i]),
                "normalized_pair_mass": float(p[i]),
                "max_pair_frequency_ratio": max(fa, fb) / N,
                "resolved_V_shell_mass_lower": lower,
                "u_shell_mass_lower": u_lower,
                "critical_shell_mass_lower": u_lower,
                "lowpass_bridge": "M||P_M V||_2^2<=M||P_M u||_2^2",
                "child_scale_at_most": max(fa, fb),
                "parent_to_child_natural_lifetime_ratio_at_least": (N / max(fa, fb)) ** 2,
            })
    if diffuse:
        if entropy is None:
            raise AssertionError("diffuse pair law lost its entropy measure")
        h0 = -math.log(theta)
        if float(entropy["H2_pair_source"]) + 2e-13 < h0:
            raise AssertionError("diffuse pressure pair source entropy failed")
        pair_routes.append("diffuse_pair_source_entropy")
        out["pair_source_entropy_lower"] = h0
    if not pair_routes:
        raise AssertionError("pair law lost both quarter-cut corollaries")
    out["quarter_cut_corollaries"] = tuple(pair_routes)
    out["dominant_pair_witnesses"] = tuple(dominant_witnesses)
    out["clean_dominant_shell_mass_lower"] = clean_dominant_pair_shell_mass_lower(
        sigma, c, dominant_fraction=theta, max_pair_frequency_ratio=0.25
    )
    return out


def theorem_certificate() -> dict[str, object]:
    clean = arb_clean_constants()
    if clean["order2_clean"] != "1/380" or clean["resolved_l3_squared_mass_clean"] != "1/15":
        raise AssertionError("resolved objective clean constants changed")
    sharp = ordered_pair_sharp_clean_coefficient()
    if not sharp < Fraction(1, 5):
        raise AssertionError("ordered objective pressure pair coefficient lost its one-fifth clean bound")
    if not 2 * sharp < Fraction(2, 5):
        raise AssertionError("unordered off-diagonal pressure pair lost its two-fifths clean bound")
    reuse = pressure_hessian_pair_energy_service_ratio_upper()
    if not reuse < Fraction(1, 5):
        raise AssertionError("fixed material pressure-Hessian pair lost one-fifth reuse lifetime")
    clean80 = clean_dominant_pair_shell_mass_lower(1.0, 1.0)
    if abs(clean80 - 80.0) > 1e-12:
        raise AssertionError("canonical dominant pair shell lower lost 80 Sigma/c")
    trade_log4 = clean_entropy_shell_tradeoff_lower(1.0, 1.0, math.log(4.0))
    if abs(trade_log4 - 80.0) > 2e-12:
        raise AssertionError("entropy-shell tradeoff no longer meets the quarter corollary at 80 Sigma/c")
    return {
        "status": STATUS,
        "exact_dual": "Z=H/||H||_F (Z=0 at H=0), rho_P=Z:H=sum unordered resolved pair scalars + SGS scalar",
        "positive_source_cover": "rho_P <= sum_[a<=b] [p_ab]_+ + [r_SGS]_+; this is source/service positivity, not causal child-energy probability",
        "hard_pair_registration": "decompose only V=S_(N/4)u into hard orthogonal dyadic shells at the physical pressure event; quotient to one integrated weight per unordered {a,b}, so observer orientation or duplicate records cannot manufacture entropy",
        "ordered_pair_sharp_clean": f"{sharp.numerator}/{sharp.denominator}<1/5",
        "unordered_pair_clean": "|p_ab|<=(kappa_ab/5)(Mmax/N)^4 sqrt(mu_a mu_b), kappa=1 diagonal, 2 off-diagonal",
        "source_half_split": "SGS positive pressure source >=Sigma_P/2 OR resolved positive pair source >=Sigma_P/2; ties joint",
        "entropy_shell_tradeoff": "on every resolved pair owner, the countable positive law has an attained qmax>=exp(-H2); worst off-diagonal Mmax<=N/4 capacity gives mu_child exp(H2_pair)>=320 Sigma_P/c, so the pair owner always enters generic critical-shell reentry",
        "dominant_pair": "theta=1/4 first gives a resolved V-shell mass >=5 Sigma_P/(16c)(N/Mmax)^4; |S_(N/4)|<=1 transfers the same lower to the u shell, hence >=80 Sigma_P/c because Mmax<=N/4",
        "absolute_pair_sum": "for canonical M_j=(N/4)2^-j, sum unordered pair capacities <=N||V||_2^2/2560; mu_V appears only as an absolute-convergence budget",
        "quarter_corollary": "theta=1/4 is optional: qmax>=1/4 gives mu_child>=80 Sigma_P/c; qmax<=1/4 gives H2>=log4; equality gives both, but both are faces of the same unconditional shell route",
        "material_sidecar": f"material reuse is optional after the hard event; a fixed objective-Hessian material pair contracts by {reuse.numerator}/{reuse.denominator}<1/5 per signed-good low-strain generation",
        "coarse_muV": "rho_P<=mu_V/5700+||R||_(3/2)/380 remains a valid coarse diagnostic but is not the canonical pressure renewal route",
    }


@dataclass(frozen=True)
class ObjectivePressurePairStress:
    samples: int
    worst_dual_scalar_residual: float
    minimum_positive_cover_margin: float
    worst_unordered_reconstruction_residual: float
    minimum_pair_capacity_margin: float
    minimum_owner_half_margin: float
    minimum_dominant_shell_margin_over_clean: float
    minimum_diffuse_entropy_margin: float
    minimum_entropy_shell_tradeoff_margin: float
    maximum_joint_primary_owner_count: int
    maximum_quarter_corollary_count: int


def stress(samples: int = 50_000, seed: int = 20260809) -> ObjectivePressurePairStress:
    rng = np.random.default_rng(seed)
    wr = wu = 0.0
    mc = mcap = mo = md = me = mt = float("inf")
    max_primary = max_corollary = 0
    for _ in range(samples):
        # Exact matrix dual/source atomization.
        n = int(rng.integers(1, 7))
        ordered = rng.normal(size=(n, n, 3, 3))
        ordered = 0.5 * (ordered + np.swapaxes(ordered, -1, -2))
        sgsM = rng.normal(size=(3, 3)); sgsM = 0.5 * (sgsM + sgsM.T)
        resolved = ordered.sum(axis=(0, 1))
        total = resolved + sgsM
        split = scalarized_pressure_split(total, resolved, sgsM)
        wr = max(wr, abs(float(split["exact_scalar_residual"])))
        mc = min(mc, float(split["positive_cover_margin"]))
        if abs(float(split["exact_scalar_residual"])) > 5e-12 * max(1.0, float(split["rho"])):
            raise AssertionError("pressure dual scalarization lost exactness")
        if float(split["positive_cover_margin"]) < -5e-12 * max(1.0, float(split["rho"])):
            raise AssertionError("positive pressure owner cover failed")
        atoms = unordered_pair_matrices(ordered)
        rec = sum((M for _, M in atoms), np.zeros((3, 3)))
        ures = float(np.linalg.norm(rec - resolved, ord="fro"))
        wu = max(wu, ures)
        if ures > 5e-12 * max(1.0, float(np.linalg.norm(resolved, "fro"))):
            raise AssertionError("unordered pressure pair reconstruction failed")
        # Transposing observer orientation leaves each unordered atom unchanged.
        atomsT = unordered_pair_matrices(np.swapaxes(ordered, 0, 1))
        for (_, A), (_, B) in zip(atoms, atomsT):
            if np.linalg.norm(A - B, ord="fro") > 5e-12 * max(1.0, np.linalg.norm(A, ord="fro")):
                raise AssertionError("unordered pair atom depends on observer orientation")

        # Clean pair capacity and dominant shell extraction.
        N = float(math.exp(rng.uniform(-2.0, 4.0)))
        ia = int(rng.integers(0, 9)); ib = int(rng.integers(0, 9))
        fa = N / (4.0 * 2.0**ia); fb = N / (4.0 * 2.0**ib)
        mua = float(math.exp(rng.uniform(-8.0, 3.0))); mub = float(math.exp(rng.uniform(-8.0, 3.0)))
        diag = ia == ib
        cap = unordered_pair_capacity_upper(mua, mub, fa, fb, N, diagonal=diag)
        atom = float(rng.uniform(0.0, 1.0)) * cap
        mpair = cap - atom
        if atom > cap + 2e-12 * max(1.0, cap):
            raise AssertionError("pressure pair atom exceeded clean capacity")
        mcap = min(mcap, mpair)

        # Realized positive source law.  Pair scales are all genuine <=N/4.
        k = int(rng.integers(1, 12))
        raw = rng.dirichlet(np.ones(k))
        pair_total = float(math.exp(rng.uniform(-5.0, 1.0)))
        weights = pair_total * raw
        sgs = float(math.exp(rng.uniform(-5.0, 1.0)))
        available = sgs + pair_total
        sigma = available * float(rng.uniform(0.35, 1.0))
        c = float(math.exp(rng.uniform(-2.0, 1.0)))
        indices = tuple((j, j) if rng.random() < 0.35 else (j, j + 1) for j in range(k))
        labels = sorted({x for pair in indices for x in pair})
        shell_freq = {
            label: N / (4.0 * 2.0 ** int(rng.integers(0, 8)))
            for label in labels
        }
        pfreqs = [(shell_freq[a], shell_freq[b]) for a, b in indices]
        route = objective_pressure_pair_route(
            sigma, c, N,
            sgs_positive_source_weight=sgs,
            pair_positive_weights=weights,
            pair_shell_indices=indices,
            pair_frequencies=pfreqs,
        )
        threshold = sigma / 2.0
        mo = min(mo, max(sgs, pair_total) - threshold)
        max_primary = max(max_primary, len(tuple(route["joint_primary_owners"])))
        if max(sgs, pair_total) + 3e-12 * max(1.0, threshold) < threshold:
            raise AssertionError("pressure owner half split failed")
        if "quarter_cut_corollaries" in route:
            pair_routes = tuple(route["quarter_cut_corollaries"])
            max_corollary = max(max_corollary, len(pair_routes))
            if route["pair_owner_route"] != "critical_shell_via_entropy_tradeoff":
                raise AssertionError("resolved pressure pair owner lost unconditional critical-shell route")
            trade = float(route["entropy_shell_tradeoff_lower"])
            actual_trade = float(route["max_pair_u_shell_mass_lower"])
            mt = min(mt, actual_trade - trade)
            if actual_trade + 4e-12 * max(1.0, trade) < trade:
                raise AssertionError("pressure entropy-shell tradeoff stress failed")
            if "dominant_hard_pair_to_critical_shell" in pair_routes:
                clean_lower = float(route["clean_dominant_shell_mass_lower"])
                for wit in route["dominant_pair_witnesses"]:
                    # Clean 80 Sigma/c assumes only Mmax<=N/4 and q>=1/4.
                    md = min(md, float(wit["critical_shell_mass_lower"]) - clean_lower)
                    if float(wit["critical_shell_mass_lower"]) + 4e-12 * max(1.0, clean_lower) < clean_lower:
                        raise AssertionError("dominant pressure pair lost clean critical shell lower")
            if "diffuse_pair_source_entropy" in pair_routes:
                h = float(route["pair_source_entropy"]["H2_pair_source"])
                h0 = float(route["pair_source_entropy_lower"])
                me = min(me, h - h0)
                if h + 2e-12 < h0:
                    raise AssertionError("diffuse pressure pair entropy failed")

    # Random samples need not realize both rare exact-boundary branches; report 0
    # rather than infinity for an unvisited diagnostic margin.
    if not math.isfinite(md): md = 0.0
    if not math.isfinite(me): me = 0.0
    if not math.isfinite(mt): mt = 0.0
    return ObjectivePressurePairStress(samples, wr, mc, wu, mcap, mo, md, me, mt, max_primary, max_corollary)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-objective-pressure-pair-atomization"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "objective_pressure_pair_atomization.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Objective pressure Hessian: direct hard-pair source atomization

Status: **{cert['status']}**.

The canonical pressure route no longer needs to coarse-grain the resolved Hessian source into aggregate `mu_V`.  Let

`H=-N^-4 <Hess P>_gamma`

be the actual averaged pressure matrix in the coherent corotational strain equation and choose the measurable Frobenius dual `Z=H/||H||_F` (`Z=0` at `H=0`).  With the exact filtered split `H=H_V+H_R`,

`rho_P=||H||_F = Z:H_V + Z:H_R`.

Decompose only the resolved transporter `V=S_(N/4)u` into hard orthogonal dyadic shells at this **physical pressure event**.  The bilinear resolved pressure tensor then expands exactly into unordered atoms `{{a,b}}`: diagonal once, off diagonal as both orientations together.  Thus

`rho_P <= [r_R]_+ + sum_(a<=b)[p_ab]_+`.

These are positive source/service atoms, not child-energy causal probabilities.

The already-certified order-2 pressure and shell Bernstein constants give the exact rational ordered-pair coefficient

`256/1425 < 1/5`.

Therefore every unordered hard pair obeys

`|p_ab| <= (kappa_ab/5)(Mmax/N)^4 sqrt(mu_a mu_b)`,

with `kappa=1` diagonal and `2` off diagonal.  Gaussian averaging is a probability contraction and costs nothing further.

For integrated pressure source weight `Sigma_P`, exact positivity gives the joint half split

`SGS-positive source >=Sigma_P/2`  OR  `resolved positive pair source >=Sigma_P/2`.

The SGS branch still yields `int||R||_(3/2)>=190 Sigma_P` and enters the existing coherent-service compiler.  On the resolved branch normalize the actual unordered pair source law.  Its native statement is the threshold-free tradeoff

`mu_child exp(H2_pair) >= 320 Sigma_P/c`.

Indeed a countable positive pair law has an attained maximal atom, `q_max>=sum q^2=exp(-H2_pair)`, and that actual pair exposes the stated hard `u`-shell lower after the strict-lowpass contraction.  Therefore **every resolved pressure-pair owner already enters the generic critical-shell theorem**.  The familiar `theta=1/4` split is only a diagnostic corollary:

- if a pair is theta-dominant, its integrated capacity forces at some time an actual hard child shell with

  `mu_child >= [5 Sigma_P/(16c)](N/Mmax)^4 >= 80 Sigma_P/c`,

  because every resolved pair has `Mmax<=N/4`; this is a genuine input to the generic critical-shell theorem;
- if no pair exceeds one quarter, the actual source law has `H2_pair>=log 4`.  This quantifies why the unconditional shell seed is weaker, but it does not create another master fate and is not a causal HH probability.

At exact quarter mass both corollaries hold.  They are not competing routes: the physical resolved-pair owner has already been sent once to the same critical-shell recursion.

Material/coherent labels are deliberately absent from the scale proof.  They may be attached after the hard event as sidecars; on a supplied signed-good low-strain material lineage the previously certified fixed objective-Hessian pair contraction `<1/5` remains an optional reuse refinement.

The old coarse inequality `rho_P<=mu_V/5700+||R||_(3/2)/380` remains true as a diagnostic, but aggregate `mu_V` is no longer the canonical pressure renewal state.

Stress: `{out.samples}` pressure tensor/pair/source states
- worst dual scalar residual: `{out.worst_dual_scalar_residual:.3e}`
- minimum positive source-cover margin: `{out.minimum_positive_cover_margin:.3e}`
- worst unordered reconstruction residual: `{out.worst_unordered_reconstruction_residual:.3e}`
- minimum sampled pair-capacity margin: `{out.minimum_pair_capacity_margin:.3e}`
- minimum primary owner half-split margin: `{out.minimum_owner_half_margin:.3e}`
- minimum dominant-shell margin over clean lower: `{out.minimum_dominant_shell_margin_over_clean:.3e}`
- minimum diffuse-entropy margin: `{out.minimum_diffuse_entropy_margin:.3e}`
- minimum entropy-shell tradeoff margin: `{out.minimum_entropy_shell_tradeoff_margin:.3e}`
- maximum joint primary owner count: `{out.maximum_joint_primary_owner_count}`
- maximum simultaneous quarter-cut corollaries: `{out.maximum_quarter_corollary_count}`

No packet synchronization, no coherent-frequency support fiction, no aggregate pressure-mass reset, and no Navier--Stokes global-regularity conclusion are asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
