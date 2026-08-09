from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.common_slice_coefficient_registration import CONTINUING_COEFFICIENT_FRACTION
from src.bargmann_root_cell_registration import (
    optimal_bargmann_fraction,
    optimal_bargmann_radius,
    unit_grid_cells_intersecting_ball_upper,
)
from src.dual_gaussian_root_registration import dual_probe_critical_mass_lower
from src.smooth_symbol_freezing import sharp_young_constant_3d
from src.weighted_causal_reuse import entropy, layered_reuse_information, parent_pushforward, random_layered_maps
from src.renyi_causal_reuse import layered_collision_reuse
from src.physical_pair_weighted_productivity import physical_log_productivity_constant, variable_productivity_root_log_lower

SHELL_LOG_HALFWIDTH = 2.0 / 25.0
SHELL_LOWER_AXIS = 2.0 / 3.0
LOW_STRAIN_ACTION = 1.0 / 30.0
DEFAULT_LOG_COV_RADIUS = 0.4
DEFAULT_REGISTRATION_FRACTION = CONTINUING_COEFFICIENT_FRACTION
DEFAULT_DUHAMEL_GENERATED_FRACTION = 1.0 / 2.0
DEFAULT_SYMMETRIZED_SOURCE_FACTOR = 2.0
ROOT_SCALE_DILATION = 25.0 / 24.0


def l2_normalized_dual_l32_scaled_upper(
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
) -> float:
    """Upper for N^(-1/2)||phi||_(3/2) of an L2-normalized dual Gaussian.

    If h is the L3-normalized Gaussian dual and r is the physical geometric
    radius of its profile covariance, exact Gaussian integration gives

      ||h||_(3/2) / ||h||_2 = (4/3) pi^(1/4) r^(-1/2).

    A log-covariance representative within delta has
      r_rep >= r exp(-sqrt(3) delta/6),
    while the shell uncertainty gives N r >= lower_axis_constant.
    """
    d = float(log_cov_radius)
    a = float(lower_axis_constant)
    if d < 0 or a <= 0:
        raise ValueError("valid covariance radius and shell axis required")
    radius_factor = a * math.exp(-math.sqrt(3.0) * d / 6.0)
    return (4.0 / 3.0) * math.pi ** 0.25 / math.sqrt(radius_factor)


def adjoint_vector_growth_upper(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
    carrier_ratio_upper: float = math.exp(SHELL_LOG_HALFWIDTH),
) -> float:
    """Scale-independent backward adjoint growth on a child slab.

    For G=-S_perp-nu|k|^2 I, the backward dual solves
      psi_dot=S_perp psi+nu|k|^2 psi.
    If T=c N^-2 and |k|<=R_k N, then
      ||psi|| <= exp(K+nu c R_k^2)||psi_terminal||.
    """
    c = float(scaled_lifetime)
    nu = float(viscosity)
    K = float(strain_action)
    R = float(carrier_ratio_upper)
    if c <= 0 or nu < 0 or K < 0 or R <= 0:
        raise ValueError("valid adjoint growth data required")
    return math.exp(K + nu * c * R * R)


def parent_l32_productivity_lower(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    generated_fraction: float = DEFAULT_DUHAMEL_GENERATED_FRACTION,
    symmetrized_source_factor: float = DEFAULT_SYMMETRIZED_SOURCE_FACTOR,
    carrier_ratio_upper: float = math.exp(SHELL_LOG_HALFWIDTH),
) -> float:
    """Lambda_l32 in ||p1||_(3/2)||p2||_(3/2) >= Lambda_l32 alpha_c.

    Here alpha_c=sqrt(N)|c_terminal| is the scale-critical selected child
    coefficient.  Pairing the symmetrized high-high source against the moving
    L2-normalized Gaussian probe gives

      |F_HH,coeff| <= s R_k N A3 ||phi||_(3/2) a1 a2.

    Integrating over T=cN^-2 and including the coefficient-space adjoint growth
    removes all powers of N.  On the Duhamel-generated branch
      |I_HH| >= generated_fraction |c_terminal|.
    """
    c = float(scaled_lifetime)
    g = float(generated_fraction)
    s = float(symmetrized_source_factor)
    R = float(carrier_ratio_upper)
    if c <= 0 or not (0 < g <= 1) or s <= 0 or R <= 0:
        raise ValueError("valid productivity data required")
    A3 = sharp_young_constant_3d()
    Cphi = l2_normalized_dual_l32_scaled_upper(log_cov_radius)
    Cadj = adjoint_vector_growth_upper(c, viscosity, LOW_STRAIN_ACTION, R)
    return g / (s * R * A3 * Cphi * Cadj * c)


def registered_coefficient_productivity_lower(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    registration_fraction: float = DEFAULT_REGISTRATION_FRACTION,
) -> float:
    """Legacy pointwise Duhamel productivity candidate.

    This stronger pointwise statement remains a useful diagnostic special case,
    but it is no longer the preferred causal premise.  The master-facing theorem
    uses physical_pair_weighted_productivity, which proves the required
    transfer-weighted logarithmic product lower directly under dT.
    """
    q = float(registration_fraction)
    if not (0 < q <= 1):
        raise ValueError("registration fraction must lie in (0,1]")
    eta_dual = dual_probe_critical_mass_lower(0.01, log_cov_radius)
    return q * q * eta_dual * parent_l32_productivity_lower(
        scaled_lifetime=scaled_lifetime,
        viscosity=viscosity,
        log_cov_radius=log_cov_radius,
    )


def anchor_coefficient_energy_fraction() -> float:
    """beta in N E_anchor >= beta alpha^2 for alpha=sqrt(N)|coherent coeff|."""
    R = optimal_bargmann_radius()
    return optimal_bargmann_fraction() / unit_grid_cells_intersecting_ball_upper(R)


def expected_log_amplitude(weights: Sequence[float], amplitudes: Sequence[float]) -> float:
    w = np.asarray(weights, float)
    a = np.asarray(amplitudes, float)
    if w.shape != a.shape or np.any(w < 0) or w.sum() <= 0 or np.any(a <= 0):
        raise ValueError("positive matching law/amplitudes required")
    w = w / w.sum()
    return float(np.dot(w, np.log(a)))


def one_layer_log_product_lower(child_log_amplitude: float, productivity: float) -> float:
    """Two-slot baseline turns the physical-weighted log-product lower into one log recursion."""
    if productivity <= 0:
        raise ValueError("positive productivity required")
    return 0.5 * math.log(productivity) + 0.5 * float(child_log_amplitude)


def root_expected_log_lower(depth: int, productivity: float, terminal_amplitude: float) -> float:
    if depth < 0 or productivity <= 0 or terminal_amplitude <= 0:
        raise ValueError("valid depth/productivity/terminal amplitude required")
    if depth == 0:
        return math.log(terminal_amplitude)
    q = 2.0 ** (-depth)
    return (1.0 - q) * math.log(productivity) + q * math.log(terminal_amplitude)


def root_expected_log_lower_variable(productivities: Sequence[float], terminal_amplitude: float) -> float:
    """Physical-weighted log recursion for layer-dependent Lambda_j.

    If ell_j >= .5 log Lambda_j + .5 ell_(j+1), then
      ell_0 >= sum_j 2^(-j-1) log Lambda_j + 2^(-L) log alpha_L.
    """
    if terminal_amplitude <= 0:
        raise ValueError("positive terminal amplitude required")
    lam = np.asarray(productivities, float)
    if lam.ndim != 1 or np.any(lam <= 0):
        raise ValueError("positive one-dimensional productivities required")
    return variable_productivity_root_log_lower(lam, math.log(terminal_amplitude))


def entropy_energy_amplitude_upper(
    amplitudes: Sequence[float],
    weights: Sequence[float],
) -> float:
    """Exact log-sum inequality H(p)+2 E_p log alpha <= log sum alpha^2."""
    a = np.asarray(amplitudes, float)
    w = np.asarray(weights, float)
    if a.shape != w.shape or np.any(a <= 0) or np.any(w < 0) or w.sum() <= 0:
        raise ValueError("positive matching amplitudes/probabilities required")
    w = w / w.sum()
    return math.log(float(np.dot(a, a))) - 2.0 * expected_log_amplitude(w, a)


def root_entropy_upper_from_energy(
    depth: int,
    base_frequency: float,
    global_energy: float,
    frame_budget: float,
    productivity: float,
    terminal_amplitude: float,
    anchor_fraction: float | None = None,
) -> float:
    """Root Shannon entropy upper without any per-root mass floor.

    Distinct root anchor cells obey N_r E_r >= beta alpha_r^2.  With
    N_r<=N_base(25/24)^L and total anchor energy <=P E_global,
      sum alpha_r^2 <= P E_global N_max / beta.
    Combine with the log-sum inequality and the multiplicative amplitude
    recursion.
    """
    if depth < 0 or min(base_frequency, global_energy, frame_budget, productivity, terminal_amplitude) <= 0:
        raise ValueError("positive root-budget data required")
    beta = anchor_coefficient_energy_fraction() if anchor_fraction is None else float(anchor_fraction)
    if beta <= 0:
        raise ValueError("positive anchor fraction required")
    nmax = base_frequency * ROOT_SCALE_DILATION**depth
    ell = root_expected_log_lower(depth, productivity, terminal_amplitude)
    return math.log(frame_budget * global_energy * nmax / beta) - 2.0 * ell


def reuse_information_lower_without_mass_floor(
    depth: int,
    base_frequency: float,
    global_energy: float,
    frame_budget: float,
    productivity: float,
    terminal_amplitude: float,
    anchor_fraction: float | None = None,
) -> dict[str, float | str | bool]:
    """Shannon and Renyi reuse lower driven by amplitude productivity, not root mass floors."""
    Hup = root_entropy_upper_from_energy(
        depth, base_frequency, global_energy, frame_budget, productivity,
        terminal_amplitude, anchor_fraction,
    )
    if Hup < -2e-12:
        return {
            "branch": "infeasible_energy_amplitude_data",
            "root_entropy_upper": Hup,
            "feasible": False,
        }
    Hup = max(0.0, Hup)
    shannon = depth * math.log(2.0) - Hup
    # H_2(root)<=H_1(root), so L log2+log Q_root >= L log2-H_1(root).
    renyi = shannon
    rich_line = depth * math.log(4.0 / 3.0)
    return {
        "branch": "amplitude_entropy_reuse_forced" if shannon > rich_line else "depth_not_yet_forced",
        "root_entropy_upper": Hup,
        "shannon_reuse_lower": shannon,
        "renyi_action_lower": renyi,
        "rich_line": rich_line,
        "one_renyi_theta_gt_one_third": bool(renyi > rich_line),
        "feasible": True,
    }


def closed_form_reuse_lower(
    depth: int,
    base_frequency: float,
    global_energy: float,
    frame_budget: float,
    productivity: float,
    terminal_amplitude: float,
    anchor_fraction: float | None = None,
) -> float:
    """The explicit lower before truncating the entropy upper at zero."""
    beta = anchor_coefficient_energy_fraction() if anchor_fraction is None else float(anchor_fraction)
    q = 2.0 ** (-depth)
    return (
        depth * math.log(48.0 / 25.0)
        - math.log(frame_budget * global_energy * base_frequency / beta)
        + 2.0 * (1.0 - q) * math.log(productivity)
        + 2.0 * q * math.log(terminal_amplitude)
    )


def _build_amplitude_tree(
    maps: Sequence[np.ndarray], productivity: float, terminal_amplitude: float,
    rng: np.random.Generator,
) -> list[np.ndarray]:
    """Adversarial wildly-unbalanced amplitudes satisfying every product gate.

    If several causal slots merge to one parent label, assign that label the
    maximum amplitude demanded by any incident slot.  This can only improve all
    product inequalities and models a reused strong reservoir without inventing
    a mass floor.
    """
    amps: list[np.ndarray] = [None] * (len(maps) + 1)  # type: ignore[list-item]
    amps[-1] = np.array([float(terminal_amplitude)])
    for j in range(len(maps) - 1, -1, -1):
        pairs = np.asarray(maps[j], int)
        child = amps[j + 1]
        npar = int(pairs.max()) + 1
        parent = np.full(npar, 1e-300, float)
        for c, ac in enumerate(child):
            # Allow enormous imbalance; only the product is fixed by physics.
            skew = float(rng.uniform(-10.0, 10.0))
            base = math.sqrt(productivity * ac)
            req0 = base * math.exp(skew)
            req1 = base * math.exp(-skew)
            parent[pairs[c, 0]] = max(parent[pairs[c, 0]], req0)
            parent[pairs[c, 1]] = max(parent[pairs[c, 1]], req1)
        amps[j] = parent
    return amps


@dataclass(frozen=True)
class AmplitudeEntropyStress:
    samples: int
    worst_layer_log_margin: float
    minimum_logsum_margin: float
    minimum_shannon_lower_margin: float
    minimum_renyi_lower_margin: float
    maximum_amplitude_ratio: float
    branch_counts: dict[str, int]


def stress(samples: int = 20_000, seed: int = 20260809) -> AmplitudeEntropyStress:
    rng = np.random.default_rng(seed)
    wl = 0.0
    ml = ms = mr = float("inf")
    maxratio = 1.0
    branches: dict[str, int] = {}
    beta = anchor_coefficient_energy_fraction()
    for _ in range(samples):
        L = int(rng.integers(1, 13))
        maps = random_layered_maps(rng, L)
        lam = float(math.exp(rng.uniform(-6.0, 0.5)))
        terminal = float(math.exp(rng.uniform(-8.0, 4.0)))
        amps = _build_amplitude_tree(maps, lam, terminal, rng)

        # Build the exact transfer-weighted laws from terminal to roots.
        w = np.ones(1)
        laws: list[np.ndarray] = [None] * (L + 1)  # type: ignore[list-item]
        laws[L] = w.copy()
        for j in range(L - 1, -1, -1):
            pairs = maps[j]
            child_ell = expected_log_amplitude(w, amps[j + 1])
            slot_parent = []
            slot_weight = []
            for c, wc in enumerate(w / w.sum()):
                for s in range(2):
                    slot_parent.append(amps[j][pairs[c, s]])
                    slot_weight.append(0.5 * wc)
            slot_ell = expected_log_amplitude(slot_weight, slot_parent)
            bound = one_layer_log_product_lower(child_ell, lam)
            margin = slot_ell - bound
            wl = max(wl, max(0.0, -margin))
            if margin < -3e-12:
                raise AssertionError("amplitude log-product recursion failed")
            w = parent_pushforward(w, pairs, len(amps[j]))
            laws[j] = w.copy()

        root = amps[0]
        rootw = laws[0]
        maxratio = max(maxratio, float(root.max() / root.min()))
        lhs = entropy(rootw) + 2.0 * expected_log_amplitude(rootw, root)
        rhs = math.log(float(np.dot(root, root)))
        ml = min(ml, rhs - lhs)
        if lhs > rhs + 3e-12:
            raise AssertionError("energy-amplitude entropy log-sum failed")

        # Choose a physical energy budget just above the minimum anchor budget.
        N0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        Nmax = N0 * ROOT_SCALE_DILATION**L
        P = float(rng.uniform(1.0, 5.0))
        Emin = beta * float(np.dot(root, root)) / (P * Nmax)
        E = Emin * float(rng.uniform(1.01, 4.0))
        lower = reuse_information_lower_without_mass_floor(L, N0, E, P, lam, terminal, beta)
        if not bool(lower["feasible"]):
            raise AssertionError("constructed feasible tree classified infeasible")

        sh = layered_reuse_information(maps)
        re = layered_collision_reuse(maps)
        actual_s = float(sh["total_reuse_information"])
        actual_r = float(re["total_action"])
        s_lo = float(lower["shannon_reuse_lower"])
        r_lo = float(lower["renyi_action_lower"])
        ms = min(ms, actual_s - s_lo)
        mr = min(mr, actual_r - r_lo)
        if s_lo > actual_s + 5e-11:
            raise AssertionError("amplitude-entropy Shannon lower exceeded exact reuse")
        if r_lo > actual_r + 5e-11:
            raise AssertionError("amplitude-entropy Renyi lower exceeded exact action")
        b = str(lower["branch"])
        branches[b] = branches.get(b, 0) + 1

    return AmplitudeEntropyStress(samples, wl, ml, ms, mr, maxratio, branches)


def theorem_certificate() -> dict[str, object]:
    beta = anchor_coefficient_energy_fraction()
    lam = physical_log_productivity_constant(1)
    return {
        "status": "EXACT_AMPLITUDE_ENTROPY_REUSE_TELESCOPE__PHYSICAL_TRANSFER_WEIGHTED_PRODUCTIVITY_SUPPLIED",
        "productivity": "preferred premise is E_dT log(alpha_p1 alpha_p2)>=E_dT log(alpha_child)+log Lambda_j under actual positive child-energy work; no pointwise pair law is required",
        "default_productivity": f"physical pair-work/KL theorem gives Lambda_1={lam:.12g} at one retained pair cell; Lambda_j may vary with finite cell refinement",
        "log_recursion": "ell_parent >= (1/2)log Lambda_j+(1/2)ell_child under the same physical transfer law",
        "depth_solution": "variable Lambda_j enter with geometric weights 2^(-j-1); polynomial pair-cell refinement therefore changes only a finite offset",
        "anchor": f"Bargmann/Moyal gives N E_anchor >= beta alpha^2 with beta={beta:.12g}",
        "entropy": "H(root)+2 E_root log alpha <= log(sum alpha^2) <= log(P E_global N_max/beta)",
        "shannon": "sum R_j=L log2-H(root) retains the linear log(48/25) slope without a per-root mass floor",
        "renyi": "H2(root)<=H1(root), so the Renyi action has the same lower bound and routes a rich layer to existing pair/entropy/cycle currencies",
        "amplitude_imbalance": "arbitrarily unbalanced parent amplitudes are absorbed by transfer-weighted logarithmic productivity; no small/large-parent currency is introduced",
        "duhamel_role": "Duhamel remains an exact support/adjoint identity only; parent productivity is derived directly from physical work by sharp Young plus KL positivity",
        "physical_weights": "all layer expectations are under actual positive child-energy work dT",
        "continuum_status": "outer moving roles and hard-event/smooth-envelope registration are supplied by companion theorems; the remaining continuum task is the recursive first-stop assembly that retains physical Young-good cells and routes every failed phase/registration/source event once to an existing cause",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=20_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-amplitude-entropy-causal-reuse"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "amplitude_entropy_causal_reuse.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    lam = physical_log_productivity_constant(1)
    beta = anchor_coefficient_energy_fraction()
    md = f"""# Amplitude-entropy causal reuse: remove the false uniform root-mass hypothesis

Status: **{cert['status']}**.

The critical amplitude of a selected coherent coefficient is

`alpha=sqrt(N)|<u,phi>|`.

The preferred productivity premise is now derived directly under the actual positive child-energy work.  Sharp Young bounds each physical pair-work density, and KL positivity against normalized time x hard-pair-cell reference gives

`E_dT log(alpha_p1 alpha_p2) >= E_dT log(alpha_child) + log Lambda_j`.

For one retained pair cell at the displayed default `c=1`, including dual-Gaussian marking and the conservative `1/4` common-slice factor on each parent, the physical-work theorem gives `Lambda_1={lam:.12g}`.  For `M_j` hard pair cells, `Lambda_j=Lambda_1/M_j`.  Polynomial refinement has a finite geometrically discounted `sum 2^(-j) log M_j`, so it changes only the offset.

Give each physical child event its two structural parent slots with the existing free `1/2` baseline.  The transfer-weighted logarithmic inequality gives

`E log alpha_parent >= (1/2) log Lambda_j + (1/2) E log alpha_child`,

hence after depth `L`

`ell_root >= sum_(j=0)^(L-1) 2^(-j-1) log Lambda_j + 2^-L log alpha_terminal`.

When all `Lambda_j=Lambda`, this reduces to the older constant-productivity formula `(1-2^-L)log Lambda + 2^-L log alpha_terminal`.

This is the key point: **arbitrary amplitude imbalance is harmless at the multiplicative level**.  One parent may be exponentially smaller than the other; their log average still sees the exact product.

Bargmann submean plus the canonical Moyal partition gives, for a root anchor cell,

`N_r E_r >= beta alpha_r^2`, `beta={beta:.12g}`.

If `w_r` is the physical transfer-weighted root law, the log-sum inequality gives exactly

`H(w_root)+2 E_w log alpha_r <= log sum_r alpha_r^2`.

Using the common-slice energy budget and `N_root<=N_base(25/24)^L`,

`H(w_root) <= log[P E_global N_base (25/24)^L / beta] - 2 ell_root`.

Therefore the exact Shannon telescope

`sum R_j=L log2-H(w_root)`

obeys

`sum R_j >= L log(48/25) - log(P E_global N_base/beta) + sum_(j=0)^(L-1) 2^(-j) log Lambda_j + 2^(1-L)log alpha_terminal`.

The physical pair theorem gives `Lambda_j=Lambda_1/M_j`; polynomial `M_j` makes the weighted logarithmic sum a finite offset.  Hence the coefficient of `L` is still `log(48/25)`.  No assumption `N E_root>=eta` for every root is used.

Moreover `H2(root)<=H1(root)`, so the exact Renyi action `sum log(1+theta_j)=L log2+log Q_root` has the **same lower bound**.  Once that lower exceeds `L log(4/3)`, some layer has `theta_j>1/3` and the existing parent-slot pair / component-entropy / same-ancestry-cycle routing applies unchanged.

Stress: `{out.samples}` random layered causal DAGs with parent amplitude ratios as large as `{out.maximum_amplitude_ratio:.3e}`
- worst one-layer log-product violation: `{out.worst_layer_log_margin:.3e}`
- minimum log-sum entropy margin: `{out.minimum_logsum_margin:.3e}`
- minimum exact Shannon-minus-lower margin: `{out.minimum_shannon_lower_margin:.3e}`
- minimum exact Renyi-minus-lower margin: `{out.minimum_renyi_lower_margin:.3e}`
- branches: `{out.branch_counts}`

This theorem changes the frontier.  The old causal root-count argument needed an absolute critical mass per root; that is incompatible with Young homogeneity.  The correct observable is the **physical-transfer-weighted logarithm of the multiplicative critical coefficient**.  The physical pair-work/KL theorem supplies it directly under `dT`, common-slice first stopping registers continuing parents, and the outer-role/event-registration theorems now supply the exact PDE carrier.  What remains is the continuum first-stop constructor that applies these alternatives measurably on every recursive physical-transfer block and routes every failed event once.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
