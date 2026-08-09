from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.dual_gaussian_root_registration import normalized_dual_probe_critical_mass_lower
from src.smooth_symbol_freezing import sharp_young_constant_3d

LOW_STRAIN_ACTION = 1.0 / 30.0
INHERITED_ENERGY_FRACTION = 1.0 / 5.0
RESIDUAL_WORK_FRACTION = 1.0 / 5.0
GENERATED_WORK_FRACTION = 8.0 / 15.0
COMMON_SLICE_REGISTRATION = 1.0 / 4.0
DEFAULT_SHELL_RADIUS = math.exp(2.0 / 25.0)
DEFAULT_SCALED_LIFETIME = 1.0


def shell_l32_energy_constant(shell_radius: float = DEFAULT_SHELL_RADIUS) -> float:
    """C_Omega in ||fhat||_(3/2) <= C_Omega sqrt(N)||fhat||_2.

    It is enough that the role support lie in the ball |xi|<=R N.  Holder on
    that finite Fourier volume gives |B_R N|^(1/6).
    """
    R = float(shell_radius)
    if R <= 0 or not math.isfinite(R):
        raise ValueError("positive finite shell radius required")
    return ((4.0 * math.pi / 3.0) * R**3) ** (1.0 / 6.0)


def generated_energy_sup_factor(
    strain_action: float = LOW_STRAIN_ACTION,
    inherited_fraction: float = INHERITED_ENERGY_FRACTION,
    residual_fraction: float = RESIDUAL_WORK_FRACTION,
    generated_fraction: float = GENERATED_WORK_FRACTION,
) -> float:
    """C_E with sup_t E(t)<=C_E W_HH on the generated branch.

    Energy Gronwall gives E(t)<=exp(2K)(E0+W_HH+W_R).  If
      E0<=a E1, W_R<=b E1, W_HH>=g E1,
    then E0+W_R <= (a+b)W_HH/g.
    """
    K = float(strain_action)
    a = float(inherited_fraction)
    b = float(residual_fraction)
    g = float(generated_fraction)
    if K < 0 or min(a, b) < 0 or g <= 0:
        raise ValueError("valid generated-branch constants required")
    return math.exp(2.0 * K) * (1.0 + (a + b) / g)


def physical_work_capacity_constant(child_wave_ratio: float = DEFAULT_SHELL_RADIUS) -> float:
    """C_Y in r_pair <= C_Y N a_c a_1 a_2.

    Symmetrizing the two parent orders costs 2 in the source and physical energy
    work contributes another factor 2.  Sharp Young supplies A_3.
    """
    R = float(child_wave_ratio)
    if R <= 0:
        raise ValueError("positive child wave ratio required")
    return 4.0 * R * sharp_young_constant_3d()


def physical_log_productivity_constant(
    pair_cells: int,
    scaled_lifetime: float = DEFAULT_SCALED_LIFETIME,
    shell_radius: float = DEFAULT_SHELL_RADIUS,
    child_wave_ratio: float = DEFAULT_SHELL_RADIUS,
) -> float:
    """Lambda_M for transfer-weighted log parent-productivity before log averaging.

    On the generated branch, physical positive pair work dT has total W and
    density r_e(t) satisfying
      r_e <= C_Y N a_c a_1 a_2.
    Relative entropy against uniform normalized time x pair-cell reference gives
      E_T log(a_1 a_2) >= log alpha_c + log Lambda_role,
    with no Duhamel weight.  Dual-Gaussian event registration and the conservative
    1/4 common-slice factor on each parent multiply Lambda_role by eta_dual/16.
    """
    M = int(pair_cells)
    c = float(scaled_lifetime)
    if M <= 0 or c <= 0:
        raise ValueError("positive pair-cell count and scaled lifetime required")
    C_y = physical_work_capacity_constant(child_wave_ratio)
    C_o = shell_l32_energy_constant(shell_radius)
    C_e = generated_energy_sup_factor()
    eta = normalized_dual_probe_critical_mass_lower()
    role = math.sqrt(GENERATED_WORK_FRACTION) / (c * M * C_y * C_o * math.sqrt(C_e))
    return (COMMON_SLICE_REGISTRATION**2) * eta * role


def transfer_weighted_log_product_lower(
    *,
    total_positive_work: float,
    final_energy: float,
    scale: float,
    duration: float,
    pair_cells: int,
    child_amplitudes: np.ndarray,
    parent_products: np.ndarray,
    physical_work_masses: np.ndarray,
    capacity_constant: float,
    shell_energy_constant: float,
    energy_sup_factor: float,
) -> dict[str, float]:
    """Finite quadrature form of the KL proof.

    The arrays represent bins of normalized time x pair-cell reference.  Each bin
    has equal reference mass.  `physical_work_masses` are actual positive work
    masses, sum to W.  The implied rate uses bin width T/#timebins and obeys the
    pointwise physical Young capacity inequality by assumption.
    """
    W = float(total_positive_work)
    E1 = float(final_energy)
    N = float(scale)
    T = float(duration)
    M = int(pair_cells)
    ac = np.asarray(child_amplitudes, float)
    ap = np.asarray(parent_products, float)
    wm = np.asarray(physical_work_masses, float)
    if W <= 0 or E1 <= 0 or N <= 0 or T <= 0 or M <= 0:
        raise ValueError("positive physical data required")
    if ac.shape != ap.shape or ac.shape != wm.shape or ac.ndim != 2 or ac.shape[1] != M:
        raise ValueError("arrays must have shape (time_bins,pair_cells)")
    if np.any(ac <= 0) or np.any(ap <= 0) or np.any(wm < 0):
        raise ValueError("positive amplitudes/products and nonnegative work required")
    if abs(float(wm.sum()) - W) > 2e-11 * max(1.0, W):
        raise ValueError("physical work masses must sum to W")
    nt = ac.shape[0]
    dt = T / nt
    rates = wm / dt
    cap = capacity_constant * N * ac * ap
    if np.any(rates > cap + 3e-11 * np.maximum(1.0, cap)):
        raise ValueError("physical pair work violates the sharp-Young capacity premise")
    if np.max(ac) > shell_energy_constant * math.sqrt(N * energy_sup_factor * W) * (1.0 + 3e-12):
        raise ValueError("child shell/energy upper premise violated")
    if W < GENERATED_WORK_FRACTION * E1 - 3e-12 * max(1.0, E1):
        raise ValueError("generated-work lower premise violated")

    prob = wm / W
    positive = prob > 0
    lhs = float(np.sum(prob[positive] * np.log(ap[positive])))
    # Uniform pair x normalized-time reference has mass 1/(nt M).  Its KL
    # density is p=nt*M*prob for each quadrature bin.
    kl = float(np.sum(prob[positive] * np.log((nt * M) * prob[positive])))
    if kl < -3e-12:
        raise AssertionError("KL positivity failed")
    alpha_upper = math.sqrt(N * E1)
    lower = (
        math.log(alpha_upper)
        + math.log(math.sqrt(GENERATED_WORK_FRACTION) / (
            (T * N * N) * M * capacity_constant * shell_energy_constant * math.sqrt(energy_sup_factor)
        ))
    )
    # T*N^2 is the scaled lifetime c.  alpha_actual may be below alpha_upper,
    # so using alpha_upper gives the strongest universal child-coefficient target.
    margin = lhs - lower
    return {
        "weighted_log_parent_product": lhs,
        "lower_bound_against_max_child_coefficient": lower,
        "margin": margin,
        "kl_to_uniform_time_pair_reference": kl,
    }


def tempered_pair_cell_penalty_upper(pair_prefactor: float, cell_power: float) -> float:
    """Bound sum 2^(-j-1) log M_j for M_j<=M0 (j+3)^p.

    Since j+3<=3*2^j, the weighted sum is at most log M0+p log 6.
    """
    M0 = float(pair_prefactor)
    p = float(cell_power)
    if M0 < 1 or p < 0:
        raise ValueError("M0>=1 and p>=0 required")
    return math.log(M0) + p * math.log(6.0)


def variable_productivity_root_log_lower(
    productivities: np.ndarray,
    terminal_log_amplitude: float,
) -> float:
    """Solve ell_j >= .5 log Lambda_j + .5 ell_(j+1)."""
    lam = np.asarray(productivities, float)
    if lam.ndim != 1 or np.any(lam <= 0):
        raise ValueError("positive one-dimensional productivities required")
    ell = float(terminal_log_amplitude)
    for x in lam[::-1]:
        ell = 0.5 * math.log(float(x)) + 0.5 * ell
    return ell


def theorem_certificate() -> dict[str, object]:
    lam1 = physical_log_productivity_constant(1)
    C_e = generated_energy_sup_factor()
    return {
        "status": "EXACT_PHYSICAL_TRANSFER_WEIGHTED_LOG_PRODUCTIVITY__NO_DUHAMEL_PAIR_WEIGHT_IDENTIFICATION",
        "work_density": "for each hard parent-pair cell, actual positive child-energy work density satisfies r_e<=C_Y N a_c a_1 a_2 by the same sharp Young trilinear bound",
        "kl_step": "normalize dT/W and compare with uniform normalized physical time x finite pair-cell reference; KL>=0 converts the physical work density directly into E_dT log(a_1 a_2)",
        "energy_sup": f"on K<=1/30, E0,W_R<=E1/5 and W_HH>=8E1/15 imply sup E <= C_E W_HH with C_E={C_e:.12g}<2",
        "shell": "finite selected Fourier shell gives a_c(t)<=C_Omega sqrt(N E(t)) by Holder; no adjoint Gaussian persistence is used",
        "registered_law": "after dual-Gaussian event marking and 1/4 common-slice registration per parent, E_dT log(alpha_p1 alpha_p2)>=log(alpha_c)+log Lambda_M",
        "default_one_pair_constant": f"for c=1 and default shell constants, Lambda_(M=1)={lam1:.12g}; for M cells Lambda_M=Lambda_1/M",
        "refinement": "if M_j<=M0(j+3)^p, the binary recursion pays at most log M0+p log 6 in its entire geometrically discounted cell-count offset",
        "causal_weights": "all expectations are under actual positive physical child-energy work dT; Duhamel is not used to select or reweight parent pairs",
        "continuum_status": "the remaining PDE assembly is to verify the sharp-Young work-density premise on the retained signed-good hard event cells and apply the existing first-stop rules to any event that fails complex Young/common-slice registration; no dGamma-to-dT pair bridge remains",
    }


@dataclass(frozen=True)
class ProductivityStress:
    samples: int
    minimum_kl_margin: float
    minimum_log_product_margin: float
    minimum_variable_recursion_margin: float
    maximum_tempered_penalty_ratio: float


def stress(samples: int = 50_000, seed: int = 20260809) -> ProductivityStress:
    rng = np.random.default_rng(seed)
    mkl = mlp = mrec = float("inf")
    mpen = 0.0
    C_y = physical_work_capacity_constant()
    C_o = shell_l32_energy_constant()
    C_e = generated_energy_sup_factor()
    for _ in range(samples):
        M = int(rng.integers(1, 20))
        nt = int(rng.integers(2, 30))
        N = float(math.exp(rng.uniform(0.0, 8.0)))
        c = float(rng.uniform(0.2, 2.0))
        T = c / (N*N)
        E1 = float(math.exp(rng.uniform(-5.0, 5.0)))
        W = float(rng.uniform(GENERATED_WORK_FRACTION, 3.0)) * E1
        acmax = C_o * math.sqrt(N * C_e * W)
        ac = np.exp(rng.uniform(math.log(acmax)-4.0, math.log(acmax), size=(nt,M)))
        # Choose arbitrary parent products, then physical rates below capacity.
        ap = np.exp(rng.uniform(-3.0, 4.0, size=(nt,M))) * math.sqrt(N*W)
        capacity = C_y * N * ac * ap
        raw = rng.random((nt,M)) * capacity
        # Rescale rates so total work is exactly W while keeping below capacity.
        current = float(raw.sum() * (T/nt))
        if current < W:
            # Increase parent capacity uniformly enough to make W admissible.
            factor = W / max(current, 1e-300) * 1.2
            ap *= factor
            capacity *= factor
            raw *= factor
            current *= factor
        rate_scale = W / current
        rates = raw * rate_scale
        wm = rates * (T/nt)
        out = transfer_weighted_log_product_lower(
            total_positive_work=W, final_energy=E1, scale=N, duration=T,
            pair_cells=M, child_amplitudes=ac, parent_products=ap,
            physical_work_masses=wm, capacity_constant=C_y,
            shell_energy_constant=C_o, energy_sup_factor=C_e,
        )
        mkl = min(mkl, float(out["kl_to_uniform_time_pair_reference"]))
        mlp = min(mlp, float(out["margin"]))
        if float(out["margin"]) < -2e-10:
            raise AssertionError("physical transfer-weighted log-product lower failed")

        L = int(rng.integers(1, 30))
        M0 = float(rng.uniform(1.0, 50.0))
        power = float(rng.uniform(0.0, 16.0))
        Ms = M0 * (np.arange(L)+3.0)**power
        lam0 = physical_log_productivity_constant(1, c)
        lams = lam0 / Ms
        ellT = float(rng.uniform(-10.0, 5.0))
        exact = variable_productivity_root_log_lower(lams, ellT)
        weights = 2.0 ** (-(np.arange(L)+1.0))
        closed = float(np.sum(weights*np.log(lams)) + (2.0**(-L))*ellT)
        mrec = min(mrec, exact-closed)
        if abs(exact-closed) > 2e-12 * max(1.0, abs(exact), abs(closed)):
            raise AssertionError("variable-productivity binary recursion failed")
        actual_pen = float(np.sum(weights*np.log(Ms)))
        upper_pen = tempered_pair_cell_penalty_upper(M0,power)
        mpen = max(mpen, actual_pen/max(upper_pen,1e-30))
        if actual_pen > upper_pen + 2e-12*max(1.0,upper_pen):
            raise AssertionError("tempered cell-count penalty bound failed")
    return ProductivityStress(samples,mkl,mlp,mrec,mpen)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--outdir",type=Path,default=Path("results-physical-pair-weighted-productivity"))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); out=stress(args.samples)
    (args.outdir/"physical_pair_weighted_productivity.json").write_text(json.dumps({"certificate":cert,"stress":asdict(out)},indent=2),encoding="utf-8")
    md=f"""# Physical pair-weighted amplitude productivity\n\nStatus: **{cert['status']}**.\n\nThe parent-product law should be averaged under the same positive child-energy work that drives causal Shannon/Renyi.  No Duhamel pair reweighting is necessary.\n\nLet `r_e(t) dt` be actual positive work in one hard parent-pair cell.  Sharp Young gives\n\n`r_e(t) <= C_Y N a_c(t) a_1,e(t) a_2,e(t)`.\n\nNormalize the physical measure by its total `W` and compare it with uniform normalized physical time times the `M` hard pair cells.  Nonnegativity of relative entropy gives\n\n`E_(dT/W) log(a_1 a_2) >= log[W/(T M C_Y N)] - E_(dT/W) log a_c`.\n\nOn the generated low-strain branch, `E0,W_R<=E1/5`, `W>=8E1/15`, and energy Gronwall imply `sup E <= C_E W` with `C_E={generated_energy_sup_factor():.12g}<2`.  Finite shell volume gives `a_c<=C_Omega sqrt(N E)`.  Since `T=cN^-2` and every L2-normalized terminal coefficient satisfies `alpha_c<=sqrt(N E1)`, all scale powers cancel:\n\n`E_(dT/W) log(a_1 a_2) >= log alpha_c + log Lambda_role,M`.\n\nDual-Gaussian event marking and the conservative common-slice `1/4` factor on each parent give\n\n`E_(dT/W) log(alpha_p1 alpha_p2) >= log alpha_c + log Lambda_M`,\n\nwith default `Lambda_1={physical_log_productivity_constant(1):.12g}` and `Lambda_M=Lambda_1/M`.\n\nThe `1/M` factor is harmless under the actual symbol-freezing refinement.  If `M_j<=M0(j+3)^p`, binary log recursion weights its depth-j constant by `2^(-j-1)`, and\n\n`sum_j 2^(-j-1) log M_j <= log M0 + p log 6 < infinity`.\n\nThus polynomially refining physical pair cells changes only a finite amplitude-entropy offset; it cannot change the linear reuse slope.  This is precisely the homogeneity expected from a quadratic cascade.\n\nStress: `{out.samples}` random physical-work/time/pair laws and variable-depth recursions\n- minimum KL margin: `{out.minimum_kl_margin:.3e}`\n- minimum physical log-product margin: `{out.minimum_log_product_margin:.3e}`\n- minimum variable-recursion identity margin: `{out.minimum_variable_recursion_margin:.3e}`\n- maximum sampled tempered-penalty/upper ratio: `{out.maximum_tempered_penalty_ratio:.6f}`\n\nThe remaining continuum assembly is now local: on the retained signed-good hard event cells, verify the sharp-Young work-density bound with the same physical normalization already used by the SGS transfer theorem, and stop at the existing transfer/phase/relink/source branches whenever complex Young or common-slice registration fails.  Duhamel remains a support/adjoint identity but no Duhamel-to-physical **pair-weight** theorem is required.  No global-regularity claim is made.\n"""
    (args.outdir/"summary.md").write_text(md,encoding="utf-8"); print(md)


if __name__=="__main__":
    main()
