from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from src.affine_sgs_boundary_ledger import sgs_increment_cubic_upper
from src.atomic_component_entropy import collision_chain
from src.coherent_averaged_strain_source import (
    EXTENDED_ASPECT,
    SHELL_LOWER_AXIS,
    coherent_local_source_weight_upper,
)
from src.coherent_increment_service import cubic_to_square_threshold
from src.fresh_service_scale_reentry import (
    STATUS as FRESH_SCALE_REENTRY_STATUS,
    fresh_service_scale_route,
)
from src.high_frequency_dissipation_reentry import (
    STATUS as HIGH_FREQUENCY_REENTRY_STATUS,
    canonical_square_lp_tail_comparison_constant,
    classify_high_tail_energy_owners,
    lp_high_clean_reentry_thresholds,
    physical_tail_dissipation_lower_from_lp,
)

from src.critical_shell_service_reentry import (
    critical_shell_bounded_service_lower,
    critical_shell_integrated_service_lower,
    dissipation_supplier_shell_mass_threshold,
)
from src.objective_pressure_pair_atomization import (
    STATUS as PRESSURE_PAIR_STATUS,
    clean_entropy_shell_tradeoff_lower,
    objective_pressure_pair_route,
)
from src.pressure_reservoir_sync import pressure_hessian_pair_energy_service_ratio_upper
from src.resolved_objective_strain_collision import sgs_gradient_stress_lower

STATUS = (
    "EXACT_COHERENT_OBJECTIVE_SOURCE_OWNER_COMPILER__LOCAL_DV_AND_VISCOSITY_TO_CRITICAL_SHELL__"
    "SGS_FRESH_SERVICE_TO_REFINEMENT_INVARIANT_SCALE_SHELL__"
    "PRESSURE_TO_SGS_OR_ENTROPY_WEIGHTED_CRITICAL_SHELL__"
    "CELL_DOMINANCE_AND_AGGREGATE_MUV_DIAGNOSTIC_ONLY__NO_PACKET_SYNCHRONIZATION"
)

OWNER_NAMES = ("local_dv", "pressure", "sgs", "viscous")
DEFAULT_DOMINANT_FRACTION = 0.25
DEFAULT_ANCESTRY_ALPHA = 0.5


def objective_owner_weight_threshold(objective_variation_action: float, scaled_lifetime: float) -> float:
    """One of four physical owner classes carries at least A_obj/(4c).

    For the coherent averaged transporter, group the source terms by physical owner:

      local_dv = bar-A quadratic + coherent Reynolds/transport corrections,
      pressure = averaged filtered pressure Hessian,
      sgs = averaged differentiated filtered SGS stress,
      viscous = averaged viscosity.

    In scaled time the positive owner weights obey

        Sigma_local + Sigma_P + Sigma_R + Sigma_nu >= A_obj/c.

    This is a norm-triangle partition only; no artificial scalar clock is introduced.
    """
    A = float(objective_variation_action)
    c = float(scaled_lifetime)
    if A <= 0 or c <= 0 or not math.isfinite(A + c):
        raise ValueError("positive finite objective action and lifetime required")
    return A / (4.0 * c)


def objective_sgs_square_service_per_source(
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> float:
    """Exact linear coefficient Y_R >= C_Y rho_R for resolved objective SGS source.

    The clean resolved source collision gives ||R||_(3/2)>=380 rho_R. Germano gives

      Q >= ||R||_(3/2)^(3/2) / [(1+g1)^(3/2) g1^(1/2)],

    while coherent square service is

      Y=(Q/g1)^(2/3)/(C_LP C_B)^2.

    The 3/2 and 2/3 powers cancel exactly, leaving

      C_Y = 380 / [g1(1+g1)(C_LP C_B)^2].

    No packet radius or temporal persistence enters this coefficient.
    """
    g1 = float(filter_l1)
    clp = float(lp_constant)
    cb = float(bernstein_constant)
    if g1 < 1.0 or clp <= 0 or cb <= 0 or not all(math.isfinite(x) for x in (g1, clp, cb)):
        raise ValueError("require g1>=1 and positive finite LP/Bernstein constants")
    return 380.0 / (g1 * (1.0 + g1) * (clp * cb) ** 2)


def objective_sgs_integrated_square_service_lower(
    source_weight: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> float:
    sigma = float(source_weight)
    if sigma <= 0 or not math.isfinite(sigma):
        raise ValueError("positive finite SGS source weight required")
    return objective_sgs_square_service_per_source(filter_l1, lp_constant, bernstein_constant) * sigma


def local_dv_reentry(
    source_weight: float,
    scaled_lifetime: float,
    viscosity: float,
    *,
    aspect: float = EXTENDED_ASPECT,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
) -> dict[str, float | str]:
    """Coherent local quadratic/Reynolds owner -> resolved D_V -> critical shell.

    Existing Gaussian source calculus proves Sigma_local <= C_local D_V. Therefore
    a positive local owner weight supplies D_V>=Sigma_local/C_local and can enter the
    generic critical-shell theorem. It remains recursive scale-critical currency,
    never an additive finite reset.
    """
    sigma = float(source_weight)
    c = float(scaled_lifetime)
    nu = float(viscosity)
    if sigma <= 0 or c <= 0 or nu < 0 or not all(math.isfinite(x) for x in (sigma, c, nu)):
        raise ValueError("valid local source/lifetime/viscosity required")
    C_local = coherent_local_source_weight_upper(1.0, aspect, lower_axis_constant)
    D0 = sigma / C_local
    mu0 = dissipation_supplier_shell_mass_threshold(D0, c)
    Y_shell = critical_shell_bounded_service_lower(mu0, c, nu)
    return {
        "owner": "local_dv",
        "local_source_per_DV": C_local,
        "resolved_DV_lower": D0,
        "critical_shell_mass_lower": mu0,
        "own_scale_service_lower_on_full_survivor": Y_shell,
        "master_semantics": "RECURSE_CRITICAL",
    }


def viscous_dv_reentry(
    source_weight: float,
    scaled_lifetime: float,
    viscosity: float,
) -> dict[str, float | str]:
    """Integrated viscous owner -> resolved D_V by Cauchy -> critical shell.

    Pointwise rho_nu <= nu sqrt(d_V)/1500. If int rho_nu >= Sigma on a scaled
    interval of length c, then

      D_V >= (1500 Sigma/nu)^2 / c.
    """
    sigma = float(source_weight)
    c = float(scaled_lifetime)
    nu = float(viscosity)
    if sigma <= 0 or c <= 0 or nu <= 0 or not all(math.isfinite(x) for x in (sigma, c, nu)):
        raise ValueError("positive finite viscous source/lifetime/viscosity required")
    D0 = (1500.0 * sigma / nu) ** 2 / c
    mu0 = dissipation_supplier_shell_mass_threshold(D0, c)
    Y_shell = critical_shell_bounded_service_lower(mu0, c, nu)
    return {
        "owner": "viscous",
        "resolved_DV_lower": D0,
        "critical_shell_mass_lower": mu0,
        "own_scale_service_lower_on_full_survivor": Y_shell,
        "master_semantics": "RECURSE_CRITICAL",
    }


def objective_sgs_episode_thresholds(
    source_weight: float,
    scaled_lifetime: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    *,
    dominant_fraction: float = DEFAULT_DOMINANT_FRACTION,
    ancestry_alpha: float = DEFAULT_ANCESTRY_ALPHA,
) -> dict[str, float | str]:
    """Aggregate positive-measure thresholds for the direct objective-SGS route.

    Let Y_tot=C_Y Sigma_R. The low/high estimate is linear after integration:

      S_low >= Y_tot - 2 D_high.

    Hence either D_high>=Y_tot/4, or S_low>=Y_tot/2. On the latter, once old
    integrated capacity is <=Y_tot/8, either Xi>=Y_tot/8 or fresh service>=Y_tot/4.
    Canonically, that fresh positive measure is pushed to the fixed LP band index
    and enters the refinement-invariant fresh-scale theorem.  The historical
    theta-dominant coherent-cell / entropy / cycle thresholds below are retained
    only for optional fine ancestry bookkeeping and backward-compatible diagnostics.
    """
    c = float(scaled_lifetime)
    theta = float(dominant_fraction)
    alpha = float(ancestry_alpha)
    if c <= 0 or not math.isfinite(c) or not (0 < theta < 1) or not (0 < alpha < 1):
        raise ValueError("valid lifetime/dominance/ancestry thresholds required")
    Y = objective_sgs_integrated_square_service_lower(source_weight, filter_l1, lp_constant, bernstein_constant)
    h_atomic = -math.log(theta)
    h_anc = alpha * h_atomic
    pair = theta**alpha - theta
    return {
        "owner": "sgs",
        "integrated_forced_square_service": Y,
        "high_frequency_dissipation_threshold": Y / 4.0,
        "integrated_low_service_lower_if_no_high": Y / 2.0,
        "old_pool_integrated_capacity_threshold": Y / 8.0,
        "selected_interface_Xi_threshold": Y / 8.0,
        "fresh_service_lower": Y / 4.0,
        "canonical_fresh_route": FRESH_SCALE_REENTRY_STATUS,
        "fresh_cell_dominance_is_canonical_renewal": "NO",
        "dominant_pair_mass_occupation_lower": theta * Y / 8.0,
        "dominant_whole_shell_mass_occupation_lower": theta * Y / 16.0,
        "dominant_peak_whole_shell_mass_lower": theta * Y / (16.0 * c),
        "atomic_entropy_lower": h_atomic,
        "ancestry_entropy_lower": h_anc,
        "same_ancestry_pair_mass_lower": pair,
        "high_frequency_dissipation_is_resolved_DV": "NO",
        "master_semantics": "high-tail / reservoir-capacity / Xi keep existing owners; fresh positive service canonically RECURSE_CRITICAL via scale pushforward; cell entropy/cycle is optional sideledger",
    }


def objective_sgs_aggregate_route(
    source_weight: float,
    scaled_lifetime: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    *,
    high_frequency_dissipation: float,
    old_pool_integrated_capacity: float,
    old_old_integrated_service: float,
    selected_interface_integrated_service: float,
    new_edge_integrated_services: Sequence[float],
    ancestry_labels: Sequence[object] | None = None,
    dominant_fraction: float = DEFAULT_DOMINANT_FRACTION,
    ancestry_alpha: float = DEFAULT_ANCESTRY_ALPHA,
) -> dict[str, object]:
    """Legacy coherent-cell realization retained for fine ancestry bookkeeping.

    High-frequency/old-pool/interface exits remain physically valid.  Once the
    fresh branch is reached, however, the canonical renewal theorem is now
    `objective_sgs_aggregate_scale_route`; the cell argmax/entropy/cycle below is
    representation-dependent and must not be used as the canonical renewal fate.
    This function remains for backward-compatible ancestry diagnostics only.
    """
    th = objective_sgs_episode_thresholds(
        source_weight,
        scaled_lifetime,
        filter_l1,
        lp_constant,
        bernstein_constant,
        dominant_fraction=dominant_fraction,
        ancestry_alpha=ancestry_alpha,
    )
    Y = float(th["integrated_forced_square_service"])
    d = float(high_frequency_dissipation)
    oldcap = float(old_pool_integrated_capacity)
    old = float(old_old_integrated_service)
    xi = float(selected_interface_integrated_service)
    vals = (d, oldcap, old, xi)
    if any(v < 0 or not math.isfinite(v) for v in vals):
        raise ValueError("finite nonnegative aggregate service data required")
    w = np.asarray(tuple(float(x) for x in new_edge_integrated_services), dtype=float)
    if w.ndim != 1 or np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("finite nonnegative fresh-edge services required")

    if d >= Y / 4.0:
        return {
            "branch": "high_frequency_dissipation",
            "branch_value": d,
            "threshold": Y / 4.0,
            "resolved_DV_supplier": "NO",
            "next_owner_interface": HIGH_FREQUENCY_REENTRY_STATUS,
            "master_semantics": "RECURSE_CRITICAL_WITH_HIGH_FREQUENCY_OWNER",
        }

    low_lower = max(0.0, Y - 2.0 * d)
    total = old + xi + float(w.sum())
    tol = 2e-12 * max(1.0, Y, total)
    if total + tol < low_lower:
        raise ValueError("integrated edge law does not realize the forced low service")
    if oldcap > Y / 8.0:
        return {
            "branch": "old_pool_not_yet_eroded",
            "branch_value": oldcap,
            "threshold": Y / 8.0,
            "master_semantics": "RECURSE_CRITICAL_RESERVOIR_CAPACITY",
        }
    if old > oldcap + tol:
        raise ValueError("old-old service exceeds certified old-pool capacity")
    if xi >= Y / 8.0:
        return {
            "branch": "selected_interface_Xi",
            "branch_value": xi,
            "threshold": Y / 8.0,
            "master_semantics": "TRANSFER_COST",
        }

    fresh = float(w.sum())
    if fresh + tol < Y / 4.0:
        raise AssertionError("fresh integrated coherent-service lower failed")
    if fresh <= 0:
        raise AssertionError("positive fresh service required")
    p = w / fresh
    imax = int(np.argmax(p))
    pmax = float(p[imax])
    theta = float(dominant_fraction)
    alpha = float(ancestry_alpha)
    if pmax >= theta:
        edge = float(w[imax])
        pair_occ = edge / 2.0
        whole_occ = pair_occ / 2.0
        clean_pair = theta * Y / 8.0
        clean_whole = theta * Y / 16.0
        peak = whole_occ / float(scaled_lifetime)
        clean_peak = clean_whole / float(scaled_lifetime)
        if pair_occ + tol < clean_pair or whole_occ + tol < clean_whole or peak + tol < clean_peak:
            raise AssertionError("dominant fresh shell supplier lower failed")
        return {
            "branch": "dominant_fresh_critical_shell",
            "fresh_service": fresh,
            "dominant_edge_service": edge,
            "pair_mass_occupation": pair_occ,
            "whole_shell_mass_occupation": whole_occ,
            "peak_whole_shell_mass_lower": peak,
            "clean_peak_whole_shell_mass_lower": clean_peak,
            "canonical_renewal_fate": False,
            "canonical_replacement": FRESH_SCALE_REENTRY_STATUS,
            "master_semantics": "SIDELEDGER_ONLY__LEGACY_CELL_CLUSTER",
        }

    q = float(np.dot(p, p))
    h = -math.log(q)
    h0 = -math.log(theta)
    if h + 1e-13 < h0:
        raise AssertionError("integrated fresh-service collision entropy failed")
    if ancestry_labels is None:
        return {
            "branch": "fresh_service_collision_entropy",
            "H_atomic": h,
            "entropy_lower": h0,
            "canonical_renewal_fate": False,
            "canonical_replacement": FRESH_SCALE_REENTRY_STATUS,
            "master_semantics": "SIDELEDGER_ONLY__LEGACY_CELL_ENTROPY",
        }
    if len(ancestry_labels) != len(w):
        raise ValueError("ancestry label length mismatch")
    chain = collision_chain(p, ancestry_labels)
    if chain["h_ancestry"] >= alpha * h0 - 1e-13:
        return {
            "branch": "fresh_service_Bellman_entropy",
            "H_atomic": h,
            "H_ancestry": chain["h_ancestry"],
            "ancestry_entropy_lower": alpha * h0,
            "canonical_renewal_fate": False,
            "canonical_replacement": FRESH_SCALE_REENTRY_STATUS,
            "master_semantics": "SIDELEDGER_ONLY__LEGACY_CELL_ANCESTRY_ENTROPY",
        }
    pair_lower = theta**alpha - theta
    if chain["hidden_pair_mass"] + 2e-13 < pair_lower:
        raise AssertionError("integrated fresh-service ancestry pair bound failed")
    return {
        "branch": "fresh_service_same_ancestry_cycle",
        "H_atomic": h,
        "H_ancestry": chain["h_ancestry"],
        "hidden_pair_mass": chain["hidden_pair_mass"],
        "hidden_pair_lower": pair_lower,
        "canonical_renewal_fate": False,
        "canonical_replacement": FRESH_SCALE_REENTRY_STATUS,
        "master_semantics": "SIDELEDGER_ONLY__LEGACY_CELL_CYCLE",
    }



def objective_sgs_aggregate_scale_route(
    source_weight: float,
    scaled_lifetime: float,
    viscosity: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    *,
    block_frequency: float,
    high_frequency_dissipation: float,
    old_pool_integrated_capacity: float,
    old_old_integrated_service: float,
    selected_interface_integrated_service: float,
    fresh_band_integrated_services: Mapping[int, float],
) -> dict[str, object]:
    """Canonical realized objective-SGS owner set after quotienting cell refinement.

    The actual low-service law is validated first.  High-frequency dissipation,
    uneroded old-pool capacity, selected interface Xi, and the fresh scale law are
    then registered independently from the same physical data.  All owners whose
    physical thresholds are met are retained jointly; there is no theorem-name or
    branch-order priority.  On the fresh owner, only the pushforward to the fixed
    LP band index enters renewal, never a coherent-cell argmax.
    """
    th = objective_sgs_episode_thresholds(
        source_weight,
        scaled_lifetime,
        filter_l1,
        lp_constant,
        bernstein_constant,
    )
    Y = float(th["integrated_forced_square_service"])
    c = float(scaled_lifetime)
    N = float(block_frequency)
    nu = float(viscosity)
    d = float(high_frequency_dissipation)
    oldcap = float(old_pool_integrated_capacity)
    old = float(old_old_integrated_service)
    xi = float(selected_interface_integrated_service)
    if c <= 0 or N <= 0 or nu < 0 or any(v < 0 or not math.isfinite(v) for v in (d, oldcap, old, xi)):
        raise ValueError("valid lifetime/frequency/viscosity and finite nonnegative service data required")
    if not all(math.isfinite(x) for x in (c, N, nu)):
        raise ValueError("finite lifetime/frequency/viscosity required")
    fresh_law = {int(j): float(v) for j, v in fresh_band_integrated_services.items()}
    if any(j > 0 or v < 0 or not math.isfinite(v) for j, v in fresh_law.items()):
        raise ValueError("finite nonnegative fresh low/base band law j<=0 required")
    fresh = sum(fresh_law.values())

    low_lower = max(0.0, Y - 2.0 * d)
    total = old + xi + fresh
    service_tol = 2e-12 * max(1.0, Y, total, low_lower)
    if old > oldcap + service_tol:
        raise ValueError("old-old service exceeds certified old-pool capacity")
    if total + service_tol < low_lower:
        raise ValueError("band-pushforward service law does not realize the forced low service")

    owners: list[str] = []
    routes: dict[str, object] = {}

    if d >= Y / 4.0:
        name = "high_frequency_dissipation"
        owners.append(name)
        routes[name] = {
            "branch_value": d,
            "threshold": Y / 4.0,
            "resolved_DV_supplier": "NO",
            "next_owner_interface": HIGH_FREQUENCY_REENTRY_STATUS,
            "master_semantics": "RECURSE_CRITICAL_WITH_HIGH_FREQUENCY_OWNER",
        }

    if oldcap > Y / 8.0:
        name = "old_pool_not_yet_eroded"
        owners.append(name)
        routes[name] = {
            "branch_value": oldcap,
            "threshold": Y / 8.0,
            "master_semantics": "RECURSE_CRITICAL_RESERVOIR_CAPACITY",
        }

    if xi >= Y / 8.0:
        name = "selected_interface_Xi"
        owners.append(name)
        routes[name] = {
            "branch_value": xi,
            "threshold": Y / 8.0,
            "master_semantics": "TRANSFER_COST",
        }

    if fresh >= Y / 4.0:
        scale = fresh_service_scale_route(
            Y,
            c,
            N,
            fresh_law,
            viscosity=nu,
        )
        name = "fresh_scale_critical_shell"
        owners.append(name)
        routes[name] = {
            "fresh_service": fresh,
            "fresh_band_law": fresh_law,
            "scale_route": scale,
            "critical_shell_mass_lower": float(scale["hard_shell_mass_lower"]),
            "H_inf_scale": float(scale["H_inf_scale"]),
            "H2_scale": float(scale["H2_scale"]),
            "next_owner_interface": FRESH_SCALE_REENTRY_STATUS,
            "coherent_cell_argmax_used": False,
            "cell_ancestry_sideledger_optional": True,
            "master_semantics": "RECURSE_CRITICAL_VIA_REFINEMENT_INVARIANT_SCALE_SHELL",
        }

    if not owners:
        raise AssertionError("complete objective-SGS service law reached no physical owner")

    return {
        "integrated_forced_square_service": Y,
        "low_service_lower": low_lower,
        "realized_low_service": total,
        "fresh_service": fresh,
        "joint_primary_owners": tuple(owners),
        "routes": routes,
        "coherent_cell_priority_used": False,
        "master_semantics": "JOINT_NATIVE_OWNERS__NO_LEXICOGRAPHIC_PRIORITY",
    }

def objective_sgs_high_frequency_physical_reentry(
    high_frequency_dissipation: float,
    viscosity: float,
    inherited_scaled_tail_energy: float,
    positive_scaled_tail_work: float,
    *,
    lp_to_physical_tail_lower: float | None = None,
) -> dict[str, object]:
    """Route the SGS high-frequency service exit through its native tail-energy law.

    `D_high` remains the smooth-LP high-frequency normalized enstrophy.  It
    reaches the orthogonal hard-tail energy theorem only through a certified
    comparison `D_tail>=c_LP D_high`.  If no constant is supplied, use the
    canonical square-normalized smooth dyadic choice `c_LP=1/4`.
    """
    c_lp = (
        canonical_square_lp_tail_comparison_constant()
        if lp_to_physical_tail_lower is None
        else float(lp_to_physical_tail_lower)
    )
    D_tail_lower = physical_tail_dissipation_lower_from_lp(
        high_frequency_dissipation, c_lp
    )
    gate = classify_high_tail_energy_owners(
        D_tail_lower,
        viscosity,
        inherited_scaled_tail_energy,
        positive_scaled_tail_work,
    )
    return {
        "owner": "sgs_high_frequency_dissipation",
        "energy_gate": gate,
        "clean_thresholds": lp_high_clean_reentry_thresholds(
            high_frequency_dissipation, c_lp, viscosity
        ),
        "lp_to_physical_tail_lower": c_lp,
        "physical_tail_dissipation_lower": D_tail_lower,
        "next_theorem_status": HIGH_FREQUENCY_REENTRY_STATUS,
        "resolved_DV_supplier": False,
        "master_semantics": "RECURSE_CRITICAL",
    }


def pressure_source_alternatives(
    source_weight: float,
    scaled_lifetime: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> dict[str, float | str | bool]:
    """Legacy coarse pressure estimate retained as a diagnostic only.

    The bound `rho_P<=mu_V/5700+||R||_(3/2)/380` remains mathematically valid,
    but aggregate `mu_V` is no longer the canonical pressure renewal state.  The
    master-facing compiler uses `pressure_canonical_interface` and the realized
    positive SGS/pair source law instead.
    """
    sigma = float(source_weight)
    c = float(scaled_lifetime)
    if sigma <= 0 or c <= 0 or not math.isfinite(sigma + c):
        raise ValueError("positive finite pressure source weight/lifetime required")
    ratio = float(pressure_hessian_pair_energy_service_ratio_upper())
    effective_sgs = sigma / 2.0
    stress_service = objective_sgs_integrated_square_service_lower(
        effective_sgs, filter_l1, lp_constant, bernstein_constant
    )
    return {
        "owner": "pressure_diagnostic",
        "resolved_lowpass_mass_occupation_lower": 2850.0 * sigma,
        "resolved_lowpass_peak_mass_lower": 2850.0 * sigma / c,
        "stress_l32_occupation_lower": 190.0 * sigma,
        "effective_sgs_source_weight_if_stress_branch": effective_sgs,
        "integrated_square_service_if_stress_branch": stress_service,
        "fixed_pair_service_ratio_upper": ratio,
        "fixed_pair_total_future_multiplier_upper": 1.0 / (1.0 - ratio),
        "resolved_mass_is_generic_critical_shell": "NO",
        "canonical_pressure_route": False,
        "master_semantics": "DIAGNOSTIC_ONLY",
    }


def pressure_canonical_interface(
    source_weight: float,
    scaled_lifetime: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> dict[str, object]:
    """Symbolic master-facing interface for one pressure owner.

    The actual averaged pressure tensor is scalarized by its event Frobenius dual.
    Its positive source law has only two physical owners: SGS pressure source and
    the resolved unordered hard-pair law.  One carries at least `Sigma_P/2`; ties
    are joint.  The pair owner is already a generic critical-shell supplier through

      mu_child exp(H2_pair) >= 320 Sigma_P/c.

    No aggregate low-pass mass appears in this canonical interface.
    """
    sigma = float(source_weight)
    c = float(scaled_lifetime)
    if sigma <= 0 or c <= 0 or not math.isfinite(sigma + c):
        raise ValueError("positive finite pressure source weight/lifetime required")
    half = sigma / 2.0
    return {
        "owner": "pressure",
        "positive_source_owner_threshold": half,
        "sgs_stress_occupation_lower_if_sgs_owner": 380.0 * half,
        "effective_sgs_source_weight_if_sgs_owner": half,
        "integrated_square_service_lower_if_sgs_owner": objective_sgs_integrated_square_service_lower(
            half, filter_l1, lp_constant, bernstein_constant
        ),
        "pair_entropy_shell_tradeoff": "mu_child exp(H2_pair)>=320 Sigma_P/c",
        "pair_quarter_shell_corollary_lower": 80.0 * sigma / c,
        "pair_quarter_entropy_corollary_lower": math.log(4.0),
        "pair_theorem_status": PRESSURE_PAIR_STATUS,
        "realized_positive_source_law_required": True,
        "aggregate_muV_canonical_route": False,
        "coarse_muV_estimate_available_as_diagnostic": True,
        "master_semantics": "PRESSURE_SGS_SERVICE_OR_PAIR_CRITICAL_SHELL",
    }


def realized_pressure_source_route(
    source_weight: float,
    scaled_lifetime: float,
    viscosity: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    *,
    block_frequency: float,
    sgs_positive_source_weight: float,
    pair_positive_weights: Sequence[float],
    pair_shell_indices: Sequence[tuple[int, int]],
    pair_frequencies: Sequence[tuple[float, float]],
    pair_dominance_fraction: float = DEFAULT_DOMINANT_FRACTION,
) -> dict[str, object]:
    """Compile the **actual** positive pressure source law into native owners.

    This function does not use the coarse `mu_V` split.  It delegates the same
    Frobenius-dual positive SGS/pair law to the pressure-pair theorem.  If SGS is
    an owner, its actual positive source weight enters the direct coherent-service
    map.  If the pair law is an owner, it always exposes a genuine `u` hard shell;
    the existing generic shell theorem then supplies a conditional own-scale
    service lower on a full no-hit natural survivor.
    """
    sigma = float(source_weight)
    c = float(scaled_lifetime)
    nu = float(viscosity)
    if sigma <= 0 or c <= 0 or nu < 0 or not all(math.isfinite(x) for x in (sigma, c, nu)):
        raise ValueError("valid pressure source/lifetime/viscosity required")
    pair_route = objective_pressure_pair_route(
        sigma,
        c,
        block_frequency,
        sgs_positive_source_weight=sgs_positive_source_weight,
        pair_positive_weights=pair_positive_weights,
        pair_shell_indices=pair_shell_indices,
        pair_frequencies=pair_frequencies,
        dominant_fraction=pair_dominance_fraction,
    )
    owners = tuple(pair_route["joint_primary_owners"])
    routes: dict[str, object] = {}
    if "sgs_pressure_source" in owners:
        sgs_weight = float(sgs_positive_source_weight)
        routes["sgs_pressure_source"] = {
            "actual_positive_source_weight": sgs_weight,
            "stress_l32_occupation_lower": 380.0 * sgs_weight,
            "integrated_forced_square_service_lower": objective_sgs_integrated_square_service_lower(
                sgs_weight, filter_l1, lp_constant, bernstein_constant
            ),
            "next_owner": "coherent_service",
            "master_semantics": "RECURSE_CRITICAL_VIA_COHERENT_SERVICE",
        }
    if "resolved_pressure_pair_law" in owners:
        mu0 = float(pair_route["critical_shell_mass_lower"])
        h2 = float(pair_route["pair_source_entropy"]["H2_pair_source"])
        clean_mu = clean_entropy_shell_tradeoff_lower(sigma, c, h2)
        tol = 8e-13 * max(1.0, mu0, clean_mu)
        if mu0 + tol < clean_mu:
            raise AssertionError("compiler weakened the certified pressure entropy-shell lower")
        y_shell = critical_shell_bounded_service_lower(mu0, c, nu)
        s_shell = critical_shell_integrated_service_lower(mu0, c, nu)
        exp_h2 = math.exp(h2)
        weighted_mu = mu0 * exp_h2
        weighted_y = y_shell * exp_h2
        weighted_s = s_shell * exp_h2
        clean_weighted_mu = 320.0 * sigma / c
        clean_weighted_y = critical_shell_bounded_service_lower(clean_weighted_mu, c, nu)
        clean_weighted_s = critical_shell_integrated_service_lower(clean_weighted_mu, c, nu)
        service_tol = 1e-12 * max(1.0, weighted_y, weighted_s, clean_weighted_y, clean_weighted_s)
        if weighted_mu + service_tol < clean_weighted_mu or weighted_y + service_tol < clean_weighted_y or weighted_s + service_tol < clean_weighted_s:
            raise AssertionError("pressure compiler lost service-complexity conjugacy")
        routes["resolved_pressure_pair_law"] = {
            "pair_source_entropy": h2,
            "critical_shell_mass_lower": mu0,
            "entropy_weighted_critical_shell_mass_lower": weighted_mu,
            "clean_entropy_weighted_critical_shell_mass_lower": clean_weighted_mu,
            "full_survivor_own_scale_service_lower": y_shell,
            "full_survivor_integrated_service_lower": s_shell,
            "entropy_weighted_full_survivor_service_lower": weighted_y,
            "entropy_weighted_full_survivor_integrated_service_lower": weighted_s,
            "clean_entropy_weighted_full_survivor_service_lower": clean_weighted_y,
            "clean_entropy_weighted_full_survivor_integrated_service_lower": clean_weighted_s,
            "next_owner": "generic_critical_shell_first_stop",
            "master_semantics": "RECURSE_CRITICAL_VIA_GENERIC_SHELL",
        }
    if not routes:
        raise AssertionError("pressure positive source law produced no physical owner")
    return {
        "owner": "pressure",
        "joint_primary_owners": owners,
        "routes": routes,
        "pressure_pair_source_law": pair_route,
        "aggregate_muV_used": False,
        "additive_reset_created": False,
        "master_semantics": "JOINT_RECURSIVE_PRESSURE_OWNERS",
    }


def compile_objective_source_owners(
    objective_variation_action: float,
    scaled_lifetime: float,
    owner_weights: Mapping[str, float],
    *,
    viscosity: float,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
    aspect: float = EXTENDED_ASPECT,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
) -> dict[str, object]:
    """Non-lexicographic source-owner compiler for one coherent objective event.

    All owner weights are positive norm/source measures in the same scaled units.
    Every owner meeting the clean pigeonhole threshold is returned; simultaneous
    causes remain a joint set. The compiler never promotes critical D_V to a finite
    reset and never invents a packet synchronization interface.
    """
    if set(owner_weights) != set(OWNER_NAMES):
        raise ValueError(f"owner weights must be exactly {OWNER_NAMES}")
    w = {k: float(owner_weights[k]) for k in OWNER_NAMES}
    if any(x < 0 or not math.isfinite(x) for x in w.values()):
        raise ValueError("finite nonnegative owner weights required")
    A = float(objective_variation_action)
    c = float(scaled_lifetime)
    sigma_star = objective_owner_weight_threshold(A, c)
    required_sum = A / c
    tol = 8e-13 * max(1.0, required_sum, sum(w.values()))
    if sum(w.values()) + tol < required_sum:
        raise ValueError("owner weights do not cover the objective source action")
    qualifying = tuple(k for k in OWNER_NAMES if w[k] + tol >= sigma_star)
    if not qualifying:
        raise AssertionError("four-owner source pigeonhole failed")

    routes: dict[str, object] = {}
    for k in qualifying:
        if k == "local_dv":
            routes[k] = local_dv_reentry(w[k], c, viscosity, aspect=aspect, lower_axis_constant=lower_axis_constant)
        elif k == "viscous":
            routes[k] = viscous_dv_reentry(w[k], c, viscosity)
        elif k == "sgs":
            routes[k] = objective_sgs_episode_thresholds(w[k], c, filter_l1, lp_constant, bernstein_constant)
        elif k == "pressure":
            routes[k] = pressure_canonical_interface(w[k], c, filter_l1, lp_constant, bernstein_constant)
    return {
        "objective_variation_action": A,
        "scaled_lifetime": c,
        "owner_threshold": sigma_star,
        "joint_owners": qualifying,
        "routes": routes,
        "additive_reset_created": False,
        "packet_synchronization_created": False,
        "status": STATUS,
    }


def theorem_certificate() -> dict[str, object]:
    ratio = pressure_hessian_pair_energy_service_ratio_upper()
    if not ratio < 1 / 5:
        raise AssertionError("objective pressure-Hessian fixed-pair service ratio lost its one-fifth life")
    return {
        "status": STATUS,
        "owner_split": "local_DV / pressure / SGS / viscosity; all threshold ties retained jointly",
        "local_owner": "Sigma_local<=C_local D_V -> generic critical-shell reentry",
        "viscous_owner": "int rho_nu>=Sigma -> D_V>=(1500 Sigma/nu)^2/c -> generic critical-shell reentry",
        "sgs_owner": "rho_R -> ||R||_(3/2) -> Q^(3/2) -> Y^(2/3), giving exact linear C_Y rho_R",
        "sgs_clean_route": "D_high>=Y/4 OR oldcap>Y/8 OR Xi>=Y/8 OR fresh NN band pushforward -> hard critical shell",
        "sgs_fresh_scale_route": "fresh F>=Y/4, p_j=F_j/F on the fixed LP band index -> mu_hard exp(H_inf_scale)>=Y/(24c); no coherent-cell dominance needed for renewal",
        "sgs_cell_sideledger": "coherent-cell dominance/entropy/ancestry cycle remains optional fine accounting only and cannot change the canonical fresh renewal fate",
        "sgs_joint_owner_rule": "high-tail / old-capacity / Xi / fresh-scale conditions are read independently from the realized law and all satisfied owners are retained jointly; no branch-order priority",
        "sgs_high_frequency_owner": "smooth-LP D_high enters the orthogonal hard-tail energy theorem only through a certified D_tail>=c_LP D_high comparison, then routes to inherited critical shell OR actual positive HH/resolved-interface regeneration; never resolved D_V",
        "pressure_owner": "actual Frobenius-dual positive source law: SGS>=Sigma_P/2 OR resolved unordered pair law>=Sigma_P/2; ties joint",
        "pressure_pair_route": "every resolved pair owner satisfies mu_child exp(H2_pair)>=320 Sigma_P/c and enters generic critical-shell reentry",
        "pressure_service_conjugacy": "on a full no-hit shell corridor, exp(H2_pair) times own-scale/integrated service is at least the generic-shell service generated by mass 320 Sigma_P/c",
        "pressure_sgs_route": "actual positive SGS pressure weight r gives int||R||_(3/2)>=380 r and direct coherent service",
        "pressure_legacy_muV": "rho_P<=mu_V/5700+||R||_(3/2)/380 remains diagnostic only and is absent from the canonical compiler state",
        "pressure_hessian_pair_ratio": f"{ratio.numerator}/{ratio.denominator}<1/5",
        "pressure_hessian_total_future_pair_capacity": "<5/4 generation-0 pair capacity on supplied signed-good low-strain lineage; optional material-reuse refinement only",
        "forbidden_identifications": (
            "aggregate pressure mu_V is not a canonical renewal state; pressure H2 is not a causal child-energy probability; fresh-scale H_inf/H2 are not causal child-energy probabilities; coherent-cell entropy is not a canonical fresh renewal fate; high-frequency SGS dissipation is not resolved D_V"
        ),
        "master_rule": "all D_V/shell/source/service outputs remain recursive scale-critical owners; no additive finite reset is created",
    }


@dataclass(frozen=True)
class ObjectiveSourceCompilerStress:
    samples: int
    worst_sgs_closed_form_relative_residual: float
    minimum_owner_pigeonhole_margin: float
    minimum_pressure_diagnostic_split_identity_margin: float
    minimum_pressure_pair_entropy_shell_margin: float
    minimum_pressure_pair_full_survivor_service_margin: float
    minimum_fresh_scale_shell_margin: float
    minimum_fresh_scale_service_conjugacy_margin: float
    minimum_local_dv_identity_margin: float
    minimum_viscous_cauchy_identity_margin: float
    maximum_joint_owner_count: int
    maximum_joint_pressure_owner_count: int
    fresh_cell_argmax_regressions: int


def stress(samples: int = 50_000, seed: int = 20260809) -> ObjectiveSourceCompilerStress:
    rng = np.random.default_rng(seed)
    ws = 0.0
    mo = mp = mpp = mps = mfs = mfss = ml = mv = float("inf")
    max_joint = max_pressure_joint = 0
    fresh_cell_argmax_regressions = 0
    for _ in range(samples):
        g1 = float(rng.uniform(1.0, 3.0))
        clp = float(rng.uniform(0.8, 3.0))
        cb = float(rng.uniform(0.8, 2.5))
        rho = float(10.0 ** rng.uniform(-8.0, 1.0))
        Cinc = sgs_increment_cubic_upper(g1, 1.0)
        rnorm = sgs_gradient_stress_lower(rho)
        q = rnorm**1.5 / Cinc
        y_comp = cubic_to_square_threshold(q, g1, clp, cb)
        y_closed = objective_sgs_square_service_per_source(g1, clp, cb) * rho
        ws = max(ws, abs(y_comp - y_closed) / max(1e-300, abs(y_closed)))
        if abs(y_comp - y_closed) > 5e-12 * max(1.0, abs(y_closed)):
            raise AssertionError("objective SGS 3/2-to-2/3 cancellation failed")

        c = float(rng.uniform(0.04, 2.0))
        nu = float(rng.uniform(0.05, 3.0))
        weights = rng.lognormal(mean=-3.0, sigma=1.2, size=4)
        total = float(weights.sum())
        A = c * total * float(rng.uniform(0.2, 1.0))
        owner_map = dict(zip(OWNER_NAMES, (float(x) for x in weights)))
        out = compile_objective_source_owners(
            A,
            c,
            owner_map,
            viscosity=nu,
            filter_l1=g1,
            lp_constant=clp,
            bernstein_constant=cb,
        )
        sigma_star = float(out["owner_threshold"])
        qset = tuple(out["joint_owners"])
        max_joint = max(max_joint, len(qset))
        best = max(owner_map.values())
        mo = min(mo, best - sigma_star)
        if not qset or best + 2e-12 < sigma_star:
            raise AssertionError("source owner pigeonhole failed")

        sigma = float(rng.lognormal(-3.0, 1.0))
        # Legacy coarse estimate remains algebraically correct but diagnostic only.
        pdiag = pressure_source_alternatives(sigma, c, g1, clp, cb)
        pmass = float(pdiag["resolved_lowpass_mass_occupation_lower"])
        pstress = float(pdiag["stress_l32_occupation_lower"])
        mp = min(mp, pmass / 5700.0 - sigma / 2.0, pstress / 380.0 - sigma / 2.0)
        if pdiag["canonical_pressure_route"] is not False:
            raise AssertionError("legacy pressure mu_V estimate re-entered the canonical compiler")
        if abs(pmass / 5700.0 - sigma / 2.0) > 1e-12 * max(1.0, sigma):
            raise AssertionError("pressure diagnostic mass split identity failed")
        if abs(pstress / 380.0 - sigma / 2.0) > 1e-12 * max(1.0, sigma):
            raise AssertionError("pressure diagnostic stress split identity failed")

        # Realized positive pressure source law: exact SGS/pair cover with one
        # physical frequency per hard shell label.
        frac = float(rng.uniform(0.0, 1.0))
        sgs_w = frac * sigma
        pair_total = (1.0 - frac) * sigma
        npair = int(rng.integers(1, 7))
        raw = rng.random(npair)
        raw /= float(raw.sum())
        pair_w = (pair_total * raw).tolist()
        Np = float(math.exp(rng.uniform(-1.0, 4.0)))
        indices = [(j, j) for j in range(npair)]
        freqs = [(Np / (4.0 * 2.0**j), Np / (4.0 * 2.0**j)) for j in range(npair)]
        pr = realized_pressure_source_route(
            sigma,
            c,
            nu,
            g1,
            clp,
            cb,
            block_frequency=Np,
            sgs_positive_source_weight=sgs_w,
            pair_positive_weights=pair_w,
            pair_shell_indices=indices,
            pair_frequencies=freqs,
        )
        powners = tuple(pr["joint_primary_owners"])
        max_pressure_joint = max(max_pressure_joint, len(powners))
        if not powners or pr["aggregate_muV_used"]:
            raise AssertionError("canonical pressure law lost its native positive owners")
        if "resolved_pressure_pair_law" in powners:
            rr = pr["routes"]["resolved_pressure_pair_law"]
            h2 = float(rr["pair_source_entropy"])
            mu = float(rr["critical_shell_mass_lower"])
            trade = 320.0 * sigma / c
            mpp = min(mpp, mu * math.exp(h2) - trade)
            if mu * math.exp(h2) + 3e-11 * max(1.0, trade) < trade:
                raise AssertionError("pressure compiler lost entropy-shell tradeoff")
            weighted_serv = float(rr["entropy_weighted_full_survivor_integrated_service_lower"])
            clean_weighted_serv = float(rr["clean_entropy_weighted_full_survivor_integrated_service_lower"])
            mps = min(mps, weighted_serv - clean_weighted_serv)
            if weighted_serv + 3e-12 * max(1.0, clean_weighted_serv) < clean_weighted_serv:
                raise AssertionError("pressure compiler lost entropy-weighted full-survivor service lower")

        # Canonical fresh SGS route: construct a fresh band law at exactly the
        # certified Y/4 threshold after keeping high/old/Xi below their faces.
        sgs_sigma = float(rng.lognormal(-3.0, 1.0))
        sth = objective_sgs_episode_thresholds(sgs_sigma, c, g1, clp, cb)
        Ysgs = float(sth["integrated_forced_square_service"])
        nband = int(rng.integers(1, 9))
        raw_band = rng.random(nband)
        raw_band /= float(raw_band.sum())
        fresh_total = Ysgs * float(rng.uniform(0.25, 0.55))
        labels_band = list(range(-nband + 1, 1))
        band_law = {j: float(p * fresh_total) for j, p in zip(labels_band, raw_band)}
        # Choose the other low-service owners so the aggregate cover is realized
        # while staying strictly below their stop thresholds.
        Dhigh = float(rng.uniform(0.0, 0.24)) * Ysgs
        oldcap_s = float(rng.uniform(0.0, 0.12)) * Ysgs
        old_s = float(rng.uniform(0.0, 1.0)) * oldcap_s
        xi_s = float(rng.uniform(0.0, 0.12)) * Ysgs
        low_need = max(0.0, Ysgs - 2.0 * Dhigh)
        current = old_s + xi_s + fresh_total
        if current < low_need:
            # Put exactly the missing physical low service into the fresh band law;
            # scaling all band atoms preserves its concentration coordinates.
            extra = low_need - current
            factor = (fresh_total + extra) / fresh_total
            band_law = {j: v * factor for j, v in band_law.items()}
            fresh_total += extra
        Nsgs = float(math.exp(rng.uniform(-1.0, 4.0)))
        sr = objective_sgs_aggregate_scale_route(
            sgs_sigma,
            c,
            nu,
            g1,
            clp,
            cb,
            block_frequency=Nsgs,
            high_frequency_dissipation=Dhigh,
            old_pool_integrated_capacity=oldcap_s,
            old_old_integrated_service=old_s,
            selected_interface_integrated_service=xi_s,
            fresh_band_integrated_services=band_law,
        )
        if tuple(sr["joint_primary_owners"]) != ("fresh_scale_critical_shell",):
            raise AssertionError("canonical fresh SGS stress state did not retain the unique fresh-scale owner")
        fresh_route = sr["routes"]["fresh_scale_critical_shell"]
        if fresh_route["coherent_cell_argmax_used"] or sr["coherent_cell_priority_used"]:
            fresh_cell_argmax_regressions += 1
            raise AssertionError("canonical fresh SGS route regressed to coherent-cell selection")
        fs = fresh_route["scale_route"]
        pmax_scale = float(fs["p_max"])
        mu_scale = float(fs["hard_shell_mass_lower"])
        clean_scale = pmax_scale * Ysgs / (24.0 * c)
        mfs = min(mfs, mu_scale - clean_scale)
        if mu_scale + 3e-12 * max(1.0, clean_scale) < clean_scale:
            raise AssertionError("compiler lost fresh scale hard-shell lower")
        weighted_s = float(fs["H_inf_weighted_full_survivor_integrated_service_lower"])
        clean_weighted_s = float(fs["clean_H_inf_weighted_full_survivor_integrated_service_lower"])
        mfss = min(mfss, weighted_s - clean_weighted_s)
        if weighted_s + 3e-12 * max(1.0, clean_weighted_s) < clean_weighted_s:
            raise AssertionError("compiler lost fresh scale service-concentration conjugacy")

        loc = local_dv_reentry(sigma, c, nu)
        C_local = float(loc["local_source_per_DV"])
        Dloc = float(loc["resolved_DV_lower"])
        ml = min(ml, C_local * Dloc - sigma)
        if abs(C_local * Dloc - sigma) > 2e-12 * max(1.0, sigma):
            raise AssertionError("local source-to-DV identity failed")

        vis = viscous_dv_reentry(sigma, c, nu)
        Dv = float(vis["resolved_DV_lower"])
        expect = (1500.0 * sigma / nu) ** 2 / c
        mv = min(mv, Dv - expect)
        if abs(Dv - expect) > 2e-12 * max(1.0, expect):
            raise AssertionError("viscous integrated Cauchy route failed")

    if not math.isfinite(mpp):
        mpp = 0.0
    if not math.isfinite(mps):
        mps = 0.0
    if not math.isfinite(mfs):
        mfs = 0.0
    if not math.isfinite(mfss):
        mfss = 0.0
    return ObjectiveSourceCompilerStress(
        samples, ws, mo, mp, mpp, mps, mfs, mfss, ml, mv,
        max_joint, max_pressure_joint, fresh_cell_argmax_regressions
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-objective-source-routing-compiler"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "objective_source_routing_compiler.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Coherent objective-source owner compiler

Status: **{cert['status']}**.

The source branch of the coherent service-or-flat gate does not need a new packet object.  Group the exact coherent averaged strain source by its physical owner:

`local_DV + pressure + SGS + viscosity`.

If `A_obj` is the objective variation action on a scaled lifetime `c`, their positive scaled weights satisfy `sum Sigma_r >= A_obj/c`, so at least one owner has `Sigma_r>=A_obj/(4c)`.  Exact ties are retained jointly; there is no theorem-name priority.

The local coherent quadratic/Reynolds owner is already bounded by `C_local D_V`, and the viscous owner obeys `rho_nu<=nu sqrt(d_V)/1500`; both therefore feed the generic critical-shell-to-own-scale-service theorem as **recursive** resolved dissipation, never as an additive reset.

For the objective SGS owner, the clean collision `||R||_(3/2)>=380 rho_R`, Germano `3/2` power and coherent-service `2/3` power cancel exactly:

`Y_R >= C_Y rho_R`,

`C_Y = 380/[g1(1+g1)(C_LP C_B)^2]`.

Thus integrated SGS source weight produces integrated coherent square service with no persistence hypothesis and no affine-radius packet.  High-frequency dissipation, old-pool capacity and selected-interface `Xi` keep their existing physical owners.  On the fresh NN branch, the canonical compiler now **first quotients coherent-cell refinement** and pushes the actual positive service measure to the fixed LP band index.  If `F_j` are the fresh band weights and `F>=Y/4`, the certified scale law gives

`mu_hard exp(H_inf^scale) >= Y/(24c)`,

so every fresh law enters the generic critical-shell first-stop theorem without a coherent-cell argmax.  Cell dominance/entropy/cycle remains optional ancestry sideledger only.  High-frequency dissipation is **not** renamed resolved `D_V`.

For pressure, the coarse estimate `rho_P<=mu_V/5700+||R||_(3/2)/380` is retained only as a diagnostic.  The canonical compiler now uses the actual Frobenius-dual positive source law

`rho_P <= [r_SGS]_+ + sum_(a<=b)[p_ab]_+`.

Thus positive SGS pressure source or the resolved unordered pair law carries at least `Sigma_P/2`, with exact ties joint.  The SGS owner uses its **actual** positive source weight `r` and gives `int||R||_(3/2)>=380r`, hence direct coherent service.  Every resolved pair owner satisfies

`mu_child exp(H2_pair) >= 320 Sigma_P/c`

and therefore enters the generic critical-shell first-stop theorem.  On a full no-hit natural survivor the compiler records the corresponding own-scale service, but it does not turn that conditional service into an unconditional event.  `H2_pair` is only the logarithmic weakening of the shell seed; it is neither a causal child-energy probability nor a separate stop.

The old fixed-material-pair ratio `{pressure_hessian_pair_energy_service_ratio_upper().numerator}/{pressure_hessian_pair_energy_service_ratio_upper().denominator}<1/5` remains a valid optional reuse refinement after material sidecars are attached; it is no longer the pressure renewal entrance.  The separate H1 pressure-third source retains its own `<1/3` theorem.

Stress: `{out.samples}` source-owner states
- worst SGS closed-form relative residual: `{out.worst_sgs_closed_form_relative_residual:.3e}`
- minimum owner-pigeonhole margin: `{out.minimum_owner_pigeonhole_margin:.3e}`
- minimum pressure diagnostic split margin: `{out.minimum_pressure_diagnostic_split_identity_margin:.3e}`
- minimum pressure entropy-shell margin: `{out.minimum_pressure_pair_entropy_shell_margin:.3e}`
- minimum pressure full-survivor service registration margin: `{out.minimum_pressure_pair_full_survivor_service_margin:.3e}`
- minimum fresh-scale shell margin: `{out.minimum_fresh_scale_shell_margin:.3e}`
- minimum fresh-scale full-survivor service-conjugacy margin: `{out.minimum_fresh_scale_service_conjugacy_margin:.3e}`
- minimum local-DV identity margin: `{out.minimum_local_dv_identity_margin:.3e}`
- minimum viscous-Cauchy identity margin: `{out.minimum_viscous_cauchy_identity_margin:.3e}`
- maximum sampled joint owner count: `{out.maximum_joint_owner_count}`
- maximum sampled joint pressure owner count: `{out.maximum_joint_pressure_owner_count}`
- fresh coherent-cell argmax regressions: `{out.fresh_cell_argmax_regressions}`

The resulting architecture is source-native: `local/viscous -> D_V -> critical shell`, `SGS -> coherent service -> high-tail / old-pool / Xi / refinement-invariant fresh scale shell`, `pressure -> actual SGS service OR entropy-weighted critical shell`.  No packet synchronization theorem and no uniform finite resource are inserted.  Final continuum master assembly and supplier-specific signed-good scale geometry remain separate.  No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
