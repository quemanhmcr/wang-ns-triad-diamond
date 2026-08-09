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
from src.high_frequency_dissipation_reentry import (
    STATUS as HIGH_FREQUENCY_REENTRY_STATUS,
    classify_high_tail_energy_owners,
    high_tail_clean_reentry_thresholds,
)

from src.critical_shell_service_reentry import (
    critical_shell_bounded_service_lower,
    dissipation_supplier_shell_mass_threshold,
)
from src.pressure_reservoir_sync import pair_energy_service_ratio_upper
from src.resolved_objective_strain_collision import sgs_gradient_stress_lower

STATUS = (
    "EXACT_COHERENT_OBJECTIVE_SOURCE_OWNER_COMPILER__LOCAL_DV_AND_VISCOSITY_TO_CRITICAL_SHELL__"
    "SGS_TO_COHERENT_SERVICE__PRESSURE_TO_SGS_OR_RESERVOIR__NO_PACKET_SYNCHRONIZATION"
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
    A theta-dominant fresh edge gives pair critical-mass occupation >=theta Y_tot/8,
    whole-shell occupation >=theta Y_tot/16 and therefore a pointwise shell event
    >=theta Y_tot/(16c). Otherwise the existing collision chain pays ancestry
    entropy or same-ancestry pair/cycle mass.
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
        "dominant_pair_mass_occupation_lower": theta * Y / 8.0,
        "dominant_whole_shell_mass_occupation_lower": theta * Y / 16.0,
        "dominant_peak_whole_shell_mass_lower": theta * Y / (16.0 * c),
        "atomic_entropy_lower": h_atomic,
        "ancestry_entropy_lower": h_anc,
        "same_ancestry_pair_mass_lower": pair,
        "high_frequency_dissipation_is_resolved_DV": "NO",
        "master_semantics": "RECURSE_CRITICAL / TRANSFER_COST according to realized service branch",
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
) -> dict[str, float | str]:
    """Route one realized integrated objective-SGS positive service law.

    This is the aggregate analogue of coherent_service_route. It preserves native
    owners and explicitly refuses the false conversion D_high -> resolved D_V.
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
            "master_semantics": "RECURSE_CRITICAL_VIA_GENERIC_SHELL",
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
            "master_semantics": "TRANSFER_COST",
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
            "master_semantics": "TRANSFER_COST",
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
        "master_semantics": "TRANSFER_COST",
    }


def objective_sgs_high_frequency_physical_reentry(
    high_frequency_dissipation: float,
    viscosity: float,
    inherited_scaled_tail_energy: float,
    positive_scaled_tail_work: float,
) -> dict[str, object]:
    """Route the SGS high-frequency service exit through its native tail-energy law.

    `D_high` remains high-frequency normalized enstrophy.  The companion theorem
    uses hard-tail Navier--Stokes energy to force inherited critical shell energy
    or actual positive nonlinear regeneration work; it never relabels the input
    resolved `D_V`.
    """
    gate = classify_high_tail_energy_owners(
        high_frequency_dissipation,
        viscosity,
        inherited_scaled_tail_energy,
        positive_scaled_tail_work,
    )
    return {
        "owner": "sgs_high_frequency_dissipation",
        "energy_gate": gate,
        "clean_thresholds": high_tail_clean_reentry_thresholds(
            high_frequency_dissipation, viscosity
        ),
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
) -> dict[str, float | str]:
    """Integrated filtered-pressure owner -> low-pass reservoir OR SGS service.

    From rho_P <= mu_V/5700 + ||R||_(3/2)/380, a pressure owner weight Sigma_P
    forces at least one of

      int mu_V >= 2850 Sigma_P,
      int ||R||_(3/2) >= 190 Sigma_P.

    The stress branch is exactly an effective objective-SGS source weight Sigma_P/2.
    The low-pass mass branch is deliberately NOT converted into a critical shell.
    It remains a pressure-reservoir occupation; on a supplied signed-good low-strain
    lineage each fixed material pair has service ratio <1/3 and total future capacity
    <3/2 of its generation-zero coefficient, so persistence requires pair relink,
    component entropy/cycle, leaving low strain, or the SGS branch.
    """
    sigma = float(source_weight)
    c = float(scaled_lifetime)
    if sigma <= 0 or c <= 0 or not math.isfinite(sigma + c):
        raise ValueError("positive finite pressure source weight/lifetime required")
    ratio = float(pair_energy_service_ratio_upper())
    effective_sgs = sigma / 2.0
    stress_service = objective_sgs_integrated_square_service_lower(
        effective_sgs, filter_l1, lp_constant, bernstein_constant
    )
    return {
        "owner": "pressure",
        "resolved_lowpass_mass_occupation_lower": 2850.0 * sigma,
        "resolved_lowpass_peak_mass_lower": 2850.0 * sigma / c,
        "stress_l32_occupation_lower": 190.0 * sigma,
        "effective_sgs_source_weight_if_stress_branch": effective_sgs,
        "integrated_square_service_if_stress_branch": stress_service,
        "fixed_pair_service_ratio_upper": ratio,
        "fixed_pair_total_future_multiplier_upper": 1.0 / (1.0 - ratio),
        "resolved_mass_is_generic_critical_shell": "NO",
        "master_semantics": "PRESSURE_RESERVOIR_OR_SGS_SERVICE",
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
            routes[k] = pressure_source_alternatives(w[k], c, filter_l1, lp_constant, bernstein_constant)
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
    ratio = pair_energy_service_ratio_upper()
    if not ratio < 1 / 3:
        raise AssertionError("pressure fixed-pair service ratio lost its one-third life")
    return {
        "status": STATUS,
        "owner_split": "local_DV / pressure / SGS / viscosity; all threshold ties retained jointly",
        "local_owner": "Sigma_local<=C_local D_V -> generic critical-shell reentry",
        "viscous_owner": "int rho_nu>=Sigma -> D_V>=(1500 Sigma/nu)^2/c -> generic critical-shell reentry",
        "sgs_owner": "rho_R -> ||R||_(3/2) -> Q^(3/2) -> Y^(2/3), giving exact linear C_Y rho_R",
        "sgs_clean_route": "D_high>=Y/4 OR oldcap>Y/8 OR Xi>=Y/8 OR fresh shell/entropy/cycle",
        "sgs_high_frequency_owner": "D_high is handed to the hard-tail energy theorem: inherited critical shell OR actual positive HH/resolved-interface regeneration; never resolved D_V",
        "pressure_owner": "int mu_V>=2850 Sigma_P OR effective SGS weight>=Sigma_P/2",
        "pressure_pair_ratio": f"{ratio.numerator}/{ratio.denominator}<1/3",
        "forbidden_identifications": (
            "pressure low-pass mass is not generic critical-shell mass; high-frequency SGS dissipation is not resolved D_V"
        ),
        "master_rule": "all D_V/shell/source/reservoir outputs remain recursive scale-critical owners; no additive finite reset is created",
    }


@dataclass(frozen=True)
class ObjectiveSourceCompilerStress:
    samples: int
    worst_sgs_closed_form_relative_residual: float
    minimum_owner_pigeonhole_margin: float
    minimum_pressure_split_identity_margin: float
    minimum_local_dv_identity_margin: float
    minimum_viscous_cauchy_identity_margin: float
    maximum_joint_owner_count: int


def stress(samples: int = 50_000, seed: int = 20260809) -> ObjectiveSourceCompilerStress:
    rng = np.random.default_rng(seed)
    ws = 0.0
    mo = mp = ml = mv = float("inf")
    max_joint = 0
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
        p = pressure_source_alternatives(sigma, c, g1, clp, cb)
        pmass = float(p["resolved_lowpass_mass_occupation_lower"])
        pstress = float(p["stress_l32_occupation_lower"])
        mp = min(mp, pmass / 5700.0 - sigma / 2.0, pstress / 380.0 - sigma / 2.0)
        if abs(pmass / 5700.0 - sigma / 2.0) > 1e-12 * max(1.0, sigma):
            raise AssertionError("pressure mass split identity failed")
        if abs(pstress / 380.0 - sigma / 2.0) > 1e-12 * max(1.0, sigma):
            raise AssertionError("pressure stress split identity failed")

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

    return ObjectiveSourceCompilerStress(samples, ws, mo, mp, ml, mv, max_joint)


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

Thus integrated SGS source weight produces integrated coherent square service with no persistence hypothesis and no affine-radius packet.  The positive service law routes exactly to high-frequency dissipation, old-pool capacity, selected-interface `Xi`, a dominant fresh critical shell, ancestry entropy, or same-ancestry cycle.  High-frequency dissipation is **not** renamed resolved `D_V`.

For pressure,

`rho_P <= mu_V/5700 + ||R||_(3/2)/380`

integrates to the honest alternative

`int mu_V >= 2850 Sigma_P`  OR  `int ||R||_(3/2) >= 190 Sigma_P`.

The stress alternative is exactly an effective SGS source weight `Sigma_P/2` and enters the same coherent-service route.  The resolved low-pass mass alternative is **not** promoted to a generic critical shell.  It remains pressure-reservoir occupation.  On a supplied signed-good low-strain lineage, each fixed materially reused low-low pair has pressure-service ratio `{pair_energy_service_ratio_upper().numerator}/{pair_energy_service_ratio_upper().denominator}<1/3`, so persistent pressure service must relink pairs, fragment into component entropy/cycle, leave low strain, or use the SGS branch.

Stress: `{out.samples}` source-owner states
- worst SGS closed-form relative residual: `{out.worst_sgs_closed_form_relative_residual:.3e}`
- minimum owner-pigeonhole margin: `{out.minimum_owner_pigeonhole_margin:.3e}`
- minimum pressure split identity margin: `{out.minimum_pressure_split_identity_margin:.3e}`
- minimum local-DV identity margin: `{out.minimum_local_dv_identity_margin:.3e}`
- minimum viscous-Cauchy identity margin: `{out.minimum_viscous_cauchy_identity_margin:.3e}`
- maximum sampled joint owner count: `{out.maximum_joint_owner_count}`

The resulting architecture is source-native: `local/viscous -> D_V -> critical shell`, `SGS -> coherent service`, `pressure -> SGS service or low-frequency reservoir`.  No packet synchronization theorem and no uniform finite resource are inserted.  Final continuum master assembly and supplier-specific signed-good scale geometry remain separate.  No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
