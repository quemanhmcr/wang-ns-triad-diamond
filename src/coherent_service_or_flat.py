from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.coherent_affine_projection import (
    EXTENDED_ASPECT,
    coherent_deformation_to_dissipation_constant,
    gaussian_core_nonaffine_forcing_upper,
    intrinsic_carrier_upper,
    normalized_dissipation_from_coherent_deformation,
)
from src.high_strain_dissipation_collision import normalized_dissipation_lower as high_strain_dissipation_lower
from src.service_or_flat_rigidity import (
    CURVATURE_DENOM,
    LOW_STRAIN_ACTION,
    PHASE_HOLONOMY_FLAT,
    affine_strain_flatness_upper,
    hodge_flatness_upper,
)

DEFAULT_RADIUS_CAP = 4.0
DEFAULT_SHELL_CARRIER_RATIO = math.exp(2.0 / 25.0)


def coherent_nonaffine_coefficient(qmax: float) -> float:
    if qmax < 0 or not math.isfinite(qmax):
        raise ValueError("finite nonnegative qmax required")
    return 1.0 + qmax / math.sqrt(2.0) + math.sqrt(7.0) / 2.0


def coherent_flat_thresholds(
    tau: float,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    carrier_ratio: float = DEFAULT_SHELL_CARRIER_RATIO,
    aspect_cap: float = EXTENDED_ASPECT,
) -> dict[str, float]:
    """Thresholds for the whole-eddy coherent service-or-flat gate.

    The transfer threshold is kept identical to the existing conservative packet
    assembly so Hodge/strain bounds are strictly stronger than tau/3.  The third
    flatness channel is now the complete non-affine Gaussian-core forcing action
    C(qmax) I_K, not a center-Hessian H1/H3 persistence quantity.
    """
    if not (0.0 < tau <= 0.1) or radius_cap <= 0 or carrier_ratio <= 0 or aspect_cap < 1:
        raise ValueError("invalid coherent flat threshold data")
    delta = tau * tau / CURVATURE_DENOM
    qmax = intrinsic_carrier_upper(radius_cap, aspect_cap, carrier_ratio)
    cnon = coherent_nonaffine_coefficient(qmax)
    return {
        "flatness_target": tau,
        "block_transfer_deficit": delta,
        "objective_strain_variation_action": tau / 60.0,
        "low_strain_action": float(LOW_STRAIN_ACTION),
        "aspect_threshold": aspect_cap,
        "radius_cap": radius_cap,
        "carrier_ratio": carrier_ratio,
        "intrinsic_carrier_upper": qmax,
        "nonaffine_forcing_coefficient": cnon,
        "coherent_deformation_action": tau / (3.0 * cnon),
        "phase_holonomy": float(PHASE_HOLONOMY_FLAT),
    }


@dataclass(frozen=True)
class FixedTransferLossGate:
    """The transfer-deficit channel of the physical whole-block gate.

    This certificate contains no causal mass.  It only says whether the supplied
    physical block deficit crosses the already-certified transfer threshold.
    A caller that owns an actual causal sublaw must bind *that same law* to this
    certificate before handing it to the single-charge compiler.
    """

    tau: float
    avg_transfer_deficit: float
    threshold: float
    triggered: bool
    cause: str | None

    def __post_init__(self) -> None:
        vals = (self.tau, self.avg_transfer_deficit, self.threshold)
        if not all(math.isfinite(float(v)) for v in vals):
            raise ValueError("finite fixed-transfer gate data required")
        if not (0.0 < self.tau <= 0.1) or self.avg_transfer_deficit < 0.0 or self.threshold <= 0.0:
            raise ValueError("invalid fixed-transfer gate data")
        expected = self.avg_transfer_deficit >= self.threshold
        if self.triggered != expected:
            raise AssertionError("fixed-transfer gate trigger disagrees with the physical deficit")
        expected_cause = "physical_transfer_cost" if expected else None
        if self.cause != expected_cause:
            raise AssertionError("fixed-transfer gate cause disagrees with the physical threshold crossing")


def fixed_transfer_loss_gate(*, tau: float, avg_transfer_deficit: float) -> FixedTransferLossGate:
    """Expose the exact transfer channel used by ``coherent_service_or_flat_gate``."""
    deficit = float(avg_transfer_deficit)
    if not math.isfinite(deficit) or deficit < 0.0:
        raise ValueError("finite nonnegative physical block transfer deficit required")
    th = coherent_flat_thresholds(float(tau))
    threshold = float(th["block_transfer_deficit"])
    triggered = deficit >= threshold
    return FixedTransferLossGate(
        tau=float(tau),
        avg_transfer_deficit=deficit,
        threshold=threshold,
        triggered=triggered,
        cause="physical_transfer_cost" if triggered else None,
    )


def coherent_connection_flatness_upper(
    avg_transfer_deficit: float,
    objective_variation_action: float,
    coherent_deformation_action: float,
    qmax: float,
) -> dict[str, float]:
    if min(avg_transfer_deficit, objective_variation_action, coherent_deformation_action) < 0:
        raise ValueError("nonnegative coherent flatness data required")
    h = hodge_flatness_upper(avg_transfer_deficit)
    s = affine_strain_flatness_upper(avg_transfer_deficit, objective_variation_action)
    n = coherent_nonaffine_coefficient(qmax) * coherent_deformation_action
    return {
        "hodge_rms": h,
        "nonconformal_strain_number": s,
        "coherent_nonaffine_connection_action": n,
        "coherent_kelvin_connection_flatness": h + s + n,
    }


def coherent_deformation_dissipation_threshold(
    tau: float,
    scaled_lifetime: float,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    carrier_ratio: float = DEFAULT_SHELL_CARRIER_RATIO,
    aspect_cap: float = EXTENDED_ASPECT,
) -> float:
    if scaled_lifetime <= 0:
        raise ValueError("positive scaled lifetime required")
    th = coherent_flat_thresholds(tau, radius_cap, carrier_ratio, aspect_cap)
    return normalized_dissipation_from_coherent_deformation(
        th["coherent_deformation_action"], scaled_lifetime, aspect_cap
    )


def coherent_service_or_flat_gate(
    *,
    tau: float,
    avg_transfer_deficit: float,
    objective_variation_action: float,
    total_strain_action: float,
    coherent_deformation_action: float,
    aspect: float,
    scale_radius: float,
    has_predecessor: bool,
    scaled_lifetime: float,
    phase_holonomy: float = 0.0,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    carrier_ratio: float = DEFAULT_SHELL_CARRIER_RATIO,
    aspect_cap: float = EXTENDED_ASPECT,
) -> dict[str, object]:
    """Whole-eddy service-or-flat gate without choosing an arbitrary primary cause.

    Several physical observables may cross their thresholds simultaneously.  The
    gate therefore returns *all* triggered causal roots.  The single-charge
    compiler, using first physical stopping/provenance, chooses the primary root.
    Only an empty trigger set is allowed to declare the block coherent-flat.
    """
    vals = (
        avg_transfer_deficit,
        objective_variation_action,
        total_strain_action,
        coherent_deformation_action,
        phase_holonomy,
    )
    if min(vals) < 0 or aspect < 1 or scale_radius <= 0 or scaled_lifetime <= 0:
        raise ValueError("invalid coherent block data")
    th = coherent_flat_thresholds(tau, radius_cap, carrier_ratio, aspect_cap)
    roots: list[dict[str, object]] = []
    transfer_gate = fixed_transfer_loss_gate(tau=tau, avg_transfer_deficit=avg_transfer_deficit)
    if not math.isclose(transfer_gate.threshold, th["block_transfer_deficit"], rel_tol=2e-15, abs_tol=0.0):
        raise AssertionError("whole-block gate and fixed-transfer channel thresholds diverged")
    if transfer_gate.triggered:
        roots.append({
            "cause": transfer_gate.cause,
            "threshold": transfer_gate.threshold,
            "value": transfer_gate.avg_transfer_deficit,
        })
    if objective_variation_action >= th["objective_strain_variation_action"]:
        roots.append({
            "cause": "coherent_averaged_strain_source",
            "threshold": th["objective_strain_variation_action"],
            "value": objective_variation_action,
        })
    if total_strain_action > th["low_strain_action"]:
        roots.append({
            "cause": "high_strain_critical_dissipation",
            "threshold": th["low_strain_action"],
            "value": total_strain_action,
            "normalized_dissipation_lower": high_strain_dissipation_lower(total_strain_action, scaled_lifetime),
        })
    if aspect > th["aspect_threshold"]:
        roots.append({
            "cause": "inherited_high_aspect" if has_predecessor else "fresh_high_aspect",
            "threshold": th["aspect_threshold"],
            "value": aspect,
        })
    if scale_radius > th["radius_cap"]:
        roots.append({
            "cause": "large_affine_radius_ancestry",
            "threshold": th["radius_cap"],
            "value": scale_radius,
            "critical_mass_lower": 0.3 * th["radius_cap"],
        })
    if phase_holonomy >= th["phase_holonomy"]:
        roots.append({
            "cause": "helical_phase_holonomy",
            "threshold": th["phase_holonomy"],
            "value": phase_holonomy,
        })
    if coherent_deformation_action >= th["coherent_deformation_action"]:
        roots.append({
            "cause": "coherent_deformation_critical_dissipation",
            "threshold": th["coherent_deformation_action"],
            "value": coherent_deformation_action,
            "normalized_dissipation_lower": normalized_dissipation_from_coherent_deformation(
                coherent_deformation_action, scaled_lifetime, aspect_cap
            ),
        })

    if roots:
        return {
            "status": "named_physical_causes",
            "triggered_causes": tuple(roots),
            "primary_selected": False,
            "primary_rule": "delegate to first-causal single-charge compiler",
        }

    flat = coherent_connection_flatness_upper(
        avg_transfer_deficit,
        objective_variation_action,
        coherent_deformation_action,
        th["intrinsic_carrier_upper"],
    )
    tol = 3e-13
    if flat["hodge_rms"] > tau / 3.0 + tol:
        raise AssertionError(("coherent Hodge threshold design failed", flat, th))
    if flat["nonconformal_strain_number"] > tau / 3.0 + tol:
        raise AssertionError(("coherent strain threshold design failed", flat, th))
    if flat["coherent_nonaffine_connection_action"] > tau / 3.0 + tol:
        raise AssertionError(("coherent nonaffine threshold design failed", flat, th))
    return {
        "status": "coherent_kelvin_extremal_flat",
        "triggered_causes": (),
        "tau": tau,
        "phase_flat": True,
        "profile_persistence_required": False,
        "eventwise_profile_marking": True,
        **flat,
    }


def theorem_certificate() -> dict[str, object]:
    th = coherent_flat_thresholds(0.01)
    dmin = coherent_deformation_dissipation_threshold(0.01, 1.0)
    return {
        "status": "EXACT_COHERENT_SERVICE_OR_FLAT_ASSEMBLY_GATE__EVENT_ROLE_REGISTRATION_SUPPLIED",
        "flatness": "sqrt(E_H_phys)+(dT)_nonconf+C(qmax) I_K <= tau",
        "nonaffine_coefficient": f"C(qmax)=1+qmax/sqrt2+sqrt7/2; default qmax={th['intrinsic_carrier_upper']:.12g}",
        "deformation_threshold": f"I_K<tau/(3C(qmax)); at tau=1/100 default threshold={th['coherent_deformation_action']:.12g}",
        "large_deformation": f"otherwise D_V>={dmin:.12g}/c for the threshold event (displayed at c=1)",
        "large_radius": "N r_g>s0 is sticky affine-radius ancestry with critical local mass >=(3/10)s0, not a finite reset",
        "profile_rule": "Christ Gaussian is an eventwise analysis mark; no frozen-profile persistence is an input to the coherent gate",
        "single_charge_rule": "the gate returns all threshold-crossing physical roots; its transfer channel is the reusable FixedTransferLossGate, and it never picks a primary by theorem-name order before the first-causal compiler",
        "H1H3_status": "existing H1/H3 theorems remain sharper optional diagnostics, not required to control the full non-affine moving Gaussian core",
        "continuum_status": "event-anchored hard roles and smooth PDE envelopes are supplied by companion theorems; remaining continuum work is the single recursive first-stop constructor which applies this gate and the physical pair-productivity route without duplicate charging",
    }


@dataclass(frozen=True)
class CoherentServiceFlatStress:
    samples: int
    minimum_hodge_flat_margin: float
    minimum_strain_flat_margin: float
    minimum_nonaffine_flat_margin: float
    minimum_dissipation_threshold: float
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260808) -> CoherentServiceFlatStress:
    rng = np.random.default_rng(seed)
    mh = ms = mn = md = float("inf")
    counts: dict[str, int] = {}
    tau = 0.01
    th = coherent_flat_thresholds(tau)
    for _ in range(samples):
        mode = int(rng.integers(0, 8))
        deficit = float(rng.uniform(0.0, 0.95)) * th["block_transfer_deficit"]
        obj = float(rng.uniform(0.0, 0.95)) * th["objective_strain_variation_action"]
        strain = float(rng.uniform(0.0, float(LOW_STRAIN_ACTION)))
        IK = float(rng.uniform(0.0, 0.95)) * th["coherent_deformation_action"]
        aspect = float(rng.uniform(1.0, EXTENDED_ASPECT))
        radius = float(rng.uniform(0.7, th["radius_cap"]))
        phase = float(rng.uniform(0.0, 0.95)) * th["phase_holonomy"]
        predecessor = bool(rng.integers(0, 2))
        c = float(rng.uniform(0.05, 2.0))
        if mode == 0:
            deficit = 1.02 * th["block_transfer_deficit"]
        elif mode == 1:
            obj = 1.02 * th["objective_strain_variation_action"]
        elif mode == 2:
            strain = 1.02 * th["low_strain_action"]
        elif mode == 3:
            aspect = 1.02 * th["aspect_threshold"]
        elif mode == 4:
            radius = 1.02 * th["radius_cap"]
        elif mode == 5:
            phase = 1.02 * th["phase_holonomy"]
        elif mode == 6:
            IK = 1.02 * th["coherent_deformation_action"]
        out = coherent_service_or_flat_gate(
            tau=tau,
            avg_transfer_deficit=deficit,
            objective_variation_action=obj,
            total_strain_action=strain,
            coherent_deformation_action=IK,
            aspect=aspect,
            scale_radius=radius,
            has_predecessor=predecessor,
            scaled_lifetime=c,
            phase_holonomy=phase,
        )
        status = str(out["status"])
        counts[status] = counts.get(status, 0) + 1
        if status == "coherent_kelvin_extremal_flat":
            mh = min(mh, tau / 3.0 - float(out["hodge_rms"]))
            ms = min(ms, tau / 3.0 - float(out["nonconformal_strain_number"]))
            mn = min(mn, tau / 3.0 - float(out["coherent_nonaffine_connection_action"]))
        else:
            causes = tuple(out["triggered_causes"])
            for root in causes:
                if root["cause"] == "coherent_deformation_critical_dissipation":
                    d = float(root["normalized_dissipation_lower"])
                    md = min(md, d)
                    if d <= 0:
                        raise AssertionError("large coherent deformation did not force positive critical dissipation")
    for x, name in ((mh, "Hodge"), (ms, "strain"), (mn, "nonaffine")):
        if not math.isfinite(x) or x < -5e-13:
            raise AssertionError(f"coherent flat {name} margin failed")
    if not math.isfinite(md) or md <= 0:
        raise AssertionError("coherent deformation branch was not exercised")
    return CoherentServiceFlatStress(samples, mh, ms, mn, md, counts)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-coherent-service-or-flat"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "coherent_service_or_flat.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    th = coherent_flat_thresholds(0.01)
    md = f"""# Coherent service-or-flat assembly gate\n\nStatus: **{cert['status']}**.\n\nThe master-facing assembly no longer needs a frozen Christ Gaussian to persist through a packet lifetime.  On each near-extremal physical-transfer event, use the certified Gaussian only as an analysis/shape mark; move the analysis eddy by the whole-eddy coherent affine regression.  The exact resolved-role theorems then identify three gauge-quotiented flatness observables:\n\n`F_coh = sqrt(E_H^phys) + (dT)_nonconf + C(qmax) I_K`.\n\nFor the default scale-matched branch `s=N r_g<=4`, transition aspect `kappa<=567/500`, and outer carrier `|k|/N<=exp(2/25)`,\n\n`qmax={th['intrinsic_carrier_upper']:.12g}`,\n`C(qmax)={th['nonaffine_forcing_coefficient']:.12g}`.\n\nAt `tau=1/100`, the conservative existing transfer gate makes the Hodge and nonconformal-strain pieces each `<tau/3`.  Choose\n\n`I_K < tau/[3 C(qmax)] = {th['coherent_deformation_action']:.12g}`.\n\nThen the **entire** spatial non-affine Gaussian-core connection action is `<tau/3`, so `F_coh<=tau`.  No H1/H3 profile-persistence hypothesis is needed for this assembly statement.\n\nIf `I_K` crosses that threshold, coherent deformation forces positive critical `D_V`; if strain is large, use the existing high-strain `D_V`; if averaged objective strain variation is large, use the coherent averaged pressure/SGS/viscous source calculus; if `N r_g>4`, route to sticky large-radius ancestry; high aspect, phase holonomy and direct transfer loss keep their existing destinations.\n\nThus the coherent gate is\n\n`transfer/source/high-strain/aspect/radius/phase/coherent-D_V  OR  coherent-Kelvin-flat`.\n\nThe original H1/H3 theorems are not discarded: they remain sharper transfer-facing diagnostics for low Hermite curvature.  They are simply no longer a logical persistence assumption needed to keep one frozen Gaussian packet alive.\n\nStress: `{out.samples}` branch states\n- minimum Hodge flat margin: `{out.minimum_hodge_flat_margin:.3e}`\n- minimum strain flat margin: `{out.minimum_strain_flat_margin:.3e}`\n- minimum nonaffine flat margin: `{out.minimum_nonaffine_flat_margin:.3e}`\n- minimum sampled coherent-deformation `D_V` lower: `{out.minimum_dissipation_threshold:.3e}`\n- branches: `{out.branch_counts}`\n\nThe remaining continuum bridge is now specifically **eventwise Gaussian marking and selected-role registration**: prove measurably that every efficient recursive physical-transfer event receives a near-Gaussian affine mark and that changing marks between common causal slices is only common affine gauge, existing covariance/symbol `Xi`, or a genuine material relink/strain/radius event.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
