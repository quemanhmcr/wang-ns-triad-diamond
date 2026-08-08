from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import numpy as np

from src.asynchronous_duhamel_sync import choose_heavy_half, initial_parent_span_ratio

LOW_STRAIN_ACTION = Fraction(1, 30)
INHERIT_ENERGY_FRACTION = Fraction(1, 5)
RESIDUAL_WORK_FRACTION = Fraction(1, 5)
PHYSICAL_HH_WORK_FRACTION = Fraction(8, 15)


def positive_source_work(c: np.ndarray, forcing: np.ndarray) -> float:
    """Positive instantaneous energy work 2[Re<c,F>]_+ for one child role."""
    c = np.asarray(c, complex)
    forcing = np.asarray(forcing, complex)
    if c.shape != forcing.shape or c.ndim != 1:
        raise ValueError("matching one-dimensional child/forcing vectors required")
    return 2.0 * max(0.0, float(np.real(np.vdot(c, forcing))))


def adjoint_projection_work_decomposition(
    c: np.ndarray,
    psi: np.ndarray,
    forcing: np.ndarray,
) -> dict[str, float]:
    """Exact physical-work split into the adjoint-response direction and its orthogonal child residual.

    Let z=<psi,c>, n=||psi||^2 and c_parallel=psi*z/n.  Then for every forcing atom

      2 Re<c,F>
      = 2/n Re(conj(z)<psi,F>) + 2 Re<c-c_parallel,F>.

    The first term is the correct energy lift of the adjoint response atom.  Raw
    dGamma omits the state-dependent factor and therefore is not itself an energy
    transfer measure.
    """
    c = np.asarray(c, complex)
    psi = np.asarray(psi, complex)
    forcing = np.asarray(forcing, complex)
    if not (c.shape == psi.shape == forcing.shape) or c.ndim != 1:
        raise ValueError("matching one-dimensional vectors required")
    n = float(np.real(np.vdot(psi, psi)))
    if n <= 0:
        raise ValueError("nonzero adjoint response required")
    z = complex(np.vdot(psi, c))
    source_atom = complex(np.vdot(psi, forcing))
    c_parallel = psi * (z / n)
    total = 2.0 * float(np.real(np.vdot(c, forcing)))
    response = 2.0 * float(np.real(np.conj(z) * source_atom)) / n
    orthogonal = 2.0 * float(np.real(np.vdot(c - c_parallel, forcing)))
    return {
        "physical_work": total,
        "response_energy_work": response,
        "orthogonal_child_work": orthogonal,
        "residual": total - response - orthogonal,
    }


def physical_hh_work_lower_bound(
    *,
    terminal_energy: float,
    initial_energy: float,
    residual_positive_work: float,
    strain_action: float,
) -> float:
    """Conservative lower bound from the selected-child energy inequality.

    For c_dot=Gc+F_HH+R and sym(G)<=||S|| I (viscosity is nonpositive),

      E_1 <= exp(2K)(E_0+W_HH^+ + W_R^+).

    Using exp(-2K)>=1-2K gives the algebraic lower bound below.  The theorem is
    used only for K<=1/30, where 1-2K is positive.
    """
    vals = (terminal_energy, initial_energy, residual_positive_work, strain_action)
    if any((not math.isfinite(v) or v < 0) for v in vals):
        raise ValueError("finite nonnegative energy/work/action data required")
    if terminal_energy <= 0:
        raise ValueError("positive terminal child energy required")
    if strain_action > 0.5:
        raise ValueError("linear exponential lower bound requires K<=1/2")
    return (1.0 - 2.0 * strain_action) * terminal_energy - initial_energy - residual_positive_work


def route_physical_energy_causality(
    *,
    terminal_energy: float,
    initial_energy: float,
    residual_positive_work: float,
    strain_action: float,
) -> dict[str, float | str]:
    """Low-strain physical inherit/residual/generate gate using actual child-energy work.

    At K>1/30 the existing high-strain theorem owns the block.  Otherwise an
    initial-energy fraction >=1/5 is sticky material inheritance; classified
    residual positive work >=1/5 E_1 delegates to its existing source/interface
    destination.  If neither fires, actual positive designated high-high child
    work is at least 8/15 E_1.
    """
    if terminal_energy <= 0:
        raise ValueError("positive terminal child energy required")
    if min(initial_energy, residual_positive_work, strain_action) < 0:
        raise ValueError("nonnegative data required")
    if strain_action > float(LOW_STRAIN_ACTION):
        return {
            "branch": "high_strain_critical_dissipation",
            "value": float(strain_action),
            "threshold": float(LOW_STRAIN_ACTION),
        }
    if initial_energy >= float(INHERIT_ENERGY_FRACTION) * terminal_energy:
        return {
            "branch": "material_energy_inheritance",
            "value": float(initial_energy),
            "threshold": float(INHERIT_ENERGY_FRACTION) * terminal_energy,
        }
    if residual_positive_work >= float(RESIDUAL_WORK_FRACTION) * terminal_energy:
        return {
            "branch": "classified_residual_physical_work",
            "value": float(residual_positive_work),
            "threshold": float(RESIDUAL_WORK_FRACTION) * terminal_energy,
        }
    lower = physical_hh_work_lower_bound(
        terminal_energy=terminal_energy,
        initial_energy=initial_energy,
        residual_positive_work=residual_positive_work,
        strain_action=strain_action,
    )
    clean = float(PHYSICAL_HH_WORK_FRACTION) * terminal_energy
    tol = 3e-15 * max(1.0, terminal_energy)
    if lower + tol < clean:
        raise AssertionError("low-strain energy gate failed to force the clean 8/15 high-high work fraction")
    return {
        "branch": "physical_high_high_transfer_generation",
        "physical_hh_work_lower": float(lower),
        "clean_threshold": clean,
    }


def heavy_half_physical_transfer(
    times: Sequence[float],
    positive_work_weights: Sequence[float],
    slab_start: float,
    slab_end: float,
) -> dict[str, float | int]:
    """Synchronize using actual positive HH child-transfer weights, not dGamma weights."""
    out = choose_heavy_half(times, positive_work_weights, slab_start, slab_end)
    total = float(out["total"])
    if total <= 0:
        raise ValueError("positive physical high-high transfer mass required")
    if float(out["mass"]) + 2e-14 * max(1.0, total) < 0.5 * total:
        raise AssertionError("physical-transfer half-slab pigeonhole failed")
    out = dict(out)
    out["normalized_parent_span_upper"] = float(initial_parent_span_ratio())
    return out


def flat_scalar_measure_countermodel(t: float = 0.5) -> dict[str, float]:
    """Exact flat scalar obstruction to K_Gamma=K_phys.

    For c'=1, c(0)=0, G=R=0 on [0,1], terminal adjoint psi=1.  Normalized
    amplitude generation is dGamma=dt, while normalized physical child-energy
    transfer is dT=2t dt.  Hence Gamma([0,t])=t and T([0,t])=t^2.
    """
    if not 0.0 <= t <= 1.0:
        raise ValueError("t must lie in [0,1]")
    gamma_cdf = float(t)
    physical_cdf = float(t * t)
    return {
        "time": float(t),
        "duhamel_cdf": gamma_cdf,
        "physical_transfer_cdf": physical_cdf,
        "cdf_gap": gamma_cdf - physical_cdf,
        "energy_lift_density": 2.0 * float(t),
    }


@dataclass(frozen=True)
class PhysicalEnergyCausalStress:
    samples: int
    worst_projection_identity_residual: float
    minimum_generation_margin: float
    minimum_half_mass_margin: float
    countermodel_half_cdf_gap: float
    branch_counts: dict[str, int]


def stress(samples: int = 50_000, seed: int = 20260808) -> PhysicalEnergyCausalStress:
    rng = np.random.default_rng(seed)
    worst_proj = 0.0
    min_gen = float("inf")
    min_half = float("inf")
    counts: dict[str, int] = {}

    for _ in range(samples):
        n = int(rng.integers(1, 7))
        c = rng.normal(size=n) + 1j * rng.normal(size=n)
        psi = rng.normal(size=n) + 1j * rng.normal(size=n)
        F = rng.normal(size=n) + 1j * rng.normal(size=n)
        out = adjoint_projection_work_decomposition(c, psi, F)
        scale = max(1.0, abs(float(out["physical_work"])), abs(float(out["response_energy_work"])), abs(float(out["orthogonal_child_work"])))
        resid = abs(float(out["residual"]))
        worst_proj = max(worst_proj, resid)
        if resid > 2e-12 * scale:
            raise AssertionError("adjoint response/full physical-work decomposition failed")

        E1 = float(rng.lognormal(mean=0.0, sigma=1.0))
        K = float(rng.uniform(0.0, 0.05))
        mode = int(rng.integers(0, 4))
        E0 = float(rng.uniform(0.0, 0.19)) * E1
        WR = float(rng.uniform(0.0, 0.19)) * E1
        if mode == 0:
            K = float(rng.uniform(float(LOW_STRAIN_ACTION) + 1e-5, 0.05))
        elif mode == 1:
            K = float(rng.uniform(0.0, float(LOW_STRAIN_ACTION)))
            E0 = float(rng.uniform(0.2, 0.8)) * E1
        elif mode == 2:
            K = float(rng.uniform(0.0, float(LOW_STRAIN_ACTION)))
            WR = float(rng.uniform(0.2, 0.8)) * E1
        else:
            K = float(rng.uniform(0.0, float(LOW_STRAIN_ACTION)))
        route = route_physical_energy_causality(
            terminal_energy=E1,
            initial_energy=E0,
            residual_positive_work=WR,
            strain_action=K,
        )
        b = str(route["branch"])
        counts[b] = counts.get(b, 0) + 1
        if b == "physical_high_high_transfer_generation":
            margin = float(route["physical_hh_work_lower"]) - float(route["clean_threshold"])
            min_gen = min(min_gen, margin)
            if margin < -2e-13 * max(1.0, E1):
                raise AssertionError("physical high-high work lower bound lost the clean threshold")

        m = int(rng.integers(2, 80))
        t0 = float(rng.uniform(-2.0, 2.0))
        T = float(rng.lognormal(mean=-1.0, sigma=1.0))
        times = t0 + rng.random(m) * T
        weights = rng.lognormal(mean=-1.0, sigma=1.0, size=m)
        hh = heavy_half_physical_transfer(times, weights, t0, t0 + T)
        hmargin = float(hh["mass"]) - 0.5 * float(hh["total"])
        min_half = min(min_half, hmargin)
        if hmargin < -2e-12 * max(1.0, float(hh["total"])):
            raise AssertionError("physical transfer failed half-slab synchronization")

    if not math.isfinite(min_gen):
        min_gen = 0.0
    cm = flat_scalar_measure_countermodel(0.5)
    if not math.isclose(float(cm["duhamel_cdf"]), 0.5) or not math.isclose(float(cm["physical_transfer_cdf"]), 0.25):
        raise AssertionError("flat scalar measure countermodel changed")
    return PhysicalEnergyCausalStress(
        samples=samples,
        worst_projection_identity_residual=worst_proj,
        minimum_generation_margin=min_gen,
        minimum_half_mass_margin=min_half,
        countermodel_half_cdf_gap=float(cm["cdf_gap"]),
        branch_counts=counts,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_PHYSICAL_ENERGY_CAUSAL_GATE_AND_MEASURE_AGNOSTIC_SYNC__CONTINUUM_EXTRACTION_CONDITIONAL",
        "countermodel": "even for G=R=0, c'=1, c(0)=0: normalized dGamma=dt but normalized physical dT=2t dt, so raw kernels are not equal",
        "adjoint_energy_lift": "2 Re<c,F> = 2||psi||^-2 Re(conj(<psi,c>)<psi,F>) + 2 Re<c_perp,F>",
        "energy_gate": "E1 <= exp(2K)(E0+W_HH^+ + W_R^+)",
        "clean_low_strain_route": "K<=1/30, E0<E1/5, W_R^+<E1/5 => W_HH^+ >= 8E1/15",
        "physical_causal_measure": "dT_HH=2[Re<c,F_HH,alpha>]_+ dt on same-time parent-pair events",
        "synchronization": "the half-slab and parabolic cone use only positivity/support, so they apply to dT_HH directly",
        "renyi_rule": "Shannon/Renyi parent-slot weights are taken from normalized dT_HH, not from raw dGamma",
        "duhamel_role": "adjoint Duhamel remains an amplitude/support diagnostic; it is not required to equal the physical transfer kernel",
        "continuum_status": "still requires actual recursive Navier-Stokes block extraction with the exact selected coefficient equation and canonical material labels",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-physical-energy-causal-bridge"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    payload = {"certificate": cert, "stress": asdict(out), "flat_scalar_countermodel": flat_scalar_measure_countermodel(0.5)}
    (args.outdir / "physical_energy_causal_bridge.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Physical-energy causal bridge\n\nStatus: **{cert['status']}**.\n\nRaw adjoint amplitude generation and physical child-energy transfer are not the same measure: the exact flat scalar model `c'=1`, `c(0)=0`, `G=R=0` has normalized `dGamma=dt` but normalized `dT=2t dt`; at half-time their cumulative masses are `1/2` and `1/4`.\n\nThe correct exact algebra is the adjoint-response projection identity\n\n`2 Re<c,F> = 2/||psi||^2 Re(conj(<psi,c>) <psi,F>) + 2 Re<c_perp,F>`.\n\nMore importantly, the selected-child energy balance itself supplies the physical causal law.  On `K=int||S||<=1/30`,\n\n`E1 <= exp(2K)(E0+W_HH^+ + W_R^+)`.\n\nUsing `exp(-2K)>=1-2K>=14/15`, if inherited energy and classified residual positive work are each below `E1/5`, then\n\n`W_HH^+ >= 8 E1/15`.\n\nThus the generated branch carries a definite amount of **actual positive high-high child-energy work**.  Decompose `F_HH` into same-time parent-pair atoms and use\n\n`dT_HH=2[Re<c,F_HH,alpha>]_+ dt`\n\nas the causal probability law.  The asynchronous half-slab/parabolic synchronization argument depends only on positivity and event support, so it applies to this physical transfer law directly.  Shannon/Renyi therefore receive the same physical weights used by Hodge/flat transfer bookkeeping, without identifying raw `dGamma` with `dT`.\n\nThe adjoint Duhamel theorem remains useful to quotient common Kelvin transport and diagnose amplitude inheritance/residuals, but raw `dGamma` is no longer a required master-facing transfer kernel.\n\nStress: `{out.samples}` synthetic algebra/synchronization states\n- worst adjoint-projection work residual: `{out.worst_projection_identity_residual:.3e}`\n- minimum physical-HH generation margin: `{out.minimum_generation_margin:.3e}`\n- minimum physical-transfer half-slab margin: `{out.minimum_half_mass_margin:.3e}`\n- exact flat countermodel CDF gap at half-time: `{out.countermodel_half_cdf_gap:.6f}`\n- branches: `{out.branch_counts}`\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
