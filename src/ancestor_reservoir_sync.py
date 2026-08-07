from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.linalg import expm

from src.single_edge_certificate import RSTAR_HI

SIGNED_GOOD_GAP = Fraction(1, 80)
CLEAN_GENERATION_PROGRESS = Fraction(8, 5)
LOW_STRAIN_ACTION = Fraction(1, 30)
CLEAN_RESERVOIR_GROWTH = Fraction(21, 20)


def generation_progress_lower_float() -> float:
    """Signed-good child/top-parent scale progress exp(gamma*-1/80)."""
    return math.exp(-float(SIGNED_GOOD_GAP)) / float(RSTAR_HI)


def low_strain_frequency_growth_upper_float() -> float:
    """Kelvin support growth on int ||S|| dt <=1/30."""
    return math.exp(float(LOW_STRAIN_ACTION))


def critical_mass_service_ratio_upper() -> Fraction:
    """Per-unit-critical-mass service coefficient ratio between generations.

    A low band at frequency M contributes to the LP increment square with
      (M/N)^2 mu.
    Under M1/M0<=21/20 and N1/N0>=8/5 the coefficient contracts by this factor.
    """
    return (CLEAN_RESERVOIR_GROWTH / CLEAN_GENERATION_PROGRESS) ** 2


def physical_energy_service_ratio_upper() -> Fraction:
    """Per-unit-physical-energy service coefficient ratio.

    Since mu=M E, low-band increment-square service is proportional to M^3 E/N^2.
    """
    return CLEAN_RESERVOIR_GROWTH**3 / CLEAN_GENERATION_PROGRESS**2


def low_band_service_from_critical_mass(
    reservoir_frequency: float,
    filter_frequency: float,
    critical_mass: float,
    beta_filter_radius: float = 1.0,
) -> float:
    if min(reservoir_frequency, filter_frequency, beta_filter_radius) <= 0 or critical_mass < 0:
        raise ValueError("invalid service data")
    return (beta_filter_radius * reservoir_frequency / filter_frequency) ** 2 * critical_mass


def low_band_service_from_physical_energy(
    reservoir_frequency: float,
    filter_frequency: float,
    physical_energy: float,
    beta_filter_radius: float = 1.0,
) -> float:
    if min(reservoir_frequency, filter_frequency, beta_filter_radius) <= 0 or physical_energy < 0:
        raise ValueError("invalid service data")
    return beta_filter_radius**2 * reservoir_frequency**3 * physical_energy / filter_frequency**2


def amortized_service_capacity_upper(
    generation: int,
    reservoir_frequency0: float,
    filter_frequency0: float,
    global_energy_cap: float,
    beta_filter_radius: float = 1.0,
) -> float:
    """Maximum low-band service of one materially reused low-strain reservoir.

    We allow the reservoir to own the entire global L2 energy at every service time,
    so this bound does not assume monotonicity of the reservoir's own amplitude.
    """
    if generation < 0:
        raise ValueError("generation must be nonnegative")
    base = low_band_service_from_physical_energy(
        reservoir_frequency0, filter_frequency0, global_energy_cap, beta_filter_radius
    )
    ratio = float(physical_energy_service_ratio_upper())
    return base * ratio**generation


def total_amortized_service_upper(
    reservoir_frequency0: float,
    filter_frequency0: float,
    global_energy_cap: float,
    beta_filter_radius: float = 1.0,
) -> float:
    base = low_band_service_from_physical_energy(
        reservoir_frequency0, filter_frequency0, global_energy_cap, beta_filter_radius
    )
    r = float(physical_energy_service_ratio_upper())
    return base / (1.0 - r)


def max_uniform_service_generations(
    service_threshold: float,
    reservoir_frequency0: float,
    filter_frequency0: float,
    global_energy_cap: float,
    beta_filter_radius: float = 1.0,
) -> int:
    """Largest q for which the universal capacity upper can still exceed threshold."""
    if service_threshold <= 0:
        raise ValueError("positive service threshold required")
    base = low_band_service_from_physical_energy(
        reservoir_frequency0, filter_frequency0, global_energy_cap, beta_filter_radius
    )
    if base < service_threshold:
        return -1
    r = float(physical_energy_service_ratio_upper())
    return int(math.floor(math.log(service_threshold / base) / math.log(r)))


def kelvin_covector(L: np.ndarray, k: np.ndarray) -> np.ndarray:
    L = np.asarray(L, float)
    k = np.asarray(k, float)
    return L.T @ k


def affine_material_step(A: np.ndarray, dt: float, L: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Exact constant-A material/Kelvin step: L1=e^{A dt}L0, k1=e^{-A^T dt}k0."""
    A = np.asarray(A, float)
    L = np.asarray(L, float)
    k = np.asarray(k, float)
    return expm(A * dt) @ L, expm(-A.T * dt) @ k


def reservoir_frequency_growth_bound(strain_action: float) -> float:
    if strain_action < 0:
        raise ValueError("nonnegative strain action required")
    return math.exp(strain_action)


def arb_reservoir_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160
    rhi = arb(RSTAR_HI.numerator) / RSTAR_HI.denominator
    progress = (-arb(1) / 80).exp() / rhi
    growth = (arb(1) / 30).exp()
    if not (progress > arb(8) / 5):
        raise AssertionError(f"signed-good scale progress did not exceed 8/5: {progress}")
    if not (growth < arb(21) / 20):
        raise AssertionError(f"low-strain Kelvin growth did not stay below 21/20: {growth}")

    cm = critical_mass_service_ratio_upper()
    pe = physical_energy_service_ratio_upper()
    if not (cm < Fraction(1, 2)):
        raise AssertionError("critical-mass service coefficient did not halve")
    if not (pe < Fraction(1, 2)):
        raise AssertionError("physical-energy service coefficient did not halve")
    if not (Fraction(1, 1) / (Fraction(1, 1) - pe) < 2):
        raise AssertionError("geometric service budget was not below two base capacities")
    return {
        "signed_good_progress_ball": str(progress),
        "clean_generation_progress": ">8/5",
        "low_strain_growth_ball": str(growth),
        "clean_reservoir_growth": "<21/20",
        "critical_mass_service_ratio": f"{cm.numerator}/{cm.denominator}",
        "critical_mass_service_clean": "<1/2",
        "physical_energy_service_ratio": f"{pe.numerator}/{pe.denominator}",
        "physical_energy_service_clean": "<1/2",
        "total_one_reservoir_budget": "<2 times its generation-0 energy-capacity coefficient",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class ReservoirSyncStress:
    samples: int
    worst_kelvin_covector_residual: float
    minimum_frequency_growth_margin: float
    maximum_critical_service_ratio: float
    maximum_energy_service_ratio: float
    minimum_amortized_margin: float


def stress(samples: int = 50_000, seed: int = 20260808) -> ReservoirSyncStress:
    rng = np.random.default_rng(seed)
    wq = 0.0
    mg = float("inf")
    mc = 0.0
    me = 0.0
    ma = float("inf")
    clean_progress = float(CLEAN_GENERATION_PROGRESS)
    clean_growth = float(CLEAN_RESERVOIR_GROWTH)
    cm_clean = float(critical_mass_service_ratio_upper())
    pe_clean = float(physical_energy_service_ratio_upper())

    for _ in range(samples):
        A = rng.normal(size=(3, 3))
        A -= np.trace(A) / 3.0 * np.eye(3)
        S = 0.5 * (A + A.T)
        sn = float(np.linalg.norm(S, 2))
        dt = float(rng.uniform(0.0, 1.0))
        if sn * dt > float(LOW_STRAIN_ACTION):
            dt *= float(LOW_STRAIN_ACTION) / max(sn * dt, 1e-30)
        L0 = rng.normal(size=(3, 3))
        while abs(np.linalg.det(L0)) < 0.1:
            L0 = rng.normal(size=(3, 3))
        k0 = rng.normal(size=3)
        if np.linalg.norm(k0) < 1e-12:
            k0[0] = 1.0
        L1, k1 = affine_material_step(A, dt, L0, k0)
        q0 = kelvin_covector(L0, k0)
        q1 = kelvin_covector(L1, k1)
        wq = max(wq, float(np.linalg.norm(q1 - q0)) / max(1.0, float(np.linalg.norm(q0))))
        actual_growth = float(np.linalg.norm(k1) / np.linalg.norm(k0))
        theoretical = math.exp(sn * dt)
        mg = min(mg, theoretical - actual_growth)
        if actual_growth > theoretical * (1.0 + 2e-12):
            raise AssertionError("Kelvin frequency growth exceeded symmetric-strain action")

        # Adversarial step ratios inside the certified clean bounds.
        mratio = float(rng.uniform(0.2, clean_growth * (1.0 - 1e-10)))
        nratio = float(rng.uniform(clean_progress * (1.0 + 1e-10), 3.0))
        cmratio = (mratio / nratio) ** 2
        eratio = mratio**3 / nratio**2
        mc = max(mc, cmratio)
        me = max(me, eratio)
        if cmratio >= 0.5 + 1e-12 or eratio >= 0.5 + 1e-12:
            raise AssertionError("one-generation service half-life failed")

        m0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        n0 = float(math.exp(rng.uniform(-2.0, 2.0)))
        E = float(math.exp(rng.uniform(-3.0, 2.0)))
        beta = float(rng.uniform(0.2, 2.0))
        q = int(rng.integers(0, 20))
        actual_upper = amortized_service_capacity_upper(q, m0, n0, E, beta)
        half_upper = low_band_service_from_physical_energy(m0, n0, E, beta) * (0.5**q)
        ma = min(ma, half_upper - actual_upper)
        if actual_upper > half_upper + 2e-12 * max(1.0, half_upper):
            raise AssertionError("clean half-life envelope failed")

    if wq > 5e-11:
        raise AssertionError("Kelvin grain covector invariance lost")
    return ReservoirSyncStress(samples, wq, mg, mc, me, ma)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-ancestor-reservoir-sync"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_reservoir_certificate()
    out = stress(args.samples)
    data = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "ancestor_reservoir_sync.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    md = f"""# Ancestor reservoir synchronization: Kelvin/Liouville spectral half-life

Status: **{cert['status']}**.

On the signed-good physical transfer core, consecutive lineage scales satisfy

`N_(q+1)/N_q > exp(gamma_*-1/80) > 8/5`.

For one materially reused low-frequency reservoir transported by the same affine/Kelvin background,

`d log|xi|/dt <= ||S||_op`.

Thus on the existing low-strain branch `int||S||dt<=1/30`, every carrier/support frequency grows by less than `exp(1/30)<21/20`.

The LP increment observable is Galilean/sweeping neutral: a low band at `M` contributes square-function service proportional to `(M/N)^2 mu_M`. Therefore its **per-unit-critical-mass** service coefficient contracts between consecutive generations by at most

`(21/32)^2 = 441/1024 < 1/2`.

Since `mu_M=M E_M`, the **per-unit-physical-energy** coefficient is `M^3/N^2` and contracts by

`(21/20)^3(5/8)^2 = 231525/512000 < 1/2`.

Even allowing this one reservoir to own the entire conserved global energy at every service time, its generation-q service capacity is bounded by a geometric half-life and its total service over all future generations is less than twice its generation-0 energy-capacity coefficient.

The material identity is not a label convention: for `Ldot=A L`, `kdot=-A^T k`, the grain covector `q=L^T k` is exactly invariant. Nonlinear creation of a new band/carrier is therefore a relinking/source event, not free reuse.

Stress: `{out.samples}`
- worst Kelvin covector residual: `{out.worst_kelvin_covector_residual:.3e}`
- minimum Kelvin frequency-growth margin: `{out.minimum_frequency_growth_margin:.3e}`
- maximum critical-mass service ratio: `{out.maximum_critical_service_ratio:.9f}`
- maximum physical-energy service ratio: `{out.maximum_energy_service_ratio:.9f}`
- minimum clean half-life envelope margin: `{out.minimum_amortized_margin:.3e}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
