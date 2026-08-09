from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.coherent_increment_service import periodic_increment_covariance_residual
from src.affine_coherent_moyal import discrete_moyal_residual
from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower
from src.high_strain_resolved_ancestor import high_strain_ancestor_mass_threshold

RESOLVED_RADIUS_RATIO = 1.0 / 4.0
HEAT_TIME_FACTOR = 1.0 / 2.0  # heat time theta_N = 1/(2 N^2)
HEAT_SUPPORT_XMAX = 1.0 / 32.0
HEAT_DEFECT_LOWER = math.exp(-HEAT_SUPPORT_XMAX)


def heat_time(child_frequency: float) -> float:
    N = float(child_frequency)
    if N <= 0 or not math.isfinite(N):
        raise ValueError("positive finite child frequency required")
    return HEAT_TIME_FACTOR / (N * N)


def heat_increment_probability_density(radius: float, child_frequency: float) -> float:
    """3D heat kernel at theta_N=1/(2N^2).

    Its characteristic function is exp(-|xi|^2/(2N^2)).  The kernel is not a
    new causal law; it is a positive probe at the intrinsic parabolic diffusion
    scale of the selected NS block.
    """
    r = float(radius)
    N = float(child_frequency)
    if r < 0 or N <= 0 or not math.isfinite(r + N):
        raise ValueError("nonnegative radius and positive finite frequency required")
    return (N / math.sqrt(2.0 * math.pi)) ** 3 * math.exp(-0.5 * (N * r) ** 2)


def heat_increment_multiplier(frequency_radius: float, child_frequency: float) -> float:
    """N^2 times the exact mean-square increment multiplier.

    If H_N is the heat kernel above,
      N^2 int H_N(r)||delta_r f||_2^2 dr
    has Fourier multiplier 2N^2(1-exp(-|xi|^2/(2N^2))).
    """
    k = float(frequency_radius)
    N = float(child_frequency)
    if k < 0 or N <= 0 or not math.isfinite(k + N):
        raise ValueError("nonnegative finite frequency and positive child frequency required")
    x = 0.5 * (k / N) ** 2
    return 2.0 * N * N * (-math.expm1(-x))


def heat_to_gradient_ratio(frequency_radius: float, child_frequency: float) -> float:
    k = float(frequency_radius)
    if k == 0:
        return 1.0
    return heat_increment_multiplier(k, child_frequency) / (k * k)


def resolved_heat_gradient_bounds(child_frequency: float) -> dict[str, float]:
    """Sharp elementary support comparison for |xi|<=N/4.

    For x=|xi|^2/(2N^2)<=1/32,
      e^(-1/32) <= (1-e^-x)/x <= 1.
    """
    N = float(child_frequency)
    if N <= 0 or not math.isfinite(N):
        raise ValueError("positive finite child frequency required")
    edge = heat_to_gradient_ratio(RESOLVED_RADIUS_RATIO * N, N)
    if edge + 1e-15 < HEAT_DEFECT_LOWER or edge > 1.0 + 1e-15:
        raise AssertionError("resolved heat/gradient support comparison failed")
    return {
        "lower": HEAT_DEFECT_LOWER,
        "upper": 1.0,
        "edge_ratio": edge,
        "support_xmax": HEAT_SUPPORT_XMAX,
    }


def spectral_heat_service(
    child_frequency: float,
    frequency_radii: np.ndarray,
    resolved_spectral_energies: np.ndarray,
) -> dict[str, float]:
    """Exact finite spectral heat-defect/gradient comparison.

    `resolved_spectral_energies` are nonnegative spectral energy masses.  This
    is a quadrature model of Plancherel, not a packet synthesis.
    """
    N = float(child_frequency)
    k = np.asarray(frequency_radii, float)
    e = np.asarray(resolved_spectral_energies, float)
    if k.shape != e.shape or k.ndim != 1 or len(k) == 0:
        raise ValueError("matching nonempty one-dimensional spectral arrays required")
    if np.any(~np.isfinite(k)) or np.any(~np.isfinite(e)) or np.any(k < 0) or np.any(e < 0):
        raise ValueError("finite nonnegative spectral data required")
    if N <= 0 or np.any(k > RESOLVED_RADIUS_RATIO * N * (1.0 + 1e-13)):
        raise ValueError("spectral support must lie in B_(N/4)")
    grad = float(np.dot(k * k, e))
    mult = np.array([heat_increment_multiplier(x, N) for x in k], float)
    heat = float(np.dot(mult, e))
    return {
        "gradient_energy": grad,
        "heat_increment_energy": heat,
        "lower_margin": heat - HEAT_DEFECT_LOWER * grad,
        "upper_margin": grad - heat,
    }


def integrated_heat_service_lower_from_dissipation(
    normalized_dissipation: float,
) -> float:
    """For S=N^3 int dt int H_N(r)||delta_r V||^2 dr, S>=e^-1/32 D_V."""
    D = float(normalized_dissipation)
    if D < 0 or not math.isfinite(D):
        raise ValueError("finite nonnegative normalized dissipation required")
    return HEAT_DEFECT_LOWER * D


def high_strain_heat_service_lower(scaled_lifetime: float) -> float:
    return integrated_heat_service_lower_from_dissipation(clean_high_strain_dissipation_lower(float(scaled_lifetime)))


def critical_ancestor_heat_service_fraction_lower() -> float:
    """Fraction of total heat service retained on high-strain critical shell atoms.

    The shell theorem gives D(G)>=D/2.  Pointwise heat service satisfies
    S(G)>=e^-1/32 D(G), while S(total)<=D(total).  Therefore
      S(G)/S(total) >= e^-1/32 / 2.
    """
    return 0.5 * HEAT_DEFECT_LOWER


def retained_heat_service_bounds(
    total_dissipation: float,
    good_ancestor_dissipation: float,
) -> dict[str, float]:
    D = float(total_dissipation)
    G = float(good_ancestor_dissipation)
    if D < 0 or G < 0 or G > D or not math.isfinite(D + G):
        raise ValueError("0<=good<=total finite dissipation required")
    good_heat_lower = HEAT_DEFECT_LOWER * G
    total_heat_upper = D
    return {
        "good_heat_lower": good_heat_lower,
        "total_heat_upper": total_heat_upper,
        "fraction_lower": good_heat_lower / total_heat_upper if D > 0 else 0.0,
    }


def theorem_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    Dstar = clean_high_strain_dissipation_lower(c)
    mustar = high_strain_ancestor_mass_threshold(c)
    Sstar = high_strain_heat_service_lower(c)
    frac = critical_ancestor_heat_service_fraction_lower()
    return {
        "status": "EXACT_HIGH_STRAIN_TO_HEAT_INCREMENT_COHERENT_SERVICE__CRITICAL_RESOLVED_ANCESTOR_FRACTION_RETAINED__OLD_POOL_ROUTING_REMAINS",
        "heat_probe": "H_N is the NS heat kernel at theta_N=1/(2N^2), with Fourier characteristic exp(-|xi|^2/(2N^2))",
        "exact_increment_identity": "int H_N(r)||delta_r V||_2^2 dr = int 2(1-exp(-|xi|^2/(2N^2)))|Vhat(xi)|^2 dxi",
        "support_comparison": f"on supp Vhat subset B_(N/4), exp(-1/32)||grad V||_2^2 <= N^2 E_H||delta_r V||_2^2 <= ||grad V||_2^2; lower={HEAT_DEFECT_LOWER:.12g}",
        "integrated_service": f"S_heat=N^3 int dt E_H||delta_r V||_2^2 >= exp(-1/32)D_V; high strain gives S_heat>={Sstar:.12g}",
        "critical_ancestor_fraction": f"combining D_V(good)>=D_V/2 with pointwise heat/dissipation comparison puts at least exp(-1/32)/2={frac:.12g} of total heat-increment service on shell-time atoms with mu_j>=mu_*={mustar:.12g}",
        "moyal": "for every r and normalized coherent window, Moyal disintegrates ||delta_r P_j V||_2^2 into positive phase-space cell energy exactly",
        "translation": "V_g(delta_r f)(X,k)=exp(-ik.r)V_g f(X-r,k)-V_g f(X,k), so each cell atom is an actual coherent spatial edge separated by the physical heat displacement r",
        "causal_role": "heat increments are a positive diagnostic/service law dominated by and comparable to resolved dissipation; they do not replace HH transfer law or create a reset currency",
        "scope": "this supplies spatial/coherent disintegration for the high-strain ancestor route; old/new material pool routing and universal slab renewal remain",
    }


@dataclass(frozen=True)
class HeatIncrementStress:
    samples: int
    minimum_spectral_lower_margin: float
    minimum_spectral_upper_margin: float
    minimum_retained_fraction_margin: float
    worst_discrete_moyal_relative_residual: float
    worst_increment_covariance_residual: float


def stress(samples: int = 50_000, seed: int = 20260809) -> HeatIncrementStress:
    rng = np.random.default_rng(seed)
    ml = mu = mf = float("inf")
    wm = wc = 0.0
    target_frac = critical_ancestor_heat_service_fraction_lower()

    for i in range(samples):
        N = float(math.exp(rng.uniform(-2.0, 8.0)))
        n = int(rng.integers(2, 100))
        k = rng.random(n) * (N / 4.0)
        e = rng.lognormal(mean=0.0, sigma=2.0, size=n)
        out = spectral_heat_service(N, k, e)
        scale = max(1.0, float(out["gradient_energy"]))
        ml = min(ml, float(out["lower_margin"]))
        mu = min(mu, float(out["upper_margin"]))
        if float(out["lower_margin"]) < -5e-12 * scale or float(out["upper_margin"]) < -5e-12 * scale:
            raise AssertionError("heat increment left resolved gradient comparison")

        D = float(math.exp(rng.uniform(-3.0, 5.0)))
        G = float(rng.uniform(0.5, 1.0)) * D
        rb = retained_heat_service_bounds(D, G)
        margin = float(rb["fraction_lower"]) - target_frac
        mf = min(mf, margin)
        if margin < -2e-14:
            raise AssertionError("critical ancestor heat-service fraction fell below e^-1/32/2")

        # Representative exact coherent identities only; the 50k loop above
        # already stresses the analytic spectral theorem cheaply.
        if i < min(samples, 2500):
            nd = int(rng.integers(8, 36))
            f = rng.normal(size=nd) + 1j * rng.normal(size=nd)
            g = rng.normal(size=nd) + 1j * rng.normal(size=nd)
            g /= np.linalg.norm(g)
            mres = discrete_moyal_residual(f, g)
            wm = max(wm, abs(mres) / max(1.0, float(np.vdot(f, f).real)))
            if abs(mres) > 2e-10 * max(1.0, float(np.vdot(f, f).real)):
                raise AssertionError("Moyal coherent service partition failed")
            shift = int(rng.integers(-nd // 3, nd // 3 + 1))
            cres = periodic_increment_covariance_residual(f, g, shift)
            wc = max(wc, cres / max(1.0, float(np.linalg.norm(f))))
            if cres > 2e-10 * max(1.0, float(np.linalg.norm(f))):
                raise AssertionError("coherent heat-increment edge covariance failed")

    return HeatIncrementStress(samples, ml, mu, mf, wm, wc)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-high-strain-heat-increment-service"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "high_strain_heat_increment_service.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# High strain becomes heat-increment coherent service\n\nStatus: **{cert['status']}**.\n\nUse no arbitrary spatial packet.  Let `H_N` be the three-dimensional Navier--Stokes heat kernel at the block's intrinsic parabolic time `theta_N=1/(2N^2)`.  Equivalently `r` is centered Gaussian with covariance `N^-2 I`, so\n\n`E exp(-i xi.r)=exp(-|xi|^2/(2N^2))`.\n\nPlancherel gives the exact heat-defect identity\n\n`E_H ||delta_r V||_2^2 = int 2(1-exp(-|xi|^2/(2N^2))) |Vhat(xi)|^2 dxi`.\n\nBecause the strict transporter has `supp Vhat subset B_(N/4)`, put `x=|xi|^2/(2N^2)<=1/32`.  The elementary integral identity `1-e^(-x)=int_0^x e^(-s)ds` gives\n\n`e^(-1/32) x <= 1-e^(-x) <= x`.\n\nHence pointwise on the full resolved support\n\n`e^(-1/32)||grad V||_2^2 <= N^2 E_H||delta_r V||_2^2 <= ||grad V||_2^2`.\n\nAfter integrating a child lifetime, the positive heat-increment service\n\n`S_heat=N^3 int dt E_H||delta_r V||_2^2`\n\nsatisfies\n\n`e^(-1/32) D_V <= S_heat <= D_V`.\n\nThus the existing high-strain gate `D_V>=32 pi^2/(75c)` immediately forces a fixed positive **physical-space increment service**.  The Gaussian is not a replacement causal law: it is the heat semigroup already intrinsic to the viscous PDE, used only to expose where resolved gradient activity lives.\n\nThe preceding resolved-ancestor theorem says at least half the actual `D_V` law lies on shell-time atoms with critical mass `mu_j>=32pi^2/(75c^2)`.  Since the heat multiplier is pointwise between `e^(-1/32)` and `1` times the gradient multiplier, the same good shell-time set carries at least\n\n`(1/2)e^(-1/32) = {critical_ancestor_heat_service_fraction_lower():.12g}`\n\nof the **entire heat-increment service law**.  The frequency ancestor and spatial service are therefore simultaneous physical marks, not two unrelated pigeonholes.\n\nFor each good shell, each heat displacement `r`, and any normalized affine coherent window, Moyal gives an exact positive phase-space disintegration of `||delta_r P_jV||_2^2`.  Translation covariance\n\n`V_g(delta_r f)(X,k)=e^(-ik.r)V_gf(X-r,k)-V_gf(X,k)`\n\nshows that every cell atom is a real coherent edge between neighborhoods separated by the actual Brownian/heat displacement.  No coherent cell is selected by argmax and no global shell mass is divided among a guessed packet count.\n\nStress: `{out.samples}` spectral heat-defect/ancestor states and representative coherent identities\n- minimum spectral lower margin: `{out.minimum_spectral_lower_margin:.3e}`\n- minimum spectral upper margin: `{out.minimum_spectral_upper_margin:.3e}`\n- minimum retained-fraction margin above `e^(-1/32)/2`: `{out.minimum_retained_fraction_margin:.3e}`\n- worst relative discrete Moyal residual: `{out.worst_discrete_moyal_relative_residual:.3e}`\n- worst relative coherent increment-covariance residual: `{out.worst_increment_covariance_residual:.3e}`\n\nThis supplies the missing **spatial/coherent entrance** for the high-strain resolved-ancestor law.  What remains is material routing of these positive coherent edges: old--old service must enter the existing reservoir half-life, old--new edges must be physical relink/interface, and genuinely new--new edges must create coherent ancestry/service.  That routing must be proved for this dissipation-seeded measure rather than assumed from the SGS source theorem.  `D_V` remains an `O(1/N)` physical cost, not an additive reset.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
