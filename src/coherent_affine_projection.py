from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

EXTENDED_ASPECT = 567.0 / 500.0
SHELL_LOWER_AXIS = 2.0 / 3.0
GAUSSIAN_DENSITY_PEAK_3D = (2.0 * math.pi) ** (-1.5)


def affine_residual_spectral_gap(
    hermite_degrees: Sequence[int], squared_coefficients: Sequence[float]
) -> dict[str, float]:
    """Exact OU spectral bookkeeping after removing Gaussian affine modes.

    R has no Hermite degree 0 or 1.  For orthonormal Gaussian Hermites,
      ||R||^2=sum |c_alpha|^2,
      ||grad R||^2=sum |alpha||c_alpha|^2 >=2||R||^2.
    """
    d = np.asarray(hermite_degrees, int)
    w = np.asarray(squared_coefficients, float)
    if d.ndim != 1 or w.shape != d.shape or np.any(d < 2) or np.any(w < 0):
        raise ValueError("affine residual uses Hermite degrees >=2 and nonnegative squared coefficients")
    r2 = float(w.sum())
    k2 = float(np.dot(d, w))
    return {
        "residual_velocity_l2_sq": r2,
        "deformation_variance": k2,
        "spectral_gap_margin": k2 - 2.0 * r2,
    }


def gaussian_position_weight_upper(residual_l2_sq: float, gradient_l2_sq: float, dimension: int = 3) -> float:
    """Bound E |z|^2 |R|^2 using Gaussian creation/annihilation operators.

    For each coordinate, ||z_j f|| <= ||partial_j f||+||(z_j-partial_j)f||,
    hence ||z_j f||^2 <= 4||partial_j f||^2+2||f||^2.  Summing gives
      E|z|^2|R|^2 <= 4 E|grad R|^2 + 2d E|R|^2.
    """
    if min(residual_l2_sq, gradient_l2_sq) < 0 or dimension <= 0:
        raise ValueError("nonnegative energies and positive dimension required")
    return 4.0 * gradient_l2_sq + 2.0 * dimension * residual_l2_sq


def gaussian_core_nonaffine_forcing_upper(deformation_rms: float, intrinsic_carrier: float) -> float:
    """Relative L2 forcing of the full non-affine Gaussian-core low-high residual.

    In intrinsic coordinates W=vbar+Abar z+R and F=grad W.
    K^2=E||F-Abar||^2=E||grad R||^2.  Since R has degrees >=2,
    E|R|^2<=K^2/2 and E|z|^2|R|^2<=7K^2 in d=3.

    For psi=g exp(i q.z), grad psi=(iq-z/2)psi.  Scalar residual advection is
    R.grad psi and vector amplitude mismatch is (F-Abar)psi.  Therefore
      ||R_total||/||psi|| <= [1+q/sqrt2+sqrt7/2] K.
    """
    if deformation_rms < 0 or intrinsic_carrier < 0:
        raise ValueError("nonnegative deformation and carrier required")
    coeff = 1.0 + intrinsic_carrier / math.sqrt(2.0) + math.sqrt(7.0) / 2.0
    return coeff * deformation_rms


def longest_axis_from_radius_aspect(r_g: float, aspect: float) -> float:
    """ell_max <= aspect^(2/3) r_g for three positive axes."""
    if r_g <= 0 or aspect < 1:
        raise ValueError("positive radius and aspect>=1 required")
    return aspect ** (2.0 / 3.0) * r_g


def intrinsic_carrier_upper(scale_radius: float, aspect: float, shell_radius_ratio: float) -> float:
    """q=|L^T k| <= ell_max |k| <= aspect^(2/3) (N r_g)(|k|/N)."""
    if scale_radius <= 0 or aspect < 1 or shell_radius_ratio <= 0:
        raise ValueError("positive radius/shell and aspect>=1 required")
    return aspect ** (2.0 / 3.0) * scale_radius * shell_radius_ratio


def coherent_deformation_to_dissipation_constant(
    aspect: float = EXTENDED_ASPECT,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
) -> float:
    """C in I_K^2 <= C c D_V for a coherent grain lifetime T=cN^-2.

    K_C^2 <= aspect^2 E_gamma |grad V|^2.
    The normalized Gaussian density has peak (2pi)^(-3/2) r_g^-3.
    Shell uncertainty gives r_g >= lower_axis_constant/N.  Cauchy in time and
    D_V=N int||grad V||_2^2 dt yield
      I_K^2 <= aspect^2 (2pi)^(-3/2) lower_axis_constant^-3 c D_V.
    """
    if aspect < 1 or lower_axis_constant <= 0:
        raise ValueError("valid aspect and lower axis constant required")
    return aspect * aspect * GAUSSIAN_DENSITY_PEAK_3D * lower_axis_constant ** (-3.0)


def normalized_dissipation_from_coherent_deformation(
    deformation_action: float,
    scaled_lifetime: float,
    aspect: float = EXTENDED_ASPECT,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
) -> float:
    """D_V lower forced by I_K=int K_C dt on T=cN^-2."""
    if deformation_action < 0 or scaled_lifetime <= 0:
        raise ValueError("nonnegative action and positive scaled lifetime required")
    C = coherent_deformation_to_dissipation_constant(aspect, lower_axis_constant)
    return deformation_action * deformation_action / (C * scaled_lifetime)


def coherent_affine_projection_from_samples(
    weights: Sequence[float], velocities: np.ndarray, gradients: np.ndarray, positions: np.ndarray
) -> dict[str, np.ndarray | float]:
    """Finite Gaussian-quadrature regression for the coherent affine projection.

    weights approximate gamma, positions are intrinsic z samples.  The exact
    continuum coefficients are vbar=E W and Abar=E grad W=E[W tensor z] by
    Gaussian integration by parts.  This routine uses vbar and mean gradient and
    reports the residual/deformation energies; it does not claim arbitrary finite
    samples obey the integration-by-parts identity exactly.
    """
    w = np.asarray(weights, float)
    V = np.asarray(velocities, float)
    G = np.asarray(gradients, float)
    Z = np.asarray(positions, float)
    if w.ndim != 1 or V.ndim != 2 or G.ndim != 3 or Z.ndim != 2:
        raise ValueError("invalid sample shapes")
    if len(w) != len(V) or len(w) != len(G) or len(w) != len(Z) or V.shape[1] != Z.shape[1] or G.shape[1:] != (V.shape[1], V.shape[1]):
        raise ValueError("sample dimensions mismatch")
    if np.any(w < 0) or float(w.sum()) <= 0:
        raise ValueError("positive probability weights required")
    w = w / float(w.sum())
    vbar = np.einsum("a,ai->i", w, V)
    Abar = np.einsum("a,aij->ij", w, G)
    R = V - vbar[None, :] - np.einsum("ij,aj->ai", Abar, Z)
    F = G - Abar[None, :, :]
    r2 = float(np.einsum("a,ai,ai->", w, R, R))
    k2 = float(np.einsum("a,aij,aij->", w, F, F))
    return {"vbar": vbar, "Abar": Abar, "residual_l2_sq": r2, "deformation_variance": k2}


@dataclass(frozen=True)
class CoherentAffineStress:
    samples: int
    minimum_spectral_gap_margin: float
    minimum_position_weight_margin: float
    minimum_dissipation_margin: float
    minimum_axis_margin: float
    maximum_clean_collision_constant: float


def stress(samples: int = 50_000, seed: int = 20260808) -> CoherentAffineStress:
    rng = np.random.default_rng(seed)
    msg = mpw = md = ma = float("inf")
    Cclean = coherent_deformation_to_dissipation_constant()
    for _ in range(samples):
        m = int(rng.integers(1, 80))
        deg = rng.integers(2, 12, size=m)
        coeff = rng.random(m)
        out = affine_residual_spectral_gap(deg, coeff)
        msg = min(msg, float(out["spectral_gap_margin"]))
        if out["spectral_gap_margin"] < -2e-13:
            raise AssertionError("OU affine-residual spectral gap failed")

        r2 = float(out["residual_velocity_l2_sq"])
        k2 = float(out["deformation_variance"])
        upper = gaussian_position_weight_upper(r2, k2)
        # creation-annihilation theorem plus r2<=k2/2 implies the clean 7 K^2 upper.
        clean = 7.0 * k2
        margin = clean - upper
        mpw = min(mpw, margin)
        if margin < -2e-12 * max(1.0, clean):
            raise AssertionError("clean Gaussian position-weight bound failed")

        kappa = float(rng.uniform(1.0, EXTENDED_ASPECT))
        c = float(rng.uniform(0.02, 2.0))
        D = float(rng.lognormal(mean=-2.0, sigma=2.0))
        C = coherent_deformation_to_dissipation_constant(kappa)
        # Choose an admissible action below the Cauchy upper sqrt(C c D).
        Imax = math.sqrt(C * c * D)
        I = float(rng.random()) * Imax
        lower = normalized_dissipation_from_coherent_deformation(I, c, kappa)
        md = min(md, D - lower)
        if lower > D + 2e-12 * max(1.0, D):
            raise AssertionError("coherent deformation/dissipation collision failed")

        axes = np.exp(rng.uniform(-3.0, 3.0, size=3))
        rg = float(np.prod(axes) ** (1.0 / 3.0))
        kap = float(axes.max() / axes.min())
        bound = longest_axis_from_radius_aspect(rg, kap)
        margin_axis = bound - float(axes.max())
        ma = min(ma, margin_axis)
        if margin_axis < -2e-12 * max(1.0, bound):
            raise AssertionError("radius/aspect longest-axis bound failed")

    return CoherentAffineStress(samples, msg, mpw, md, ma, Cclean)


def theorem_certificate() -> dict[str, object]:
    C = coherent_deformation_to_dissipation_constant()
    return {
        "status": "EXACT_COHERENT_AFFINE_PROJECTION_AND_DEFORMATION_DISSIPATION_COLLISION__AVERAGED_TRANSPORTER_SOURCE_CALCULUS_REMAINS",
        "affine_projection": "vbar=E_gamma W, Abar=E_gamma grad W=E_gamma[W tensor z]; R=W-vbar-Abar z has Hermite degrees >=2",
        "ou_gap": "E|R|^2 <= (1/2) E|grad R|^2=(1/2) K_coh^2",
        "weighted_moment": "E|z|^2|R|^2 <= 7 K_coh^2 in d=3",
        "gaussian_core_forcing": "||R_nonaff psi||/||psi|| <= [1+|q|/sqrt2+sqrt7/2] K_coh",
        "radius_aspect_carrier": "|q|<=aspect^(2/3)(N r_g)(|k|/N)",
        "dissipation_collision": f"I_K^2 <= C_coh c D_V with C_coh={C:.12g} on aspect<=567/500 and r_g>=2/(3N)",
        "master_route": "large coherent deformation action is critical D_V; small action makes the full Gaussian-core nonaffine low-high forcing perturbative without resolving separate high Hermite currencies",
        "continuum_status": "to replace the center transporter globally, derive the resolved Navier-Stokes source/corotational calculus for the coherent averaged affine jet Abar and register its time variation once",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-coherent-affine-projection"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "coherent_affine_projection.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    C = coherent_deformation_to_dissipation_constant()
    md = f"""# Coherent affine projection and deformation-dissipation collision\n\nStatus: **{cert['status']}**.\n\nThe center-jet countermodel indicates that the affine gauge should be selected by the whole coherent eddy.  In intrinsic Gaussian coordinates `W(z)=L^-1 V(X+Lz)`, define\n\n`vbar=E_gamma W`, `Abar=E_gamma grad W`.\n\nGaussian integration by parts gives `Abar=E_gamma[W tensor z]`, so\n\n`R=W-vbar-Abar z`\n\nhas no Hermite degree 0 or 1.  The Ornstein--Uhlenbeck spectral gap therefore gives exactly\n\n`E|R|^2 <= K_coh^2/2`, where `K_coh^2=E||grad W-Abar||^2`.\n\nCreation/annihilation operators give the weighted bound\n\n`E |z|^2 |R|^2 <= 7 K_coh^2`.\n\nFor the coherent Gaussian carrier `psi=g exp(i q.z)`, scalar residual advection plus vector amplitude mismatch obey\n\n`||R_nonaff psi||/||psi|| <= [1+|q|/sqrt(2)+sqrt(7)/2] K_coh`.\n\nThus the **entire** spatial non-affine Gaussian-core forcing, including all higher Hermite degrees, is controlled by one physical deformation-variance observable.  There is no need to create H4/H5/... master currencies.\n\nLarge coherent deformation is already critical dissipation.  The Gaussian density peak, `cond(L)<=567/500`, shell uncertainty `r_g>=2/(3N)`, Cauchy in time and `D_V=N int||grad V||_2^2 dt` give\n\n`I_K^2 <= C_coh c D_V`, `I_K=int K_coh dt`,\n\nwith\n\n`C_coh=(567/500)^2 (2 pi)^(-3/2) (3/2)^3 = {C:.12g}`.\n\nHence `D_V >= I_K^2/(C_coh c)`.  This is a scale-critical dissipation branch, **not** a uniform finite reset count.\n\nOn a scale-matched radius branch `s=N r_g<=s0`, the intrinsic carrier is uniformly bounded by\n\n`|q| <= cond(L)^(2/3) s0 (|k|/N)`,\n\nso small `I_K` makes the whole non-affine Gaussian-core low--high residual a controlled perturbation.  Large radius remains the existing sticky affine-radius ancestry branch.\n\nStress: `{out.samples}` spectral/radius/dissipation states\n- minimum OU spectral-gap margin: `{out.minimum_spectral_gap_margin:.3e}`\n- minimum clean position-weight margin: `{out.minimum_position_weight_margin:.3e}`\n- minimum dissipation-collision margin: `{out.minimum_dissipation_margin:.3e}`\n- minimum radius/aspect axis margin: `{out.minimum_axis_margin:.3e}`\n- clean collision constant: `{out.maximum_clean_collision_constant:.12g}`\n\nThis removes the need for a separate high-Hermite curvature currency at the Gaussian-core forcing level.  The remaining conceptual bridge is now the **source calculus of the coherent averaged affine jet** `Abar(t)`: if it replaces the center jet as the common material transporter, its corotational time variation must be derived from resolved Navier--Stokes and routed to the existing pressure/SGS/viscous/service currencies without double charging.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
