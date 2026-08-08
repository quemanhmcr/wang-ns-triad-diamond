from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.affine_coherent_bessel import CLEAN_BESSEL
from src.coherent_affine_projection import EXTENDED_ASPECT, SHELL_LOWER_AXIS

P_YOUNG = 1.5
Q_DUAL = 3.0
DEFAULT_PROFILE_ERROR = 0.01
DEFAULT_LOG_COV_RADIUS = 0.4
DEFAULT_RADIUS_CAP = 4.0
ROOT_GRID_COLORS = 5**6
SHELL_LOG_HALFWIDTH = 2.0 / 25.0
SCALE_BIN_WIDTH = 2.0 / 25.0
SCALE_COLORS = 4


def dual_gaussian_pairing_lower(log_cov_radius: float) -> float:
    """Uniform L^(3/2)-L^3 pairing of a Gaussian and quantized dual Gaussian.

    Write a normalized frequency Gaussian as G_C=A_C exp(-xi.C.xi/2), and the
    L^3-normalized dual attached to covariance representative D as
    h_D=B_D exp(-xi.D.xi/4).  If every eigenvalue of log(C^-1/2 D C^-1/2) has
    absolute value <=delta, then

      <G_C,h_D> = prod_i [(3/2) r_i^(1/3)/(1+r_i/2)]^(1/2).

    Frobenius log-distance <=delta implies the eigenvalue bound used here.
    """
    d = float(log_cov_radius)
    if d < 0 or not math.isfinite(d):
        raise ValueError("finite nonnegative log-covariance radius required")
    one = []
    for sign in (-1.0, 1.0):
        r = math.exp(sign * d)
        one.append(math.sqrt(1.5 * r ** (1.0 / 3.0) / (1.0 + 0.5 * r)))
    return min(one) ** 3


def dual_probe_l2_norm_sq(profile_radius_rep: float) -> float:
    """L2 norm squared of the L3-normalized Gaussian dual probe.

    For profile physical geometric radius r_g (the profile covariance, not the
    half-covariance of the L2 probe itself), ||h||_2^2=3 sqrt(pi)/(2 r_g).
    """
    if profile_radius_rep <= 0 or not math.isfinite(profile_radius_rep):
        raise ValueError("positive finite representative profile radius required")
    return 3.0 * math.sqrt(math.pi) / (2.0 * profile_radius_rep)


def representative_radius_ratio_lower(log_cov_radius: float) -> float:
    """If d_log(Sigma_rep,Sigma)<=delta in 3D, r_rep/r >=e^(-sqrt3 delta/6)."""
    if log_cov_radius < 0:
        raise ValueError("nonnegative log-covariance radius required")
    return math.exp(-math.sqrt(3.0) * log_cov_radius / 6.0)


def dual_probe_critical_mass_lower(
    profile_error: float = DEFAULT_PROFILE_ERROR,
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    scale_radius_lower: float = SHELL_LOWER_AXIS,
) -> float:
    """Lower N |<f,phi_rep>|^2 for an actual near-Gaussian Young role.

    ||f-G||_(3/2)<=eps, ||G||_(3/2)=1.  h_rep is L3 normalized, so
      |<f,h_rep>| >= pairing(G,h_rep)-eps.
    Normalize phi=h_rep/||h_rep||_2.  The shell gives N r_g>=s_min and covariance
    quantization gives r_rep>=r_g exp(-sqrt3 delta/6).  Therefore

      N|<f,phi>|^2 >= [2/(3sqrt(pi))] (P_delta-eps)^2
                       s_min exp(-sqrt3 delta/6).
    """
    eps = float(profile_error)
    d = float(log_cov_radius)
    smin = float(scale_radius_lower)
    if not (0.0 <= eps < 1.0) or d < 0 or smin <= 0:
        raise ValueError("invalid dual-probe registration data")
    pairing = dual_gaussian_pairing_lower(d)
    if pairing <= eps:
        return 0.0
    return (
        2.0
        / (3.0 * math.sqrt(math.pi))
        * (pairing - eps) ** 2
        * smin
        * representative_radius_ratio_lower(d)
    )


def scale_bin_index(scale: float, bin_width: float = SCALE_BIN_WIDTH) -> int:
    if scale <= 0 or bin_width <= 0:
        raise ValueError("positive scale/bin width required")
    return math.floor(math.log(scale) / bin_width)


def scale_bin_reference(bin_index: int, bin_width: float = SCALE_BIN_WIDTH) -> float:
    if bin_width <= 0:
        raise ValueError("positive bin width required")
    return math.exp((bin_index + 0.5) * bin_width)


def scale_bin_shell_union(bin_index: int, bin_width: float = SCALE_BIN_WIDTH, shell_halfwidth: float = SHELL_LOG_HALFWIDTH) -> tuple[float, float]:
    """Frequency support union for all scales N in one logarithmic bin."""
    if bin_width <= 0 or shell_halfwidth < 0:
        raise ValueError("valid bin/shell widths required")
    lo = math.exp(bin_index * bin_width - shell_halfwidth)
    hi = math.exp((bin_index + 1) * bin_width + shell_halfwidth)
    return lo, hi


def same_color_scale_shells_are_disjoint(
    color_period: int = SCALE_COLORS,
    bin_width: float = SCALE_BIN_WIDTH,
    shell_halfwidth: float = SHELL_LOG_HALFWIDTH,
) -> bool:
    """Bins whose indices differ by color_period have disjoint shell unions."""
    if color_period <= 0:
        raise ValueError("positive scale color period required")
    return (color_period - 1) * bin_width > 2.0 * shell_halfwidth


def normalized_covariance_eigenvalue_bounds(
    radius_cap: float = DEFAULT_RADIUS_CAP,
    aspect_cap: float = EXTENDED_ASPECT,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
    bin_width: float = SCALE_BIN_WIDTH,
) -> tuple[float, float]:
    """Eigenvalue bounds for N_b^2 Sigma inside one logarithmic scale bin.

    If N/N_b lies in [e^-h/2,e^h/2], then (N_b/N)^2 lies in [e^-h,e^h].
    Combine this with the shell/radius/aspect bounds on N^2 Sigma.
    """
    if radius_cap <= 0 or aspect_cap < 1 or lower_axis_constant <= 0 or bin_width <= 0:
        raise ValueError("invalid radius/aspect/bin data")
    eig_lo = math.exp(-bin_width) * lower_axis_constant**2
    eig_hi = math.exp(bin_width) * aspect_cap ** (4.0 / 3.0) * radius_cap**2
    return eig_lo, eig_hi


def frobenius_radius_for_affine_log_radius(log_cov_radius: float, eigenvalue_lower: float) -> float:
    """Euclidean SPD tolerance ensuring affine-invariant log distance <=delta.

    If A,B>=mI and ||A-B||_F<=eps, then
      E=A^-1/2(B-A)A^-1/2, ||E||_F<=eps/m.
    For eps/m=x<1, ||log(I+E)||_F<=x/(1-x).
    Choosing x=delta/(1+delta) makes this <=delta.
    """
    d = float(log_cov_radius)
    m = float(eigenvalue_lower)
    if d <= 0 or m <= 0:
        raise ValueError("positive metric radius and eigenvalue lower bound required")
    return m * d / (1.0 + d)


def affine_log_distance_upper_from_frobenius(frobenius_distance: float, eigenvalue_lower: float) -> float:
    eps = float(frobenius_distance)
    m = float(eigenvalue_lower)
    if eps < 0 or m <= 0 or eps >= m:
        raise ValueError("require 0<=eps<m")
    x = eps / m
    return x / (1.0 - x)


def covariance_cover_number_upper(
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    aspect_cap: float = EXTENDED_ASPECT,
    lower_axis_constant: float = SHELL_LOWER_AXIS,
) -> int:
    """Crude finite cover whose cells have affine-invariant SPD radius <=delta.

    The compact SPD set lies in the six-dimensional Frobenius ball of radius
    sqrt(3) M.  A maximal eps-separated subset has cardinality at most
    (1+2R/eps)^6.  Taking eps=m delta/(1+delta) guarantees every representative
    is within affine-invariant log distance delta by the preceding lemma.
    """
    d = float(log_cov_radius)
    if d <= 0:
        raise ValueError("positive covariance quantization radius required")
    m, M = normalized_covariance_eigenvalue_bounds(radius_cap, aspect_cap, lower_axis_constant)
    eps = frobenius_radius_for_affine_log_radius(d, m)
    R = math.sqrt(3.0) * M
    return int(math.ceil((1.0 + 2.0 * R / eps) ** 6))


def phase_space_color_count(separation: int = 4, dimension: int = 6) -> int:
    """Color unit cells in one representative dual-probe phase coordinate.

    Residues modulo separation+1 in each coordinate suffice even after allowing
    points to move within their unit cells.  This coloring is used only for the
    analysis Bessel budget; it is not asserted to be the canonical material label.
    """
    if separation <= 0 or dimension <= 0:
        raise ValueError("positive separation/dimension required")
    return (separation + 1) ** dimension


def effective_root_frame_budget(
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    aspect_cap: float = EXTENDED_ASPECT,
) -> float:
    """Finite Bessel budget after covariance binning and canonical zeta coloring."""
    if not same_color_scale_shells_are_disjoint():
        raise AssertionError("default logarithmic scale coloring does not separate outer shells")
    bins = covariance_cover_number_upper(log_cov_radius, radius_cap, aspect_cap)
    return float(SCALE_COLORS * bins * phase_space_color_count(4, 6)) * float(CLEAN_BESSEL)


def registered_root_count_upper(
    global_energy: float,
    root_scale_upper: float,
    profile_error: float = DEFAULT_PROFILE_ERROR,
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    aspect_cap: float = EXTENDED_ASPECT,
) -> float:
    """Count distinct registered root cells using only actual analysis coefficients."""
    if global_energy < 0 or root_scale_upper <= 0:
        raise ValueError("valid global energy/root scale required")
    eta = dual_probe_critical_mass_lower(profile_error, log_cov_radius)
    if eta <= 0:
        raise ValueError("dual-probe mass lower is not positive")
    P = effective_root_frame_budget(log_cov_radius, radius_cap, aspect_cap)
    return P * global_energy * root_scale_upper / eta


def renyi_action_lower_with_registered_roots(
    depth: int,
    global_energy: float,
    base_scale: float,
    profile_error: float = DEFAULT_PROFILE_ERROR,
    log_cov_radius: float = DEFAULT_LOG_COV_RADIUS,
    radius_cap: float = DEFAULT_RADIUS_CAP,
    aspect_cap: float = EXTENDED_ASPECT,
) -> float:
    """Existing binary/Renyi root action with the dual-probe root budget inserted."""
    if depth < 0 or global_energy <= 0 or base_scale <= 0:
        raise ValueError("valid depth/energy/base scale required")
    root_scale_upper = base_scale * (25.0 / 24.0) ** depth
    n0 = registered_root_count_upper(
        global_energy,
        root_scale_upper,
        profile_error,
        log_cov_radius,
        radius_cap,
        aspect_cap,
    )
    return depth * math.log(2.0) - math.log(max(n0, 1.0))


def divergence_free_projection_pairing_residual(u: np.ndarray, phi: np.ndarray, P: np.ndarray) -> float:
    """Finite-dimensional model of <u,P phi>=<u,phi> for Pu=u, P=P*=P^2."""
    u = np.asarray(u, complex)
    phi = np.asarray(phi, complex)
    P = np.asarray(P, complex)
    if P.ndim != 2 or P.shape[0] != P.shape[1] or u.shape != phi.shape or u.shape != (P.shape[0],):
        raise ValueError("dimension mismatch")
    return float(abs(np.vdot(u, P @ phi) - np.vdot(u, phi)))


@dataclass(frozen=True)
class DualGaussianRootStress:
    samples: int
    minimum_pairing_margin: float
    minimum_root_mass_margin_over_one_fifth: float
    minimum_cover_margin: float
    worst_projection_pairing_residual: float
    covariance_bins: int
    effective_frame_budget: float


def stress(samples: int = 50_000, seed: int = 20260808) -> DualGaussianRootStress:
    rng = np.random.default_rng(seed)
    mp = mm = mc = float("inf")
    wp = 0.0
    d = DEFAULT_LOG_COV_RADIUS
    lower = dual_gaussian_pairing_lower(d)
    eta = dual_probe_critical_mass_lower()
    if eta <= 0.2:
        raise AssertionError("default dual-probe root mass does not beat 1/5")
    for _ in range(samples):
        # Exact determinant formula for a diagonal covariance mismatch inside the log ball.
        lam = rng.uniform(-d, d, size=3)
        r = np.exp(lam)
        exact = float(np.prod(np.sqrt(1.5 * r ** (1.0 / 3.0) / (1.0 + 0.5 * r))))
        mp = min(mp, exact - lower)
        if exact + 2e-14 < lower:
            raise AssertionError("dual Gaussian pairing lower failed")

        eps = float(rng.uniform(0.0, DEFAULT_PROFILE_ERROR))
        dd = float(rng.uniform(0.0, DEFAULT_LOG_COV_RADIUS))
        mass = dual_probe_critical_mass_lower(eps, dd)
        # The default endpoint is the weakest over the sampled rectangle.
        mm = min(mm, mass - eta)
        if mass + 2e-14 < eta:
            raise AssertionError("dual-probe mass monotonic endpoint lower failed")

        # Four logarithmic scale colors make same-color outer shell subspaces disjoint.
        if not same_color_scale_shells_are_disjoint():
            raise AssertionError("scale shell coloring failed")
        b = int(rng.integers(-30, 30))
        lo1, hi1 = scale_bin_shell_union(b)
        lo2, hi2 = scale_bin_shell_union(b + SCALE_COLORS)
        if not (hi1 < lo2):
            raise AssertionError("same-color scale shell unions overlap")

        # The Euclidean SPD net radius rigorously implies affine-log radius <=d.
        m_eig, _ = normalized_covariance_eigenvalue_bounds()
        eps_cov = frobenius_radius_for_affine_log_radius(d, m_eig)
        aff_upper = affine_log_distance_upper_from_frobenius(eps_cov, m_eig)
        margin = d - aff_upper
        mc = min(mc, margin)
        if margin < -3e-14:
            raise AssertionError("Frobenius covariance cover did not imply affine-log radius")

        # Projection pairing on an actual divergence-free vector.
        n = int(rng.integers(2, 8))
        rank = int(rng.integers(1, n + 1))
        Qp, _ = np.linalg.qr(rng.normal(size=(n, rank)) + 1j * rng.normal(size=(n, rank)))
        P = Qp @ Qp.conj().T
        coeff = rng.normal(size=rank) + 1j * rng.normal(size=rank)
        u = Qp @ coeff
        phi = rng.normal(size=n) + 1j * rng.normal(size=n)
        res = divergence_free_projection_pairing_residual(u, phi, P)
        wp = max(wp, res)
        if res > 2e-11 * max(1.0, np.linalg.norm(u) * np.linalg.norm(phi)):
            raise AssertionError("Leray-projected dual probe changed divergence-free coefficient")

    return DualGaussianRootStress(
        samples,
        mp,
        mm,
        mc,
        wp,
        covariance_cover_number_upper(),
        effective_root_frame_budget(),
    )


def theorem_certificate() -> dict[str, object]:
    eta = dual_probe_critical_mass_lower()
    bins = covariance_cover_number_upper()
    P = effective_root_frame_budget()
    return {
        "status": "EXACT_DUAL_GAUSSIAN_ANALYSIS_ROOT_QUANTUM__TRANSFER_CELL_ALIGNMENT_REMAINS",
        "duality": "||f-G||_(3/2)<=eps and L3-normalized Gaussian dual h imply |<f,h>|>=<G,h>-eps",
        "covariance_pairing": f"d_log<=0.4 gives Gaussian primal-dual pairing >={dual_gaussian_pairing_lower(0.4):.12g}",
        "root_quantum": f"with eps=1/100 and N r_g>2/3, quantized dual probe gives N|<f,phi>|^2>={eta:.12g}>1/5",
        "leray": "for divergence-free u, replacing a vector Gaussian probe by its Leray projection preserves the coefficient exactly",
        "covariance_cover": f"inside one log-scale bin, transition-aspect 2/3<Nr_g<=4 covariances rescaled by the bin reference need at most {bins} affine-log-SPD bins of radius 0.4",
        "scale_coloring": "log-scale bins of width 2/25 split into 4 colors; same-color outer shell projectors are frequency-disjoint and hence orthogonal",
        "phase_space_coloring": "within each scale/covariance representative, unit cells in that dual-probe phase coordinate split into 5^6 colors; one color is 4-separated for the Bessel budget only",
        "bessel": f"within one scale bin/covariance/color the equal-covariance analysis budget is 25/4; orthogonality across same-color scale bins makes the total uniform P_eff={P:.12g}",
        "causal_effect": "the huge P_eff changes only the finite logarithmic root offset; the L log(48/25) reuse slope is unchanged",
        "important_scope": "this is an actual analysis coefficient of the selected orthogonal Fourier/helical role; project the dual probe by the same outer role projector so no other field component can cancel it. It is not a synthesis coefficient or an L2-closeness claim",
        "continuum_status": "remaining registration is transfer-cell/material-label alignment: associate the complex Gaussian parent mark with the actual transfer-selected causal root cell or route misaligned work to existing physical currencies",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-dual-gaussian-root-registration"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "dual_gaussian_root_registration.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    eta = dual_probe_critical_mass_lower()
    md = f"""# Dual-Gaussian root registration: Christ proximity gives actual analysis energy\n\nStatus: **{cert['status']}**.\n\nA phase-aligned complex Gaussian near-profile should not be inserted into the velocity as a fictitious synthesis component.  Use it only to design a dual Gaussian **analysis probe**.  Christ supplies the magnitude Gaussian; the present theorem is conditional on the separate phase/polarization control being strong enough to lift that mark to complex `L^(3/2)` proximity.\n\nLet `||f||_(3/2)=||G||_(3/2)=1`, `||f-G||_(3/2)<=eps`.  The exact L3-dual Gaussian is `h_G=|G|^(-1/2)G`, with `||h_G||_3=1` and `<G,h_G>=1`.  Quantize the profile covariance to a representative within log-SPD radius `delta`, and use the corresponding L3-normalized dual `h_rep`.  If `delta<=0.4`, exact Gaussian integration gives\n\n`<G,h_rep> >= {dual_gaussian_pairing_lower(0.4):.12g}`.\n\nTherefore `|<f,h_rep>| >= <G,h_rep>-eps`.  After L2-normalizing `phi=h_rep/||h_rep||_2`,\n\n`||h_rep||_2^2 = 3 sqrt(pi)/(2 r_g,rep)`.\n\nThe shell lower axis and covariance quantization give the scale-independent actual coefficient bound\n\n`N |<f,phi>|^2 >= {eta:.12g}`\n\nat `eps=1/100`, `delta=0.4`.  In particular it is **strictly larger than 1/5**, the clean critical root quantum used by the causal/Renyi modules.\n\nThis is a coefficient of the actual selected role.  If that role is obtained by an exact self-adjoint outer Fourier/helical projector `Q`, use the probe `Q phi`; then `<u,Q phi>=<Qu,phi>` exactly and `||Q phi||_2<=1`, so normalization cannot reduce the coefficient.  For divergence-free roles, Leray projection is likewise coefficient-preserving.  No L2 closeness of the Christ remainder is asserted or needed.\n\nVariable root scale also does not destroy the energy count.  Put `log N` in bins of width `2/25` and color the bins modulo `4`.  Because the outer role shell halfwidth is `2/25`, two distinct bins of one color have disjoint physical Fourier support, so their exact outer role projectors are orthogonal.  Within one scale bin use its reference `N_b`; the rescaled covariances `N_b^2 Sigma` lie in a fixed compact subset of six-dimensional `Sym(3)`.  A crude Frobenius net, chosen fine enough to guarantee affine-invariant log-SPD radius `delta=0.4`, uses at most `{out.covariance_bins}` bins.  For each fixed scale/covariance representative, unit cells in that representative dual-probe phase coordinate may be colored with `5^6={ROOT_GRID_COLORS}` colors so cells of one color are 4-separated.  This auxiliary coloring is only an analysis-budget device; it is not identified with the canonical material label.  Inside one scale bin/covariance/color the exact affine coherent Bessel theorem gives analysis budget `25/4`; orthogonality across same-color scale bins prevents any factor growing with causal depth.  Thus all registered roots have one finite effective budget\n\n`P_eff <= {out.effective_frame_budget:.12g}`.\n\nThis constant is intentionally huge but **scale independent**.  In the causal root estimate it enters only through `log P_eff`, so it changes the finite depth offset and not the linear reuse slope `log(48/25)`.\n\nThe theorem therefore closes a subtle gap **once complex phase-aligned `L^(3/2)` Gaussian proximity is available**: such profile information produces a scale-critical quantum in the actual `L2` energy analysis, via duality, without pretending `L^(3/2)` closeness implies `L2` closeness.\n\nStress: `{out.samples}` covariance/pairing/Leray states\n- minimum Gaussian pairing margin: `{out.minimum_pairing_margin:.3e}`\n- minimum root-mass margin above the default endpoint: `{out.minimum_root_mass_margin_over_one_fifth:.3e}`\n- minimum covariance-cover radius margin: `{out.minimum_cover_margin:.3e}`\n- worst Leray coefficient residual: `{out.worst_projection_pairing_residual:.3e}`\n\nTwo registration issues remain.  First, the current inverse-Young ledger directly supplies Gaussian proximity for magnitudes; the separate phase/polarization theorems must be assembled into the complex phase-aligned proximity assumed here.  Second,  A causal parent slot already has a transfer-selected material coherent label.  The dual probe above is centered at the Christ Gaussian mark.  The final assembly must prove that a fixed transfer-weighted fraction of roots are aligned with that mark (or else the misaligned physical work is already cross/relink/backscatter/service currency).  Measurable selection of a Christ mark is not the hard part; this **transfer-cell alignment** is.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
