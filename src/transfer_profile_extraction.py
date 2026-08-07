from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


def trilinear_replacement_loss(epsilon: float) -> float:
    """Normalized transfer loss when each of three unit inputs moves by eps.

    For a trilinear form with sharp operator norm one, if ||f_j||=1 and
    ||f_j-G_j||<=eps, telescoping gives
      eps + (1+eps)eps + (1+eps)^2 eps.
    """
    if epsilon < 0:
        raise ValueError("epsilon must be nonnegative")
    return 3.0 * epsilon + 3.0 * epsilon * epsilon + epsilon ** 3


def gaussian_profile_transfer_lower(block_efficiency: float, gaussian_distance: float) -> float:
    if not (0.0 <= block_efficiency <= 1.0 + 1e-12):
        raise ValueError("block efficiency must lie in [0,1]")
    return block_efficiency - trilinear_replacement_loss(gaussian_distance)


def shell_critical_l2_mass_lower(lp_mass_inside_shell: float, shell_volume_coefficient: float) -> float:
    """Lower bound for N ||f||_2^2 from shell L^{3/2} mass.

    If Omega_N has volume <= C_Omega N^3, Holder gives
      ||f 1_Omega||_{3/2} <= |Omega|^{1/6} ||f||_2,
    hence N ||f||_2^2 >= C_Omega^{-1/3} ||f 1_Omega||_{3/2}^2.
    Plancherel constants are understood in the unitary Fourier convention.
    """
    if lp_mass_inside_shell < 0 or shell_volume_coefficient <= 0:
        raise ValueError("invalid shell mass parameters")
    return lp_mass_inside_shell ** 2 / (shell_volume_coefficient ** (1.0 / 3.0))


@dataclass(frozen=True)
class OneShotProfileCertificate:
    block_efficiency: float
    gaussian_distance: float
    replacement_loss: float
    gaussian_transfer_lower: float
    gaussian_lp_mass_inside_shell_lower: float
    gaussian_critical_l2_mass_lower: float


def one_shot_profile_certificate(
    block_efficiency: float,
    gaussian_distance: float,
    shell_volume_coefficient: float,
) -> OneShotProfileCertificate:
    """Deterministic algebra after applying Christ's inverse Young theorem.

    If each normalized role f_j is supported in its frequency block Omega_j and
    Christ supplies ||f_j-G_j||_{3/2}<=eps, then G_j has at least 1-eps of
    L^{3/2} norm inside Omega_j.  This routine records the exact transfer and
    critical-L2 consequences; the inverse-Young modulus delta_Christ(eps) is an
    external analytic theorem, not numerically invented here.
    """
    eps = gaussian_distance
    if not (0 <= eps < 1):
        raise ValueError("gaussian distance must be in [0,1)")
    loss = trilinear_replacement_loss(eps)
    transfer = block_efficiency - loss
    lp_inside = max(0.0, 1.0 - eps)
    critical = shell_critical_l2_mass_lower(lp_inside, shell_volume_coefficient)
    return OneShotProfileCertificate(block_efficiency, eps, loss, transfer, lp_inside, critical)


def finite_energy_shell_lp_upper(critical_l2_mass: float, shell_volume_coefficient: float) -> float:
    """Reverse presentation of the same Holder bridge.

    ||f||_{3/2} <= C_Omega^{1/6} (N||f||_2^2)^{1/2}.
    """
    if critical_l2_mass < 0 or shell_volume_coefficient <= 0:
        raise ValueError("invalid mass parameters")
    return shell_volume_coefficient ** (1.0 / 6.0) * math.sqrt(critical_l2_mass)




def arb_narrow_shell_mass_certificate() -> dict[str, str]:
    """Rigorous critical-mass bridge for the certified log shell |log(|xi|/N)|<=2/25.

    At one-percent Christ profile distance, Holder on the full spherical shell
    gives a universal lower bound N||G||_2^2 > 3/4 in the unitary Fourier
    convention.
    """
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover - exercised in Actions
        raise RuntimeError("python-flint is required for the rigorous shell-mass certificate") from exc
    ctx.prec = 160
    sigma = arb(2) / 25
    eps = arb(1) / 100
    pi = arb.pi()
    c_omega = (4 * pi / 3) * ((3 * sigma).exp() - (-3 * sigma).exp())
    mass = (1 - eps) ** 2 / c_omega.root(3)
    if not (mass > arb(3) / 4):
        raise AssertionError(f"narrow-shell critical mass certificate failed: {mass}")
    return {
        "shell_halfwidth": "2/25",
        "profile_distance": "1/100",
        "shell_volume_coefficient_ball": str(c_omega),
        "critical_l2_mass_ball": str(mass),
        "critical_l2_mass_lower_bound": "3/4",
        "status": "CERTIFIED",
    }


def stress(samples: int = 100_000, seed: int = 20260807) -> dict[str, float]:
    import numpy as np

    rng = np.random.default_rng(seed)
    worst_transfer_margin = float("inf")
    worst_mass = float("inf")
    for _ in range(samples):
        eps = float(rng.uniform(1e-8, 0.03))
        # Choose block deficit below the scale of the profile perturbation.
        block_eff = 1.0 - float(rng.uniform(0.0, 0.02))
        cvol = float(rng.uniform(0.5, 50.0))
        cert = one_shot_profile_certificate(block_eff, eps, cvol)
        # Direct telescoping polynomial identity.
        direct = eps + (1 + eps) * eps + (1 + eps) ** 2 * eps
        if abs(direct - cert.replacement_loss) > 1e-14:
            raise AssertionError("trilinear replacement identity failed")
        # Holder critical-mass bridge must invert exactly.
        lp_back = finite_energy_shell_lp_upper(cert.gaussian_critical_l2_mass_lower, cvol)
        if lp_back + 1e-12 < cert.gaussian_lp_mass_inside_shell_lower:
            raise AssertionError("critical mass / Lp bridge failed")
        worst_transfer_margin = min(worst_transfer_margin, cert.gaussian_transfer_lower)
        worst_mass = min(worst_mass, cert.gaussian_critical_l2_mass_lower)
    return {
        "samples": samples,
        "minimum_gaussian_transfer_lower": worst_transfer_margin,
        "minimum_critical_l2_mass_lower": worst_mass,
        "replacement_loss_at_1_percent": trilinear_replacement_loss(0.01),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=100_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-transfer-profile"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_narrow_shell_mass_certificate()
    result = stress(args.samples)
    payload = {"certificate": cert, "stress": result}
    (args.outdir / "transfer_profile_extraction.json").write_text(json.dumps(payload, indent=2))
    md = f"""# One-shot transfer-preserving Gaussian profile algebra\n\nStatus: **{cert['status']}** for the narrow-shell critical-mass bridge.\n\n- certified shell halfwidth: `{cert['shell_halfwidth']}`\n- certified 1% profile critical-L2 mass: `> {cert['critical_l2_mass_lower_bound']}`\n- shell-volume enclosure: `{cert['shell_volume_coefficient_ball']}`\n- critical-mass enclosure: `{cert['critical_l2_mass_ball']}`\n- random parameter checks: `{result['samples']}`\n- exact replacement loss at 1% profile distance: `{result['replacement_loss_at_1_percent']:.9f}`\n- minimum Gaussian-transfer lower bound in stress: `{result['minimum_gaussian_transfer_lower']:.9f}`\n- minimum shell-critical-L2 lower bound in stress: `{result['minimum_critical_l2_mass_lower']:.9e}`\n\nThe existence of the Gaussian approximation is supplied analytically by Christ's\nnear-extremizer theorem for Young convolution.  This workflow verifies only the\nnew deterministic transfer/mass consequences and does not pretend to numerically\ncertify Christ's theorem.\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
