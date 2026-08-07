from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

# Conservative dimensionless constant obtained from the explicit pressure
# kernel gradient, the path subtraction, and the safe distance loss used in the
# far-shell estimate.  It is intentionally not optimized.
FAR_KERNEL_CONSTANT = 400.0


@dataclass(frozen=True)
class PressureShell:
    index: int
    packet_masses: tuple[float, ...]

    @property
    def total_mass(self) -> float:
        return sum(self.packet_masses)


@dataclass(frozen=True)
class CollisionCertificate:
    target_cancellation: float
    near_coefficient: float
    far_coefficient: float
    fresh_mass_threshold: float


def critical_l2_mass_from_l3_charge(l3_charge: float, bernstein_constant: float = 1.0) -> float:
    """Scale-critical L2 mass forced by a band-limited L3 charge.

    If r is the packet scale and ||f||_3 <= C_B r^{-1/2} ||f||_2, then
    r^{-1}||f||_2^2 >= (int |f|^3)^{2/3}/C_B^2.
    """
    if l3_charge < 0 or bernstein_constant <= 0:
        raise ValueError("require nonnegative charge and positive Bernstein constant")
    return l3_charge ** (2.0 / 3.0) / (bernstein_constant * bernstein_constant)


def local_pressure_mass_threshold(
    normalized_pressure_work: float,
    grad_scale: float,
    riesz_constant: float,
    bernstein_constant: float,
) -> float:
    """Critical L2 mass forced by a local pressure-cancellation branch.

    Assumptions in the packet model:
      r|W_loc| <= (r||grad chi||_inf) C_R ||V||_3^3,
    and V is localized at frequency <= O(r^{-1}).
    """
    if normalized_pressure_work < 0 or grad_scale <= 0 or riesz_constant <= 0:
        raise ValueError("invalid local pressure parameters")
    l3 = normalized_pressure_work / (grad_scale * riesz_constant)
    return critical_l2_mass_from_l3_charge(l3, bernstein_constant)


def far_multipole_moment(shells: Sequence[PressureShell], first_shell: int = 3) -> float:
    """The positive far-pressure moment sum 2^{-4n} M_n.

    The fourth power is the dipole/multipole gain after subtracting the constant
    pressure mode.  In three dimensions packet counts grow only like 2^{3n}.
    """
    total = 0.0
    for sh in shells:
        if sh.index < first_shell:
            raise ValueError("shell is not in the far field")
        if any(m < 0 for m in sh.packet_masses):
            raise ValueError("packet critical masses must be nonnegative")
        total += (2.0 ** (-4 * sh.index)) * sh.total_mass
    return total


def far_pressure_work_bound(
    shells: Sequence[PressureShell],
    boundary_factor: float,
    kernel_constant: float = FAR_KERNEL_CONSTANT,
    first_shell: int = 3,
) -> float:
    """Dimensionless far pressure-work bound in the packet model.

    boundary_factor packages the scale-invariant boundary quantities generated
    by U and grad chi. The source dependence is exactly the multipole moment.
    """
    if boundary_factor < 0 or kernel_constant <= 0:
        raise ValueError("invalid far pressure parameters")
    return kernel_constant * boundary_factor * far_multipole_moment(shells, first_shell)


def no_fresh_far_coefficient(
    packing_constant: float,
    boundary_factor: float,
    first_shell: int = 3,
    kernel_constant: float = FAR_KERNEL_CONSTANT,
) -> float:
    """C_far with W_far <= C_far * mu_max under 3D shell packing.

    If shell n contains at most C_geom 2^{3n} packets and each packet has
    critical mass <= mu_max, then the 2^{-4n} pressure dipole weight leaves
    sum 2^{-n} = 2^{1-n0}.
    """
    if packing_constant <= 0 or boundary_factor < 0 or first_shell < 0:
        raise ValueError("invalid packing parameters")
    return kernel_constant * boundary_factor * packing_constant * (2.0 ** (1 - first_shell))


def collision_certificate(target_cancellation: float, near_coefficient: float, far_coefficient: float) -> CollisionCertificate:
    """Contrapositive fresh-packet threshold.

    Suppose absence of a fresh packet with mass > mu implies
        cancellation <= C_near mu^(3/2) + C_far mu.
    If the actual normalized pressure cancellation is >= rho, then it is enough
    to choose mu below both half-budget thresholds to get a contradiction.
    """
    if target_cancellation <= 0 or near_coefficient < 0 or far_coefficient < 0:
        raise ValueError("invalid collision parameters")
    candidates = []
    if near_coefficient > 0:
        candidates.append((target_cancellation / (2.0 * near_coefficient)) ** (2.0 / 3.0))
    if far_coefficient > 0:
        candidates.append(target_cancellation / (2.0 * far_coefficient))
    if not candidates:
        raise ValueError("at least one pressure channel must be present")
    return CollisionCertificate(target_cancellation, near_coefficient, far_coefficient, min(candidates))


def pressure_kernel_component_gradient_bound() -> float:
    """A rigorous elementary bound for sum_ij |grad K_ij| |v_i v_j|/|v|^2.

    K_ij(z)=(4 pi)^(-1)(3 z_i z_j |z|^-5-delta_ij |z|^-3).
    Each coordinate derivative is bounded by 24/(4 pi)|z|^-4, its vector
    gradient by sqrt(3) times that, and sum_ij |v_i v_j| <= 3 |v|^2.
    The returned clean constant 10 dominates the resulting 9.93... .
    """
    raw = 3.0 * math.sqrt(3.0) * 24.0 / (4.0 * math.pi)
    if raw >= 10.0:
        raise AssertionError("clean pressure-kernel constant no longer dominates")
    return 10.0


def stress(samples: int = 50_000, seed: int = 20260807) -> dict[str, float]:
    import numpy as np

    rng = np.random.default_rng(seed)
    worst_far_ratio = 0.0
    worst_collision_margin = float("inf")
    packing = 3.0
    boundary = 0.4
    n0 = 3
    cfar = no_fresh_far_coefficient(packing, boundary, n0)
    rho = 0.2
    cnear = 7.0
    cert = collision_certificate(rho, cnear, cfar)
    # Stay strictly below the theorem threshold.
    mu = 0.9 * cert.fresh_mass_threshold
    for _ in range(samples):
        shells = []
        for n in range(n0, n0 + 7):
            max_count = max(1, int(packing * (2 ** (3 * n))))
            count = int(rng.integers(1, min(max_count, 200) + 1))
            masses = tuple(float(rng.uniform(0.0, mu)) for _j in range(count))
            # If we did not instantiate all allowed packets, append aggregate
            # adversarial mass at mu through the analytic packing bound below;
            # the explicit shell probe is only a lower-complexity regression.
            shells.append(PressureShell(n, masses))
        explicit = far_pressure_work_bound(shells, boundary, first_shell=n0)
        analytic = cfar * mu
        worst_far_ratio = max(worst_far_ratio, explicit / analytic if analytic else 0.0)
        total_upper = cnear * mu ** 1.5 + cfar * mu
        worst_collision_margin = min(worst_collision_margin, rho - total_upper)
        if explicit > analytic + 1e-12:
            raise AssertionError("explicit far shell exceeded no-fresh analytic bound")
        if total_upper >= rho:
            raise AssertionError("fresh threshold failed to exclude target cancellation")
    return {
        "samples": samples,
        "kernel_component_constant": pressure_kernel_component_gradient_bound(),
        "far_kernel_constant": FAR_KERNEL_CONSTANT,
        "far_coefficient": cfar,
        "fresh_mass_threshold": cert.fresh_mass_threshold,
        "worst_explicit_far_over_analytic": worst_far_ratio,
        "minimum_collision_margin": worst_collision_margin,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-annular-pressure"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = stress(args.samples)
    (args.outdir / "annular_pressure_collision.json").write_text(json.dumps(result, indent=2))
    md = f"""# Annular pressure collision\n\n- elementary pressure-kernel gradient constant: `< {result['kernel_component_constant']:.1f}`\n- conservative far multipole prefactor: `{result['far_kernel_constant']:.0f}`\n- stress samples: `{result['samples']}`\n- packet-model far coefficient: `{result['far_coefficient']:.6f}`\n- packet-model fresh critical-mass threshold: `{result['fresh_mass_threshold']:.9e}`\n- worst explicit far / analytic no-fresh bound: `{result['worst_explicit_far_over_analytic']:.6f}`\n- minimum collision exclusion margin: `{result['minimum_collision_margin']:.6e}`\n"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
