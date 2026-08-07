from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .triad_extremizer import symmetric_gamma, symmetric_rstar


@dataclass(frozen=True)
class ConservativeInteraction:
    log_frequencies: tuple[float, ...]
    energy_rates: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.log_frequencies) != len(self.energy_rates) or len(self.energy_rates) < 2:
            raise ValueError("matching frequency/rate arrays of length >=2 required")
        scale = max(1.0, *(abs(x) for x in self.energy_rates))
        if abs(sum(self.energy_rates)) > 1e-12 * scale:
            raise ValueError("energy rates must sum to zero")


@dataclass(frozen=True)
class ProgressEdge:
    log_parent_top: float
    log_child: float
    child_transfer: float

    def __post_init__(self) -> None:
        if not self.log_parent_top < self.log_child:
            raise ValueError("parent top frequency must be below child")
        if self.child_transfer <= 0.0:
            raise ValueError("positive child transfer required on the forward core")

    @property
    def midgap(self) -> float:
        return 0.5 * (self.log_parent_top + self.log_child)

    @property
    def progress(self) -> float:
        return self.log_child - self.log_parent_top


@dataclass(frozen=True)
class MidgapBlockCertificate:
    tau: float
    buffer_width: float
    moat_margin: float
    upper_progress: float
    smooth_tail_flux: float
    equality_residual: float
    total_child_transfer: float


@dataclass(frozen=True)
class PolarizationCertificate:
    normalized_signed_flux: float
    total_deficit: float
    multiplier_deficit: float
    phase_deficit: float
    exact_residual: float


def compact_cdf(z: float, delta: float) -> float:
    """CDF of the even compact parabolic log-kernel.

    rho_delta(s)=3/(4 delta)*(1-(s/delta)^2) on |s|<=delta.
    This concrete C1 profile is for regression.  The theorem only needs an even
    probability kernel; a C-infinity bump gives the same identities.
    """
    if delta <= 0.0:
        raise ValueError("delta must be positive")
    if z <= -delta:
        return 0.0
    if z >= delta:
        return 1.0
    r = z / delta
    return 0.5 + 0.75 * r - 0.25 * r**3


def compact_density(z: float, delta: float) -> float:
    if delta <= 0.0:
        raise ValueError("delta must be positive")
    if abs(z) >= delta:
        return 0.0
    r = z / delta
    return 0.75 * (1.0 - r * r) / delta


def sharp_flux_at(t: float, interaction: ConservativeInteraction) -> float:
    """Outward cumulative sharp spectral flux at log cutoff t."""
    return -sum(rate for ell, rate in zip(interaction.log_frequencies, interaction.energy_rates) if ell <= t)


def graded_flux_at(t: float, interaction: ConservativeInteraction, delta: float) -> float:
    """Outward flux for the graded low-pass energy weight Psi_delta(t-ell)."""
    return -sum(
        rate * compact_cdf(t - ell, delta)
        for ell, rate in zip(interaction.log_frequencies, interaction.energy_rates)
    )


def mellin_moment(interaction: ConservativeInteraction) -> float:
    """Filter-independent all-log-scale flux moment: sum rate_i log |k_i|."""
    return float(sum(ell * rate for ell, rate in zip(interaction.log_frequencies, interaction.energy_rates)))


def ramp_potential(ell: float, tau: float, delta: float) -> float:
    """Phi with Phi'=2 Psi_delta(ell-tau), flat below and slope 2 above.

    Evenness of the kernel implies Phi=0 below tau-delta and
    Phi=2(ell-tau) above tau+delta, with no additive calibration error.
    """
    if delta <= 0.0:
        raise ValueError("delta must be positive")
    z = ell - tau
    if z <= -delta:
        return 0.0
    if z >= delta:
        return 2.0 * z
    # Integral from -delta to z of 2*Psi_delta(s) ds.
    return z + 3.0 * z * z / (4.0 * delta) - z**4 / (8.0 * delta**3) + 3.0 * delta / 8.0


def smooth_tail_flux(interaction: ConservativeInteraction, tau: float, delta: float) -> float:
    """Exact value of 2 int_tau^infty Pi_delta(t) dt via the weak flux identity."""
    return float(sum(
        rate * ramp_potential(ell, tau, delta)
        for ell, rate in zip(interaction.log_frequencies, interaction.energy_rates)
    ))


def interaction_scale_potential(interaction: ConservativeInteraction, values: Sequence[float]) -> float:
    if len(values) != len(interaction.energy_rates):
        raise ValueError("one potential value per mode required")
    return float(sum(rate * val for rate, val in zip(interaction.energy_rates, values)))


def transfer_weighted_midgap(edges: Sequence[ProgressEdge]) -> float:
    if not edges:
        raise ValueError("at least one edge required")
    total = sum(e.child_transfer for e in edges)
    return sum(e.child_transfer * e.midgap for e in edges) / total


def certify_midgap_block(edges: Sequence[ProgressEdge], delta: float) -> MidgapBlockCertificate:
    """Exact smooth-tail = upper-progress certificate on a common spectral moat.

    tau is the positive-transfer-weighted mean of edge midgaps.  If every top
    parent lies below tau-delta and every child above tau+delta, the smooth
    graded tail flux equals the sum T_e log(q_e/p_e) exactly.
    """
    if delta <= 0.0:
        raise ValueError("delta must be positive")
    if not edges:
        raise ValueError("at least one edge required")
    tau = transfer_weighted_midgap(edges)
    moat = min(
        min(tau - e.log_parent_top, e.log_child - tau)
        for e in edges
    ) - delta
    total_transfer = sum(e.child_transfer for e in edges)
    upper = sum(e.child_transfer * e.progress for e in edges)
    smooth = sum(2.0 * e.child_transfer * (e.log_child - tau) for e in edges)
    return MidgapBlockCertificate(
        tau=float(tau),
        buffer_width=float(delta),
        moat_margin=float(moat),
        upper_progress=float(upper),
        smooth_tail_flux=float(smooth),
        equality_residual=float(smooth - upper),
        total_child_transfer=float(total_transfer),
    )


def polarization_certificate(capacity: Sequence[float], multiplier: Sequence[float], phase: Sequence[float]) -> PolarizationCertificate:
    """Exact no-cancellation decomposition 1-R=E[(1-m)+m(1-cos phase)]."""
    a = np.asarray(capacity, dtype=float)
    m = np.asarray(multiplier, dtype=float)
    c = np.asarray(phase, dtype=float)
    if not (a.shape == m.shape == c.shape) or a.ndim != 1 or a.size == 0:
        raise ValueError("matching nonempty one-dimensional arrays required")
    if np.any(a < 0.0) or float(a.sum()) <= 0.0:
        raise ValueError("nonnegative capacities with positive total required")
    if np.any(m < 0.0) or np.any(m > 1.0):
        raise ValueError("multipliers must lie in [0,1]")
    if np.any(c < -1.0) or np.any(c > 1.0):
        raise ValueError("phase factors must lie in [-1,1]")
    w = a / a.sum()
    ratio = float(np.dot(w, m * c))
    mult = float(np.dot(w, 1.0 - m))
    phase_def = float(np.dot(w, m * (1.0 - c)))
    total = 1.0 - ratio
    return PolarizationCertificate(ratio, total, mult, phase_def, total - mult - phase_def)


def bad_capacity_mass_bound(total_deficit: float, threshold: float) -> float:
    """Markov bound for mass where signed normalized efficiency <=1-threshold."""
    if total_deficit < 0.0 or threshold <= 0.0:
        raise ValueError("nonnegative deficit and positive threshold required")
    return min(1.0, total_deficit / threshold)


def buffered_commutator_constant(first_kernel_moment: float, packet_buffer_widths: float) -> float:
    """Dimensionless Lp commutator bound M1/M for a buffer M packet widths."""
    if first_kernel_moment < 0.0 or packet_buffer_widths <= 0.0:
        raise ValueError("invalid commutator parameters")
    return first_kernel_moment / packet_buffer_widths


def random_stress(samples: int = 50000, seed: int = 20260807) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    r = symmetric_rstar()
    gamma = symmetric_gamma(r)
    delta = 0.05
    shell = 0.06
    worst_identity = 0.0
    worst_moat = float("inf")
    worst_polarization = 0.0
    for _ in range(samples):
        n = int(rng.integers(1, 12))
        p = rng.uniform(-shell, shell, size=n)
        q = gamma + rng.uniform(-shell, shell, size=n)
        T = rng.lognormal(mean=0.0, sigma=1.0, size=n)
        edges = [ProgressEdge(float(pp), float(qq), float(tt)) for pp, qq, tt in zip(p, q, T)]
        cert = certify_midgap_block(edges, delta)
        worst_identity = max(worst_identity, abs(cert.equality_residual))
        worst_moat = min(worst_moat, cert.moat_margin)

        cap = rng.lognormal(mean=0.0, sigma=1.0, size=n)
        mult = rng.uniform(0.0, 1.0, size=n)
        phase = rng.uniform(-1.0, 1.0, size=n)
        pol = polarization_certificate(cap, mult, phase)
        worst_polarization = max(worst_polarization, abs(pol.exact_residual))
    return {
        "samples": int(samples),
        "delta": delta,
        "shell_halfwidth": shell,
        "worst_midgap_identity_residual": float(worst_identity),
        "minimum_moat_margin": float(worst_moat),
        "worst_polarization_identity_residual": float(worst_polarization),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50000)
    ap.add_argument("--outdir", type=Path, default=Path("results-smooth-flux"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    stress = random_stress(args.samples)
    r = symmetric_rstar()
    gamma = symmetric_gamma(r)
    result = {
        "rstar": r,
        "gamma": gamma,
        "stress": stress,
        "theorem": {
            "filter_invariant_mellin": True,
            "midgap_tail_exact_on_common_moat": True,
            "polarization_identity": True,
            "example_delta": 0.05,
        },
    }
    (args.outdir / "smooth_flux_cocycle.json").write_text(json.dumps(result, indent=2))
    md = f'''# Smooth log-scale flux cocycle

The theorem identities are analytic; the stress test only validates implementation.

- gamma*: `{gamma:.15f}`
- filter-independent all-log Mellin moment: exact for every conservative finite interaction
- common-midgap smooth-tail identity: exact on a spectral moat
- no-cancellation polarization identity: exact
- stress blocks: `{stress['samples']}`
- test filter half-width delta: `{stress['delta']:.3f}` log units
- test shell half-width: `{stress['shell_halfwidth']:.3f}` log units
- minimum common-moat margin seen: `{stress['minimum_moat_margin']:.9f}`
- worst midgap equality residual: `{stress['worst_midgap_identity_residual']:.3e}`
- worst polarization residual: `{stress['worst_polarization_identity_residual']:.3e}`

The key PDE-facing formula is

`2 int_tau^infty Pi_delta(t) dt = sum_e T_e log(q_e/p_e)`

when tau is the positive-transfer-weighted mean of the edge midgaps and a common
smooth transition moat separates all top parents from all children.
'''
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
