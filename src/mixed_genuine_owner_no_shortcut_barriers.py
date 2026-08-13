from __future__ import annotations

import math
from dataclasses import dataclass

from src.critical_shell_service_reentry import theorem_certificate as critical_shell_certificate
from src.fresh_service_scale_reentry import theorem_certificate as fresh_scale_certificate
from src.physical_multicurrency_master import theorem_certificate as multicurrency_certificate
from src.radial_spectral_crossing_layer_cake import theorem_certificate as radial_certificate


STATUS = (
    "DRAFT_MIXED_GENUINE_OWNER_NO_SHORTCUT_BARRIERS__"
    "CRITICAL_NE_NOT_FINITE_RESET__RADIAL_TRAFFIC_NOT_BUDGET__"
    "GENERIC_SHELL_REGISTRATION_NOT_PROGRESS__FRESH_SCALE_HAS_NO_DIRECTION"
)


@dataclass(frozen=True)
class CriticalEnergyGeometricChain:
    base_frequency: float
    critical_mass: float
    ratio: float
    depth: int
    shell_frequencies: tuple[float, ...]
    shell_energies: tuple[float, ...]
    total_shell_energy: float
    infinite_chain_energy_upper: float

    def __post_init__(self) -> None:
        if not (self.base_frequency > 0 and self.critical_mass > 0 and self.ratio > 1):
            raise ValueError("positive base/critical mass and geometric UV ratio >1 required")
        if self.depth <= 0:
            raise ValueError("positive chain depth required")
        if len(self.shell_frequencies) != self.depth or len(self.shell_energies) != self.depth:
            raise ValueError("chain coordinates must match depth")
        if self.total_shell_energy > self.infinite_chain_energy_upper + 1e-13 * max(1.0, self.infinite_chain_energy_upper):
            raise AssertionError("finite prefix exceeded geometric infinite-chain energy upper")


def critical_energy_geometric_chain(
    base_frequency: float,
    critical_mass: float,
    depth: int,
    ratio: float = 2.0,
) -> CriticalEnergyGeometricChain:
    """Anti-reset scaling model: NE=mu at every geometrically increasing shell.

    This is deliberately **not** an existence claim for a Navier--Stokes
    trajectory. It proves only that the global L2 energy budget, by itself,
    cannot count scale-critical shell events.  If N_j=N_0 r^j and N_j E_j=mu,
    then E_j=mu/N_j and the entire infinite geometric allocation has finite
    total energy mu/N_0 * r/(r-1).
    """
    N0 = float(base_frequency)
    mu = float(critical_mass)
    r = float(ratio)
    L = int(depth)
    if not (math.isfinite(N0) and math.isfinite(mu) and math.isfinite(r)):
        raise ValueError("finite chain parameters required")
    if N0 <= 0 or mu <= 0 or r <= 1 or L <= 0:
        raise ValueError("positive base/critical mass, r>1 and depth>0 required")
    freqs = tuple(N0 * (r**j) for j in range(L))
    energies = tuple(mu / N for N in freqs)
    total = math.fsum(energies)
    upper = (mu / N0) * r / (r - 1.0)
    return CriticalEnergyGeometricChain(
        base_frequency=N0,
        critical_mass=mu,
        ratio=r,
        depth=L,
        shell_frequencies=freqs,
        shell_energies=energies,
        total_shell_energy=total,
        infinite_chain_energy_upper=upper,
    )


def certified_no_shortcut_barriers() -> dict[str, object]:
    multi = multicurrency_certificate()
    radial = radial_certificate()
    shell = critical_shell_certificate()
    fresh = fresh_scale_certificate()

    physics = str(multi.get("physics", ""))
    if "critical NE" not in physics or "D_V" not in physics:
        raise AssertionError("multicurrency anti-reset theorem no longer exposes the critical NE/D_V barrier")
    if radial.get("gross_crossing_budget") is not False:
        raise AssertionError("radial crossing theorem unexpectedly declared gross traffic a finite budget")
    if "not fabricated" not in str(shell.get("scale_scope", "")):
        raise AssertionError("generic shell theorem no longer states its no-progress scope")
    fresh_scope = str(fresh.get("scale_scope", ""))
    if "2N" not in fresh_scope or "not asserted" not in fresh_scope:
        raise AssertionError("fresh service theorem no longer exposes its direction-free scale scope")

    return {
        "status": STATUS,
        "critical_reset_forbidden": True,
        "critical_reset_reason": physics,
        "gross_radial_budget_forbidden": True,
        "gross_radial_budget_source": radial.get("crossing_semantics"),
        "generic_shell_progress_forbidden": True,
        "generic_shell_scale_scope": shell.get("scale_scope"),
        "fresh_scale_direction_forbidden": True,
        "fresh_scale_scope": fresh.get("scale_scope"),
        "consequence": "a mixed-owner closure must use a new native cross-owner law or typed ledger; it cannot be obtained by scalarizing critical mass, normalized dissipation, radial traffic, generic shell registration, or fresh-scale selection",
        "claims_global_regularity": False,
    }


def theorem_certificate() -> dict[str, object]:
    barriers = certified_no_shortcut_barriers()
    chain = critical_energy_geometric_chain(1.0, 1.0, 64, 2.0)
    return {
        **barriers,
        "critical_energy_scaling_counterexample": "for N_j=N0 r^j and N_j E_j=mu, sum_j E_j=(mu/N0) r/(r-1)<infinity; this is a budget-scaling anti-theorem, not an NS trajectory-existence statement",
        "sample_geometric_chain_total_energy": chain.total_shell_energy,
        "sample_geometric_chain_infinite_upper": chain.infinite_chain_energy_upper,
        "scope": "these are proof-obstruction guards for the draft native normal form; they do not prove or disprove existence of an infinite NS event path and make no global-regularity claim",
    }
