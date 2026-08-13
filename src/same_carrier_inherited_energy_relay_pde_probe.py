from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from src.helical_mode_set_energy_continuity_pde_probe import (
    HelicalModeContinuityPDEProbe,
    run_probe as run_mode_continuity_probe,
)

STATUS = (
    "ACTUAL_FOURIER_GALERKIN_NS_SAME_PHYSICAL_HELICAL_CARRIER_STOCK__"
    "INHERITED_ENDPOINT_ENERGY_IS_STOCK_NOT_NEW_WORK__"
    "GROSS_CANONICAL_DW_PLUS_MINUS_AND_VISCOSITY_CONTINUITY"
)


@dataclass(frozen=True)
class ActualNSInheritedStockObservation:
    resolution: int
    child_mode: tuple[int, int, int]
    helicity: int
    initial_energy: float
    final_energy: float
    inherited_fraction: float
    inheritance_gate_margin: float
    integrated_positive_work: float
    integrated_negative_work: float
    viscous_dissipation: float
    interval_continuity_native_residual: float
    global_energy_balance_relative_residual: float

    def __post_init__(self) -> None:
        vals = (
            self.initial_energy,
            self.final_energy,
            self.inherited_fraction,
            self.integrated_positive_work,
            self.integrated_negative_work,
            self.viscous_dissipation,
            self.interval_continuity_native_residual,
            self.global_energy_balance_relative_residual,
        )
        if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in vals):
            raise ValueError("finite nonnegative actual NS stock/work data required")
        if self.initial_energy <= 0.0 or self.final_energy <= 0.0:
            raise AssertionError("actual NS selected physical carrier lost positive endpoint stock")
        if abs(self.inherited_fraction - self.initial_energy / self.final_energy) > 2e-13 * max(1.0, self.inherited_fraction):
            raise AssertionError("actual NS inherited fraction changed from endpoint physical stock")
        if self.inheritance_gate_margin < -2e-12 * max(1.0, self.initial_energy, self.final_energy):
            raise AssertionError("actual NS stock failed the E0>=E1/5 inheritance face")


@dataclass(frozen=True)
class SameCarrierInheritedEnergyPDEProbe:
    status: str
    upstream_mode_continuity: HelicalModeContinuityPDEProbe
    observations: tuple[ActualNSInheritedStockObservation, ...]
    minimum_inheritance_gate_margin: float
    maximum_mode_continuity_native_residual: float
    maximum_global_energy_balance_relative_residual: float

    def __post_init__(self) -> None:
        if self.status != STATUS or not self.observations:
            raise ValueError("actual NS same-carrier stock referee provenance mismatch")
        if self.minimum_inheritance_gate_margin < -2e-12:
            raise AssertionError("actual NS inheritance face was not realized")
        if self.maximum_mode_continuity_native_residual > 5e-5:
            raise AssertionError("actual NS physical carrier continuity left certified tolerance")
        if self.maximum_global_energy_balance_relative_residual > 5e-5:
            raise AssertionError("actual NS global energy balance left certified tolerance")


def run_probe(
    *,
    resolutions: Sequence[int] = (24, 28),
    cutoff: int = 7,
    steps: int = 32,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    phase_sign: int = 1,
) -> SameCarrierInheritedEnergyPDEProbe:
    upstream = run_mode_continuity_probe(
        resolutions=tuple(int(v) for v in resolutions),
        cutoff=int(cutoff),
        steps=int(steps),
        viscosity=float(viscosity),
        amplitude=float(amplitude),
        duration=float(duration),
        phase_sign=int(phase_sign),
    )
    obs = tuple(
        ActualNSInheritedStockObservation(
            resolution=r.resolution,
            child_mode=r.child_mode,
            helicity=r.helicity,
            initial_energy=r.initial_energy,
            final_energy=r.final_energy,
            inherited_fraction=r.initial_energy / r.final_energy,
            inheritance_gate_margin=r.initial_energy - 0.2 * r.final_energy,
            integrated_positive_work=r.integrated_positive_work,
            integrated_negative_work=r.integrated_negative_work,
            viscous_dissipation=r.viscous_dissipation,
            interval_continuity_native_residual=r.interval_continuity_native_residual,
            global_energy_balance_relative_residual=r.global_energy_balance_relative_residual,
        )
        for r in upstream.runs
    )
    return SameCarrierInheritedEnergyPDEProbe(
        status=STATUS,
        upstream_mode_continuity=upstream,
        observations=obs,
        minimum_inheritance_gate_margin=min(o.inheritance_gate_margin for o in obs),
        maximum_mode_continuity_native_residual=max(o.interval_continuity_native_residual for o in obs),
        maximum_global_energy_balance_relative_residual=max(o.global_energy_balance_relative_residual for o in obs),
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=STATUS)
    ap.add_argument("--resolutions", type=int, nargs="+", default=(24, 28))
    ap.add_argument("--cutoff", type=int, default=7)
    ap.add_argument("--steps", type=int, default=32)
    ap.add_argument("--viscosity", type=float, default=0.03)
    ap.add_argument("--amplitude", type=float, default=1.0)
    ap.add_argument("--duration", type=float, default=0.001)
    ap.add_argument("--phase-sign", type=int, default=1)
    ap.add_argument("--outdir", type=Path, default=Path("results-same-carrier-inherited-energy-relay-ns"))
    args = ap.parse_args()
    out = run_probe(
        resolutions=args.resolutions,
        cutoff=args.cutoff,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        phase_sign=args.phase_sign,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "same_carrier_inherited_energy_relay_pde.json").write_text(
        json.dumps(asdict(out), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "# Actual Fourier--Galerkin Navier--Stokes referee: inherited carrier stock",
        "",
        f"Status: **{STATUS}**.",
        "",
        "The referee follows one fixed physical Fourier--helical mode through an actual dealiased Galerkin Navier--Stokes trajectory.  Its endpoint energy is persistent stock on that physical carrier; gross canonical dW+ / dW- and viscosity satisfy the already-certified between-time continuity law.  This finite-PDE referee checks the stock ontology only.  The smooth same-Q/same-psi no-hit condition remains the exact typed theorem and is not replaced by a numerical surrogate.",
        "",
        f"- minimum E0-E1/5 inheritance margin: `{out.minimum_inheritance_gate_margin:.12g}`",
        f"- maximum modal continuity native residual: `{out.maximum_mode_continuity_native_residual:.3e}`",
        f"- maximum global energy-balance residual: `{out.maximum_global_energy_balance_relative_residual:.3e}`",
    ]
    for o in out.observations:
        lines.extend([
            "",
            f"## resolution {o.resolution}",
            f"- physical carrier: mode `{o.child_mode}`, helicity `{o.helicity}`",
            f"- initial/final stock: `{o.initial_energy:.12g}` / `{o.final_energy:.12g}`",
            f"- inherited stock fraction E0/E1: `{o.inherited_fraction:.12g}`",
            f"- gross positive/negative nonlinear work: `{o.integrated_positive_work:.12g}` / `{o.integrated_negative_work:.12g}`",
            f"- viscous dissipation: `{o.viscous_dissipation:.12g}`",
        ])
    lines.extend([
        "",
        "No FIFO/LIFO pairing of earlier deposits to later withdrawals is introduced, and no hard interaction cell is treated as a between-time wallet.  No global-regularity claim is made.",
    ])
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
