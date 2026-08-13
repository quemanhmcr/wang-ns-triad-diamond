from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

from src.continuum_helical_edge_measure_registration import unitary_fourier_convolution_factor
from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel
from src.hard_tail_true_upward_supply import hard_tail_upward_supply_split
from src.hard_tail_true_upward_supply_pde_probe import (
    HardTailTrueUpwardSupplyPDEProbe,
    _selected_closed_triad,
    _spectral_geometry,
    adversarial_mixed_fate_initial_state,
    run_probe as run_upstream_probe,
)
from src.pure_uv_true_upward_natural_window import (
    PURE_UV_PARENT_UPPER_RATIO,
    PURE_UV_SHELL_INDEX,
    STATUS as THEOREM_STATUS,
    pure_uv_first_shell_law,
)

STATUS = (
    "ORTHOGONAL_FOURIER_GALERKIN_NS_PURE_UV_TRUE_UPWARD_FIRST_SHELL__"
    "ACTUAL_CANONICAL_RECIPIENT_CAUSE__M_EQUALS_2N__H_EQUALS_U_BY_STRICT_SUPPORT__"
    "NO_OUTPUT_SCALE_SELECTION"
)


@dataclass(frozen=True)
class ActualPureUVLawObservation:
    resolution: int
    cutoff: int
    radial_boundary: float
    pure_uv_physical_work: float
    pure_uv_common_work: float
    recipient_charge_count: int
    donor_sidecar_count: int
    recipient_shell_index: int
    recipient_shell_scale: float
    p_scale: float
    minimum_parent_above_quarter_shell_margin: float
    maximum_parent_to_shell_ratio: float

    def __post_init__(self) -> None:
        if self.resolution < 16 or self.cutoff <= 0 or self.radial_boundary <= 0.0:
            raise ValueError("valid Fourier-Galerkin pure-UV observation geometry required")
        if self.pure_uv_physical_work <= 0.0 or self.pure_uv_common_work <= 0.0:
            raise ValueError("actual pure-UV observation lost positive canonical work")
        if self.recipient_charge_count <= 0 or self.donor_sidecar_count <= 0:
            raise ValueError("actual pure-UV observation lost recipient/donor provenance")
        if self.recipient_shell_index != PURE_UV_SHELL_INDEX:
            raise AssertionError("actual NS pure-UV observation left the first shell")
        if abs(self.recipient_shell_scale - 2.0 * self.radial_boundary) > 5.0e-12 * self.recipient_shell_scale:
            raise AssertionError("actual NS pure-UV shell is not M=2N")
        if self.p_scale != 1.0:
            raise AssertionError("actual NS pure-UV law acquired an artificial output-scale selection")
        if self.minimum_parent_above_quarter_shell_margin <= 0.0:
            raise AssertionError("actual NS pure-UV parent touched resolved cutoff support")
        if self.maximum_parent_to_shell_ratio > PURE_UV_PARENT_UPPER_RATIO + 5.0e-12:
            raise AssertionError("actual NS pure-UV parent escaped the 3M/2 corridor")


@dataclass(frozen=True)
class PureUVNaturalWindowPDEProbe:
    status: str
    theorem_status: str
    upstream_probe: HardTailTrueUpwardSupplyPDEProbe
    observations: tuple[ActualPureUVLawObservation, ...]
    maximum_common_work_representation_relative_spread: float
    minimum_cutoff_support_margin: float
    maximum_scale_probability_residual: float
    coexistence_with_resolved_contact_observed: bool

    def __post_init__(self) -> None:
        if self.status != STATUS or self.theorem_status != THEOREM_STATUS:
            raise ValueError("pure-UV PDE probe status provenance mismatch")
        if not self.observations:
            raise ValueError("nonempty actual pure-UV observations required")
        if self.maximum_common_work_representation_relative_spread > 5.0e-8:
            raise AssertionError("same finite-cutoff pure-UV common work changed across FFT representations")
        if self.minimum_cutoff_support_margin <= 0.0:
            raise AssertionError("pure-UV parent support no longer forces h=u")
        if self.maximum_scale_probability_residual > 5.0e-15:
            raise AssertionError("pure-UV actual PDE law acquired p_scale != 1")
        if not self.coexistence_with_resolved_contact_observed:
            raise AssertionError("actual NS fixture lost simultaneous pure-UV/contact upward phenomena")


def _actual_law_observation(*, resolution: int, cutoff: int, amplitude: float, radial_boundary: float) -> ActualPureUVLawObservation:
    n = int(resolution)
    k, k2, dealias, actual_cutoff = _spectral_geometry(n, int(cutoff))
    if actual_cutoff != int(cutoff):
        raise AssertionError("pure-UV law observation changed requested Galerkin cutoff")
    state = adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=float(amplitude))
    triad = _selected_closed_triad(state)
    kernel = cyclic_triad_measure_kernel(
        triad, quotient_measure_mass=1.0 / unitary_fourier_convolution_factor()
    )
    split = hard_tail_upward_supply_split(triad, kernel, boundary=float(radial_boundary))
    law = pure_uv_first_shell_law(split)
    pure = tuple(a for a in split.atoms if a.pure_uv_hh_by_support)
    M = law.recipient_shell_scale
    margin = min(min(a.interaction_parent_radii) - 0.25 * M for a in pure)
    return ActualPureUVLawObservation(
        resolution=n,
        cutoff=int(cutoff),
        radial_boundary=float(radial_boundary),
        pure_uv_physical_work=law.total_canonical_positive_mass,
        pure_uv_common_work=law.total_common_unit_work,
        recipient_charge_count=len(law.recipient_submeasures),
        donor_sidecar_count=sum(len(c.donor_closed_mode_indices) for c in law.recipient_submeasures),
        recipient_shell_index=PURE_UV_SHELL_INDEX,
        recipient_shell_scale=M,
        p_scale=law.p_scale,
        minimum_parent_above_quarter_shell_margin=margin,
        maximum_parent_to_shell_ratio=max(a.comparable_parent_upper_ratio for a in pure),
    )


def _relative_spread(values: Sequence[float]) -> float:
    vals = tuple(float(v) for v in values)
    if not vals:
        return 0.0
    scale = max(max(abs(v) for v in vals), 1.0e-300)
    return (max(vals) - min(vals)) / scale


def run_probe(
    *,
    main_resolutions: Sequence[int] = (24, 28),
    deep_resolutions: Sequence[int] = (20, 24),
    main_cutoff: int = 7,
    deep_cutoff: int = 2,
    main_steps: int = 48,
    deep_steps: int = 16,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    main_duration: float = 0.001,
    deep_duration: float = 0.0002,
    closed_tail_resolutions: Sequence[int] = (20, 24),
    closed_tail_steps: int = 24,
    closed_tail_duration: float = 0.0003,
) -> PureUVNaturalWindowPDEProbe:
    upstream = run_upstream_probe(
        main_resolutions=tuple(int(v) for v in main_resolutions),
        deep_resolutions=tuple(int(v) for v in deep_resolutions),
        main_cutoff=int(main_cutoff),
        deep_cutoff=int(deep_cutoff),
        main_steps=int(main_steps),
        deep_steps=int(deep_steps),
        viscosity=float(viscosity),
        amplitude=float(amplitude),
        main_duration=float(main_duration),
        deep_duration=float(deep_duration),
        closed_tail_resolutions=tuple(int(v) for v in closed_tail_resolutions),
        closed_tail_steps=int(closed_tail_steps),
        closed_tail_duration=float(closed_tail_duration),
    )
    obs = tuple(
        _actual_law_observation(
            resolution=int(n), cutoff=int(main_cutoff), amplitude=float(amplitude), radial_boundary=8.0
        )
        for n in main_resolutions
    )
    spread = _relative_spread([o.pure_uv_common_work for o in obs])
    min_margin = min(o.minimum_parent_above_quarter_shell_margin for o in obs)
    p_residual = max(abs(o.p_scale - 1.0) for o in obs)
    coexist = all(o.pure_uv_work > 0.0 and o.resolved_contact_work > 0.0 for o in upstream.selected_pure_support)
    return PureUVNaturalWindowPDEProbe(
        status=STATUS,
        theorem_status=THEOREM_STATUS,
        upstream_probe=upstream,
        observations=obs,
        maximum_common_work_representation_relative_spread=spread,
        minimum_cutoff_support_margin=min_margin,
        maximum_scale_probability_residual=p_residual,
        coexistence_with_resolved_contact_observed=coexist,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--main-resolutions", type=int, nargs="+", default=(24, 28))
    ap.add_argument("--deep-resolutions", type=int, nargs="+", default=(20, 24))
    ap.add_argument("--main-cutoff", type=int, default=7)
    ap.add_argument("--deep-cutoff", type=int, default=2)
    ap.add_argument("--main-steps", type=int, default=48)
    ap.add_argument("--deep-steps", type=int, default=16)
    ap.add_argument("--viscosity", type=float, default=0.03)
    ap.add_argument("--amplitude", type=float, default=1.0)
    ap.add_argument("--main-duration", type=float, default=0.001)
    ap.add_argument("--deep-duration", type=float, default=0.0002)
    ap.add_argument("--closed-tail-resolutions", type=int, nargs="+", default=(20, 24))
    ap.add_argument("--closed-tail-steps", type=int, default=24)
    ap.add_argument("--closed-tail-duration", type=float, default=0.0003)
    ap.add_argument("--outdir", type=Path, default=Path("results-pure-uv-natural-window-ns"))
    args = ap.parse_args()
    out = run_probe(
        main_resolutions=args.main_resolutions,
        deep_resolutions=args.deep_resolutions,
        main_cutoff=args.main_cutoff,
        deep_cutoff=args.deep_cutoff,
        main_steps=args.main_steps,
        deep_steps=args.deep_steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        main_duration=args.main_duration,
        deep_duration=args.deep_duration,
        closed_tail_resolutions=args.closed_tail_resolutions,
        closed_tail_steps=args.closed_tail_steps,
        closed_tail_duration=args.closed_tail_duration,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "pure_uv_true_upward_natural_window_pde.json").write_text(
        json.dumps(asdict(out), indent=2, sort_keys=True) + "\n"
    )
    lines = [
        "# Actual Fourier--Galerkin Navier--Stokes referee: pure-UV true-upward natural-window binding",
        "",
        f"Status: **{STATUS}**.",
        "",
        "The referee reuses the certified full Galerkin hard-tail probe rather than copying its solver.  The same actual NS state simultaneously carries pure-UV and resolved-contact upward submeasures; only the pure branch is read here.  Its same-time donor sidecars are coalesced to recipient canonical submeasures, and the strict parent support is checked directly against M/4.",
        "",
        f"- pure common-work cross-FFT relative spread: `{out.maximum_common_work_representation_relative_spread:.3e}`",
        f"- minimum physical parent margin above M/4: `{out.minimum_cutoff_support_margin:.12g}`",
        f"- maximum p_scale residual: `{out.maximum_scale_probability_residual:.3e}`",
        f"- pure/contact coexistence observed: `{out.coexistence_with_resolved_contact_observed}`",
        f"- upstream pure-support representation residual: `{out.upstream_probe.maximum_pure_support_work_representation_native_residual:.3e}`",
    ]
    for o in out.observations:
        lines.extend([
            "",
            f"## resolution {o.resolution}",
            f"- pure physical/common work: `{o.pure_uv_physical_work:.12g}` / `{o.pure_uv_common_work:.12g}`",
            f"- recipient charges / donor sidecars: `{o.recipient_charge_count}` / `{o.donor_sidecar_count}`",
            f"- shell index/scale: `{o.recipient_shell_index}` / `{o.recipient_shell_scale:.12g}`",
            f"- parent upper ratio: `{o.maximum_parent_to_shell_ratio:.12g}`",
        ])
    lines.extend([
        "",
        "This is a finite-PDE referee for the support/cause identities.  The continuum total-variation/Young capacity remains the analytic theorem and is not replaced by finite-grid numerics.  No global-regularity claim is made.",
    ])
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(asdict(out), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
