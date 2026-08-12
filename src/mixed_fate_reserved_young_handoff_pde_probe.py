from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.canonical_positive_edge_work_routing import (
    route_canonical_positive_edge_work,
    single_hard_role_map,
)
from src.continuum_helical_edge_measure_pde_probe import (
    _index,
    _leray_dealias,
    _nonlinear_term,
    _rk4_step,
    _snapshot_with_ledger,
    _spectral_geometry,
)
from src.helical import coupling_g, helical_basis
from src.mixed_fate_reserved_young_handoff import inherited_negative_work

STATUS = (
    "EVOLVED_DEALIASED_FOURIER_GALERKIN_NS_MIXED_FATE_HARD_CELL_AUDIT__"
    "ACTUAL_GOOD_BAD_NEGATIVE_EDGE_WORK__FULL_SIGNED_CELL_RESERVATION_IDENTITIES"
)

CHILD = (7, 6, 5)
GOOD_PARENTS = ((5, 0, 4), (2, 6, 1))
BAD_PARENTS = ((7, 7, 7), (0, -1, -2))
NEGATIVE_PARENTS = ((6, 7, 7), (1, -1, -2))
HELICITIES = (1, -1, 1)


def _unit_phase_for_child_work(
    x: tuple[int, int, int],
    y: tuple[int, int, int],
    *,
    desired_phase_alignment: float,
) -> complex:
    sx, sy, sz = HELICITIES
    xv = np.asarray(x, dtype=float)
    yv = np.asarray(y, dtype=float)
    zv = np.asarray(CHILD, dtype=float)
    g = coupling_g(xv, yv, -zv, sx, sy, sz)
    if abs(g) == 0.0:
        raise AssertionError("adversarial NS mixed-fate pair has zero Waleffe coupling")
    signed_frequency = sx * float(np.linalg.norm(xv)) - sy * float(np.linalg.norm(yv))
    sign_a = 1.0 if signed_frequency >= 0.0 else -1.0
    c = float(desired_phase_alignment)
    if c not in (-1.0, 1.0):
        raise ValueError("probe uses exact aligned or anti-aligned child work")
    # Child and y-parent amplitudes have unit phase.  The helical registration
    # defines c=sign(a) Re(conj(az) conj(g) ax ay); this phase makes c exact.
    return sign_a * c * g / abs(g)


def _set_real_helical_mode(
    coeff: np.ndarray,
    wavevector: tuple[int, int, int],
    helicity: int,
    amplitude: complex,
) -> None:
    n = int(coeff.shape[1])
    k = np.asarray(wavevector, dtype=float)
    value = complex(amplitude) * helical_basis(k, int(helicity))
    idx = _index(wavevector, n)
    neg = tuple(-v for v in wavevector)
    nidx = _index(neg, n)
    if np.linalg.norm(coeff[(slice(None),) + idx]) > 1.0e-14:
        raise AssertionError("adversarial Galerkin mode was assigned twice")
    if np.linalg.norm(coeff[(slice(None),) + nidx]) > 1.0e-14:
        raise AssertionError("adversarial conjugate Galerkin mode was assigned twice")
    coeff[(slice(None),) + idx] = value
    coeff[(slice(None),) + nidx] = np.conjugate(value)


def adversarial_mixed_fate_initial_state(
    resolution: int,
    k: np.ndarray,
    k2: np.ndarray,
    dealias: np.ndarray,
    *,
    amplitude: float = 1.0,
) -> np.ndarray:
    """Sparse real divergence-free Fourier polynomial with three fates at one child."""
    n = int(resolution)
    amp = float(amplitude)
    if not math.isfinite(amp) or amp <= 0.0:
        raise ValueError("positive finite adversarial NS amplitude required")
    coeff = np.zeros((3, n, n, n), dtype=complex)
    sx, sy, sz = HELICITIES

    _set_real_helical_mode(coeff, CHILD, sz, 1.0)

    gx, gy = GOOD_PARENTS
    _set_real_helical_mode(coeff, gx, sx, _unit_phase_for_child_work(gx, gy, desired_phase_alignment=1.0))
    _set_real_helical_mode(coeff, gy, sy, 1.0)

    bx, by = BAD_PARENTS
    _set_real_helical_mode(coeff, bx, sx, 0.35 * _unit_phase_for_child_work(bx, by, desired_phase_alignment=1.0))
    _set_real_helical_mode(coeff, by, sy, 0.35)

    nx, ny = NEGATIVE_PARENTS
    _set_real_helical_mode(coeff, nx, sx, 0.25 * _unit_phase_for_child_work(nx, ny, desired_phase_alignment=-1.0))
    _set_real_helical_mode(coeff, ny, sy, 0.25)

    coeff *= amp
    state = coeff * float(n**3)
    return _leray_dealias(state, k, k2, dealias)


@dataclass(frozen=True)
class MixedFateReservedYoungPDEProbe:
    status: str
    resolution: int
    cutoff: int
    steps: int
    snapshots: int
    mixed_fate_snapshots: int
    snapshots_with_good_work: int
    snapshots_with_bad_positive_work: int
    snapshots_with_negative_work: int
    initial_maximum_good_signed_efficiency: float
    initial_mixed_good_work: float
    worst_signed_ns_reconstruction_relative: float
    worst_hahn_pushforward_identity_relative: float
    worst_reservation_identity_relative: float
    maximum_bad_assistance_fraction_of_total_variation_upper: float
    stage_zero_first_time_failures: int
    geometry_good_marking_promotions: int


def run_probe(
    *,
    resolution: int = 24,
    steps: int = 16,
    viscosity: float = 0.03,
    amplitude: float = 1.0,
    duration: float = 0.001,
    snapshot_count: int = 3,
    tau: float = 0.1,
) -> MixedFateReservedYoungPDEProbe:
    n = int(resolution)
    count = int(steps)
    snaps = int(snapshot_count)
    if n < 24 or n % 2:
        raise ValueError("mixed-fate NS probe requires an even grid at least 24")
    if count < 16 or snaps < 2 or snaps > count + 1:
        raise ValueError("mixed-fate NS probe requires at least sixteen RK4 steps and two snapshots")
    k, k2, dealias, cutoff = _spectral_geometry(n, 7)
    state = adversarial_mixed_fate_initial_state(n, k, k2, dealias, amplitude=amplitude)
    dt = float(duration) / count
    sample_indices = tuple(sorted({round(j * count / (snaps - 1)) for j in range(snaps)}))

    mixed_snapshots = 0
    good_snapshots = 0
    bad_snapshots = 0
    negative_snapshots = 0
    initial_max_good = 0.0
    initial_mixed_good = 0.0
    worst_signed = 0.0
    worst_hahn = 0.0
    worst_reservation = 0.0
    max_bad_tv = 0.0
    first_time_failures = 0
    marking_promotions = 0

    for step in range(count + 1):
        nonlinear = _nonlinear_term(state, k, k2, dealias)
        if step in sample_indices:
            row, ledger = _snapshot_with_ledger(
                state,
                k,
                k2,
                dealias,
                cutoff,
                child_mode=CHILD,
                nonlinear_hat=nonlinear,
            )
            routing = route_canonical_positive_edge_work(
                ledger,
                tau=float(tau),
                mode_roles=single_hard_role_map(ledger),
            )
            good_snapshots += int(routing.good_positive_work > 0.0)
            bad_snapshots += int(routing.bad_positive_work > 0.0)
            negative_snapshots += int(ledger.negative_edge_work > 0.0)
            marking_promotions += int(routing.young_eligible.marking_good or routing.young_eligible.young_certified)
            if routing.bad_route is not None:
                first_time_failures += int(routing.bad_route.transfer_partition.first_time is not None)

            signed_scale = max(
                abs(float(row["actual_child_work"])),
                abs(ledger.signed_direct_work),
                ledger.positive_edge_work + ledger.negative_edge_work,
                1.0e-300,
            )
            worst_signed = max(
                worst_signed,
                abs(float(row["actual_child_work"]) - ledger.signed_direct_work) / signed_scale,
            )

            found_mixed = False
            for cell in routing.hard_cell_compression.cells:
                g = float(cell.inherited_good_positive_work)
                b = float(cell.inherited_bad_positive_work)
                nwork = inherited_negative_work(cell)
                T = float(cell.signed_work)
                scale = max(g + b + nwork, abs(T), 1.0e-300)
                worst_hahn = max(worst_hahn, abs(T - g - b + nwork) / scale)
                reservation = T - b
                worst_reservation = max(worst_reservation, abs(reservation - g + nwork) / scale)
                variation_upper = g + b + nwork
                if variation_upper > 0.0:
                    max_bad_tv = max(max_bad_tv, b / variation_upper)
                if g > 0.0 and b > 0.0:
                    found_mixed = True
            mixed_snapshots += int(found_mixed)

            if step == 0:
                if routing.good_positive_work <= 0.0:
                    raise AssertionError("engineered actual NS state failed to produce geometry-good positive work")
                if routing.bad_positive_work <= 0.0:
                    raise AssertionError("engineered actual NS state failed to produce terminal bad-positive work")
                if ledger.negative_edge_work <= 0.0:
                    raise AssertionError("engineered actual NS state failed to produce negative signed work")
                if not found_mixed:
                    raise AssertionError("engineered actual NS state failed to expose the mixed-fate hard-cell seam")
                initial_max_good = max(edge.signed_efficiency for edge in routing.good_support)
                initial_mixed_good = routing.young_eligible.unresolved_mixed_positive_work
                if not initial_max_good > 1.0 - 1.0e-4:
                    raise AssertionError("engineered near-extremal integer triad missed the canonical good threshold")
                if not initial_mixed_good > 0.0:
                    raise AssertionError("mixed-fate good dW+ was not retained as unresolved causal mass")

        if step < count:
            state = _rk4_step(state, dt, float(viscosity), k, k2, dealias)

    if mixed_snapshots == 0 or good_snapshots == 0 or bad_snapshots == 0 or negative_snapshots == 0:
        raise AssertionError("actual NS mixed-fate audit failed to exercise all signed-work fates")
    if first_time_failures or marking_promotions:
        raise AssertionError("actual NS mixed-fate audit invented a stage-zero clock or geometry-to-Young promotion")
    if max(worst_signed, worst_hahn, worst_reservation) > 4.0e-8:
        raise AssertionError("actual NS mixed-fate work provenance identity lost numerical accuracy")

    return MixedFateReservedYoungPDEProbe(
        status=STATUS,
        resolution=n,
        cutoff=cutoff,
        steps=count,
        snapshots=len(sample_indices),
        mixed_fate_snapshots=mixed_snapshots,
        snapshots_with_good_work=good_snapshots,
        snapshots_with_bad_positive_work=bad_snapshots,
        snapshots_with_negative_work=negative_snapshots,
        initial_maximum_good_signed_efficiency=initial_max_good,
        initial_mixed_good_work=initial_mixed_good,
        worst_signed_ns_reconstruction_relative=worst_signed,
        worst_hahn_pushforward_identity_relative=worst_hahn,
        worst_reservation_identity_relative=worst_reservation,
        maximum_bad_assistance_fraction_of_total_variation_upper=max_bad_tv,
        stage_zero_first_time_failures=first_time_failures,
        geometry_good_marking_promotions=marking_promotions,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resolution", type=int, default=24)
    parser.add_argument("--steps", type=int, default=16)
    parser.add_argument("--viscosity", type=float, default=0.03)
    parser.add_argument("--amplitude", type=float, default=1.0)
    parser.add_argument("--duration", type=float, default=0.001)
    parser.add_argument("--snapshots", type=int, default=3)
    parser.add_argument("--tau", type=float, default=0.1)
    parser.add_argument("--outdir", type=Path, default=Path("results-mixed-fate-reserved-young-handoff-pde"))
    args = parser.parse_args()
    out = run_probe(
        resolution=args.resolution,
        steps=args.steps,
        viscosity=args.viscosity,
        amplitude=args.amplitude,
        duration=args.duration,
        snapshot_count=args.snapshots,
        tau=args.tau,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "probe.json").write_text(json.dumps(asdict(out), indent=2, sort_keys=True) + "\n")
    summary = f"""# Actual Galerkin NS mixed-fate signed-cell audit\n\nStatus: **{STATUS}**.\n\nA real divergence-free Fourier polynomial is evolved by the repository's 2/3-dealiased Galerkin Navier--Stokes RK4 path.  At the selected physical child `z={CHILD}`, the initial state contains a near-extremal signed-good parent pair plus positive nonforward and negative-work parent pairs.  They are deliberately compressed into the same deterministic hard product cell only to stress the already-fixed edge-space Hahn provenance.\n\n- resolution/cutoff: `{out.resolution}` / `{out.cutoff}`\n- steps/snapshots: `{out.steps}` / `{out.snapshots}`\n- mixed-fate snapshots: `{out.mixed_fate_snapshots}`\n- snapshots with good/bad/negative work: `{out.snapshots_with_good_work}` / `{out.snapshots_with_bad_positive_work}` / `{out.snapshots_with_negative_work}`\n- initial maximum good signed efficiency: `{out.initial_maximum_good_signed_efficiency:.12g}`\n- initial unresolved mixed good dW+ mass: `{out.initial_mixed_good_work:.12g}`\n- worst actual-NS signed reconstruction relative residual: `{out.worst_signed_ns_reconstruction_relative:.3e}`\n- worst `T=g+b-n` hard-cell relative residual: `{out.worst_hahn_pushforward_identity_relative:.3e}`\n- worst `T-b=g-n` reservation-certificate relative residual: `{out.worst_reservation_identity_relative:.3e}`\n- stage-zero first-time failures: `{out.stage_zero_first_time_failures}`\n- geometry-good marking promotions: `{out.geometry_good_marking_promotions}`\n\nThe total-variation upper reported by the probe is only a sign/provenance diagnostic and is **not** substituted for the sharp continuous Young bound.  The analytic handoff theorem remains the statement about the actual full signed hard-cell trilinear form and a separately certified Young norm upper.  No global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
