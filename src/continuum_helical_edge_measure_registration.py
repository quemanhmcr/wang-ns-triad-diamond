from __future__ import annotations

import argparse
import itertools
import json
import math
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.coherent_service_or_flat import (
    coherent_flat_thresholds,
    coherent_service_or_flat_gate,
)
from src.helical import coupling_g, helical_basis, stable_norm3
from src.helical_physical_edge_registration import (
    STATUS as HELICAL_EDGE_STATUS,
)
from src.helical_physical_edge_registration import (
    HelicalPhysicalEdgeRegistration,
    _positive_product as _edge_positive_product,
    _signed_product as _edge_signed_product,
    leray_project,
    register_helical_physical_edge,
)
from src.physical_pair_weighted_productivity import physical_work_capacity_constant
from src.physical_transfer_defect_moat import (
    CLEAN_CHANGE_OF_MEASURE,
    GOOD_THRESHOLD,
    PHYSICAL_DEFECT_MEAN_FACTOR,
    physical_good_core_defect_mean_upper,
)
from src.single_edge_certificate import float_jstar
from src.smooth_flux_cocycle import (
    child_transfer_density_condition_number,
    near_extremal_gap_radius,
    polarization_certificate,
    positive_core_mass_lower_bound,
)
from src.triad_extremizer import symmetric_gamma, symmetric_rstar

STATUS = (
    "EXACT_CONTINUUM_HELICAL_EDGE_MEASURE_REGISTRATION__UNITARY_FOURIER__"
    "UNORDERED_PARENT_QUOTIENT__EXACT_HELICITY_RECONSTRUCTION__SIGNED_BEFORE_HAHN__"
    "NATIVE_CAPACITY_POLARIZATION__PHYSICAL_GOOD_CORE_CHANGE_OF_MEASURE"
)

UNITARY_FOURIER_CONVOLUTION_FACTOR = (2.0 * math.pi) ** (-1.5)
GOOD_CORE_ETA = float(GOOD_THRESHOLD)
LOW_COST_DEFICIT_CEILING = 1.0 / 20_000.0


def _complex_norm3(value: np.ndarray) -> float:
    q = np.asarray(value, dtype=complex)
    if q.shape != (3,) or np.any(~np.isfinite(q.real)) or np.any(~np.isfinite(q.imag)):
        raise ValueError("finite complex three-vector required")
    return float(math.hypot(*(abs(complex(x)) for x in q)))


def _finite_complex_scalar(value: complex, name: str) -> complex:
    z = complex(value)
    if not (math.isfinite(z.real) and math.isfinite(z.imag)):
        raise ValueError(f"{name} must be finite")
    return z


def _finite_sum(values: Sequence[float], name: str) -> float:
    vals = tuple(float(value) for value in values)
    if not all(math.isfinite(value) for value in vals):
        raise ValueError(f"{name} terms must be finite")
    try:
        out = math.fsum(vals)
    except OverflowError as exc:
        raise ValueError(f"{name} left the finite native range") from exc
    if not math.isfinite(out):
        raise ValueError(f"{name} left the finite native range")
    return out


def _relative_gap(actual: complex, expected: complex, *, scale: float | None = None) -> float:
    a = _finite_complex_scalar(actual, "actual native quantity")
    b = _finite_complex_scalar(expected, "expected native quantity")
    gap = abs(a - b)
    if not math.isfinite(gap):
        raise ValueError("native comparison gap must be finite")
    if scale is None:
        native = max(abs(a), abs(b))
    else:
        raw_scale = float(scale)
        if not math.isfinite(raw_scale):
            raise ValueError("native comparison scale must be finite")
        native = abs(raw_scale)
    if native == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    out = float(gap / native)
    if not math.isfinite(out):
        raise ValueError("native relative gap must be finite")
    return out


def _require_native_equal(
    name: str,
    actual: complex,
    expected: complex,
    *,
    scale: float | None = None,
    relative_tolerance: float = 5e-10,
) -> None:
    if _relative_gap(actual, expected, scale=scale) > relative_tolerance:
        raise AssertionError(f"{name} failed its native-scale identity")


def _relative_vector_gap(actual: np.ndarray, expected: np.ndarray, *, scale: float | None = None) -> float:
    a = _cvec3(actual, "actual vector")
    b = _cvec3(expected, "expected vector")
    gap = _complex_norm3(a - b)
    if scale is None:
        native = max(_complex_norm3(a), _complex_norm3(b))
    else:
        raw_scale = float(scale)
        if not math.isfinite(raw_scale):
            raise ValueError("native vector comparison scale must be finite")
        native = abs(raw_scale)
    if native == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    out = gap / native
    if not math.isfinite(out):
        raise ValueError("native vector relative gap must be finite")
    return out


def _physical_triad(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, float]:
    nx, ny, nz = map(stable_norm3, (x, y, z))
    if min(nx, ny, nz) == 0.0:
        raise ValueError("nonzero parent and child wavevectors required")
    if stable_norm3(x + y - z) > 2e-12 * max(nx, ny, nz):
        raise ValueError("physical triad requires z=x+y")
    return nx, ny, nz


def unitary_fourier_convolution_factor() -> float:
    """Product/convolution coefficient for the repository's unitary R^3 Fourier transform."""
    return UNITARY_FOURIER_CONVOLUTION_FACTOR


def unitary_sharp_young_physical_work_upper(child_wave_ratio: float) -> float:
    """Sharper unitary-Fourier version of the existing clean physical-work upper.

    ``physical_work_capacity_constant`` intentionally uses the clean larger bound
    ``4 R A3``.  The exact unitary Fourier product contributes the additional
    factor ``(2pi)^(-3/2)<1``.  Keeping these two constants distinct prevents an
    analytic upper bound from being mistaken for the native physical multiplier.
    """
    r = float(child_wave_ratio)
    if r <= 0.0 or not math.isfinite(r):
        raise ValueError("positive finite child wave ratio required")
    return unitary_fourier_convolution_factor() * physical_work_capacity_constant(r)


def _vec3(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (3,) or np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must be a finite real three-vector")
    return out


def _cvec3(value: np.ndarray, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=complex)
    if out.shape != (3,) or np.any(~np.isfinite(out.real)) or np.any(~np.isfinite(out.imag)):
        raise ValueError(f"{name} must be a finite complex three-vector")
    return out


def divergence_relative_residual(k: np.ndarray, value: np.ndarray) -> float:
    q = _vec3(k, "wavevector")
    v = _cvec3(value, "Fourier vector")
    nk = stable_norm3(q)
    if nk == 0.0:
        raise ValueError("nonzero wavevector required")
    nv = _complex_norm3(v)
    if nv == 0.0:
        return 0.0
    qhat = q / nk
    return float(abs(np.dot(qhat, v)) / nv)


def _require_divergence_free(k: np.ndarray, value: np.ndarray, name: str) -> np.ndarray:
    v = _cvec3(value, name)
    if divergence_relative_residual(k, v) > 2e-10:
        raise ValueError(f"{name} must be divergence free at its wavevector")
    return v


def helical_coefficients(k: np.ndarray, value: np.ndarray) -> dict[int, complex]:
    """Exact eventwise two-sector helical coordinates of a divergence-free Fourier vector."""
    q = _vec3(k, "wavevector")
    v = _require_divergence_free(q, value, "Fourier vector")
    if float(np.linalg.norm(q)) <= 0.0:
        raise ValueError("nonzero wavevector required")
    return {s: complex(np.vdot(helical_basis(q, s), v)) for s in (-1, 1)}


def helical_reconstruction(k: np.ndarray, value: np.ndarray) -> tuple[np.ndarray, float]:
    q = _vec3(k, "wavevector")
    v = _require_divergence_free(q, value, "Fourier vector")
    coeff = helical_coefficients(q, v)
    reconstructed = sum((coeff[s] * helical_basis(q, s) for s in (-1, 1)), np.zeros(3, complex))
    nv = _complex_norm3(v)
    gap = _complex_norm3(reconstructed - v)
    residual = 0.0 if nv == 0.0 and gap == 0.0 else (math.inf if nv == 0.0 else gap / nv)
    if residual > 3e-10:
        raise AssertionError("helical sectors failed to reconstruct a divergence-free Fourier vector")
    return reconstructed, residual


def ordered_parent_curl_source(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
) -> np.ndarray:
    """One ordered Fourier integrand ``P_z(u_x x omega_y)`` before parent quotient."""
    x = _vec3(x, "x")
    y = _vec3(y, "y")
    z = _vec3(z, "z")
    ux = _require_divergence_free(x, ux, "ux")
    uy = _require_divergence_free(y, uy, "uy")
    _physical_triad(x, y, z)
    omega_y = 1j * np.cross(y, uy)
    return leray_project(z, np.cross(ux, omega_y))


def unordered_parent_curl_source_vector(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
) -> np.ndarray:
    """Orientation-free two-parent orbit source on the quotient ``{x,y}``.

    The quotient measure is analytically ``lambda_z^unord=(1/2)(pi_z)_# dx``.
    Its orbit integrand is the sum of the two ordered parent terms.  No
    lexicographic representative is selected.
    """
    xy = ordered_parent_curl_source(x, y, z, ux, uy)
    yx = ordered_parent_curl_source(y, x, z, uy, ux)
    return xy + yx


def direct_vector_child_work(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
) -> float:
    z = _vec3(z, "z")
    uz = _require_divergence_free(z, uz, "uz")
    source = unordered_parent_curl_source_vector(x, y, z, ux, uy)
    return 2.0 * float(np.real(np.vdot(uz, source)))


@dataclass(frozen=True)
class ContinuumModalEdgeAtom:
    """One helicity-resolved edge with only its quotient base-measure mass added.

    All physical work/capacity/geometry/phase data are derived from the already
    certified one-edge registration object.  A caller cannot inject a desired
    capacity, work, multiplier or phase at this layer.
    """

    registration: HelicalPhysicalEdgeRegistration
    quotient_measure_mass: float

    def __post_init__(self) -> None:
        q = float(self.quotient_measure_mass)
        if not math.isfinite(q) or q < 0.0:
            raise ValueError("nonnegative finite unordered quotient-measure mass required")
        if self.registration.young_norm_used_as_capacity or self.registration.duhamel_weight_used_as_causal_law:
            raise ValueError("continuum atom inherited a forbidden capacity/causal replacement")
        if not self.registration.unordered_parent_orientation:
            raise ValueError("continuum atom must already be parent-orientation quotiented")

    @property
    def base_factor(self) -> float:
        return _edge_positive_product(
            (unitary_fourier_convolution_factor(), float(self.quotient_measure_mass)),
            "continuum quotient base factor",
        )

    @property
    def signed_work_mass(self) -> float:
        return _edge_signed_product(
            (self.base_factor, self.registration.signed_child_energy_work),
            "continuum signed work mass",
        )

    @property
    def capacity_mass(self) -> float:
        return _edge_positive_product(
            (self.base_factor, self.registration.native_modal_capacity),
            "continuum capacity mass",
        )

    @property
    def signed_progress_mass(self) -> float:
        return _edge_signed_product(
            (self.base_factor, self.registration.signed_upper_progress_work),
            "continuum signed progress mass",
        )

    @property
    def multiplier(self) -> float:
        return self.registration.normalized_multiplier

    @property
    def phase(self) -> float:
        return self.registration.phase_alignment

    @property
    def signed_efficiency(self) -> float:
        return self.multiplier * self.phase

    @property
    def scale_progress(self) -> float:
        return self.registration.scale_progress


@dataclass(frozen=True)
class ContinuumFiberRegistration:
    quotient_measure_mass: float
    direct_signed_work_density: float
    modal_signed_work_density: float
    direct_signed_progress_density: float
    modal_signed_progress_density: float
    ordered_quotient_source_residual: float
    parent_swap_residual: float
    helical_reconstruction_residual: float
    signed_work_reconstruction_residual: float
    signed_progress_reconstruction_residual: float
    modal_atoms: tuple[ContinuumModalEdgeAtom, ...]

    def __post_init__(self) -> None:
        if not math.isfinite(self.quotient_measure_mass) or self.quotient_measure_mass < 0.0:
            raise ValueError("nonnegative finite quotient mass required")
        vals = (
            self.direct_signed_work_density,
            self.modal_signed_work_density,
            self.direct_signed_progress_density,
            self.modal_signed_progress_density,
            self.ordered_quotient_source_residual,
            self.parent_swap_residual,
            self.helical_reconstruction_residual,
            self.signed_work_reconstruction_residual,
            self.signed_progress_reconstruction_residual,
        )
        if not all(math.isfinite(v) for v in vals):
            raise ValueError("finite continuum fiber certificate required")
        if len(self.modal_atoms) != 8:
            raise ValueError("all eight helical interaction sectors must be retained at the event")
        if any(v < 0.0 for v in (
            self.ordered_quotient_source_residual,
            self.parent_swap_residual,
            self.helical_reconstruction_residual,
        )):
            raise ValueError("continuum fiber provenance residuals must be nonnegative")
        q = float(self.quotient_measure_mass)
        for atom in self.modal_atoms:
            _require_native_equal(
                "continuum atom quotient-measure mass binding",
                atom.quotient_measure_mass,
                q,
                scale=max(abs(atom.quotient_measure_mass), abs(q)),
            )
        modal_work = float(sum(a.registration.signed_child_energy_work for a in self.modal_atoms))
        work_scale = max(
            abs(self.direct_signed_work_density),
            sum(abs(a.registration.signed_child_energy_work) for a in self.modal_atoms),
        )
        _require_native_equal(
            "continuum fiber modal work binding",
            self.modal_signed_work_density,
            modal_work,
            scale=work_scale,
        )
        _require_native_equal(
            "continuum fiber direct/modal work reconstruction",
            self.direct_signed_work_density,
            modal_work,
            scale=work_scale,
        )
        _require_native_equal(
            "continuum fiber stored work residual",
            self.signed_work_reconstruction_residual,
            self.direct_signed_work_density - self.modal_signed_work_density,
            scale=work_scale,
        )
        modal_progress = float(sum(a.registration.signed_upper_progress_work for a in self.modal_atoms))
        progress_scale = max(
            abs(self.direct_signed_progress_density),
            sum(abs(a.registration.signed_upper_progress_work) for a in self.modal_atoms),
        )
        _require_native_equal(
            "continuum fiber modal progress binding",
            self.modal_signed_progress_density,
            modal_progress,
            scale=progress_scale,
        )
        _require_native_equal(
            "continuum fiber direct/modal progress reconstruction",
            self.direct_signed_progress_density,
            modal_progress,
            scale=progress_scale,
        )
        _require_native_equal(
            "continuum fiber stored progress residual",
            self.signed_progress_reconstruction_residual,
            self.direct_signed_progress_density - self.modal_signed_progress_density,
            scale=progress_scale,
        )


def register_continuum_triad_fiber(
    *,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    uz: np.ndarray,
    quotient_measure_mass: float,
) -> ContinuumFiberRegistration:
    """Resolve one arbitrary divergence-free continuum triad fiber into 8 physical helical edges."""
    x = _vec3(x, "x")
    y = _vec3(y, "y")
    z = _vec3(z, "z")
    ux = _require_divergence_free(x, ux, "ux")
    uy = _require_divergence_free(y, uy, "uy")
    uz = _require_divergence_free(z, uz, "uz")
    qmass = float(quotient_measure_mass)
    if not math.isfinite(qmass) or qmass < 0.0:
        raise ValueError("nonnegative finite quotient-measure mass required")
    nx, ny, nz = _physical_triad(x, y, z)

    _, rx = helical_reconstruction(x, ux)
    _, ry = helical_reconstruction(y, uy)
    _, rz = helical_reconstruction(z, uz)
    hres = max(rx, ry, rz)

    source_xy = ordered_parent_curl_source(x, y, z, ux, uy)
    source_yx = ordered_parent_curl_source(y, x, z, uy, ux)
    source_unordered = unordered_parent_curl_source_vector(x, y, z, ux, uy)
    source_scale = max(_complex_norm3(source_unordered), _complex_norm3(source_xy) + _complex_norm3(source_yx))
    source_res = _relative_vector_gap(source_unordered, source_xy + source_yx, scale=source_scale)
    if source_res > 3e-11:
        raise AssertionError("unordered quotient orbit failed ordered convolution reconstruction")

    swapped = unordered_parent_curl_source_vector(y, x, z, uy, ux)
    swap_res = _relative_vector_gap(swapped, source_unordered)
    if swap_res > 3e-11:
        raise AssertionError("unordered parent quotient depended on parent orientation")

    cx = helical_coefficients(x, ux)
    cy = helical_coefficients(y, uy)
    cz = helical_coefficients(z, uz)
    atoms: list[ContinuumModalEdgeAtom] = []
    for sx, sy, sz in itertools.product((-1, 1), repeat=3):
        reg = register_helical_physical_edge(
            x=x,
            y=y,
            z=z,
            sx=sx,
            sy=sy,
            sz=sz,
            ax=cx[sx],
            ay=cy[sy],
            az=cz[sz],
        )
        atoms.append(ContinuumModalEdgeAtom(reg, qmass))

    direct_work = 2.0 * float(np.real(np.vdot(uz, source_unordered)))
    modal_work = float(sum(a.registration.signed_child_energy_work for a in atoms))
    work_res = direct_work - modal_work
    work_scale = max(abs(direct_work), sum(abs(a.registration.signed_child_energy_work) for a in atoms))
    _require_native_equal(
        "eight-helicity direct vector child-work reconstruction",
        direct_work,
        modal_work,
        scale=work_scale,
        relative_tolerance=5e-10,
    )

    progress = max(
        0.0,
        math.log(nz / max(nx, ny)),
    )
    direct_progress = direct_work * progress
    modal_progress = float(sum(a.registration.signed_upper_progress_work for a in atoms))
    progress_res = direct_progress - modal_progress
    progress_scale = max(abs(direct_progress), sum(abs(a.registration.signed_upper_progress_work) for a in atoms))
    _require_native_equal(
        "eight-helicity upper-progress reconstruction",
        direct_progress,
        modal_progress,
        scale=progress_scale,
        relative_tolerance=5e-10,
    )

    return ContinuumFiberRegistration(
        quotient_measure_mass=qmass,
        direct_signed_work_density=direct_work,
        modal_signed_work_density=modal_work,
        direct_signed_progress_density=direct_progress,
        modal_signed_progress_density=modal_progress,
        ordered_quotient_source_residual=source_res,
        parent_swap_residual=swap_res,
        helical_reconstruction_residual=hres,
        signed_work_reconstruction_residual=work_res,
        signed_progress_reconstruction_residual=progress_res,
        modal_atoms=tuple(atoms),
    )


@dataclass(frozen=True)
class ContinuumEdgeMeasureLedger:
    fibers: int
    modal_edges: int
    physical_fibers: tuple[ContinuumFiberRegistration, ...]
    quotient_measure_mass: float
    signed_direct_work: float
    signed_modal_work: float
    positive_edge_work: float
    negative_edge_work: float
    aggregate_positive_work: float
    fiber_positive_work: float
    positive_dominance_over_aggregate: float
    positive_dominance_over_fibers: float
    positive_forward_work: float
    positive_nonforward_work: float
    capacity_mass: float
    signed_direct_progress: float
    signed_registered_progress: float
    normalized_signed_flux: float
    block_transfer_deficit: float
    multiplier_deficit: float
    phase_deficit: float
    polarization_residual: float
    good_core_capacity_mass: float
    good_core_positive_work: float
    good_core_capacity_fraction: float
    good_core_physical_to_capacity_rn_min: float | None
    good_core_physical_to_capacity_rn_max: float | None
    direct_work_reconstruction_residual: float
    direct_progress_reconstruction_residual: float
    capacity_is_causal_law: bool = False
    causal_law: str = "positive_Hahn_part_of_actual_child_energy_edge_work"
    parent_orientation_chosen: bool = False

    def __post_init__(self) -> None:
        if self.fibers <= 0 or self.modal_edges <= 0:
            raise ValueError("positive continuum fiber and modal-edge counts required")
        if len(self.physical_fibers) != self.fibers or sum(len(f.modal_atoms) for f in self.physical_fibers) != self.modal_edges:
            raise ValueError("continuum ledger physical-fiber provenance count mismatch")
        numeric = (
            self.quotient_measure_mass,
            self.signed_direct_work,
            self.signed_modal_work,
            self.positive_edge_work,
            self.negative_edge_work,
            self.aggregate_positive_work,
            self.fiber_positive_work,
            self.positive_dominance_over_aggregate,
            self.positive_dominance_over_fibers,
            self.positive_forward_work,
            self.positive_nonforward_work,
            self.capacity_mass,
            self.signed_direct_progress,
            self.signed_registered_progress,
            self.normalized_signed_flux,
            self.block_transfer_deficit,
            self.multiplier_deficit,
            self.phase_deficit,
            self.polarization_residual,
            self.good_core_capacity_mass,
            self.good_core_positive_work,
            self.good_core_capacity_fraction,
            self.direct_work_reconstruction_residual,
            self.direct_progress_reconstruction_residual,
        )
        if not all(math.isfinite(float(value)) for value in numeric):
            raise ValueError("finite continuum edge-measure ledger required")
        for value in (self.good_core_physical_to_capacity_rn_min, self.good_core_physical_to_capacity_rn_max):
            if value is not None and not math.isfinite(float(value)):
                raise ValueError("finite Radon-Nikodym provenance required when present")
        if self.quotient_measure_mass < 0.0:
            raise ValueError("nonnegative finite quotient-measure mass required")
        if self.capacity_mass <= 0.0:
            raise ValueError("positive continuum capacity mass required")
        if self.positive_edge_work < 0.0 or self.negative_edge_work < 0.0:
            raise ValueError("Hahn masses must be nonnegative")
        if self.capacity_is_causal_law or self.parent_orientation_chosen:
            raise ValueError("capacity/orientation bookkeeping was promoted to physical causality")
        if self.block_transfer_deficit < -5e-10 or self.block_transfer_deficit > 2.0 + 5e-10:
            raise AssertionError("normalized signed capacity deficit left [0,2]")


def continuum_edge_measure_ledger(fibers: Sequence[ContinuumFiberRegistration]) -> ContinuumEdgeMeasureLedger:
    """Integrate the exact signed helical edge law against a positive quotient base measure.

    This finite sequence is a quadrature/atomic representation of the underlying
    Radon measure.  The theorem statement is measure-level: signed reconstruction
    comes first, then Hahn splitting.  In particular positive edge mass dominates
    aggregate positive work but is not identified with it under cancellation.
    """
    fs = tuple(fibers)
    if not fs:
        raise ValueError("at least one continuum triad fiber required")
    atoms = tuple(a for f in fs for a in f.modal_atoms)
    cf = unitary_fourier_convolution_factor()

    signed_direct = _finite_sum(
        (
            _edge_signed_product(
                (cf, f.quotient_measure_mass, f.direct_signed_work_density),
                "continuum direct signed-work fiber mass",
            )
            for f in fs
        ),
        "continuum direct signed work",
    )
    signed_modal = _finite_sum((a.signed_work_mass for a in atoms), "continuum modal signed work")
    direct_work_res = signed_direct - signed_modal
    signed_work_scale = max(
        abs(signed_direct),
        _finite_sum((abs(a.signed_work_mass) for a in atoms), "continuum absolute modal work scale"),
    )
    _require_native_equal(
        "continuum signed helical edge measure direct NS work",
        signed_direct,
        signed_modal,
        scale=signed_work_scale,
        relative_tolerance=7e-10,
    )

    positive = _finite_sum((max(a.signed_work_mass, 0.0) for a in atoms), "continuum positive Hahn mass")
    negative = _finite_sum((max(-a.signed_work_mass, 0.0) for a in atoms), "continuum negative Hahn mass")
    hahn_res = (positive - negative) - signed_modal
    _require_native_equal(
        "continuum Hahn split signed edge reconstruction",
        positive - negative,
        signed_modal,
        scale=positive + negative,
        relative_tolerance=7e-11,
    )
    aggregate_positive = max(0.0, signed_direct)
    fiber_positive = _finite_sum(
        (
            _edge_positive_product(
                (cf, f.quotient_measure_mass, max(f.direct_signed_work_density, 0.0)),
                "continuum fiber-positive work mass",
            )
            for f in fs
        ),
        "continuum fiber-positive work",
    )
    positive_scale = max(positive, fiber_positive, aggregate_positive)
    positive_slack = 7e-11 * positive_scale
    if positive + positive_slack < fiber_positive or fiber_positive + positive_slack < aggregate_positive:
        raise AssertionError("positive edge Hahn mass failed physical positive-work dominance")

    capacity = _finite_sum((a.capacity_mass for a in atoms), "continuum total capacity mass")
    if capacity <= 0.0:
        raise ValueError("positive modal capacity required")
    direct_progress = _finite_sum(
        (
            _edge_signed_product(
                (cf, f.quotient_measure_mass, f.direct_signed_progress_density),
                "continuum direct signed-progress fiber mass",
            )
            for f in fs
        ),
        "continuum direct signed progress",
    )
    registered_progress = _finite_sum(
        (a.signed_progress_mass for a in atoms), "continuum registered signed progress"
    )
    direct_progress_res = direct_progress - registered_progress
    progress_scale = max(
        abs(direct_progress),
        _finite_sum((abs(a.signed_progress_mass) for a in atoms), "continuum absolute progress scale"),
        _edge_positive_product((capacity, float_jstar()), "continuum J-star capacity scale"),
    )
    _require_native_equal(
        "continuum A*J*c direct upper-progress work",
        direct_progress,
        registered_progress,
        scale=progress_scale,
        relative_tolerance=7e-10,
    )

    capacities = np.asarray([a.capacity_mass for a in atoms], dtype=float)
    multipliers = np.asarray([a.multiplier for a in atoms], dtype=float)
    phases = np.asarray([a.phase for a in atoms], dtype=float)
    pol = polarization_certificate(capacities, multipliers, phases)
    jstar = float_jstar()
    ratio_raw = registered_progress / (jstar * capacity)
    if ratio_raw < -1.0 - 5e-10 or ratio_raw > 1.0 + 5e-10:
        raise AssertionError("continuum signed progress exceeded the global capacity envelope")
    ratio = max(-1.0, min(1.0, ratio_raw))
    if abs(ratio - pol.normalized_signed_flux) > 7e-11:
        raise AssertionError("measure-level polarization ratio disagreed with actual A*J*c progress")
    deficit = 1.0 - ratio

    positive_forward = _finite_sum(
        (max(a.signed_work_mass, 0.0) for a in atoms if a.scale_progress > 0.0),
        "continuum positive forward work",
    )
    positive_nonforward = _finite_sum(
        (max(a.signed_work_mass, 0.0) for a in atoms if a.scale_progress <= 0.0),
        "continuum positive nonforward work",
    )

    good = [a for a in atoms if a.capacity_mass > 0.0 and a.signed_efficiency > 1.0 - GOOD_CORE_ETA]
    good_capacity = _finite_sum((a.capacity_mass for a in good), "continuum good-core capacity mass")
    good_work = _finite_sum((a.signed_work_mass for a in good), "continuum good-core physical work")
    rn_min: float | None = None
    rn_max: float | None = None
    if good_capacity > 0.0:
        if good_work <= 0.0:
            raise AssertionError("signed-good capacity core did not carry positive physical work")
        rn_values = [
            (a.signed_work_mass / good_work) / (a.capacity_mass / good_capacity)
            for a in good
        ]
        rn_min = min(rn_values)
        rn_max = max(rn_values)

    return ContinuumEdgeMeasureLedger(
        fibers=len(fs),
        modal_edges=len(atoms),
        physical_fibers=fs,
        quotient_measure_mass=_finite_sum((f.quotient_measure_mass for f in fs), "continuum quotient-measure mass"),
        signed_direct_work=signed_direct,
        signed_modal_work=signed_modal,
        positive_edge_work=positive,
        negative_edge_work=negative,
        aggregate_positive_work=aggregate_positive,
        fiber_positive_work=fiber_positive,
        positive_dominance_over_aggregate=positive - aggregate_positive,
        positive_dominance_over_fibers=positive - fiber_positive,
        positive_forward_work=positive_forward,
        positive_nonforward_work=positive_nonforward,
        capacity_mass=capacity,
        signed_direct_progress=direct_progress,
        signed_registered_progress=registered_progress,
        normalized_signed_flux=ratio,
        block_transfer_deficit=deficit,
        multiplier_deficit=pol.multiplier_deficit,
        phase_deficit=pol.phase_deficit,
        polarization_residual=pol.exact_residual,
        good_core_capacity_mass=good_capacity,
        good_core_positive_work=good_work,
        good_core_capacity_fraction=good_capacity / capacity,
        good_core_physical_to_capacity_rn_min=rn_min,
        good_core_physical_to_capacity_rn_max=rn_max,
        direct_work_reconstruction_residual=direct_work_res,
        direct_progress_reconstruction_residual=direct_progress_res,
    )


def _replay_physical_ledger(ledger: ContinuumEdgeMeasureLedger) -> ContinuumEdgeMeasureLedger:
    """Recompute summary observables from the bound physical fiber law.

    A typed summary is not continuation/provenance authority.  Downstream gates
    must replay the actual registered fibers so ``dataclasses.replace`` or an
    equivalent forged summary cannot manufacture low deficit or a good core.
    """
    replayed = continuum_edge_measure_ledger(tuple(ledger.physical_fibers))
    numeric_fields = (
        "quotient_measure_mass",
        "signed_direct_work",
        "signed_modal_work",
        "positive_edge_work",
        "negative_edge_work",
        "aggregate_positive_work",
        "fiber_positive_work",
        "positive_forward_work",
        "positive_nonforward_work",
        "capacity_mass",
        "signed_direct_progress",
        "signed_registered_progress",
        "normalized_signed_flux",
        "block_transfer_deficit",
        "multiplier_deficit",
        "phase_deficit",
        "polarization_residual",
        "good_core_capacity_mass",
        "good_core_positive_work",
        "good_core_capacity_fraction",
    )
    for name in numeric_fields:
        _require_native_equal(
            f"continuum ledger replay field {name}",
            getattr(ledger, name),
            getattr(replayed, name),
            relative_tolerance=8e-10,
        )
    for name in ("good_core_physical_to_capacity_rn_min", "good_core_physical_to_capacity_rn_max"):
        actual = getattr(ledger, name)
        expected = getattr(replayed, name)
        if actual is None or expected is None:
            if actual is not expected:
                raise AssertionError(f"continuum ledger replay field {name} lost physical provenance")
        else:
            _require_native_equal(
                f"continuum ledger replay field {name}",
                actual,
                expected,
                relative_tolerance=8e-10,
            )
    if ledger.capacity_is_causal_law != replayed.capacity_is_causal_law or ledger.causal_law != replayed.causal_law or ledger.parent_orientation_chosen != replayed.parent_orientation_chosen:
        raise AssertionError("continuum ledger replay changed causal/orientation provenance")
    return replayed


@dataclass(frozen=True)
class PhysicalGoodCoreCertificate:
    block_transfer_deficit: float
    eta: float
    certified_capacity_fraction_lower: float
    realized_capacity_fraction: float
    gap_radius_upper: float
    density_condition_number: float
    clean_normalized_rn_lower: float
    clean_normalized_rn_upper: float
    realized_normalized_rn_min: float
    realized_normalized_rn_max: float
    physical_hodge_defect_mean_upper: float


def signed_good_core_physical_law(ledger: ContinuumEdgeMeasureLedger) -> PhysicalGoodCoreCertificate:
    """Convert the same low-deficit capacity law to actual positive child work on its good core."""
    ledger = _replay_physical_ledger(ledger)
    eps = float(ledger.block_transfer_deficit)
    if not (0.0 <= eps < LOW_COST_DEFICIT_CEILING):
        raise ValueError("physical good-core change of measure requires deficit <1/20000")
    lower_mass = positive_core_mass_lower_bound(eps, GOOD_CORE_ETA)
    if lower_mass <= 0.5:
        raise AssertionError("low-cost continuum block lost the certified half-capacity good core")
    if ledger.good_core_capacity_fraction + 5e-10 < lower_mass:
        raise AssertionError("realized continuum good core violated Markov capacity lower bound")
    if ledger.good_core_physical_to_capacity_rn_min is None or ledger.good_core_physical_to_capacity_rn_max is None:
        raise AssertionError("low-cost continuum block produced no realized physical signed-good core")

    gamma = symmetric_gamma(symmetric_rstar())
    gap = near_extremal_gap_radius(GOOD_CORE_ETA)
    cond = child_transfer_density_condition_number(GOOD_CORE_ETA, gamma)
    if not cond < CLEAN_CHANGE_OF_MEASURE:
        raise AssertionError("single-edge physical/capacity condition no longer lies below 53/50")
    clean_lo = 1.0 / CLEAN_CHANGE_OF_MEASURE
    clean_hi = CLEAN_CHANGE_OF_MEASURE
    rn_lo = float(ledger.good_core_physical_to_capacity_rn_min)
    rn_hi = float(ledger.good_core_physical_to_capacity_rn_max)
    if rn_lo < clean_lo - 5e-9 or rn_hi > clean_hi + 5e-9:
        raise AssertionError("realized normalized physical/capacity Radon-Nikodym ratio left the clean 50/53..53/50 band")
    mean_upper = physical_good_core_defect_mean_upper(eps)
    if abs(mean_upper - PHYSICAL_DEFECT_MEAN_FACTOR * eps) > 2e-15 * max(1.0, mean_upper):
        raise AssertionError("physical defect mean upper changed from the certified 106/25 law")

    return PhysicalGoodCoreCertificate(
        block_transfer_deficit=eps,
        eta=GOOD_CORE_ETA,
        certified_capacity_fraction_lower=lower_mass,
        realized_capacity_fraction=ledger.good_core_capacity_fraction,
        gap_radius_upper=gap,
        density_condition_number=cond,
        clean_normalized_rn_lower=clean_lo,
        clean_normalized_rn_upper=clean_hi,
        realized_normalized_rn_min=rn_lo,
        realized_normalized_rn_max=rn_hi,
        physical_hodge_defect_mean_upper=mean_upper,
    )


def edge_measure_to_service_or_flat(
    ledger: ContinuumEdgeMeasureLedger,
    *,
    tau: float,
    objective_variation_action: float,
    total_strain_action: float,
    coherent_deformation_action: float,
    aspect: float,
    scale_radius: float,
    has_predecessor: bool,
    scaled_lifetime: float,
    phase_holonomy: float = 0.0,
) -> dict[str, object]:
    """Feed the coherent gate the deficit computed from the actual continuum edge ledger.

    The caller is deliberately not allowed to supply ``avg_transfer_deficit``.
    """
    ledger = _replay_physical_ledger(ledger)
    out = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=ledger.block_transfer_deficit,
        objective_variation_action=objective_variation_action,
        total_strain_action=total_strain_action,
        coherent_deformation_action=coherent_deformation_action,
        aspect=aspect,
        scale_radius=scale_radius,
        has_predecessor=has_predecessor,
        scaled_lifetime=scaled_lifetime,
        phase_holonomy=phase_holonomy,
    )
    return {
        **out,
        "edge_measure_transfer_deficit": ledger.block_transfer_deficit,
        "transfer_deficit_source": "actual_continuum_signed_A_J_c_edge_measure",
        "capacity_used_as_causal_law": False,
    }


def theorem_certificate() -> dict[str, object]:
    tau = 0.01
    service_delta = float(coherent_flat_thresholds(tau)["block_transfer_deficit"])
    gamma = symmetric_gamma(symmetric_rstar())
    gap = near_extremal_gap_radius(GOOD_CORE_ETA)
    cond = child_transfer_density_condition_number(GOOD_CORE_ETA, gamma)
    clean_upper = physical_work_capacity_constant(1.0)
    unitary_upper = unitary_sharp_young_physical_work_upper(1.0)
    if not unitary_upper < clean_upper:
        raise AssertionError("unitary Fourier physical-work upper should be strictly sharper than the clean 4RA3 bound")
    if not service_delta < LOW_COST_DEFICIT_CEILING:
        raise AssertionError("coherent transfer gate no longer lies inside the physical good-core regime")
    return {
        "status": STATUS,
        "upstream_one_edge_status": HELICAL_EDGE_STATUS,
        "unitary_fourier": "fhat=(2pi)^(-3/2) integral exp(-ix.k) f(x)dx, so product convolution carries C_F=(2pi)^(-3/2)",
        "parent_quotient": "for fixed child z, pi_z(x)={x,z-x} and lambda_z^unord=(1/2)(pi_z)_# dx; the orbit integrand is the sum of both ordered parent terms, so no orientation selector is physical",
        "helicity": "arbitrary divergence-free parent/child Fourier vectors resolve into the two orthogonal helical sectors at the event; summing all 8 (sx,sy,sz) work channels reconstructs direct vector NS work exactly",
        "signed_measure": "dW=C_F T_e d(lambda_unord) is constructed signed before Hahn; W_plus-W_minus=W while W_plus >= [W]_+ under cancellation",
        "capacity_measure": "dA=C_F A_e d(lambda_unord), A_e=4|z||a_x a_y a_z|; dA is a positive reference measure, never the causal child-work law",
        "progress_measure": "dF=g_scale dW=J_e c_e dA and R=F/(J_* A); 1-R=E_A[(1-m)+m(1-c)] exactly",
        "nonforward": "positive child work with g_scale=J=0 remains in dW^+ and is not relabeled signed-good; it contributes full multiplier deficit to the capacity reference",
        "good_core": f"for eta={GOOD_CORE_ETA:g}, low deficit gives capacity core r=m c>1-eta; |g-gamma*|<={gap:.12g} and normalized physical/capacity laws differ by condition {cond:.12g}<53/50",
        "moat_handoff": "on deficit <1/20000, same physical good-core law supplies E_phys[D]<=(106/25) epsilon to physical_transfer_defect_moat",
        "service_handoff": "coherent_service_or_flat receives epsilon=1-F/(J_*A) only from this physical ledger; the caller cannot inject an avg_transfer_deficit scalar",
        "normalization_distinction": f"native unitary Young work upper at R=1 is C_F*(4A3)={unitary_upper:.12g}, while existing clean productivity upper 4A3={clean_upper:.12g} deliberately dominates it",
        "service_default_deficit_threshold": service_delta,
        "young_distinction": "epsilon=1-F/(J_*A) is the edge geometry/phase signed-efficiency deficit relative to actual modal capacity; it is not the separate Young norm-saturation deficit, which remains downstream",
        "scope": "this closes continuum signed edge-measure registration and its low-deficit physical-law handoff; it does not prove that every generic HH block is low deficit or terminate nonforward/high-deficit physical events",
    }


@dataclass(frozen=True)
class ContinuumEdgeMeasureStress:
    samples: int
    core_blocks: int
    worst_helical_reconstruction_residual: float
    worst_ordered_quotient_residual: float
    worst_parent_swap_residual: float
    worst_direct_modal_work_relative_residual: float
    worst_direct_modal_progress_relative_residual: float
    minimum_hahn_aggregate_margin: float
    minimum_hahn_fiber_margin: float
    worst_polarization_residual: float
    maximum_signed_flux_ratio: float
    minimum_signed_flux_ratio: float
    positive_nonforward_blocks: int
    minimum_good_core_capacity_margin: float
    minimum_good_core_rn_lower_margin: float
    minimum_good_core_rn_upper_margin: float
    minimum_service_transfer_cost_margin: float
    low_deficit_flat_blocks: int


def _relative(a: float, b: float) -> float:
    return _relative_gap(a, b)


def _random_divergence_free(rng: np.random.Generator, k: np.ndarray) -> np.ndarray:
    for _ in range(20):
        raw = rng.normal(size=3) + 1j * rng.normal(size=3)
        v = leray_project(k, raw)
        if np.linalg.norm(v) > 1e-5:
            return v
    raise AssertionError("failed to sample a nonzero divergence-free Fourier vector")


def _symmetric_extremal_fiber(
    *,
    quotient_measure_mass: float = 1.0,
    u_residual: float = 0.0,
    v_residual: float = 0.0,
) -> ContinuumFiberRegistration:
    rstar = symmetric_rstar()
    gamma = symmetric_gamma(rstar)
    nx = math.exp(-(gamma + v_residual) - 0.5 * u_residual)
    ny = math.exp(-(gamma + v_residual) + 0.5 * u_residual)
    xx = 0.5 * (1.0 + nx * nx - ny * ny)
    yy2 = nx * nx - xx * xx
    if yy2 <= 0.0:
        raise ValueError("residuals left the physical triangle domain")
    yy = math.sqrt(yy2)
    x = np.array([xx, yy, 0.0])
    y = np.array([1.0 - xx, -yy, 0.0])
    z = x + y
    sx, sy, sz = 1, -1, 1
    ax = 1.0 + 0.0j
    ay = 1.0 + 0.0j
    g = coupling_g(x, y, -z, sx, sy, sz)
    signed_frequency = sx * nx - sy * ny
    target_sign = 1.0 if signed_frequency >= 0.0 else -1.0
    az = target_sign * np.exp(-1j * np.angle(g))
    ux = ax * helical_basis(x, sx)
    uy = ay * helical_basis(y, sy)
    uz = az * helical_basis(z, sz)
    return register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=ux,
        uy=uy,
        uz=uz,
        quotient_measure_mass=quotient_measure_mass,
    )


def _nonforward_positive_fiber(quotient_measure_mass: float = 1.0) -> ContinuumFiberRegistration:
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([-0.8, 0.6, 0.0])
    z = x + y
    sx, sy, sz = 1, -1, 1
    ax = 1.0 + 0.0j
    ay = 1.0 + 0.0j
    g = coupling_g(x, y, -z, sx, sy, sz)
    signed_frequency = sx * np.linalg.norm(x) - sy * np.linalg.norm(y)
    target_sign = 1.0 if signed_frequency >= 0.0 else -1.0
    az = target_sign * np.exp(-1j * np.angle(g))
    return register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=ax * helical_basis(x, sx),
        uy=ay * helical_basis(y, sy),
        uz=az * helical_basis(z, sz),
        quotient_measure_mass=quotient_measure_mass,
    )


def stress(samples: int = 3_000, core_blocks: int = 100, seed: int = 20260811) -> ContinuumEdgeMeasureStress:
    if samples <= 0 or core_blocks <= 0:
        raise ValueError("positive stress sizes required")
    rng = np.random.default_rng(seed)
    wh = wo = ws = ww = wp = wpol = 0.0
    mha = mhf = mc = mrl = mru = msvc = float("inf")
    maxr = -math.inf
    minr = math.inf
    nnon = 0
    flat_blocks = 0
    batch: list[ContinuumFiberRegistration] = []

    for i in range(samples):
        while True:
            x = rng.normal(size=3)
            y = rng.normal(size=3)
            z = x + y
            if min(np.linalg.norm(x), np.linalg.norm(y), np.linalg.norm(z)) > 0.2 and np.linalg.norm(np.cross(x, y)) > 0.08:
                break
        ux = _random_divergence_free(rng, x)
        uy = _random_divergence_free(rng, y)
        uz = _random_divergence_free(rng, z)
        q = float(np.exp(rng.uniform(-2.0, 2.0)))
        f = register_continuum_triad_fiber(x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, quotient_measure_mass=q)
        wh = max(wh, f.helical_reconstruction_residual)
        wo = max(wo, f.ordered_quotient_source_residual)
        ws = max(ws, f.parent_swap_residual)
        ww = max(ww, _relative(f.direct_signed_work_density, f.modal_signed_work_density))
        wp = max(wp, _relative(f.direct_signed_progress_density, f.modal_signed_progress_density))
        batch.append(f)
        if len(batch) == 5 or i == samples - 1:
            led = continuum_edge_measure_ledger(tuple(batch))
            mha = min(mha, led.positive_dominance_over_aggregate)
            mhf = min(mhf, led.positive_dominance_over_fibers)
            wpol = max(wpol, abs(led.polarization_residual))
            maxr = max(maxr, led.normalized_signed_flux)
            minr = min(minr, led.normalized_signed_flux)
            if led.positive_nonforward_work > 0.0:
                nnon += 1
            batch.clear()

    # Deterministic physical nonforward positive work must survive with zero progress.
    non = continuum_edge_measure_ledger((_nonforward_positive_fiber(),))
    if non.positive_nonforward_work <= 0.0 or non.positive_forward_work > 5e-12:
        raise AssertionError("nonforward positive physical work was lost or mislabeled forward")
    service_non = edge_measure_to_service_or_flat(
        non,
        tau=0.01,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=True,
        scaled_lifetime=1.0,
    )
    roots = tuple(service_non.get("triggered_causes", ()))
    transfer_root = next((r for r in roots if r["cause"] == "physical_transfer_cost"), None)
    if transfer_root is None:
        raise AssertionError("native nonforward capacity deficit did not enter physical transfer-cost branch")
    msvc = min(msvc, float(transfer_root["value"]) - float(transfer_root["threshold"]))

    # Low-deficit blocks use actual helical fibers near the certified extremizer.
    for _ in range(core_blocks):
        fs: list[ContinuumFiberRegistration] = []
        for _ in range(5):
            u = float(rng.uniform(0.0, 2e-6))
            v = float(rng.uniform(-2e-6, 2e-6))
            q = float(np.exp(rng.uniform(-1.0, 1.0)))
            fs.append(_symmetric_extremal_fiber(quotient_measure_mass=q, u_residual=u, v_residual=v))
        led = continuum_edge_measure_ledger(tuple(fs))
        if led.block_transfer_deficit >= LOW_COST_DEFICIT_CEILING:
            raise AssertionError("near-extremal physical block left low-deficit good-core regime")
        core = signed_good_core_physical_law(led)
        mc = min(mc, core.realized_capacity_fraction - core.certified_capacity_fraction_lower)
        mrl = min(mrl, core.realized_normalized_rn_min - core.clean_normalized_rn_lower)
        mru = min(mru, core.clean_normalized_rn_upper - core.realized_normalized_rn_max)

    # Exact symmetric edge should feed the coherent gate without caller-supplied transfer deficit.
    exact = continuum_edge_measure_ledger((_symmetric_extremal_fiber(),))
    flat = edge_measure_to_service_or_flat(
        exact,
        tau=0.01,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=True,
        scaled_lifetime=1.0,
    )
    if flat["status"] == "coherent_kelvin_extremal_flat":
        flat_blocks += 1
    else:
        # Floating equality may retain a tiny transfer-cost root; it must be no
        # larger than numerical tolerance around the exact zero-deficit state.
        causes = tuple(flat.get("triggered_causes", ()))
        tr = next((r for r in causes if r["cause"] == "physical_transfer_cost"), None)
        if tr is None or float(tr["value"]) > 5e-10:
            raise AssertionError("exact extremal physical edge did not reach flat gate modulo numerical equality tolerance")

    return ContinuumEdgeMeasureStress(
        samples=samples,
        core_blocks=core_blocks,
        worst_helical_reconstruction_residual=wh,
        worst_ordered_quotient_residual=wo,
        worst_parent_swap_residual=ws,
        worst_direct_modal_work_relative_residual=ww,
        worst_direct_modal_progress_relative_residual=wp,
        minimum_hahn_aggregate_margin=mha,
        minimum_hahn_fiber_margin=mhf,
        worst_polarization_residual=wpol,
        maximum_signed_flux_ratio=maxr,
        minimum_signed_flux_ratio=minr,
        positive_nonforward_blocks=nnon,
        minimum_good_core_capacity_margin=mc,
        minimum_good_core_rn_lower_margin=mrl,
        minimum_good_core_rn_upper_margin=mru,
        minimum_service_transfer_cost_margin=msvc,
        low_deficit_flat_blocks=flat_blocks,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=3_000)
    ap.add_argument("--core-blocks", type=int, default=100)
    ap.add_argument("--outdir", type=Path, default=Path("results-continuum-helical-edge-measure"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples, args.core_blocks)
    (args.outdir / "continuum_helical_edge_measure_registration.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Continuum helical edge-measure registration: signed NS work before Hahn

Status: **{cert['status']}**.

The unitary Fourier convention contributes `C_F=(2pi)^(-3/2)`.  For fixed child
`z`, quotient the ordered parent involution `x <-> z-x` without choosing an
orientation:

`lambda_z^unord = (1/2)(pi_z)_# dx`, `pi_z(x)={{x,z-x}}`.

The orbit source is the sum of the two parent orders.  Arbitrary divergence-free
Fourier vectors then resolve at the event into all eight helical sectors and the
modal works reconstruct direct vector Navier--Stokes work exactly.

On each certified one-edge registration define

`dW = C_F T_e d lambda_unord`,
`dA = C_F A_e d lambda_unord`,
`dF = g_e dW = J_e c_e dA`.

`dW` is signed physical child-energy work. `dA` is a positive interaction-capacity
reference and is **not** causal.  Only after exact signed reconstruction is the
Hahn law formed; therefore `W_plus-W_minus=W` and `W_plus >= [W]_+`, with no false
identity between aggregate positive work and pointwise positive edge mass.

The actual block transfer deficit is computed internally as

`epsilon = 1 - F/(J_* A)`.

The coherent service-or-flat adapter accepts this ledger and never accepts a
caller-supplied transfer deficit.  If the deficit is low, the same continuum law
supplies the signed-good core.  At `eta=1e-4`, the certified normalized physical
child-work/capacity Radon--Nikodym ratio lies in `[50/53,53/50]`, so the existing
physical transfer-defect moat receives its genuine premise and
`E_phys[D] <= (106/25) epsilon`.

Positive nonforward work remains physical: it sits in `dW^+` while `J=0` and is
not relabeled signed-good.

Stress:
- arbitrary continuum fibers: `{out.samples}`
- near-extremal physical core blocks: `{out.core_blocks}`
- worst helical reconstruction residual: `{out.worst_helical_reconstruction_residual:.3e}`
- worst ordered/unordered quotient residual: `{out.worst_ordered_quotient_residual:.3e}`
- worst parent-swap residual: `{out.worst_parent_swap_residual:.3e}`
- worst direct/modal work relative residual: `{out.worst_direct_modal_work_relative_residual:.3e}`
- worst direct/modal progress relative residual: `{out.worst_direct_modal_progress_relative_residual:.3e}`
- minimum Hahn margin over aggregate positive work: `{out.minimum_hahn_aggregate_margin:.3e}`
- minimum Hahn margin over fiber-positive work: `{out.minimum_hahn_fiber_margin:.3e}`
- worst polarization residual: `{out.worst_polarization_residual:.3e}`
- signed flux ratio range: `[{out.minimum_signed_flux_ratio:.9f},{out.maximum_signed_flux_ratio:.9f}]`
- random blocks carrying positive nonforward edge work: `{out.positive_nonforward_blocks}`
- minimum good-core capacity Markov margin: `{out.minimum_good_core_capacity_margin:.3e}`
- minimum good-core RN lower margin: `{out.minimum_good_core_rn_lower_margin:.3e}`
- minimum good-core RN upper margin: `{out.minimum_good_core_rn_upper_margin:.3e}`
- minimum nonforward transfer-cost margin: `{out.minimum_service_transfer_cost_margin:.3e}`
- exact-extremal coherent-flat blocks: `{out.low_deficit_flat_blocks}`

This theorem closes the continuum **measure-registration** seam.  It does not say
every generic HH block is low deficit, and it does not terminate nonforward or
high-deficit physical events.  No Navier--Stokes global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
