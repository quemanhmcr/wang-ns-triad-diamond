"""Adversarial native-scale and provenance checks for continuum helical edge measures.

These tests deliberately challenge theorem boundaries rather than numerical
fixtures.  They are intended to run in GitHub Actions on the independent audit
lane; local development should remain static/basic only.
"""

from dataclasses import replace
import math

import numpy as np
import pytest

from src.continuum_helical_edge_measure_registration import (
    continuum_edge_measure_ledger,
    edge_measure_to_service_or_flat,
    register_continuum_triad_fiber,
    signed_good_core_physical_law,
)
from src.helical import coupling_g, helical_basis


def _edge(scale: float = 1.0):
    x = scale * np.asarray((0.52, 0.31, 0.0))
    y = scale * np.asarray((0.52, -0.31, 0.0))
    return x, y, x + y


def _fiber(*, wave_scale: float = 1.0, amplitude: float = 1.0, quotient_mass: float = 1.0):
    x, y, z = _edge(wave_scale)
    sx, sy, sz = 1, -1, 1
    g = coupling_g(x, y, -z, sx, sy, sz)
    signed_frequency = sx * np.linalg.norm(x) - sy * np.linalg.norm(y)
    target_sign = 1.0 if signed_frequency >= 0.0 else -1.0
    az = target_sign * np.exp(-1j * np.angle(g))
    return register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=amplitude * helical_basis(x, sx),
        uy=amplitude * helical_basis(y, sy),
        uz=amplitude * az * helical_basis(z, sz),
        quotient_measure_mass=quotient_mass,
    )


def _service(ledger):
    return edge_measure_to_service_or_flat(
        ledger,
        tau=0.01,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=True,
        scaled_lifetime=1.0,
    )


def test_continuum_registration_is_covariant_under_tiny_uniform_wavevector_dilation():
    base = continuum_edge_measure_ledger((_fiber(),))
    scale = 1.0e-120
    tiny = continuum_edge_measure_ledger((_fiber(wave_scale=scale),))

    assert tiny.block_transfer_deficit == pytest.approx(base.block_transfer_deficit, abs=2.0e-11)
    assert tiny.normalized_signed_flux == pytest.approx(base.normalized_signed_flux, abs=2.0e-11)
    assert tiny.signed_direct_work / base.signed_direct_work == pytest.approx(scale, rel=3.0e-11)
    assert tiny.capacity_mass / base.capacity_mass == pytest.approx(scale, rel=3.0e-11)


def test_continuum_registration_rejects_foreign_tiny_child_at_native_scale():
    x, y, z = _edge(1.0e-12)
    foreign = z + np.asarray((1.5e-12, 0.0, 0.0))
    with pytest.raises(ValueError, match=r"z=x\+y|triad|parent pair"):
        register_continuum_triad_fiber(
            x=x,
            y=y,
            z=foreign,
            ux=helical_basis(x, 1),
            uy=helical_basis(y, -1),
            uz=helical_basis(foreign, 1),
            quotient_measure_mass=1.0,
        )


def test_tiny_nondivergencefree_vector_cannot_hide_behind_absolute_unit_floor():
    x, y, z = _edge(1.0e-12)
    # The parallel component is order-one relative to this tiny Fourier vector.
    # Multiplying the modal amplitude by the same native scale must not make a
    # physically longitudinal mode look divergence free merely because |k.u|<1.
    bad_ux = 1.0e-12 * x / np.linalg.norm(x)
    with pytest.raises(ValueError, match="divergence free"):
        register_continuum_triad_fiber(
            x=x,
            y=y,
            z=z,
            ux=bad_ux.astype(complex),
            uy=1.0e-12 * helical_basis(y, -1),
            uz=1.0e-12 * helical_basis(z, 1),
            quotient_measure_mass=1.0,
        )


def test_fiber_certificate_cannot_rebind_quotient_mass_without_rebinding_atoms():
    fiber = _fiber(amplitude=1.0e-6, quotient_mass=1.0)
    forged = replace(fiber, quotient_measure_mass=2.0)
    with pytest.raises((ValueError, AssertionError), match="quotient|measure|mass|reconstruct|provenance"):
        continuum_edge_measure_ledger((forged,))


def test_fiber_certificate_cannot_forge_native_direct_work_density():
    fiber = _fiber(amplitude=1.0e-7)
    native = abs(fiber.modal_signed_work_density)
    assert native > 0.0
    forged = replace(
        fiber,
        direct_signed_work_density=fiber.direct_signed_work_density + 0.25 * native,
        signed_work_reconstruction_residual=0.0,
    )
    with pytest.raises((ValueError, AssertionError), match="work|reconstruct|provenance|fiber"):
        continuum_edge_measure_ledger((forged,))


def test_high_deficit_ledger_cannot_forge_a_low_deficit_good_core_certificate():
    x = np.asarray((1.0, 0.0, 0.0))
    y = np.asarray((-0.8, 0.6, 0.0))
    z = x + y
    physical = register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=helical_basis(x, 1),
        uy=helical_basis(y, -1),
        uz=1j * helical_basis(z, 1),
        quotient_measure_mass=1.0,
    )
    ledger = continuum_edge_measure_ledger((physical,))
    assert ledger.block_transfer_deficit > 0.5

    forged = replace(
        ledger,
        block_transfer_deficit=1.0e-8,
        good_core_capacity_fraction=1.0,
        good_core_physical_to_capacity_rn_min=1.0,
        good_core_physical_to_capacity_rn_max=1.0,
    )
    with pytest.raises((ValueError, AssertionError), match="ledger|deficit|physical|provenance|replay"):
        signed_good_core_physical_law(forged)


def test_service_adapter_replays_physical_ledger_instead_of_trusting_forged_deficit():
    x = np.asarray((1.0, 0.0, 0.0))
    y = np.asarray((-0.8, 0.6, 0.0))
    z = x + y
    physical = register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=helical_basis(x, 1),
        uy=helical_basis(y, -1),
        uz=1j * helical_basis(z, 1),
        quotient_measure_mass=1.0,
    )
    ledger = continuum_edge_measure_ledger((physical,))
    assert ledger.block_transfer_deficit > 0.5
    forged = replace(ledger, block_transfer_deficit=0.0)

    with pytest.raises((ValueError, AssertionError), match="ledger|deficit|physical|provenance|replay"):
        _service(forged)


def test_nonzero_continuum_capacity_must_not_silently_underflow_to_zero():
    with pytest.raises(ValueError, match="capacity|native.*range|underflow"):
        _fiber(amplitude=1.0e-120)
