import math

import pytest

from src.signed_good_generated_epoch_two_layer_pde_probe import (
    MIDDLE_CHILD,
    STATUS,
    TOP_CHILD,
    simulate_signed_good_two_layer_galerkin,
)


def test_naive_two_triad_geometry_is_rejected_when_actual_ns_work_routes_elsewhere():
    assert STATUS.startswith("DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES")
    assert MIDDLE_CHILD == (6, 0, 0)
    assert TOP_CHILD == (8, 4, 4)
    with pytest.raises(AssertionError, match="top exact triad.*classified_residual_physical_work"):
        simulate_signed_good_two_layer_galerkin(
            resolution=28,
            steps=64,
            viscosity=0.02,
            amplitude=96.0,
            scaled_lifetime=0.05,
            middle_seed_weight=5.0e-4,
            partner_polarization_angle=0.0,
        )


def test_two_exact_signed_good_triads_telescope_on_one_navier_stokes_orbit():
    out = simulate_signed_good_two_layer_galerkin(
        resolution=28,
        steps=64,
        viscosity=0.02,
        amplitude=96.0,
        scaled_lifetime=0.05,
        middle_seed_weight=5.0e-4,
        partner_polarization_angle=math.pi / 4.0,
    )
    assert out.partner_polarization_angle == pytest.approx(math.pi / 4.0)
    assert 3.0 / 5.0 < out.middle_parent_child_ratio < 5.0 / 8.0
    assert 3.0 / 5.0 < out.top_parent_child_ratio < 5.0 / 8.0
    assert out.middle_initial_energy_fraction < 0.2
    assert out.top_initial_energy_fraction < 0.2
    assert out.middle_actual_positive_hh_work >= out.middle_energy_gate_lower * (1.0 - 2.0e-5)
    assert out.top_actual_positive_hh_work >= out.top_energy_gate_lower * (1.0 - 2.0e-5)
    assert out.middle_heavy_half_fraction >= 0.5
    assert out.top_heavy_half_fraction >= 0.5
    assert out.first_common_reference_time > 0.0
    assert out.last_common_reference_time <= 0.0
    assert out.epoch_layer_count == 2
    assert out.epoch_hits_initial_boundary is True
