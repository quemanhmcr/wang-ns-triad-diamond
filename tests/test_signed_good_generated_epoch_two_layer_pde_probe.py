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
