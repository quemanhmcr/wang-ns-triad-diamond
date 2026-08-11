from src.helical_physical_edge_pde_probe import (
    STATUS,
    simulate_helical_edge_on_galerkin_ns,
)


def test_evolved_ns_parent_pair_is_exactly_the_sum_of_its_helical_edges():
    out = simulate_helical_edge_on_galerkin_ns(
        resolution=24,
        steps=32,
        viscosity=0.02,
        amplitude=64.0,
        scaled_lifetime=0.05,
    )
    assert STATUS.startswith("EVOLVED_DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES")
    assert out.positive_physical_work_samples > 0
    assert out.minimum_actual_source_norm > 0.0
    assert out.worst_source_registration_residual < 2.0e-9
    assert out.worst_full_pair_work_registration_residual < 2.0e-9
    assert out.worst_registered_upper_residual < 2.0e-9
    assert out.global_energy_balance_relative_residual < 2.0e-5
    assert out.maximum_normalized_multiplier <= 1.0 + 5.0e-10
