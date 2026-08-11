from src.continuum_helical_edge_measure_pde_probe import (
    STATUS,
    simulate_continuum_edge_measure_on_galerkin_ns,
)


def test_actual_galerkin_ns_convolution_is_the_registered_unordered_edge_measure():
    out = simulate_continuum_edge_measure_on_galerkin_ns(
        resolution=20,
        steps=24,
        viscosity=0.03,
        amplitude=4.0,
        duration=0.006,
        snapshot_count=3,
    )
    assert STATUS.startswith("EVOLVED_DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES")
    assert out.unordered_pairs > 0
    assert out.modal_edges == 8 * out.unordered_pairs
    assert out.nonzero_source_snapshots == out.snapshots
    assert out.worst_ordered_source_residual < 3.0e-8
    assert out.worst_unordered_source_residual < 3.0e-8
    assert out.worst_ordered_unordered_residual < 3.0e-8
    assert out.worst_signed_work_residual < 3.0e-8
    assert out.worst_signed_modal_work_residual < 3.0e-8
    assert out.worst_progress_residual < 3.0e-8
    assert out.worst_hahn_residual < 3.0e-8
    assert out.global_energy_balance_relative_residual < 5.0e-5
    assert out.maximum_global_nonlinear_work_relative_rate < 5.0e-10
    assert out.maximum_divergence_relative_to_initial_l2 < 5.0e-11


def test_same_galerkin_system_is_invariant_under_fft_grid_embedding():
    from src.continuum_helical_edge_measure_pde_probe import run_probe

    out = run_probe(
        (20, 24),
        steps=16,
        viscosity=0.03,
        amplitude=4.0,
        duration=0.003,
        snapshot_count=3,
    )
    assert out.common_cutoff == 5
    assert out.common_final_child_energy_resolution_spread < 5.0e-8
    assert out.common_integrated_child_work_resolution_spread < 5.0e-8
    # Native truncations are different PDEs; their spread is recorded, not gated.
    assert out.native_final_child_energy_resolution_spread >= 0.0
    assert out.native_integrated_child_work_resolution_spread >= 0.0
