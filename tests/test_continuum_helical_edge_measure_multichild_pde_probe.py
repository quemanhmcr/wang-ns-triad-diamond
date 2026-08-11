from src.continuum_helical_edge_measure_multichild_pde_probe import (
    MULTI_CHILD_MODES,
    STATUS,
    simulate_multichild_edge_measure_on_galerkin_ns,
)


def test_multiple_children_share_one_actual_ns_orbit_and_one_outer_work_law():
    out = simulate_multichild_edge_measure_on_galerkin_ns(
        resolution=20,
        spectral_cutoff=5,
        child_modes=MULTI_CHILD_MODES,
        steps=12,
        duration=0.002,
        snapshot_count=2,
    )
    assert STATUS.startswith("EVOLVED_DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES")
    assert out.child_snapshot_registrations == 2 * len(MULTI_CHILD_MODES)
    assert out.minimum_unordered_pairs_per_child > 0
    assert out.minimum_modal_edges_per_child == 8 * out.minimum_unordered_pairs_per_child
    assert out.worst_child_source_residual < 3.0e-8
    assert out.worst_child_work_residual < 3.0e-8
    assert out.worst_child_progress_residual < 3.0e-8
    assert out.worst_child_hahn_residual < 3.0e-8
    assert out.worst_joint_work_residual < 3.0e-8
    assert out.worst_joint_modal_work_residual < 3.0e-8
    assert out.worst_joint_hahn_residual < 3.0e-8
    assert out.worst_joint_progress_residual < 3.0e-8
    assert out.global_energy_balance_relative_residual < 5.0e-5
    assert out.maximum_divergence_relative_to_initial_l2 < 5.0e-11
