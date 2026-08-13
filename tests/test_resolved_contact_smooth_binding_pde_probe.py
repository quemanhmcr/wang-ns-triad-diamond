from src.resolved_contact_smooth_binding_pde_probe import run_probe


def test_actual_galerkin_ns_referee_keeps_borderline_and_strict_deep_physics_distinct():
    out = run_probe(resolutions=(24,), steps=8, duration=5.0e-5, amplitude=0.8)
    assert len(out.boundary_contact_runs) == 1
    assert len(out.interior_transition_runs) == 1
    assert len(out.strict_deep_runs) == 1
    border = out.boundary_contact_runs[0]
    transition = out.interior_transition_runs[0]
    strict = out.strict_deep_runs[0]
    assert border.target_snapshots == border.total_snapshots
    assert border.maximum_mixed_fraction <= 5e-12
    assert border.minimum_transition_hh_fraction >= 1.0 - 5e-12
    assert transition.target_snapshots == transition.total_snapshots
    assert 1e-4 < transition.minimum_mixed_fraction < 1.0 - 1e-4
    assert 1e-4 < transition.minimum_transition_hh_fraction < 1.0 - 1e-4
    assert transition.worst_ks_to_mixed_edge_native_residual <= 5e-10
    assert transition.worst_ks_signed_identity_native_residual <= 5e-10
    assert transition.worst_existing_owner_work_split_native_residual <= 5e-10
    assert transition.worst_existing_skew_antisymmetry_native_residual <= 5e-10
    assert transition.worst_existing_strain_symmetry_native_residual <= 5e-10
    assert transition.worst_existing_skew_closure_balance_native_residual <= 5e-10
    assert transition.maximum_existing_skew_donor_path_length <= 1
    assert strict.target_snapshots == strict.total_snapshots
    assert strict.minimum_mixed_fraction >= 1.0 - 5e-12
    assert strict.maximum_transition_hh_fraction <= 5e-12
    assert strict.worst_ks_to_mixed_edge_native_residual <= 5e-10
    assert strict.worst_ks_signed_identity_native_residual <= 5e-10
    assert strict.worst_ks_skew_pair_residual <= 5e-10
    assert strict.worst_ks_strain_pair_residual <= 5e-10
    assert strict.worst_canonical_ks_mass_residual <= 5e-10
    assert strict.worst_existing_owner_work_split_native_residual <= 5e-10
    assert strict.worst_existing_skew_antisymmetry_native_residual <= 5e-10
    assert strict.worst_existing_strain_symmetry_native_residual <= 5e-10
    assert strict.worst_existing_skew_closure_balance_native_residual <= 5e-10
    assert strict.maximum_existing_skew_donor_path_length <= 1
    assert (transition.skew_bound_snapshots + transition.strain_bound_snapshots) > 0
    assert (strict.skew_bound_snapshots + strict.strain_bound_snapshots) > 0
    assert out.maximum_representation_spread <= 5e-10
