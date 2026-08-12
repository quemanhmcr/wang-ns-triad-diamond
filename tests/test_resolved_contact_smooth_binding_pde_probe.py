from src.resolved_contact_smooth_binding_pde_probe import run_probe


def test_actual_galerkin_ns_referee_keeps_borderline_and_strict_deep_physics_distinct():
    out = run_probe(resolutions=(24,), steps=8, duration=5.0e-5, amplitude=0.8)
    assert len(out.borderline_runs) == 1
    assert len(out.strict_deep_runs) == 1
    border = out.borderline_runs[0]
    strict = out.strict_deep_runs[0]
    assert border.target_snapshots == border.total_snapshots
    assert border.maximum_mixed_fraction <= 5e-12
    assert border.minimum_transition_hh_fraction >= 1.0 - 5e-12
    assert strict.target_snapshots == strict.total_snapshots
    assert strict.minimum_mixed_fraction >= 1.0 - 5e-12
    assert strict.maximum_transition_hh_fraction <= 5e-12
    assert strict.worst_ks_to_mixed_edge_native_residual <= 5e-10
    assert strict.worst_ks_signed_identity_native_residual <= 5e-10
    assert strict.worst_ks_skew_pair_residual <= 5e-10
    assert strict.worst_ks_strain_pair_residual <= 5e-10
    assert strict.worst_canonical_ks_mass_residual <= 5e-10
    assert out.maximum_representation_spread <= 5e-10
