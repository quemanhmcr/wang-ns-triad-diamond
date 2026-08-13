from src.resolved_contact_native_binding_pde_probe import run_probe


def test_actual_galerkin_NS_keeps_canonical_cause_fixed_while_cutoff_repartitions_Vh_and_hh():
    out = run_probe(
        resolutions=(24,),
        boundary_cutoff=2,
        interior_cutoff=3,
        steps=8,
        duration=5.0e-5,
        viscosity=0.05,
        amplitude=0.8,
    )
    assert out.maximum_representation_spread < 5e-10
    boundary = out.boundary_contact_runs[0]
    interior = out.interior_contact_runs[0]
    assert boundary.maximum_transition_profile_q < 5e-12
    assert boundary.maximum_plateau_profile_q < 5e-12
    assert 0.0 < interior.minimum_transition_profile_q < 1.0
    assert interior.maximum_transition_profile_q < 1.0
    assert interior.minimum_plateau_profile_q > 1.0 - 5e-12
    assert interior.worst_cutoff_repartition_gauge_residual < 3e-12
    assert interior.worst_ks_signed_identity_native_residual < 5e-10
    assert interior.worst_canonical_ks_positive_cover_defect < 5e-10
