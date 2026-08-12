from src.hard_tail_true_upward_supply_pde_probe import run_probe


def test_actual_ns_keeps_full_tail_reading_separate_from_pure_and_deep_selected_support():
    out=run_probe(
        main_resolutions=(24,),
        deep_resolutions=(20,),
        main_cutoff=7,
        deep_cutoff=2,
        main_steps=16,
        deep_steps=8,
        viscosity=0.03,
        amplitude=0.8,
        main_duration=0.00025,
        deep_duration=0.00008,
        closed_tail_resolutions=(20,),
        closed_tail_steps=8,
        closed_tail_duration=0.00008,
    )
    assert out.radial_tail_probe.runs[0].tail_interval_balance_native_residual < 5e-5
    pure=out.selected_pure_support[0]
    assert pure.upward_work>0.0
    assert pure.pure_uv_work>0.0
    assert pure.all_pure_atoms_first_shell
    assert pure.all_energy_donors_are_interaction_parents
    assert pure.maximum_parent_to_shell_ratio<=1.5+5e-12
    deep=out.deep_contact_runs[0]
    assert deep.initial_deep_upward_work>0.0
    assert deep.initial_pure_uv_work==0.0
    assert deep.snapshots_with_deep_upward>=1
    assert deep.maximum_deep_donor_to_quarter_shell_excess<=5e-12
    full=out.closed_triad_tail_runs[0]
    assert full.tail_continuity_native_residual < 5e-5
    assert full.worst_full_cyclic_tail_reconstruction_native_residual < 5e-8
    assert full.worst_full_cyclic_boundary_divergence_native_residual < 5e-8
    assert full.inherited_owner or full.true_upward_owner
    assert out.resolved_contact_declared_interface_owner is False
