from src.pure_uv_true_upward_natural_window_pde_probe import run_probe


def test_actual_galerkin_pure_uv_law_is_first_shell_single_charged_and_representation_native():
    out = run_probe(
        main_resolutions=(24,),
        deep_resolutions=(20,),
        main_steps=16,
        deep_steps=8,
        main_duration=0.00025,
        deep_duration=0.00008,
        closed_tail_resolutions=(20,),
        closed_tail_steps=8,
        closed_tail_duration=0.0001,
    )
    assert out.maximum_common_work_representation_relative_spread <= 5e-8
    assert out.minimum_cutoff_support_margin > 0.0
    assert out.maximum_scale_probability_residual == 0.0
    assert out.coexistence_with_resolved_contact_observed
    assert all(o.recipient_shell_index == 1 and o.p_scale == 1.0 for o in out.observations)
    assert out.upstream_probe.maximum_pure_support_work_representation_native_residual <= 5e-8
