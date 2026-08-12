from src.cyclic_helical_triad_donor_kernel_pde_probe import STATUS, run_probe


def test_evolved_actual_galerkin_ns_preserves_cyclic_donor_kernel_and_signed_good_side_recipient():
    out = run_probe(resolutions=(24,), cutoff=7, steps=16, duration=0.0005, snapshot_count=3)
    assert out.status == STATUS
    assert out.common_cutoff == 7
    assert len(out.runs) == 1
    run = out.runs[0]
    assert run.worst_cyclic_energy_conservation_relative < 4e-8
    assert run.worst_cyclic_coupling_native_residual < 5e-12
    assert run.worst_measure_donor_marginal_relative < 4e-8
    assert run.worst_measure_recipient_marginal_relative < 4e-8
    assert run.global_energy_balance_relative_residual < 5e-5
    assert run.maximum_global_nonlinear_work_relative_rate < 5e-10
    assert run.maximum_divergence_relative_to_initial_l2 < 5e-11
    assert run.initial_signed_good_efficiency > 1.0 - 1.0e-4
    assert 0.3 < run.initial_side_to_child_ratio < 1.0 / 3.0
    assert 0.75 < run.initial_child_to_donor_ratio < 10.0 / 13.0
    assert 3.0 / 13.0 < run.initial_side_to_donor_ratio < 0.25
    assert run.initial_side_forward_ratio < 1.0
    assert run.initial_side_geometric_multiplier == 0.0
    assert run.initial_donor_count == 1
    assert run.initial_recipient_count == 2
