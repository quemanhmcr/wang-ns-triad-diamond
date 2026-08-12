import pytest

from src.cyclic_hard_cell_single_charge_quotient_pde_probe import STATUS, run_probe


def test_actual_ns_positive_phase_exposes_one_donor_two_recipient_good_bad_single_charge():
    out = run_probe(
        resolutions=(24,), cutoff=7, steps=16, duration=0.0005, snapshot_count=3, phase_sign=1
    )
    assert out.status == STATUS
    run = out.runs[0]
    assert run.initial_donor_count == 1
    assert run.initial_recipient_count == 2
    assert run.initial_good_recipient_mass > 0.0
    assert run.initial_bad_recipient_mass > 0.0
    assert run.initial_coarse_self_loop_fraction == pytest.approx(1.0)
    assert max(
        run.worst_balance_native_residual,
        run.worst_donor_marginal_native_residual,
        run.worst_recipient_marginal_native_residual,
        run.worst_fate_partition_native_residual,
        run.worst_restricted_pushforward_native_residual,
        run.worst_coarse_total_native_residual,
        run.worst_coarse_fate_native_residual,
    ) < 5e-8


def test_actual_ns_phase_reversal_exposes_two_donors_recombining_to_one_recipient_charge():
    out = run_probe(
        resolutions=(24,), cutoff=7, steps=16, duration=0.0005, snapshot_count=3, phase_sign=-1
    )
    run = out.runs[0]
    assert run.initial_donor_count == 2
    assert run.initial_recipient_count == 1
    assert run.initial_overlapping_recipient_charge_count == 1
    assert run.initial_coarse_self_loop_fraction == pytest.approx(1.0)


def test_same_cutoff_actual_ns_hard_cell_charges_are_fft_representation_invariant():
    out = run_probe(
        resolutions=(24, 28), cutoff=7, steps=16, duration=0.0005, snapshot_count=3, phase_sign=1
    )
    assert out.maximum_total_positive_work_representation_native_residual < 5e-8
    assert out.maximum_good_recipient_work_representation_native_residual < 5e-8
    assert out.maximum_bad_recipient_work_representation_native_residual < 5e-8
    assert {run.cutoff for run in out.runs} == {7}


def test_invalid_phase_sign_is_rejected_before_physical_audit():
    with pytest.raises(ValueError, match="phase_sign"):
        run_probe(resolutions=(24,), steps=16, snapshot_count=3, phase_sign=0)
