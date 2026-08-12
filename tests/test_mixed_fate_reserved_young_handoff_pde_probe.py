from src.mixed_fate_reserved_young_handoff_pde_probe import STATUS, run_probe


def test_evolved_actual_galerkin_ns_exposes_mixed_good_bad_negative_hard_cell_without_rehahn():
    out = run_probe(resolution=24, steps=16, duration=0.0005, snapshot_count=2, tau=0.1)
    assert out.status == STATUS
    assert out.mixed_fate_snapshots > 0
    assert out.snapshots_with_good_work > 0
    assert out.snapshots_with_bad_positive_work > 0
    assert out.snapshots_with_negative_work > 0
    assert out.initial_maximum_good_signed_efficiency > 1.0 - 1.0e-4
    assert out.initial_mixed_good_work > 0.0
    assert out.worst_signed_ns_reconstruction_relative < 4e-8
    assert out.worst_hahn_pushforward_identity_relative < 4e-8
    assert out.worst_reservation_identity_relative < 4e-8
    assert out.stage_zero_first_time_failures == 0
    assert out.geometry_good_marking_promotions == 0
