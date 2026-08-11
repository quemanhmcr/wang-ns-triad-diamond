from src.canonical_positive_edge_work_routing_pde_probe import STATUS, run_probe


def test_actual_galerkin_ns_routes_the_same_canonical_positive_edge_law():
    out = run_probe(resolution=20, steps=16, duration=0.0025, snapshot_count=3, tau=0.1)
    assert out.status == STATUS
    assert out.positive_work_snapshots > 0
    assert out.bad_work_snapshots > 0
    assert out.worst_signed_ns_reconstruction_relative < 3e-8
    assert out.worst_positive_mass_reconstruction_relative < 3e-9
    assert out.worst_hard_pushforward_relative < 3e-9
    assert 0.0 <= out.maximum_exact_role_mixed_good_fraction <= 1.0
    assert 0.0 <= out.maximum_coarsened_mixed_good_fraction <= 1.0 + 1e-10
    assert out.stage_zero_first_time_failures == 0
    assert out.geometry_good_marking_promotions == 0
