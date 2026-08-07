from src.cycle_hodge import (
    butterfly_rigidity_certificate,
    flat_butterfly_motif,
    nonflat_reuse_motif,
    optimal_geometry,
    triad_residual_identity,
    weighted_hodge_energy,
)


def test_triad_residual_identity():
    row = triad_residual_identity(0.2, -0.1, 0.8, 0.3)
    assert abs(row["lhs"] - row["rhs"]) < 1e-12


def test_nonflat_cycle_has_exact_energy():
    gamma = 0.4928152853421352
    row = weighted_hodge_energy(nonflat_reuse_motif(), gamma)
    assert row["cycle_rank"] == 1
    assert abs(row["energy"] - gamma * gamma / 5.0) < 1e-10
    assert abs(row["energy"] - row["dual_energy"]) < 1e-10
    assert row["dual_stationarity"] < 1e-10
    assert not row["flat"]


def test_equal_transfer_nonflat_motif_has_certified_positive_block_cost():
    gamma = 0.4928152853421352
    # Three triads have normalized transfer weights 1/3; each contributes two
    # Hodge arcs with the same inherited conductance.
    row = weighted_hodge_energy(nonflat_reuse_motif(), gamma, weights=[1.0 / 3.0] * 6)
    assert abs(row["energy"] - gamma * gamma / 15.0) < 1e-11
    block_log_cost = 0.5 * row["energy"]
    assert abs(block_log_cost - gamma * gamma / 30.0) < 1e-11
    assert block_log_cost > 0.008


def test_butterfly_is_scale_flat_despite_cycle():
    gamma = 0.4928152853421352
    row = weighted_hodge_energy(flat_butterfly_motif(), gamma)
    assert row["cycle_rank"] == 1
    assert row["energy"] < 1e-20
    levels = row["levels"]
    assert abs((levels["m"] - levels["a"]) - 1.0) < 1e-10
    assert abs((levels["n"] - levels["a"]) - 1.0) < 1e-10
    assert abs((levels["d"] - levels["m"]) - 1.0) < 1e-10


def test_exact_butterfly_returns_central_direction():
    geom = optimal_geometry(0.6109041015867660)
    row = butterfly_rigidity_certificate(geom["c"])
    assert abs(row["b_dot_c"] - row["required_b_dot_c"]) < 1e-12
    assert abs(row["m_dot_n"] - geom["c"]) < 1e-12
    assert row["d_error"] < 1e-12
    assert abs(row["internal_angle_gap"] - 0.5 * geom["theta"]) < 1e-12


def test_planar_exact_erosion():
    from src.cycle_hodge import planar_erosion_step
    theta = 1.0
    # Pairs (0,1), (1,2), (2,3) give children 0.5,1.5,2.5.
    row = planar_erosion_step([0.0, 1.0, 2.0, 3.0], theta)
    assert row["actual_exact_next_diameter"] == 2.0
    assert row["next_diameter_bound"] == 2.0


def test_planar_fresh_span_telescopes():
    from src.cycle_hodge import planar_fresh_span_lower_bound
    theta = 1.2
    # Maintaining the same diameter for 5 exact generations costs at least 5 theta.
    cost = planar_fresh_span_lower_bound(3.0, 3.0, 5, theta)
    assert abs(cost - 6.0) < 1e-12


def test_near_planar_erosion_bound():
    from src.cycle_hodge import planar_erosion_step
    row = planar_erosion_step([0.0, 3.0], 1.0, fresh_angles=[-0.2, 3.4], pair_tolerance=0.1, midpoint_error=0.02)
    assert abs(row["fresh_expansion"] - 0.6) < 1e-12
    assert abs(row["next_diameter_bound"] - (3.6 - 0.9 + 0.04)) < 1e-12
