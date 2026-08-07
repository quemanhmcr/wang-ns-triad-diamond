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
    gamma = 0.4928152849469459
    row = weighted_hodge_energy(nonflat_reuse_motif(), gamma)
    assert row["cycle_rank"] == 1
    assert abs(row["energy"] - gamma * gamma / 5.0) < 1e-10
    assert abs(row["energy"] - row["dual_energy"]) < 1e-10
    assert row["dual_stationarity"] < 1e-10
    assert not row["flat"]


def test_butterfly_is_scale_flat_despite_cycle():
    gamma = 0.4928152849469459
    row = weighted_hodge_energy(flat_butterfly_motif(), gamma)
    assert row["cycle_rank"] == 1
    assert row["energy"] < 1e-20
    levels = row["levels"]
    assert abs((levels["m"] - levels["a"]) - 1.0) < 1e-10
    assert abs((levels["n"] - levels["a"]) - 1.0) < 1e-10
    assert abs((levels["d"] - levels["m"]) - 1.0) < 1e-10


def test_exact_butterfly_returns_central_direction():
    geom = optimal_geometry(0.6109041018281888)
    row = butterfly_rigidity_certificate(geom["c"])
    assert abs(row["b_dot_c"] - row["required_b_dot_c"]) < 1e-12
    assert abs(row["m_dot_n"] - geom["c"]) < 1e-12
    assert row["d_error"] < 1e-12
    assert abs(row["internal_angle_gap"] - 0.5 * geom["theta"]) < 1e-12
