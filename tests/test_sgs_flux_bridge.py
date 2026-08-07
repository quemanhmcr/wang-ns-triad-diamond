from src.sgs_flux_bridge import divergence_free_random_field, sgs_flux_identity


def test_periodic_sgs_flux_equals_resolved_nonlinear_work():
    u = divergence_free_random_field(8, 31)
    row = sgs_flux_identity(u, 2.7)
    scale = max(abs(row["mean_sgs_flux"]), abs(row["mean_resolved_nonlinear_work"]), 1e-12)
    assert abs(row["sgs_vs_resolved_error"]) <= 2e-10 * scale + 2e-12
    assert abs(row["projection_work_error"]) <= 2e-10 * scale + 2e-12
    assert row["relative_divergence"] < 1e-12
