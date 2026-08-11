from src.full_natural_service_corridor_pde_probe import STATUS, run_physical_pde_probe


def test_dealiased_navier_stokes_corridor_keeps_q2_cover_and_service_on_one_interval():
    result = run_physical_pde_probe(resolutions=(12, 16), steps=24, duration=0.001)
    assert result.status == STATUS
    assert result.scaled_lifetime == result.duration * result.carrier_frequency**2
    assert result.final_carrier_energy_resolution_spread < 2.0e-2
    assert all(run.maximum_absolute_carrier_nonlinear_work > 1.0e-8 for run in result.runs)
    assert all(run.integrated_same_interval_bounded_service > 0 for run in result.runs)
