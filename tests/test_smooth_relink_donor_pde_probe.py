from src.smooth_relink_donor_pde_probe import (
    STATUS,
    run_smooth_relink_physical_pde_probe,
)


def test_actual_navier_stokes_linearized_work_enters_finite_smooth_donor_route():
    result = run_smooth_relink_physical_pde_probe(
        resolutions=(12,),
        steps=4,
        duration=0.01,
        viscosity=0.05,
        carrier_frequency=4.0,
    )
    assert result.status == STATUS
    assert result.final_positive_relink_resolution_spread == 0.0
    run = result.runs[0]
    assert run.relink_owner_snapshots > 0
    assert run.maximum_positive_relink_work > 0.0
    assert run.maximum_shortest_donor_path_length <= 2
    assert run.master_route_failures == 0
