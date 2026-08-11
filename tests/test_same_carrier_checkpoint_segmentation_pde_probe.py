from src.same_carrier_checkpoint_segmentation_pde_probe import (
    STATUS,
    run_same_carrier_physical_pde_probe,
)


def test_actual_navier_stokes_path_is_invariant_under_same_carrier_checkpoint_cuts():
    result = run_same_carrier_physical_pde_probe(
        resolutions=(20,),
        steps=8,
        duration=0.015625,
        carrier_frequency=4.0,
        natural_window_count=4,
    )
    assert result.status == STATUS
    assert result.terminal_amplitude_resolution_spread == 0.0
    run = result.runs[0]
    assert run.maximum_segmentation_first_time_residual == 0.0
    assert run.maximum_absolute_hh_impulse > 0.0
    assert run.maximum_absolute_residual_impulse > 0.0
    assert run.maximum_absolute_imaginary_impulse > 0.0
    assert run.fixed_natural_windows_before_t0 == 4
    assert run.fixed_window_interior_zeno_possible is False
