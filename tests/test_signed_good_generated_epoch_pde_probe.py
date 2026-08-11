from src.signed_good_generated_epoch_pde_probe import (
    PARENT_P,
    STATUS,
    simulate_signed_good_triad_galerkin,
)


def test_actual_navier_stokes_triad_supplies_the_signed_good_physical_work_law():
    out = simulate_signed_good_triad_galerkin(
        resolution=24,
        steps=48,
        viscosity=0.02,
        amplitude=64.0,
        scaled_lifetime=0.05,
    )
    assert STATUS.startswith("DEALIASED_FOURIER_GALERKIN_NAVIER_STOKES")
    assert PARENT_P == (3, 2, 0)
    assert 3.0 / 5.0 < out.signed_good_parent_child_ratio < 5.0 / 8.0
    assert out.route_branch == "physical_high_high_transfer_generation"
    assert out.actual_positive_hh_work >= out.energy_gate_hh_work_lower * (1.0 - 2.0e-6)
    assert out.selected_heavy_half_fraction >= 0.5
    assert out.residual_positive_work_to_final_child_energy < 0.2
    assert out.low_pass_strain_action_upper < 1.0e-10
    assert out.final_child_energy > out.initial_child_energy
    assert out.epoch_hits_initial_boundary is True
