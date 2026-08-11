from src.high_strain_descending_epoch_pde_probe import (
    simulate_high_strain_epoch_galerkin,
)


def test_actual_galerkin_high_strain_event_enters_physical_descending_epoch_route():
    out = simulate_high_strain_epoch_galerkin(
        resolution=12,
        steps=64,
        viscosity=0.05,
        amplitude=256.0,
        child_frequency=16.0,
        scaled_lifetime=1.0,
    )
    assert out.root_strain_action >= out.high_strain_action_threshold
    assert out.root_normalized_resolved_dissipation >= out.high_strain_dissipation_lower
    assert out.collision_relative_margin >= -5.0e-3
    assert out.global_reservoir_relative_margin >= -2.0e-10
    assert out.retained_critical_fraction >= 0.5 - 5.0e-3
    assert out.selected_shell_critical_dissipation > 0.0
    assert out.renewal_frequency / out.child_frequency <= 3.0 / 16.0 * (1.0 + 1.0e-11)
    assert out.descendant_high_strain_excluded_by_spectral_gap is True
    assert out.maximum_descendant_resolved_gradient_relative_to_root <= 1.0e-24
    assert out.observed_epoch_steps == 1
