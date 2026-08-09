import math

import numpy as np
import pytest

from src.critical_annular_carrier_service_reentry import (
    persistent_carrier_critical_mass_lower,
    uniform_bounded_square_service_lower,
)
from src.critical_shell_service_reentry import (
    critical_shell_backward_first_hit,
    critical_shell_bounded_service_lower,
    critical_shell_natural_outcome,
    critical_shell_persistent_carrier_mass_lower,
    critical_shell_survivor_coefficient_mass_lower,
    critical_shell_terminal_mass_lower,
    dissipation_supplier_retained_fraction_lower,
    dissipation_supplier_shell_mass_threshold,
    fresh_dominant_service_shell_mass_lower,
    theorem_certificate,
    two_cell_cluster_to_whole_shell_mass_lower,
)
from src.high_strain_resolved_ancestor import high_strain_ancestor_mass_threshold
from src.nn_critical_heat_carrier_seed import LOW_STRAIN_ACTION, renewal_scale


def test_high_strain_is_exact_specialization_of_generic_shell_theorem():
    for c, nu in ((0.5, 0.0), (1.0, 1.0), (2.0, 1.5)):
        mu = high_strain_ancestor_mass_threshold(c)
        assert math.isclose(
            critical_shell_persistent_carrier_mass_lower(mu, c, nu),
            persistent_carrier_critical_mass_lower(c, nu),
            rel_tol=2e-14,
        )
        assert math.isclose(
            critical_shell_bounded_service_lower(mu, c, nu),
            uniform_bounded_square_service_lower(c, nu),
            rel_tol=2e-14,
        )
        assert math.isclose(
            critical_shell_survivor_coefficient_mass_lower(mu),
            math.pi**2 / (50.0 * c * c),
            rel_tol=2e-14,
        )


def test_generic_D0_supplier_retains_half_actual_resolved_dissipation():
    for c, D0, D in ((1.0, 2.0, 2.0), (0.4, 0.3, 1.2), (2.0, 5.0, 20.0)):
        assert math.isclose(dissipation_supplier_shell_mass_threshold(D0, c), D0 / c)
        assert dissipation_supplier_retained_fraction_lower(D, D0, c) >= 0.5 - 1e-14


def test_dominant_fresh_two_cell_cluster_supplies_whole_shell_without_packet_floor():
    Y = 8.0
    pair_mass = (1.0 / 4.0) * Y / 8.0
    whole = two_cell_cluster_to_whole_shell_mass_lower(pair_mass)
    assert whole == Y / 64.0
    assert fresh_dominant_service_shell_mass_lower(Y) == Y / 64.0


def test_generic_corridor_has_only_three_native_monitors_and_records_horizon():
    amp = 2.0
    ell = np.linspace(0.0, 1.0, 5)
    hit = critical_shell_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1 / 60, 5),
        residual_impulse_abs=np.linspace(0.0, 0.1 * amp, 5),
        hh_impulse_abs=np.linspace(0.0, 0.2 * amp, 5),
    )
    assert hit["first_elapsed"] is None
    assert hit["observed_elapsed_end"] == 1.0
    assert set(hit["individual_debuts"]) == {
        "high_strain_critical_dissipation",
        "classified_role_interface_impulse",
        "hh_regeneration_impulse",
    }


def test_incomplete_monitor_horizon_cannot_certify_survivor_or_boundary():
    c = 1.0
    M = 4.0
    A = renewal_scale(M)
    T = c / A**2
    mu0 = 1.0
    amp = math.sqrt(critical_shell_terminal_mass_lower(mu0) / A)
    short = np.linspace(0.0, T / 2.0, 4)
    hit = critical_shell_backward_first_hit(
        short,
        terminal_amplitude=amp,
        strain_action=np.zeros(4),
        residual_impulse_abs=np.zeros(4),
        hh_impulse_abs=np.zeros(4),
    )
    with pytest.raises(ValueError, match="do not cover"):
        critical_shell_natural_outcome(
            event_time=2 * T,
            renewal_frequency=A,
            shell_critical_mass_lower=mu0,
            scaled_lifetime=c,
            viscosity=1.0,
            terminal_coefficient=amp,
            endpoint_coefficient=amp,
            hh_impulse=0j,
            residual_interface_impulse=0j,
            first_hit=hit,
        )


def test_first_hit_threshold_amplitude_must_match_registered_terminal_coefficient():
    c = 1.0
    M = 4.0
    A = renewal_scale(M)
    T = c / A**2
    mu0 = 1.0
    amp = math.sqrt(critical_shell_terminal_mass_lower(mu0) / A)
    ell = np.linspace(0.0, T, 4)
    hit = critical_shell_backward_first_hit(
        ell,
        terminal_amplitude=2.0 * amp,
        strain_action=np.zeros(4),
        residual_impulse_abs=np.zeros(4),
        hh_impulse_abs=np.zeros(4),
    )
    with pytest.raises(ValueError, match="thresholds do not match"):
        critical_shell_natural_outcome(
            event_time=2 * T,
            renewal_frequency=A,
            shell_critical_mass_lower=mu0,
            scaled_lifetime=c,
            viscosity=1.0,
            terminal_coefficient=amp,
            endpoint_coefficient=amp,
            hh_impulse=0j,
            residual_interface_impulse=0j,
            first_hit=hit,
        )


def test_full_generic_shell_corridor_creates_positive_own_scale_service_before_materiality():
    c = 1.0
    M = 4.0
    A = renewal_scale(M)
    T = c / A**2
    mu0 = 2.0
    amp = math.sqrt(1.2 * critical_shell_terminal_mass_lower(mu0) / A)
    ell = np.linspace(0.0, T, 5)
    hit = critical_shell_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1 / 60, 5),
        residual_impulse_abs=np.linspace(0.0, 0.1 * amp, 5),
        hh_impulse_abs=np.linspace(0.0, 0.2 * amp, 5),
    )
    ir = 0.1 * amp
    ih = 0.2j * amp
    out = critical_shell_natural_outcome(
        event_time=2 * T,
        renewal_frequency=A,
        shell_critical_mass_lower=mu0,
        scaled_lifetime=c,
        viscosity=1.0,
        terminal_coefficient=amp,
        endpoint_coefficient=amp - ir - ih,
        hh_impulse=ih,
        residual_interface_impulse=ir,
        first_hit=hit,
    )
    assert out["classification"] == "full_natural_own_scale_service"
    assert out["clean_retained_coefficient_mass_lower"] == critical_shell_survivor_coefficient_mass_lower(mu0)
    assert out["uniform_square_service_lower"] > 0
    assert out["materiality_assigned"] == "only_after_service_via_actual_Moyal_endpoints"


def test_certificate_keeps_DV_sampling_noncausal_and_scale_progress_supplier_specific():
    cert = theorem_certificate()
    assert "diagnostic" in cert["dissipation_supplier"]
    assert "never HH causal probabilities" in cert["dissipation_supplier"]
    assert "supplier-specific signed-good progress" in cert["scale_scope"]
    assert "pressure/source routing" in cert["scope"]
