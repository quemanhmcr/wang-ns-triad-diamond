"""Adversarial audit guards for the full-natural service corridor quotient.

These tests exercise the representation boundaries which ordinary positive
fixtures do not see: parabolic time rescaling, physical scale provenance,
positive-service provenance, and arbitrarily small critical masses.
"""

import math

import numpy as np
import pytest

from src.critical_shell_service_reentry import (
    critical_shell_backward_first_hit,
    critical_shell_natural_outcome,
    critical_shell_terminal_mass_lower,
)
from src.full_natural_service_corridor_quotient import (
    FULL_NATURAL_SERVICE_WITNESS,
    FullNaturalServiceCorridor,
    endpoint_comparable_hard_shell_cover,
    endpoint_hard_shell_cover_from_full_natural_outcome,
    material_partition_is_same_corridor_measure,
    realized_endpoint_hard_shell_witnesses,
)
from src.high_strain_critical_carrier_reentry import (
    critical_seed_backward_first_hit,
    critical_seed_natural_outcome,
)
from src.nn_critical_heat_carrier_seed import renewal_carrier_critical_mass_lower, renewal_scale


def _actual_full_natural_outcome(*, parent_shell_frequency: float = 4.0) -> tuple[dict[str, object], float]:
    c = 1.0
    M = float(parent_shell_frequency)
    A = renewal_scale(M)
    T = c / A**2
    mu0 = 2.0
    amp = math.sqrt(1.2 * critical_shell_terminal_mass_lower(mu0) / A)
    elapsed = np.linspace(0.0, T, 5)
    first_hit = critical_shell_backward_first_hit(
        elapsed,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1.0 / 60.0, 5),
        residual_impulse_abs=np.linspace(0.0, 0.1 * amp, 5),
        hh_impulse_abs=np.linspace(0.0, 0.2 * amp, 5),
    )
    residual = 0.1 * amp
    hh = 0.2j * amp
    outcome = critical_shell_natural_outcome(
        event_time=2.0 * T,
        renewal_frequency=A,
        shell_critical_mass_lower=mu0,
        scaled_lifetime=c,
        viscosity=1.0,
        terminal_coefficient=amp,
        endpoint_coefficient=amp - residual - hh,
        hh_impulse=hh,
        residual_interface_impulse=residual,
        first_hit=first_hit,
    )
    assert outcome["classification"] == FULL_NATURAL_SERVICE_WITNESS
    return outcome, A


@pytest.mark.parametrize("parabolic_scale", [1.0, 1.0e4, 1.0e8])
def test_corridor_identity_rejects_the_same_relative_error_at_every_parabolic_scale(parabolic_scale: float):
    """An exact corridor guard cannot acquire an absolute one-second tolerance."""
    A = parabolic_scale
    expected = 1.0 / A**2
    with pytest.raises(ValueError, match="completed natural corridor"):
        FullNaturalServiceCorridor(
            terminal_time=2.0 * expected,
            endpoint_time=1.5 * expected,
            renewal_frequency=A,
            scaled_lifetime=1.0,
            uniform_service_lower=1.0,
            integrated_service_lower=1.0,
            endpoint_carrier_critical_mass_lower=1.0,
        )


def test_endpoint_cover_cannot_rebind_a_certified_carrier_to_a_foreign_shell_scale():
    outcome, A = _actual_full_natural_outcome()
    true_parent = A / 0.75
    endpoint_hard_shell_cover_from_full_natural_outcome(
        outcome,
        parent_shell_frequency=true_parent,
    )
    with pytest.raises(ValueError, match="frequency|scale|provenance"):
        endpoint_hard_shell_cover_from_full_natural_outcome(
            outcome,
            parent_shell_frequency=100.0 * true_parent,
        )


def test_positive_corridor_service_cannot_be_relabelled_as_a_zero_edge_measure():
    corridor = FullNaturalServiceCorridor(
        terminal_time=2.0,
        endpoint_time=1.0,
        renewal_frequency=1.0,
        scaled_lifetime=1.0,
        uniform_service_lower=0.5,
        integrated_service_lower=0.5,
        endpoint_carrier_critical_mass_lower=0.5,
    )
    with pytest.raises(ValueError, match="positive service|service lower|same measure"):
        material_partition_is_same_corridor_measure(
            corridor,
            edge_weights=[0.0, 0.0],
            old_here=[True, False],
            old_neighbor=[True, False],
        )


def test_positive_tiny_carrier_cannot_be_realized_by_two_zero_hard_shells():
    cover = endpoint_comparable_hard_shell_cover(
        parent_shell_frequency=1.0,
        endpoint_carrier_critical_mass=1.0e-18,
    )
    with pytest.raises(ValueError, match="do not realize"):
        realized_endpoint_hard_shell_witnesses(cover, (0.0, 0.0))


def test_generic_uv_corridor_cannot_be_certified_from_half_of_its_monitor_history():
    A = 1.0e8
    c = 1.0
    T = c / A**2
    mu0 = 2.0
    amp = math.sqrt(1.2 * critical_shell_terminal_mass_lower(mu0) / A)
    elapsed = np.linspace(0.0, 0.5 * T, 4)
    first_hit = critical_shell_backward_first_hit(
        elapsed,
        terminal_amplitude=amp,
        strain_action=np.zeros(4),
        residual_impulse_abs=np.zeros(4),
        hh_impulse_abs=np.zeros(4),
    )
    with pytest.raises(ValueError, match="do not cover"):
        critical_shell_natural_outcome(
            event_time=2.0 * T,
            renewal_frequency=A,
            shell_critical_mass_lower=mu0,
            scaled_lifetime=c,
            viscosity=1.0,
            terminal_coefficient=amp,
            endpoint_coefficient=amp,
            hh_impulse=0j,
            residual_interface_impulse=0j,
            first_hit=first_hit,
        )


def test_high_strain_uv_corridor_requires_a_full_monitor_horizon():
    A = 1.0e8
    c = 1.0
    T = c / A**2
    amp = math.sqrt(1.1 * renewal_carrier_critical_mass_lower(c) / A)
    elapsed = np.linspace(0.0, 0.5 * T, 4)
    first_hit = critical_seed_backward_first_hit(
        elapsed,
        terminal_amplitude=amp,
        strain_action=np.zeros(4),
        residual_impulse_abs=np.zeros(4),
        hh_impulse_abs=np.zeros(4),
    )
    with pytest.raises(ValueError, match="do not cover"):
        critical_seed_natural_outcome(
            event_time=2.0 * T,
            renewal_frequency=A,
            scaled_lifetime=c,
            viscosity=1.0,
            terminal_coefficient=amp,
            endpoint_coefficient=amp,
            hh_impulse=0j,
            residual_interface_impulse=0j,
            first_hit=first_hit,
        )


def test_high_strain_monitor_thresholds_are_bound_to_the_actual_terminal_amplitude():
    A = 3.0
    c = 1.0
    T = c / A**2
    amp = math.sqrt(1.1 * renewal_carrier_critical_mass_lower(c) / A)
    elapsed = np.linspace(0.0, T, 4)
    first_hit = critical_seed_backward_first_hit(
        elapsed,
        terminal_amplitude=10.0 * amp,
        strain_action=np.zeros(4),
        residual_impulse_abs=np.zeros(4),
        hh_impulse_abs=np.zeros(4),
    )
    with pytest.raises(ValueError, match="thresholds do not match"):
        critical_seed_natural_outcome(
            event_time=2.0 * T,
            renewal_frequency=A,
            scaled_lifetime=c,
            viscosity=1.0,
            terminal_coefficient=amp,
            endpoint_coefficient=amp,
            hh_impulse=0j,
            residual_interface_impulse=0j,
            first_hit=first_hit,
        )
