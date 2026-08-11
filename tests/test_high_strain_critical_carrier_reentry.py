import math

import numpy as np
import pytest

from src.common_slice_coefficient_registration import (
    HH_COEFFICIENT_OBSTRUCTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
)
from src.high_strain_critical_carrier_reentry import (
    CriticalDissipationAtom,
    critical_seed_backward_first_hit,
    critical_seed_natural_outcome,
    pushforward_critical_dissipation_law,
    theorem_certificate,
)
from src.full_natural_service_corridor_quotient import endpoint_hard_shell_cover_from_full_natural_outcome
from src.high_strain_resolved_ancestor import high_strain_ancestor_mass_threshold
from src.nn_critical_heat_carrier_seed import renewal_carrier_critical_mass_lower


def _source_seed(*, scaled_lifetime: float, renewal_frequency: float, event_time: float, terminal_mass: float):
    M = renewal_frequency / 0.75
    return pushforward_critical_dissipation_law(
        (
            CriticalDissipationAtom(
                mass=1.0,
                child_frequency=4.0 * M,
                shell_upper_frequency=M,
                shell_energy_u=terminal_mass / renewal_frequency,
                time=event_time,
            ),
        ),
        scaled_lifetime=scaled_lifetime,
    )[0]


def test_actual_critical_dissipation_law_pushes_to_carrier_seeds_without_nn():
    c = 1.0
    N = 64.0
    mu = high_strain_ancestor_mass_threshold(c)
    atoms = (
        CriticalDissipationAtom(2.0, N, N / 4, 1.2 * mu / (N / 4), 0.1),
        CriticalDissipationAtom(3.0, N, N / 8, 2.0 * mu / (N / 8), 0.2),
    )
    seeds = pushforward_critical_dissipation_law(atoms, scaled_lifetime=c)
    assert math.isclose(sum(x.normalized_dissipation_weight for x in seeds), 1.0)
    assert all(x.renewal_critical_mass >= renewal_carrier_critical_mass_lower(c) for x in seeds)
    assert all(x.natural_lifetime_ratio >= 256 / 9 for x in seeds)


def test_generic_corridor_has_three_native_monitors_and_no_material_monitor():
    amp = 2.0
    ell = np.linspace(0.0, 1.0, 5)
    hit = critical_seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1 / 60, 5),
        residual_impulse_abs=np.linspace(0.0, 0.1 * amp, 5),
        hh_impulse_abs=np.linspace(0.0, 0.2 * amp, 5),
    )
    assert hit["first_elapsed"] is None
    assert set(hit["individual_debuts"]) == {
        "high_strain_critical_dissipation",
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
        HH_COEFFICIENT_OBSTRUCTION,
    }
    assert hit["requires_physical_energy_reentry"] is False
    assert hit["coefficient_impulses_used_as_work"] is False


def test_high_strain_hh_coefficient_hit_requests_energy_reentry_without_work_assignment():
    amp = 2.0
    hit = critical_seed_backward_first_hit(
        np.linspace(0.0, 1.0, 5),
        terminal_amplitude=amp,
        strain_action=np.zeros(5),
        residual_impulse_abs=np.zeros(5),
        hh_impulse_abs=np.linspace(0.0, 1.2, 5),
    )
    assert hit["joint_first_stops"] == (HH_COEFFICIENT_OBSTRUCTION,)
    assert hit["requires_physical_energy_reentry"] is True
    assert hit["coefficient_impulses_used_as_work"] is False


def test_full_generic_critical_corridor_creates_own_scale_service_without_nn():
    c = 1.0
    A = 3.0
    T = c / A**2
    terminal_mass = renewal_carrier_critical_mass_lower(c)
    amp = math.sqrt(terminal_mass / A)
    seed = _source_seed(
        scaled_lifetime=c,
        renewal_frequency=A,
        event_time=2.0 * T,
        terminal_mass=terminal_mass,
    )
    ell = np.linspace(0.0, T, 5)
    hit = critical_seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1 / 60, 5),
        residual_impulse_abs=np.linspace(0.0, 0.1 * amp, 5),
        hh_impulse_abs=np.linspace(0.0, 0.2 * amp, 5),
    )
    ir = 0.1 * amp
    ih = 0.2j * amp
    out = critical_seed_natural_outcome(
        source_seed=seed,
        event_time=2.0 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
        viscosity=1.0,
        terminal_coefficient=amp,
        endpoint_coefficient=amp - ir - ih,
        hh_impulse=ih,
        residual_interface_impulse=ir,
        first_hit=hit,
    )
    assert out["classification"] == "full_natural_own_scale_service"
    assert out["uniform_square_service_lower"] > 0
    assert out["integrated_bounded_heat_service_lower"] > 0
    assert out["nn_seed_required"] is False
    assert out["materiality_assigned"] == "only_after_service_via_exact_Moyal_OO_ON_NN"
    assert out["service_same_corridor_witness"] is True
    assert out["service_adds_recursion_depth"] is False
    assert out["physical_time_drop"] == pytest.approx(T)
    assert out["corridor_endpoint_time"] == pytest.approx(T)
    cover = endpoint_hard_shell_cover_from_full_natural_outcome(out, parent_shell_frequency=A / 0.75)
    assert cover["candidate_ratios_to_parent"] == pytest.approx((0.75, 1.5))
    assert cover["guaranteed_max_hard_shell_critical_mass_lower"] == pytest.approx((2.0/3.0) * out["endpoint_carrier_critical_mass_lower"])


def test_high_strain_or_hh_hit_stays_named_recursive_not_service():
    c = 1.0
    A = 2.0
    T = c / A**2
    terminal_mass = 1.1 * renewal_carrier_critical_mass_lower(c)
    amp = math.sqrt(terminal_mass / A)
    seed = _source_seed(
        scaled_lifetime=c,
        renewal_frequency=A,
        event_time=2.0 * T,
        terminal_mass=terminal_mass,
    )
    ell = np.linspace(0.0, T, 5)
    hit = critical_seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1 / 20, 5),
        residual_impulse_abs=np.zeros(5),
        hh_impulse_abs=np.zeros(5),
    )
    out = critical_seed_natural_outcome(
        source_seed=seed,
        event_time=2 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
        viscosity=1.0,
        terminal_coefficient=amp,
        endpoint_coefficient=amp,
        hh_impulse=0j,
        residual_interface_impulse=0j,
        first_hit=hit,
    )
    assert out["classification"] == "named_first_stop"
    assert "high_strain_critical_dissipation" in out["joint_first_stops"]


def test_certificate_explicitly_bypasses_nn_only_for_renewal_entrance():
    cert = theorem_certificate()
    assert "NN_NOT_REQUIRED_FOR_RENEWAL_ENTRANCE" in cert["status"]
    assert "no longer requires" in cert["architectural_shortcut"]
    assert "remain valid refinements" in cert["architectural_shortcut"]
    assert "source/SGS" in cert["scope"]
