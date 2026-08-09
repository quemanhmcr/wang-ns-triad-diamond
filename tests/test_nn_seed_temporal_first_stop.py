import math

import numpy as np

from src.nn_critical_heat_carrier_seed import renewal_carrier_critical_mass_lower
from src.nn_seed_temporal_first_stop import (
    backward_natural_endpoint,
    inherited_seed_critical_mass_lower,
    natural_corridor_outcome,
    renewed_natural_duration,
    seed_backward_first_hit,
    theorem_certificate,
)


def test_clean_inherited_seed_mass_is_pi_squared_over_50c2():
    for c in (0.5, 1.0, 2.3):
        assert math.isclose(inherited_seed_critical_mass_lower(c), math.pi**2 / (50 * c**2))


def test_backward_natural_endpoint_truncates_at_initial_surface():
    A = 2.0
    c = 1.0
    T = renewed_natural_duration(A, c)
    out = backward_natural_endpoint(0.4 * T, A, c)
    assert out["hits_initial_boundary"] is True
    assert out["backward_endpoint"] == 0.0
    out2 = backward_natural_endpoint(2.0 * T, A, c)
    assert out2["hits_initial_boundary"] is False
    assert math.isclose(float(out2["elapsed_available"]), T)


def test_no_hit_full_natural_corridor_keeps_quarter_coefficient_and_critical_mass():
    c = 1.0
    A = 3.0
    T = renewed_natural_duration(A, c)
    terminal_mass = renewal_carrier_critical_mass_lower(c)
    amp = math.sqrt(terminal_mass / A)
    ell = np.linspace(0.0, T, 5)
    IR = np.linspace(0.0, 0.1 * amp, 5)
    IH = np.linspace(0.0, 0.2 * amp, 5)
    hit = seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.linspace(0.0, 1 / 60, 5),
        residual_impulse_abs=IR,
        hh_impulse_abs=IH,
        material_boundary_distance=np.linspace(1.0, 0.6, 5),
    )
    ir = 0.1 * amp
    ih = 0.2j * amp
    zs = amp - ir - ih
    out = natural_corridor_outcome(
        event_time=2.0 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
        terminal_coefficient=amp,
        endpoint_coefficient=zs,
        hh_impulse=ih,
        residual_interface_impulse=ir,
        first_hit=hit,
    )
    assert out["classification"] == "full_natural_corridor_survivor"
    assert float(out["endpoint_amplitude"]) >= amp / 4
    assert float(out["endpoint_critical_mass"]) >= inherited_seed_critical_mass_lower(c) - 1e-12
    assert out["nn_endpoint_witness_survives"] is True
    assert out["whole_carrier_declared_nn"] is False


def test_material_boundary_contact_is_named_first_stop():
    amp = 2.0
    ell = np.linspace(0.0, 1.0, 5)
    hit = seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.zeros(5),
        residual_impulse_abs=np.zeros(5),
        hh_impulse_abs=np.zeros(5),
        material_boundary_distance=np.linspace(1.0, 0.0, 5),
    )
    assert "material_boundary_contact" in hit["joint_causes"]


def test_initial_boundary_is_absorbing_even_without_dynamic_hit():
    c = 1.0
    A = 2.0
    T = renewed_natural_duration(A, c)
    amp = math.sqrt(renewal_carrier_critical_mass_lower(c) / A)
    ell = np.linspace(0.0, 0.5 * T, 5)
    hit = seed_backward_first_hit(
        ell,
        terminal_amplitude=amp,
        strain_action=np.zeros(5),
        residual_impulse_abs=np.zeros(5),
        hh_impulse_abs=np.zeros(5),
        material_boundary_distance=np.ones(5),
    )
    out = natural_corridor_outcome(
        event_time=0.5 * T,
        renewal_frequency=A,
        scaled_lifetime=c,
        terminal_coefficient=amp,
        endpoint_coefficient=amp,
        hh_impulse=0j,
        residual_interface_impulse=0j,
        first_hit=hit,
    )
    assert out["classification"] == "initial_boundary_root"


def test_certificate_keeps_interface_and_hh_impulses_out_of_physical_work_law():
    cert = theorem_certificate()
    assert "not promoted to physical work" in cert["interface_scope"]
    assert "physical-energy causal gate" in cert["hh_scope"]
    assert "NN material witness" in cert["scope"]
    assert "efficiency/service renewal" in cert["scope"]
