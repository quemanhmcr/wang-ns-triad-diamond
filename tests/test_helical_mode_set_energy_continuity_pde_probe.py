import pytest

from src.helical_mode_set_energy_continuity_pde_probe import STATUS, run_probe
from src.mixed_fate_reserved_young_handoff_pde_probe import HELICITIES


def test_actual_ns_helical_mode_uses_gross_canonical_edge_hahn_work_in_interval_balance():
    out = run_probe(resolutions=(24,), cutoff=7, steps=16, duration=0.0005, helicity=HELICITIES[2])
    assert out.status == STATUS
    run = out.runs[0]
    assert run.initial_energy > 0.0
    assert run.positive_work_steps > 0
    assert run.negative_work_steps > 0
    assert run.interval_continuity_native_residual < 5e-5
    assert run.worst_instantaneous_signed_reconstruction_native_residual < 5e-9


def test_actual_ns_phase_reversal_still_obeys_same_stock_work_viscosity_continuity():
    out = run_probe(
        resolutions=(24,), cutoff=7, steps=16, duration=0.0005,
        helicity=HELICITIES[2], phase_sign=-1,
    )
    run = out.runs[0]
    assert run.initial_energy > 0.0
    assert run.interval_continuity_native_residual < 5e-5
    assert run.worst_instantaneous_signed_reconstruction_native_residual < 5e-9


def test_opposite_child_helicity_is_a_real_mode_even_when_initial_stock_is_zero_or_tiny():
    out = run_probe(
        resolutions=(24,), cutoff=7, steps=16, duration=0.0005,
        helicity=-HELICITIES[2],
    )
    run = out.runs[0]
    assert run.final_energy >= 0.0
    assert run.interval_continuity_native_residual < 5e-5


def test_same_cutoff_helical_mode_continuity_is_fft_representation_invariant():
    out = run_probe(resolutions=(24, 28), cutoff=7, steps=16, duration=0.0005, helicity=HELICITIES[2])
    assert out.maximum_initial_energy_representation_relative_residual < 5e-7
    assert out.maximum_final_energy_representation_relative_residual < 5e-7
    assert out.maximum_integrated_positive_work_representation_relative_residual < 5e-7
    assert out.maximum_integrated_negative_work_representation_relative_residual < 5e-7
    assert out.maximum_viscous_dissipation_representation_relative_residual < 5e-7


def test_invalid_phase_sign_is_rejected():
    with pytest.raises(ValueError, match="phase_sign"):
        run_probe(resolutions=(24,), steps=16, phase_sign=0)
