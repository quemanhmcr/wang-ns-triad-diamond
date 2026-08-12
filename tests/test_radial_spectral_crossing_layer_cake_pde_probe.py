from __future__ import annotations

import pytest

from src.radial_spectral_crossing_layer_cake_pde_probe import STATUS, run_probe


def test_actual_ns_radial_tail_and_selected_upward_crossing():
    out = run_probe(resolutions=(24,), cutoff=7, steps=16, duration=0.0005, radial_boundary=8.0, phase_sign=1)
    run = out.runs[0]
    assert out.status == STATUS
    assert run.tail_interval_balance_native_residual < 5e-5
    assert run.initial_selected_upward_crossing > 0.0
    assert run.initial_selected_downward_crossing == pytest.approx(0.0, abs=1e-12)
    assert run.worst_selected_radial_divergence_native_residual < 1e-9
    assert run.worst_selected_truncated_layer_cake_native_residual < 1e-9


def test_sign_reversed_actual_ns_has_initial_downward_crossing():
    out = run_probe(resolutions=(24,), cutoff=7, steps=16, duration=0.0005, radial_boundary=8.0, phase_sign=-1)
    run = out.runs[0]
    assert run.initial_selected_downward_crossing > 0.0
    assert run.initial_selected_upward_crossing == pytest.approx(0.0, abs=1e-12)
    assert run.tail_interval_balance_native_residual < 5e-5


def test_same_cutoff_radial_ns_is_fft_representation_invariant_on_native_envelopes():
    out = run_probe(resolutions=(24, 28), cutoff=7, steps=16, duration=0.0005, radial_boundary=8.0, phase_sign=1)
    assert out.tail_representation_native_scale > 0.0
    assert out.selected_crossing_representation_native_scale > 0.0
    assert out.maximum_initial_tail_energy_representation_native_residual < 5e-7
    assert out.maximum_final_tail_energy_representation_native_residual < 5e-7
    assert out.maximum_integrated_tail_work_representation_native_residual < 5e-7
    assert out.maximum_tail_viscosity_representation_native_residual < 5e-7
    assert out.maximum_integrated_selected_upward_representation_native_residual < 5e-7
    assert out.maximum_integrated_selected_downward_representation_native_residual < 5e-7
