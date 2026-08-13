from src.same_carrier_inherited_energy_relay_pde_probe import run_probe


def test_actual_ns_endpoint_stock_realizes_inheritance_face_without_new_temporal_matching():
    out = run_probe(resolutions=(24,), steps=16, duration=0.00025)
    assert out.observations
    assert out.minimum_inheritance_gate_margin > 0.0
    assert out.maximum_mode_continuity_native_residual <= 5e-5
    assert out.maximum_global_energy_balance_relative_residual <= 5e-5
    o = out.observations[0]
    assert o.inherited_fraction >= 0.2
    assert o.initial_energy > 0.0 and o.final_energy > 0.0
