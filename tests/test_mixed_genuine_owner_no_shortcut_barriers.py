import math

from src.mixed_genuine_owner_no_shortcut_barriers import (
    certified_no_shortcut_barriers,
    critical_energy_geometric_chain,
    theorem_certificate,
)


def test_geometric_critical_shell_allocation_has_finite_energy_budget():
    out = critical_energy_geometric_chain(8.0, 3.0, 40, 2.0)
    assert out.shell_frequencies[1] == 16.0
    assert math.isclose(out.shell_energies[0], 3.0 / 8.0)
    assert out.total_shell_energy < out.infinite_chain_energy_upper
    assert math.isclose(out.infinite_chain_energy_upper, 0.75)


def test_existing_certified_theorems_forbid_scalar_shortcuts():
    out = certified_no_shortcut_barriers()
    assert out["critical_reset_forbidden"]
    assert out["gross_radial_budget_forbidden"]
    assert out["generic_shell_progress_forbidden"]
    assert out["fresh_scale_direction_forbidden"]
    assert not out["claims_global_regularity"]


def test_barrier_certificate_states_budget_model_is_not_ns_existence_claim():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "not an NS trajectory-existence statement" in cert["critical_energy_scaling_counterexample"]
    assert "do not prove or disprove" in cert["scope"]
