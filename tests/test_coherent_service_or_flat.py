import math

from src.coherent_affine_projection import coherent_deformation_to_dissipation_constant
from src.coherent_service_or_flat import (
    coherent_service_or_flat_gate,
    coherent_deformation_dissipation_threshold,
    coherent_flat_thresholds,
    coherent_nonaffine_coefficient,
)


def _base():
    tau = 0.01
    th = coherent_flat_thresholds(tau)
    return tau, th


def test_default_coherent_flat_branch_has_three_tau_over_three_pieces():
    tau, th = _base()
    out = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=0.9 * th["block_transfer_deficit"],
        objective_variation_action=0.9 * th["objective_strain_variation_action"],
        total_strain_action=0.9 * th["low_strain_action"],
        coherent_deformation_action=0.9 * th["coherent_deformation_action"],
        aspect=1.05,
        scale_radius=3.0,
        has_predecessor=True,
        scaled_lifetime=1.0,
        phase_holonomy=0.1,
    )
    assert out["status"] == "coherent_kelvin_extremal_flat"
    assert out["triggered_causes"] == ()
    assert out["profile_persistence_required"] is False
    assert out["hodge_rms"] < tau / 3
    assert out["nonconformal_strain_number"] < tau / 3
    assert out["coherent_nonaffine_connection_action"] < tau / 3


def test_coherent_deformation_threshold_routes_to_critical_dissipation():
    tau, th = _base()
    out = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=0.0,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=1.01 * th["coherent_deformation_action"],
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=True,
        scaled_lifetime=0.5,
    )
    assert out["status"] == "named_physical_causes"
    root = next(r for r in out["triggered_causes"] if r["cause"] == "coherent_deformation_critical_dissipation")
    assert root["normalized_dissipation_lower"] > 0
    assert math.isclose(
        coherent_deformation_dissipation_threshold(tau, 0.5),
        th["coherent_deformation_action"] ** 2 /
        (0.5 * coherent_deformation_to_dissipation_constant()),
    )


def test_large_radius_is_sticky_ancestry_not_uniform_reset():
    tau, th = _base()
    out = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=0.0,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.01 * th["radius_cap"],
        has_predecessor=True,
        scaled_lifetime=1.0,
    )
    root = next(r for r in out["triggered_causes"] if r["cause"] == "large_affine_radius_ancestry")
    assert math.isclose(root["critical_mass_lower"], 0.3 * th["radius_cap"])


def test_nonaffine_coefficient_is_exact_gaussian_core_bound():
    q = 2.0
    assert math.isclose(coherent_nonaffine_coefficient(q), 1 + q / math.sqrt(2) + math.sqrt(7) / 2)


def test_averaged_high_strain_inherits_existing_critical_dissipation_collision():
    tau, th = _base()
    out = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=0.0,
        objective_variation_action=0.0,
        total_strain_action=1.01 * th["low_strain_action"],
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=True,
        scaled_lifetime=1.0,
    )
    root = next(r for r in out["triggered_causes"] if r["cause"] == "high_strain_critical_dissipation")
    assert root["normalized_dissipation_lower"] > 0


def test_gate_reports_simultaneous_physical_causes_without_lexicographic_primary():
    tau, th = _base()
    out = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=1.1 * th["block_transfer_deficit"],
        objective_variation_action=1.1 * th["objective_strain_variation_action"],
        total_strain_action=1.1 * th["low_strain_action"],
        coherent_deformation_action=1.1 * th["coherent_deformation_action"],
        aspect=1.1 * th["aspect_threshold"],
        scale_radius=1.1 * th["radius_cap"],
        has_predecessor=False,
        scaled_lifetime=1.0,
        phase_holonomy=1.1 * th["phase_holonomy"],
    )
    assert out["status"] == "named_physical_causes"
    assert out["primary_selected"] is False
    assert len(out["triggered_causes"]) == 7
