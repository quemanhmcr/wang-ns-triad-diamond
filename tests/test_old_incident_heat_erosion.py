import math

import pytest

from src.old_incident_heat_erosion import (
    CLEAN_HEAT_POOL_RATIO,
    band_addressed_material_partition,
    canonical_nn_critical_thresholds,
    first_nn_critical_generation,
    forced_nn_service_lower,
    old_incident_heat_capacity_upper,
    positive_measure_intersection_lower,
    theorem_certificate,
)
from src.high_strain_heat_increment_service import high_strain_heat_service_lower


def test_old_incident_is_exactly_oo_plus_on_and_below_old_shell_law():
    out = band_addressed_material_partition(
        [1.0, 2.0, 3.0, 4.0],
        [True, True, True, False],
        [True, True, False, False],
        [True, False, False, False],
    )
    assert math.isclose(out["old_old"], 1.0)
    assert math.isclose(out["old_new_interface"], 2.0)
    assert math.isclose(out["old_incident"], 3.0)
    assert math.isclose(out["old_shell_service"], 6.0)
    assert math.isclose(out["new_new"], 7.0)
    assert math.isclose(out["ownership_partition_residual"], 0.0)
    assert math.isclose(out["incident_identity_residual"], 0.0)
    assert out["old_shell_capacity_margin"] >= 0.0


def test_old_endpoint_outside_old_shell_is_rejected():
    with pytest.raises(ValueError):
        band_addressed_material_partition([1.0], [False], [True], [False])


def test_old_incident_capacity_inherits_heat_specific_geometric_ratio():
    c0 = old_incident_heat_capacity_upper(
        generation=0,
        initial_low_cut_ratio=0.25,
        initial_block_frequency=64.0,
        frame_energy_bound=1.0,
        global_energy=3.0,
    )
    c1 = old_incident_heat_capacity_upper(
        generation=1,
        initial_low_cut_ratio=0.25,
        initial_block_frequency=64.0,
        frame_energy_bound=1.0,
        global_energy=3.0,
    )
    assert math.isclose(c1 / c0, float(CLEAN_HEAT_POOL_RATIO))


def test_canonical_nn_critical_overlap_is_half_of_existing_good_fraction():
    th = canonical_nn_critical_thresholds()
    g = th["critical_heat_fraction"]
    eps = th["old_incident_fraction_target"]
    assert math.isclose(eps, g / 2.0)
    assert math.isclose(th["nn_critical_intersection_fraction_lower"], g / 2.0)
    assert th["new_new_fraction_lower"] > 0.75


def test_positive_measure_intersection_is_sharp_inclusion_exclusion():
    assert math.isclose(positive_measure_intersection_lower(10.0, 8.0, 4.0), 2.0)
    assert math.isclose(positive_measure_intersection_lower(10.0, 5.0, 4.0), 0.0)


def test_finite_age_forces_nn_critical_seed():
    c = 1.0
    Sstar = high_strain_heat_service_lower(c)
    eps = canonical_nn_critical_thresholds()["old_incident_fraction_target"]
    C0 = 100.0 * Sstar
    q = first_nn_critical_generation(scaled_lifetime=c, initial_old_capacity=C0)
    r = float(CLEAN_HEAT_POOL_RATIO)
    assert C0 * r**q <= eps * Sstar
    assert q == 0 or C0 * r ** (q - 1) > eps * Sstar


def test_nn_remainder_uses_total_positive_heat_law():
    assert math.isclose(forced_nn_service_lower(5.0, 1.25), 3.75)


def test_certificate_does_not_claim_cell_mass_or_universal_renewal():
    cert = theorem_certificate()
    assert "OO+ON" in cert["old_incident"]
    assert "does not give a per-cell mass floor" in cert["scope"]
    assert "universal epoch/slab renewal" in cert["scope"]
