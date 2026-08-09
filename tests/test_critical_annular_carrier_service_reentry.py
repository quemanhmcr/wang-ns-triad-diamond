import math

import numpy as np
import pytest

from src.critical_annular_carrier_service_reentry import (
    BOUNDED_HEAT_RADIUS,
    arb_bounded_heat_probe_certificate,
    bounded_heat_defect_fraction_lower,
    gaussian_3d_tail_probability,
    heat_defect_fraction_lower,
    integrated_bounded_heat_service_lower,
    material_service_partition,
    no_hit_prefix_amplitude_lower,
    persistent_carrier_critical_mass_lower,
    renewed_analysis_probe_growth_upper,
    service_epoch_reentry_certificate,
    theorem_certificate,
    transported_annular_support_ratios,
    uniform_bounded_square_service_lower,
)
from src.nn_seed_temporal_first_stop import inherited_seed_critical_mass_lower


def test_every_strict_no_hit_prefix_keeps_more_than_quarter_terminal_amplitude():
    amp = 7.0
    for ir_frac, hh_frac in ((0.0, 0.0), (0.24, 0.49), (0.10, 0.20)):
        low = no_hit_prefix_amplitude_lower(amp, ir_frac * amp, hh_frac * amp)
        assert low > amp / 4


def test_viscous_adjoint_cost_is_scale_independent_and_kept_in_carrier_mass():
    c = 1.0
    nu = 1.0
    J = renewed_analysis_probe_growth_upper(c, nu)
    assert J > 1.0
    expected = inherited_seed_critical_mass_lower(c) / J**2
    assert math.isclose(persistent_carrier_critical_mass_lower(c, nu), expected, rel_tol=1e-14)


def test_transport_keeps_annular_lower_edge_above_low_low_radius():
    lo, hi = transported_annular_support_ratios()
    assert lo > 0.5
    assert hi > 1.5


def test_radius_three_retains_clean_half_of_annular_heat_lower():
    pytest.importorskip("flint")
    cert = arb_bounded_heat_probe_certificate()
    assert cert["status"].startswith("ARB_CERTIFIED_RADIUS_3")
    full = heat_defect_fraction_lower()
    bounded = bounded_heat_defect_fraction_lower()
    assert math.isclose(bounded, 0.5 * full, rel_tol=1e-14)
    assert full - 4.0 * gaussian_3d_tail_probability(BOUNDED_HEAT_RADIUS) > bounded


def test_uniform_bounded_service_is_strictly_positive_and_integrates_with_natural_lifetime():
    for c, nu in ((0.5, 0.0), (1.0, 1.0), (2.0, 1.5)):
        y = uniform_bounded_square_service_lower(c, nu)
        assert y > 0
        assert math.isclose(integrated_bounded_heat_service_lower(c, nu), c * y, rel_tol=1e-14)


def test_materiality_is_reread_from_the_new_positive_service_law():
    w = np.array([1.0, 2.0, 3.0, 4.0])
    old0 = np.array([True, True, False, False])
    old1 = np.array([True, False, True, False])
    out = material_service_partition(w, old0, old1)
    assert out["old_old"] == 1.0
    assert out["old_new_interface"] == 5.0
    assert out["new_new"] == 4.0
    assert out["total"] == 10.0
    assert out["partition_residual"] == 0.0


def test_uniform_service_threshold_enters_existing_finite_old_pool_stopping_epoch():
    y = uniform_bounded_square_service_lower(1.0, 1.0)
    out = service_epoch_reentry_certificate(1.0, 1000.0 * y, 1.0)
    assert out["renewed_service_threshold"] == y
    assert out["bounded_displacement_radius_over_A"] == 3.0
    assert out["first_forced_generation"] >= 1
    assert out["old_capacity_at_forced_generation"] <= out["target_old_capacity"]


def test_certificate_does_not_promote_service_to_hh_efficiency_or_whole_carrier_nn():
    pytest.importorskip("flint")
    cert = theorem_certificate()
    assert "not near-extremal HH transfer efficiency" in cert["no_efficiency_overclaim"]
    assert "not whole-carrier ownership" in cert["material_reentry"]
    assert "universal source/relink" in cert["scope"]
