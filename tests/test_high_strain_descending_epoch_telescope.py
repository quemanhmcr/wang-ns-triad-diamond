import math

import pytest

from src.high_strain_descending_epoch_telescope import (
    HIGH_STRAIN_RENEWAL_RATIO_UPPER,
    STATUS,
    HighStrainRenewalStep,
    high_strain_epoch_telescope,
    kinetic_energy_gradient_dissipation_upper,
    physical_high_strain_frequency_floor,
    theorem_certificate,
)
from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower
from src.nn_critical_heat_carrier_seed import renewal_scale


def _step(N: float, D: float, ratio: float = 3.0 / 16.0) -> HighStrainRenewalStep:
    M = (ratio / (3.0 / 4.0)) * N
    return HighStrainRenewalStep(N, M, renewal_scale(M), D)


def test_exact_high_strain_ancestor_renewal_descends_by_at_most_three_sixteenths():
    N = 64.0
    row = _step(N, 2.0)
    assert row.ancestor_shell_frequency == pytest.approx(N / 4.0)
    assert row.renewal_frequency == pytest.approx((3.0 / 16.0) * N)
    assert row.renewal_frequency / row.child_frequency == pytest.approx(HIGH_STRAIN_RENEWAL_RATIO_UPPER)


def test_complete_time_overlap_is_allowed_because_scale_weights_sum_geometrically():
    c = 1.0
    Dstar = clean_high_strain_dissipation_lower(c)
    G = 12.0
    N0 = 100.0 * Dstar / G
    rows = []
    N = N0
    for _ in range(3):
        rows.append(_step(N, 1.05 * Dstar))
        N = rows[-1].renewal_frequency

    out = high_strain_epoch_telescope(rows, total_gradient_dissipation=G, scaled_lifetime=c)
    assert out.arbitrary_time_overlap_allowed is True
    assert out.normalized_dissipation_used_as_global_reset is False
    assert out.step_count == 3
    assert out.normalized_dissipation_sum <= out.normalized_dissipation_capacity_upper
    assert sum(x.child_frequency for x in rows) <= out.geometric_frequency_sum_upper


def test_high_strain_event_has_physical_frequency_floor_from_global_gradient_reservoir():
    c = 0.8
    G = 7.0
    floor = physical_high_strain_frequency_floor(G, c)
    Dstar = clean_high_strain_dissipation_lower(c)
    assert floor == pytest.approx(Dstar / G)

    row = _step(1.1 * floor, Dstar)
    out = high_strain_epoch_telescope((row,), total_gradient_dissipation=G, scaled_lifetime=c)
    assert out.last_child_frequency >= out.physical_frequency_floor
    assert out.physical_frequency_floor == pytest.approx(floor)


def test_unforced_ns_energy_inequality_converts_initial_energy_to_gradient_reservoir_upper():
    E0 = 9.0
    nu = 0.75
    Gupper = kinetic_energy_gradient_dissipation_upper(E0, nu)
    assert Gupper == pytest.approx(E0 / (2.0 * nu))
    with pytest.raises(ValueError, match="positive viscosity"):
        kinetic_energy_gradient_dissipation_upper(E0, 0.0)


def test_nonconsecutive_restart_cannot_be_hidden_inside_one_high_strain_epoch():
    c = 1.0
    Dstar = clean_high_strain_dissipation_lower(c)
    G = 20.0
    first = _step(20.0, Dstar)
    # This is individually a valid high-strain step, but it is not the carrier
    # produced by the previous physical ancestor/renewal law.
    second = _step(0.9 * first.child_frequency, Dstar)
    with pytest.raises(ValueError, match="actual renewed carrier scale consecutively"):
        high_strain_epoch_telescope((first, second), total_gradient_dissipation=G, scaled_lifetime=c)


def test_normalized_dissipation_above_global_reservoir_capacity_is_rejected():
    c = 1.0
    G = 2.0
    N = 5.0
    row = _step(N, N * G + 1.0)
    with pytest.raises(ValueError, match="exceeds the supplied global gradient reservoir"):
        high_strain_epoch_telescope((row,), total_gradient_dissipation=G, scaled_lifetime=c)


def test_certificate_closes_only_eventually_pure_high_strain_recurrence():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "N_next/N<=3/16" in cert["physical_renewal_scale"]
    assert "overlap arbitrarily" in cert["global_reservoir"]
    assert "D_j>=D_*" in cert["frequency_floor"]
    assert "not promoted" in cert["anti_reset"]
    assert "non-high-strain epoch breakers" in cert["master_consequence"]
    assert "mixed recurrence" in cert["scope"]
    assert "no Navier-Stokes global-regularity claim" in cert["scope"]


def test_frequency_floor_count_is_finite_on_exact_geometric_descent():
    c = 1.0
    G = 10.0
    Dstar = clean_high_strain_dissipation_lower(c)
    floor = Dstar / G
    N = 100.0 * floor
    rows = []
    while N >= floor and len(rows) < 20:
        rows.append(_step(N, Dstar))
        N = rows[-1].renewal_frequency
    out = high_strain_epoch_telescope(rows, total_gradient_dissipation=G, scaled_lifetime=c)
    assert math.isfinite(out.certified_count_upper)
    assert out.step_count <= out.certified_count_upper
