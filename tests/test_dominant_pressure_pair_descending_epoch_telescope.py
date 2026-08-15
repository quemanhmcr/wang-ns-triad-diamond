import math
import pytest

from src.dominant_pressure_pair_descending_epoch_telescope import (
    DominantPressurePairRenewalStep,
    dominant_pressure_pair_child_mass_floor,
    dominant_pressure_pair_epoch_telescope,
    dominant_pressure_pair_physical_frequency_floor,
    dominant_pressure_pair_source_floor_from_objective_stop,
    theorem_certificate,
)


def _step(parent: float, child: float, sigma: float = 1.0, c: float = 1.0):
    mu = dominant_pressure_pair_child_mass_floor(sigma, c)
    return DominantPressurePairRenewalStep(
        parent_frequency=parent,
        child_frequency=child,
        pressure_source_weight=sigma,
        normalized_pair_mass=0.25,
        child_critical_mass=mu,
        scaled_lifetime=c,
    )


def test_objective_stop_floor_gives_twenty_A_over_c_squared_child_mass():
    A = 3.0
    c = 2.0
    sigma = dominant_pressure_pair_source_floor_from_objective_stop(A, c)
    mu = dominant_pressure_pair_child_mass_floor(sigma, c)
    assert sigma == pytest.approx(A / (4.0 * c))
    assert mu == pytest.approx(20.0 * A / (c * c))


def test_physical_frequency_floor_comes_from_actual_shell_mass_and_global_energy():
    sigma = 2.0
    c = 1.0
    E = 40.0
    mu = dominant_pressure_pair_child_mass_floor(sigma, c)
    out = dominant_pressure_pair_physical_frequency_floor(E, sigma, c)
    assert mu == pytest.approx(160.0)
    assert out == pytest.approx(4.0)


def test_four_quarter_descents_close_at_global_energy_frequency_floor():
    # sigma=1,c=1 gives mu*=80.  E=10 gives N_min=8.
    # N0=2048 permits at most log_4(2048/8)=4 transitions.
    rows = (
        _step(2048.0, 512.0),
        _step(512.0, 128.0),
        _step(128.0, 32.0),
        _step(32.0, 8.0),
    )
    out = dominant_pressure_pair_epoch_telescope(
        rows,
        global_energy_upper=10.0,
        pressure_source_weight_floor=1.0,
        scaled_lifetime=1.0,
    )
    assert out.transition_count == 4
    assert out.maximum_transition_count == 4
    assert out.physical_frequency_floor == pytest.approx(8.0)
    assert out.final_child_frequency == pytest.approx(8.0)
    assert out.maximum_observed_scale_ratio == pytest.approx(0.25)
    assert not out.pressure_entropy_used_as_cost
    assert not out.critical_mass_used_as_additive_reset


def test_fifth_quarter_descent_is_physically_incompatible_with_same_mass_floor():
    rows = (
        _step(2048.0, 512.0),
        _step(512.0, 128.0),
        _step(128.0, 32.0),
        _step(32.0, 8.0),
        _step(8.0, 2.0),
    )
    with pytest.raises(ValueError):
        dominant_pressure_pair_epoch_telescope(
            rows,
            global_energy_upper=10.0,
            pressure_source_weight_floor=1.0,
            scaled_lifetime=1.0,
        )


def test_nonconsecutive_child_registration_is_rejected():
    rows = (_step(512.0, 128.0), _step(64.0, 16.0))
    with pytest.raises(ValueError, match="selected child shell"):
        dominant_pressure_pair_epoch_telescope(
            rows,
            global_energy_upper=100.0,
            pressure_source_weight_floor=1.0,
            scaled_lifetime=1.0,
        )


def test_diffuse_pair_is_not_admitted_to_dominant_epoch():
    with pytest.raises(ValueError, match="quarter-dominant"):
        DominantPressurePairRenewalStep(
            parent_frequency=100.0,
            child_frequency=25.0,
            pressure_source_weight=1.0,
            normalized_pair_mass=0.2,
            child_critical_mass=80.0,
            scaled_lifetime=1.0,
        )


def test_generic_shell_three_quarters_registration_cannot_enter_pressure_epoch():
    with pytest.raises(ValueError, match="N_next<=N/4"):
        _step(100.0, 75.0)


def test_exact_power_of_four_count_has_no_off_by_one():
    sigma = 1.0
    c = 1.0
    E = 1.0
    # mu*=80, Nmin=80, root/Nmin=4^3.
    rows = (
        _step(5120.0, 1280.0, sigma, c),
        _step(1280.0, 320.0, sigma, c),
        _step(320.0, 80.0, sigma, c),
    )
    out = dominant_pressure_pair_epoch_telescope(
        rows,
        global_energy_upper=E,
        pressure_source_weight_floor=sigma,
        scaled_lifetime=c,
    )
    assert out.maximum_transition_count == 3
    assert math.isclose(out.final_child_frequency, 80.0)


def test_certificate_keeps_diffuse_and_mixed_recurrence_open():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "diffuse pressure H2 is not used as a cost" in cert["pressure_entropy"]
    assert "cross-family source/strain/HH recurrence" in cert["scope"]
