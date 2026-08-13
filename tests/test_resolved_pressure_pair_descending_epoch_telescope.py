import math
import pytest

from src.resolved_pressure_pair_descending_epoch_telescope import (
    ResolvedPressurePairRenewalStep,
    certify_resolved_pressure_pair_owner,
    resolved_pressure_pair_epoch_telescope,
    resolved_pressure_pair_parent_frequency_floor,
    resolved_pressure_pair_source_floor_from_objective_stop,
    theorem_certificate,
)


def _owner(N: float, *, sigma: float = 1.0, c: float = 1.0, atoms: int = 5):
    # Equal positive pair atoms make qmax=1/atoms.  atoms=5 is diffuse (<1/4),
    # demonstrating that the theorem does not use the quarter cut.
    weights = [sigma / atoms] * atoms
    indices = [(j, j) for j in range(atoms)]
    freqs = [(N / 4.0, N / 4.0)] * atoms
    return certify_resolved_pressure_pair_owner(
        sigma,
        c,
        N,
        sgs_positive_source_weight=0.0,
        pair_positive_weights=weights,
        pair_shell_indices=indices,
        pair_frequencies=freqs,
    )


def _step(N: float, *, sigma: float = 1.0, c: float = 1.0, atoms: int = 5):
    owner = _owner(N, sigma=sigma, c=c, atoms=atoms)
    return ResolvedPressurePairRenewalStep(
        owner=owner,
        child_frequency=N / 4.0,
        child_critical_mass=owner.selected_pair_shell_mass_lower,
    )


def test_parent_floor_from_actual_pair_source_and_absolute_capacity():
    # sigma=1,c=1,E=10 -> N_min=128.
    out = resolved_pressure_pair_parent_frequency_floor(10.0, 1.0, 1.0)
    assert out == pytest.approx(128.0)


def test_objective_stop_floor_gives_320_A_over_c_squared_energy_parent_floor():
    A = 2.0
    c = 4.0
    E = 5.0
    sigma = resolved_pressure_pair_source_floor_from_objective_stop(A, c)
    floor = resolved_pressure_pair_parent_frequency_floor(E, sigma, c)
    assert sigma == pytest.approx(A / (4.0 * c))
    assert floor == pytest.approx(320.0 * A / (c * c * E))


def test_diffuse_pressure_pair_epoch_is_finite_without_entropy_cost():
    # qmax=1/5<1/4 on every row.  N_min=128 and root/Nmin=16=4^2,
    # so at most three consecutive resolved-pair events: 2048,512,128.
    rows = (_step(2048.0), _step(512.0), _step(128.0))
    out = resolved_pressure_pair_epoch_telescope(
        rows,
        global_energy_upper=10.0,
        pressure_source_weight_floor=1.0,
        scaled_lifetime=1.0,
    )
    assert out.event_count == 3
    assert out.maximum_event_count == 3
    assert out.parent_frequency_floor == pytest.approx(128.0)
    assert out.minimum_observed_pair_mass == pytest.approx(0.2)
    assert out.maximum_observed_pair_entropy > math.log(4.0)
    assert not out.dominant_cut_used
    assert not out.pressure_entropy_used_as_cost
    assert not out.pair_capacity_summed_across_events


def test_fourth_diffuse_pressure_event_fails_below_local_parent_floor():
    rows = (_step(2048.0), _step(512.0), _step(128.0), _step(32.0))
    with pytest.raises(ValueError, match="capacity envelope|frequency floor"):
        resolved_pressure_pair_epoch_telescope(
            rows,
            global_energy_upper=10.0,
            pressure_source_weight_floor=1.0,
            scaled_lifetime=1.0,
        )


def test_dominant_pair_is_included_as_special_case_without_using_cut():
    rows = (_step(512.0, atoms=2), _step(128.0, atoms=2))
    out = resolved_pressure_pair_epoch_telescope(
        rows,
        global_energy_upper=20.0,
        pressure_source_weight_floor=1.0,
        scaled_lifetime=1.0,
    )
    assert out.minimum_observed_pair_mass == pytest.approx(0.5)
    assert not out.dominant_cut_used


def test_sgs_only_pressure_owner_cannot_enter_resolved_pair_epoch():
    with pytest.raises(ValueError, match="did not certify a resolved positive pair owner"):
        certify_resolved_pressure_pair_owner(
            1.0,
            1.0,
            512.0,
            sgs_positive_source_weight=1.0,
            pair_positive_weights=[0.1],
            pair_shell_indices=[(0, 0)],
            pair_frequencies=[(128.0, 128.0)],
        )


def test_nonconsecutive_selected_child_path_is_rejected():
    first = _step(512.0)
    second = _step(64.0)
    with pytest.raises(ValueError, match="selected physical child shell"):
        resolved_pressure_pair_epoch_telescope(
            (first, second),
            global_energy_upper=100.0,
            pressure_source_weight_floor=1.0,
            scaled_lifetime=1.0,
        )


def test_certificate_explicitly_keeps_capacity_noncausal_and_fresh_sgs_open():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "never normalized, charged, summed across events" in cert["capacity_semantics"]
    assert "dominant and diffuse" in cert["entropy"]
    assert "pressure-SGS/fresh-SGS" in cert["scope"]
