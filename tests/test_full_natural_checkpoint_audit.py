"""Adversarial audit of the full-natural checkpoint boundary.

These tests target native PDE provenance.  A no-hit corridor may be quotiented as
an event vertex, but its actual parent shell, carrier scale, elapsed parabolic
time, endpoint shell state, and successor registration may not be rebound by an
observer merely because some floating-point numbers are close.
"""

import math

import pytest

from src.full_natural_checkpoint_quotient import (
    FullNaturalCheckpoint,
    checkpoint_chain_ledger,
    checkpoint_from_full_natural_outcome,
    checkpoint_reregistration,
    checkpoint_transition_from_full_natural_outcome,
    geometric_uv_checkpoint_time,
)
from src.full_natural_service_corridor_quotient import (
    FULL_NATURAL_SERVICE_WITNESS,
    RENEWAL_TO_PARENT_SHELL_RATIO,
    TWO_HARD_SHELL_COVER_FACTOR,
)


def _native_duration(A: float, c: float) -> float:
    """Evaluate c/A^2 without first forming the potentially overflowing A^2."""
    return (math.sqrt(c) / A) ** 2


def _full_outcome(
    M: float,
    c: float,
    *,
    mu: float = 2.0,
    terminal_time: float | None = None,
    physical_time_drop: float | None = None,
    parent_shell_critical_mass_lower: float | None = None,
) -> dict[str, object]:
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    natural = _native_duration(A, c)
    drop = natural if physical_time_drop is None else float(physical_time_drop)
    t = 8.0 * natural if terminal_time is None else float(terminal_time)
    parent_mass = mu if parent_shell_critical_mass_lower is None else float(parent_shell_critical_mass_lower)
    return {
        "classification": FULL_NATURAL_SERVICE_WITNESS,
        "joint_first_stops": (),
        "required_elapsed": drop,
        "observed_elapsed_end": drop,
        "corridor_terminal_time": t,
        "corridor_endpoint_time": t - drop,
        "corridor_endpoint_elapsed_from_terminal": drop,
        "physical_time_drop": drop,
        "renewal_frequency": A,
        "scaled_lifetime": c,
        "parent_shell_frequency": M,
        "parent_shell_critical_mass_lower": parent_mass,
        "service_same_corridor_witness": True,
        "service_adds_recursion_depth": False,
        "uniform_square_service_lower": 0.2,
        "integrated_bounded_heat_service_lower": 0.2 * c,
        "endpoint_carrier_critical_mass_lower": mu,
        "requires_physical_energy_reentry": False,
        "coefficient_impulses_used_as_work": False,
    }


def test_checkpoint_rejects_half_length_uv_corridor_in_native_units():
    M, c = 1.0e8, 1.0
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    natural = _native_duration(A, c)
    bad = _full_outcome(M, c, physical_time_drop=0.5 * natural)

    with pytest.raises(ValueError, match="corridor|elapsed|provenance|physical|interval"):
        checkpoint_from_full_natural_outcome(
            bad,
            parent_shell_frequency=M,
            scaled_lifetime=c,
        )


def test_checkpoint_cannot_rebind_parent_and_lifetime_at_the_same_clock_length():
    native_M, native_c = 8.0, 1.0
    outcome = _full_outcome(native_M, native_c)

    # c/A^2 is unchanged under (M,c) -> (2M,4c), so a clock-only adapter
    # cannot notice that it has attached the PDE state to a foreign shell.
    with pytest.raises(ValueError, match="parent|renewal|lifetime|provenance"):
        checkpoint_from_full_natural_outcome(
            outcome,
            parent_shell_frequency=2.0 * native_M,
            scaled_lifetime=4.0 * native_c,
        )


def test_checkpoint_replays_the_complete_monitor_horizon_instead_of_flags():
    M, c = 8.0, 1.0
    bad = _full_outcome(M, c)
    bad["observed_elapsed_end"] = 0.5 * float(bad["required_elapsed"])

    with pytest.raises(ValueError, match="horizon|incomplete|corridor"):
        checkpoint_from_full_natural_outcome(
            bad,
            parent_shell_frequency=M,
            scaled_lifetime=c,
        )


def test_deep_uv_checkpoint_keeps_local_elapsed_when_global_endpoint_rounds():
    M, c = 1.0e10, 1.0
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    natural = _native_duration(A, c)
    outcome = _full_outcome(M, c, terminal_time=1.0)
    assert float(outcome["corridor_endpoint_time"]) == float(outcome["corridor_terminal_time"])

    checkpoint = checkpoint_from_full_natural_outcome(
        outcome,
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    assert checkpoint.physical_time_drop == pytest.approx(natural, rel=2e-12, abs=0.0)
    assert checkpoint.endpoint_time == checkpoint.terminal_time


def test_deep_uv_typed_chain_telescopes_when_all_global_clocks_round_together():
    M, c = 1.0e10, 1.0
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c, terminal_time=1.0),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    endpoint_masses = (2.0, 0.8)
    successor_outcome = _full_outcome(
        0.75 * M,
        c,
        terminal_time=first.endpoint_time,
        parent_shell_critical_mass_lower=endpoint_masses[0],
    )
    transition = checkpoint_transition_from_full_natural_outcome(
        first,
        endpoint_masses,
        successor_outcome,
    )
    ledger = checkpoint_chain_ledger((transition,))

    assert transition.successor_checkpoint.endpoint_time == first.terminal_time
    assert ledger["physical_time_drop"] > 0.0
    assert ledger["endpoint_time_drop"] == 0.0
    assert ledger["absolute_clock_residual_diagnostic"] == pytest.approx(
        ledger["physical_time_drop"], rel=2e-12, abs=0.0
    )
    assert ledger["time_telescope_residual"] == 0.0


def test_tiny_parent_scale_cannot_accept_a_foreign_corridor_scale():
    M, c = 1.0e-13, 1.0
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    foreign_A = 2.0 * A
    drop = _native_duration(foreign_A, c)

    with pytest.raises(ValueError, match="corridor scale|renewal scale"):
        FullNaturalCheckpoint(
            terminal_time=4.0 * drop,
            physical_time_drop=drop,
            parent_shell_frequency=M,
            parent_shell_critical_mass_lower=2.0,
            corridor_frequency=foreign_A,
            scaled_lifetime=c,
            endpoint_carrier_critical_mass_lower=2.0,
            endpoint_shell_candidates=(A, 2.0 * A),
        )


def test_tiny_parent_scale_cannot_certify_zero_endpoint_shells():
    M, c = 1.0e-13, 1.0
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    drop = _native_duration(A, c)

    with pytest.raises(ValueError, match="candidate|shell|positive"):
        FullNaturalCheckpoint(
            terminal_time=4.0 * drop,
            physical_time_drop=drop,
            parent_shell_frequency=M,
            parent_shell_critical_mass_lower=2.0,
            corridor_frequency=A,
            scaled_lifetime=c,
            endpoint_carrier_critical_mass_lower=2.0,
            endpoint_shell_candidates=(0.0, 0.0),
        )


def test_tiny_carrier_does_not_promote_zero_mass_as_a_joint_witness():
    M, c, mu = 8.0, 1.0, 1.0e-20
    checkpoint = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c, mu=mu),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    lower = TWO_HARD_SHELL_COVER_FACTOR * mu
    reread = checkpoint_reregistration(checkpoint, (lower, 0.0))

    assert reread["joint_endpoint_witness_frequencies"] == pytest.approx(
        (RENEWAL_TO_PARENT_SHELL_RATIO * M,)
    )
    assert reread["joint_endpoint_witness_critical_masses"] == pytest.approx((lower,))


def test_chain_rejects_an_unrelated_contiguous_parent_shell():
    M, c = 8.0, 1.0
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    unrelated_M = 100.0
    second = checkpoint_from_full_natural_outcome(
        _full_outcome(unrelated_M, c, terminal_time=first.endpoint_time),
        parent_shell_frequency=unrelated_M,
        scaled_lifetime=c,
    )

    with pytest.raises(ValueError, match="scale|witness|transition|provenance"):
        checkpoint_chain_ledger((first, second))


def test_chain_cannot_take_the_losing_cover_branch_without_its_actual_mass():
    M, c = 8.0, 1.0
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    reread = checkpoint_reregistration(first, (2.0, 0.8))
    assert reread["joint_endpoint_witness_frequencies"] == pytest.approx((0.75 * M,))

    losing_M = 1.5 * M
    second = checkpoint_from_full_natural_outcome(
        _full_outcome(losing_M, c, terminal_time=first.endpoint_time),
        parent_shell_frequency=losing_M,
        scaled_lifetime=c,
    )

    # The raw checkpoint pair contains no typed record carrying `reread` and its
    # winning mass into the next producer, so it must fail closed.
    with pytest.raises(ValueError, match="witness|transition|mass|provenance"):
        checkpoint_chain_ledger((first, second))


def test_typed_transition_reuses_the_state_selected_frequency_and_mass():
    M, c = 8.0, 1.0
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    endpoint_masses = (2.0, 0.8)
    winner_M, winner_mass = 0.75 * M, endpoint_masses[0]
    successor_outcome = _full_outcome(
        winner_M,
        c,
        terminal_time=first.endpoint_time,
        parent_shell_critical_mass_lower=winner_mass,
    )

    transition = checkpoint_transition_from_full_natural_outcome(
        first,
        endpoint_masses,
        successor_outcome,
    )
    assert transition.successor_checkpoint.parent_shell_frequency == pytest.approx(winner_M)
    assert transition.successor_checkpoint.parent_shell_critical_mass_lower == pytest.approx(winner_mass)
    ledger = checkpoint_chain_ledger((transition,))
    assert ledger["checkpoints"] == 2
    assert ledger["certified_transitions"] == 1
    assert ledger["recursive_events_added"] == 0


def test_typed_transition_rejects_a_rebound_winner_mass():
    M, c = 8.0, 1.0
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    endpoint_masses = (2.0, 0.8)
    bad_successor = _full_outcome(
        0.75 * M,
        c,
        terminal_time=first.endpoint_time,
        parent_shell_critical_mass_lower=1.0,
    )

    with pytest.raises(ValueError, match="frequency and mass|witness.*mass"):
        checkpoint_transition_from_full_natural_outcome(
            first,
            endpoint_masses,
            bad_successor,
        )


def test_typed_transition_rejects_a_foreign_endpoint_time_token():
    M, c = 8.0, 1.0
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    endpoint_masses = (2.0, 0.8)
    foreign_time = math.nextafter(first.endpoint_time, math.inf)
    bad_successor = _full_outcome(
        0.75 * M,
        c,
        terminal_time=foreign_time,
        parent_shell_critical_mass_lower=endpoint_masses[0],
    )

    with pytest.raises(ValueError, match="endpoint time token"):
        checkpoint_transition_from_full_natural_outcome(
            first,
            endpoint_masses,
            bad_successor,
        )


def test_chain_contiguity_has_no_one_unit_absolute_time_floor():
    M, c = 1.0e8, 1.0
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    first_drop = _native_duration(A, c)
    first = checkpoint_from_full_natural_outcome(
        _full_outcome(M, c, terminal_time=100.0 * first_drop),
        parent_shell_frequency=M,
        scaled_lifetime=c,
    )
    next_M = A
    mismatched_terminal = first.endpoint_time + 0.5 * first_drop
    second = checkpoint_from_full_natural_outcome(
        _full_outcome(next_M, c, terminal_time=mismatched_terminal),
        parent_shell_frequency=next_M,
        scaled_lifetime=c,
    )

    with pytest.raises(ValueError, match="contiguous|time|transition"):
        checkpoint_chain_ledger((first, second))


def test_geometric_uv_clock_stays_positive_when_A_squared_would_overflow():
    M, c, ratio = 1.0e160, 1.0, 1.5
    A = RENEWAL_TO_PARENT_SHELL_RATIO * M
    first = _native_duration(A, c)
    expected = first / (-math.expm1(-2.0 * math.log(ratio)))
    actual = geometric_uv_checkpoint_time(M, c, ratio)

    assert first > 0.0
    assert actual > 0.0
    assert actual == pytest.approx(expected, rel=5e-4, abs=0.0)
