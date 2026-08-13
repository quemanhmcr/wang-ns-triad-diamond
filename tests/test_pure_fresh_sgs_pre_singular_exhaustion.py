import pytest

from src.descending_fresh_sgs_scale_epoch_telescope import certify_fresh_sgs_scale_owner
from src.pure_fresh_sgs_pre_singular_exhaustion import (
    FixedSmoothSGSFilterEnvelope,
    PureFreshSGSFirstStopStep,
    doubling_fresh_hard_shell_mass_lower,
    fresh_doubling_is_excluded,
    fresh_scale_active_band_count_upper,
    fresh_scale_exponential_moment_upper_from_owner_floor,
    fresh_sgs_physical_slab_duration_lower,
    max_scale_atom_lower_from_exponential_moment,
    pure_fresh_sgs_pre_singular_exhaustion,
    resolved_sgs_source_rate_energy_upper,
    theorem_certificate,
)


def _filter():
    # The theorem is parametric in the two fixed smooth-kernel norms.  Unit values
    # are a regression normalization, not a claim about the repository cutoff.
    return FixedSmoothSGSFilterEnvelope(1.0, 1.0)


def _owner(N: float, j: int, *, Y: float = 1.0, c: float = 1.0):
    return certify_fresh_sgs_scale_owner(Y, c, N, {j: Y / 4.0})


def _step(N: float, j: int, child: float, earlier: float, later: float, *, sigma: float = 1.0):
    owner = _owner(N, j)
    return PureFreshSGSFirstStopStep(
        owner=owner,
        child_frequency=child,
        child_critical_mass=owner.selected_hard_shell_mass_lower,
        sgs_source_weight=sigma,
        earlier_time=earlier,
        later_time=later,
    )


def test_weighted_moment_forces_a_top_scale_atom_without_entropy_cost():
    K = fresh_scale_exponential_moment_upper_from_owner_floor(
        8.0,
        global_energy_upper=1.0,
        forced_square_service_floor=2.0,
        scaled_lifetime_upper=0.5,
    )
    assert K == pytest.approx(32.0)
    assert fresh_scale_active_band_count_upper(K) == 6
    assert max_scale_atom_lower_from_exponential_moment(K) == pytest.approx(1.0 / 12.0)


def test_doubling_lower_eventually_beats_compact_pre_singular_h1_upper():
    low = doubling_fresh_hard_shell_mass_lower(
        2.0**30,
        global_energy_upper=1.0,
        forced_square_service_floor=1.0,
        scaled_lifetime_upper=1.0,
    )
    assert low > 0.0
    assert fresh_doubling_is_excluded(
        2.0**30,
        global_energy_upper=1.0,
        forced_square_service_floor=1.0,
        scaled_lifetime_upper=1.0,
        pre_singular_h1_seminorm_sq_upper=1.0,
    )


def test_fixed_filter_energy_envelope_gives_positive_source_time_floor():
    C = resolved_sgs_source_rate_energy_upper(4.0, 2.0, _filter())
    assert C > 0.0
    dt = fresh_sgs_physical_slab_duration_lower(
        4.0,
        sgs_source_weight_floor=0.5,
        global_energy_upper=2.0,
        filter_envelope=_filter(),
    )
    assert dt > 0.0


def test_mixed_down_same_up_fresh_word_is_physically_exhausted_on_smooth_interval():
    # Use deliberately long slabs so only topology/scale/source guards are tested.
    # Geometry: 4 -> 2 (j=-1 lower), 2 -> 2 (j=0 lower), 2 -> 4 (j=0 upper).
    a = _step(4.0, -1, 2.0, 2.0, 3.0)
    b = _step(2.0, 0, 2.0, 1.0, 2.0)
    c = _step(2.0, 0, 4.0, 0.0, 1.0)
    out = pure_fresh_sgs_pre_singular_exhaustion(
        (a, b, c),
        global_energy_upper=8.0,
        pre_singular_h1_seminorm_sq_upper=100.0,
        forced_square_service_floor=1.0,
        sgs_source_weight_floor=1.0,
        scaled_lifetime_upper=1.0,
        filter_envelope=_filter(),
    )
    assert out.event_count == 3
    assert out.maximum_event_count_from_physical_time >= 3
    assert out.reachable_frequency_upper >= out.maximum_parent_frequency


def test_nonconsecutive_observer_bins_cannot_masquerade_as_recursive_first_stop_word():
    a = _step(4.0, -1, 2.0, 2.0, 3.0)
    b = _step(2.0, 0, 2.0, 0.5, 1.5)
    with pytest.raises(ValueError, match="meet at their physical endpoint"):
        pure_fresh_sgs_pre_singular_exhaustion(
            (a, b),
            global_energy_upper=8.0,
            pre_singular_h1_seminorm_sq_upper=100.0,
            forced_square_service_floor=1.0,
            sgs_source_weight_floor=1.0,
            scaled_lifetime_upper=1.0,
            filter_envelope=_filter(),
        )


def test_certificate_keeps_h1_local_and_entropy_noncausal():
    cert = theorem_certificate()
    assert "fixed compact pre-singular" in cert["doubling_exclusion"]
    assert "not charged" in cert["deterministic_concentration"]
    assert "no global H1" in cert["forbidden_shortcuts"]
