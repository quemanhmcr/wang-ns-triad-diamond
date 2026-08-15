import pytest

from src.descending_fresh_sgs_scale_epoch_telescope import (
    DescendingFreshSGSRenewalStep,
    certify_fresh_sgs_scale_owner,
    descending_fresh_sgs_epoch_telescope,
    fresh_sgs_forced_service_floor_from_objective_stop,
    fresh_sgs_parent_frequency_floor,
    theorem_certificate,
)


def _owner(N: float, *, Y: float = 1.0, c: float = 1.0, j: int = -2):
    # Concentrate the fresh law on one canonical band; this is enough to test
    # the typed scale route.  j=-2 gives candidates N/4 and N/2.
    return certify_fresh_sgs_scale_owner(Y, c, N, {j: Y / 4.0})


def _step(N: float, *, Y: float = 1.0, c: float = 1.0, j: int = -2, choose_upper: bool = True):
    owner = _owner(N, Y=Y, c=c, j=j)
    child = max(owner.hard_shell_candidates) if choose_upper else min(owner.hard_shell_candidates)
    return DescendingFreshSGSRenewalStep(
        owner=owner,
        child_frequency=child,
        child_critical_mass=owner.selected_hard_shell_mass_lower,
    )


def test_fresh_parent_floor_is_Y_over_16cE():
    assert fresh_sgs_parent_frequency_floor(10.0, 8.0, 2.0) == pytest.approx(8.0 / 320.0)


def test_objective_sgs_source_floor_composes_exact_linear_CY():
    A = 2.0
    c = 1.5
    g1 = 1.0
    clp = 2.0
    cb = 3.0
    Y = fresh_sgs_forced_service_floor_from_objective_stop(A, c, g1, clp, cb)
    CY = 380.0 / (g1 * (1.0 + g1) * (clp * cb) ** 2)
    assert Y == pytest.approx(CY * A / (4.0 * c))


def test_strict_half_scale_fresh_epoch_is_finite():
    # Y=1,c=1,E=1/16 gives N_min=1.  Root 8 permits 4 events at parents 8,4,2,1.
    rows = (_step(8.0), _step(4.0), _step(2.0), _step(1.0))
    out = descending_fresh_sgs_epoch_telescope(
        rows,
        global_energy_upper=1.0 / 16.0,
        forced_square_service_floor=1.0,
        scaled_lifetime=1.0,
    )
    assert out.parent_frequency_floor == pytest.approx(1.0)
    assert out.event_count == 4
    assert out.maximum_event_count == 4
    assert out.maximum_observed_scale_ratio == pytest.approx(0.5)


def test_fifth_half_scale_event_fails_below_parent_floor():
    rows = (_step(8.0), _step(4.0), _step(2.0), _step(1.0), _step(0.5))
    with pytest.raises(ValueError, match="energy envelope|frequency floor"):
        descending_fresh_sgs_epoch_telescope(
            rows,
            global_energy_upper=1.0 / 16.0,
            forced_square_service_floor=1.0,
            scaled_lifetime=1.0,
        )


def test_top_band_non_descending_fresh_child_is_rejected_not_scalarized():
    owner = _owner(8.0, j=0)
    child = max(owner.hard_shell_candidates)  # 2N
    with pytest.raises(ValueError, match="strictly descending"):
        DescendingFreshSGSRenewalStep(
            owner=owner,
            child_frequency=child,
            child_critical_mass=owner.selected_hard_shell_mass_lower,
        )


def test_j_minus_one_upper_child_N_is_also_excluded():
    owner = _owner(8.0, j=-1)
    child = max(owner.hard_shell_candidates)  # N
    with pytest.raises(ValueError, match="strictly descending"):
        DescendingFreshSGSRenewalStep(
            owner=owner,
            child_frequency=child,
            child_critical_mass=owner.selected_hard_shell_mass_lower,
        )


def test_scope_names_top_band_survivor_explicitly():
    cert = theorem_certificate()
    assert cert["status"].startswith("DRAFT_")
    assert "j=0 or j=-1" in cert["top_band_survivor"]
    assert "top-band fresh renewal" in cert["scope"]
