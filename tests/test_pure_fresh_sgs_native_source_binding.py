import pytest

from src.descending_fresh_sgs_scale_epoch_telescope import certify_fresh_sgs_scale_owner
from src.pure_fresh_sgs_native_source_binding import (
    NativePureFreshSGSFirstStopStep,
    bind_fresh_owner_to_native_sgs_source,
    native_pure_fresh_sgs_pre_singular_exhaustion,
    pure_fresh_source_floors_from_objective_action,
    theorem_certificate,
)
from src.pure_fresh_sgs_pre_singular_exhaustion import FixedSmoothSGSFilterEnvelope

CB = 190.0**0.5  # with g1=C_LP=1 this makes C_Y=1


def _bound_owner(N: float, j: int, sigma: float = 1.0):
    owner = certify_fresh_sgs_scale_owner(sigma, 1.0, N, {j: sigma / 4.0})
    return bind_fresh_owner_to_native_sgs_source(
        owner, sigma, filter_l1=1.0, lp_constant=1.0, bernstein_constant=CB
    )


def test_objective_action_supplies_matched_uniform_source_and_service_floors():
    sigma, Y = pure_fresh_source_floors_from_objective_action(
        4.0, 2.0, filter_l1=1.0, lp_constant=1.0, bernstein_constant=CB
    )
    assert sigma == pytest.approx(0.5)
    assert Y == pytest.approx(0.5)


def test_unrelated_source_scalar_fails_closed():
    owner = certify_fresh_sgs_scale_owner(1.0, 1.0, 4.0, {0: 0.25})
    with pytest.raises(ValueError, match="not the canonical image"):
        bind_fresh_owner_to_native_sgs_source(
            owner, 2.0, filter_l1=1.0, lp_constant=1.0, bernstein_constant=CB
        )


def test_native_wrapper_routes_only_bound_source_service_steps():
    a_bind = _bound_owner(4.0, -1)
    b_bind = _bound_owner(2.0, 0)
    a = NativePureFreshSGSFirstStopStep(
        a_bind, 2.0, a_bind.owner.selected_hard_shell_mass_lower, 100.0, 200.0
    )
    b = NativePureFreshSGSFirstStopStep(
        b_bind, 2.0, b_bind.owner.selected_hard_shell_mass_lower, 0.0, 100.0
    )
    out = native_pure_fresh_sgs_pre_singular_exhaustion(
        (a, b),
        global_energy_upper=8.0,
        pre_singular_h1_seminorm_sq_upper=100.0,
        forced_square_service_floor=1.0,
        sgs_source_weight_floor=1.0,
        scaled_lifetime_upper=1.0,
        filter_envelope=FixedSmoothSGSFilterEnvelope(1.0, 1.0),
    )
    assert out.event_count == 2


def test_certificate_declares_fail_closed_binding():
    cert = theorem_certificate()
    assert "Y=C_Y Sigma_R" in cert["single_law"]
    assert "unrelated" in cert["fail_closed"]
