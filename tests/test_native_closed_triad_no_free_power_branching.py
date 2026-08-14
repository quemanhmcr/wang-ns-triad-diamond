import math

import pytest

from src.cyclic_helical_triad_donor_kernel import signed_good_integer_triad
from src.native_closed_triad_no_free_power_branching import (
    POWER_CONTINUATION_UPPER,
    POWER_LOG_COST_LOWER,
    SIDE_DIVERSION_LOWER,
    STATUS,
    closed_triad_current_law,
    theorem_certificate,
    weighted_signed_good_restriction,
)


def test_one_scalar_closed_triad_current_has_both_native_null_laws():
    out = closed_triad_current_law((5.0, 3.0, 4.0), (1, -1, 1), 2.75)
    assert math.fsum(out.works) == pytest.approx(0.0, abs=1e-13)
    assert math.fsum(a * t for a, t in zip(out.signed_frequencies, out.works)) == pytest.approx(0.0, abs=1e-13)
    assert not out.temporal_matching_used


def test_same_weight_restriction_preserves_cyclic_good_branching_and_single_charge():
    _triad, cert = signed_good_integer_triad()
    out = weighted_signed_good_restriction(
        (cert, cert, cert, cert),
        (0.125, 3.0, 0.0, 7.25),
    )
    assert out.donor_work == pytest.approx(out.good_work + out.side_work, rel=6e-10)
    assert out.continuation_ratio < float(POWER_CONTINUATION_UPPER)
    assert out.side_diversion_ratio > float(SIDE_DIVERSION_LOWER)
    assert out.local_log_cost > POWER_LOG_COST_LOWER
    assert out.side_is_existing_transfer_work_loss
    assert not out.canonical_good_cause_replaced
    assert not out.side_separately_charged
    assert not out.between_time_stock_claimed
    assert not out.temporal_matching_used
    assert not out.later_hahn_used
    assert not out.master_composition_certified


def test_zero_weight_family_cannot_mint_a_power_event():
    _triad, cert = signed_good_integer_triad()
    with pytest.raises(ValueError):
        weighted_signed_good_restriction((cert,), (0.0,))


def test_theorem_certificate_is_explicitly_local_not_master_closed():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "10 D_q/13" in cert["same_weight_radon_restriction"]
    assert "3 D_q/13" in cert["same_weight_radon_restriction"]
    assert "not between-time stock" in cert["forbidden_interpretation"]
    assert "still requires a theorem" in cert["open_master_seam"]
    assert cert["global_regularity_claimed"] is False
