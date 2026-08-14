import pytest

from src.intrinsic_generated_donor_time_telescope import (
    BACKWARD_SUM_COEFF,
    DONOR_CHILD_RATIO_UPPER,
    HALF_CHILD_TO_DONOR_SPAN_MAX,
    LIFETIME_GROWTH_MIN,
    ONE_STEP_BACKSHIFT_MIN,
    GeneratedEnergyDonorStep,
    generated_energy_donor_epoch_telescope,
    theorem_certificate,
)


def _step(child, donor, start, end, c=1.0):
    return GeneratedEnergyDonorStep(
        child_frequency=child,
        donor_frequency=donor,
        scaled_lifetime=c,
        work_support_start=start,
        work_support_end=end,
        physical_hh_work_mass=1.0,
        physical_hh_work_total=1.5,
        physical_hh_work_lower=1.2,
    )


def test_constants_are_derived_only_from_upper_donor_ratio_and_heavy_half():
    assert DONOR_CHILD_RATIO_UPPER == 5/8
    assert float(LIFETIME_GROWTH_MIN) == pytest.approx(64/25)
    assert HALF_CHILD_TO_DONOR_SPAN_MAX == 25/128
    assert float(ONE_STEP_BACKSHIFT_MIN) == pytest.approx(6859/16000)
    assert float(BACKWARD_SUM_COEFF) == pytest.approx(6859/24960)
    cert=theorem_certificate()
    assert "no lower donor/child ratio" in cert["scope"]
    assert "no log-progress J" in cert["scope"]


def test_one_step_needs_no_three_fifths_lower_ratio():
    row=_step(child=10.0,donor=4.0,start=1.0,end=1.005)
    assert row.donor_child_ratio == pytest.approx(0.4)
    assert row.donor_child_ratio < 3/5
    assert row.normalized_donor_span < float(HALF_CHILD_TO_DONOR_SPAN_MAX)


def test_two_layer_physical_donor_lineage_has_intrinsic_backward_shift():
    # T_d0=1/16.  Put the next support at the right edge of the previous common
    # interval to exercise the worst allowed asynchronous placement.
    first=_step(child=10.0,donor=4.0,start=1.0,end=1.005)
    s0=first.common_reference_time
    second_child=4.0
    second_donor=1.5
    width=0.020
    end=first.work_support_end
    start=end-width
    assert start >= s0
    second=_step(child=second_child,donor=second_donor,start=start,end=end)
    out=generated_energy_donor_epoch_telescope((first,second))
    assert out.layer_count == 2
    assert out.minimum_lifetime_growth > float(LIFETIME_GROWTH_MIN)
    assert out.cumulative_reference_backshift >= out.minimum_cumulative_backshift - 1e-12


def test_ratio_at_or_above_five_eighths_fails_closed():
    with pytest.raises(ValueError,match="5/8"):
        _step(child=8.0,donor=5.0,start=1.0,end=1.001)
