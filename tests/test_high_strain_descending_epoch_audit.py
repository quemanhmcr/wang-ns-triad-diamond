import math

import pytest

from src.high_strain_descending_epoch_telescope import (
    HIGH_STRAIN_RENEWAL_RATIO_UPPER,
    HighStrainEpochCertificate,
    HighStrainRenewalStep,
    _frequency_floor_count_upper,
    high_strain_epoch_telescope,
)
from src.high_strain_dissipation_collision import clean_high_strain_dissipation_lower
from src.nn_critical_heat_carrier_seed import RENEWAL_SCALE_FACTOR, renewal_scale


def _physical_step(
    child_frequency: float,
    normalized_dissipation: float,
    ratio: float = HIGH_STRAIN_RENEWAL_RATIO_UPPER,
) -> HighStrainRenewalStep:
    ancestor = (ratio / RENEWAL_SCALE_FACTOR) * child_frequency
    return HighStrainRenewalStep(
        child_frequency,
        ancestor,
        renewal_scale(ancestor),
        normalized_dissipation,
    )


def test_native_scale_ancestor_cannot_ascend_inside_absolute_unit_tolerance():
    native = 1.0e-120
    ancestor = 2.0 * native
    with pytest.raises(ValueError, match="M<=N/4"):
        HighStrainRenewalStep(
            native,
            ancestor,
            renewal_scale(ancestor),
            native,
        )


def test_native_scale_subthreshold_dissipation_is_not_promoted_to_high_strain():
    scaled_lifetime = 1.0e120
    threshold = clean_high_strain_dissipation_lower(scaled_lifetime)
    row = _physical_step(4.0 * threshold, 0.5 * threshold)
    with pytest.raises(ValueError, match="below the physical high-strain"):
        high_strain_epoch_telescope(
            (row,),
            total_gradient_dissipation=1.0,
            scaled_lifetime=scaled_lifetime,
        )


def test_native_scale_global_gradient_reservoir_cannot_be_overdrawn():
    scaled_lifetime = 1.0e120
    threshold = clean_high_strain_dissipation_lower(scaled_lifetime)
    row = _physical_step(2.0 * threshold, 3.0 * threshold)
    with pytest.raises(ValueError, match="exceeds the supplied global gradient reservoir"):
        high_strain_epoch_telescope(
            (row,),
            total_gradient_dissipation=1.0,
            scaled_lifetime=scaled_lifetime,
        )


def test_native_scale_foreign_child_cannot_cross_consecutive_carrier_binding():
    scaled_lifetime = 1.0e120
    threshold = clean_high_strain_dissipation_lower(scaled_lifetime)
    first = _physical_step(100.0 * threshold, threshold)
    foreign_child = 2.0 * first.renewal_frequency
    second = _physical_step(foreign_child, threshold)
    with pytest.raises(ValueError, match="actual renewed carrier scale consecutively"):
        high_strain_epoch_telescope(
            (first, second),
            total_gradient_dissipation=1.0,
            scaled_lifetime=scaled_lifetime,
        )


def test_epoch_certificate_rejects_nan_observed_scale_ratio():
    with pytest.raises(ValueError, match="scale ratio"):
        HighStrainEpochCertificate(
            step_count=1,
            root_frequency=1.0,
            last_child_frequency=1.0,
            total_gradient_dissipation=1.0,
            high_strain_dissipation_lower=1.0,
            physical_frequency_floor=1.0,
            geometric_frequency_sum_upper=1.25,
            normalized_dissipation_sum=1.0,
            normalized_dissipation_capacity_upper=1.25,
            frequency_floor_count_upper=1,
            weighted_capacity_count_upper=1,
            certified_count_upper=1,
            maximum_observed_scale_ratio=float("nan"),
        )


def test_frequency_floor_count_uses_log_geometry_without_ratio_underflow():
    count = _frequency_floor_count_upper(1.0e300, 1.0e-300)
    expected_coordinate = (
        math.log(1.0e-300) - math.log(1.0e300)
    ) / math.log(HIGH_STRAIN_RENEWAL_RATIO_UPPER)
    assert count == math.floor(expected_coordinate + 8.0e-12) + 1
    assert count > 800
