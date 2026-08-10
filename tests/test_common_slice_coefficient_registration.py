import math

from src.common_slice_coefficient_registration import (
    CONTINUING_PRODUCT_FRACTION,
    HH_COEFFICIENT_OBSTRUCTION,
    ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
    common_slice_natural_window_margin,
    continuing_pair_product_lower,
    registration_first_stop,
)


def test_async_cone_leaves_clean_common_slice_window_margin():
    assert math.isclose(common_slice_natural_window_margin(), 67.0 / 195.0)


def test_no_source_stop_forces_one_quarter_inheritance():
    z_event = 1.0 + 0.0j
    i_hh = 0.2 + 0.0j
    i_r = 0.1 + 0.0j
    z_slice = z_event - i_hh - i_r
    out = registration_first_stop(z_event, z_slice, i_hh, i_r)
    assert out["branch"] == "registered_material_inheritance"
    assert out["continuing"] is True
    assert float(out["slice_amplitude"]) >= 0.25


def test_large_hh_impulse_is_an_energy_reentry_obstruction_not_work():
    z_event = 1.0 + 0.0j
    i_hh = 0.6 + 0.0j
    i_r = 0j
    out = registration_first_stop(z_event, z_event - i_hh, i_hh, i_r)
    assert out["branch"] == "hh_coefficient_obstruction_stop"
    assert out["continuing"] is False
    assert out["requires_physical_energy_reentry"] is True
    assert out["coefficient_impulses_used_as_work"] is False


def test_large_residual_locates_energy_reentry_before_source_ownership():
    z_event = 1.0 + 0.0j
    i_hh = 0j
    i_r = 0.3 + 0.0j
    out = registration_first_stop(z_event, z_event - i_r, i_hh, i_r)
    assert out["branch"] == "role_interface_coefficient_obstruction_stop"
    assert out["continuing"] is False
    assert out["requires_physical_energy_reentry"] is True
    assert out["coefficient_impulses_used_as_work"] is False


def test_two_continuing_parents_keep_one_sixteenth_product():
    assert math.isclose(CONTINUING_PRODUCT_FRACTION, 1.0 / 16.0)
    assert math.isclose(continuing_pair_product_lower(8.0), 0.5)


def test_simultaneous_first_stops_are_returned_without_primary_priority():
    z_event = 1.0 + 0.0j
    i_hh = 0.6 + 0.0j
    i_r = 0.3 + 0.0j
    z_slice = z_event - i_hh - i_r
    out = registration_first_stop(z_event, z_slice, i_hh, i_r, material_relink=True)
    assert out["branch"] == "multiple_first_stops_before_common_slice"
    assert set(out["first_stops"]) == {
        "material_relink",
        ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,
        HH_COEFFICIENT_OBSTRUCTION,
    }
    assert out["continuing"] is False
    assert out["primary_selected"] is False
