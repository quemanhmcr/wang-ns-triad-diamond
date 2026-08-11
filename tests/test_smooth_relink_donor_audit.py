import pytest

from src.continuum_master_event_quotient import energy_reentry_master_route
from src.smooth_quadratic_carrier_interface import (
    RELINK_OWNER,
    STRAIN_OWNER,
    GaugeQuotientedInterfaceWork,
    positive_smooth_interface_split,
)
from src.smooth_relink_donor_quotient import (
    SmoothRelinkDonorCertificate,
    smooth_relink_donor_quotient,
)


BASE_PAIR_MATRIX = (
    (0.0, 3.0, -1.0),
    (-3.0, 0.0, 2.0),
    (1.0, -2.0, 0.0),
)


def _pure_relink_work(scale: float, *, pair_scale: float = 1.0) -> GaugeQuotientedInterfaceWork:
    relink = tuple(scale * value for value in (2.0, -1.0, -1.0))
    pair = tuple(
        tuple(scale * pair_scale * value for value in row)
        for row in BASE_PAIR_MATRIX
    )
    return GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=relink,
        signed_physical_relink_atoms=relink,
        signed_existing_strain_atoms=(0.0, 0.0, 0.0),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=pair,
    )


def _smooth_reentry(
    work: GaugeQuotientedInterfaceWork,
    owners: tuple[str, ...],
) -> dict[str, object]:
    return {
        "branch": "smooth_interface_physical_work",
        "joint_interface_owners": owners,
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
        "gauge_quotiented_interface_work_certificate": work,
    }


def test_smooth_relink_donor_law_is_covariant_at_tiny_native_work_scale():
    reference = smooth_relink_donor_quotient(_pure_relink_work(1.0))
    tiny = smooth_relink_donor_quotient(_pure_relink_work(1.0e-120))
    huge = smooth_relink_donor_quotient(_pure_relink_work(1.0e120))

    for rescaled in (tiny, huge):
        assert rescaled["recipient_roles"] == reference["recipient_roles"] == (0,)
        assert rescaled["terminal_negative_net_donor_roles"] == reference[
            "terminal_negative_net_donor_roles"
        ]
        assert rescaled["maximum_shortest_donor_path_length"] == reference[
            "maximum_shortest_donor_path_length"
        ]
    assert positive_smooth_interface_split(_pure_relink_work(1.0e-120))[
        "joint_physical_owners"
    ] == (RELINK_OWNER,)


def test_native_scale_pair_flux_deficit_cannot_hide_behind_an_observer_unit_floor():
    # The relink rows are twice the row sums of the supplied pair law.  At this
    # native scale the old max(1, ...) tolerances accepted both the false binding
    # and a 25 percent deficit in actual incoming donor flux.
    unbound = _pure_relink_work(1.0e-12, pair_scale=0.5)
    with pytest.raises(ValueError, match="row law|incoming donor flux"):
        smooth_relink_donor_quotient(unbound)


def test_master_replays_physical_split_instead_of_trusting_claimed_owner_labels():
    work = _pure_relink_work(1.0)
    forged = _smooth_reentry(work, (STRAIN_OWNER,))
    with pytest.raises(TypeError, match="owner|replay|relink"):
        energy_reentry_master_route(
            "positive native smooth interface work",
            2.0,
            forged,
        )


def test_master_binds_supplied_mass_to_actual_positive_native_interface_work():
    work = _pure_relink_work(1.0)
    with pytest.raises(ValueError, match="mass|positive native"):
        energy_reentry_master_route(
            "positive native smooth interface work",
            200.0,
            _smooth_reentry(work, (RELINK_OWNER,)),
        )


@pytest.mark.parametrize("bad_mass", [0.0, -1.0, float("nan"), float("inf")])
def test_master_rejects_nonphysical_smooth_interface_route_mass(bad_mass: float):
    work = _pure_relink_work(1.0)
    with pytest.raises(ValueError, match="mass|positive|finite"):
        energy_reentry_master_route(
            "positive native smooth interface work",
            bad_mass,
            _smooth_reentry(work, (RELINK_OWNER,)),
        )


def test_mixed_relink_strain_routes_actual_strain_component_not_full_native_mass():
    relink = (2.0, -1.0, -1.0)
    strain = (2.0, 0.0, 0.0)
    native = tuple(a + b for a, b in zip(relink, strain, strict=True))
    work = GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=native,
        signed_physical_relink_atoms=relink,
        signed_existing_strain_atoms=strain,
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=BASE_PAIR_MATRIX,
    )
    route = energy_reentry_master_route(
        "positive native smooth interface work",
        4.0,
        _smooth_reentry(work, (RELINK_OWNER, STRAIN_OWNER)),
    )

    assert route.owner_bundle is not None
    assert route.owner_bundle.owners == (STRAIN_OWNER,)
    assert route.owner_bundle.physical_measure == "positive existing smooth strain work"
    assert route.owner_bundle.mass == pytest.approx(2.0)
    assert route.mass == pytest.approx(4.0)


def test_typed_donor_certificate_rejects_residuals_large_in_its_native_scale():
    with pytest.raises(ValueError, match="residual|native"):
        SmoothRelinkDonorCertificate(
            relink_owner=RELINK_OWNER,
            recipient_roles=(0,),
            terminal_negative_net_donor_roles=(1,),
            maximum_shortest_donor_path_length=1,
            role_count=2,
            positive_relink_work=1.0e-30,
            recipient_positive_incoming_flux=1.0e-30,
            pair_antisymmetry_residual=1.0e-12,
            row_binding_residual=0.0,
            total_relink_work_residual=0.0,
        )


@pytest.mark.parametrize(
    ("residual_field", "bad_value"),
    [
        ("gauge_transport_operator_residual", float("nan")),
        ("gauge_transport_operator_residual", float("inf")),
        ("gauge_transport_operator_residual", -1.0),
        ("skew_decomposition_residual", float("nan")),
        ("skew_decomposition_residual", float("inf")),
        ("skew_decomposition_residual", -1.0),
    ],
)
def test_gauge_quotiented_work_rejects_invalid_provenance_residuals(
    residual_field: str,
    bad_value: float,
):
    values = {
        "signed_native_interface_atoms": (2.0, -1.0, -1.0),
        "signed_physical_relink_atoms": (2.0, -1.0, -1.0),
        "signed_existing_strain_atoms": (0.0, 0.0, 0.0),
        "gauge_transport_operator_residual": 0.0,
        "skew_decomposition_residual": 0.0,
        "signed_physical_relink_pair_matrix": BASE_PAIR_MATRIX,
    }
    values[residual_field] = bad_value
    with pytest.raises(ValueError, match="gauge|skew|residual|finite|nonnegative"):
        GaugeQuotientedInterfaceWork(**values)
