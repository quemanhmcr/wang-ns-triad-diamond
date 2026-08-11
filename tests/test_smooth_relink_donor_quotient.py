import pytest

from src.smooth_quadratic_carrier_interface import RELINK_OWNER, GaugeQuotientedInterfaceWork
from src.smooth_relink_donor_quotient import (
    SMOOTH_RELINK_SAME_EVENT_RELAY,
    STATUS,
    SmoothRelinkDonorCertificate,
    smooth_relink_donor_quotient,
    theorem_certificate,
)


def _pure_relink_work() -> GaugeQuotientedInterfaceWork:
    T = (
        (0.0, 3.0, -1.0),
        (-3.0, 0.0, 2.0),
        (1.0, -2.0, 0.0),
    )
    relink = (2.0, -1.0, -1.0)
    return GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=relink,
        signed_physical_relink_atoms=relink,
        signed_existing_strain_atoms=(0.0, 0.0, 0.0),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=T,
    )


def test_positive_smooth_relink_is_finite_same_event_donor_flux_not_generation():
    out = smooth_relink_donor_quotient(_pure_relink_work())
    assert out["relay_kind"] == SMOOTH_RELINK_SAME_EVENT_RELAY
    assert out["recipient_roles"] == (0,)
    assert set(out["terminal_negative_net_donor_roles"]) == {1, 2}
    assert out["positive_relink_work"] == pytest.approx(2.0)
    assert out["recipient_positive_incoming_flux"] >= out["positive_relink_work"]
    assert out["maximum_shortest_donor_path_length"] <= 2
    assert out["same_physical_event"] is True
    assert out["same_physical_time"] is True
    assert out["new_causal_charge_created"] is False
    assert out["recursive_generation_created"] is False
    assert out["scale_progress_created"] is False
    assert out["hard_smooth_measure_identification_used"] is False
    assert isinstance(out["certificate"], SmoothRelinkDonorCertificate)


def test_smooth_donor_quotient_rejects_unbound_pair_matrix_even_if_atoms_sum_to_zero():
    work = _pure_relink_work()
    bad = GaugeQuotientedInterfaceWork(
        work.signed_native_interface_atoms,
        work.signed_physical_relink_atoms,
        work.signed_existing_strain_atoms,
        work.gauge_transport_operator_residual,
        work.skew_decomposition_residual,
        (
            (0.0, 3.25, -1.0),
            (-3.0, 0.0, 2.0),
            (1.0, -2.0, 0.0),
        ),
    )
    with pytest.raises(ValueError, match="exact conservative K_phys row law"):
        smooth_relink_donor_quotient(bad)


def test_relink_atoms_without_bound_pair_matrix_cannot_be_promoted_to_donor_relay():
    legacy = GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=(2.0, -1.0, -1.0),
        signed_physical_relink_atoms=(2.0, -1.0, -1.0),
        signed_existing_strain_atoms=(0.0, 0.0, 0.0),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
    )
    with pytest.raises(ValueError, match="bound square K_phys pair matrix"):
        smooth_relink_donor_quotient(legacy)


def test_certificate_keeps_smooth_and_hard_measures_distinct_and_closes_only_relink_depth():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "remain distinct physical disintegrations" in cert["measure_separation"]
    assert "physical event/time" in cert["event_semantics"]
    assert "no causal charge, recursive generation, or scale progress" in cert["event_semantics"]
    assert "not a recursive generation owner" in cert["master_boundary"]
    assert "does not terminate genuine strain" in cert["scope"]
    assert RELINK_OWNER in cert["master_boundary"]
