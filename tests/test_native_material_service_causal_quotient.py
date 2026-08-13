import numpy as np
import pytest

from src.common_slice_coefficient_registration import registration_first_stop
from src.fresh_service_scale_reentry import fresh_service_scale_route, pushforward_fresh_edges_to_bands
from src.native_material_service_causal_quotient import (
    MATERIAL_FRESH_SERVICE_PROVENANCE,
    MATERIAL_INTERFACE_SERVICE_PROVENANCE,
    RAW_MATERIAL_CAUSE_LABELS,
    smooth_role_subset_flux_from_kphys_relink,
    material_ownership_rereading_anti_theorem,
    positive_material_service_causal_quotient,
    require_native_service_cause_hits,
    require_native_service_owner_labels,
)
from src.physical_branch_compiler import CauseHit, PhysicalCause
from src.smooth_quadratic_carrier_interface import GaugeQuotientedInterfaceWork


def _exact_relink_work(T: np.ndarray) -> GaugeQuotientedInterfaceWork:
    T = np.asarray(T, float)
    relink = T.sum(axis=1)
    return GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=tuple(float(x) for x in relink),
        signed_physical_relink_atoms=tuple(float(x) for x in relink),
        signed_existing_strain_atoms=tuple(0.0 for _ in relink),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=tuple(tuple(float(x) for x in row) for row in T),
    )


def test_kphys_smooth_role_subset_work_is_exact_boundary_flux_not_source():
    T = np.array(
        [
            [0.0, 2.0, -5.0, 1.0],
            [-2.0, 0.0, 3.0, 4.0],
            [5.0, -3.0, 0.0, -7.0],
            [-1.0, -4.0, 7.0, 0.0],
        ]
    )
    out = smooth_role_subset_flux_from_kphys_relink(_exact_relink_work(T), (0, 1))
    assert out.selected_signed_relink_work == pytest.approx(out.signed_boundary_flux_into_selected)
    assert out.signed_boundary_flux_into_selected == pytest.approx(
        out.positive_boundary_inflow - out.positive_boundary_outflow
    )
    assert out.total_relink_work_residual == pytest.approx(0.0)
    assert out.native_relink_strain_split_residual == pytest.approx(0.0)
    assert out.subset_divergence_residual == pytest.approx(0.0)
    assert out.recursive_generation_created is False
    assert out.physical_source_created is False


def test_smooth_role_subset_flux_rejects_interface_certificate_that_breaks_native_relink_strain_split():
    T = np.array([[0.0, 2.0], [-2.0, 0.0]])
    bad = GaugeQuotientedInterfaceWork(
        signed_native_interface_atoms=(3.0, -2.0),
        signed_physical_relink_atoms=(2.0, -2.0),
        signed_existing_strain_atoms=(0.0, 0.0),
        gauge_transport_operator_residual=0.0,
        skew_decomposition_residual=0.0,
        signed_physical_relink_pair_matrix=((0.0, 2.0), (-2.0, 0.0)),
    )
    with pytest.raises(ValueError, match=r"native=relink\+strain"):
        smooth_role_subset_flux_from_kphys_relink(bad, (0,))



def test_internal_kphys_circulation_cancels_inside_smooth_role_subset():
    # The very large 0<->1 circulation is internal to O={0,1}; it must disappear
    # from the subset balance. Only O--O^c pair flux may remain.
    T = np.array(
        [
            [0.0, 1e6, 3.0],
            [-1e6, 0.0, -2.0],
            [-3.0, 2.0, 0.0],
        ]
    )
    out = smooth_role_subset_flux_from_kphys_relink(_exact_relink_work(T), (0, 1))
    assert out.signed_boundary_flux_into_selected == pytest.approx(1.0)
    assert out.selected_signed_relink_work == pytest.approx(1.0)
    assert out.recursive_generation_created is False


def test_oo_on_nn_partition_is_downstream_of_one_positive_service_law():
    weights = np.array([1.0, 2.0, 4.0, 8.0])
    q = positive_material_service_causal_quotient(
        service_measure="positive_sgs_increment_service",
        native_owner=PhysicalCause.RESOLVED_SOURCE.value,
        edge_weights=weights,
        old_here=(True, True, False, False),
        old_neighbor=(True, False, True, False),
    )
    assert q.total_service == pytest.approx(15.0)
    assert q.old_old_service == pytest.approx(1.0)
    assert q.old_new_interface_service == pytest.approx(6.0)
    assert q.new_new_service == pytest.approx(8.0)
    assert q.material_provenance == (
        MATERIAL_INTERFACE_SERVICE_PROVENANCE,
        MATERIAL_FRESH_SERVICE_PROVENANCE,
    )
    assert q.native_owner == PhysicalCause.RESOLVED_SOURCE.value
    assert q.service_created_by_material_partition is False
    assert q.recursive_generation_created_by_material_partition is False


def test_same_service_law_can_have_different_material_partition_without_new_generation():
    weights = np.array([2.0, 3.0, 5.0, 7.0, 11.0, 13.0])
    out = material_ownership_rereading_anti_theorem(
        service_measure="fixed_positive_physical_service",
        native_owner=PhysicalCause.HIGH_STRAIN_DISSIPATION.value,
        edge_weights=weights,
        first_old_here=(True, True, True, False, False, False),
        first_old_neighbor=(False, True, False, True, False, False),
        second_old_here=(True, False, True, False, True, False),
        second_old_neighbor=(True, False, True, False, True, False),
    )
    a = out["first"]
    b = out["second"]
    assert out["same_underlying_service"] is True
    assert out["same_native_owner"] is True
    assert a.total_service == pytest.approx(b.total_service)
    assert out["ownership_partition_l1_change"] > 0.0
    assert out["material_generation_created"] is False


def test_fixed_positive_law_needs_no_invented_supplier_to_prove_material_noncreation():
    q = positive_material_service_causal_quotient(
        service_measure="actual_state_increment_square",
        native_owner=None,
        edge_weights=(1.0, 3.0, 2.0),
        old_here=(True, False, False),
        old_neighbor=(False, False, True),
    )
    assert q.native_owner is None
    assert q.total_service == pytest.approx(6.0)
    assert q.service_created_by_material_partition is False
    assert q.recursive_generation_created_by_material_partition is False


@pytest.mark.parametrize("raw", sorted(RAW_MATERIAL_CAUSE_LABELS))
def test_raw_material_names_are_locators_not_native_service_owners(raw: str):
    with pytest.raises(TypeError, match="not an independent Navier-Stokes generation owner"):
        require_native_service_owner_labels((raw,))


def test_common_slice_material_exit_remains_a_locator_until_native_resolution():
    out = registration_first_stop(
        1.0 + 0.0j,
        1.0 + 0.0j,
        0.0 + 0.0j,
        0.0 + 0.0j,
        material_relink=True,
    )
    assert out["continuing"] is False
    assert out["first_stops"] == (PhysicalCause.MATERIAL_RELINK.value,)
    with pytest.raises(TypeError, match="not an independent Navier-Stokes generation owner"):
        require_native_service_owner_labels(out["first_stops"])


def test_fresh_nn_service_reenters_through_physical_scale_not_material_generation():
    weights = (1.0, 2.0, 1.0, 4.0)
    bands = (0, -1, -2, -1)
    fresh_law = pushforward_fresh_edges_to_bands(
        weights,
        bands,
        old_here=(False, False, False, False),
        old_neighbor=(False, False, False, False),
    )
    route = fresh_service_scale_route(
        integrated_square_service_threshold=32.0,
        scaled_lifetime=1.0,
        block_frequency=8.0,
        fresh_band_weights=fresh_law,
        viscosity=0.03,
    )
    assert route["fresh_service"] == pytest.approx(sum(weights))
    assert route["next_owner"] == "generic_critical_shell_first_stop"
    assert route["master_semantics"] == "RECURSE_CRITICAL_VIA_GENERIC_SHELL"
    assert "whole hard u-shell is not declared fresh material" in route["material_semantics"]



def test_raw_material_causehit_is_rejected_but_native_source_hit_is_preserved():
    with pytest.raises(TypeError, match="carrier/material-state locator"):
        require_native_service_cause_hits((CauseHit(0.2, PhysicalCause.MATERIAL_RELINK),))
    source = CauseHit(0.2, PhysicalCause.RESOLVED_SOURCE, 3.0, "actual SGS/source service")
    assert require_native_service_cause_hits((source,)) == (source,)
