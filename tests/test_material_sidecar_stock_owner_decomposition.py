from dataclasses import replace

import pytest

from src.material_label_carrier_quotient import (
    MATERIAL_MEMBERSHIP_EVENT,
    SELECTED_FAMILY_EVENT,
    carrier_registration_with_material_sidecars,
    selected_family_switch_sidecar,
)
from src.material_sidecar_stock_owner_decomposition import (
    MEMBERSHIP_PROVENANCE_CURRENCY,
    SELECTED_FAMILY_MOYAL_CURRENCY,
    MaterialSidecarCharge,
    material_sidecar_stock_decomposition,
    same_state_selected_family_switch_anti_theorem,
    selected_family_service_no_escape_binding,
    theorem_certificate,
)
from src.same_carrier_checkpoint_segmentation_quotient import partition_same_carrier_path
from src.same_carrier_inherited_energy_relay import same_carrier_inherited_energy_relay


def _segments():
    return partition_same_carrier_path(
        carrier_id="sidecar-stock",
        terminal_amplitude=1.0,
        elapsed_times=(0.0, 0.1, 0.2, 0.3, 0.4),
        strain_action=(0.0, 0.004, 0.008, 0.012, 0.016),
        residual_impulse_abs=(0.0, 0.03, 0.05, 0.04, 0.06),
        hh_impulse_abs=(0.0, 0.08, 0.12, 0.10, 0.16),
        checkpoint_indices=(2,),
    )


def _stock(material):
    return same_carrier_inherited_energy_relay(
        _segments(),
        initial_time=1.6,
        terminal_time=2.0,
        initial_energy=0.31,
        terminal_energy=1.0,
        residual_positive_work=0.12,
        strain_action=0.016,
        material_registration=material,
        initial_endpoint_is_non_event_carrier_slice=True,
    )


def _registration(*, membership, family, switch_energy=0.0):
    z = 1.0 + 0.0j
    ih = 0.08 + 0.02j
    ii = 0.03 - 0.01j
    return carrier_registration_with_material_sidecars(
        z,
        z - ih - ii,
        ih,
        ii,
        intrinsic_material_membership_change=membership,
        selected_family_change=family,
        selected_family_switch_energy=switch_energy,
        same_smooth_role=True,
        same_analysis_probe=True,
    )


def test_certificate_keeps_stock_membership_and_selected_family_in_different_currencies():
    cert = theorem_certificate()
    assert "zero sidecar energy" in cert["membership"]
    assert "exact Moyal" in cert["selected_family"]
    assert "not smooth K_phys" in cert["service_binding"]
    assert "never added to dW" in cert["currency_separation"]
    assert "zero generation depth" in cert["recurrence_scope"]
    assert "separate physical service/source evidence" in cert["recurrence_scope"]
    assert cert["later_hahn_used"] is False
    assert cert["claims_global_regularity"] is False


def test_membership_only_inherited_stock_is_one_stock_charge_plus_zero_charge_provenance():
    stock = _stock(_registration(membership=True, family=False))
    assert stock.selected_family_switch_energy == 0.0
    out = material_sidecar_stock_decomposition(stock)
    assert out.inherited_stock_mass == pytest.approx(stock.initial_energy)
    assert out.stock_charge_count == 1
    assert out.membership_only_or_no_sidecar
    assert not out.requires_selected_family_ancestry_routing
    assert len(out.charges) == 1
    charge = out.charges[0]
    assert charge.event == MATERIAL_MEMBERSHIP_EVENT
    assert charge.currency == MEMBERSHIP_PROVENANCE_CURRENCY
    assert charge.charge == 0.0
    assert charge.recursion_classification == "zero_charge_provenance"


def test_selected_family_switch_binds_stored_charge_to_exact_moyal_symmetric_difference_energy():
    energies = [0.4, 0.9, 1.7, 0.2]
    switch = selected_family_switch_sidecar(energies, {0, 1}, {1, 2})
    R = switch["symmetric_difference_energy"]
    stock = _stock(_registration(membership=False, family=True, switch_energy=R))
    assert stock.material_sidecars == (SELECTED_FAMILY_EVENT,)
    assert stock.selected_family_switch_energy == pytest.approx(R)
    out = material_sidecar_stock_decomposition(stock, selected_family_switch_certificate=switch)
    assert out.requires_selected_family_ancestry_routing
    assert out.inherited_stock_mass == pytest.approx(0.31)
    assert len(out.charges) == 1
    charge = out.charges[0]
    assert charge.event == SELECTED_FAMILY_EVENT
    assert charge.currency == SELECTED_FAMILY_MOYAL_CURRENCY
    assert charge.charge == pytest.approx(R)
    assert charge.recursion_classification == "zero_generation_depth_selected_service_boundary"
    assert not charge.physical_work
    assert not charge.carrier_stock
    assert not charge.smooth_k_phys_relink


def test_selected_family_without_exact_switch_certificate_and_charge_transplant_both_fail_closed():
    energies = [0.5, 1.0, 2.0]
    switch = selected_family_switch_sidecar(energies, {0}, {1})
    R = switch["symmetric_difference_energy"]
    stock = _stock(_registration(membership=True, family=True, switch_energy=R))
    with pytest.raises(TypeError, match="requires exact Moyal"):
        material_sidecar_stock_decomposition(stock)

    other = selected_family_switch_sidecar([1.0, 4.0, 8.0], {0}, {2})
    with pytest.raises(TypeError, match="lost the exact selected-family Moyal charge"):
        material_sidecar_stock_decomposition(stock, selected_family_switch_certificate=other)


def test_positive_moyal_switch_charge_can_exist_with_identical_coherent_state_and_zero_increments():
    energies = [0.4, 1.1, 2.3, 0.7]
    anti = same_state_selected_family_switch_anti_theorem(energies, {0, 1}, {1, 2, 3})
    assert anti.symmetric_difference_energy > 0.0
    assert anti.coherent_ledger_relink_energy == pytest.approx(anti.symmetric_difference_energy)
    assert anti.positive_increment_work == 0.0
    assert anti.negative_increment_work == 0.0
    assert anti.identical_state_energy_residual == 0.0
    assert not anti.generation_event_inferred
    assert not anti.physical_work_inferred


def test_exact_selected_cell_service_law_reads_same_R_switch_but_not_as_Kphys_or_work():
    energies = [0.3, 0.8, 1.5]
    switch = selected_family_switch_sidecar(energies, {0}, {1, 2})
    R = switch["symmetric_difference_energy"]
    stock = _stock(_registration(membership=False, family=True, switch_energy=R))
    out = material_sidecar_stock_decomposition(stock, selected_family_switch_certificate=switch)

    # Saturate P+ = E_final + P- + R with R as the unique >= P+/3 branch.
    pminus = 0.1 * R
    efinal = 0.1 * R
    pplus = pminus + efinal + R
    binding = selected_family_service_no_escape_binding(
        out,
        positive_selected_service=pplus,
        negative_selected_service=pminus,
        final_selected_energy=efinal,
    )
    assert binding.branch == "relink_symmetric_difference"
    assert binding.branch_value == pytest.approx(R)
    assert not binding.identified_with_physical_work
    assert not binding.identified_with_carrier_stock
    assert not binding.identified_with_smooth_k_phys


def test_membership_and_family_sidecars_do_not_clone_or_reweight_inherited_stock():
    switch = selected_family_switch_sidecar([0.2, 0.6, 1.1, 0.4], {0, 2}, {1, 3})
    R = switch["symmetric_difference_energy"]
    stock = _stock(_registration(membership=True, family=True, switch_energy=R))
    out = material_sidecar_stock_decomposition(stock, selected_family_switch_certificate=switch)
    assert out.stock_charge_count == 1
    assert out.inherited_stock_mass == pytest.approx(stock.initial_energy)
    assert {c.event for c in out.charges} == {MATERIAL_MEMBERSHIP_EVENT, SELECTED_FAMILY_EVENT}
    with pytest.raises(ValueError):
        replace(out, stock_cloned_per_sidecar=True)
    with pytest.raises(ValueError):
        replace(out, sidecar_charge_added_to_physical_work=True)
    with pytest.raises(ValueError):
        replace(out, k_phys_identification_used=True)


def test_sidecar_charge_type_rejects_work_stock_and_kphys_aliases():
    with pytest.raises(ValueError):
        MaterialSidecarCharge(
            MATERIAL_MEMBERSHIP_EVENT,
            MEMBERSHIP_PROVENANCE_CURRENCY,
            0.0,
            "zero_charge_provenance",
            physical_work=True,
        )
    with pytest.raises(ValueError):
        MaterialSidecarCharge(
            SELECTED_FAMILY_EVENT,
            SELECTED_FAMILY_MOYAL_CURRENCY,
            1.0,
            "zero_generation_depth_selected_service_boundary",
            smooth_k_phys_relink=True,
        )
