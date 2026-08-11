import math

import pytest

from src.continuum_master_event_quotient import (
    DiagnosticConcentration,
    EventDisposition,
    NativePathLedger,
    RecursiveEventState,
    STATUS,
    SupplierKind,
    bounded_scale_full_survivor_steps_to_boundary,
    canonical_owner_bundle,
    compose_native_ledgers,
    forbid_diagnostic_as_causal_action,
    full_natural_survivor_endpoint,
    geometric_uv_natural_time_sum,
    log_scale_telescope,
    master_escape_dichotomy,
    physical_time_telescope,
    supplier_scale_certificate,
    theorem_certificate,
    trace_full_natural_survivors,
    validate_supplier_scale,
    zero_charge_owner_relay,
    zero_charge_witness_relay,
)


def test_zero_charge_relays_cannot_duplicate_one_physical_law():
    base = canonical_owner_bundle("actual positive SGS source law", 3.5, ("sgs", "sgs"))
    out = zero_charge_owner_relay(base, ("coherent_service", "critical_shell", "coherent_service"))
    assert out.mass == 3.5
    assert out.owners == ("coherent_service", "critical_shell", "sgs")
    assert sum([out.mass]) == 3.5


def test_witness_relay_can_change_observable_units_without_creating_a_second_charge():
    relay = zero_charge_witness_relay(
        "actual resolved pressure-pair source law",
        "critical u-shell mass",
        "pair capacity plus strict-lowpass contraction",
    )
    assert not relay.causal_charge_created
    assert not relay.diagnostic_probability_created
    assert relay.upstream_physical_measure != relay.downstream_state_observable


def test_event_state_has_physical_time_and_t0_is_absorbing():
    state = RecursiveEventState(0.7, 4.0, "actual positive HH work", ("HH_regeneration",))
    assert not state.absorbing
    root = RecursiveEventState(
        0.0,
        4.0,
        "initial data",
        ("t=0",),
        EventDisposition.ABSORBING_INITIAL_BOUNDARY,
    )
    assert root.absorbing
    with pytest.raises(ValueError):
        RecursiveEventState(0.1, 4.0, "initial data", ("t=0",), EventDisposition.ABSORBING_INITIAL_BOUNDARY)


def test_supplier_scale_progress_is_exactly_branch_specific():
    validate_supplier_scale(SupplierKind.GENERATED_SIGNED_GOOD_HH, 10.0, 6.1)
    validate_supplier_scale(SupplierKind.RESOLVED_DISSIPATION, 10.0, 2.5)
    validate_supplier_scale(SupplierKind.PRESSURE_PAIR, 10.0, 1.0)
    validate_supplier_scale(SupplierKind.FRESH_SGS_SCALE, 10.0, 20.0)
    validate_supplier_scale(SupplierKind.FRESH_SGS_SCALE, 10.0, 0.1)
    validate_supplier_scale(SupplierKind.HIGH_TAIL, 10.0, 20.0)
    with pytest.raises(ValueError):
        validate_supplier_scale(SupplierKind.GENERATED_SIGNED_GOOD_HH, 10.0, 6.5)
    with pytest.raises(ValueError):
        validate_supplier_scale(SupplierKind.RESOLVED_DISSIPATION, 10.0, 3.0)
    with pytest.raises(ValueError):
        validate_supplier_scale(SupplierKind.HIGH_TAIL, 10.0, 19.0)
    assert supplier_scale_certificate(SupplierKind.FRESH_SGS_SCALE).directional_progress == "no_directional_progress_supplied"
    assert "no_supplier_relative" in supplier_scale_certificate(SupplierKind.GENERIC_CRITICAL_SHELL).directional_progress


def test_physical_time_and_log_scale_telescope_without_common_clock():
    tout = physical_time_telescope([1.0, 0.8, 0.31, 0.2])
    assert tout["sum_physical_time_drops"] == pytest.approx(0.8)
    assert tout["residual"] == pytest.approx(0.0, abs=1e-15)
    sout = log_scale_telescope([4.0, 2.0, 8.0, 3.0])
    assert sout["endpoint_log_scale_change"] == pytest.approx(math.log(3.0 / 4.0))
    assert sout["residual"] == pytest.approx(0.0, abs=1e-15)


def test_full_natural_survivor_uses_its_own_physical_duration_and_t0_truncates():
    out = full_natural_survivor_endpoint(1.0, 2.0, 1.0)
    assert out["requested_duration"] == pytest.approx(0.25)
    assert out["end_time"] == pytest.approx(0.75)
    assert not out["hits_initial_boundary"]
    assert out["disposition"] == "full_natural_analysis_checkpoint"
    assert out["physical_event_created"] is False
    assert out["recursion_edges_added"] == 0
    root = full_natural_survivor_endpoint(0.2, 2.0, 1.0)
    assert root["end_time"] == 0.0
    assert root["hits_initial_boundary"]
    assert root["disposition"] == EventDisposition.ABSORBING_INITIAL_BOUNDARY.value


def test_bounded_scale_full_survivor_tail_must_hit_initial_boundary():
    t0, mbar, c = 1.1, 3.0, 0.7
    k = bounded_scale_full_survivor_steps_to_boundary(t0, mbar, c)
    assert k == math.ceil(t0 * mbar * mbar / c)
    # Using the largest allowed scale gives the shortest admissible natural windows,
    # hence it is the hardest bounded-scale case for reaching t=0.
    out = trace_full_natural_survivors(t0, [mbar] * k, c)
    assert out["hits_initial_boundary"]
    assert out["final_time"] == 0.0
    assert out["physical_time_telescope_residual"] == pytest.approx(0.0, abs=1e-14)
    assert out["recursive_events_added"] == 0


def test_uv_geometric_survivors_have_finite_total_physical_natural_time():
    total = geometric_uv_natural_time_sum(2.0, 1.0, 2.0)
    assert total == pytest.approx((1.0 / 4.0) / (1.0 - 1.0 / 4.0))
    assert total == pytest.approx(1.0 / 3.0)
    assert math.isfinite(total)


def test_typed_ledger_composes_without_scalar_exchange_rates():
    a = NativePathLedger(0.2, math.log(2.0), 0.3, 0.0, 0.01, (("energy", 2.0),))
    b = NativePathLedger(0.1, -math.log(4.0), 0.0, 0.5, 0.02, (("energy", 1.0), ("enstrophy", 4.0)))
    out = compose_native_ledgers((a, b))
    assert out.physical_time_drop == pytest.approx(0.3)
    assert out.log_scale_change == pytest.approx(-math.log(2.0))
    assert out.multiplicative_transfer_cost == pytest.approx(0.3)
    assert out.causal_reuse_action == pytest.approx(0.5)
    assert out.global_resources == (("energy", 3.0), ("enstrophy", 4.0))


def test_concentration_coordinates_cannot_be_promoted_to_causal_entropy():
    h = DiagnosticConcentration("H_inf_time", math.log(4.0), "positive comparable HH work")
    with pytest.raises(TypeError, match="cannot be charged as causal entropy"):
        forbid_diagnostic_as_causal_action(h)


def test_certificate_states_exact_remaining_escape_dichotomy_without_regularity_claim():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "same-carrier checkpoint segmentation" in cert["infinite_escape_dichotomy"]
    assert "infinite recursive EVENT path" in cert["infinite_escape_dichotomy"]
    assert "not a canonical physical lineage" in cert["uv_obstruction"]
    assert "cannot replace the event-anchored smooth carrier" in cert["same_carrier_checkpoint_segmentation"]
    assert "not a global no-escape" in cert["scope"]
    assert "not causal Shannon/Renyi action" in cert["diagnostic_separation"]
    dich = master_escape_dichotomy()
    assert "cannot replace the event-anchored smooth carrier" in dich["proof"]
    assert "continuation of the same carrier" in dich["proof"]


def test_unrouted_coefficient_obstruction_cannot_enter_canonical_master_state():
    from src.common_slice_coefficient_registration import ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION

    with pytest.raises(TypeError, match="unrouted coefficient obstruction"):
        RecursiveEventState(
            0.6,
            4.0,
            "coefficient first-stop locator",
            (ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,),
        )
    with pytest.raises(TypeError, match="unrouted coefficient obstruction"):
        canonical_owner_bundle(
            "coefficient first-stop locator",
            1.0,
            (ROLE_INTERFACE_COEFFICIENT_OBSTRUCTION,),
        )


def test_actual_energy_reentry_adapter_creates_only_physical_owner_bundle():
    from src.continuum_master_event_quotient import owner_bundle_from_energy_reentry

    reentry = {
        "branch": "physical_high_high_transfer_generation",
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
    }
    out = owner_bundle_from_energy_reentry("actual positive q^2-weighted HH work", 2.5, reentry)
    assert out.mass == pytest.approx(2.5)
    assert out.owners == ("physical_high_high_transfer_generation",)

    with pytest.raises(TypeError, match="coefficient impulse"):
        owner_bundle_from_energy_reentry(
            "bad locator",
            1.0,
            {
                "branch": "physical_high_high_transfer_generation",
                "coefficient_impulse_used_as_physical_work": True,
                "observer_partition_motion_charged_as_physics": False,
            },
        )


def test_full_natural_service_witness_cannot_be_promoted_to_separate_owner_event():
    from src.full_natural_service_corridor_quotient import FULL_NATURAL_SERVICE_WITNESS

    with pytest.raises(TypeError, match="same-corridor physical witness"):
        RecursiveEventState(
            0.5,
            4.0,
            "completed full-natural service law",
            (FULL_NATURAL_SERVICE_WITNESS,),
        )
    with pytest.raises(TypeError, match="same-corridor physical witness"):
        canonical_owner_bundle(
            "completed full-natural service law",
            1.0,
            (FULL_NATURAL_SERVICE_WITNESS,),
        )


def test_full_natural_checkpoint_cannot_enter_recursive_event_state_or_owner_bundle():
    from src.full_natural_checkpoint_quotient import FULL_NATURAL_CHECKPOINT

    with pytest.raises(TypeError, match="analysis re-registration"):
        canonical_owner_bundle(
            "no-hit horizon checkpoint",
            1.0,
            (FULL_NATURAL_CHECKPOINT,),
        )
    with pytest.raises(TypeError, match="analysis checkpoint"):
        RecursiveEventState(
            0.5,
            4.0,
            "no-hit horizon",
            (),
            EventDisposition.FULL_NATURAL_SURVIVOR,
        )


def _smooth_relink_reentry(*, with_strain: bool = False):
    from src.smooth_quadratic_carrier_interface import (
        RELINK_OWNER,
        STRAIN_OWNER,
        GaugeQuotientedInterfaceWork,
    )

    T = (
        (0.0, 3.0, -1.0),
        (-3.0, 0.0, 2.0),
        (1.0, -2.0, 0.0),
    )
    relink = (2.0, -1.0, -1.0)
    strain = (2.0, 0.0, 0.0) if with_strain else (0.0, 0.0, 0.0)
    native = tuple(a + b for a, b in zip(relink, strain, strict=True))
    work = GaugeQuotientedInterfaceWork(
        native,
        relink,
        strain,
        0.0,
        0.0,
        T,
    )
    owners = (RELINK_OWNER, STRAIN_OWNER) if with_strain else (RELINK_OWNER,)
    return {
        "branch": "smooth_interface_physical_work",
        "joint_interface_owners": owners,
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
        "gauge_quotiented_interface_work_certificate": work,
    }


def test_raw_smooth_relink_cannot_enter_recursive_owner_bundle():
    from src.smooth_quadratic_carrier_interface import RELINK_OWNER

    with pytest.raises(TypeError, match="same-event donor relay"):
        canonical_owner_bundle("smooth K_phys relink", 1.0, (RELINK_OWNER,))


def test_pure_smooth_relink_energy_reentry_becomes_zero_depth_same_event_relay():
    from src.continuum_master_event_quotient import energy_reentry_master_route
    from src.smooth_relink_donor_quotient import SMOOTH_RELINK_SAME_EVENT_RELAY

    route = energy_reentry_master_route("positive native smooth interface work", 2.0, _smooth_relink_reentry())
    assert route.owner_bundle is None
    assert route.recursive_event_created is False
    assert route.same_event_relays == (SMOOTH_RELINK_SAME_EVENT_RELAY,)
    assert route.smooth_relink_donor_certificate is not None
    assert route.smooth_relink_donor_certificate.recursive_generation_created is False


def test_smooth_relink_plus_strain_keeps_only_strain_as_recursive_owner():
    from src.continuum_master_event_quotient import energy_reentry_master_route
    from src.smooth_quadratic_carrier_interface import STRAIN_OWNER
    from src.smooth_relink_donor_quotient import SMOOTH_RELINK_SAME_EVENT_RELAY

    route = energy_reentry_master_route(
        "positive native smooth interface work",
        4.0,
        _smooth_relink_reentry(with_strain=True),
    )
    assert route.owner_bundle is not None
    assert route.owner_bundle.owners == (STRAIN_OWNER,)
    assert route.recursive_event_created is True
    assert route.same_event_relays == (SMOOTH_RELINK_SAME_EVENT_RELAY,)


def test_legacy_owner_bundle_wrapper_rejects_pure_relink_zero_depth_route():
    from src.continuum_master_event_quotient import owner_bundle_from_energy_reentry

    with pytest.raises(TypeError, match="same-event donor relay"):
        owner_bundle_from_energy_reentry(
            "positive native smooth interface work",
            2.0,
            _smooth_relink_reentry(),
        )


def test_master_certificate_closes_eventually_pure_high_strain_but_not_mixed_recurrence():
    cert = theorem_certificate()
    hs = cert["high_strain_descending_epoch"]
    assert "genuine recursive owner" in hs
    assert "N_next/N<=3/16" in hs
    assert "arbitrary interval overlap" in hs
    assert "without promoting D_V to a global reset" in hs
    dich = master_escape_dichotomy()
    assert "eventually consisting only of high-strain" in dich["statement"]
    assert "mixed recurrence" in dich["remaining_physics"]
    assert "eventually-pure high-strain recurrence" in cert["scope"]
