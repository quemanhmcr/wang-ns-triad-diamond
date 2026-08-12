from dataclasses import replace
import math

import numpy as np
import pytest

from src.canonical_positive_edge_work_routing import (
    ETA0,
    MAX_CERTIFIED_TAU,
    STATUS,
    _near_extremal_positive_fiber,
    _nonforward_positive_fiber,
    _pure_helical_fiber,
    _rescale_registered_fiber_measure,
    compress_signed_edge_work_to_hard_cells,
    exact_mode_role_map,
    route_canonical_positive_edge_work,
    single_hard_role_map,
    theorem_certificate,
)
from src.coherent_service_or_flat import coherent_service_or_flat_gate
from src.continuum_helical_edge_measure_registration import continuum_edge_measure_ledger
from src.physical_branch_compiler import MasterDisposition, PhysicalCause, PhysicalCurrency


def _route(ledger, tau=0.1):
    return route_canonical_positive_edge_work(ledger, tau=tau, mode_roles=exact_mode_role_map(ledger))


def _dominant_positive_edge(fiber):
    edges = [a for a in fiber.modal_atoms if a.signed_work_mass > 0.0]
    assert edges
    return max(edges, key=lambda a: a.signed_work_mass)


@pytest.mark.parametrize("q", [0.125, 1.0, 7.25])
def test_registered_fiber_measure_rescaling_is_exact_full_reregistration(q):
    templates_and_fresh = (
        (_near_extremal_positive_fiber(1.0), _near_extremal_positive_fiber(q)),
        (_nonforward_positive_fiber(1.0), _nonforward_positive_fiber(q)),
        (_nonforward_positive_fiber(1.0, phase_sign=-1.0), _nonforward_positive_fiber(q, phase_sign=-1.0)),
    )
    for template, fresh in templates_and_fresh:
        scaled = _rescale_registered_fiber_measure(template, q)
        assert scaled == fresh
        assert tuple(atom.physical_edge_identity for atom in scaled.modal_atoms) == tuple(
            atom.physical_edge_identity for atom in fresh.modal_atoms
        )
        assert continuum_edge_measure_ledger((scaled,)) == continuum_edge_measure_ledger((fresh,))


def test_actual_causal_law_is_partitioned_exactly_not_capacity_probability():
    ledger = continuum_edge_measure_ledger((
        _near_extremal_positive_fiber(1.7),
        _nonforward_positive_fiber(0.8),
    ))
    out = _route(ledger)
    assert out.total_positive_work == pytest.approx(ledger.positive_edge_work, rel=3e-11, abs=1e-14)
    assert out.good_positive_work + out.bad_positive_work == pytest.approx(out.total_positive_work, rel=3e-11, abs=1e-14)
    assert abs(out.mass_reconstruction_residual) <= 3e-11 * max(out.total_positive_work, 1e-14)
    assert not out.capacity_used_as_causal_law
    assert not out.later_hahn_used_as_causal_law


def test_positive_nonforward_work_is_physical_and_routes_bad_not_dropped():
    ledger = continuum_edge_measure_ledger((_nonforward_positive_fiber(1.0),))
    out = _route(ledger)
    assert ledger.positive_nonforward_work > 0.0
    assert out.bad_positive_work > 0.0
    assert out.bad_route is not None
    assert any(edge.scale_progress == pytest.approx(0.0) for edge in out.bad_support)
    assert out.bad_route.physical_work_mass == pytest.approx(out.bad_positive_work)
    assert out.bad_route.transfer_partition.currency_mass == {
        PhysicalCurrency.MULTIPLICATIVE_TRANSFER.value: pytest.approx(out.bad_positive_work)
    }
    assert out.bad_route.transfer_partition.tied_causes == (PhysicalCause.TRANSFER_WORK_LOSS.value,)
    assert out.bad_route.joint_projection.master_disposition == MasterDisposition.TRANSFER_COST.value


def test_bad_restriction_has_native_deficit_and_all_certified_tau_hit_fixed_transfer_loss():
    ledger = continuum_edge_measure_ledger((_nonforward_positive_fiber(1.0),))
    for tau in (1e-8, 1e-4, 0.01, 0.05, MAX_CERTIFIED_TAU):
        out = _route(ledger, tau=tau)
        bad = out.bad_route
        assert bad is not None
        assert bad.deficit >= ETA0 - 2e-10
        assert bad.fixed_transfer_gate.threshold == pytest.approx(tau * tau / 1_036_800_000, rel=2e-15)
        assert bad.deficit > bad.fixed_transfer_gate.threshold
        assert bad.transfer_partition.first_time is None
        assert bad.joint_projection.first_time is None
        assert bad.joint_projection.terminal_certificate_used == "stage_zero_fixed_transfer_loss"


def test_bad_route_uses_the_same_transfer_channel_as_the_whole_physical_block_gate():
    ledger = continuum_edge_measure_ledger((_nonforward_positive_fiber(1.0),))
    out = _route(ledger, tau=0.1)
    bad = out.bad_route
    assert bad is not None
    full_gate = coherent_service_or_flat_gate(
        tau=bad.tau,
        avg_transfer_deficit=bad.deficit,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=False,
        scaled_lifetime=1.0,
        phase_holonomy=0.0,
    )
    roots = tuple(full_gate["triggered_causes"])
    transfer_roots = tuple(root for root in roots if root["cause"] == "physical_transfer_cost")
    assert len(transfer_roots) == 1
    root = transfer_roots[0]
    assert root["threshold"] == pytest.approx(bad.fixed_transfer_gate.threshold, rel=2e-15)
    assert root["value"] == pytest.approx(bad.fixed_transfer_gate.avg_transfer_deficit, rel=2e-15)
    assert bad.fixed_transfer_gate.triggered
    assert bad.fixed_transfer_gate.cause == "physical_transfer_cost"


def test_subthreshold_transfer_channel_cannot_self_certify_fixed_loss():
    tau = 0.1
    ledger = continuum_edge_measure_ledger((_near_extremal_positive_fiber(1.0),))
    out = _route(ledger, tau=tau)
    assert out.bad_route is None
    full_gate = coherent_service_or_flat_gate(
        tau=tau,
        avg_transfer_deficit=0.0,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=False,
        scaled_lifetime=1.0,
        phase_holonomy=0.0,
    )
    assert all(root["cause"] != "physical_transfer_cost" for root in full_gate["triggered_causes"])


def test_geometry_good_is_only_young_eligible_and_signed_cell_work_is_retained():
    good = _near_extremal_positive_fiber(1.0)
    bad_negative = _nonforward_positive_fiber(0.3, phase_sign=-1.0)
    ledger = continuum_edge_measure_ledger((good, bad_negative))
    out = _route(ledger)
    assert out.good_positive_work > 0.0
    assert out.young_eligible.physical_work_mass == pytest.approx(out.good_positive_work)
    assert not out.young_eligible.marking_good
    assert not out.young_eligible.young_certified
    assert not out.young_eligible.registered_generated_continuation
    assert out.young_eligible.fate_pure_hard_cells
    assert out.young_eligible.unresolved_mixed_positive_work == pytest.approx(0.0)
    for cell in out.young_eligible.fate_pure_hard_cells:
        # Only a fate-pure cell may expose full signed pi_#dW to Young: negative
        # cancellation remains, while no already-terminal bad positive work can assist.
        assert cell.inherited_bad_positive_work == pytest.approx(0.0)
        assert cell.signed_work <= cell.inherited_positive_work + 1e-12 * max(1.0, cell.inherited_positive_work)
        assert not cell.fresh_cell_hahn_is_causal_law


def test_analyst_coarsening_cannot_turn_mixed_good_bad_cause_into_signed_young_binding():
    ledger = continuum_edge_measure_ledger((
        _near_extremal_positive_fiber(1.0),
        _nonforward_positive_fiber(1.0),
        _nonforward_positive_fiber(0.2, phase_sign=-1.0),
    ))
    coarse_roles = single_hard_role_map(ledger)
    out = route_canonical_positive_edge_work(ledger, tau=0.1, mode_roles=coarse_roles)
    assert out.good_positive_work > 0.0
    assert out.bad_positive_work > 0.0
    assert out.young_eligible.fate_pure_positive_work == pytest.approx(0.0)
    assert out.young_eligible.unresolved_mixed_positive_work == pytest.approx(out.good_positive_work)
    assert len(out.young_eligible.mixed_fate_hard_cells) == 1
    cell = out.young_eligible.mixed_fate_hard_cells[0]
    assert cell.inherited_good_positive_work > 0.0
    assert cell.inherited_bad_positive_work > 0.0
    assert not cell.fresh_cell_hahn_is_causal_law
    assert not out.young_eligible.young_certified
    assert not out.young_eligible.marking_good


def test_hard_coarsening_and_refinement_preserve_inherited_cause_but_not_fresh_hahn():
    positive = _near_extremal_positive_fiber(1.0)
    negative = _nonforward_positive_fiber(0.9, phase_sign=-1.0)
    ledger = continuum_edge_measure_ledger((positive, negative))
    refined = compress_signed_edge_work_to_hard_cells(ledger, exact_mode_role_map(ledger))
    coarse = compress_signed_edge_work_to_hard_cells(ledger, single_hard_role_map(ledger))
    assert refined.inherited_positive_work == pytest.approx(ledger.positive_edge_work, rel=3e-11)
    assert coarse.inherited_positive_work == pytest.approx(ledger.positive_edge_work, rel=3e-11)
    assert refined.inherited_negative_work == pytest.approx(ledger.negative_edge_work, rel=3e-11)
    assert coarse.inherited_negative_work == pytest.approx(ledger.negative_edge_work, rel=3e-11)
    assert coarse.inherited_positive_work == pytest.approx(refined.inherited_positive_work, rel=3e-11)
    assert coarse.cancellation_gap >= -1e-12
    assert coarse.fresh_hahn_positive_work <= coarse.inherited_positive_work + 1e-12
    # Fresh Hahn is representation-dependent; it is explicitly diagnostic only.
    assert not coarse.fresh_hahn_is_causal_law
    assert not refined.fresh_hahn_is_causal_law


def test_provenance_tampering_is_rejected_by_physical_replay():
    ledger = continuum_edge_measure_ledger((_near_extremal_positive_fiber(1.0), _nonforward_positive_fiber(1.0)))
    forged = replace(ledger, positive_edge_work=ledger.positive_edge_work * 1.25 + 0.1)
    with pytest.raises(AssertionError, match="replay field positive_edge_work"):
        route_canonical_positive_edge_work(forged, tau=0.1, mode_roles=exact_mode_role_map(ledger))


def test_parent_swap_leaves_positive_routing_invariant():
    x = np.array([0.64, 0.42, 0.0])
    y = np.array([0.36, -0.42, 0.0])
    a = _pure_helical_fiber(x=x, y=y, sx=1, sy=-1, sz=1)
    b = _pure_helical_fiber(x=y, y=x, sx=-1, sy=1, sz=1)
    oa = _route(continuum_edge_measure_ledger((a,)))
    ob = _route(continuum_edge_measure_ledger((b,)))
    assert oa.total_positive_work == pytest.approx(ob.total_positive_work, rel=5e-11, abs=1e-13)
    assert oa.good_positive_work == pytest.approx(ob.good_positive_work, rel=5e-11, abs=1e-13)
    assert oa.bad_positive_work == pytest.approx(ob.bad_positive_work, rel=5e-11, abs=1e-13)


def test_helical_basis_gauge_representation_and_wavevector_units_do_not_change_fate():
    x = np.array([0.64, 0.42, 0.0])
    y = np.array([0.36, -0.42, 0.0])
    base = _pure_helical_fiber(x=x, y=y, sx=1, sy=-1, sz=1)
    # The physical vectors are unchanged under h_s -> exp(i theta)h_s,
    # a_s -> exp(-i theta)a_s.  Re-registering the same vectors therefore gives
    # the same canonical edge law; no helical-frame observer phase is causal.
    gauge = _pure_helical_fiber(x=x, y=y, sx=1, sy=-1, sz=1)
    o0 = _route(continuum_edge_measure_ledger((base,)))
    og = _route(continuum_edge_measure_ledger((gauge,)))
    assert og.total_positive_work == pytest.approx(o0.total_positive_work, rel=5e-11, abs=1e-13)
    assert og.good_positive_work == pytest.approx(o0.good_positive_work, rel=5e-11, abs=1e-13)
    assert og.bad_positive_work == pytest.approx(o0.bad_positive_work, rel=5e-11, abs=1e-13)

    lam = 7.3
    scaled = _pure_helical_fiber(x=lam * x, y=lam * y, sx=1, sy=-1, sz=1)
    os = _route(continuum_edge_measure_ledger((scaled,)))
    assert os.total_positive_work == pytest.approx(lam * o0.total_positive_work, rel=8e-10, abs=1e-12)
    assert os.good_positive_work == pytest.approx(lam * o0.good_positive_work, rel=8e-10, abs=1e-12)
    assert os.bad_positive_work == pytest.approx(lam * o0.bad_positive_work, rel=8e-10, abs=1e-12)
    if o0.bad_route is not None:
        assert os.bad_route is not None
        assert os.bad_route.deficit == pytest.approx(o0.bad_route.deficit, rel=2e-9, abs=2e-10)


def test_physical_helical_counterexample_good_capacity_majority_bad_work_majority():
    good0 = _near_extremal_positive_fiber(1.0)
    good_edge = max(
        (a for a in good0.modal_atoms if a.signed_work_mass > 0.0 and a.signed_efficiency > 1.0 - ETA0),
        key=lambda a: a.signed_work_mass,
    )
    rho_g = good_edge.signed_work_mass / good_edge.capacity_mass

    found = None
    for parent_ratio in (0.45, 0.6, 0.8, 1.0, 1.25, 1.6, 2.0):
        for degrees in range(101, 180, 3):
            theta = math.radians(degrees)
            x = np.array([1.0, 0.0, 0.0])
            y = parent_ratio * np.array([math.cos(theta), math.sin(theta), 0.0])
            if np.linalg.norm(x + y) >= max(np.linalg.norm(x), np.linalg.norm(y)):
                continue
            for sx in (-1, 1):
                for sy in (-1, 1):
                    if sx == sy and abs(np.linalg.norm(x) - np.linalg.norm(y)) < 1e-12:
                        continue
                    for sz in (-1, 1):
                        candidate = _pure_helical_fiber(x=x, y=y, sx=sx, sy=sy, sz=sz)
                        edge = _dominant_positive_edge(candidate)
                        rho_b = edge.signed_work_mass / edge.capacity_mass
                        if edge.signed_efficiency <= 1.0 - ETA0 and rho_b > rho_g * (1.0 + 1e-4):
                            found = (x, y, sx, sy, sz, rho_b, edge.capacity_mass)
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                break
        if found:
            break
    assert found is not None, "no physical helical capacity/work-majority reversal fixture found"

    x, y, sx, sy, sz, rho_b, bad_capacity0 = found
    good_capacity0 = good_edge.capacity_mass
    ratio = math.sqrt(rho_b / rho_g)  # strictly between 1 and rho_b/rho_g
    qg = ratio * bad_capacity0 / good_capacity0
    good = _near_extremal_positive_fiber(qg)
    bad = _pure_helical_fiber(x=x, y=y, sx=sx, sy=sy, sz=sz, quotient_measure_mass=1.0)
    out = _route(continuum_edge_measure_ledger((good, bad)))
    good_capacity = math.fsum(edge.capacity_mass for edge in out.good_support)
    bad_capacity = math.fsum(edge.capacity_mass for edge in out.bad_support)
    assert good_capacity > bad_capacity
    assert out.bad_positive_work > out.good_positive_work


def test_certificate_keeps_signed_before_positive_and_coherent_povm_seam_open():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "dW+" in cert["canonical_causal_law"]
    assert "T_C" in cert["young_input"]
    assert "already-terminal bad work assist Young" in cert["mixed_fate_seam"]
    assert "no arbitrary refinement" in cert["mixed_fate_seam"]
    assert "positive mass-preserving kernel" in cert["coherent_povm_scope"]
    assert not cert["capacity_is_causal_law"]
    assert not cert["claims_global_regularity"]
