import math

import numpy as np
import pytest

from src.high_tail_binary_work_reentry import (
    STATUS,
    binary_hh_common_work_law,
    common_unit_regeneration_owners,
    high_tail_hh_binary_reentry,
    own_scale_reweighting_counterexample,
    theorem_certificate,
)


def _atoms(net_positive: float) -> np.ndarray:
    # Signed work = net_positive, while atomic positive mass is larger.
    return np.array([[[net_positive + 3.0, -1.0], [-1.0, -1.0]]], dtype=float)


def test_common_unit_regeneration_owner_has_clean_half_tail_lower():
    D = 2.0
    nu = 0.5
    W = 1.2
    shell = [0.7, 0.6]
    out = common_unit_regeneration_owners(D, nu, W, shell, 0.8, 0.6)
    assert out["scaled_shell_positive_work_sum"] == pytest.approx(1.3)
    assert out["owner_threshold"] == pytest.approx(0.65)
    assert out["clean_owner_threshold"] == pytest.approx(0.5)
    assert out["joint_owners"] == ("positive_HH_regeneration",)
    assert out["causal_probability_reweighted_by_shell_scale"] is False


def test_common_unit_exact_tie_keeps_hh_and_interface_jointly():
    out = common_unit_regeneration_owners(
        1.0, 1.0, 1.0, [0.5, 0.5], 0.5, 0.5
    )
    assert set(out["joint_owners"]) == {
        "positive_HH_regeneration",
        "positive_resolved_cross_interface",
    }


def test_event_resolved_hahn_atoms_dominate_positive_HH_and_define_binary_law():
    atoms = [_atoms(2.0), _atoms(1.0)]
    out = binary_hh_common_work_law([1, 3], atoms)
    assert out["aggregate_positive_HH_common_work"] == pytest.approx(3.0)
    assert out["binary_positive_common_work"] == pytest.approx(9.0)
    assert out["atomic_positive_dominance_margin"] == pytest.approx(6.0)
    assert sum(e["probability"] for e in out["events"]) == pytest.approx(1.0)
    assert out["causal_probability_uses_common_N_work"] is True
    assert out["causal_probability_uses_Mj_reweighting"] is False


def test_common_work_unit_rescaling_does_not_change_binary_probabilities():
    atoms = [_atoms(2.0), _atoms(1.0)]
    a = binary_hh_common_work_law([1, 4], atoms)
    b = binary_hh_common_work_law([1, 4], [17.0 * x for x in atoms])
    pa = [e["probability"] for e in a["events"]]
    pb = [e["probability"] for e in b["events"]]
    assert pa == pytest.approx(pb)


def test_HH_primary_owner_supplies_clean_binary_work_but_not_productivity():
    D = 1.0
    nu = 1.0
    W = 1.0
    shell = [0.6, 0.4]
    H = 0.75
    I = 0.25
    atoms = [_atoms(0.45), _atoms(0.30)]
    out = high_tail_hh_binary_reentry(D, nu, W, shell, H, I, [1, 2], atoms)
    assert out["HH_primary"] is True
    assert out["clean_binary_positive_common_work_if_HH_owner"] == pytest.approx(0.5)
    binary = out["binary_HH_law"]
    assert binary["binary_positive_common_work"] >= 0.5
    assert binary["H_inf_weighted_selected_own_scale_binary_work"] >= 1.0
    assert out["productivity_energy_gate_supplied"] is False
    assert out["Young_near_extremality_supplied"] is False
    assert out["parent_child_scale_locality_supplied"] is False
    assert out["master_semantics"] == "JOINT_PHYSICAL_WORK_OWNERS__INTERFACE_QUOTIENT_BEFORE_RECURSION"
    assert "RESOLVED_INTERFACE_DONOR_QUOTIENT" in out["next_owner_if_interface"]


def test_scale_reweighting_counterexample_changes_cause_probabilities():
    out = own_scale_reweighting_counterexample()
    assert out["common_unit_probabilities"] == pytest.approx((0.5, 0.5))
    assert out["maximum_probability_distortion"] > 0.49
    assert out["own_scale_reweighted_probabilities"][1] > 0.99


def test_certificate_forbids_scale_reweighting_and_keeps_locality_open():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "N times actual positive child-energy work" in cert["common_unit"]
    assert "must not redefine probabilities" in cert["anti_reweight"]
    assert "KL/log-productivity remains conditional" in cert["productivity_scope"]
    assert "nonlocal K>>M geometry remains" in cert["locality_scope"]
    assert "same-event conservative donor tracing" in cert["interface_continuation"]
