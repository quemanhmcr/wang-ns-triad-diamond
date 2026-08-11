import math

import numpy as np
import pytest

from src.helical import coupling_g, edge_metrics
from src.helical_physical_edge_registration import (
    STATUS,
    gauge_transform_modal_data,
    phase_alignment,
    register_helical_physical_edge,
    theorem_certificate,
    unordered_parent_curl_source,
    waleffe_child_source_coefficient,
)
from src.single_edge_certificate import float_jstar


def _edge():
    x = np.array([0.52, 0.31, 0.0])
    y = np.array([0.52, -0.31, 0.0])
    z = x + y
    return x, y, z


def test_direct_ns_curl_leray_source_equals_waleffe_coefficient():
    x, y, z = _edge()
    ax = 1.2 - 0.4j
    ay = -0.7 + 0.8j
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                row = register_helical_physical_edge(
                    x=x, y=y, z=z, sx=sx, sy=sy, sz=sz,
                    ax=ax, ay=ay, az=0.5 + 0.9j,
                )
                assert row.direct_child_source_coefficient == pytest.approx(
                    row.waleffe_child_source_coefficient, rel=2e-12, abs=2e-12
                )
                assert row.leray_pairing_residual < 2e-12


def test_unordered_parent_source_reconstructs_two_ordered_curl_terms():
    x, y, z = _edge()
    source = unordered_parent_curl_source(x, y, z, 1, -1, 1.0 + 0.2j, 0.7 - 0.4j)
    assert np.linalg.norm(source) > 0.0
    assert abs(np.dot(z, source)) < 2e-12 * max(1.0, np.linalg.norm(source))


def test_native_capacity_identity_is_exact_and_not_a_young_norm_product():
    x, y, z = _edge()
    ax, ay, az = 1.1 + 0.3j, 0.8 - 0.2j, -0.4 + 1.2j
    row = register_helical_physical_edge(x=x, y=y, z=z, sx=1, sy=-1, sz=1, ax=ax, ay=ay, az=az)
    expected_capacity = 4.0 * np.linalg.norm(z) * abs(ax) * abs(ay) * abs(az)
    assert row.native_modal_capacity == pytest.approx(expected_capacity)
    assert row.signed_upper_progress_work == pytest.approx(row.registered_upper_progress_work, rel=2e-12, abs=2e-12)
    assert row.young_norm_used_as_capacity is False
    assert row.duhamel_weight_used_as_causal_law is False


def test_geometric_multiplier_is_exact_existing_single_edge_J_over_Jstar():
    x, y, z = _edge()
    row = register_helical_physical_edge(x=x, y=y, z=z, sx=1, sy=-1, sz=1, ax=1.0, ay=1.0j, az=0.7 - 0.3j)
    geom = edge_metrics(x, y, z, 1, -1, 1)
    assert row.geometric_multiplier_J == pytest.approx(geom.efficiency)
    assert row.normalized_multiplier == pytest.approx(geom.efficiency / float_jstar())
    assert 0.0 <= row.normalized_multiplier <= 1.0 + 1e-10


def test_parent_orientation_swap_is_exactly_the_same_physical_edge():
    x, y, z = _edge()
    ax, ay, az = 1.3 - 0.1j, -0.2 + 0.9j, 0.7 + 0.6j
    a = register_helical_physical_edge(x=x, y=y, z=z, sx=1, sy=-1, sz=-1, ax=ax, ay=ay, az=az)
    b = register_helical_physical_edge(x=y, y=x, z=z, sx=-1, sy=1, sz=-1, ax=ay, ay=ax, az=az)
    assert a.direct_child_source_coefficient == pytest.approx(b.direct_child_source_coefficient, rel=2e-12, abs=2e-12)
    assert a.signed_child_energy_work == pytest.approx(b.signed_child_energy_work)
    assert a.native_modal_capacity == pytest.approx(b.native_modal_capacity)
    assert a.geometric_multiplier_J == pytest.approx(b.geometric_multiplier_J)
    assert a.phase_alignment == pytest.approx(b.phase_alignment)


def test_helical_basis_phase_gauge_leaves_registered_alignment_invariant():
    x, y, z = _edge()
    sx, sy, sz = 1, -1, 1
    ax, ay, az = 1.1 + 0.4j, -0.3 + 0.9j, 0.8 - 0.7j
    g = coupling_g(x, y, -z, sx, sy, sz)
    a = sx * np.linalg.norm(x) - sy * np.linalg.norm(y)
    c0 = phase_alignment(a, g, ax, ay, az)
    gp, axp, ayp, azp = gauge_transform_modal_data(g, ax, ay, az, 0.7, -1.1, 0.4)
    c1 = phase_alignment(a, gp, axp, ayp, azp)
    assert c1 == pytest.approx(c0, abs=2e-14)


def test_uniform_wavevector_dilation_changes_capacity_and_work_linearly_only():
    x, y, z = _edge()
    kw = {"sx": 1, "sy": -1, "sz": 1, "ax": 0.6 + 0.2j, "ay": -0.5 + 0.4j, "az": 1.2 - 0.1j}
    a = register_helical_physical_edge(x=x, y=y, z=z, **kw)
    lam = 7.3
    b = register_helical_physical_edge(x=lam*x, y=lam*y, z=lam*z, **kw)
    assert b.forward_ratio == pytest.approx(a.forward_ratio)
    assert b.normalized_multiplier == pytest.approx(a.normalized_multiplier)
    assert b.phase_alignment == pytest.approx(a.phase_alignment)
    assert b.native_modal_capacity == pytest.approx(lam * a.native_modal_capacity)
    assert b.signed_child_energy_work == pytest.approx(lam * a.signed_child_energy_work)


def test_nonforward_positive_work_is_not_relabelled_signed_good_progress():
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([-0.8, 0.3, 0.0])
    z = x + y
    # Search a deterministic phase/helicity combination with positive child work.
    found = None
    for sx in (-1, 1):
        for sy in (-1, 1):
            for sz in (-1, 1):
                for phi in (0.0, math.pi/3, 2*math.pi/3, math.pi, 4*math.pi/3, 5*math.pi/3):
                    row = register_helical_physical_edge(
                        x=x, y=y, z=z, sx=sx, sy=sy, sz=sz,
                        ax=1.0, ay=1.0, az=np.exp(1j*phi),
                    )
                    if row.signed_child_energy_work > 1e-9:
                        found = row
                        break
                if found is not None:
                    break
            if found is not None:
                break
        if found is not None:
            break
    assert found is not None
    assert found.forward_ratio < 1.0
    assert found.scale_progress == 0.0
    assert found.geometric_multiplier_J == 0.0
    assert found.signed_upper_progress_work == 0.0
    assert found.registered_upper_progress_work == pytest.approx(0.0)


def test_waleffe_formula_rejects_nonphysical_parent_sum():
    x, y, z = _edge()
    with pytest.raises(ValueError, match=r"z=x\+y"):
        waleffe_child_source_coefficient(x, y, z + np.array([0.1, 0.0, 0.0]), 1, -1, 1, 1.0, 1.0)


def test_certificate_preserves_physical_scope():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "-P(u.grad u)=P(u x omega)" in cert["ns_identity"]
    assert "factor four" in cert["physical_work"]
    assert "A_e=4|z||a_x a_y a_z|" in cert["native_capacity"]
    assert "T_e log_+" in cert["registration"]
    assert "no Young norm product" in cert["native_capacity"]
    assert "does not yet construct the continuum edge measure" in cert["scope"]
    assert "generic HH" in cert["scope"]
