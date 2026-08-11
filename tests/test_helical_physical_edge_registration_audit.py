"""Adversarial native-scale and typed-certificate checks for physical edges."""

from dataclasses import replace

import numpy as np
import pytest

from src.helical_physical_edge_registration import register_helical_physical_edge


def _edge(scale: float = 1.0):
    x = scale * np.asarray((0.52, 0.31, 0.0))
    y = scale * np.asarray((0.52, -0.31, 0.0))
    return x, y, x + y


def _registration(scale: float = 1.0):
    x, y, z = _edge(scale)
    return register_helical_physical_edge(
        x=x,
        y=y,
        z=z,
        sx=1,
        sy=-1,
        sz=1,
        ax=1.1 + 0.3j,
        ay=0.8 - 0.2j,
        az=-0.4 + 1.2j,
    )


def test_uniform_wavevector_dilation_survives_at_tiny_native_scale():
    base = _registration()
    scale = 1.0e-120
    tiny = _registration(scale)
    assert tiny.forward_ratio == pytest.approx(base.forward_ratio, rel=2.0e-12)
    assert tiny.normalized_multiplier == pytest.approx(base.normalized_multiplier, rel=2.0e-12)
    assert tiny.phase_alignment == pytest.approx(base.phase_alignment, abs=2.0e-12)
    assert tiny.signed_child_energy_work / base.signed_child_energy_work == pytest.approx(
        scale, rel=2.0e-12
    )


def test_absolute_triad_tolerance_cannot_accept_a_foreign_tiny_child():
    x, y, z = _edge(1.0e-12)
    foreign = z + np.asarray((1.5e-12, 0.0, 0.0))
    with pytest.raises(ValueError, match=r"z=x\+y|parent pair"):
        register_helical_physical_edge(
            x=x,
            y=y,
            z=foreign,
            sx=1,
            sy=-1,
            sz=1,
            ax=1.0,
            ay=1.0,
            az=1.0j,
        )


def test_typed_edge_cannot_forge_the_physical_registration_identity():
    row = _registration()
    with pytest.raises((ValueError, AssertionError), match="identity|registration"):
        replace(
            row,
            registered_upper_progress_work=row.registered_upper_progress_work
            + 0.25 * row.native_modal_capacity * row.global_multiplier_Jstar,
            upper_progress_identity_residual=0.0,
        )


def test_typed_edge_cannot_forge_direct_waleffe_equality():
    row = _registration()
    with pytest.raises((ValueError, AssertionError), match="Waleffe|source|coefficient"):
        replace(
            row,
            waleffe_child_source_coefficient=row.waleffe_child_source_coefficient + 1.0,
            direct_vs_waleffe_residual=0.0,
        )


def test_nonzero_modal_capacity_must_not_silently_underflow_to_zero():
    x, y, z = _edge()
    with pytest.raises(ValueError, match="capacity|native.*range|underflow"):
        register_helical_physical_edge(
            x=x,
            y=y,
            z=z,
            sx=1,
            sy=-1,
            sz=1,
            ax=1.0e-120,
            ay=1.0e-120,
            az=1.0e-120,
        )
