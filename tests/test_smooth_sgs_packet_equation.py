import numpy as np
from src.smooth_sgs_packet_equation import (
    selected_parent_lower_ratio, low_low_output_radius,
    moving_multiplier_residual, transported_gaussian_symbol,
    source_taxonomy_is_disjoint,
)
from src.single_edge_certificate import float_rstar


def test_good_core_parents_clear_three_fifths():
    assert selected_parent_lower_ratio(float_rstar()) > 3/5


def test_quarter_transport_lowpass_cannot_make_selected_role():
    assert low_low_output_radius(1/4) == 1/2
    assert 1/2 < 3/5


def test_affine_transport_heisenberg_residual():
    A=np.array([[.2,.3,0.],[-.1,-.4,.2],[.1,.2,.2]])
    A-=np.trace(A)/3*np.eye(3)
    xi=np.array([.7,-.4,1.1]); C=np.diag([1.2,.8,.5])
    _,grad,dt=transported_gaussian_symbol(A,.03,xi,C)
    assert abs(moving_multiplier_residual(A,xi,grad,dt))<1e-12


def test_micro_macro_source_ledgers_are_disjoint():
    assert source_taxonomy_is_disjoint()
