import math

import numpy as np

from src.resolved_cutoff_repartition_relay import (
    cutoff_repartition_residual,
    renewed_lifetime_ratio,
    renewed_low_low_gap,
    renewed_parent_role_lower_relative_to_parent,
)


def test_cutoff_change_is_exact_repartition_of_same_quadratic_ns_term():
    rng=np.random.default_rng(7)
    n=4
    T=rng.normal(size=(n,n,n))+1j*rng.normal(size=(n,n,n))
    u=rng.normal(size=n)+1j*rng.normal(size=n)
    V0=rng.normal(size=n)+1j*rng.normal(size=n)
    V1=rng.normal(size=n)+1j*rng.normal(size=n)
    Q=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
    out=cutoff_repartition_residual(T,u,V0,V1,Q)
    assert out['old_to_new']<1e-11
    assert out['old_to_full_ns']<1e-11
    assert out['new_to_full_ns']<1e-11


def test_parent_scale_renewal_keeps_smooth_carrier_above_new_lowlow_output():
    lower=renewed_parent_role_lower_relative_to_parent(5/8,1/30,1/30)
    assert lower>(1/2)
    assert renewed_low_low_gap(5/8,1/30,1/30)>0
    assert math.isclose(lower,(22/25)*math.exp(-1/15))


def test_parent_parabolic_lifetime_window_is_the_existing_signed_good_window():
    assert math.isclose(renewed_lifetime_ratio(5/8),64/25)
    assert math.isclose(renewed_lifetime_ratio(3/5),25/9)
