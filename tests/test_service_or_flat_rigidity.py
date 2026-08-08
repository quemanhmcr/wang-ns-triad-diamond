from fractions import Fraction
from src.service_or_flat_rigidity import (
    CURVATURE_DENOM, arb_threshold_certificate, classify_service_or_flat,
    clean_thresholds, connection_flatness_upper,
)


def test_tau_one_percent_exact_delta():
    assert Fraction(1,CURVATURE_DENOM*10000)==Fraction(1,10_368_000_000_000)


def test_clean_flat_block_is_tau_flat():
    tau=.01; th=clean_thresholds(tau); d=th['block_transfer_deficit']*.5
    out=classify_service_or_flat(tau,d,d,d,.5*th['h3_source_impulse'],.5*th['h1_source_impulse'],.5*th['objective_strain_variation_action'],.5*th['low_strain_action'],1.02,True)
    assert out['branch']=='kelvin_extremal_flat'
    assert out['kelvin_connection_flatness']<=tau


def test_high_aspect_routes_to_inherited_or_fresh():
    tau=.01; th=clean_thresholds(tau); z=0.0
    a=classify_service_or_flat(tau,z,z,z,z,z,z,.01,1.2,True); b=classify_service_or_flat(tau,z,z,z,z,z,z,.01,1.2,False)
    assert a['branch']=='inherited_high_aspect'; assert b['branch']=='fresh_high_aspect'


def test_transfer_threshold_routes_before_flatness():
    tau=.01; th=clean_thresholds(tau)
    out=classify_service_or_flat(tau,th['block_transfer_deficit'],0,0,0,0,0,.01,1.0,True)
    assert out['branch']=='physical_transfer_cost'


def test_arb_optional():
    import pytest; pytest.importorskip('flint'); c=arb_threshold_certificate(); assert c['status'].startswith('CERTIFIED')
