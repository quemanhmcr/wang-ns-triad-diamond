from fractions import Fraction
from src.coherent_service_stopping import (
    clean_old_pool_ratio, epoch_certificate, first_forced_cost_generation,
)


def test_clean_ratio_is_below_half():
    assert clean_old_pool_ratio()==Fraction(231525,512000)
    assert clean_old_pool_ratio()<Fraction(1,2)


def test_immediate_stop_when_old_pool_already_small():
    assert first_forced_cost_generation(1.0,.1)==0


def test_stopping_generation_is_minimal():
    Y=.2; C0=10.; q=first_forced_cost_generation(Y,C0); r=float(clean_old_pool_ratio()); target=Y/8
    assert C0*r**q<=target
    assert q==0 or C0*r**(q-1)>target


def test_epoch_certificate_names_clean_routes():
    c=epoch_certificate(1.0,8.0)
    assert 'Y/32' in c['forced_alternatives']
    assert c['first_forced_generation']>0
