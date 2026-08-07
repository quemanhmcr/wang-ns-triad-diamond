import math
from src.h1_covariant_source import (
    arb_source_constants_certificate,
    clean_interaction_variation_upper,
    h1_dominant_strain_or_source,
)


def test_clean_source_bound():
    assert clean_interaction_variation_upper(2.0, 3.0) == 166.0


def test_h1_dominant_source_or_strain():
    I1=1.0; IB=2.0; T=1.0
    out=h1_dominant_strain_or_source(I1,IB,T,I1/44,0.0)
    assert out['branch']=='curvature_source'
    out=h1_dominant_strain_or_source(I1,IB,T,0.0,1/2376)
    assert out['branch']=='base_strain_action'


def test_arb_constants_optional():
    import pytest
    pytest.importorskip('flint')
    c=arb_source_constants_certificate(); assert c['clean_connection_coefficient']=='54'
