import math
import numpy as np
from src.sideband_coherence_daughter import (
    arb_sideband_capacity_certificate,
    continuous_variation_margin,
    first_duhamel_dichotomy,
    h3_forcing_vector,
)
from src.hermite_helicity_ledger import sym3


def test_continuous_variation_examples():
    # Constant forcing saturates through the daughter term.
    assert continuous_variation_margin(2.0,2.0,0.0,1.0)==0.0
    # A perfectly reversing forcing can have zero impulse but must pay variation.
    assert continuous_variation_margin(1.0,0.0,2.0,1.0)==0.0


def test_duhamel_dichotomy_branches():
    assert first_duhamel_dichotomy(1.0,0.7,0.8,1.0)['branch']=='coherent_daughter'
    assert first_duhamel_dichotomy(1.0,0.1,1.8,1.0)['branch']=='dephasing_source'


def test_h3_forcing_norm_isometry():
    rng=np.random.default_rng(4); T=sym3(rng.normal(size=(3,3,3)))
    f=h3_forcing_vector(T)
    assert math.isclose(np.dot(f,f),3/8*np.sum(T*T),rel_tol=1e-13,abs_tol=1e-13)


def test_arb_sideband_constants():
    cert=arb_sideband_capacity_certificate()
    assert cert['H3_clean_reverse_constant']=='1/160'
    assert cert['H1_clean_reverse_constant']=='1/16'
