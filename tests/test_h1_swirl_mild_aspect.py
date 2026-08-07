import math
from src.h1_swirl_mild_aspect import arb_isotropic_hook_certificate, mild_aspect_lower


def test_clean_mild_aspect_constant():
    assert mild_aspect_lower(21/20)>1/25


def test_arb_hook_certificate_optional():
    import pytest
    pytest.importorskip('flint')
    c=arb_isotropic_hook_certificate(); assert c['isotropic_Qpol_lower']=='1/10'; assert c['mild_Qpol_lower']=='1/25'
