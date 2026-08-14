import math
import numpy as np
import pytest

from src.critical_unique_donor_pushforward import (
    critical_donor_scale_measure,
    critical_efficient_recipient_donor_pushforward,
    theorem_certificate,
)
from src.curl_spectral_curvature_balance import sharp_critical_frontier_geometry
from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel, register_closed_helical_triad
from src.helical import coupling_g


def _sharp_triad(q=1.0):
    out=sharp_critical_frontier_geometry(); D=out.same_helicity_parent_ratio; S=out.opposite_helicity_parent_ratio
    c=(1-D*D-S*S)/(2*D*S); ss=math.sqrt(1-c*c)
    x=np.array([D,0.,0.]); y=np.array([S*c,S*ss,0.]); z=x+y
    g=coupling_g(x,y,-z,1,-1,1)
    triad=register_closed_helical_triad(wavevectors=(-z,x,y),helicities=(1,1,-1),amplitudes=(1.,g/abs(g),1.))
    recipient=triad.slot_for_edge_child_wavevector(z).closed_mode_index
    return triad,cyclic_triad_measure_kernel(triad,quotient_measure_mass=q),recipient


def test_sharp_critical_recipient_pushes_to_one_lower_scale_energy_donor_without_rehahn():
    triad,kernel,i=_sharp_triad(2.0)
    out=critical_efficient_recipient_donor_pushforward(triad,kernel,recipient_closed_mode_index=i)
    assert out.critical_efficiency > out.critical_efficiency_threshold
    assert out.donor_child_ratio < 5/8
    assert out.canonical_recipient_positive_mass == pytest.approx(out.pushed_donor_mass,rel=3e-10)
    assert not out.later_hahn_used
    assert not out.argmax_scale_selector_used


def test_full_donor_scale_distribution_is_mass_preserving_not_argmax_selected():
    rows=[]
    for q in (0.3,1.0,2.7):
        triad,kernel,i=_sharp_triad(q)
        rows.append(critical_efficient_recipient_donor_pushforward(triad,kernel,recipient_closed_mode_index=i))
    law=critical_donor_scale_measure(rows)
    assert law.total_recipient_mass == pytest.approx(law.total_donor_pushforward_mass,rel=3e-10)
    assert law.maximum_donor_child_ratio < 5/8
    assert len(law.atoms)==3
    assert not law.argmax_scale_selector_used


def test_certificate_keeps_between_time_continuity_as_real_remaining_problem():
    cert=theorem_certificate()
    assert "no argmax" in cert["distribution"]
    assert "between-time continuation" in cert["scope"]
