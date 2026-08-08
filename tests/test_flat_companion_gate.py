import math
from src.flat_companion_gate import classify_flat_companion, cap_mass_lower_from_barycenter


def test_low_barycenter_pays_entropy():
    out=classify_flat_companion(1e-6,.9,.999)
    assert out['branch']=='old_parent_collision_entropy'
    assert out['entropy_lower']>1/200


def test_two_concentrated_marginals_force_separated_cores():
    H=(.01/3)**2
    out=classify_flat_companion(H,.995,.997)
    assert out['branch']=='two_trackable_parent_cores'
    assert out['barycenter_direction_angle_lower']>1
    assert out['mass_each_cap_lower']>=7/9
    assert out['cap_gap_chord_lower']>1/3


def test_markov_cap_mass_exact_clean_value():
    assert abs(cap_mass_lower_from_barycenter(.99,.3)-7/9)<1e-14
