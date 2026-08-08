import numpy as np
from src.adjoint_kelvin_duhamel import (
    adjoint_pairing_derivative_residual, exact_piecewise_duhamel,
    inherit_or_generate_route, phase_aligned_positive_generation,
    signed_good_lifetime_ratio_bounds, common_parent_natural_window_lower,
)


def test_adjoint_pairing_cancels_any_linear_generator():
    G=np.array([[1.,2.],[-3.,.5]],complex); c=np.array([1+2j,-.5j]); psi=np.array([.3-.1j,2+1j]); F=np.array([.2j,-1.])
    assert abs(adjoint_pairing_derivative_residual(G,c,psi,F))<1e-13


def test_exact_piecewise_duhamel():
    G=[np.array([[-1.,.2],[.2,-.5]],complex)]
    out=exact_piecewise_duhamel(np.array([1.+0j,.2j]),np.array([.5+.1j,1.-.2j]),G,[.3],[np.array([.2+.3j,-.1j])],[np.array([.01j,.02])])
    assert abs(out['duhamel_residual'])<1e-13


def test_clean_inherit_generate_triangle():
    # exact z1=z0+hh+r, inherited and residual below A/4 forces hh>A/2
    z0=.1+0j; r=.1j; hh=1-z0-r
    out=inherit_or_generate_route(1+0j,z0,r,hh)
    assert out['branch']=='high_high_generation'
    assert out['value']>=.5


def test_phase_aligned_positive_generation():
    z=[1+1j,-.4+.2j,.2-.1j]
    out=phase_aligned_positive_generation(z)
    assert out['positive_mass']+1e-14>=abs(sum(z))
    assert abs(np.sum(out['weights'])-1)<1e-14


def test_signed_good_lifetime_and_half_slab_overlap():
    lo,hi=signed_good_lifetime_ratio_bounds()
    assert lo.numerator==64 and lo.denominator==25
    assert hi.numerator==25 and hi.denominator==9
    assert abs(common_parent_natural_window_lower(1.0,.5)-103/50)<1e-14
