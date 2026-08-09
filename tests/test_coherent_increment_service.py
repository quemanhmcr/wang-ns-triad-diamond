import math
import numpy as np
from src.coherent_increment_service import (
    coherent_service_route, cubic_to_square_threshold, discrete_cell_increment_energies,
    exact_certificate, low_square_service_lower, periodic_increment_covariance_residual,
)


def test_cubic_square_homogeneity():
    y=cubic_to_square_threshold(8.0,1.0,2.0,1.0)
    assert math.isclose(y,1.0)


def test_translation_covariance_discrete():
    rng=np.random.default_rng(4); n=21
    f=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g)
    assert periodic_increment_covariance_residual(f,g,3)<1e-11


def test_local_cell_increment_capacity():
    rng=np.random.default_rng(5); n=17
    f=rng.normal(size=n)+1j*rng.normal(size=n); g=rng.normal(size=n)+1j*rng.normal(size=n); g/=np.linalg.norm(g)
    labels=np.indices((n,n)).sum(axis=0)%4
    inc,here,nbr=discrete_cell_increment_energies(f,g,2,labels)
    assert np.all(inc<=2*(here+nbr)+1e-11)


def test_clean_interface_route():
    Y=1.0; d=.1
    assert low_square_service_lower(Y,d)>=.5
    r=coherent_service_route(Y,d,.1,.08,.2,[.55],['A'])
    assert r['branch']=='selected_interface_Xi'


def test_clean_dominant_new_cluster_mass():
    Y=1.; d=.1
    r=coherent_service_route(Y,d,.1,.05,.05,[.6,.05,.05],['A','B','C'])
    assert r['branch']=='dominant_new_coherent_cluster'
    assert r['coherent_critical_mass_lower']>=1/32


def test_entropy_or_cycle_constants():
    Y=1.; d=.1
    w=[.10]*8
    r=coherent_service_route(Y,d,.05,.02,.02,w,list(range(8)))
    assert r['branch']=='new_service_Bellman_entropy'
    assert r['H_ancestry']>=math.log(2)-1e-12


def test_canonical_lp_registration_keeps_observable_and_pde_tail_distinct():
    cert=exact_certificate()
    assert 'square-normalized' in cert['canonical_lp_registration']
    assert 'D_tail>=D_high/4' in cert['canonical_lp_registration']
    assert 'smooth-LP service observable' in cert['observer_pde_separation']
    assert 'orthogonal PDE dissipation currency' in cert['observer_pde_separation']
