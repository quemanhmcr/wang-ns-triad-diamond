import itertools
import math

import numpy as np

from src.helical import (
    check_helical_eigenvector,
    coupling_g,
    coupling_magnitude_closed,
    diamond_metrics,
    edge_metrics,
)


def test_helical_eigenvectors_and_reality_relation():
    for k in [np.array([1.2, -0.7, 2.3]), np.array([-2.0, 0.4, 0.9])]:
        for s in (-1, 1):
            assert check_helical_eigenvector(k, s)


def test_coupling_matches_closed_magnitude():
    k = np.array([1.0, 0.0, 0.0])
    p = np.array([-0.2, 1.1, 0.0])
    q = -(k+p)
    nk, np_, nq = map(np.linalg.norm, (k,p,q))
    for sk,sp,sq in itertools.product((-1,1), repeat=3):
        direct = abs(coupling_g(k,p,q,sk,sp,sq))
        closed = coupling_magnitude_closed(nk,np_,nq,sk,sp,sq)
        assert math.isclose(direct, closed, rel_tol=3e-10, abs_tol=3e-10)


def test_energy_and_helicity_coefficients_cancel():
    k,p,q = 0.8,1.1,1.5
    for sk,sp,sq in itertools.product((-1,1), repeat=3):
        ak = sp*p-sq*q
        ap = sq*q-sk*k
        aq = sk*k-sp*p
        assert abs(ak+ap+aq) < 1e-12
        assert abs(sk*k*ak + sp*p*ap + sq*q*aq) < 1e-12


def test_symmetric_reference_near_reported_value():
    x=0.610904
    q=np.array([1.0,0.0,0.0])
    theta=math.acos(1/(2*x*x)-1)
    a=x*np.array([math.cos(theta/2), math.sin(theta/2),0])
    b=x*np.array([math.cos(theta/2),-math.sin(theta/2),0])
    z=a+b
    e=edge_metrics(a,b,z,1,-1,1)
    assert abs(e.efficiency-0.100110) < 5e-4


def test_diamond_phase_frustration_is_finite():
    a=np.array([1.,0.,0.])
    b=np.array([0.2,1.1,0.])
    c=np.array([0.3,-0.1,1.2])
    m=diamond_metrics(a,b,c,(1,-1,1,-1,1,-1))
    assert 0 <= m["phase_frustration"] <= math.pi


def test_exact_scale_holonomy_identity():
    gamma=0.493
    la,lb,lc,lm=0.17,-0.31,0.44,1.08
    r1=la-lb
    r2=lm-(la+lb)/2-gamma
    r3=lm-lc
    r4=lb-lc
    assert abs((r2-r3+0.5*r1+r4)+gamma) < 1e-12
