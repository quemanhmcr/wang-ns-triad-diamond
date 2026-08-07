import math
import numpy as np

from src.multicommodity_hodge import (
    DirectedEdge, hodge_energy, multicommodity_rayleigh,
    synchronization_certificate, integer_gauge_mismatch_bound,
    tree_pair_resistance_identity,
)


def test_single_cycle_rayleigh_exact_when_cycle_spans_cycle_space():
    edges=[DirectedEdge(0,1,1.0,1.0),DirectedEdge(1,2,1.0,1.0),DirectedEdge(2,0,1.0,1.0)]
    z=np.array([1.0,1.0,1.0])
    cert=multicommodity_rayleigh(3,edges,[z],np.array([1.0]))
    assert abs(cert['margin']) < 1e-10
    assert cert['hodge_energy'] > 0


def test_gradient_target_zero_hodge():
    phi=np.array([0.0,1.0,3.0])
    edges=[DirectedEdge(0,1,1.2,1.0),DirectedEdge(1,2,0.7,2.0),DirectedEdge(0,2,2.0,3.0)]
    assert hodge_energy(3,edges)['energy'] < 1e-20


def test_synchronization_equal_shift_zero_variance():
    gamma=0.5
    nv=4
    old_pairs=[(0,2),(1,2),(2,3)]
    new_pairs=[(0,2),(1,2),(2,3)]
    h_old=gamma*np.array([0,0,1,2.],float)
    h_new=h_old+gamma
    cert=synchronization_certificate(nv,old_pairs,new_pairs,h_old,h_new,[0,1],np.array([0.4,0.6]))
    assert cert['shift_variance'] < 1e-20
    assert cert['hodge_energy'] < 1e-20


def test_integer_mismatch_probability_bound():
    gamma=0.7
    p=np.array([0.5,0.3,0.2])
    d=gamma*np.array([0,1,1],float)
    mean=(p*d).sum(); var=(p*(d-mean)**2).sum()
    exact=0.0
    for i in range(3):
        for j in range(3):
            if d[i] != d[j]: exact += p[i]*p[j]
    assert exact <= integer_gauge_mismatch_bound(var,gamma)+1e-12


def test_tree_pair_resistance_identity():
    out=tree_pair_resistance_identity(4,[(0,1,1.0),(1,2,2.0),(1,3,0.5)],np.array([0.4,0.3,0.2,0.1]))
    assert abs(out['margin']) < 1e-12
