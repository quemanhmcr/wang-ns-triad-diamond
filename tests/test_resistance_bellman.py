import math
import numpy as np
from src.resistance_bellman import Edge, normalize_edges, tree_resistance_matrix, poisson_certificate, component_collision

def test_two_vertex_exact():
    edges=normalize_edges([Edge(0,1,2.0)])
    p=np.array([0.5,0.5]); R=tree_resistance_matrix(2,edges); lam=1.0
    cert=poisson_certificate(edges,p,R,lam)
    assert cert['edgeQ'] <= cert['edgeQbound']+1e-12
    # direct expectation of same-component collision under independent Poisson cut
    qcut=1-math.exp(-1/(edges[0].c*lam))
    q_keep=1-qcut
    q0,_=component_collision(2,edges,0,p)
    q1,_=component_collision(2,edges,1,p)
    assert abs(q_keep*q0+qcut*q1-cert['expQ'])<1e-12

def test_edge_entropy_bound_random():
    rng=np.random.default_rng(3)
    for _ in range(1000):
        w=rng.dirichlet(np.ones(rng.integers(2,30)))
        edges=[Edge(i,i+1,float(x)) for i,x in enumerate(w)]
        lam=float(10**rng.uniform(0,3))
        rho=sum(min(e.c,1/lam) for e in edges)
        q=sum(e.c**2 for e in edges)
        assert q <= 1-rho+1/lam+1e-12

def test_simultaneous_witness_small_path():
    edges=normalize_edges([Edge(0,1,1),Edge(1,2,0.2),Edge(2,3,1)])
    p=np.ones(4)/4; R=tree_resistance_matrix(4,edges); cert=poisson_certificate(edges,p,R,5.0)
    found=False
    for mask in range(8):
        q,c=component_collision(4,edges,mask,p)
        if q<=cert['qbound']+1e-12 and c<=cert['cbound']+1e-12: found=True
    assert found
