import math
import numpy as np
from src.sideband_pair_ancestry import (
    ancestry_weights_equal_edge_class_mass, collision_chain,
    pair_rescue_ancestry_route, same_ancestry_endpoint_law, split_pair_rescue,
)


def test_endpoint_ancestry_equals_edge_class_mass():
    edges=np.array([[0,1],[1,2],[2,3],[0,3]])
    weights=np.array([2.,1.,3.,4.]); anc=np.array([0,0,1,1])
    assert ancestry_weights_equal_edge_class_mass(weights,edges,anc)<1e-14
    w,WA,_=same_ancestry_endpoint_law(weights,edges,anc)
    c=collision_chain(w,WA); assert c['Q_ancestry']>=c['Q_atomic']


def test_cross_branch():
    edges=np.array([[0,1],[0,3],[2,1],[2,3]]); weights=np.ones(4); anc=np.array([0,1,0,1])
    out=pair_rescue_ancestry_route(weights,edges,anc); assert out['branch']=='cross_Xi'
