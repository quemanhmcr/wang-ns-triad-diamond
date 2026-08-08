import math
import numpy as np
from src.renyi_causal_reuse import layer_collision_reuse,layered_collision_reuse,rich_layer_route


def test_injective_binary_slot_baseline():
    w=np.array([.6,.4]); pairs=np.array([[0,1],[2,3]])
    r=layer_collision_reuse(w,pairs,4)
    assert abs(r['q_parent']-.5*r['q_child'])<1e-14
    assert abs(r['hidden_parent_slot_pair_mass'])<1e-14
    assert abs(r['theta'])<1e-14


def test_reused_parent_increases_collision():
    w=np.array([.5,.5]); pairs=np.array([[0,1],[0,1]])
    r=layer_collision_reuse(w,pairs,2)
    assert r['hidden_parent_slot_pair_mass']>0
    assert r['theta']>0


def test_layered_action_telescope():
    maps=[np.array([[0,1],[1,2]]),np.array([[0,1]])]
    out=layered_collision_reuse(maps)
    assert abs(out['total_action']-(2*math.log(2)+math.log(out['root_collision'])))<1e-13


def test_rich_low_entropy_gives_parent_slot_pair():
    w=np.array([1.0]); pairs=np.array([[0,0]])
    out=rich_layer_route(w,pairs)
    assert out['branch']=='weighted_parent_slot_reuse_pair'
    assert out['hidden_pair_mass']>1/12
