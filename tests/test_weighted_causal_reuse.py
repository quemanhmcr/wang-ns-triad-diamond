import math
import numpy as np
from src.weighted_causal_reuse import layer_reuse_information,layered_reuse_information,random_surjective_map,rich_layer_threshold


def test_injective_parent_map_has_zero_reuse_information():
    w=np.array([.7,.3]); pairs=np.array([[0,1],[2,3]])
    out=layer_reuse_information(w,pairs,4)
    assert abs(out['reuse_information'])<1e-14


def test_complete_parent_reuse_has_positive_information():
    w=np.array([.5,.5]); pairs=np.array([[0,1],[0,1]])
    out=layer_reuse_information(w,pairs,2)
    assert out['reuse_information']>0


def test_layered_telescope():
    maps=[np.array([[0,1],[1,2]]),np.array([[0,1]])]
    out=layered_reuse_information(maps)
    assert abs(out['total_reuse_information']-(2*math.log(2)-out['root_entropy']))<1e-13


def test_clean_threshold_is_log_four_thirds():
    assert abs(rich_layer_threshold()-math.log(4/3))<1e-15
