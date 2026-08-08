import math
from src.causal_binary_ancestry import (
    binary_product_ratio,binary_reuse_action,causal_cycle_identity,quarter_reuse_depth,reuse_fractions,reuse_savings,
)


def test_full_binary_tree_has_zero_reuse_and_cycle_rank():
    n=[8,4,2,1]
    assert reuse_savings(n)==[0,0,0]
    assert causal_cycle_identity(n)==(0,0)
    assert binary_product_ratio(n)==1
    assert binary_reuse_action(n)==0


def test_layer_reuse_cycle_identity():
    n=[4,3,2,1]
    beta,s=causal_cycle_identity(n)
    assert beta==s==3
    assert abs(binary_reuse_action(n)-(3*math.log(2)-math.log(4)))<1e-14


def test_product_identity():
    n=[3,2,1]
    rho=reuse_fractions(n)
    assert abs((1-rho[0])*(1-rho[1])-3/4)<1e-14


def test_quarter_stop_is_finite():
    assert quarter_reuse_depth(1.0,1.0)>0
