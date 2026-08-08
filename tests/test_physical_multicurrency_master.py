import math

from src.physical_multicurrency_master import (
    additive_reset_count_upper,
    multicurrency_log_efficiency_lower,
    transfer_costly_count_lower,
)


def test_additive_resource_count_is_single_charge_sum():
    assert math.isclose(additive_reset_count_upper({'E':(3.,.5),'D':(2.,.25)}),14.)


def test_transfer_count_formula_on_exact_boundary():
    # L=10, NA=2, NT=3, NF=5 and 5*.2=(3+2+1)*.2-0.2 is admissible with Z=0 here.
    low=transfer_costly_count_lower(10,2,.2,.2,0.)
    assert low<=3+1e-14


def test_multicurrency_rate_is_original_master_rate():
    out=multicurrency_log_efficiency_lower(
        depth=100,budgets={'energy':(2.,.5)},c0=.01,kappa0=.17,potential_reset=.2,potential_error=.1,xi=.05
    )
    assert math.isclose(out['asymptotic_rate'],.01*.17/(.17+.2))
    assert out['transfer_costly_lower']>0
