import math
from src.coherent_sgs_episode import (
    coherent_sgs_episode_costs, coherent_square_service_per_sgs_source,
    source_homogeneity_residual,
)


def test_source_to_coherent_service_is_linear():
    args=(1.2,1.1,1.5,1.3)
    c=coherent_square_service_per_sgs_source(*args)
    assert c>0
    assert abs(source_homogeneity_residual(.037,*args))<1e-12


def test_source_weighted_constants():
    out=coherent_sgs_episode_costs(.1,.5,1.2,1.1,1.5,1.3,old_pool_capacity=.001)
    cy=out['coherent_square_service_per_source']; sig=out['total_source_weight']
    assert math.isclose(out['high_frequency_dissipation'],cy*sig/16)
    assert math.isclose(out['selected_interface_Xi'],cy*sig/32)
    assert math.isclose(out['integrated_new_coherent_critical_mass'],cy*sig/128)
    assert math.isclose(out['peak_new_coherent_critical_mass']*.5,out['integrated_new_coherent_critical_mass'])
