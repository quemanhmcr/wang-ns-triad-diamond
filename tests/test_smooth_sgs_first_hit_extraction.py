import math

import numpy as np

from src.smooth_sgs_first_hit_extraction import (
    PhysicalPathMonitor,
    ThresholdTopology,
    circle_geodesic_holonomy,
    first_physical_corridor_exit,
    moyal_cell_energy_rate_identity,
    registration_no_hit_exhaustion,
    rescale_monitor_units,
    superlevel_debut_piecewise_linear,
    transported_moyal_cell_rate_upper,
)


def test_closed_and_strict_superlevel_debuts_distinguish_boundary_plateau():
    t=(0.0,1.0,2.0,3.0)
    v=(0.0,1.0,1.0,2.0)
    assert superlevel_debut_piecewise_linear(t,v,1.0,ThresholdTopology.CLOSED)==1.0
    assert superlevel_debut_piecewise_linear(t,v,1.0,ThresholdTopology.STRICT)==2.0


def test_joint_first_exit_keeps_exact_tie_without_priority():
    t=(0.0,0.5,1.0)
    a=PhysicalPathMonitor("strain",1.0,(0.0,1.0,2.0),ThresholdTopology.CLOSED)
    b=PhysicalPathMonitor("source",7.0,(6.0,7.0,8.0),ThresholdTopology.CLOSED)
    out=first_physical_corridor_exit(t,(b,a))
    assert out.first_time==0.5
    assert out.joint_causes==("source","strain")


def test_independent_change_of_units_leaves_first_exit_unchanged():
    t=(0.0,0.5,1.0)
    a=PhysicalPathMonitor("a",2.0,(0.0,2.0,4.0))
    b=PhysicalPathMonitor("b",3.0,(0.0,1.0,4.0))
    ref=first_physical_corridor_exit(t,(a,b))
    changed=first_physical_corridor_exit(t,(rescale_monitor_units(a,1e9),rescale_monitor_units(b,1e-7)))
    assert ref.first_time==changed.first_time
    assert ref.joint_causes==changed.joint_causes


def test_material_moyal_cell_energy_rate_is_exactly_cauchy_schwarz_controlled():
    F=np.array([1+2j,-3+1j,2-1j],complex)
    dF=np.array([2-1j,1+4j,-2+3j],complex)
    out=moyal_cell_energy_rate_identity(F,dF,np.array([True,False,True]))
    expected=2*float(np.real(np.vdot(F[[0,2]],dF[[0,2]])))
    assert math.isclose(out["cell_energy_derivative"],expected)
    assert out["margin"]>=-1e-13


def test_transport_moyal_rate_bound_has_correct_product_structure():
    upper=transported_moyal_cell_rate_upper(cell_energy=4.0,u_rate_l2=3.0,u_l2=5.0,window_rate_l2=0.2)
    assert math.isclose(upper,2*2*(3+1))


def test_circle_geodesic_has_no_principal_arg_branch_jump():
    eps=1e-7
    a=complex(math.cos(math.pi-eps),math.sin(math.pi-eps))
    b=complex(math.cos(-math.pi+eps),math.sin(-math.pi+eps))
    assert math.isclose(circle_geodesic_holonomy(a),circle_geodesic_holonomy(b),abs_tol=1e-12)


def test_no_backward_obstruction_is_registered_generated_survivor():
    z_event=1+0j
    i_r=.1j
    i_hh=.2+0j
    z_slice=z_event-i_r-i_hh
    out=registration_no_hit_exhaustion(z_event,z_slice,i_hh,i_r)
    assert out["classification"]=="registered_generated_survivor"
    assert out["registered_amplitude_lower"]>=0.25


def test_backward_material_relink_is_named_stop_not_survivor():
    z_event=1+0j
    i_r=0j
    i_hh=0j
    z_slice=1+0j
    out=registration_no_hit_exhaustion(z_event,z_slice,i_hh,i_r,material_relink=True)
    assert out["classification"]=="named_backward_physical_stop"
    assert "material_relink" in out["stop_causes"]
