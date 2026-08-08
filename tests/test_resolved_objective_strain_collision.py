import math
import numpy as np

from src.resolved_objective_strain_collision import (
    l32_derivative_to_linf_constant,
    pressure_hessian_clean_routes,
    quadratic_source_enstrophy_lower,
    resolved_corotational_strain_rhs,
    resolved_l3_squared_mass_constant,
    resolved_material_gradient_rhs,
    sgs_gradient_stress_lower,
    sym,
    viscous_source_enstrophy_lower,
)


def test_resolved_corotational_identity():
    A=np.array([[1.,2.,0.],[-1.,-.5,.2],[.1,.3,-.5]])
    Hp=sym(np.array([[.2,1.,0.],[0.,-.1,.4],[.2,.1,-.1]]))
    GR=np.arange(9,dtype=float).reshape(3,3)/10
    DA=np.flipud(GR); nu=.7
    rhs=resolved_material_gradient_rhs(A,Hp,GR,DA,nu)
    S=sym(A); O=.5*(A-A.T)
    direct=sym(rhs)+S@O-O@S
    assert np.linalg.norm(direct-resolved_corotational_strain_rhs(A,Hp,GR,DA,nu))<1e-13


def test_clean_order2_and_l3_constants_are_below_targets():
    assert l32_derivative_to_linf_constant(2,.5)<1/380
    assert resolved_l3_squared_mass_constant(.25)<1/15


def test_clean_routes():
    rho=.003
    assert math.isclose(quadratic_source_enstrophy_lower(rho),96*math.pi**2*rho)
    assert math.isclose(sgs_gradient_stress_lower(rho),380*rho)
    r=pressure_hessian_clean_routes(rho)
    assert r['resolved_critical_mass']==2850*rho
    assert r['stress_l32']==190*rho
    assert viscous_source_enstrophy_lower(rho,.5)>0
