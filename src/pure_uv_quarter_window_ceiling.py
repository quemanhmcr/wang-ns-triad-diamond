from __future__ import annotations
import math
STATUS='DRAFT_PURE_UV_QUARTER_WINDOW_PRE_SINGULAR_CEILING__PTIME_AT_LEAST_ONE_QUARTER__H2_TAIL_EXCLUDES_HIGH_DOUBLING'
def pure_uv_time_fraction_lower()->float: return 0.25
def pure_uv_shell_mass_lower(parent_frequency:float,physical_tail_dissipation_lower:float,viscosity:float,global_energy:float,scaled_lifetime:float)->float:
    N=float(parent_frequency); D=float(physical_tail_dissipation_lower); nu=float(viscosity); E=float(global_energy); c=float(scaled_lifetime)
    if min(N,D,nu,E,c)<=0 or not all(math.isfinite(x) for x in (N,D,nu,E,c)): raise ValueError('positive finite physical inputs required')
    root=pure_uv_time_fraction_lower()*nu*D/(9.0*c*math.sqrt(math.pi)*N*E)
    return root*root
def presingular_h2_shell_mass_upper(parent_frequency:float,h2_seminorm_sq_upper:float)->float:
    N=float(parent_frequency); H=float(h2_seminorm_sq_upper)
    if min(N,H)<=0 or not all(math.isfinite(x) for x in (N,H)): raise ValueError('positive finite N,H2 bound required')
    return 2.0*H/(N**3)
def pure_uv_parent_frequency_ceiling(physical_tail_dissipation_lower:float,viscosity:float,global_energy:float,scaled_lifetime:float,h2_seminorm_sq_upper:float)->float:
    D=float(physical_tail_dissipation_lower); nu=float(viscosity); E=float(global_energy); c=float(scaled_lifetime); H=float(h2_seminorm_sq_upper)
    if min(D,nu,E,c,H)<=0 or not all(math.isfinite(x) for x in (D,nu,E,c,H)): raise ValueError('positive finite ceiling data required')
    C=nu*D/(36.0*c*math.sqrt(math.pi)*E)
    return 2.0*H/(C*C)
def theorem_certificate():
    return {'status':STATUS,'measure_law':'M=2N gives T_child=T_parent/4; four quarter windows cover the parent slab, so the supremal sliding positive-work fraction p_time>=1/4','mass_lower':'sqrt(mu)>=nu D_tail/(36 c sqrt(pi) N E_global)','smooth_tail':'finite pre-singular H2 gives mu_(2N)<=2||u||_Hdot2^2/N^3','ceiling':'N<=2592 pi c^2 E_global^2 H2/(nu^2 D_tail^2)','scope':'eventwise pure-UV ceiling; deep resolved-contact upward work remains separate'}
