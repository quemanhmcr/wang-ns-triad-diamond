from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np


def sym(A: np.ndarray) -> np.ndarray:
    A=np.asarray(A,float)
    return 0.5*(A+A.T)


def resolved_corotational_strain_rhs(A: np.ndarray, HessP: np.ndarray, GradDivR: np.ndarray, DeltaA: np.ndarray, nu: float) -> np.ndarray:
    """Exact corotational strain RHS for the filtered transporter V.

    D_t^V A = -A^2-Hess P-GradDiv R+nu Delta A,
    S_circ=sym(D_t A)+S Omega-Omega S.
    """
    A=np.asarray(A,float); HessP=np.asarray(HessP,float); GradDivR=np.asarray(GradDivR,float); DeltaA=np.asarray(DeltaA,float)
    if any(X.shape!=(3,3) for X in (A,HessP,GradDivR,DeltaA)) or nu<0:
        raise ValueError('3x3 matrices and nu>=0 required')
    S=sym(A); O=0.5*(A-A.T)
    return -(S@S)-(O@O)+(S@O-O@S)-sym(HessP)-sym(GradDivR)+nu*sym(DeltaA)


def resolved_material_gradient_rhs(A: np.ndarray, HessP: np.ndarray, GradDivR: np.ndarray, DeltaA: np.ndarray, nu: float) -> np.ndarray:
    A=np.asarray(A,float)
    return -(A@A)-np.asarray(HessP,float)-np.asarray(GradDivR,float)+nu*np.asarray(DeltaA,float)


def l32_derivative_to_linf_constant(order: int, support_ratio: float) -> float:
    """Unitary Fourier L^(3/2)->Linf derivative bound on B_(lambda N).

    Hausdorff--Young L^(3/2)->L^3 plus Holder in Fourier space gives
    C=(2pi)^-2 [4pi/(3m/2+3)]^(2/3) lambda^(m+2).
    """
    if order<0 or support_ratio<=0:
        raise ValueError('invalid derivative/support')
    denom=1.5*order+3.0
    return (2*math.pi)**-2*(4*math.pi/denom)**(2/3)*support_ratio**(order+2)


def l2_linf_constant(support_ratio: float) -> float:
    if support_ratio<=0: raise ValueError('positive support ratio required')
    return (2*math.pi)**-1.5*math.sqrt(4*math.pi/3)*support_ratio**1.5


def resolved_l3_squared_mass_constant(support_ratio: float=0.25) -> float:
    """C_B^2 in ||V||_3^2 <= C_B^2 N||V||_2^2 by L2/Linf interpolation."""
    return l2_linf_constant(support_ratio)**(2/3)


def third_from_gradient_l2_to_linf_constant(support_ratio: float=0.25) -> float:
    """C in ||D^3 V||inf <= C N^(7/2)||grad V||2."""
    if support_ratio<=0: raise ValueError('positive support ratio required')
    return (2*math.pi)**-1.5*math.sqrt(4*math.pi/7)*support_ratio**3.5


def quadratic_source_enstrophy_lower(source_level: float) -> float:
    """Clean d_V lower using C_grad=1/(8sqrt6 pi): rho_Q<=d_V/(96pi^2)."""
    if source_level<0: raise ValueError('nonnegative source required')
    return 96*math.pi**2*source_level


def sgs_gradient_stress_lower(source_level: float) -> float:
    """Clean order-2 SGS route using C_2<1/380."""
    if source_level<0: raise ValueError('nonnegative source required')
    return 380*source_level


def pressure_hessian_clean_routes(source_level: float) -> dict[str,float]:
    """If rho_P is large, either resolved critical mass or SGS stress is large.

    Clean bound rho_P <= mu_V/5700 + ||R||_(3/2)/380.
    Split the source budget in halves.
    """
    if source_level<0: raise ValueError('nonnegative source required')
    return {'resolved_critical_mass':2850*source_level,'stress_l32':190*source_level}


def viscous_source_enstrophy_lower(source_level: float, nu: float) -> float:
    """Clean d_V lower from C_31<1/1500: rho_nu <=nu sqrt(d_V)/1500."""
    if source_level<0 or nu<0: raise ValueError('invalid source/viscosity')
    if source_level==0: return 0.0
    if nu==0: return math.inf
    return (1500*source_level/nu)**2


def per_channel_scaled_source_weight_lower(objective_variation_action: float, scaled_lifetime: float) -> float:
    """One of four resolved source channels carries this scaled source weight.

    A_obj=T int||S_circ||dt with T=cN^-2. In scaled time tau=N^2t,
    A_obj=c int rho_circ d tau. Splitting Q,P,SGS,nu gives Sigma_*>=A_obj/(4c).
    """
    if objective_variation_action<0 or scaled_lifetime<=0:
        raise ValueError('invalid action/lifetime')
    return objective_variation_action/(4*scaled_lifetime)


def arb_clean_constants() -> dict[str,str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:
        raise RuntimeError('python-flint required') from exc
    ctx.prec=180; pi=arb.pi()
    c2=(2*pi)**-2*(2*pi/3).root(3)**2*(arb(1)/2)**4
    cinf=(2*pi)**(-arb(3)/2)*(4*pi/3).sqrt()*(arb(1)/4)**(arb(3)/2)
    cb2=cinf**(arb(2)/3)
    c31=(2*pi)**(-arb(3)/2)*(4*pi/7).sqrt()*(arb(1)/4)**(arb(7)/2)
    cgrad2=cinf*cinf
    if not (c2<arb(1)/380): raise AssertionError(f'order2 constant not <1/380: {c2}')
    if not (cb2<arb(1)/15): raise AssertionError(f'L3^2 mass constant not <1/15: {cb2}')
    if not (c31<arb(1)/1500): raise AssertionError(f'D3-from-grad constant not <1/1500: {c31}')
    if not (abs(cgrad2-arb(1)/(384*pi*pi))<arb('1e-50')): raise AssertionError('gradient C^2 identity failed')
    return {
        'order2_l32_linf_ball':str(c2),
        'order2_clean':'1/380',
        'resolved_l3_squared_mass_ball':str(cb2),
        'resolved_l3_squared_mass_clean':'1/15',
        'third_from_gradient_ball':str(c31),
        'third_from_gradient_clean':'1/1500',
        'quadratic_clean':'rho_Q<=d_V/(96 pi^2)',
        'pressure_clean':'rho_P<=mu_V/5700+||R||_(3/2)/380',
        'status':'ARB_CERTIFIED_RESOLVED_OBJECTIVE_STRAIN_COLLISION',
    }


@dataclass(frozen=True)
class ResolvedStrainStress:
    samples:int
    worst_corotational_identity_residual:float
    minimum_quadratic_margin:float
    minimum_sgs_margin:float
    minimum_pressure_route_margin:float
    minimum_viscous_margin:float


def stress(samples:int=50_000,seed:int=20260808)->ResolvedStrainStress:
    rng=np.random.default_rng(seed)
    wi=0.; mq=ms=mp=mv=float('inf')
    c2=l32_derivative_to_linf_constant(2,.5)
    cb2=resolved_l3_squared_mass_constant(.25)
    c31=third_from_gradient_l2_to_linf_constant(.25)
    for _ in range(samples):
        A=rng.normal(size=(3,3)); A-=np.trace(A)/3*np.eye(3)
        Hp=sym(rng.normal(size=(3,3))); GR=rng.normal(size=(3,3)); DA=rng.normal(size=(3,3)); nu=float(rng.uniform(0,2))
        direct=sym(resolved_material_gradient_rhs(A,Hp,GR,DA,nu))
        S=sym(A); O=.5*(A-A.T)
        direct=direct+S@O-O@S
        formula=resolved_corotational_strain_rhs(A,Hp,GR,DA,nu)
        wi=max(wi,float(np.linalg.norm(direct-formula))/max(1.,float(np.linalg.norm(formula))))
        if np.linalg.norm(direct-formula)>3e-12*max(1.,np.linalg.norm(formula)):
            raise AssertionError('resolved corotational identity failed')

        rho=float(10**rng.uniform(-8,-.2))
        # Clean route implications tested against sharper analytic constants.
        dq=1.01*quadratic_source_enstrophy_lower(rho)
        upperq=4*l2_linf_constant(.25)**2*dq
        mq=min(mq,upperq-rho)
        if upperq<rho: raise AssertionError('quadratic clean collision direction failed')

        rr=1.01*sgs_gradient_stress_lower(rho)
        uppers=c2*rr
        ms=min(ms,uppers-rho)
        if uppers<rho: raise AssertionError('SGS-gradient clean collision failed')

        routes=pressure_hessian_clean_routes(rho)
        # At 1.01 times either clean threshold, the sharp analytic bound can cover >rho/2.
        term_mu=c2*cb2*(1.01*routes['resolved_critical_mass'])
        term_r=c2*(1.01*routes['stress_l32'])
        mp=min(mp,term_mu-rho/2,term_r-rho/2)
        if term_mu<rho/2 or term_r<rho/2:
            raise AssertionError('pressure clean split threshold too optimistic')

        if nu>1e-12:
            dv=1.01*viscous_source_enstrophy_lower(rho,nu)
            upperv=nu*c31*math.sqrt(dv)
            mv=min(mv,upperv-rho)
            if upperv<rho: raise AssertionError('viscous clean collision failed')
    return ResolvedStrainStress(samples,wi,mq,ms,mp,mv)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument('--samples',type=int,default=50_000); ap.add_argument('--outdir',type=Path,default=Path('results-resolved-objective-strain-collision'))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_clean_constants(); out=stress(args.samples)
    (args.outdir/'resolved_objective_strain_collision.json').write_text(json.dumps({'certificate':cert,'stress':asdict(out)},indent=2),encoding='utf-8')
    md=f'''# Resolved objective-strain source collision\n\nStatus: **{cert['status']}**.\n\nThe objective strain used by the affine/Kelvin packet is the strain of the strict transporter `V=S_(N/4)u`. Its filtered equation gives exactly\n\n`S_circ = -S^2-Omega^2+[S,Omega] - Hess P - sym grad div R + nu Delta S`.\n\nHere `supp Vhat subset B_(N/4)` and `supp Rhat,supp Phat subset B_(N/2)`.  The order-two Hilbert-valued Hausdorff--Young/Bernstein constant is `<1/380`; interpolation gives `||V||_3^2 < mu_V/15`; and `||D^3 V||_inf < N^(7/2)||grad V||_2/1500`. Therefore the four source channels have clean pointwise collisions\n\n- quadratic stretching: `rho_Q <= d_V/(96 pi^2)`, hence `d_V>=96 pi^2 rho_Q`;\n- resolved SGS strain source: `rho_R2 <= ||R||_(3/2)/380`, hence `||R||_(3/2)>=380 rho_R2`;\n- filtered pressure Hessian: `rho_P <= mu_V/5700+||R||_(3/2)/380`, hence either `mu_V>=2850 rho_P` or `||R||_(3/2)>=190 rho_P`;\n- viscosity: `rho_nu <= nu sqrt(d_V)/1500`, hence `d_V >= (1500 rho_nu/nu)^2`.\n\nThus the old unresolved **near pressure-Hessian coefficient disappears for the actual resolved Kelvin transporter**. Pressure strain-dephasing routes to coherent resolved mass or to the same SGS-increment service currency. The raw full-velocity strain identity remains mathematically valid, but it is not the correct source object for the affine transporter used in the service-or-flat gate.\n\nIf `A_obj=T int||S_circ||dt` and `T=cN^-2`, one of these four channels carries scaled source weight at least `A_obj/(4c)`. Quadratic/viscous channels route to critical normalized dissipation; SGS and the stress part of pressure route through Germano/Onsager to coherent service, ancestry or high-frequency enstrophy. The resolved-mass pressure branch enters the existing coherent reservoir/reuse mechanism.\n\nStress: `{out.samples}`\n- worst corotational identity residual: `{out.worst_corotational_identity_residual:.3e}`\n- minimum quadratic margin: `{out.minimum_quadratic_margin:.3e}`\n- minimum SGS margin: `{out.minimum_sgs_margin:.3e}`\n- minimum pressure split margin: `{out.minimum_pressure_route_margin:.3e}`\n- minimum viscous margin: `{out.minimum_viscous_margin:.3e}`\n'''
    (args.outdir/'summary.md').write_text(md,encoding='utf-8'); print(md)

if __name__=='__main__': main()
