from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.affine_gaussian_forcing import transform_hessian, whitened_velocity_hessian


def normalized_material_coordinate_rhs(L: np.ndarray, Ux: np.ndarray, U0: np.ndarray, A0: np.ndarray, y: np.ndarray) -> np.ndarray:
    """D_t z for z=L^-1 y, Xdot=U0, Ldot=A0 L and y=x-X."""
    L=np.asarray(L,float); y=np.asarray(y,float)
    return np.linalg.solve(L, np.asarray(Ux,float)-np.asarray(U0,float)-np.asarray(A0,float)@y)


def quadratic_velocity(U0: np.ndarray, A0: np.ndarray, H: np.ndarray, y: np.ndarray) -> np.ndarray:
    y=np.asarray(y,float)
    return np.asarray(U0,float)+np.asarray(A0,float)@y+0.5*np.einsum('ijk,j,k->i',np.asarray(H,float),y,y)


def affine_curvature_norm(L: np.ndarray, H: np.ndarray) -> float:
    return float(np.linalg.norm(whitened_velocity_hessian(L,H)))


def material_window_leakage_bound(kappa_aff: float, M: float, grad_base: float=1.0, support_radius: float=2.0) -> float:
    """|D_t chi(z/M)| <= (grad_base support_radius^2/2) kappa_aff M."""
    if min(kappa_aff,M,grad_base,support_radius)<0 or M==0:
        raise ValueError('nonnegative parameters and M>0 required')
    return 0.5*grad_base*support_radius**2*kappa_aff*M


def physical_window_gradient_bound(L: np.ndarray, M: float, grad_base: float=1.0) -> float:
    if M<=0: raise ValueError('M>0 required')
    return grad_base*float(np.linalg.norm(np.linalg.inv(np.asarray(L,float)),2))/M


def clean_shell_gradient_bound(N: float, M: float, grad_base: float=1.0) -> float:
    """N^-1 ||grad chi|| <= (3/2) grad_base/M from lmin>2/(3N)."""
    if N<=0 or M<=0: raise ValueError('positive N,M')
    return 1.5*grad_base/M


def convolution_commutator_bound(filter_first_moment: float, N: float, grad_chi: float, f_l2: float=1.0) -> float:
    """Young/MVT bound for [chi,G_N*]f with dimensionless first moment m1(G)."""
    if min(filter_first_moment,N,grad_chi,f_l2)<0 or N==0: raise ValueError('invalid parameters')
    return filter_first_moment*grad_chi*f_l2/N


def clean_affine_commutator_bound(filter_first_moment: float, M: float, grad_base: float=1.0, f_l2: float=1.0) -> float:
    """Shell-axis consequence: <=(3/2)m1(G) C_chi M^-1 ||f||_2."""
    if M<=0: raise ValueError('M>0 required')
    return 1.5*filter_first_moment*grad_base*f_l2/M


def affine_balance_optimum(a: float, b: float, kappa_aff: float) -> tuple[float,float]:
    if min(a,b,kappa_aff)<=0: raise ValueError('positive a,b,kappa_aff')
    M=math.sqrt(a/(b*kappa_aff))
    return M,2.0*math.sqrt(a*b*kappa_aff)


@dataclass(frozen=True)
class WindowStress:
    samples:int
    worst_material_identity_residual:float
    worst_taylor_leakage_ratio:float
    worst_shell_gradient_ratio:float
    worst_affine_curvature_residual:float
    worst_optimizer_residual:float
    extreme_condition_number:float
    extreme_leakage_ratio:float


def stress(samples:int=50_000,seed:int=20260807)->WindowStress:
    rng=np.random.default_rng(seed)
    wi=wc=wg=wa=wo=0.0
    for _ in range(samples):
        # SPD grain frame with broad aspect range.
        Q,_=np.linalg.qr(rng.normal(size=(3,3))); axes=np.exp(rng.uniform(-4,4,size=3));L=Q@np.diag(axes)@Q.T
        U0=rng.normal(size=3);A0=rng.normal(size=(3,3));H=rng.normal(size=(3,3,3));H=.5*(H+H.swapaxes(1,2))
        M=float(rng.uniform(.5,20.));R=float(rng.uniform(.5,2.)); z=rng.normal(size=3);z*=R*M/max(np.linalg.norm(z),1e-15);y=L@z
        Ux=quadratic_velocity(U0,A0,H,y)
        rhs=normalized_material_coordinate_rhs(L,Ux,U0,A0,y)
        B=whitened_velocity_hessian(L,H)
        expected=.5*np.einsum('abc,b,c->a',B,z,z)
        wi=max(wi,float(np.linalg.norm(rhs-expected))/max(1.,float(np.linalg.norm(expected))))
        # Any normalized base-window gradient with norm <=grad_base produces this leakage.
        g=rng.normal(size=3);g/=np.linalg.norm(g); grad_base=float(rng.uniform(.1,2.)); g*=grad_base/M
        actual=abs(float(np.dot(g,rhs)))
        bound=material_window_leakage_bound(float(np.linalg.norm(B)),M,grad_base,R)
        if actual>bound+2e-9: raise AssertionError('affine Taylor leakage bound failed')
        if bound>1e-14: wc=max(wc,actual/bound)

        N=10**float(rng.uniform(-2,4)); lmin=(2.0/(3.0*N))*float(rng.uniform(1.0001,5.0)); Q2,_=np.linalg.qr(rng.normal(size=(3,3))); a2=np.array([lmin,lmin*rng.uniform(1,30),lmin*rng.uniform(1,50)]);L2=Q2@np.diag(a2)@Q2.T
        actualg=physical_window_gradient_bound(L2,M,grad_base)/N
        clean=clean_shell_gradient_bound(N,M,grad_base)
        if actualg>clean+1e-9: raise AssertionError('shell gradient bound failed')
        if clean>1e-14: wg=max(wg,actualg/clean)

        S=Q@np.diag(np.exp(rng.uniform(-3,3,size=3)))@Q.T
        Bp=whitened_velocity_hessian(S@L,transform_hessian(S,H))
        wa=max(wa,float(np.linalg.norm(Bp-B))/max(1.,float(np.linalg.norm(B))))

        aa=float(rng.uniform(.1,4));bb=float(rng.uniform(.1,4));kap=10**float(rng.uniform(-7,-1));Mo,Eo=affine_balance_optimum(aa,bb,kap);direct=aa/Mo+bb*kap*Mo;wo=max(wo,abs(direct-Eo))

    # Extreme common affine squeeze: intrinsic leakage ratio stays bounded.
    A=1e6; S=np.diag([A,A**-.5,A**-.5]); L=np.diag([1.3,.9,1.1]);H=rng.normal(size=(3,3,3));H=.5*(H+H.swapaxes(1,2));B=whitened_velocity_hessian(L,H);Bp=whitened_velocity_hessian(S@L,transform_hessian(S,H))
    M=8.;R=2.;z=np.array([R*M,0.,0.]);rhs=.5*np.einsum('abc,b,c->a',Bp,z,z);g=np.array([1.,0.,0.])/M;bound=material_window_leakage_bound(float(np.linalg.norm(Bp)),M,1.,R);er=abs(float(g@rhs))/max(bound,1e-15)
    return WindowStress(samples,wi,wc,wg,wa,wo,float(np.linalg.cond(S@L)),er)


def main()->None:
    ap=argparse.ArgumentParser();ap.add_argument('--samples',type=int,default=50_000);ap.add_argument('--outdir',type=Path,default=Path('results-affine-window-balance'));args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    out=stress(args.samples)
    payload={'theorem':{
        'material_coordinate':'D_t z=L^-1[U(X+Lz)-U(X)-A(X)Lz]',
        'curvature':'kappa_aff=sup ||L^-1(nabla^2U)[L,L]||',
        'window_leakage':'|D_t chi(z/M)| <= (Cgrad Rchi^2/2) kappa_aff M',
        'shell_gradient':'N^-1||grad_x chi|| <= (3/2) Cgrad/M',
        'filter_commutator':'||[chi,G_N*]f||_2 <= m1(G) N^-1||grad chi||_inf||f||_2 <= (3/2)m1(G)Cgrad M^-1||f||_2',
        'balance':'a/M+b kappa_aff M, M*=sqrt(a/(b kappa_aff))',
    },'stress':asdict(out)}
    (args.outdir/'affine_window_balance.json').write_text(json.dumps(payload,indent=2))
    md=f"""# Affine ellipsoidal moving-window balance

Transport `Xdot=U(X)` and `Ldot=A(X)L`.  Then the normalized grain coordinate
obeys the exact non-affine remainder identity

`D_t z=L^-1[U(X+Lz)-U(X)-A(X)Lz]`.

For `chi(z/M)` this gives `O(kappa_aff M)` material leakage, while the certified
shell lower axis gives the complementary `O(1/M)` physical-gradient/commutator
scale.  Thus the affine window retains the same square-root balance
`a/M+b kappa_aff M` without an aspect-ratio penalty.

- random checks: `{out.samples}`
- worst material identity residual: `{out.worst_material_identity_residual:.3e}`
- worst Taylor leakage/bound ratio: `{out.worst_taylor_leakage_ratio:.9f}`
- worst shell gradient/bound ratio: `{out.worst_shell_gradient_ratio:.9f}`
- worst affine-curvature invariance residual: `{out.worst_affine_curvature_residual:.3e}`
- worst optimizer residual: `{out.worst_optimizer_residual:.3e}`
- extreme transformed condition number: `{out.extreme_condition_number:.3e}`
- extreme leakage/bound ratio: `{out.extreme_leakage_ratio:.9f}`

This closes the geometry of an ellipsoidal moving moat and the generic physical
convolution-filter commutator coefficient.  It does not yet insert those bounds
into the full localized SGS/pressure wave-packet identity; pressure/window, `RU`,
viscous-boundary and partition-overlap terms remain continuum terms.
"""
    (args.outdir/'summary.md').write_text(md);print(md)

if __name__=='__main__':main()
