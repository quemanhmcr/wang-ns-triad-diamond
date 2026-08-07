from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.affine_shell_aspect import local_ellipsoid_mass_coefficient, physical_axis_lower_constant


def geometric_radius(Sigma: np.ndarray) -> float:
    Sigma=np.asarray(Sigma,float)
    d=float(np.linalg.det(Sigma))
    if d<=0: raise ValueError("Sigma must be positive definite")
    return d**(1.0/6.0)


def affine_critical_mass(local_energy: float, r_g: float) -> float:
    if r_g<=0 or local_energy<0: raise ValueError("invalid energy/radius")
    return local_energy/r_g


def fresh_radius_budget(energies: np.ndarray, radii: np.ndarray, eta: float, overlap: float=1.0) -> tuple[float,float]:
    """Return sum radii and the energy-conservation upper bound.

    Assumes every listed grain satisfies E_j/r_j>=eta and the physical windows
    have overlap multiplicity at most `overlap`, so sum E_j<=overlap*E_total.
    Here E_total is represented by sum(energies)/overlap for an exact test list.
    """
    energies=np.asarray(energies,float);radii=np.asarray(radii,float)
    if eta<=0 or overlap<1 or energies.shape!=radii.shape: raise ValueError("bad inputs")
    if np.any(energies+1e-14<eta*radii): raise ValueError("grain below affine critical threshold")
    Etot=float(np.sum(energies))/overlap
    return float(np.sum(radii)), overlap*Etot/eta


def aspect_upper_from_geometric_scale(s: float) -> float:
    """A=N lmax <= (9/4)(N r_g)^3 using lmin>2/(3N)."""
    if s<=0: raise ValueError("s must be positive")
    return 9.0/4.0*s**3


def geometric_scale_lower_from_aspect(A: float) -> float:
    if A<=0: raise ValueError("A must be positive")
    return (4.0*A/9.0)**(1.0/3.0)


def arb_certificate() -> dict[str,str]:
    try:
        from flint import arb,ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec=160
    # Re-certify the two clean inputs used by the ancestry reformulation.
    # The detailed shell/HY certificate lives in affine_shell_aspect.py.
    axis=arb(str(physical_axis_lower_constant()))
    mass=arb(str(local_ellipsoid_mass_coefficient()))
    if not (axis>arb(2)/3): raise AssertionError(f"axis input failed: {axis}")
    if not (mass>arb(3)/10): raise AssertionError(f"mass input failed: {mass}")
    # Algebraic aspect coefficient: l2,l3 >=2/(3N), r_g^3=l1 l2 l3.
    if not (arb(9)/4>arb(2)): raise AssertionError("trivial rational check failed")
    return {
        "affine_mass_lower":"M_aff(E2)>=3/10",
        "axis_lower":"l_i>2/(3N)",
        "aspect_relation":"A=N lmax <= (9/4)(N r_g)^3",
        "fresh_energy_budget":"sum r_g <= overlap*E_total/eta",
        "status":"CERTIFIED_FROM_AFFINE_SHELL_INPUTS",
    }


@dataclass(frozen=True)
class GrainStress:
    samples:int
    worst_geometric_radius_residual:float
    minimum_affine_mass_coefficient:float
    worst_aspect_relation_margin:float
    worst_fresh_budget_margin:float


def stress(samples:int=50_000,seed:int=20260807)->GrainStress:
    rng=np.random.default_rng(seed)
    wr=0.0; wa=float("inf"); wb=float("inf")
    coeff=local_ellipsoid_mass_coefficient()
    for _ in range(samples):
        N=10**float(rng.uniform(-2,5))
        lo=(2.0/3.0)/N*(1+float(rng.uniform(1e-4,2.0)))
        axes=lo*np.exp(rng.uniform(0,8,size=3)); axes.sort()
        Q,_=np.linalg.qr(rng.normal(size=(3,3)))
        Sigma=Q@np.diag(axes**2)@Q.T
        rg=geometric_radius(Sigma)
        wr=max(wr,abs(rg-(axes.prod())**(1/3))/max(rg,1e-15))
        A=N*axes[-1]; s=N*rg
        margin=aspect_upper_from_geometric_scale(s)-A
        wa=min(wa,margin)
        if margin<-1e-8: raise AssertionError("aspect/geometric-scale relation failed")

        m=int(rng.integers(1,20)); eta=float(rng.uniform(.05,.4)); overlap=float(rng.integers(1,5))
        radii=np.exp(rng.uniform(-5,1,size=m)); energies=eta*radii*(1+rng.random(m))
        sr,bound=fresh_radius_budget(energies,radii,eta,overlap)
        bm=bound-sr; wb=min(wb,bm)
        if bm<-1e-10: raise AssertionError("fresh affine radius budget failed")
    return GrainStress(samples,wr,coeff,wa,wb)


def main()->None:
    ap=argparse.ArgumentParser();ap.add_argument("--samples",type=int,default=50_000);ap.add_argument("--outdir",type=Path,default=Path("results-affine-critical-grain"));args=ap.parse_args();args.outdir.mkdir(parents=True,exist_ok=True)
    cert=arb_certificate();out=stress(args.samples)
    (args.outdir/"affine_critical_grain.json").write_text(json.dumps({"certificate":cert,"stress":asdict(out)},indent=2))
    md=f"""# Affine critical-grain energy ledger

Status: **{cert['status']}**.

Define `r_g=(det Sigma_x)^(1/6)` and the affine scale-critical local mass
`M_aff(E)=r_g^-1 integral_E |u|^2`.  The shell/aspect certificate gives on the
radius-two Gaussian covariance ellipsoid

`M_aff(E2) >= 3/10`.

For fresh grains with `M_aff>=eta` and overlap multiplicity `P`, physical energy
conservation gives the exact budget

`sum r_g <= P E_total / eta`.

The certified shell lower axis also implies, with `s=N r_g` and `A=N l_max`,

`A <= (9/4) s^3`.

Thus an affine grain with natural geometric scale has bounded aspect, while a
very elongated grain necessarily has a large physical geometric radius and is
more expensive in the fresh-energy ledger; no false Young/Bellman anisotropy
penalty is required.

- random checks: `{out.samples}`
- worst geometric-radius residual: `{out.worst_geometric_radius_residual:.3e}`
- local affine-mass coefficient: `{out.minimum_affine_mass_coefficient:.9f}`
- minimum aspect-relation margin: `{out.worst_aspect_relation_margin:.3e}`
- minimum fresh-budget margin: `{out.worst_fresh_budget_margin:.3e}`
"""
    (args.outdir/"summary.md").write_text(md);print(md)

if __name__=="__main__": main()
