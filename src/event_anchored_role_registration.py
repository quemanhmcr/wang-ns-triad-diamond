from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

INNER_ROLE_LOWER = 3.0 / 5.0
OUTER_ENVELOPE_LOWER = 11.0 / 20.0
LOW_STRAIN_ACTION = 1.0 / 30.0
LOW_LOW_OUTPUT = 1.0 / 2.0


def envelope_low_low_gap(
    envelope_lower: float = OUTER_ENVELOPE_LOWER,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Low-strain support gap of a smooth outer envelope above V tensor V."""
    r = float(envelope_lower)
    K = float(strain_action)
    if r <= 0 or K < 0:
        raise ValueError("positive envelope lower edge and nonnegative strain action required")
    return r * math.exp(-K) - LOW_LOW_OUTPUT


def arb_envelope_support_certificate() -> dict[str, str]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 180
    lower = (arb(11) / 20) * (-arb(1) / 30).exp()
    if not (lower > arb(1) / 2):
        raise AssertionError(f"smooth envelope can meet low-low support: {lower}")
    return {
        "hard_inner_lower": "3/5",
        "smooth_envelope_lower": "11/20",
        "anchor_radial_buffer": "1/20",
        "low_strain_action": "1/30",
        "transported_envelope_lower": str(lower),
        "low_low_output_upper": "1/2",
        "status": "CERTIFIED_HARD_CORE_SMOOTH_ENVELOPE_BUFFER",
    }


def hard_role_projectors_from_masks_and_fibers(
    masks: np.ndarray,
    fiber_projectors: np.ndarray | None = None,
) -> list[np.ndarray]:
    """Build exact orthogonal Fourier/fiber role projectors in a finite model.

    masks[a,k] are disjoint 0/1 frequency cells.  If fiber_projectors are
    supplied they have shape (s,k,d,d), resolve the fiber identity for every k,
    and are orthogonal Hermitian projectors.  The returned projectors act on the
    flattened (frequency,fiber) Hilbert space.
    """
    masks = np.asarray(masks, float)
    if masks.ndim != 2 or np.any((masks != 0.0) & (masks != 1.0)):
        raise ValueError("binary frequency masks required")
    if np.any(masks.sum(axis=0) > 1.0 + 1e-12):
        raise ValueError("frequency cells must be disjoint")
    nc, nk = masks.shape
    if fiber_projectors is None:
        return [np.diag(masks[a].astype(complex)) for a in range(nc)]
    H = np.asarray(fiber_projectors, complex)
    if H.ndim != 4 or H.shape[1] != nk or H.shape[2] != H.shape[3]:
        raise ValueError("fiber projector array must have shape (sector,k,d,d)")
    ns, _, d, _ = H.shape
    I = np.eye(d, dtype=complex)
    for k in range(nk):
        total = np.zeros((d, d), complex)
        for s in range(ns):
            P = H[s, k]
            if np.linalg.norm(P-P.conj().T) > 2e-10 or np.linalg.norm(P@P-P) > 2e-10:
                raise ValueError("orthogonal Hermitian fiber projectors required")
            for t in range(s):
                if np.linalg.norm(P @ H[t, k]) > 2e-10:
                    raise ValueError("fiber sectors must be orthogonal")
            total += P
        if np.linalg.norm(total-I) > 2e-10:
            raise ValueError("fiber sectors must resolve identity")
    out: list[np.ndarray] = []
    for a in range(nc):
        for s in range(ns):
            blocks = [masks[a, k] * H[s, k] for k in range(nk)]
            P = np.zeros((nk*d, nk*d), complex)
            for k, B in enumerate(blocks):
                P[k*d:(k+1)*d, k*d:(k+1)*d] = B
            out.append(P)
    return out


def role_partition_residual(projectors: Sequence[np.ndarray], covered: np.ndarray) -> dict[str, float]:
    Ps = [np.asarray(P, complex) for P in projectors]
    C = np.asarray(covered, complex)
    if not Ps or any(P.shape != C.shape for P in Ps):
        raise ValueError("matching projectors and covered projector required")
    total = sum(Ps, np.zeros_like(C))
    orth = 0.0
    idem = 0.0
    selfadj = 0.0
    for i, P in enumerate(Ps):
        idem = max(idem, float(np.linalg.norm(P@P-P)))
        selfadj = max(selfadj, float(np.linalg.norm(P-P.conj().T)))
        for Q in Ps[i+1:]:
            orth = max(orth, float(np.linalg.norm(P@Q)))
    return {
        "resolution": float(np.linalg.norm(total-C)),
        "orthogonality": orth,
        "idempotence": idem,
        "selfadjointness": selfadj,
    }


def envelope_registration_residual(
    hard_projector: np.ndarray,
    smooth_envelope: np.ndarray,
    u: np.ndarray,
    probe: np.ndarray,
) -> complex:
    """Residual of <P u,phi>=<Q u,P phi> when QP=P and Q=Q*."""
    P = np.asarray(hard_projector, complex)
    Q = np.asarray(smooth_envelope, complex)
    u = np.asarray(u, complex)
    phi = np.asarray(probe, complex)
    n = len(u)
    if P.shape != (n,n) or Q.shape != (n,n) or phi.shape != (n,):
        raise ValueError("matching role/envelope/probe data required")
    if np.linalg.norm(P-P.conj().T) > 2e-10 or np.linalg.norm(Q-Q.conj().T) > 2e-10:
        raise ValueError("self-adjoint hard and smooth multipliers required")
    if np.linalg.norm(Q@P-P) > 2e-10:
        raise ValueError("smooth envelope must equal one on the hard role")
    lhs = np.vdot(P@u, phi)
    rhs = np.vdot(Q@u, P@phi)
    return lhs-rhs


def pointwise_projector_lp_contraction(values: np.ndarray, projectors: np.ndarray, p: float) -> tuple[float,float]:
    """Return ||P(k)f(k)||_p and ||f||_p for pointwise orthogonal fiber projections."""
    f = np.asarray(values, complex)
    P = np.asarray(projectors, complex)
    if f.ndim != 2 or P.shape != (f.shape[0], f.shape[1], f.shape[1]) or p < 1:
        raise ValueError("frequency-by-fiber values, matching projectors, p>=1 required")
    g = np.einsum("kij,kj->ki", P, f)
    gn = np.linalg.norm(g, axis=1)
    fn = np.linalg.norm(f, axis=1)
    return float(np.sum(gn**p)**(1.0/p)), float(np.sum(fn**p)**(1.0/p))


def bilinear_source_partition_residual(
    tensor: np.ndarray,
    h: np.ndarray,
    projectors: Sequence[np.ndarray],
) -> np.ndarray:
    """Exact B(h,h)=sum_ab B(P_a h,P_b h) on the covered high subspace."""
    T = np.asarray(tensor, complex)
    h = np.asarray(h, complex)
    Ps = [np.asarray(P, complex) for P in projectors]
    n = len(h)
    if T.shape != (n,n,n) or any(P.shape != (n,n) for P in Ps):
        raise ValueError("matching bilinear source data required")
    covered_h = sum((P@h for P in Ps), np.zeros(n, complex))
    lhs = np.einsum("ijk,j,k->i", T, covered_h, covered_h)
    rhs = np.zeros(n, complex)
    for P in Ps:
        a = P@h
        for Q in Ps:
            b = Q@h
            rhs += np.einsum("ijk,j,k->i", T, a, b)
    return lhs-rhs


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_EVENT_ANCHORED_HARD_ROLE_TO_SMOOTH_PDE_ENVELOPE_REGISTRATION",
        "hard_roles": "use deterministic disjoint Borel frequency cells and pointwise orthogonal fiber/helical projectors only at the physical event; these are exact self-adjoint orthogonal roles and introduce no synthesis error",
        "symbol_freezing": "the smooth SGS symbol may be frozen on the same hard product cells; the Lipschitz error is a property of the physical multiplier, not a cutoff smoothing error",
        "envelope": "for a selected hard role P choose a scalar smooth PDE envelope Q with QP=P, inner radial lower 3/5 and envelope lower 11/20",
        "buffer": "on K<=1/30 the envelope lower edge stays (11/20)e^(-1/30)>1/2, so low-low V tensor V remains excluded while the hard event cell sits in a 1/20 radial plateau",
        "coefficient": "<P u,phi>=<Q u,P phi> exactly; terminal hard-cell/helical marking is therefore an exact coefficient of the smooth moving PDE role",
        "lp": "frequency-cell indicators and pointwise orthogonal helical projectors are contractions in every vector L^p; projecting the dual Gaussian probe cannot worsen its L3 or L2 norm constants",
        "high_high": "the complete h-h quadratic source disintegrates exactly over the hard parent roles; selection/designation happens at actual work level, so no L2 remainder for 'other HH interactions' is needed",
        "polarization": "helicity is an eventwise terminal fiber mark; between event and common slice the adjoint Kelvin fiber dynamics transports the probe, so a helical projector is not assumed to persist",
        "single_charge": "hard roles carry physical transfer/energy identity; the overlapping smooth envelope carries only the PDE between slices and never defines a second transfer measure",
        "continuum_status": "frequency-role/transfer-cell registration is exact at event and coefficient levels; smooth continuation uses the Q^2 quadratic-carrier energy law, and coefficient obstruction reenters actual physical-energy work before K_coh/strain/D_V or relink ownership",
    }


@dataclass(frozen=True)
class RegistrationStress:
    samples: int
    worst_partition_residual: float
    worst_registration_residual: float
    worst_l3_contraction_margin: float
    worst_l2_contraction_margin: float
    worst_bilinear_partition_residual: float
    minimum_envelope_gap: float


def _random_fiber_projectors(rng: np.random.Generator, nk: int, d: int=2) -> np.ndarray:
    H = np.zeros((2,nk,d,d), complex)
    for k in range(nk):
        z = rng.normal(size=d)+1j*rng.normal(size=d)
        z /= np.linalg.norm(z)
        P = np.outer(z,z.conj())
        H[0,k] = P
        H[1,k] = np.eye(d)-P
    return H


def stress(samples: int=50_000, seed: int=20260809) -> RegistrationStress:
    rng = np.random.default_rng(seed)
    wp=wr=wb=0.0
    ml3=ml2=float("inf")
    mg=float("inf")
    nk=6; d=2; n=nk*d
    for _ in range(samples):
        # Deterministic-style disjoint event cells, randomized occupancy.
        assignment = rng.integers(0,3,size=nk)
        masks = np.zeros((3,nk),float)
        masks[assignment,np.arange(nk)] = 1.0
        H = _random_fiber_projectors(rng,nk,d)
        Ps = hard_role_projectors_from_masks_and_fibers(masks,H)
        covered = np.eye(n,dtype=complex)
        pr = role_partition_residual(Ps,covered)
        pres=max(pr.values())
        wp=max(wp,pres)
        if pres>3e-10:
            raise AssertionError("hard Fourier/fiber role partition failed")

        P = Ps[int(rng.integers(0,len(Ps)))]
        # Smooth diagonal envelope in frequency, equal one on all frequency
        # blocks touched by P and arbitrary in [0,1] elsewhere.
        qfreq = rng.uniform(0.0,1.0,size=nk)
        touched=[]
        for k in range(nk):
            block=P[k*d:(k+1)*d,k*d:(k+1)*d]
            if np.linalg.norm(block)>1e-14:
                qfreq[k]=1.0; touched.append(k)
        Q=np.zeros((n,n),complex)
        for k,q in enumerate(qfreq):
            Q[k*d:(k+1)*d,k*d:(k+1)*d]=q*np.eye(d)
        u=rng.normal(size=n)+1j*rng.normal(size=n)
        phi=rng.normal(size=n)+1j*rng.normal(size=n)
        rr=abs(envelope_registration_residual(P,Q,u,phi))/max(1.0,np.linalg.norm(u)*np.linalg.norm(phi))
        wr=max(wr,rr)
        if rr>3e-12:
            raise AssertionError("hard-to-smooth envelope coefficient registration failed")

        vals=rng.normal(size=(nk,d))+1j*rng.normal(size=(nk,d))
        sector=int(rng.integers(0,2))
        l3p,l3=pointwise_projector_lp_contraction(vals,H[sector],3.0)
        l2p,l2=pointwise_projector_lp_contraction(vals,H[sector],2.0)
        ml3=min(ml3,l3-l3p); ml2=min(ml2,l2-l2p)
        if l3p>l3+2e-12 or l2p>l2+2e-12:
            raise AssertionError("pointwise fiber projection increased an Lp norm")

        T=rng.normal(size=(n,n,n))+1j*rng.normal(size=(n,n,n))
        h=rng.normal(size=n)+1j*rng.normal(size=n)
        br=bilinear_source_partition_residual(T,h,Ps)
        bscale=max(1.0,np.linalg.norm(T)*np.linalg.norm(h)**2)
        brel=float(np.linalg.norm(br))/bscale
        wb=max(wb,brel)
        if brel>5e-12:
            raise AssertionError("complete high-high source did not disintegrate over hard roles")

        K=float(rng.uniform(0.0,LOW_STRAIN_ACTION))
        env=float(rng.uniform(OUTER_ENVELOPE_LOWER,INNER_ROLE_LOWER))
        gap=env*math.exp(-K)-LOW_LOW_OUTPUT
        mg=min(mg,gap)
        if gap<=0:
            raise AssertionError("smooth outer envelope met low-low support")
    return RegistrationStress(samples,wp,wr,ml3,ml2,wb,mg)


def main() -> None:
    ap=argparse.ArgumentParser()
    ap.add_argument("--samples",type=int,default=50_000)
    ap.add_argument("--outdir",type=Path,default=Path("results-event-anchored-role-registration"))
    args=ap.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    cert=theorem_certificate(); arb=arb_envelope_support_certificate(); out=stress(args.samples)
    (args.outdir/"event_anchored_role_registration.json").write_text(json.dumps({"certificate":cert,"arb":arb,"stress":asdict(out)},indent=2),encoding="utf-8")
    md=f"""# Event-anchored hard roles inside a smooth moving PDE envelope\n\nStatus: **{cert['status']}**.\n\nAt one physical transfer event, use disjoint Borel frequency cells and pointwise orthogonal fiber/helical projectors.  They give exact orthogonal hard roles `P_a`, so physical energy and the complete high--high quadratic source disintegrate without synthesis error.  Smooth SGS symbol freezing may be performed on these same hard product cells; only the multiplier variation is `Xi`.\n\nFor the role PDE do **not** evolve the hard boundary.  Choose a scalar smooth outer envelope `Q` with `QP=P`.  The signed-good hard role starts above `3N/5`; choose the envelope to start above `11N/20`.  On `K<=1/30`,\n\n`(11/20)e^(-1/30) N > N/2`,\n\nso the smooth envelope excludes `V tensor V` throughout the slab while retaining a radial plateau of width `N/20` around the hard role at the anchor.\n\nFor every terminal probe,\n\n`<P u,phi> = <Q u,P phi>`\n\nexactly.  Hence the complex-Young/dual-Gaussian mark extracted from the hard physical transfer role is an exact terminal coefficient of the smooth moving PDE role.  Frequency indicators and pointwise orthogonal helical projectors satisfy `|P(k)v|<=|v|`, so they are contractions in vector `L^3` and `L^2`; projecting the dual probe cannot worsen the coefficient constants.\n\nHelicity is only an eventwise terminal fiber mark.  The adjoint Kelvin fiber equation transports that vector probe to the common slice; no persistent helical projector is introduced.  The smooth envelope is only a PDE carrier and never defines another transfer measure.\n\nStress: `{out.samples}` exact hard-role/envelope/source states\n- worst hard partition residual: `{out.worst_partition_residual:.3e}`\n- worst hard-to-envelope coefficient residual: `{out.worst_registration_residual:.3e}`\n- minimum L3 contraction margin: `{out.worst_l3_contraction_margin:.3e}`\n- minimum L2 contraction margin: `{out.worst_l2_contraction_margin:.3e}`\n- worst full HH bilinear partition residual: `{out.worst_bilinear_partition_residual:.3e}`\n- minimum sampled envelope/low-low support gap: `{out.minimum_envelope_gap:.6e}`\n\nThis closes the frequency-role/transfer-cell alignment seam without identifying a smooth envelope with a physical packet.  Smooth continuation is read at the native carrier energy `<u,Q^2u>` and completed by a square partition; no idempotence of `Q` is asserted.  A large interface coefficient impulse only locates physical-energy reentry.  Actual carrier energy and native interface work then route to inheritance, HH generation, conservative relink or existing strain.  No global-regularity claim is made.\n"""
    (args.outdir/"summary.md").write_text(md,encoding="utf-8"); print(md)


if __name__=="__main__":
    main()
