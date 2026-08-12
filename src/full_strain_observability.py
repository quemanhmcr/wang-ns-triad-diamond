from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import numpy as np

RSTAR_LO = Fraction(61090410158, 100_000_000_000)
RSTAR_HI = Fraction(61090410160, 100_000_000_000)
_I2 = np.eye(2)
_I2.setflags(write=False)


def tracefree_2x2(M: np.ndarray) -> np.ndarray:
    M = np.asarray(M, dtype=float)
    tr = 0.5 * (float(M[0, 0]) + float(M[1, 1]))
    return np.array(
        [[float(M[0, 0]) - tr, float(M[0, 1])],
         [float(M[1, 0]), float(M[1, 1]) - tr]],
        dtype=float,
    )


def extremal_parent_directions(rstar: float) -> tuple[np.ndarray, np.ndarray]:
    c = 1.0 / (2.0 * rstar)
    s = math.sqrt(1.0 - c * c)
    return np.array([c, s, 0.0]), np.array([c, -s, 0.0])


def transverse_frame(khat: np.ndarray) -> np.ndarray:
    """For an in-plane khat, return columns (in-plane tangent, z)."""
    kh = np.asarray(khat, dtype=float)
    kh = kh / np.linalg.norm(kh)
    t = np.array([-kh[1], kh[0], 0.0])
    t /= np.linalg.norm(t)
    return np.column_stack([t, np.array([0.0, 0.0, 1.0])])


def restriction_tracefree_norm2(S: np.ndarray, E: np.ndarray) -> float:
    R = E.T @ np.asarray(S, dtype=float) @ E
    D = tracefree_2x2(R)
    return float(np.sum(D * D))


def strain_observability(S: np.ndarray, rstar: float = 0.610904101586766) -> tuple[float, float]:
    """Return Q and ||S||_F^2 for symmetric trace-free S.

    Q is the sum of trace-free restriction energies on the triad plane and the
    two parent polarization planes k_i^perp at the symmetric extremizer.
    """
    S = np.asarray(S, dtype=float)
    if np.linalg.norm(S - S.T) > 1e-10 or abs(float(np.trace(S))) > 1e-10:
        raise ValueError("S must be symmetric and trace free")
    k1, k2 = extremal_parent_directions(rstar)
    Eplane = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
    Q = restriction_tracefree_norm2(S, Eplane)
    Q += restriction_tracefree_norm2(S, transverse_frame(k1))
    Q += restriction_tracefree_norm2(S, transverse_frame(k2))
    return Q, float(np.sum(S * S))


def arb_observability_certificate() -> dict[str, str]:
    """Certify Q >= 13/20 ||S||_F^2 throughout the r* bracket.

    In coordinates
      S=[[a,b,x],[b,d,y],[x,y,-a-d]],
    put C=cos(phi)^2=1/(4 r_*^2).  Q-lambda||S||^2 is block diagonal in
    (a,d), b, x, y.  We certify positivity of the scalar blocks and the two
    Sylvester minors of the (a,d) block with Arb.
    """
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 160

    def aq(q: Fraction):
        return arb(q.numerator) / q.denominator

    r = aq(RSTAR_LO).union(aq(RSTAR_HI))
    C = 1 / (4 * r * r)
    lam = arb(13) / 20

    # Coefficients of Q-lambda*N, N=||S||_F^2.
    A = C * C - 4 * C + arb(9) / 2 - 2 * lam
    D = C * C + 2 * C + arb(3) / 2 - 2 * lam
    L = -2 * C * C + 2 * C + 3 - 2 * lam  # coefficient of a*d
    B = -4 * C * C + 4 * C + 2 - 2 * lam
    X = 4 * (1 - C) - 2 * lam
    Y = 4 * C - 2 * lam
    det_ad = A * D - (L / 2) * (L / 2)

    zero = arb(0)
    for name, val in {"A": A, "det_ad": det_ad, "B": B, "X": X, "Y": Y}.items():
        if not (val > zero):
            raise AssertionError(f"observability positivity failed for {name}: {val}")

    # The off-plane x mode gives the actual smallest generalized eigenvalue
    # 2(1-C), about 0.6602495, so 13/20 is conservative.
    dangerous = 2 * (1 - C)
    if not (dangerous > lam):
        raise AssertionError(f"dangerous shear mode failed: {dangerous}")

    return {
        "rstar_ball": str(r),
        "cos2_half_angle_ball": str(C),
        "certified_observability_lower": "13/20",
        "dangerous_offplane_shear_ratio_ball": str(dangerous),
        "ad_first_minor_ball": str(A),
        "ad_determinant_ball": str(det_ad),
        "b_margin_ball": str(B),
        "x_margin_ball": str(X),
        "y_margin_ball": str(Y),
        "combined_shape_helicity_lower": "559/2000",
        "status": "CERTIFIED",
    }


@dataclass(frozen=True)
class ObservabilityStress:
    samples: int
    worst_ratio: float
    worst_combined_margin: float


def stress(samples: int = 100_000, seed: int = 20260807) -> ObservabilityStress:
    rng = np.random.default_rng(seed)
    r = 0.610904101586766
    worst_ratio = float("inf")
    worst_combined = float("inf")
    for _ in range(samples):
        M = rng.normal(size=(3, 3))
        S = 0.5 * (M + M.T)
        S -= np.trace(S) / 3.0 * np.eye(3)
        Q, N = strain_observability(S, r)
        if N < 1e-18:
            continue
        ratio = Q / N
        worst_ratio = min(worst_ratio, ratio)

        k1, k2 = extremal_parent_directions(r)
        Eplane = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
        Dp2 = restriction_tracefree_norm2(S, Eplane)
        T1 = restriction_tracefree_norm2(S, transverse_frame(k1))
        T2 = restriction_tracefree_norm2(S, transverse_frame(k2))
        # Certified scalar-shape speed H >= 43/100 Dp2; helicity mixer
        # |zeta_i|^2=T_i/2.  Therefore H+(43/50)sum|zeta|^2
        # >= (43/100) Q >= 559/2000 ||S||^2.
        lhs_lower = (43.0 / 100.0) * Dp2 + (43.0 / 100.0) * (T1 + T2)
        margin = lhs_lower / N - 559.0 / 2000.0
        worst_combined = min(worst_combined, margin)
        if ratio < 13.0 / 20.0 - 2e-12 or margin < -2e-12:
            raise AssertionError("full strain observability violated")
    return ObservabilityStress(samples, worst_ratio, worst_combined)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=100_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-full-strain"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = arb_observability_certificate()
    out = stress(args.samples)
    (args.outdir / "full_strain_observability.json").write_text(json.dumps({"certificate": cert, "stress": out.__dict__}, indent=2))
    md = f"""# Full 3D strain observability from shape + helicity

Status: **{cert['status']}**.

- certified triad/polarization tomography: `Q >= 13/20 ||S||_F^2`
- dangerous off-plane shear generalized ratio: `{cert['dangerous_offplane_shear_ratio_ball']}`
- combined physical observable:
  `H_shape + (43/50)(|zeta_1|^2+|zeta_2|^2) >= 559/2000 ||S||_F^2`
- random trace-free strains: `{out.samples}`
- worst observed Q/||S||^2: `{out.worst_ratio:.9f}`
- minimum combined-observable margin: `{out.worst_combined_margin:.3e}`

Thus an extremal three-dimensional triad cannot hide a large incompressible
symmetric strain from both its scalar side-length geometry and the opposite-
helicity conversion of its two parent polarizations.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
