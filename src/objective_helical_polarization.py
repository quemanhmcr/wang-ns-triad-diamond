from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def split_2x2(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    M = np.asarray(M, dtype=float)
    return 0.5 * (M + M.T), 0.5 * (M - M.T)


def objective_transverse_generator(A_perp: np.ndarray) -> np.ndarray:
    """Generator in the objective transverse frame.

    If the transverse frame spin is chosen as Omega_frame=-skew(A_perp), then
    coefficient dynamics has generator -sym(A_perp).  Rigid rotation is thus a
    gauge, while symmetric strain is physical.
    """
    S, _ = split_2x2(A_perp)
    return -S


def circular_matrix_from_real_symmetric(S: np.ndarray) -> np.ndarray:
    """Represent a real symmetric 2x2 tensor in (h_+,h_-) circular basis."""
    S = np.asarray(S, dtype=float)
    hp = np.array([1.0, 1j]) / math.sqrt(2.0)
    hm = np.array([1.0, -1j]) / math.sqrt(2.0)
    H = np.column_stack([hp, hm])
    return np.conjugate(H).T @ S @ H


def objective_helical_generator(A_perp: np.ndarray, viscous_rate: float = 0.0) -> np.ndarray:
    """Exact helical generator after removing transverse frame rotation."""
    S, _ = split_2x2(A_perp)
    H = circular_matrix_from_real_symmetric(S)
    return -H - float(viscous_rate) * np.eye(2, dtype=complex)


def strain_coordinates(D: np.ndarray) -> tuple[float, float]:
    """For trace-free symmetric D, return delta,beta in [[delta,beta],[beta,-delta]]."""
    D = np.asarray(D, dtype=float)
    if np.linalg.norm(D - D.T) > 1e-12 or abs(float(np.trace(D))) > 1e-12:
        raise ValueError("D must be symmetric trace free")
    return float(D[0, 0]), float(D[0, 1])


def strain_commutator(D1: np.ndarray, D2: np.ndarray) -> np.ndarray:
    return np.asarray(D1, float) @ np.asarray(D2, float) - np.asarray(D2, float) @ np.asarray(D1, float)


def strain_area_commutator(delta1: float, beta1: float, delta2: float, beta2: float) -> np.ndarray:
    area = delta1 * beta2 - beta1 * delta2
    return 2.0 * area * np.array([[0.0, 1.0], [-1.0, 0.0]])


def coherence_second_magnus_bound(epsilon: float, dT: float) -> float:
    """Bound the norm of the second Magnus term under ||D(t)-D0||op<=eps*d.

    With ||D0||op=d,
      ||[D(t1),D(t2)]||op <= (4 eps + 2 eps^2)d^2.
    The ordered triangle and the 1/2 Magnus prefactor give
      ||Omega_2|| <= (eps + eps^2/2)(dT)^2.
    This is a bound only on the second Magnus generator, not on the full
    time-ordered exponential.
    """
    if epsilon < 0 or dT < 0:
        raise ValueError("nonnegative parameters required")
    return (epsilon + 0.5 * epsilon * epsilon) * dT * dT


def frozen_tracefree_propagator(delta: float, beta: float, t: float) -> np.ndarray:
    """Exact exp(-Dt) for D^2=(delta^2+beta^2)I."""
    D = np.array([[delta, beta], [beta, -delta]], dtype=float)
    d = math.hypot(delta, beta)
    if d == 0:
        return np.eye(2)
    return math.cosh(d * t) * np.eye(2) - math.sinh(d * t) / d * D


@dataclass(frozen=True)
class ObjectivePolarizationStress:
    samples: int
    worst_objective_generator_residual: float
    worst_circular_generator_residual: float
    worst_commutator_residual: float
    worst_frozen_propagator_residual: float
    coherence_second_magnus_at_repo_threshold: float


def stress(samples: int = 50_000, seed: int = 20260807) -> ObjectivePolarizationStress:
    rng = np.random.default_rng(seed)
    wobj = wcirc = wcomm = wprop = 0.0
    for _ in range(samples):
        A = rng.normal(size=(2, 2))
        S, W = split_2x2(A)
        # Choosing frame spin=-W makes -A-(-W)=-S.
        direct = -A + W
        wobj = max(wobj, float(np.linalg.norm(direct - objective_transverse_generator(A))))

        visc = float(rng.uniform(0.0, 2.0))
        direct_h = -circular_matrix_from_real_symmetric(S) - visc * np.eye(2)
        wcirc = max(wcirc, float(np.linalg.norm(direct_h - objective_helical_generator(A, visc))))

        d1, b1, d2, b2 = rng.normal(size=4)
        D1 = np.array([[d1, b1], [b1, -d1]])
        D2 = np.array([[d2, b2], [b2, -d2]])
        wcomm = max(wcomm, float(np.linalg.norm(strain_commutator(D1, D2) - strain_area_commutator(d1, b1, d2, b2))))

        t = float(rng.uniform(-0.6, 0.6))
        P = frozen_tracefree_propagator(d1, b1, t)
        # Analytic derivative identity follows from D^2=d^2 I; compare with eig exp.
        vals, vecs = np.linalg.eigh(D1)
        exact = vecs @ np.diag(np.exp(-vals * t)) @ vecs.T
        wprop = max(wprop, float(np.linalg.norm(P - exact)))

    return ObjectivePolarizationStress(
        samples,
        wobj,
        wcirc,
        wcomm,
        wprop,
        coherence_second_magnus_bound(1 / 20, 1 / 30),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-objective-helical"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    payload = out.__dict__ | {
        "objective_frame_law": "G_obj=-sym(A_perp)-nu|k|^2 I",
        "commutator_law": "[D1,D2]=2(delta1 beta2-beta1 delta2) J",
        "interpretation": "first order is helicity mixing; second-order noncommutativity is local geometric phase",
    }
    (args.outdir / "objective_helical_polarization.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Objective helical polarization dynamics

- exact objective transverse generator: `-sym(A_perp)` after frame spin cancels `skew(A_perp)`
- circular/helical trace-free strain is off diagonal: first-order deformation is helicity conversion, not Berry phase
- exact commutator: `[D1,D2]=2(delta1 beta2-beta1 delta2) J`
- repository coherence thresholds `eps=1/20`, `dT=1/30` give second-Magnus bound
  `{out.coherence_second_magnus_at_repo_threshold:.9e}`
- random checks: `{out.samples}`
- worst objective-generator residual: `{out.worst_objective_generator_residual:.3e}`
- worst circular-generator residual: `{out.worst_circular_generator_residual:.3e}`
- worst strain-area commutator residual: `{out.worst_commutator_residual:.3e}`
- worst frozen propagator residual: `{out.worst_frozen_propagator_residual:.3e}`

The second-Magnus number is not advertised as a bound for the full time-ordered
propagator.  It identifies the correct local geometric-phase mechanism: rotation
of polarization appears from noncommuting strain orientations, while first-order
symmetric strain is a helicity mixer.  Failure of strain-orientation coherence is
already routed to the objective-strain/source ledger.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
