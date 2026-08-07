from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from src.helical import coupling_magnitude_closed

Array = np.ndarray


def unit(v: Array) -> Array:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    if n <= 0:
        raise ValueError("zero vector")
    return v / n


def helical_with_normal(k: Array, s: int, normal: Array) -> Array:
    """Helical eigenvector in the oriented-plane gauge set by ``normal``.

    The unit normal must be perpendicular to k.  With khat=k/|k| and
    t=normal x khat,

        h_s=(t+i s normal)/sqrt(2),
        i k x h_s=s |k| h_s.

    This gauge is SO(3)-covariant when k and normal are rotated together.
    """
    if s not in (-1, 1):
        raise ValueError("helicity must be +/-1")
    kh = unit(k)
    n = unit(normal)
    if abs(float(np.dot(kh, n))) > 2e-10:
        raise ValueError("normal must be perpendicular to k")
    t = np.cross(n, kh)
    t /= np.linalg.norm(t)
    return (t + 1j * s * n) / math.sqrt(2.0)


def triad_normal(x: Array, y: Array) -> Array:
    n = np.cross(np.asarray(x, float), np.asarray(y, float))
    nn = float(np.linalg.norm(n))
    if nn <= 1e-13 * max(1.0, np.linalg.norm(x) * np.linalg.norm(y)):
        raise ValueError("degenerate triad plane")
    return n / nn


def forward_normal_coupling(x: Array, y: Array, z: Array, sx: int, sy: int, sz: int) -> complex:
    """Waleffe coupling in a reality-compatible triad-normal gauge.

    Here z=x+y is the physical child, while the closed triad is (x,y,-z).
    Parents use h_s(.,n) and the closed child uses conj(h_s(z,n)); this
    respects the physical reality identification between z and -z.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    if np.linalg.norm(x + y - z) > 1e-9 * max(1.0, np.linalg.norm(z)):
        raise ValueError("z must equal x+y")
    n = triad_normal(x, y)
    hx = helical_with_normal(x, sx, n)
    hy = helical_with_normal(y, sy, n)
    hz = helical_with_normal(z, sz, n)
    hminus_z = np.conjugate(hz)
    return -0.5 * np.dot(np.cross(np.conjugate(hx), np.conjugate(hy)), np.conjugate(hminus_z))


def coupling_phase_is_quadrature(g: complex, atol: float = 1e-11) -> bool:
    return abs(float(np.real(g))) <= atol * max(1.0, abs(g))


def signed_normal_angle(k: Array, n1: Array, n2: Array) -> float:
    """Signed angle rotating n1 to n2 about +k, in (-pi,pi]."""
    kh = unit(k)
    a = unit(n1)
    b = unit(n2)
    if max(abs(float(np.dot(kh, a))), abs(float(np.dot(kh, b)))) > 2e-10:
        raise ValueError("normals must be perpendicular to k")
    c = float(np.clip(np.dot(a, b), -1.0, 1.0))
    ss = float(np.dot(kh, np.cross(a, b)))
    return math.atan2(ss, c)


def normal_transition_phase(k: Array, s: int, n1: Array, n2: Array) -> float:
    """Exact spin-1 transition phase between two triad-normal gauges.

    h_s(k;n2)=exp(-i s psi) h_s(k;n1), where psi is the signed dihedral
    rotation of n1 to n2 about k.
    """
    return -float(s) * signed_normal_angle(k, n1, n2)


def wrap_angle(x: float) -> float:
    return float((x + math.pi) % (2.0 * math.pi) - math.pi)


def spherical_helical(theta: float, phi: float, s: int) -> Array:
    """Standard local spherical gauge h=(e_theta+i s e_phi)/sqrt(2)."""
    if s not in (-1, 1):
        raise ValueError("helicity must be +/-1")
    st, ct = math.sin(theta), math.cos(theta)
    sp, cp = math.sin(phi), math.cos(phi)
    e_theta = np.array([ct * cp, ct * sp, -st], dtype=float)
    e_phi = np.array([-sp, cp, 0.0], dtype=float)
    return (e_theta + 1j * s * e_phi) / math.sqrt(2.0)


def berry_connection_phi(theta: float, s: int) -> float:
    """Coefficient A_phi in A=i<h,dh>=s cos(theta) dphi."""
    if s not in (-1, 1):
        raise ValueError("helicity must be +/-1")
    return float(s) * math.cos(theta)


def berry_curvature_theta_phi(theta: float, s: int) -> float:
    """Coefficient F_{theta phi}=-s sin(theta)."""
    if s not in (-1, 1):
        raise ValueError("helicity must be +/-1")
    return -float(s) * math.sin(theta)


def berry_chern_number(s: int) -> int:
    """(2pi)^-1 int_{S2} F=-2s: the spin-1 helicity line has Chern number +/-2."""
    if s not in (-1, 1):
        raise ValueError("helicity must be +/-1")
    return -2 * s


def rotation_matrix(axis: Array, angle: float) -> Array:
    a = unit(axis)
    x, y, z = a
    K = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + math.sin(angle) * K + (1.0 - math.cos(angle)) * (K @ K)


def transverse_strain_helicity_matrix(delta: float, beta: float) -> Array:
    """Trace-free real transverse strain in the circular/helical basis.

    D=[[delta,beta],[beta,-delta]].  In (h_+,h_-) order its representation is
    [[0, delta-i beta],[delta+i beta,0]].
    """
    z = complex(delta, -beta)
    return np.array([[0.0 + 0.0j, z], [np.conjugate(z), 0.0 + 0.0j]])


def frozen_pure_strain_conversion(d: float, t: float) -> tuple[float, float, float]:
    """Exact helicity conversion for D=diag(d,-d), k fixed normal to the plane.

    Starting from a unit + helical amplitude, exp(-Dt) gives circular-basis
    amplitudes (cosh(dt), -sinh(dt)).  Returns plus amplitude magnitude,
    minus amplitude magnitude, and opposite-helicity energy fraction.
    """
    x = abs(float(d) * float(t))
    cp = math.cosh(x)
    sm = math.sinh(x)
    frac = sm * sm / (cp * cp + sm * sm)
    return cp, sm, frac


@dataclass(frozen=True)
class SpinTransportStress:
    samples: int
    worst_eigenvector_residual: float
    worst_rotation_covariance_residual: float
    worst_transition_residual: float
    worst_coupling_real_fraction: float
    worst_coupling_magnitude_residual: float
    worst_berry_connection_residual: float
    worst_helicity_matrix_residual: float


def stress(samples: int = 50_000, seed: int = 20260807) -> SpinTransportStress:
    rng = np.random.default_rng(seed)
    worst_eig = worst_rot = worst_trans = worst_real = worst_mag = worst_berry = worst_mix = 0.0
    for _ in range(samples):
        # random nondegenerate forward triad
        x = rng.normal(size=3)
        y = rng.normal(size=3)
        if np.linalg.norm(np.cross(x, y)) < 0.1 * np.linalg.norm(x) * np.linalg.norm(y):
            y = y + np.roll(x, 1)
        z = x + y
        if np.linalg.norm(z) < 1e-3:
            z = z + np.array([0.1, -0.2, 0.3])
            y = z - x
        n = triad_normal(x, y)
        s = int(rng.choice([-1, 1]))
        h = helical_with_normal(x, s, n)
        eig = 1j * np.cross(x, h) - s * np.linalg.norm(x) * h
        worst_eig = max(worst_eig, float(np.linalg.norm(eig)) / max(1.0, np.linalg.norm(x)))

        # SO(3) covariance of the triad-adapted gauge
        R = rotation_matrix(rng.normal(size=3), float(rng.uniform(-math.pi, math.pi)))
        hR = helical_with_normal(R @ x, s, R @ n)
        worst_rot = max(worst_rot, float(np.linalg.norm(hR - R @ h)))

        # Spin-1 normal transition law around the same carrier.
        psi = float(rng.uniform(-math.pi, math.pi))
        Rk = rotation_matrix(x, psi)
        n2 = Rk @ n
        h2 = helical_with_normal(x, s, n2)
        pred = np.exp(-1j * s * psi) * h
        worst_trans = max(worst_trans, float(np.linalg.norm(h2 - pred)))

        # Triad-normal coupling is in exact quadrature and has the closed magnitude.
        sx, sy, sz = (int(q) for q in rng.choice([-1, 1], size=3))
        g = forward_normal_coupling(x, y, z, sx, sy, sz)
        worst_real = max(worst_real, abs(float(np.real(g))) / max(1e-15, abs(g)))
        closed = coupling_magnitude_closed(np.linalg.norm(x), np.linalg.norm(y), np.linalg.norm(z), sx, sy, sz)
        worst_mag = max(worst_mag, abs(abs(g) - closed) / max(1.0, closed))

        # Finite-difference check of A_phi=s cos(theta).
        theta = float(rng.uniform(0.15, math.pi - 0.15))
        phi = float(rng.uniform(-math.pi, math.pi))
        dphi = 1e-7
        hs = spherical_helical(theta, phi, s)
        hs2 = spherical_helical(theta, phi + dphi, s)
        overlap_phase = float(np.angle(np.vdot(hs, hs2))) / dphi
        # arg <h(phi),h(phi+dphi)> = -A_phi dphi + o(dphi)
        worst_berry = max(worst_berry, abs(overlap_phase + berry_connection_phi(theta, s)))

        # Circular representation of trace-free transverse strain.
        delta, beta = rng.normal(size=2)
        D = np.array([[delta, beta], [beta, -delta]], dtype=float)
        hp = np.array([1.0, 1j]) / math.sqrt(2.0)
        hm = np.array([1.0, -1j]) / math.sqrt(2.0)
        H = np.column_stack([hp, hm])
        direct = np.conjugate(H).T @ D @ H
        formula = transverse_strain_helicity_matrix(delta, beta)
        worst_mix = max(worst_mix, float(np.linalg.norm(direct - formula)))

    return SpinTransportStress(samples, worst_eig, worst_rot, worst_trans, worst_real, worst_mag, worst_berry, worst_mix)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-helical-spin"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    payload = out.__dict__ | {
        "berry_curvature": "F_s=-s sin(theta) dtheta wedge dphi",
        "chern_numbers": {"s=+1": berry_chern_number(1), "s=-1": berry_chern_number(-1)},
        "single_triad_phase_statement": "triad-normal coupling is purely imaginary, so its phase is locally constant away from zeros",
        "frozen_pure_strain_conversion": "u_+=cosh(dt), u_-=-sinh(dt)",
    }
    (args.outdir / "helical_spin_transport.json").write_text(json.dumps(payload, indent=2))
    md = f"""# Helical spin / triad-normal transport

- exact local Berry curvature: `F_s=-s sin(theta) dtheta wedge dphi`
- exact Chern numbers: `c1(s=+1)=-2`, `c1(s=-1)=+2`
- single-triad normal gauge: coupling phase is quadrature (`+/- pi/2`) away from coupling zeros
- random 3D checks: `{out.samples}`
- worst SO(3)-covariance residual: `{out.worst_rotation_covariance_residual:.3e}`
- worst spin-1 normal-transition residual: `{out.worst_transition_residual:.3e}`
- worst coupling real/absolute ratio: `{out.worst_coupling_real_fraction:.3e}`
- worst coupling-magnitude residual: `{out.worst_coupling_magnitude_residual:.3e}`
- worst Berry-connection finite-difference residual: `{out.worst_berry_connection_residual:.3e}`
- worst transverse-strain/helicity-matrix residual: `{out.worst_helicity_matrix_residual:.3e}`

A single triad's absolute Berry phase is not charged as a cascade defect: the
triad's moving normal supplies an SO(3)-covariant gauge in which rigid rotation is
phase-free and the Waleffe coupling phase is constant.  The physical geometric
obstruction appears only when one Fourier mode is reused by incident triads with
different normals; their transition function is the spin-1 dihedral phase
`exp(-i s psi)`.
"""
    (args.outdir / "summary.md").write_text(md)
    print(md)


if __name__ == "__main__":
    main()
