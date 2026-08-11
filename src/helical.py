from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

import numpy as np

Array = np.ndarray


def _wavevector(k: Array, name: str = "wavevector") -> Array:
    q = np.asarray(k, dtype=float)
    if q.shape != (3,) or np.any(~np.isfinite(q)):
        raise ValueError(f"{name} must be a finite three-vector")
    return q


def stable_norm3(k: Array) -> float:
    """Euclidean norm without squaring away tiny or overflowing large modes."""
    q = _wavevector(k)
    return float(math.hypot(*(float(x) for x in q)))


def _canonical_sign(k: Array) -> int:
    """Return +1 if the first nonzero component is positive, else -1."""
    for x in _wavevector(k):
        if float(x) != 0.0:
            return 1 if x > 0 else -1
    raise ValueError("zero wavevector")


def helical_basis(k: Array, s: int) -> Array:
    """A deterministic unit helical vector with h_s(-k)=conj(h_s(k))."""
    k = _wavevector(k)
    if s not in (-1, 1):
        raise ValueError("helicity must be ±1")
    norm = stable_norm3(k)
    if norm == 0:
        raise ValueError("zero wavevector")

    sign = _canonical_sign(k)
    if sign < 0:
        return np.conjugate(helical_basis(-k, s))

    khat = k / norm
    refs = np.eye(3)
    ref = refs[np.argmin(np.abs(refs @ khat))]
    e1 = np.cross(ref, khat)
    e1 /= stable_norm3(e1)
    e2 = np.cross(khat, e1)
    # i k × h = s |k| h
    h = (e1 + 1j * s * e2) / math.sqrt(2.0)
    return h


def check_helical_eigenvector(k: Array, s: int, atol: float = 1e-10) -> bool:
    h = helical_basis(k, s)
    lhs = 1j * np.cross(k, h)
    rhs = s * stable_norm3(k) * h
    return bool(np.allclose(lhs, rhs, atol=atol, rtol=atol))


def coupling_g(k: Array, p: Array, q: Array, sk: int, sp: int, sq: int) -> complex:
    """Waleffe geometric coefficient; requires k+p+q=0."""
    k = _wavevector(k, "k")
    p = _wavevector(p, "p")
    q = _wavevector(q, "q")
    scale = max(stable_norm3(k), stable_norm3(p), stable_norm3(q))
    if scale == 0.0:
        raise ValueError("zero wavevector")
    if stable_norm3(k + p + q) > 2e-12 * scale:
        raise ValueError("triad does not close")
    hk = helical_basis(k, sk)
    hp = helical_basis(p, sp)
    hq = helical_basis(q, sq)
    return -0.5 * np.dot(np.cross(np.conjugate(hp), np.conjugate(hq)), np.conjugate(hk))


def triangle_area(k: float, p: float, q: float) -> float:
    x = (k + p + q) * (-k + p + q) * (k - p + q) * (k + p - q)
    return 0.25 * math.sqrt(max(0.0, x))


def coupling_magnitude_closed(k: float, p: float, q: float, sk: int, sp: int, sq: int) -> float:
    area = triangle_area(k, p, q)
    return area * abs(sk * k + sp * p + sq * q) / (2.0 * math.sqrt(2.0) * k * p * q)


@dataclass(frozen=True)
class EdgeResult:
    efficiency: float
    raw_prefactor: float
    g_abs: float
    g_phase: float
    target_phase: float
    forward_ratio: float


def edge_metrics(x: Array, y: Array, z: Array, sx: int, sy: int, sz: int) -> EdgeResult:
    """Metrics for parents x,y feeding child z=x+y.

    The closed triad is (x,y,-z). The child-energy coefficient is
    A = sx|x| - sy|y|. The phase target is chosen so dE_child/dt > 0.
    """
    x = _wavevector(x, "x")
    y = _wavevector(y, "y")
    z = _wavevector(z, "z")
    nx, ny, nz = map(stable_norm3, (x, y, z))
    scale = max(nx, ny, nz)
    if scale == 0.0:
        return EdgeResult(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    if stable_norm3(x + y - z) > 2e-12 * scale:
        raise ValueError("z must equal x+y")
    if min(nx, ny, nz) == 0.0:
        return EdgeResult(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    g = coupling_g(x, y, -z, sx, sy, sz)
    a = sx * nx - sy * ny
    forward_ratio = nz / max(nx, ny)
    progress = max(0.0, math.log(forward_ratio))
    raw = abs(a) * abs(g) / nz
    efficiency = progress * raw
    target = 0.0 if a >= 0 else math.pi
    return EdgeResult(
        efficiency=float(efficiency),
        raw_prefactor=float(raw),
        g_abs=float(abs(g)),
        g_phase=float(np.angle(g)),
        target_phase=float(target),
        forward_ratio=float(forward_ratio),
    )


def wrap_angle(x: float) -> float:
    return float((x + math.pi) % (2.0 * math.pi) - math.pi)


def diamond_metrics(a: Array, b: Array, c: Array, signs: Tuple[int, int, int, int, int, int]) -> Dict[str, object]:
    sa, sb, sc, sm, sn, sd = signs
    m = a + b
    n = b + c
    d = a + b + c
    edges = {
        "ab_m": edge_metrics(a, b, m, sa, sb, sm),
        "mc_d": edge_metrics(m, c, d, sm, sc, sd),
        "bc_n": edge_metrics(b, c, n, sb, sc, sn),
        "an_d": edge_metrics(a, n, d, sa, sn, sd),
    }
    geom_holonomy = wrap_angle(
        edges["ab_m"].g_phase + edges["mc_d"].g_phase
        - edges["bc_n"].g_phase - edges["an_d"].g_phase
    )
    target_holonomy = wrap_angle(
        edges["ab_m"].target_phase + edges["mc_d"].target_phase
        - edges["bc_n"].target_phase - edges["an_d"].target_phase
    )
    frustration = abs(wrap_angle(geom_holonomy - target_holonomy))
    vals = np.array([e.efficiency for e in edges.values()], dtype=float)
    return {
        "edges": edges,
        "min_efficiency": float(np.min(vals)),
        "geom_mean_efficiency": float(np.prod(np.maximum(vals, 0.0)) ** 0.25),
        "mean_efficiency": float(np.mean(vals)),
        "phase_frustration": float(frustration),
        "geom_holonomy": float(geom_holonomy),
        "target_holonomy": float(target_holonomy),
        "vectors": {"a": a, "b": b, "c": c, "m": m, "n": n, "d": d},
    }


def all_signs() -> Iterable[Tuple[int, int, int, int, int, int]]:
    return itertools.product((-1, 1), repeat=6)
