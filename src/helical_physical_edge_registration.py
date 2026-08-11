from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from src.helical import coupling_g, edge_metrics, helical_basis
from src.single_edge_certificate import float_jstar

STATUS = (
    "EXACT_HELICAL_PHYSICAL_EDGE_REGISTRATION__DIRECT_LERAY_CURL_EQUALS_WALEFFE__"
    "UNORDERED_PARENT_PAIR_FACTOR_FOUR__NATIVE_MODAL_CAPACITY__"
    "SIGNED_UPPER_PROGRESS_EQUALS_A_J_C"
)


def _vec3(x: np.ndarray, name: str) -> np.ndarray:
    v = np.asarray(x, dtype=float)
    if v.shape != (3,) or np.any(~np.isfinite(v)):
        raise ValueError(f"{name} must be a finite three-vector")
    return v


def _complex_scalar(x: complex, name: str) -> complex:
    z = complex(x)
    if not (math.isfinite(z.real) and math.isfinite(z.imag)):
        raise ValueError(f"{name} must be finite")
    return z


def _helicity(s: int, name: str) -> int:
    q = int(s)
    if q not in (-1, 1):
        raise ValueError(f"{name} must be ±1")
    return q


def leray_project(k: np.ndarray, value: np.ndarray) -> np.ndarray:
    """Complex Fourier Leray projection at one nonzero real wavevector."""
    q = _vec3(k, "wavevector")
    v = np.asarray(value, dtype=complex)
    if v.shape != (3,) or np.any(~np.isfinite(v.real)) or np.any(~np.isfinite(v.imag)):
        raise ValueError("finite complex three-vector required")
    n2 = float(np.dot(q, q))
    if n2 <= 0.0:
        raise ValueError("nonzero wavevector required")
    return v - q * (np.dot(q, v) / n2)


def unordered_parent_curl_source(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    sx: int,
    sy: int,
    ax: complex,
    ay: complex,
) -> np.ndarray:
    """Actual unordered two-parent NS source at child ``z=x+y``.

    For incompressible Navier--Stokes,

        -P(u·grad u) = P(u x omega).

    With one helical Fourier mode on each parent, the unordered convolution orbit
    contains the two ordered terms ``u_x x omega_y`` and ``u_y x omega_x``.
    This function returns their exact Leray-projected child source.  No Young norm,
    packet count or causal normalization enters.
    """
    x = _vec3(x, "x")
    y = _vec3(y, "y")
    z = _vec3(z, "z")
    sx = _helicity(sx, "sx")
    sy = _helicity(sy, "sy")
    ax = _complex_scalar(ax, "ax")
    ay = _complex_scalar(ay, "ay")
    scale = max(1.0, float(np.linalg.norm(x)), float(np.linalg.norm(y)), float(np.linalg.norm(z)))
    if np.linalg.norm(x + y - z) > 2e-12 * scale:
        raise ValueError("physical parent pair must satisfy z=x+y")
    nx = float(np.linalg.norm(x))
    ny = float(np.linalg.norm(y))
    nz = float(np.linalg.norm(z))
    if min(nx, ny, nz) <= 1e-14:
        raise ValueError("nonzero parent and child wavevectors required")

    hx = helical_basis(x, sx)
    hy = helical_basis(y, sy)
    ux = ax * hx
    uy = ay * hy
    wx = sx * nx * ux
    wy = sy * ny * uy
    raw = np.cross(ux, wy) + np.cross(uy, wx)
    return leray_project(z, raw)


def waleffe_child_source_coefficient(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    sx: int,
    sy: int,
    sz: int,
    ax: complex,
    ay: complex,
) -> complex:
    """Child helical coefficient of the same unordered NS source.

    The repository convention is ``g=coupling_g(x,y,-z,sx,sy,sz)``.  Reality of
    the helical frame gives

        conj(h_z)·(h_x x h_y) = -2 conj(g),

    while ``u_x x omega_y + u_y x omega_x`` contributes the signed factor
    ``sy|y|-sx|x|``.  Therefore

        <h_z,F_z> = 2 (sx|x|-sy|y|) conj(g) a_x a_y.
    """
    x = _vec3(x, "x")
    y = _vec3(y, "y")
    z = _vec3(z, "z")
    sx = _helicity(sx, "sx")
    sy = _helicity(sy, "sy")
    sz = _helicity(sz, "sz")
    ax = _complex_scalar(ax, "ax")
    ay = _complex_scalar(ay, "ay")
    scale = max(1.0, float(np.linalg.norm(x)), float(np.linalg.norm(y)), float(np.linalg.norm(z)))
    if np.linalg.norm(x + y - z) > 2e-12 * scale:
        raise ValueError("physical parent pair must satisfy z=x+y")
    nx = float(np.linalg.norm(x))
    ny = float(np.linalg.norm(y))
    if min(nx, ny, float(np.linalg.norm(z))) <= 1e-14:
        raise ValueError("nonzero parent and child wavevectors required")
    g = coupling_g(x, y, -z, sx, sy, sz)
    signed_frequency = sx * nx - sy * ny
    return 2.0 * signed_frequency * np.conjugate(g) * ax * ay


def direct_child_source_coefficient(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    sx: int,
    sy: int,
    sz: int,
    ax: complex,
    ay: complex,
) -> complex:
    source = unordered_parent_curl_source(x, y, z, sx, sy, ax, ay)
    hz = helical_basis(_vec3(z, "z"), _helicity(sz, "sz"))
    return complex(np.vdot(hz, source))


def phase_alignment(
    signed_frequency: float,
    coupling: complex,
    ax: complex,
    ay: complex,
    az: complex,
) -> float:
    """Signed physical phase/orientation efficiency ``c_e in [-1,1]``.

    The sign of ``sx|x|-sy|y|`` is geometry, while the remaining real part is the
    phase alignment of the actual modal amplitudes with the Waleffe coefficient.
    Their product is invariant under parent swap and under reciprocal helical-basis
    phase gauge.
    """
    a = float(signed_frequency)
    g = _complex_scalar(coupling, "coupling")
    ax = _complex_scalar(ax, "ax")
    ay = _complex_scalar(ay, "ay")
    az = _complex_scalar(az, "az")
    den = abs(g) * abs(ax) * abs(ay) * abs(az)
    if a == 0.0 or den == 0.0:
        return 0.0
    raw = math.copysign(1.0, a) * float(
        np.real(np.conjugate(az) * np.conjugate(g) * ax * ay)
    ) / den
    if abs(raw) > 1.0 + 5e-12:
        raise AssertionError("modal phase alignment left [-1,1]")
    return float(max(-1.0, min(1.0, raw)))


@dataclass(frozen=True)
class HelicalPhysicalEdgeRegistration:
    parent_x_frequency: float
    parent_y_frequency: float
    child_frequency: float
    parent_top_frequency: float
    forward_ratio: float
    scale_progress: float
    signed_frequency_factor: float
    coupling_abs: float
    geometric_multiplier_J: float
    global_multiplier_Jstar: float
    normalized_multiplier: float
    phase_alignment: float
    native_modal_capacity: float
    signed_child_energy_work: float
    signed_upper_progress_work: float
    registered_upper_progress_work: float
    direct_child_source_coefficient: complex
    waleffe_child_source_coefficient: complex
    direct_vs_waleffe_residual: float
    leray_pairing_residual: float
    upper_progress_identity_residual: float
    positive_forward_work: bool
    unordered_parent_orientation: bool = True
    young_norm_used_as_capacity: bool = False
    duhamel_weight_used_as_causal_law: bool = False

    def __post_init__(self) -> None:
        positive = (
            self.parent_x_frequency,
            self.parent_y_frequency,
            self.child_frequency,
            self.parent_top_frequency,
            self.global_multiplier_Jstar,
        )
        if not all(math.isfinite(x) and x > 0 for x in positive):
            raise ValueError("positive finite physical frequencies/Jstar required")
        finite = (
            self.forward_ratio,
            self.scale_progress,
            self.signed_frequency_factor,
            self.coupling_abs,
            self.geometric_multiplier_J,
            self.normalized_multiplier,
            self.phase_alignment,
            self.native_modal_capacity,
            self.signed_child_energy_work,
            self.signed_upper_progress_work,
            self.registered_upper_progress_work,
            self.direct_vs_waleffe_residual,
            self.leray_pairing_residual,
            self.upper_progress_identity_residual,
        )
        if not all(math.isfinite(x) for x in finite):
            raise ValueError("finite physical edge registration data required")
        if self.native_modal_capacity < 0 or self.geometric_multiplier_J < 0 or self.coupling_abs < 0:
            raise ValueError("capacity and geometric magnitudes must be nonnegative")
        if self.normalized_multiplier < -5e-12 or self.normalized_multiplier > 1.0 + 5e-10:
            raise AssertionError("single-edge multiplier exceeds the certified global Jstar envelope")
        if abs(self.phase_alignment) > 1.0 + 5e-12:
            raise AssertionError("phase alignment must lie in [-1,1]")
        if not self.unordered_parent_orientation or self.young_norm_used_as_capacity or self.duhamel_weight_used_as_causal_law:
            raise ValueError("physical helical edge registration used a forbidden observer/causal replacement")


def register_helical_physical_edge(
    *,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    sx: int,
    sy: int,
    sz: int,
    ax: complex,
    ay: complex,
    az: complex,
) -> HelicalPhysicalEdgeRegistration:
    """Register one actual unordered helical parent pair on its native edge law."""
    x = _vec3(x, "x")
    y = _vec3(y, "y")
    z = _vec3(z, "z")
    sx = _helicity(sx, "sx")
    sy = _helicity(sy, "sy")
    sz = _helicity(sz, "sz")
    ax = _complex_scalar(ax, "ax")
    ay = _complex_scalar(ay, "ay")
    az = _complex_scalar(az, "az")
    scale = max(1.0, float(np.linalg.norm(x)), float(np.linalg.norm(y)), float(np.linalg.norm(z)))
    if np.linalg.norm(x + y - z) > 2e-12 * scale:
        raise ValueError("physical parent pair must satisfy z=x+y")
    nx, ny, nz = map(float, map(np.linalg.norm, (x, y, z)))
    if min(nx, ny, nz) <= 1e-14:
        raise ValueError("nonzero parent and child wavevectors required")

    direct = direct_child_source_coefficient(x, y, z, sx, sy, sz, ax, ay)
    waleffe = waleffe_child_source_coefficient(x, y, z, sx, sy, sz, ax, ay)
    coeff_scale = max(1.0, abs(direct), abs(waleffe))
    coeff_res = abs(direct - waleffe)
    if coeff_res > 2e-10 * coeff_scale:
        raise AssertionError("direct Leray/curl source disagrees with Waleffe helical coefficient")

    # Pairing a divergence-free child helical vector makes Leray projection free.
    hx = helical_basis(x, sx)
    hy = helical_basis(y, sy)
    ux = ax * hx
    uy = ay * hy
    raw = np.cross(ux, sy * ny * uy) + np.cross(uy, sx * nx * ux)
    hz = helical_basis(z, sz)
    unprojected = complex(np.vdot(hz, raw))
    leray_res = abs(unprojected - direct)
    if leray_res > 2e-10 * max(1.0, abs(unprojected), abs(direct)):
        raise AssertionError("Leray projection changed a divergence-free child helical pairing")

    g = coupling_g(x, y, -z, sx, sy, sz)
    signed_frequency = sx * nx - sy * ny
    work = 2.0 * float(np.real(np.conjugate(az) * direct))
    parent_top = max(nx, ny)
    forward_ratio = nz / parent_top
    progress = max(0.0, math.log(forward_ratio))
    upper = work * progress

    # This is the native modal interaction capacity of the same physical edge.
    # Factor 2: the two parent orders in the unordered convolution orbit.
    # Factor 2: physical child energy derivative 2 Re(conj(a_z) F_z).
    capacity = 4.0 * nz * abs(ax) * abs(ay) * abs(az)
    metrics = edge_metrics(x, y, z, sx, sy, sz)
    J = float(metrics.efficiency)
    jstar = float_jstar()
    c = phase_alignment(signed_frequency, g, ax, ay, az)
    registered = capacity * J * c
    identity_res = upper - registered
    tol = 5e-10 * max(1.0, abs(upper), abs(registered), capacity * jstar)
    if abs(identity_res) > tol:
        raise AssertionError("physical upper-progress work failed the exact A*J*c registration identity")
    multiplier = J / jstar
    if multiplier > 1.0 + 5e-10:
        raise AssertionError("physical edge geometric multiplier exceeded global Jstar")

    return HelicalPhysicalEdgeRegistration(
        parent_x_frequency=nx,
        parent_y_frequency=ny,
        child_frequency=nz,
        parent_top_frequency=parent_top,
        forward_ratio=forward_ratio,
        scale_progress=progress,
        signed_frequency_factor=signed_frequency,
        coupling_abs=float(abs(g)),
        geometric_multiplier_J=J,
        global_multiplier_Jstar=jstar,
        normalized_multiplier=multiplier,
        phase_alignment=c,
        native_modal_capacity=capacity,
        signed_child_energy_work=work,
        signed_upper_progress_work=upper,
        registered_upper_progress_work=registered,
        direct_child_source_coefficient=direct,
        waleffe_child_source_coefficient=waleffe,
        direct_vs_waleffe_residual=coeff_res,
        leray_pairing_residual=leray_res,
        upper_progress_identity_residual=identity_res,
        positive_forward_work=(progress > 0.0 and work > 0.0),
    )


def gauge_transform_modal_data(
    coupling: complex,
    ax: complex,
    ay: complex,
    az: complex,
    theta_x: float,
    theta_y: float,
    theta_z: float,
) -> tuple[complex, complex, complex, complex]:
    """Reciprocal helical-basis phase gauge for the closed triad ``(x,y,-z)``."""
    g = _complex_scalar(coupling, "coupling")
    angles = tuple(float(t) for t in (theta_x, theta_y, theta_z))
    if not all(math.isfinite(t) for t in angles):
        raise ValueError("finite gauge angles required")
    tx, ty, tz = angles
    gp = np.exp(1j * (tz - tx - ty)) * g
    return (
        complex(gp),
        complex(np.exp(-1j * tx) * ax),
        complex(np.exp(-1j * ty) * ay),
        complex(np.exp(-1j * tz) * az),
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "ns_identity": "for incompressible NS, -P(u.grad u)=P(u x omega); one unordered parent orbit contributes u_x x omega_y + u_y x omega_x at z=x+y",
        "helical_projection": "with repository coupling g(x,y,-z), <h_z,F_z>=2(s_x|x|-s_y|y|)conj(g)a_xa_y; direct curl/Leray pairing is checked against this coefficient",
        "physical_work": "production child-energy convention is T_e=2 Re(conj(a_z)<h_z,F_z>), hence the unordered parent orbit carries the canonical factor four",
        "native_capacity": "A_e=4|z||a_x a_y a_z| is modal interaction capacity of the same physical Fourier edge; no Young norm product is substituted for it",
        "geometry": "J_e=log_+(|z|/p_top)|s_x|x|-s_y|y|| |g_e|/|z| and m_e=J_e/J_* with the global sign-exhausted single-edge J_*",
        "phase": "c_e=sign(s_x|x|-s_y|y|) Re(conj(a_z)conj(g)a_xa_y)/(|g a_x a_y a_z|) lies in [-1,1]",
        "registration": "exact signed identity T_e log_+(|z|/p_top)=A_e J_e c_e on the same physical event; backscatter is retained before positive Hahn restriction",
        "orientation": "swapping the two parent orientations changes both the signed frequency factor and Waleffe g sign, leaving source, work, capacity, J and c invariant",
        "gauge": "reciprocal phase changes of helical bases and modal amplitudes leave c_e and the physical edge identity invariant",
        "scale": "uniform wavevector dilation leaves J_e, m_e, c_e and forward ratio invariant while both T_e and A_e scale linearly",
        "scope": "this is one physical Fourier/helical edge registration theorem only; it does not yet construct the continuum edge measure, prove a low-deficit block, or identify every generic HH event with signed-good generation",
    }


@dataclass(frozen=True)
class HelicalPhysicalEdgeStress:
    samples: int
    worst_direct_waleffe_relative_residual: float
    worst_leray_pairing_relative_residual: float
    worst_upper_progress_relative_residual: float
    worst_parent_swap_relative_residual: float
    worst_gauge_phase_alignment_residual: float
    worst_scale_invariance_residual: float
    maximum_normalized_multiplier: float
    minimum_positive_forward_phase_alignment: float
    positive_forward_samples: int
    nonforward_samples: int


def _relative(a: complex, b: complex) -> float:
    return float(abs(a - b) / max(1.0, abs(a), abs(b)))


def stress(samples: int = 50_000, seed: int = 20260811) -> HelicalPhysicalEdgeStress:
    rng = np.random.default_rng(seed)
    wcoeff = wleray = wident = wswap = wgauge = wscale = 0.0
    maxm = 0.0
    minpos = math.inf
    npos = nnon = 0

    for _ in range(samples):
        # Generic nondegenerate physical parent vectors; z is the actual child.
        while True:
            x = rng.normal(size=3)
            y = rng.normal(size=3)
            z = x + y
            if min(np.linalg.norm(x), np.linalg.norm(y), np.linalg.norm(z)) > 0.15 and np.linalg.norm(np.cross(x, y)) > 0.08:
                break
        sx = int(rng.choice((-1, 1)))
        sy = int(rng.choice((-1, 1)))
        sz = int(rng.choice((-1, 1)))
        amps = np.exp(rng.uniform(-3.0, 3.0, size=3)) * np.exp(1j * rng.uniform(-math.pi, math.pi, size=3))
        ax, ay, az = map(complex, amps)

        row = register_helical_physical_edge(x=x, y=y, z=z, sx=sx, sy=sy, sz=sz, ax=ax, ay=ay, az=az)
        wcoeff = max(wcoeff, row.direct_vs_waleffe_residual / max(1.0, abs(row.direct_child_source_coefficient), abs(row.waleffe_child_source_coefficient)))
        wleray = max(wleray, row.leray_pairing_residual / max(1.0, abs(row.direct_child_source_coefficient)))
        wident = max(wident, _relative(row.signed_upper_progress_work, row.registered_upper_progress_work))
        maxm = max(maxm, row.normalized_multiplier)
        if row.scale_progress == 0.0:
            nnon += 1
        if row.positive_forward_work:
            npos += 1
            minpos = min(minpos, row.phase_alignment)
            if row.phase_alignment <= -2e-12:
                raise AssertionError("positive forward physical work has negative registered phase alignment")

        swap = register_helical_physical_edge(x=y, y=x, z=z, sx=sy, sy=sx, sz=sz, ax=ay, ay=ax, az=az)
        wswap = max(
            wswap,
            _relative(row.direct_child_source_coefficient, swap.direct_child_source_coefficient),
            _relative(row.signed_child_energy_work, swap.signed_child_energy_work),
            _relative(row.native_modal_capacity, swap.native_modal_capacity),
            _relative(row.geometric_multiplier_J, swap.geometric_multiplier_J),
            _relative(row.phase_alignment, swap.phase_alignment),
        )

        g = coupling_g(x, y, -z, sx, sy, sz)
        th = rng.uniform(-math.pi, math.pi, size=3)
        gp, axp, ayp, azp = gauge_transform_modal_data(g, ax, ay, az, *map(float, th))
        cp = phase_alignment(row.signed_frequency_factor, gp, axp, ayp, azp)
        wgauge = max(wgauge, abs(cp - row.phase_alignment))

        lam = float(math.exp(rng.uniform(-4.0, 4.0)))
        scaled = register_helical_physical_edge(
            x=lam * x,
            y=lam * y,
            z=lam * z,
            sx=sx,
            sy=sy,
            sz=sz,
            ax=ax,
            ay=ay,
            az=az,
        )
        wscale = max(
            wscale,
            abs(scaled.normalized_multiplier - row.normalized_multiplier),
            abs(scaled.phase_alignment - row.phase_alignment),
            abs(scaled.forward_ratio - row.forward_ratio),
            _relative(scaled.native_modal_capacity, lam * row.native_modal_capacity),
            _relative(scaled.signed_child_energy_work, lam * row.signed_child_energy_work),
        )

    if minpos == math.inf:
        minpos = 0.0
    return HelicalPhysicalEdgeStress(
        samples=samples,
        worst_direct_waleffe_relative_residual=wcoeff,
        worst_leray_pairing_relative_residual=wleray,
        worst_upper_progress_relative_residual=wident,
        worst_parent_swap_relative_residual=wswap,
        worst_gauge_phase_alignment_residual=wgauge,
        worst_scale_invariance_residual=wscale,
        maximum_normalized_multiplier=maxm,
        minimum_positive_forward_phase_alignment=minpos,
        positive_forward_samples=npos,
        nonforward_samples=nnon,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-helical-physical-edge-registration"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "helical_physical_edge_registration.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Helical physical edge registration: native NS capacity before Young

Status: **{cert['status']}**.

For one actual unordered helical parent pair `x,y` feeding `z=x+y`, the nonlinear Navier--Stokes source is read directly as

`F_z=P_z(u_x x omega_y + u_y x omega_x)`.

The repository Waleffe convention gives exactly

`<h_z,F_z>=2(s_x|x|-s_y|y|)conj(g_e)a_xa_y`.

Because physical child energy differentiates as `2 Re(conj(a_z)<h_z,F_z>)`, the signed work is

`T_e=4(s_x|x|-s_y|y|) Re[conj(a_z)conj(g_e)a_xa_y]`.

Define the native modal interaction capacity

`A_e=4|z||a_xa_ya_z|`,

the geometric upper-progress multiplier

`J_e=log_+(|z|/p_top)|s_x|x|-s_y|y|| |g_e|/|z|`,

and the signed phase/orientation alignment `c_e in [-1,1]`.  Then on the **same physical Fourier event**,

`T_e log_+(|z|/p_top)=A_e J_e c_e`.

`A_e` is not a Young norm product.  It is the available modal interaction amplitude of this triad, including the exact factor two from the unordered parent orbit and factor two from physical energy differentiation.  Young/Christ enters only downstream when one asks whether a block of these real edges is near saturation.

Stress: `{out.samples}` random physical helical triads/amplitudes
- worst direct-Leray / Waleffe coefficient relative residual: `{out.worst_direct_waleffe_relative_residual:.3e}`
- worst Leray-free child pairing residual: `{out.worst_leray_pairing_relative_residual:.3e}`
- worst `T log = A J c` relative residual: `{out.worst_upper_progress_relative_residual:.3e}`
- worst unordered-parent swap residual: `{out.worst_parent_swap_relative_residual:.3e}`
- worst helical phase-gauge alignment residual: `{out.worst_gauge_phase_alignment_residual:.3e}`
- worst uniform-wavevector scale invariance residual: `{out.worst_scale_invariance_residual:.3e}`
- maximum sampled `J/J_*`: `{out.maximum_normalized_multiplier:.12g}`
- minimum phase alignment on sampled positive forward-work edges: `{out.minimum_positive_forward_phase_alignment:.3e}`
- positive forward-work samples: `{out.positive_forward_samples}`
- nonforward samples (zero upper-progress multiplier, retained rather than mislabeled signed-good): `{out.nonforward_samples}`

This theorem registers one physical edge only.  It does **not** yet construct the continuum edge measure or claim that generic/nonforward HH is signed-good.  No Navier--Stokes global-regularity conclusion is asserted.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
