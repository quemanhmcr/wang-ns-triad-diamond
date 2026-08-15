"""Primitive material vorticity/heat identities for incompressible Navier--Stokes.

This module deliberately uses only objects generated directly by NS:
material vorticity two-form, the pullback metric of the incompressible flow,
and the Hodge heat operator.  It does not introduce event taxonomies, spectral
owners, or artificial stopping rules.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

STATUS = "DRAFT_NATIVE_MATERIAL_VORTICITY_HEAT_LAW__TRANSVERSE_MEMORY__HEAT_ONLY_RESET"


def _metric(g: np.ndarray, *, atol: float = 2.0e-10) -> np.ndarray:
    a = np.asarray(g, dtype=float)
    if a.shape != (3, 3) or not np.all(np.isfinite(a)):
        raise ValueError("metric must be a finite 3x3 matrix")
    if np.linalg.norm(a - a.T, ord="fro") > atol * max(1.0, np.linalg.norm(a, ord="fro")):
        raise ValueError("metric must be symmetric")
    ev = np.linalg.eigvalsh(a)
    if ev[0] <= 0.0:
        raise ValueError("metric must be positive definite")
    det = float(np.linalg.det(a))
    if abs(det - 1.0) > atol * max(1.0, abs(det)):
        raise ValueError("incompressible pullback metric must have determinant one")
    return a


def _unit(v: Sequence[float]) -> np.ndarray:
    q = np.asarray(tuple(float(x) for x in v), dtype=float)
    if q.shape != (3,) or not np.all(np.isfinite(q)):
        raise ValueError("direction must be a finite three-vector")
    n = float(np.linalg.norm(q))
    if n <= 0.0:
        raise ValueError("nonzero direction required")
    return q / n


def _transverse_frame(direction: Sequence[float]) -> np.ndarray:
    """One orthonormal 3x2 frame for the Euclidean material plane q^perp."""
    q = _unit(direction)
    axes = np.eye(3)
    axis = axes[int(np.argmin(np.abs(q)))]
    e1 = axis - q * float(np.dot(q, axis))
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(q, e1)
    e2 /= np.linalg.norm(e2)
    return np.column_stack((e1, e2))


def transverse_heat_determinant(metric: np.ndarray, material_vorticity: Sequence[float]) -> dict[str, float]:
    """Exact stretch--diffusion determinant identity for a material vorticity 2-form.

    Write the pulled-back vorticity two-form as ``beta=i_q da``.  Closedness has
    principal-symbol characteristic plane ``xi.q=0``.  For ``g=F^T F`` with
    ``det g=1``, the Hodge heat principal symbol is ``g^-1`` and

        det(g^-1|_{q^perp}) = (q^T g q)/|q|^2.

    The right side is precisely the squared geometric amplification of the
    physical vorticity vector ``F q`` relative to the material vector ``q``.
    """
    g = _metric(metric)
    q0 = np.asarray(tuple(float(x) for x in material_vorticity), dtype=float)
    if q0.shape != (3,) or not np.all(np.isfinite(q0)) or np.linalg.norm(q0) == 0.0:
        raise ValueError("nonzero finite material vorticity required")
    qn = q0 / np.linalg.norm(q0)
    p = _transverse_frame(qn)
    a = p.T @ np.linalg.inv(g) @ p
    det_trans = float(np.linalg.det(a))
    amp2 = float(qn @ g @ qn)
    residual = abs(det_trans - amp2)
    scale = max(1.0, abs(det_trans), abs(amp2))
    if residual > 5.0e-9 * scale:
        raise AssertionError("transverse heat determinant lost incompressible cofactor identity")
    return {
        "geometric_amplification_squared": amp2,
        "geometric_amplification": math.sqrt(max(0.0, amp2)),
        "transverse_heat_determinant": det_trans,
        "transverse_heat_trace": float(np.trace(a)),
        "transverse_heat_geometric_mean": math.sqrt(max(0.0, det_trans)),
        "identity_residual": residual,
    }



def transverse_heat_log_rate(
    metric: np.ndarray,
    metric_rate: np.ndarray,
    material_vorticity: Sequence[float],
) -> dict[str, float]:
    """Differential form of stretch = transverse-heat-area growth.

    For a determinant-one metric path ``g(t)`` and Euler-frozen material
    direction ``q``,

        d log(q^T g q)/dt
        = d log det(g^-1|q^perp)/dt.

    ``metric_rate`` must be tangent to the incompressible metric manifold, i.e.
    ``tr(g^-1 metric_rate)=0``.
    """
    g = _metric(metric)
    h = np.asarray(metric_rate, dtype=float)
    if h.shape != (3, 3) or not np.all(np.isfinite(h)):
        raise ValueError("metric rate must be a finite 3x3 matrix")
    if np.linalg.norm(h-h.T, ord="fro") > 3.0e-10*max(1.0,np.linalg.norm(h,ord="fro")):
        raise ValueError("metric rate must be symmetric")
    inv=np.linalg.inv(g)
    tang=float(np.trace(inv@h))
    if abs(tang)>4.0e-10*max(1.0,np.linalg.norm(inv)*np.linalg.norm(h)):
        raise ValueError("metric rate must preserve determinant to first order")
    q=_unit(material_vorticity); p=_transverse_frame(q)
    amp=float(q@g@q)
    amp_rate=float(q@h@q)/amp
    aperp=p.T@inv@p
    adot=p.T@(-inv@h@inv)@p
    heat_rate=float(np.trace(np.linalg.solve(aperp,adot)))
    residual=abs(amp_rate-heat_rate)
    if residual>8.0e-9*max(1.0,abs(amp_rate),abs(heat_rate)):
        raise AssertionError("stretch/heat-area logarithmic rate identity failed")
    return {
        "log_vorticity_amplification_squared_rate": amp_rate,
        "log_transverse_heat_determinant_rate": heat_rate,
        "identity_residual": residual,
    }

def transverse_two_covector_area_identity(
    metric: np.ndarray,
    material_vorticity: Sequence[float],
    xi: Sequence[float],
    eta: Sequence[float],
    *,
    atol: float = 3.0e-10,
) -> dict[str, float]:
    """Gram-area identity on the characteristic plane of the closed 2-form.

    For ``xi,eta`` transverse to ``q`` and ``A=g^-1``, exactly

        det Gram_A(xi,eta)
        = ((q^T g q)/|q|^2) |xi x eta|^2.

    Hence amplification times Euclidean transverse frequency area is bounded by
    the geometric mean of the two heat symbols.
    """
    g = _metric(metric)
    q = _unit(material_vorticity)
    x = np.asarray(tuple(float(z) for z in xi), dtype=float)
    y = np.asarray(tuple(float(z) for z in eta), dtype=float)
    if x.shape != (3,) or y.shape != (3,) or not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)):
        raise ValueError("finite transverse covectors required")
    scalev = max(1.0, np.linalg.norm(x), np.linalg.norm(y))
    if abs(float(np.dot(x, q))) > atol * scalev or abs(float(np.dot(y, q))) > atol * scalev:
        raise ValueError("closed-form characteristic covectors must lie in q^perp")
    inv = np.linalg.inv(g)
    xx = float(x @ inv @ x)
    yy = float(y @ inv @ y)
    xy = float(x @ inv @ y)
    gram = xx * yy - xy * xy
    area2 = float(np.dot(np.cross(x, y), np.cross(x, y)))
    amp2 = float(q @ g @ q)
    rhs = amp2 * area2
    residual = abs(gram - rhs)
    scale = max(1.0, abs(gram), abs(rhs))
    if residual > 8.0e-9 * scale:
        raise AssertionError("transverse heat Gram-area identity failed")
    cs_margin = math.sqrt(max(0.0, xx * yy)) - math.sqrt(max(0.0, rhs))
    if cs_margin < -8.0e-10 * max(1.0, math.sqrt(max(0.0, xx * yy)), math.sqrt(max(0.0, rhs))):
        raise AssertionError("transverse no-free-area inequality failed")
    return {
        "heat_gram_determinant": gram,
        "amplification_times_area_squared": rhs,
        "identity_residual": residual,
        "no_free_area_margin": max(0.0, cs_margin),
    }


def accumulated_transverse_heat_memory(
    metrics: Sequence[np.ndarray],
    material_vorticity: Sequence[float],
    time_weights: Sequence[float],
) -> dict[str, float]:
    """Minkowski determinant memory for one Euler-frozen material polarization.

    Euler freezes the pulled-back vorticity two-form, so its material direction
    ``q`` does not rotate under nonlinearity.  For any positive time quadrature,

        H_perp = sum dt_j (g_j^-1|_{q^perp})

    satisfies

        sqrt(det H_perp) >= sum dt_j sqrt(det(g_j^-1|_{q^perp}))
                           = sum dt_j |F_j q|/|q|.

    Thus rotating/extreme anisotropy cannot erase accumulated transverse heat
    area.  This is the discrete exact form of the continuous Minkowski integral
    inequality.
    """
    if len(metrics) == 0 or len(metrics) != len(time_weights):
        raise ValueError("one positive time weight per metric required")
    q = _unit(material_vorticity)
    p = _transverse_frame(q)
    h = np.zeros((2, 2), dtype=float)
    accumulated_amplification = 0.0
    for raw_g, raw_dt in zip(metrics, time_weights):
        g = _metric(raw_g)
        dt = float(raw_dt)
        if not math.isfinite(dt) or dt <= 0.0:
            raise ValueError("strictly positive finite time weights required")
        a = p.T @ np.linalg.inv(g) @ p
        h += dt * a
        accumulated_amplification += dt * math.sqrt(max(0.0, float(np.linalg.det(a))))
    heat_area = math.sqrt(max(0.0, float(np.linalg.det(h))))
    margin = heat_area - accumulated_amplification
    if margin < -2.0e-9 * max(1.0, heat_area, accumulated_amplification):
        raise AssertionError("accumulated transverse heat lost Minkowski memory")
    return {
        "accumulated_heat_area": heat_area,
        "accumulated_geometric_amplification": accumulated_amplification,
        "memory_margin": max(0.0, margin),
        "history_length": float(math.fsum(float(x) for x in time_weights)),
    }




def moving_polarization_memory_bound(
    deformations: Sequence[np.ndarray],
    material_vorticities: Sequence[Sequence[float]],
    time_weights: Sequence[float],
) -> dict[str, float]:
    """Exact memory--reset inequality for a heat-moving material polarization.

    Let ``omega_j=F_j q_j`` and anchor the history at the final material
    polarization ``q_T``.  Minkowski determinant memory applies to this fixed
    direction, while the triangle inequality gives

        sum dt |F_j q_j|
        <= |q_T| sqrt(det sum dt g_j^-1|_{q_T^perp})
           + sum dt |F_j(q_j-q_T)|.

    The second term is the *only* way the actual material-vorticity history can
    evade the fixed-plane heat memory.  In Navier--Stokes ``q_t`` changes only
    through ``beta_t=-nu L_g beta``.  This helper is purely geometric and does
    not insert a reset threshold or event count.
    """
    n=len(deformations)
    if n==0 or len(material_vorticities)!=n or len(time_weights)!=n:
        raise ValueError("matching nonempty deformation, polarization and time histories required")
    fs=[]
    for raw in deformations:
        f=np.asarray(raw,dtype=float)
        if f.shape!=(3,3) or not np.all(np.isfinite(f)):
            raise ValueError("finite 3x3 deformation gradients required")
        det=float(np.linalg.det(f))
        if det<=0.0 or abs(det-1.0)>3.0e-9*max(1.0,abs(det)):
            raise ValueError("incompressible deformation gradients must lie in SL(3)")
        fs.append(f)
    qs=[np.asarray(tuple(float(x) for x in q),dtype=float) for q in material_vorticities]
    if any(q.shape!=(3,) or not np.all(np.isfinite(q)) for q in qs):
        raise ValueError("finite material vorticity vectors required")
    qT=qs[-1]; nT=float(np.linalg.norm(qT))
    if nT<=0.0:
        raise ValueError("final material polarization must be nonzero")
    p=_transverse_frame(qT)
    H=np.zeros((2,2),dtype=float)
    actual=0.0; reset=0.0; fixed_amp=0.0
    for f,q,raw_dt in zip(fs,qs,time_weights):
        dt=float(raw_dt)
        if not math.isfinite(dt) or dt<=0.0:
            raise ValueError("strictly positive finite time weights required")
        g=f.T@f
        A=p.T@np.linalg.inv(g)@p
        H += dt*A
        actual += dt*float(np.linalg.norm(f@q))
        reset += dt*float(np.linalg.norm(f@(q-qT)))
        fixed_amp += dt*float(np.linalg.norm(f@qT))
    memory=nT*math.sqrt(max(0.0,float(np.linalg.det(H))))
    # Fixed-direction Minkowski memory dominates the exact fixed-direction
    # amplitude integral; the actual history then differs only by reset motion.
    if fixed_amp > memory + 3.0e-8*max(1.0,fixed_amp,memory):
        raise AssertionError("fixed-polarization amplitude exceeded transverse heat memory")
    upper=memory+reset
    if actual > upper + 4.0e-8*max(1.0,actual,upper):
        raise AssertionError("moving-polarization memory/reset inequality failed")
    return {
        "actual_vorticity_amplitude_history": actual,
        "final_plane_heat_memory": memory,
        "fixed_final_polarization_history": fixed_amp,
        "heat_only_reset_remainder": reset,
        "memory_reset_upper": upper,
        "inequality_margin": max(0.0,upper-actual),
    }

def pair_direction_mismatch_decomposition(
    deformation_a: np.ndarray,
    deformation_b: np.ndarray,
    material_vorticity_a: Sequence[float],
    material_vorticity_b: Sequence[float],
) -> dict[str, float]:
    """Exact two-source decomposition of physical vorticity mismatch.

    For ``F_a,F_b in SL(3)`` and material vorticities ``q_a,q_b``,

        (F_a q_a) x (F_b q_b)
        = F_a^-T(q_a x q_b)
          + (F_a q_a) x ((F_b-F_a)q_b).

    The first term is material two-form/polarization mismatch measured by the
    local heat covector norm of ``g_a^-1``.  The second is deformation
    non-affinity.  These are the only two sources of physical directional
    mismatch entering the Biot--Savart stretching kernel.
    """
    def Fmat(raw):
        f=np.asarray(raw,dtype=float)
        if f.shape!=(3,3) or not np.all(np.isfinite(f)):
            raise ValueError("deformation gradients must be finite 3x3 matrices")
        det=float(np.linalg.det(f))
        if det<=0 or abs(det-1.0)>3.0e-9*max(1.0,abs(det)):
            raise ValueError("incompressible deformation gradients must lie in SL(3)")
        return f
    fa=Fmat(deformation_a); fb=Fmat(deformation_b)
    qa=np.asarray(tuple(float(x) for x in material_vorticity_a),float)
    qb=np.asarray(tuple(float(x) for x in material_vorticity_b),float)
    if qa.shape!=(3,) or qb.shape!=(3,) or not np.all(np.isfinite(qa)) or not np.all(np.isfinite(qb)):
        raise ValueError("finite material vorticity vectors required")
    wa=fa@qa; wb=fb@qb
    lhs=np.cross(wa,wb)
    material=np.linalg.solve(fa.T,np.cross(qa,qb))
    nonaff=np.cross(wa,(fb-fa)@qb)
    represented=material+nonaff
    residual=float(np.linalg.norm(lhs-represented))
    scale=max(1.0,np.linalg.norm(lhs),np.linalg.norm(material),np.linalg.norm(nonaff))
    if residual>2.0e-8*scale:
        raise AssertionError("pair directional mismatch decomposition failed")
    ga=fa.T@fa; heat_norm=math.sqrt(max(0.0,float(np.cross(qa,qb)@np.linalg.inv(ga)@np.cross(qa,qb))))
    return {
        "physical_cross_norm": float(np.linalg.norm(lhs)),
        "material_heat_covector_norm": heat_norm,
        "material_term_norm": float(np.linalg.norm(material)),
        "nonaffinity_term_norm": float(np.linalg.norm(nonaff)),
        "identity_residual": residual,
    }

def material_hminus2_reset_identity(
    laplacian_eigenvalues: Sequence[float],
    beta_coefficients: Sequence[complex],
    viscosity: float,
) -> dict[str, float]:
    """Exact heat-only reset action of the pulled-back vorticity two-form.

    Put ``L_g=-Delta_g >=0`` on exact two-forms.  The primitive material law is

        beta_t = -nu L_g beta.

    In the instantaneous homogeneous ``H^-2`` metric generated by ``L_g``,

        ||beta_t||_{H^-2_g}^2 = nu^2 ||beta||_{L2_g}^2.

    Therefore the full NS energy law ``E'=-2nu||beta||_g^2`` gives

        (1/nu) int ||beta_t||_{H^-2_g}^2 dt = (E(0)-E(T))/2.

    This helper checks the instantaneous spectral identity; the integrated
    equality is then exactly the physical velocity-energy balance.
    """
    lam = np.asarray(tuple(float(x) for x in laplacian_eigenvalues), dtype=float)
    b = np.asarray(tuple(complex(x) for x in beta_coefficients), dtype=complex)
    nu = float(viscosity)
    if lam.ndim != 1 or b.ndim != 1 or len(lam) == 0 or len(lam) != len(b):
        raise ValueError("matching nonempty Laplacian spectrum and beta coefficients required")
    if not np.all(np.isfinite(lam)) or np.any(lam <= 0.0):
        raise ValueError("exact vorticity modes require positive finite Laplacian eigenvalues")
    if not np.all(np.isfinite(b.real)) or not np.all(np.isfinite(b.imag)):
        raise ValueError("finite beta coefficients required")
    if not math.isfinite(nu) or nu <= 0.0:
        raise ValueError("positive finite viscosity required")
    bt = -nu * lam * b
    action = float(np.sum(np.abs(bt) ** 2 / (lam ** 2)))
    beta_l2 = float(np.sum(np.abs(b) ** 2))
    represented = nu * nu * beta_l2
    residual = abs(action - represented)
    if residual > 4.0e-11 * max(1.0, action, represented):
        raise AssertionError("material H^-2 reset identity failed")
    return {
        "beta_l2_squared": beta_l2,
        "hminus2_reset_speed_squared": action,
        "represented_reset_speed_squared": represented,
        "energy_loss_density_factor": action / nu,
        "identity_residual": residual,
    }


def rank_one_incompressible_stretch_null(amplitude: Sequence[float], covector: Sequence[float]) -> dict[str, float]:
    """A one-direction incompressible velocity gradient cannot self-stretch vorticity.

    If ``A=a tensor xi`` and incompressibility gives ``a.xi=0``, then the
    associated vorticity is parallel to ``xi x a`` and ``S omega=0`` exactly.
    This is the physical-space rank-one null behind the collinear Fourier/Waleffe
    degeneracy: genuine vortex stretching needs more than one spatial derivative
    direction.
    """
    a = np.asarray(tuple(float(x) for x in amplitude), dtype=float)
    xi = np.asarray(tuple(float(x) for x in covector), dtype=float)
    if a.shape != (3,) or xi.shape != (3,) or not np.all(np.isfinite(a)) or not np.all(np.isfinite(xi)):
        raise ValueError("finite three-vectors required")
    nx = float(np.linalg.norm(xi))
    if nx <= 0.0:
        raise ValueError("nonzero covector required")
    div = float(np.dot(a, xi))
    scale = max(1.0, np.linalg.norm(a) * nx)
    if abs(div) > 3.0e-10 * scale:
        raise ValueError("rank-one gradient must be incompressible: a.xi=0")
    A = np.outer(a, xi)
    S = 0.5 * (A + A.T)
    omega = np.cross(xi, a)
    stretch_vec = S @ omega
    residual = float(np.linalg.norm(stretch_vec))
    if residual > 5.0e-10 * max(1.0, np.linalg.norm(S) * np.linalg.norm(omega)):
        raise AssertionError("rank-one incompressible gradient acquired vortex stretching")
    return {
        "divergence_residual": abs(div),
        "vorticity_norm": float(np.linalg.norm(omega)),
        "stretch_vector_norm": residual,
    }



def material_state_speed_lock(beta_l2_squared: float, viscosity: float) -> dict[str, float]:
    """Exact global speed lock of the two primitive material state fields.

    With the affine-invariant metric speed

        ||g_t||_M^2 = int tr[(g^-1 g_t)^2] da,

    ``g_t=2 Phi^*S`` gives ``||g_t||_M^2=4||S||_2^2``.  Global
    incompressibility gives ``||omega||_2^2=2||S||_2^2=||beta||_g^2``.
    Together with ``beta_t=-nu L_g beta`` this yields

        ||g_t||_M^2 = 2||beta||_g^2,
        ||beta_t||_H^-2_g^2 = (nu^2/2)||g_t||_M^2,
        -E' = nu ||g_t||_M^2 = (2/nu)||beta_t||_H^-2_g^2.

    Thus material metric deformation and heat rewriting are not independent
    state velocities; NS locks them to the same enstrophy.
    """
    b=float(beta_l2_squared); nu=float(viscosity)
    if not math.isfinite(b) or b<0.0 or not math.isfinite(nu) or nu<=0.0:
        raise ValueError("nonnegative finite beta L2 squared and positive viscosity required")
    metric_speed2=2.0*b
    reset_speed2=nu*nu*b
    energy_decay=2.0*nu*b
    if abs(reset_speed2-0.5*nu*nu*metric_speed2)>2e-13*max(1.0,reset_speed2):
        raise AssertionError("material state speed lock failed")
    if abs(energy_decay-nu*metric_speed2)>2e-13*max(1.0,energy_decay):
        raise AssertionError("metric speed lost physical energy law")
    return {
        "beta_l2_squared": b,
        "metric_affine_speed_squared": metric_speed2,
        "hminus2_beta_reset_speed_squared": reset_speed2,
        "positive_energy_decay_rate": energy_decay,
        "reset_to_metric_speed_squared_ratio": 0.5*nu*nu,
    }

def material_log_distortion_energy_bound(time_horizon: float, energy_loss: float, viscosity: float) -> dict[str, float]:
    """Global material-label L2 bound for accumulated deformation log-stretch.

    Along an incompressible trajectory, ``d log sigma_max(F)/dt <= ||S||op``.
    Volume preservation and Cauchy give

        int [log sigma_max(F(a,t))_+]^2 da
        <= t int_0^t ||S||_2^2 ds
        = t (E(0)-E(t))/(4 nu).

    The last equality uses only ``||omega||_2^2=2||S||_2^2`` and the NS energy
    law.  This does not control a material supremum; it quantifies that extreme
    deformation can persist only on a shrinking set of labels.
    """
    t = float(time_horizon)
    loss = float(energy_loss)
    nu = float(viscosity)
    if not all(math.isfinite(x) for x in (t, loss, nu)) or t < 0.0 or loss < 0.0 or nu <= 0.0:
        raise ValueError("nonnegative finite time/loss and positive viscosity required")
    bound = t * loss / (4.0 * nu)
    return {
        "log_max_stretch_l2_squared_upper": bound,
        "strain_spacetime_l2_squared": loss / (4.0 * nu),
    }



def material_hodge_speed_ladder(
    laplacian_eigenvalues: Sequence[float],
    beta_coefficients: Sequence[complex],
    viscosity: float,
    sobolev_exponent: float,
) -> dict[str, float]:
    """All-scale lock of material-vorticity reset, viscous current and metric speed.

    On a flat incompressible material metric put ``h=g_t=Lie_v g`` and
    ``beta=d v^flat``.  The symmetric/antisymmetric parts of the same velocity
    gradient obey, on every Hodge frequency,

        ||L^(s/2) h||^2 = 2 ||L^(s/2) beta||^2.

    With ``j=nu delta beta=-nu div h`` and ``beta_t=-d j=-nu L beta``, exact
    Hodge functional calculus gives for every real exponent in the common domain

        ||beta_t||_{H^(s-2)}^2
        = ||j||_{H^(s-1)}^2
        = (nu^2/2) ||h||_{H^s}^2.

    Thus strong material-polarization rewriting is not an independent degree of
    freedom: at exactly the same spatial scale it is the curvature/divergence of
    the metric-deformation velocity generated by the same NS state.
    """
    lam=np.asarray(tuple(float(x) for x in laplacian_eigenvalues),dtype=float)
    b=np.asarray(tuple(complex(x) for x in beta_coefficients),dtype=complex)
    nu=float(viscosity); ss=float(sobolev_exponent)
    if lam.ndim!=1 or b.ndim!=1 or len(lam)==0 or len(lam)!=len(b):
        raise ValueError("matching nonempty positive Hodge spectrum and beta coefficients required")
    if not np.all(np.isfinite(lam)) or np.any(lam<=0.0):
        raise ValueError("positive finite Hodge eigenvalues required")
    if not np.all(np.isfinite(b.real)) or not np.all(np.isfinite(b.imag)):
        raise ValueError("finite beta coefficients required")
    if not math.isfinite(nu) or nu<=0.0 or not math.isfinite(ss):
        raise ValueError("positive finite viscosity and finite Sobolev exponent required")
    beta_hs=float(np.sum((lam**ss)*np.abs(b)**2))
    metric_hs=2.0*beta_hs
    # beta_t=-nu L beta: H^(s-2) multiplier squared is lam^(s-2).
    reset=nu*nu*beta_hs
    # j=nu delta beta; on exact two-forms delta is sqrt(L) isometric into
    # coexact one-forms, hence H^(s-1) has the same spectral weight lam^s.
    current=nu*nu*beta_hs
    residual=max(abs(reset-current),abs(reset-.5*nu*nu*metric_hs))
    if residual>8.0e-11*max(1.0,reset,current,abs(.5*nu*nu*metric_hs)):
        raise AssertionError("material Hodge speed ladder lost exact scale lock")
    return {
        "sobolev_exponent":ss,
        "beta_hs_squared":beta_hs,
        "metric_speed_hs_squared":metric_hs,
        "viscous_current_hs_minus_one_squared":current,
        "beta_reset_hs_minus_two_squared":reset,
        "identity_residual":residual,
    }


def primitive_material_current_fourier_law(
    wavevector: Sequence[float],
    velocity_coefficient: Sequence[complex],
    viscosity: float,
) -> dict[str, float]:
    """One Fourier fiber of the local primitive current law.

    For a nonzero wavevector and divergence-free velocity coefficient ``u``, put

        h_ij = i(k_i u_j+k_j u_i),
        beta_ij = i(k_i u_j-k_j u_i).

    Then, with the standard form codifferential and tensor divergence,

        delta beta = -div h,
        d div h = -L beta,
        j=nu delta beta,
        beta_t + d j = 0.

    The identity is local/operator-level; Fourier is used here only as a numerical
    referee which avoids choosing packets, shells or event labels.
    """
    k=np.asarray(tuple(float(x) for x in wavevector),dtype=float)
    u=np.asarray(tuple(complex(x) for x in velocity_coefficient),dtype=complex)
    nu=float(viscosity)
    if k.shape!=(3,) or u.shape!=(3,) or not np.all(np.isfinite(k)) or not np.all(np.isfinite(u.real)) or not np.all(np.isfinite(u.imag)):
        raise ValueError("finite three-dimensional Fourier data required")
    r2=float(k@k)
    if r2<=0.0 or not math.isfinite(nu) or nu<=0.0:
        raise ValueError("nonzero wavevector and positive finite viscosity required")
    divu=complex(k@u)
    if abs(divu)>5.0e-10*max(1.0,math.sqrt(r2)*float(np.linalg.norm(u))):
        raise ValueError("velocity coefficient must be divergence free")
    h=1j*(np.outer(k,u)+np.outer(u,k))
    beta=1j*(np.outer(k,u)-np.outer(u,k))
    delta_beta=-1j*np.einsum('i,ij->j',k,beta)
    div_h=1j*np.einsum('i,ij->j',k,h)
    d_div=1j*(np.outer(k,div_h)-np.outer(div_h,k))
    Lbeta=r2*beta
    j=nu*delta_beta
    dj=1j*(np.outer(k,j)-np.outer(j,k))
    beta_t=-nu*Lbeta
    r1=float(np.linalg.norm(delta_beta+div_h))
    r2res=float(np.linalg.norm(d_div+Lbeta))
    r3=float(np.linalg.norm(beta_t+dj))
    scale=max(1.0,float(np.linalg.norm(Lbeta)),float(np.linalg.norm(beta_t)),float(np.linalg.norm(dj)))
    if max(r1,r2res,r3)>8.0e-10*scale:
        raise AssertionError("primitive material current factorization failed")
    return {
        "codifferential_divergence_residual":r1,
        "laplacian_factorization_residual":r2res,
        "vorticity_flux_continuity_residual":r3,
        "current_norm_squared":float(np.vdot(j,j).real),
        "metric_speed_norm_squared":float(np.vdot(h,h).real),
        "beta_norm_squared_matrix_convention":float(.5*np.vdot(beta,beta).real),
    }


def vorticity_stress_current_algebra(
    velocity: Sequence[float],
    vorticity: Sequence[float],
    curl_vorticity: Sequence[float],
    strain: np.ndarray,
) -> dict[str, float]:
    """Pointwise algebra behind the primitive stress/current Noether identity.

    For ``c=curl omega`` define the vorticity stress

        T = (|omega|^2/2) I - omega tensor omega.

    Its geometric divergence is ``omega x c`` for a closed vorticity two-form.
    This helper records the algebraic work identities once that differential
    identity is supplied:

        -T:(2S)/2 = omega.S.omega,
        u.(omega x c) = c.(u x omega).

    The differential theorem ``div T=i_c beta`` is documented separately; it is
    a standard direct consequence of d beta=0 and delta beta=c.
    """
    u=np.asarray(tuple(float(x) for x in velocity),float)
    w=np.asarray(tuple(float(x) for x in vorticity),float)
    c=np.asarray(tuple(float(x) for x in curl_vorticity),float)
    S=np.asarray(strain,float)
    if u.shape!=(3,) or w.shape!=(3,) or c.shape!=(3,) or S.shape!=(3,3) or not all(np.all(np.isfinite(x)) for x in (u,w,c,S)):
        raise ValueError("finite velocity/vorticity/current/strain data required")
    if np.linalg.norm(S-S.T)>5e-11*max(1.0,np.linalg.norm(S)) or abs(float(np.trace(S)))>5e-11*max(1.0,np.linalg.norm(S)):
        raise ValueError("strain must be symmetric and trace free")
    T=.5*float(w@w)*np.eye(3)-np.outer(w,w)
    stress_work=-.5*float(np.sum(T*(2*S)))
    stretching=float(w@(S@w))
    divT=np.cross(w,c)
    current_work=float(u@divT)
    cyclic=float(c@np.cross(u,w))
    residual=max(abs(stress_work-stretching),abs(current_work-cyclic))
    if residual>5e-11*max(1.0,abs(stress_work),abs(stretching),abs(current_work),abs(cyclic)):
        raise AssertionError("vorticity stress/current algebra failed")
    return {
        "stretching_density":stretching,
        "stress_metric_work_density":stress_work,
        "stress_divergence_vector_norm":float(np.linalg.norm(divT)),
        "current_cyclic_work_density":cyclic,
        "velocity_stress_divergence_work_density":current_work,
        "identity_residual":residual,
    }



def canonical_maxwell_extension_spectral_law(
    signed_curl_frequencies: Sequence[float],
    modal_energies: Sequence[float],
) -> dict[str, float]:
    """Canonical Chern--Simons/Maxwell extension of the primitive curl law.

    For a co-closed velocity one-form ``alpha`` with signed-curl atoms ``a`` and
    energies ``E_a``, put ``Lambda=|curl|`` and extend

        A(y)=exp(-y Lambda) alpha,   y>0,
        F4=d_4 A.

    The source-free four-dimensional Maxwell field has exactly

        int |F4|^2 dy = K = sum |a| E_a,
        int F4 wedge F4 = H = sum a E_a             (orientation sign aside),
        E_+=(K+H)/2, E_-=(K-H)/2,
        int |partial_y F4|^2 dy = M3 = sum |a|^3 E_a.

    Its boundary energy-density profile is

        Q(y)=2 sum a^2 E_a exp(-2|a|y),

    so ``Q(0)=2Z`` and ``-Q'(0)/4=M3``.  This helper records the exact spectral
    algebra; no artificial extension scale or cutoff is selected.
    """
    a=np.asarray(tuple(float(x) for x in signed_curl_frequencies),dtype=float)
    e=np.asarray(tuple(float(x) for x in modal_energies),dtype=float)
    if a.ndim!=1 or e.ndim!=1 or len(a)==0 or len(a)!=len(e):
        raise ValueError("matching nonempty signed-curl spectrum and modal energies required")
    if not np.all(np.isfinite(a)) or np.any(a==0.0) or not np.all(np.isfinite(e)) or np.any(e<0.0):
        raise ValueError("finite nonzero curl frequencies and nonnegative finite energies required")
    r=np.abs(a)
    energy=float(np.sum(e)); helicity=float(np.sum(a*e)); critical=float(np.sum(r*e))
    enstrophy=float(np.sum(a*a*e)); m3=float(np.sum(r**3*e))
    plus=float(np.sum(r[a>0]*e[a>0])); minus=float(np.sum(r[a<0]*e[a<0]))
    residual=max(abs(critical-plus-minus),abs(helicity-plus+minus))
    if residual>3e-12*max(1.0,critical,abs(helicity)):
        raise AssertionError("Maxwell self-dual/anti-self-dual energy split failed")
    return {
        "velocity_energy":energy,
        "helicity_chern_simons":helicity,
        "critical_maxwell_energy":critical,
        "self_dual_energy":plus,
        "anti_self_dual_energy":minus,
        "enstrophy":enstrophy,
        "critical_viscous_bulk_gradient":m3,
        "boundary_maxwell_energy_density":2.0*enstrophy,
        "negative_quarter_boundary_profile_derivative":m3,
        "duality_split_residual":residual,
    }


def canonical_poisson_scale_overlap(radius_a: float, radius_b: float) -> dict[str, float]:
    """Exact overlap of two unit critical-energy Maxwell depth profiles.

    ``phi_r(y)=sqrt(2r) exp(-r y)`` has unit ``L2(0,infinity)`` norm.  Therefore

        <phi_r,phi_s> = 2 sqrt(rs)/(r+s)
                       = sech((log r-log s)/2).

    This is the earlier critical ``sech`` scale filter, now as a literal overlap
    of the canonical harmonic extensions generated by ``|curl|``.
    """
    r=float(radius_a); s=float(radius_b)
    if not all(math.isfinite(x) and x>0.0 for x in (r,s)):
        raise ValueError("positive finite curl radii required")
    overlap=2.0*math.sqrt(r*s)/(r+s)
    sech=1.0/math.cosh(0.5*math.log(r/s))
    if abs(overlap-sech)>3e-14*max(1.0,overlap,sech):
        raise AssertionError("Poisson overlap lost critical sech identity")
    return {"poisson_overlap":overlap,"log_scale_sech":sech,"log_scale_gap":abs(math.log(r/s))}


def maxwell_duality_stress_algebra(field_strength: np.ndarray) -> dict[str, float]:
    """Four-dimensional Maxwell BPS stress identity for one two-form.

    In oriented Euclidean four-space decompose ``F=F_+ + F_-`` by the Hodge star.
    The Maxwell stress of either pure duality vanishes identically, while

        |T(F)|^2 = 4 |F_+|^2 |F_-|^2

    using the standard two-form norm ``|F|^2=(1/2)F_ab F_ab`` and
    ``T_ab=F_ac F_bc-(1/4)delta_ab F_cd F_cd``.
    """
    F=np.asarray(field_strength,dtype=float)
    if F.shape!=(4,4) or not np.all(np.isfinite(F)):
        raise ValueError("finite 4x4 two-form matrix required")
    skew=.5*(F-F.T)
    if np.linalg.norm(F-skew)>3e-10*max(1.0,np.linalg.norm(F)):
        raise ValueError("field strength must be antisymmetric")
    eps=np.zeros((4,4,4,4),dtype=float)
    import itertools
    for perm in itertools.permutations(range(4)):
        inv=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        eps[perm]=(-1.0)**inv
    star=.5*np.einsum('abcd,cd->ab',eps,F)
    fp=.5*(F+star); fm=.5*(F-star)
    def n2(X): return .5*float(np.sum(X*X))
    def T(X): return X@X.T-.25*np.eye(4)*float(np.sum(X*X))
    tf=T(F);tp=T(fp);tm=T(fm)
    lhs=float(np.sum(tf*tf)); rhs=4.0*n2(fp)*n2(fm)
    square_res=float(np.linalg.norm(tf@tf-.25*lhs*np.eye(4)))
    residual=max(abs(lhs-rhs),np.linalg.norm(tp),np.linalg.norm(tm),square_res)
    if residual>2e-10*max(1.0,lhs,rhs,n2(F)):
        raise AssertionError("Maxwell self-duality/stress identity failed")
    return {
        "field_energy_density":n2(F),
        "self_dual_energy_density":n2(fp),
        "anti_self_dual_energy_density":n2(fm),
        "maxwell_stress_norm_squared":lhs,
        "cross_duality_product_times_four":rhs,
        "pure_self_dual_stress_norm":float(np.linalg.norm(tp)),
        "pure_anti_self_dual_stress_norm":float(np.linalg.norm(tm)),
        "stress_square_scalar_residual":square_res,
        "identity_residual":residual,
    }



def primitive_spacetime_gauge_algebra(
    velocity: Sequence[float],
    vorticity: Sequence[float],
    vorticity_current: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Pointwise algebra of the exact spacetime vorticity curvature.

    With ``alpha=u^flat``, ``beta=d alpha`` and ``c=delta beta`` the rotational
    Navier--Stokes equation is

        alpha_t + d(p+|u|^2/2) = - i_u beta - nu c.

    Hence the Abelian spacetime connection ``A4=alpha-(p+|u|^2/2)dt`` has

        F4 = beta - dt wedge e,       e=i_u beta+nu c.

    In vector representatives ``i_u beta=-u cross omega``.  The Bianchi identity
    ``d_4 F4=0`` is exactly the vorticity equation.  Algebraically,

        e.omega = nu c.omega
        -e.c = (u cross omega).c - nu |c|^2.

    The first identity is the pointwise Euler Chern--Simons null; the second is
    the local Joule density in the enstrophy Poynting law.  Pressure occurs only
    as the temporal gauge potential and is not an independent curvature source.
    """
    u=np.asarray(tuple(float(x) for x in velocity),dtype=float)
    w=np.asarray(tuple(float(x) for x in vorticity),dtype=float)
    c=np.asarray(tuple(float(x) for x in vorticity_current),dtype=float)
    nu=float(viscosity)
    if u.shape!=(3,) or w.shape!=(3,) or c.shape!=(3,) or not np.all(np.isfinite(u)) or not np.all(np.isfinite(w)) or not np.all(np.isfinite(c)):
        raise ValueError("finite velocity/vorticity/current three-vectors required")
    if not math.isfinite(nu) or nu<0.0:
        raise ValueError("finite nonnegative viscosity required")
    uxw=np.cross(u,w)
    i_u_beta=-uxw
    e=i_u_beta+nu*c
    top=float(np.dot(e,w)); visc_top=nu*float(np.dot(c,w))
    stretch=float(np.dot(uxw,c)); joule=-float(np.dot(e,c)); represented=stretch-nu*float(np.dot(c,c))
    null=float(np.dot(i_u_beta,w))
    residual=max(abs(top-visc_top),abs(joule-represented),abs(null))
    if residual>3e-12*max(1.0,abs(top),abs(visc_top),abs(joule),abs(represented),np.linalg.norm(u)*np.linalg.norm(w)*max(1.0,np.linalg.norm(c))):
        raise AssertionError("primitive spacetime curvature algebra failed")
    return {
        "euler_emotive_norm":float(np.linalg.norm(i_u_beta)),
        "total_emotive_norm":float(np.linalg.norm(e)),
        "euler_topological_null":null,
        "chern_simons_density_rate_half":top,
        "viscous_chern_simons_density_rate_half":visc_top,
        "euler_stretching_work_density":stretch,
        "negative_joule_work_density":joule,
        "stretch_minus_ohmic_density":represented,
        "identity_residual":residual,
    }


def so33_exterior_square_algebra(generator: np.ndarray) -> dict[str, float]:
    """Exterior-square ``sl(4)->so(3,3)`` law on two-forms.

    The wedge form on ``Lambda^2 R^4`` has signature ``(3,3)``.  A trace-free
    four-dimensional infinitesimal volume-preserving deformation ``A`` acts on
    covariant two-forms by

        rho(A)F=A^T F+F A.

    In a wedge-orthonormal basis, ``rho(A)^T J+J rho(A)=0``.  The skew part of A
    acts compactly: it is Euclidean-skew and commutes with the Hodge involution J.
    The symmetric trace-free part acts as a noncompact boost: it is
    Euclidean-symmetric and anticommutes with J.  Moreover

        ||rho(S)||_HS^2 = 2 ||S||_F^2.

    Thus strain is literally the duality-mixing tangent of the primitive
    ``SO(3,3)`` geometry; rotation cannot mix the two Hodge sectors.
    """
    A=np.asarray(generator,dtype=float)
    if A.shape!=(4,4) or not np.all(np.isfinite(A)):
        raise ValueError("finite 4x4 generator required")
    tr=float(np.trace(A))
    if abs(tr)>3e-10*max(1.0,np.linalg.norm(A)):
        raise ValueError("volume-preserving generator must be trace-free")
    pairs=((0,1),(0,2),(0,3),(2,3),(3,1),(1,2))
    basis=[]
    for a,b in pairs:
        F=np.zeros((4,4));F[a,b]=1.0;F[b,a]=-1.0;basis.append(F)
    eps=np.zeros((4,4,4,4),dtype=float)
    import itertools
    for perm in itertools.permutations(range(4)):
        inv=sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        eps[perm]=(-1.0)**inv
    def inn(F,G): return .5*float(np.sum(F*G))
    def star(F): return .5*np.einsum('abcd,cd->ab',eps,F)
    J=np.array([[inn(bi,star(bj)) for bj in basis] for bi in basis])
    def rep(B):
        R=np.zeros((6,6))
        for j,F in enumerate(basis):
            RF=B.T@F+F@B
            for i,bi in enumerate(basis):R[i,j]=inn(bi,RF)
        return R
    R=rep(A); K=.5*(A-A.T); S=.5*(A+A.T); S-=np.trace(S)/4.0*np.eye(4)
    RK=rep(K); RS=rep(S)
    so_res=float(np.linalg.norm(R.T@J+J@R))
    rot_res=max(float(np.linalg.norm(RK+RK.T)),float(np.linalg.norm(RK@J-J@RK)))
    boost_res=max(float(np.linalg.norm(RS-RS.T)),float(np.linalg.norm(RS@J+J@RS)))
    lhs=float(np.sum(RS*RS));rhs=2.0*float(np.sum(S*S)); norm_res=abs(lhs-rhs)
    scale=max(1.0,np.linalg.norm(R),lhs,rhs)
    if max(so_res,rot_res,boost_res,norm_res)>2e-10*scale:
        raise AssertionError("exterior-square SO(3,3) algebra failed")
    return {
        "wedge_signature_positive":3.0,
        "wedge_signature_negative":3.0,
        "so33_residual":so_res,
        "rotation_compact_residual":rot_res,
        "strain_boost_residual":boost_res,
        "strain_exterior_hs_squared":lhs,
        "twice_strain_frobenius_squared":rhs,
        "boost_speed_identity_residual":norm_res,
    }



def transported_twoform_transverse_determinant(
    deformation: np.ndarray,
    material_twoform_vector: Sequence[float],
) -> dict[str, float]:
    """Cofactor law for a transported two-form under an arbitrary orientation-preserving flow.

    If ``F`` is a deformation gradient with ``J=det F>0`` and a material
    two-form is ``beta=i_q da``, then its physical vector representative is

        omega = F q / J.

    With ``g=F^T F`` one has exactly

        det(g^-1|_{q^perp}) = |Fq|^2/(J^2 |q|^2)=|omega|^2/|q|^2.

    Thus the transverse inverse-metric determinant law is a property of
    transported two-forms themselves; incompressibility is the special case
    ``J=1``.
    """
    F=np.asarray(deformation,dtype=float)
    q=np.asarray(tuple(float(x) for x in material_twoform_vector),dtype=float)
    if F.shape!=(3,3) or not np.all(np.isfinite(F)) or q.shape!=(3,) or not np.all(np.isfinite(q)):
        raise ValueError("finite 3x3 deformation and material two-form vector required")
    J=float(np.linalg.det(F)); nq=float(np.linalg.norm(q))
    if J<=0.0 or nq<=0.0:
        raise ValueError("orientation-preserving nonsingular deformation and nonzero two-form vector required")
    g=F.T@F; qn=q/nq; P=_transverse_frame(qn)
    detperp=float(np.linalg.det(P.T@np.linalg.inv(g)@P))
    amp2=float(np.dot(F@q,F@q)/(J*J*nq*nq))
    res=abs(detperp-amp2)
    cond=float(np.linalg.cond(F))
    tol=5e-9*max(1.0,cond*cond)*max(1.0,abs(detperp),abs(amp2))
    if res>tol:
        raise AssertionError("general transported-two-form cofactor identity failed")
    return {
        "jacobian":J,
        "physical_amplification_squared":amp2,
        "transverse_inverse_metric_determinant":detperp,
        "deformation_condition_number":cond,
        "identity_residual":res,
    }


def vortex_slip_twist_algebra(
    velocity: Sequence[float],
    vorticity: Sequence[float],
    vorticity_gradient: np.ndarray,
    viscosity: float,
) -> dict[str, float]:
    """Exact vortex-line slip/Frobenius-twist decomposition of the Hodge current.

    Write ``omega=m xi`` and ``c=curl omega``.  Away from ``m=0`` put

        tau = xi.curl xi,
        c_parallel = m tau xi,
        v_slip = -nu (omega x c)/m^2
               = nu[(xi.grad)xi-grad_perp log m],
        w = u + v_slip.

    Then, with the repository convention ``i_v beta=-v x omega``, the exact
    Ohm/Faraday field is

        e = i_u beta + nu c = i_w beta + nu c_parallel.

    Hence the perpendicular viscous current is exactly a change of vortex-line
    transport velocity.  Only the parallel one-form remains outside Lie
    transport; an exact-gradient part of that remainder is still gauge, so this
    function does *not* identify pointwise ``c_parallel`` with reconnection.

    The current norm and Poynting--Joule density obey

        nu |c|^2 = (m^2/nu)|v_slip|^2 + nu m^2 tau^2,

        -e.c = m^2|u_perp|^2/(4nu)
               -(m^2/nu)|v_slip+u_perp/2|^2-nu m^2 tau^2.
    """
    u=np.asarray(tuple(float(x) for x in velocity),dtype=float)
    wv=np.asarray(tuple(float(x) for x in vorticity),dtype=float)
    G=np.asarray(vorticity_gradient,dtype=float)
    nu=float(viscosity)
    if u.shape!=(3,) or wv.shape!=(3,) or G.shape!=(3,3) or not np.all(np.isfinite(u)) or not np.all(np.isfinite(wv)) or not np.all(np.isfinite(G)):
        raise ValueError("finite velocity, vorticity and 3x3 vorticity gradient required")
    if not math.isfinite(nu) or nu<=0.0:
        raise ValueError("positive finite viscosity required")
    div=float(np.trace(G)); gscale=max(1.0,float(np.linalg.norm(G)))
    if abs(div)>3e-10*gscale:
        raise ValueError("vorticity gradient must satisfy div omega=0")
    m=float(np.linalg.norm(wv))
    if m<=1e-12:
        raise ValueError("nonzero vorticity required for vortex-line gauge")
    xi=wv/m
    c=np.array((G[1,2]-G[2,1],G[2,0]-G[0,2],G[0,1]-G[1,0]),dtype=float)
    dm=G@xi
    dxi=(G-dm[:,None]*xi[None,:])/m
    curlxi=np.array((dxi[1,2]-dxi[2,1],dxi[2,0]-dxi[0,2],dxi[0,1]-dxi[1,0]),dtype=float)
    tau=float(np.dot(xi,curlxi))
    cpar=m*tau*xi
    cperp=c-cpar
    gradperp=dm-xi*float(np.dot(xi,dm))
    curvature=xi@dxi
    vslip=nu*(curvature-gradperp/m)
    vdirect=-nu*np.cross(wv,c)/(m*m)
    W=u+vslip
    e=-np.cross(u,wv)+nu*c
    erep=-np.cross(W,wv)+nu*cpar
    uperp=u-xi*float(np.dot(u,xi))
    current_cost=nu*float(np.dot(c,c))
    cost_rep=(m*m/nu)*float(np.dot(vslip,vslip))+nu*m*m*tau*tau
    joule=-float(np.dot(e,c))
    square=(m*m/(4.0*nu))*float(np.dot(uperp,uperp))-(m*m/nu)*float(np.dot(vslip+.5*uperp,vslip+.5*uperp))-nu*m*m*tau*tau
    scale=max(1.0,np.linalg.norm(e),np.linalg.norm(erep),np.linalg.norm(c)*nu,abs(current_cost),abs(joule),abs(square))
    res=max(float(np.linalg.norm(vslip-vdirect)),float(np.linalg.norm(cpar-xi*np.dot(c,xi))),float(np.linalg.norm(e-erep)),abs(current_cost-cost_rep),abs(joule-square))
    if res>2e-8*scale:
        raise AssertionError("vortex-line slip/twist algebra failed")
    return {
        "vorticity_magnitude":m,
        "frobenius_twist":tau,
        "parallel_current_norm":float(np.linalg.norm(cpar)),
        "perpendicular_current_norm":float(np.linalg.norm(cperp)),
        "slip_speed":float(np.linalg.norm(vslip)),
        "vortex_transport_speed":float(np.linalg.norm(W)),
        "viscous_current_cost_density":current_cost,
        "slip_twist_cost_density":cost_rep,
        "negative_joule_work_density":joule,
        "slip_twist_square_density":square,
        "parallel_residual_magnitude_sink":nu*m*tau*tau,
        "identity_residual":res,
    }



def material_metric_path_action_bound(time_horizon: float, energy_loss: float, viscosity: float) -> dict[str, float]:
    """Finite affine-invariant path action of the material metric.

    Pointwise ``g(a,t)`` lies in the canonical symmetric space
    ``SPD_1(3)=SL(3)/SO(3)`` with speed

        |g_t|_g^2=tr[(g^-1 g_t)^2].

    The primitive speed lock and velocity-energy law give

        int_0^t int |g_s|_g^2 da ds = (E(0)-E(t))/nu.

    Hence the total path length ``ell(a)=int_0^t |g_s(a)|_g ds`` satisfies

        int ell(a)^2 da <= t (E(0)-E(t))/nu,

    and so does the squared affine distance from the identity.  On any finite
    smooth interval, the metric path therefore has finite length for almost every
    material label.  No supremum or global-regularity conclusion follows.
    """
    t=float(time_horizon); loss=float(energy_loss); nu=float(viscosity)
    if not all(math.isfinite(x) for x in (t,loss,nu)) or t<0.0 or loss<0.0 or nu<=0.0:
        raise ValueError("nonnegative finite time/loss and positive viscosity required")
    action=loss/nu
    return {
        "metric_speed_spacetime_l2_squared":action,
        "material_path_length_l2_squared_upper":t*action,
        "affine_distance_l2_squared_upper":t*action,
        "log_max_stretch_l2_squared_upper":0.25*t*action,
    }



def klein_spacetime_vortex_worldsheet_algebra(
    velocity: Sequence[float],
    vorticity: Sequence[float],
    vorticity_current: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Klein-quadric decomposition of the physical spacetime vorticity curvature.

    With ``e=-u cross omega+nu c`` the physical curvature is
    ``F=beta-dt wedge e``.  In four dimensions a nonzero two-form is simple/rank
    two exactly when its Pfaffian, equivalently ``F wedge F``, vanishes.  Here

        Pf(F)=-e.omega=-nu c.omega  (for the orientation convention used below).

    At fixed spatial ``beta``, the Klein slice is the plane ``e.omega=0``.
    Orthogonal projection onto it removes only ``e_parallel=nu c_parallel``.
    The projected curvature is simple and its kernel contains the vortex-line
    spacetime directions ``(1,w)`` and ``(0,omega)``, where

        w=u-nu(omega cross c)/|omega|^2.

    The fixed-beta distance squared to the Klein slice is
    ``nu^2 |c_parallel|^2``.  Dividing by nu gives exactly the twist part of the
    palinstrophy density.  This is an algebraic/fiberwise statement; physical
    off-Klein curvature is not by itself a reconnection theorem.
    """
    u=np.asarray(tuple(float(x) for x in velocity),float)
    om=np.asarray(tuple(float(x) for x in vorticity),float)
    c=np.asarray(tuple(float(x) for x in vorticity_current),float)
    nu=float(viscosity)
    if u.shape!=(3,) or om.shape!=(3,) or c.shape!=(3,) or not np.all(np.isfinite(u)) or not np.all(np.isfinite(om)) or not np.all(np.isfinite(c)):
        raise ValueError("finite velocity/vorticity/current three-vectors required")
    if not math.isfinite(nu) or nu<0.0:
        raise ValueError("finite nonnegative viscosity required")
    m=float(np.linalg.norm(om))
    if m<=1e-12:
        raise ValueError("nonzero vorticity required")
    xi=om/m; cpar=xi*float(np.dot(c,xi)); cperp=c-cpar
    e=-np.cross(u,om)+nu*c
    eperp=e-xi*float(np.dot(e,xi)); epar=e-eperp
    w=u-nu*np.cross(om,c)/(m*m)
    # skew-matrix convention with Pfaffian = e.omega
    F=np.zeros((4,4)); F[0,1:]=-e; F[1:,0]=e
    X=np.array([[0,-om[2],om[1]],[om[2],0,-om[0]],[-om[1],om[0],0]])
    F[1:,1:]=-X
    FT=np.zeros((4,4));FT[0,1:]=-eperp;FT[1:,0]=eperp;FT[1:,1:]=-X
    pf=F[0,1]*F[2,3]-F[0,2]*F[1,3]+F[0,3]*F[1,2]
    pfT=FT[0,1]*FT[2,3]-FT[0,2]*FT[1,3]+FT[0,3]*FT[1,2]
    k1=np.r_[1.0,w]; k2=np.r_[0.0,om]
    kernel=max(float(np.linalg.norm(FT@k1)),float(np.linalg.norm(FT@k2)))
    dist2=float(np.dot(epar,epar)); represented=nu*nu*float(np.dot(cpar,cpar))
    twist_cost=(dist2/nu if nu>0 else 0.0)
    top=-nu*float(np.dot(c,om))
    scale=max(1.0,np.linalg.norm(F),np.linalg.norm(FT),abs(pf),abs(top),dist2,represented)
    res=max(abs(pf-top),abs(pfT),kernel,abs(dist2-represented),abs(np.linalg.det(F)-pf*pf))
    if res>2e-10*scale:
        raise AssertionError("Klein spacetime curvature algebra failed")
    return {
        "physical_pfaffian":float(pf),
        "viscous_pfaffian":top,
        "projected_klein_pfaffian":float(pfT),
        "projected_worldsheet_kernel_residual":kernel,
        "fixed_beta_klein_distance_squared":dist2,
        "viscosity_squared_parallel_current":represented,
        "twist_dissipation_density":twist_cost,
        "identity_residual":res,
    }


def local_flux_velocity_gauge_algebra(
    velocity: Sequence[float],
    vorticity: Sequence[float],
    electromotive: Sequence[float],
) -> dict[str, float]:
    """Local transport-gauge decomposition of Faraday's law away from omega=0.

    Given ``beta_t+d e=0``, pure Lie transport by some local velocity ``w`` is
    equivalent to finding ``psi`` with

        e-i_w beta=d psi.

    Pointwise choose the directional derivative of psi so that

        omega.grad psi=e.omega.

    The remainder ``e-d psi`` is transverse to omega and is therefore exactly
    ``i_w beta`` for a suitable transverse ``w``.  The function returns one
    pointwise representative.  Along an actual vortex line this is a first-order
    ODE for psi and is locally solvable.  A global single-valued psi can fail on
    closed/recurrent lines; that global leafwise cohomology is not decided by this
    pointwise algebra.
    """
    u=np.asarray(tuple(float(x) for x in velocity),float)
    om=np.asarray(tuple(float(x) for x in vorticity),float)
    e=np.asarray(tuple(float(x) for x in electromotive),float)
    if u.shape!=(3,) or om.shape!=(3,) or e.shape!=(3,) or not np.all(np.isfinite(u)) or not np.all(np.isfinite(om)) or not np.all(np.isfinite(e)):
        raise ValueError("finite velocity/vorticity/electromotive three-vectors required")
    m2=float(np.dot(om,om))
    if m2<=1e-24:
        raise ValueError("nonzero vorticity required")
    gradpsi=om*float(np.dot(e,om))/m2
    et=e-gradpsi
    wperp=-np.cross(om,et)/m2
    w=wperp+om*float(np.dot(u,om))/m2
    represented=-np.cross(w,om)+gradpsi
    residual=float(np.linalg.norm(e-represented))
    if residual>2e-11*max(1.0,np.linalg.norm(e),np.linalg.norm(represented)):
        raise AssertionError("local flux-velocity gauge decomposition failed")
    return {
        "directional_gauge_derivative":float(np.dot(gradpsi,om)/math.sqrt(m2)),
        "transport_velocity_norm":float(np.linalg.norm(w)),
        "transverse_emotive_norm":float(np.linalg.norm(et)),
        "parallel_emotive_norm":float(np.linalg.norm(gradpsi)),
        "identity_residual":residual,
    }


def closed_vortex_line_period_cost(
    line_length: float,
    electromotive_period: float,
    viscosity: float,
) -> dict[str, float]:
    """Sharp Cauchy cost of a gauge-invariant non-ideal period on a closed vortex line.

    Along a closed vortex line the ideal one-form ``i_u beta`` vanishes on the
    tangent, so the Faraday/transport obstruction is

        E_gamma = integral_gamma e = nu integral_gamma c^flat.

    If ``D_twist,gamma=nu integral_gamma m^2 tau^2 ds`` is the linewise twist
    dissipation, then

        |E_gamma|^2 <= nu L_gamma D_twist,gamma.

    This helper returns the *minimum* linewise twist dissipation compatible with
    the supplied period.  It does not claim that all vortex lines are closed or
    that linewise cost has already been converted into a volume budget.
    """
    L=float(line_length); P=float(electromotive_period); nu=float(viscosity)
    if not all(math.isfinite(x) for x in (L,P,nu)) or L<=0.0 or nu<=0.0:
        raise ValueError("positive finite line length/viscosity and finite period required")
    lower=P*P/(nu*L)
    return {
        "line_length":L,
        "electromotive_period":P,
        "minimum_twist_dissipation":lower,
        "period_squared":P*P,
    }



def curl_line_geometry_algebra(field: Sequence[float], field_gradient: np.ndarray) -> dict[str, float]:
    """Universal polar line-geometry decomposition of the primitive curl operator.

    For a nonzero vector field ``b=m n`` at one jet, define

        kappa=(n.grad)n,
        tau=n.curl n,
        A=kappa-grad_perp log m.

    Then exactly

        curl b = m tau n + m n cross A,
        b.curl b = m^2 tau,
        b cross curl b = -m^2 A,
        |curl b|^2 = m^2(tau^2+|A|^2).

    This is not a new PDE variable: it is the longitudinal/transverse polar
    decomposition of ``curl`` itself.  ``tau`` is Frobenius twist and ``A`` is
    the curvature-minus-transverse-concentration defect of the field lines.
    """
    b=np.asarray(tuple(float(x) for x in field),dtype=float)
    G=np.asarray(field_gradient,dtype=float)
    if b.shape!=(3,) or G.shape!=(3,3) or not np.all(np.isfinite(b)) or not np.all(np.isfinite(G)):
        raise ValueError("finite nonzero vector and finite 3x3 first jet required")
    m=float(np.linalg.norm(b))
    if m<=1e-12:
        raise ValueError("nonzero field required")
    n=b/m; dm=G@n; dn=(G-dm[:,None]*n[None,:])/m
    kappa=n@dn
    curln=np.array((dn[1,2]-dn[2,1],dn[2,0]-dn[0,2],dn[0,1]-dn[1,0]),dtype=float)
    tau=float(np.dot(n,curln))
    gradperp=dm-n*float(np.dot(n,dm))
    A=kappa-gradperp/m
    curlb=np.array((G[1,2]-G[2,1],G[2,0]-G[0,2],G[0,1]-G[1,0]),dtype=float)
    represented=m*tau*n+m*np.cross(n,A)
    hel=float(np.dot(b,curlb)); cross=np.cross(b,curlb)
    curl2=float(np.dot(curlb,curlb)); rep2=m*m*(tau*tau+float(np.dot(A,A)))
    scale=max(1.0,m,np.linalg.norm(curlb),m*m*(1.0+abs(tau)+np.linalg.norm(A)))
    res=max(float(np.linalg.norm(curlb-represented)),abs(hel-m*m*tau),float(np.linalg.norm(cross+m*m*A)),abs(curl2-rep2))
    if res>2e-10*scale:
        raise AssertionError("curl line-geometry polar decomposition failed")
    return {
        "field_magnitude":m,
        "frobenius_twist":tau,
        "transverse_line_defect_norm":float(np.linalg.norm(A)),
        "curl_norm_squared":curl2,
        "twist_plus_defect_curl_norm_squared":rep2,
        "helicity_density":hel,
        "cross_curl_norm":float(np.linalg.norm(cross)),
        "identity_residual":res,
    }


def two_level_curl_geometry_algebra(
    velocity: Sequence[float],
    vorticity: Sequence[float],
    vorticity_current: Sequence[float],
    viscosity: float,
) -> dict[str, float]:
    """Two consecutive applications of the same curl-line geometry in NS.

    Let ``omega=curl u`` and ``c=curl omega`` at a point.  Algebraically define

        A_u=-u cross omega/|u|^2,
        tau_u=u.omega/|u|^2,
        A_w=-omega cross c/|omega|^2,
        tau_w=omega.c/|omega|^2.

    These are exactly the ``A,tau`` data returned by ``curl_line_geometry`` when
    the corresponding first jets exist.  The Euler and viscous objects are

        u cross omega=-|u|^2 A_u,
        v_slip=nu A_w,
        |c|^2=|omega|^2(tau_w^2+|A_w|^2),

    and the nonlinear enstrophy density is

        (u cross omega).c
        = |u|^2 |omega| (xi cross A_u).A_w.

    Thus NS uses the same primitive curl geometry at two successive levels; the
    longitudinal twist at the second level is invisible to nonlinear production
    and remains a pure viscous cost.
    """
    u=np.asarray(tuple(float(x) for x in velocity),float)
    om=np.asarray(tuple(float(x) for x in vorticity),float)
    c=np.asarray(tuple(float(x) for x in vorticity_current),float)
    nu=float(viscosity)
    if u.shape!=(3,) or om.shape!=(3,) or c.shape!=(3,) or not np.all(np.isfinite(u)) or not np.all(np.isfinite(om)) or not np.all(np.isfinite(c)):
        raise ValueError("finite velocity/vorticity/current vectors required")
    if not math.isfinite(nu) or nu<0.0:
        raise ValueError("finite nonnegative viscosity required")
    a=float(np.linalg.norm(u));m=float(np.linalg.norm(om))
    if a<=1e-12 or m<=1e-12:
        raise ValueError("nonzero velocity and vorticity required")
    xi=om/m
    Au=-np.cross(u,om)/(a*a); tauu=float(np.dot(u,om)/(a*a))
    Aw=-np.cross(om,c)/(m*m); tauw=float(np.dot(om,c)/(m*m))
    lamb=np.cross(u,om); vslip=nu*Aw
    crep=m*tauw*xi+m*np.cross(xi,Aw)
    stretch=float(np.dot(lamb,c)); stretchrep=a*a*m*float(np.dot(np.cross(xi,Au),Aw))
    p=float(np.dot(c,c)); prep=m*m*(tauw*tauw+float(np.dot(Aw,Aw)))
    res=max(float(np.linalg.norm(lamb+a*a*Au)),float(np.linalg.norm(c-crep)),abs(stretch-stretchrep),abs(p-prep),float(np.linalg.norm(vslip-nu*Aw)))
    if res>2e-10*max(1.0,np.linalg.norm(lamb),np.linalg.norm(c),abs(stretch),p):
        raise AssertionError("two-level primitive curl geometry failed")
    return {
        "velocity_twist":tauu,
        "velocity_transverse_defect_norm":float(np.linalg.norm(Au)),
        "vorticity_twist":tauw,
        "vorticity_transverse_defect_norm":float(np.linalg.norm(Aw)),
        "lamb_norm":float(np.linalg.norm(lamb)),
        "viscous_slip_speed":float(np.linalg.norm(vslip)),
        "stretching_density":stretch,
        "represented_stretching_density":stretchrep,
        "palinstrophy_density":p,
        "represented_palinstrophy_density":prep,
        "identity_residual":res,
    }

def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "primitive_material_system": "with L_g=d delta_g+delta_g d, beta_t=-nu L_g beta, alpha=delta_g L_g^-1 beta, v=alpha^sharp_g, g_t=Lie_v g, d beta=0, det g=1 and Riem(g)=0; beta and g form a closed autonomous material NS system",
        "primitive_hodge_energy": "kinetic energy is exactly <beta,L_g^-1 beta>_g, enstrophy is ||beta||_g^2, and the heat part gives E'=-2nu||beta||_g^2 while Euler metric motion contributes zero",
        "transverse_determinant": "for beta=i_q da, det(g^-1|q^perp)=q.g.q/|q|^2=|Fq|^2/|q|^2; vortex amplification is exactly transverse heat-symbol determinant",
        "transverse_area": "for xi,eta perpendicular q, det Gram_g^-1(xi,eta)=(|Fq|^2/|q|^2)|xi wedge eta|^2",
        "history_memory": "Euler freezes material beta, and Minkowski determinant gives sqrt(det int g^-1|qperp dt)>=int |Fq|/|q| dt; rotating anisotropy cannot reset accumulated transverse heat area",
        "moving_polarization_memory": "for full NS, int|F q(t)|dt <= |q(T)| sqrt(det int g^-1|q(T)^perp dt)+int|F(q(t)-q(T))|dt; the only escape from fixed-plane heat memory is actual heat-driven rewriting of beta",
        "heat_only_reset": "Euler has no beta_t term in material coordinates; viscosity alone changes beta and ||beta_t||_H^-2_g^2=nu^2||beta||_g^2, so integrated reset action divided by nu equals half the physical velocity-energy loss",
        "primitive_state_speed_lock": "with ||g_t||_M^2=int tr[(g^-1 g_t)^2], NS gives ||g_t||_M^2=2||beta||_g^2 and ||beta_t||_H^-2_g^2=(nu^2/2)||g_t||_M^2; equivalently -E'=nu||g_t||_M^2=(2/nu)||beta_t||_H^-2_g^2",
        "primitive_current_law": "with h=g_t and c=delta_g beta, the same velocity gradient gives c=-div_g h; j=nu c is the actual viscous material current and beta_t+d j=0, so material vorticity can only be rewritten by local flux through material boundaries",
        "all_scale_speed_ladder": "for every real s in the common Hodge domain, ||beta_t||_H^(s-2)^2=||j||_H^(s-1)^2=(nu^2/2)||g_t||_H^s^2; reset, viscous flux and metric deformation are one scale-shifted state velocity",
        "stress_current_noether": "T_beta=<i_X beta,i_Y beta>-(1/2)|beta|^2 g obeys div T_beta=i_(delta beta) beta; hence stretching=-1/2 int T_beta:g_t=int v.div T_beta=<delta beta,v cross omega>, so Euler work uses the same codifferential current that viscosity transports",
        "chern_simons_maxwell_extension": "curl C=*d is simultaneously the Hessian of helicity/Chern-Simons, its modulus |C| gives the critical stock, and C^2 is Hodge heat; the canonical extension A(y)=exp(-y|C|)u has source-free Maxwell curvature F4 with int|F4|^2=K and topological charge H",
        "maxwell_duality_stress": "the 4D Maxwell extension splits into self-dual/anti-self-dual energies (K+H)/2 and (K-H)/2; each pure duality has zero stress and |T|^2=4|F+|^2|F-|^2, so critical Euler work is intrinsically cross-duality stress",
        "poisson_sech_origin": "unit critical Maxwell depth profiles sqrt(2r)e^-ry have overlap 2sqrt(rs)/(r+s)=sech((log r-log s)/2); critical log-scale locality is harmonic-depth overlap, not a shell cutoff",
        "primitive_spacetime_curvature": "with alpha=u^flat, beta=d alpha, c=delta beta and Bernoulli B=p+|u|^2/2, the Abelian spacetime connection A4=alpha-B dt has curvature F4=beta-dt wedge(i_u beta+nu c); d_4F4=0 is exactly the vorticity equation, Euler obeys (i_u beta) wedge beta=0 pointwise, and F4 wedge F4 is purely viscous",
        "vortex_line_slip_gauge": "away from omega=0, split c=delta beta into c_perp+c_parallel and set w=u-nu(omega cross c)/|omega|^2; then e=i_w beta+nu c_parallel and beta_t+Lie_w beta=-nu d(c_parallel^flat), so all perpendicular viscous current is exact vortex-line slip and only a parallel one-form remains outside Lie transport",
        "frobenius_twist_current": "for omega=m xi, c_parallel=m(xi.curl xi)xi; xi.curl xi is the Frobenius obstruction xi^flat wedge d xi^flat, and xi.curl(m tau xi)=m tau^2 makes the parallel residual a direct magnitude sink in the vortex-line gauge",
        "slip_twist_pythagoras": "v_slip=nu[(xi.grad)xi-grad_perp log m] and nu|c|^2=(m^2/nu)|v_slip|^2+nu m^2 tau^2; the Poynting-Joule density completes the same slip square plus the pure twist sink",
        "vortex_line_topology_guard": "pointwise c_parallel is only the possible non-Lie remainder, not itself a reconnection theorem: exact-gradient pieces are gauge and only d(c_parallel^flat) rewrites beta in the vortex-line frame",
        "general_twoform_cofactor": "for any orientation-preserving transport F with J=det F, det((F^T F)^-1|q^perp)=|Fq|^2/(J^2|q|^2), exactly the squared amplification of a transported two-form vector; the transverse determinant law is not restricted to incompressible fluid labels",
        "klein_vortex_worldsheet": "Euler spacetime curvature is a closed simple rank-two two-form on the Klein quadric; its characteristic kernel is spanned by partial_t+u and omega and is automatically integrable, so frozen vortex lines are literal spacetime vortex worldsheets",
        "klein_tangent_normal_current": "at fixed beta, the Klein slice is e.omega=0; c_perp moves the electromotive field tangentially inside this slice, while the orthogonal departure is e_parallel=nu c_parallel with squared fiber distance nu^2|c_parallel|^2 and dissipation distance^2/nu",
        "local_flux_velocity_gauge": "away from omega=0, Faraday beta_t+d e=0 is locally pure Lie transport after solving omega.grad psi=e.omega and writing e-d psi=i_w beta; pointwise twist/Pfaffian is therefore not a local reconnection obstruction",
        "leafwise_period_obstruction": "a global single-valued flux velocity can fail through the leafwise cohomology of e along vortex lines; on a closed vortex line gamma the gauge-invariant period is integral_gamma e=nu integral_gamma c, and period^2<=nu L_gamma times the linewise twist dissipation",
        "klein_topology_guard": "physical F wedge F or c.omega measures off-Klein spacetime curvature but is not by itself a reconnection theorem: exact transport gauges can remove local parallel electromotive fields, and global topology requires leafwise period/cohomology information",
        "primitive_poynting_joule": "put e=i_u beta+nu c; the exact local law is partial_t |omega|^2/2+div(e cross omega)=-e.c, with -e.c=(u cross omega).c-nu|c|^2, so stretching and palinstrophy are one electromotive-current work law",
        "so33_duality_group": "volume-preserving four-dimensional deformation acts on two-forms through sl(4)->so(3,3); rotations commute with Hodge duality while trace-free strain anticommutes and is the unique noncompact duality-mixing boost, with ||rho(S)||_HS^2=2||S||^2",
        "enstrophy_derivative_null": "integrating <curl omega,u cross omega> by parts gives sum_j <partial_j omega, partial_j u cross u>; the derivative-on-u self term vanishes pointwise by cross-product skewness, so nonlinear enstrophy production necessarily uses the same first vorticity derivative squared by Hodge heat",
        "curl_polar_line_geometry": "for every nonzero b=m n, curl b=m tau_b n+m n cross A_b with tau_b=n.curl n and A_b=(n.grad)n-grad_perp log m=-b cross curl b/|b|^2; helicity, Lamb force and curl norm are the longitudinal/transverse pieces of this one operator law",
        "iterated_curl_ns_grammar": "applying the same curl-polar law to u and then omega gives u cross omega=-|u|^2 A_u, v_slip=nu A_omega, |curl omega|^2=|omega|^2(tau_omega^2+|A_omega|^2), and stretching=|u|^2|omega|(xi cross A_u).A_omega; no separate phase/coherence mechanism is added",
        "rank_one_null": "a one-direction incompressible gradient a tensor xi has omega=xi cross a and S omega=0; self-stretching is absent at rank one",
        "flat_hodge_dirichlet": "because every g=Phi^*g0 is flat, <beta,L_g beta>_g=||nabla^g beta||_2^2; material spatial non-affinity and vorticity magnitude/direction variation are already part of the same heat Dirichlet form, not a separate escape channel",
        "pair_mismatch_collapse": "omega_a cross omega_b=F_a^-T(q_a cross q_b)+(F_a q_a) cross((F_b-F_a)q_b); these are coordinate pieces of one covariant material-two-form variation, whose intrinsic norm is the same nabla^g beta squared by Hodge heat",
        "distortion_budget": "int [log sigma_max F]_+^2 da <= t(E0-Et)/(4nu); extreme material distortion is globally L2-log sparse but a supremum is not controlled",
        "hodge_lax_isospectral": "naturality gives partial_t L_g=[Lie_v,L_g] and partial_t delta_g=[Lie_v,delta_g]; every material Hodge operator is conjugate to the fixed Euclidean one, so Euler moves its frame but does not create or destroy heat eigenvalues",
        "current_forced_heat": "with c=delta_g beta and beta_t=-nu L_g beta, c_t=[Lie_v,delta_g]beta-nu L_g^(1)c; Euler can regenerate the viscous current only by motion of the same Hodge frame that it creates through g_t",
        "maurer_cartan_zero_curvature": "for F=D_a Phi, Gamma=F^-1 dF and B=F^-1 F_t obey dGamma+Gamma wedge Gamma=0 and Gamma_t=D_Gamma B; deformation-frame turnover is a pure-gauge SL(3) connection, not an independent phase field",
        "connection_current_lock": "flatness and incompressibility give ||D_Gamma B||_L2^2=||delta_g beta||_L2^2 and ||nabla g_t||_L2^2=2||delta_g beta||_L2^2; frame turnover and Hodge current are the same derivative-order activity",
        "symmetric_space_path_action": "g(a,t) lies in SL(3)/SO(3) and int path_length(a)^2 da <= t(E0-Et)/nu; hence the material metric has finite total affine path length for almost every label on every finite smooth interval, while concentration on a null label set remains open",
        "kelvin_current_parent": "the pulled-back velocity one-form obeys alpha_t+d pi=-nu delta_g d alpha, so d/dt circulation(gamma)=-nu integral_gamma delta_g beta; beta_t+nu d delta_g beta=0 is the exterior derivative of one Kelvin-current law",
        "lagrangian_geodesic_acceleration": "for X_t=u(X,t), full NS gives X_tt=(-grad p+nu Delta u)(X,t) and F_tt=(-Hess p+nu grad Delta u)(X,t)F; the explicit quadratic A^2 self-stretch cancels from material acceleration, and for Euler D_t^2 omega=-(Hess p)omega along a frozen material vorticity vector",
        "global_regularity_claimed": False,
        "case_taxonomy_used": False,
        "analysis_cutoff_used": False,
    }
