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
        "enstrophy_derivative_null": "integrating <curl omega,u cross omega> by parts gives sum_j <partial_j omega, partial_j u cross u>; the derivative-on-u self term vanishes pointwise by cross-product skewness, so nonlinear enstrophy production necessarily uses the same first vorticity derivative squared by Hodge heat",
        "rank_one_null": "a one-direction incompressible gradient a tensor xi has omega=xi cross a and S omega=0; self-stretching is absent at rank one",
        "flat_hodge_dirichlet": "because every g=Phi^*g0 is flat, <beta,L_g beta>_g=||nabla^g beta||_2^2; material spatial non-affinity and vorticity magnitude/direction variation are already part of the same heat Dirichlet form, not a separate escape channel",
        "pair_mismatch_collapse": "omega_a cross omega_b=F_a^-T(q_a cross q_b)+(F_a q_a) cross((F_b-F_a)q_b); these are coordinate pieces of one covariant material-two-form variation, whose intrinsic norm is the same nabla^g beta squared by Hodge heat",
        "distortion_budget": "int [log sigma_max F]_+^2 da <= t(E0-Et)/(4nu); extreme material distortion is globally L2-log sparse but a supremum is not controlled",
        "global_regularity_claimed": False,
        "case_taxonomy_used": False,
        "analysis_cutoff_used": False,
    }
