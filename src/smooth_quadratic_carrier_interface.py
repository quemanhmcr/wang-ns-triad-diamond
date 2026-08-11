from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.event_anchored_role_registration import envelope_registration_residual
from src.nonaffine_role_interface_work import adjoint_split
from src.outer_moving_role_extraction import bilinear_apply, linearized_resolved
from src.physical_energy_causal_bridge import route_physical_energy_causality


STATUS = (
    "EXACT_SMOOTH_QUADRATIC_CARRIER_INTERFACE__"
    "Q2_ENERGY_LAW__COMMON_GAUGE_QUOTIENT_BEFORE_PHYSICAL_RELINK__"
    "SYMMETRIC_WORK_EXISTING_STRAIN__COEFFICIENT_OBSTRUCTION_ENERGY_REENTRY"
)

RELINK_OWNER = "conservative_smooth_role_relink"
STRAIN_OWNER = "existing_resolved_strain_deformation"


def _matching_square_operators(operators: Sequence[np.ndarray], *, name: str) -> list[np.ndarray]:
    out = [np.asarray(A, dtype=complex) for A in operators]
    if not out or any(A.ndim != 2 or A.shape[0] != A.shape[1] for A in out):
        raise ValueError(f"nonempty square {name} required")
    n = out[0].shape[0]
    if any(A.shape != (n, n) for A in out):
        raise ValueError(f"matching square {name} required")
    return out


def _selfadjoint_positive_contraction_residual(A: np.ndarray) -> tuple[float, float, float]:
    selfadjoint = float(np.linalg.norm(A - A.conj().T))
    H = 0.5 * (A + A.conj().T)
    eig = np.linalg.eigvalsh(H)
    lower = max(0.0, -float(np.min(eig)))
    upper = max(0.0, float(np.max(eig)) - 1.0)
    return selfadjoint, lower, upper


def quadratic_partition_diagnostics(
    analysis_operators: Sequence[np.ndarray],
    analysis_rates: Sequence[np.ndarray] | None = None,
) -> dict[str, float]:
    """Validate the natural smooth energy partition sum_a A_a^2=I.

    The propagated PDE objects are positive self-adjoint analysis multipliers,
    not hard projectors.  Their channel energies are ||A_a u||_2^2 and hence
    their physical effects are eta_a=A_a^2.  If rates are supplied, the exact
    differentiated partition is also checked:

        sum_a (dot A_a A_a + A_a dot A_a)=0.
    """
    As = _matching_square_operators(analysis_operators, name="analysis operators")
    n = As[0].shape[0]
    worst_selfadjoint = worst_lower = worst_upper = 0.0
    for A in As:
        sa, lo, hi = _selfadjoint_positive_contraction_residual(A)
        worst_selfadjoint = max(worst_selfadjoint, sa)
        worst_lower = max(worst_lower, lo)
        worst_upper = max(worst_upper, hi)

    total = sum((A @ A for A in As), np.zeros((n, n), dtype=complex))
    partition = float(np.linalg.norm(total - np.eye(n)))
    rate_selfadjoint = derivative = 0.0
    if analysis_rates is not None:
        dAs = _matching_square_operators(analysis_rates, name="analysis rates")
        if len(dAs) != len(As) or dAs[0].shape != As[0].shape:
            raise ValueError("one analysis rate per smooth carrier required")
        derivative_matrix = np.zeros((n, n), dtype=complex)
        for A, dA in zip(As, dAs, strict=True):
            rate_selfadjoint = max(rate_selfadjoint, float(np.linalg.norm(dA - dA.conj().T)))
            derivative_matrix += dA @ A + A @ dA
        derivative = float(np.linalg.norm(derivative_matrix))

    scale = max(1.0, float(n), float(len(As)))
    if max(worst_selfadjoint, worst_lower, worst_upper, partition) > 4e-10 * scale:
        raise ValueError("positive self-adjoint quadratic analysis partition required")
    if max(rate_selfadjoint, derivative) > 6e-10 * scale:
        raise ValueError("analysis rates must preserve the quadratic partition")
    return {
        "worst_selfadjoint_residual": worst_selfadjoint,
        "worst_positivity_defect": worst_lower,
        "worst_contraction_defect": worst_upper,
        "quadratic_partition_residual": partition,
        "worst_rate_selfadjoint_residual": rate_selfadjoint,
        "quadratic_partition_derivative_residual": derivative,
    }


def quadratic_complement(smooth_envelope: np.ndarray) -> np.ndarray:
    """Return R=(I-Q^2)^(1/2) for a positive self-adjoint contraction Q.

    Analytically the preferred construction is to choose Q=cos(theta) and
    R=sin(theta) from one smooth angle, which preserves smoothness at plateaux.
    This finite-dimensional functional-calculus helper certifies the same exact
    energy relation Q^2+R^2=I and is used for cross-module regression tests.
    """
    Q = np.asarray(smooth_envelope, dtype=complex)
    if Q.ndim != 2 or Q.shape[0] != Q.shape[1]:
        raise ValueError("square smooth envelope required")
    sa, lo, hi = _selfadjoint_positive_contraction_residual(Q)
    scale = max(1.0, float(np.linalg.norm(Q)))
    if max(sa, lo, hi) > 3e-10 * scale:
        raise ValueError("positive self-adjoint contraction required")
    H = 0.5 * (Q + Q.conj().T)
    eig, U = np.linalg.eigh(H)
    eig = np.clip(eig, 0.0, 1.0)
    R = (U * np.sqrt(np.maximum(0.0, 1.0 - eig * eig))) @ U.conj().T
    R = 0.5 * (R + R.conj().T)
    quadratic_partition_diagnostics((Q, R))
    return R


def _commutator_work(A: np.ndarray, L: np.ndarray, u: np.ndarray) -> float:
    w = A @ u
    return 2.0 * float(np.real(np.vdot(w, (L @ A - A @ L) @ u)))


def hard_linear_complement_skew_defect(
    smooth_envelope: np.ndarray,
    skew_operator: np.ndarray,
    state: np.ndarray,
) -> dict[str, float]:
    """Expose why I-Q is not the energy complement of a smooth envelope.

    For K*=-K and self-adjoint Q,

      I_Q(K)+I_(I-Q)(K)=4 Re <Q(I-Q)u,K u>,

    which need not vanish.  Replacing I-Q by the quadratic complement R with
    Q^2+R^2=I restores exact skew conservation.
    """
    Q = np.asarray(smooth_envelope, dtype=complex)
    K = np.asarray(skew_operator, dtype=complex)
    u = np.asarray(state, dtype=complex)
    n = len(u)
    if Q.shape != (n, n) or K.shape != (n, n):
        raise ValueError("matching smooth-envelope/skew/state data required")
    if float(np.linalg.norm(K + K.conj().T)) > 3e-10 * max(1.0, float(np.linalg.norm(K))):
        raise ValueError("skew-adjoint operator required")
    R = quadratic_complement(Q)
    linear = np.eye(n, dtype=complex) - Q
    defect = _commutator_work(Q, K, u) + _commutator_work(linear, K, u)
    overlap = 4.0 * float(np.real(np.vdot(Q @ linear @ u, K @ u)))
    quadratic = _commutator_work(Q, K, u) + _commutator_work(R, K, u)
    scale = max(1.0, float(np.linalg.norm(K)) * float(np.linalg.norm(u)) ** 2)
    if max(abs(defect - overlap), abs(quadratic)) > 2e-10 * scale:
        raise AssertionError("smooth-envelope complement identities failed")
    return {
        "linear_complement_skew_defect": defect,
        "overlap_defect_formula": overlap,
        "overlap_formula_residual": defect - overlap,
        "quadratic_complement_skew_residual": quadratic,
    }


@dataclass(frozen=True)
class GaugeQuotientedInterfaceWork:
    """Physical smooth-interface work after exact observer-gauge quotient.

    The signed atoms are physical-routing eligible only because the analysis roles
    have first been proved to move by the supplied common physical transport
    generator.  Arbitrary motion of an analysis partition cannot construct this
    certificate through :func:`smooth_quadratic_interface_balance`.
    """

    signed_native_interface_atoms: tuple[float, ...]
    signed_physical_relink_atoms: tuple[float, ...]
    signed_existing_strain_atoms: tuple[float, ...]
    gauge_transport_operator_residual: float
    skew_decomposition_residual: float
    signed_physical_relink_pair_matrix: tuple[tuple[float, ...], ...] = ()

    def __post_init__(self) -> None:
        residuals = (
            float(self.gauge_transport_operator_residual),
            float(self.skew_decomposition_residual),
        )
        if not all(math.isfinite(value) and value >= 0.0 for value in residuals):
            raise ValueError(
                "gauge-transport and skew-decomposition residuals must be finite and nonnegative"
            )


def smooth_quadratic_interface_balance(
    analysis_operators: Sequence[np.ndarray],
    analysis_rates: Sequence[np.ndarray],
    resolved_operator: np.ndarray,
    state: np.ndarray,
    *,
    common_gauge_skew_operator: np.ndarray,
) -> dict[str, object]:
    """Exact Q^2 interface law after quotienting common transported gauge motion.

    Let eta_a=A_a^2, sum eta_a=I and L=K+S.  A moving analysis partition is
    *not* automatically physical.  The supplied common gauge generator G must be
    skew-adjoint and must transport every actual analysis role by

        dot A_a + [G,A_a] = 0.

    This is the finite-dimensional form of transporting the smooth Fourier role
    by the common affine/Kelvin PDE motion.  Only after that exact identity is
    verified do we write

        K = G + K_phys.

    The dot-A work then cancels the G-work exactly.  The remaining K_phys work is
    conservative physical relinking, while S is the existing strain/deformation
    work.  Hence observer-chosen partition motion is quotiented before causal
    ownership is assigned.
    """
    As = _matching_square_operators(analysis_operators, name="analysis operators")
    dAs = _matching_square_operators(analysis_rates, name="analysis rates")
    if len(As) != len(dAs):
        raise ValueError("one analysis rate per smooth carrier required")
    quadratic_partition_diagnostics(As, dAs)
    u = np.asarray(state, dtype=complex)
    L = np.asarray(resolved_operator, dtype=complex)
    G = np.asarray(common_gauge_skew_operator, dtype=complex)
    n = As[0].shape[0]
    if u.shape != (n,) or L.shape != (n, n) or G.shape != (n, n):
        raise ValueError("matching resolved operator, common gauge generator and state required")

    gscale = max(1.0, float(np.linalg.norm(G)))
    if float(np.linalg.norm(G + G.conj().T)) > 3e-10 * gscale:
        raise ValueError("common gauge transport generator must be skew-adjoint")

    transport_residuals = [dA + G @ A - A @ G for A, dA in zip(As, dAs, strict=True)]
    transport_scale = max(
        1.0,
        max(float(np.linalg.norm(dA)) for dA in dAs),
        max(float(np.linalg.norm(G @ A - A @ G)) for A in As),
    )
    transport_residual = max(float(np.linalg.norm(R)) for R in transport_residuals)
    if transport_residual > 6e-10 * transport_scale:
        raise ValueError(
            "analysis partition motion is not the certified common affine/Kelvin transport; quotient observer motion before charging physics"
        )

    K, S = adjoint_split(L)
    Kphys = K - G
    skew_decomposition_residual = float(np.linalg.norm(K - G - Kphys))
    kscale = max(1.0, float(np.linalg.norm(K)), float(np.linalg.norm(G)), float(np.linalg.norm(Kphys)))
    if float(np.linalg.norm(Kphys + Kphys.conj().T)) > 4e-10 * kscale:
        raise AssertionError("gauge-quotiented residual skew operator lost skew-adjointness")

    etas = [A @ A for A in As]
    detas = [dA @ A + A @ dA for A, dA in zip(As, dAs, strict=True)]
    moving = np.array([float(np.real(np.vdot(u, deta @ u))) for deta in detas])
    gauge_skew = np.array(
        [-2.0 * float(np.real(np.vdot(eta @ u, G @ u))) for eta in etas]
    )
    relink = np.array(
        [-2.0 * float(np.real(np.vdot(eta @ u, Kphys @ u))) for eta in etas]
    )
    strain = np.array(
        [-2.0 * float(np.real(np.vdot(eta @ u, S @ u))) for eta in etas]
    )
    direct_native = np.array(
        [
            float(np.real(np.vdot(u, deta @ u)))
            - 2.0 * float(np.real(np.vdot(eta @ u, L @ u)))
            for eta, deta in zip(etas, detas, strict=True)
        ]
    )
    native = relink + strain

    outer_interface = np.empty(len(As), dtype=float)
    diagonal_resolved = np.empty(len(As), dtype=float)
    for a, (A, dA) in enumerate(zip(As, dAs, strict=True)):
        w = A @ u
        outer_interface[a] = (
            2.0 * float(np.real(np.vdot(w, dA @ u))) + _commutator_work(A, L, u)
        )
        diagonal_resolved[a] = 2.0 * float(np.real(np.vdot(w, L @ w)))

    m = len(As)
    TK = np.zeros((m, m), dtype=float)
    DS = np.zeros((m, m), dtype=float)
    for a in range(m):
        for b in range(m):
            TK[a, b] = -2.0 * float(np.real(np.vdot(etas[a] @ u, Kphys @ (etas[b] @ u))))
            DS[a, b] = -2.0 * float(np.real(np.vdot(etas[a] @ u, S @ (etas[b] @ u))))

    u2 = float(np.linalg.norm(u)) ** 2
    rate_scale = sum(
        float(np.linalg.norm(A)) * float(np.linalg.norm(dA))
        for A, dA in zip(As, dAs, strict=True)
    )
    scale = max(1.0, (float(np.linalg.norm(L)) + float(np.linalg.norm(G)) + 2.0 * rate_scale) * u2)
    gauge_cancel_residual = float(np.linalg.norm(moving + gauge_skew))
    direct_split_residual = float(np.linalg.norm(direct_native - native))
    outer_residual = float(np.linalg.norm(direct_native - (outer_interface - diagonal_resolved)))
    relink_pair_residual = float(np.linalg.norm(TK + TK.T))
    strain_pair_residual = float(np.linalg.norm(DS - DS.T))
    relink_row_residual = float(np.linalg.norm(TK.sum(axis=1) - relink))
    strain_row_residual = float(np.linalg.norm(DS.sum(axis=1) - strain))
    moving_total = float(moving.sum())
    relink_total = float(relink.sum())
    strain_total_residual = float(strain.sum() + 2.0 * np.real(np.vdot(u, S @ u)))
    native_total_residual = float(native.sum() + 2.0 * np.real(np.vdot(u, S @ u)))
    residuals = (
        gauge_cancel_residual,
        direct_split_residual,
        outer_residual,
        relink_pair_residual,
        strain_pair_residual,
        relink_row_residual,
        strain_row_residual,
        abs(moving_total),
        abs(relink_total),
        abs(strain_total_residual),
        abs(native_total_residual),
    )
    if max(residuals) > 3e-10 * scale:
        raise AssertionError("gauge-quotiented smooth quadratic interface lost its exact energy structure")

    work_certificate = GaugeQuotientedInterfaceWork(
        tuple(float(x) for x in native),
        tuple(float(x) for x in relink),
        tuple(float(x) for x in strain),
        transport_residual,
        skew_decomposition_residual,
        tuple(tuple(float(x) for x in row) for row in TK),
    )
    return {
        "moving_analysis_partition_work": moving,
        "signed_common_gauge_skew_work": gauge_skew,
        "signed_physical_relink_work": relink,
        "signed_existing_strain_work": strain,
        "signed_native_interface_work": native,
        "outer_equation_interface_work": outer_interface,
        "diagonal_resolved_role_work": diagonal_resolved,
        "physical_relink_pair_matrix": TK,
        "strain_synthesis_pair_matrix": DS,
        "identity_scale": scale,
        "gauge_transport_operator_residual": transport_residual,
        "gauge_work_cancellation_residual": gauge_cancel_residual,
        "native_gauge_quotient_residual": direct_split_residual,
        "native_outer_recombination_residual": outer_residual,
        "relink_pair_antisymmetry_residual": relink_pair_residual,
        "strain_pair_symmetry_residual": strain_pair_residual,
        "relink_pair_row_sum_residual": relink_row_residual,
        "strain_pair_row_sum_residual": strain_row_residual,
        "total_moving_analysis_work": moving_total,
        "total_physical_relink_work": relink_total,
        "global_strain_reconstruction_residual": strain_total_residual,
        "global_native_interface_reconstruction_residual": native_total_residual,
        "observer_partition_motion_charged_as_physics": False,
        "gauge_transport_verified": True,
        "work_certificate": work_certificate,
    }


def smooth_carrier_energy_identity(
    *,
    tensor: np.ndarray,
    state: np.ndarray,
    resolved_state: np.ndarray,
    smooth_envelope: np.ndarray,
    envelope_rate: np.ndarray,
    viscosity_operator: np.ndarray | None = None,
    viscosity: float = 0.0,
) -> dict[str, float]:
    """Direct Q^2-weighted Navier--Stokes carrier-energy identity.

    For u_t=-B(u,u)+nu D u, w=Q u and eta=Q^2,

      d||w||^2/dt + diss_Q
       = <u,dot eta u> - 2 Re<eta u,B(u,u)>.

    Repartitioning with V and h=u-V gives exact low-low, high-high and native
    moving-interface work.  The same native interface is obtained from the
    outer-role equation only after subtracting the diagonal L_V role work from
    the work of (dot Q+[L_V,Q])u.
    """
    u = np.asarray(state, dtype=complex)
    V = np.asarray(resolved_state, dtype=complex)
    Q = np.asarray(smooth_envelope, dtype=complex)
    dQ = np.asarray(envelope_rate, dtype=complex)
    n = len(u)
    if V.shape != (n,) or Q.shape != (n, n) or dQ.shape != (n, n):
        raise ValueError("matching smooth-carrier PDE data required")
    sa, lo, hi = _selfadjoint_positive_contraction_residual(Q)
    if max(sa, lo, hi, float(np.linalg.norm(dQ - dQ.conj().T))) > 3e-10 * max(1.0, float(n)):
        raise ValueError("positive self-adjoint smooth envelope and self-adjoint rate required")
    if viscosity_operator is None:
        D = np.zeros((n, n), dtype=complex)
    else:
        D = np.asarray(viscosity_operator, dtype=complex)
        if D.shape != (n, n):
            raise ValueError("matching viscosity operator required")
    nu = float(viscosity)
    if nu < 0 or not math.isfinite(nu):
        raise ValueError("finite nonnegative viscosity required")
    dscale = max(1.0, float(np.linalg.norm(D)))
    if float(np.linalg.norm(D - D.conj().T)) > 3e-10 * dscale:
        raise ValueError("self-adjoint viscosity generator required")
    if float(np.max(np.linalg.eigvalsh(0.5 * (D + D.conj().T)))) > 3e-10 * dscale:
        raise ValueError("nonpositive viscosity generator required")
    visc_commutator = float(np.linalg.norm(Q @ D - D @ Q))
    if visc_commutator > 3e-10 * max(1.0, float(np.linalg.norm(Q)) * dscale):
        raise ValueError("smooth Fourier envelope must commute with viscosity")

    eta = Q @ Q
    deta = dQ @ Q + Q @ dQ
    h = u - V
    Buu = bilinear_apply(tensor, u, u)
    BVV = bilinear_apply(tensor, V, V)
    Bhh = bilinear_apply(tensor, h, h)
    LVu = linearized_resolved(tensor, V, u)
    w = Q @ u
    du = -Buu + nu * (D @ u)
    dw = dQ @ u + Q @ du

    energy_rate = 2.0 * float(np.real(np.vdot(w, dw)))
    viscous_dissipation = -2.0 * nu * float(np.real(np.vdot(w, D @ w)))
    if viscous_dissipation < -3e-11 * max(1.0, abs(energy_rate)):
        raise AssertionError("nonpositive viscosity produced negative carrier dissipation")
    moving = float(np.real(np.vdot(u, deta @ u)))
    physical_nonlinear = -2.0 * float(np.real(np.vdot(eta @ u, Buu)))
    direct_rhs = moving + physical_nonlinear

    low_low = 2.0 * float(np.real(np.vdot(eta @ u, BVV)))
    high_high = -2.0 * float(np.real(np.vdot(eta @ u, Bhh)))
    native_interface = moving - 2.0 * float(np.real(np.vdot(eta @ u, LVu)))
    repartitioned_rhs = low_low + high_high + native_interface

    Lw = linearized_resolved(tensor, V, w)
    RQ = dQ @ u + Lw - Q @ LVu
    outer_interface_work = 2.0 * float(np.real(np.vdot(w, RQ)))
    diagonal_resolved_work = 2.0 * float(np.real(np.vdot(w, Lw)))
    outer_recombined_interface = outer_interface_work - diagonal_resolved_work

    lhs = energy_rate + viscous_dissipation
    scale = max(1.0, abs(lhs), abs(direct_rhs), abs(repartitioned_rhs))
    direct_residual = lhs - direct_rhs
    repartition_residual = direct_rhs - repartitioned_rhs
    outer_residual = native_interface - outer_recombined_interface
    if max(abs(direct_residual), abs(repartition_residual), abs(outer_residual)) > 2e-10 * scale:
        raise AssertionError("direct smooth-carrier energy identity failed")
    return {
        "carrier_energy_rate": energy_rate,
        "carrier_viscous_dissipation": viscous_dissipation,
        "moving_quadratic_weight_work": moving,
        "physical_nonlinear_work": physical_nonlinear,
        "low_low_work": low_low,
        "high_high_work": high_high,
        "native_interface_work": native_interface,
        "outer_equation_interface_work": outer_interface_work,
        "diagonal_resolved_role_work": diagonal_resolved_work,
        "direct_energy_identity_residual": direct_residual,
        "resolved_repartition_residual": repartition_residual,
        "outer_to_native_interface_residual": outer_residual,
        "viscosity_commutator_residual": visc_commutator,
    }


def positive_smooth_interface_split(
    interface_work: GaugeQuotientedInterfaceWork,
) -> dict[str, object]:
    """Hahn cover after observer-gauge motion has already been quotiented.

    This function deliberately does not accept raw arrays.  Causal relink/strain
    owners may be formed only from a :class:`GaugeQuotientedInterfaceWork`
    produced after exact common-transport verification.
    """
    total = np.asarray(interface_work.signed_native_interface_atoms, dtype=float)
    relink = np.asarray(interface_work.signed_physical_relink_atoms, dtype=float)
    strain = np.asarray(interface_work.signed_existing_strain_atoms, dtype=float)
    if total.ndim != 1 or len(total) == 0 or relink.shape != total.shape or strain.shape != total.shape:
        raise ValueError("matching nonempty gauge-quotiented smooth-interface work atoms required")
    if np.any(~np.isfinite(total)) or np.any(~np.isfinite(relink)) or np.any(~np.isfinite(strain)):
        raise ValueError("finite gauge-quotiented smooth-interface work atoms required")
    scale = max(float(np.max(np.abs(total))), float(np.max(np.abs(relink))), float(np.max(np.abs(strain))))
    identity_residual = float(np.max(np.abs(total - relink - strain)))
    if identity_residual > 5e-12 * scale:
        raise ValueError("native interface atoms must split exactly into physical relink plus strain")
    W = float(np.maximum(total, 0.0).sum())
    WR = float(np.maximum(relink, 0.0).sum())
    WS = float(np.maximum(strain, 0.0).sum())
    cover = WR + WS - W
    work_scale = max(W, WR, WS)
    tol = 8e-13 * work_scale
    if cover < -tol:
        raise AssertionError("positive smooth-interface work escaped physical relink+strain Hahn cover")
    owners: list[str] = []
    threshold = 0.5 * W
    if W > 0:
        if WR + tol >= threshold:
            owners.append(RELINK_OWNER)
        if WS + tol >= threshold:
            owners.append(STRAIN_OWNER)
        if not owners:
            raise AssertionError("positive smooth interface lost both native physical owners")
    return {
        "positive_native_interface_work": W,
        "positive_conservative_relink_work": WR,
        "positive_existing_strain_work": WS,
        "positive_cover_margin": cover,
        "owner_threshold": threshold,
        "joint_physical_owners": tuple(owners),
        "signed_identity_max_residual": identity_residual,
        "gauge_transport_operator_residual": interface_work.gauge_transport_operator_residual,
        "new_interface_currency_created": False,
        "observer_partition_motion_charged_as_physics": False,
        "primary_selected": False,
    }


def coefficient_obstruction_energy_reentry(
    *,
    terminal_coefficient: complex,
    terminal_probe_l2: float,
    terminal_carrier_energy: float,
    initial_carrier_energy: float,
    strain_action: float,
    coefficient_obstruction_impulse: complex,
    interface_work: GaugeQuotientedInterfaceWork,
) -> dict[str, object]:
    """Resolve a coefficient obstruction only through actual Q^2 energy/work.

    The coefficient impulse is an interval locator and is absent from every
    physical threshold.  The interface input must already carry an exact
    observer-gauge quotient certificate.  Only then may positive interface work
    become conservative relink provenance or existing-strain ownership.
    """
    probe = float(terminal_probe_l2)
    E1 = float(terminal_carrier_energy)
    E0 = float(initial_carrier_energy)
    K = float(strain_action)
    impulse = complex(coefficient_obstruction_impulse)
    z = complex(terminal_coefficient)
    vals = (probe, E1, E0, K, z.real, z.imag, impulse.real, impulse.imag)
    if not all(math.isfinite(x) for x in vals) or probe <= 0 or E1 <= 0 or min(E0, K) < 0:
        raise ValueError("finite coefficient, positive carrier energy/probe and nonnegative history required")
    coefficient_lower = abs(z) ** 2 / (probe * probe)
    tol = 8e-13 * max(1.0, E1, coefficient_lower)
    if E1 + tol < coefficient_lower:
        raise ValueError("terminal carrier energy violates its registered coefficient lower bound")

    split = positive_smooth_interface_split(interface_work)
    interface_positive = float(split["positive_native_interface_work"])
    gate = route_physical_energy_causality(
        terminal_energy=E1,
        initial_energy=E0,
        residual_positive_work=interface_positive,
        strain_action=K,
    )
    gate_branch = str(gate["branch"])
    if gate_branch == "classified_residual_physical_work":
        branch = "smooth_interface_physical_work"
        owners = tuple(split["joint_physical_owners"])
        if not owners:
            raise AssertionError("energy-selected smooth-interface work has no physical owner")
    else:
        branch = gate_branch
        owners = ()
    return {
        "branch": branch,
        "energy_gate_branch": gate_branch,
        "energy_gate": gate,
        "terminal_coefficient_energy_lower": coefficient_lower,
        "terminal_carrier_energy": E1,
        "positive_native_interface_work": interface_positive,
        "joint_interface_owners": owners,
        "coefficient_obstruction_magnitude": abs(impulse),
        "coefficient_obstruction_is_interval_locator": True,
        "coefficient_impulse_used_as_physical_work": False,
        "observer_partition_motion_charged_as_physics": False,
        "causal_weight_source": "actual_gauge_quotiented_smooth_carrier_energy_work",
        "interface_split": split,
        "gauge_quotiented_interface_work_certificate": interface_work,
    }


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "native_object": "for a smooth PDE carrier w=Q u the physical energy is ||w||_2^2=<u,Q^2u>; Q is an analysis operator and eta=Q^2 is its energy effect",
        "hard_smooth_separation": "hard orthogonal P is read only at an actual event; a smooth non-idempotent Q with QP=P is propagated between events and is completed by a square partition",
        "quadratic_partition": "positive self-adjoint smooth carriers satisfy sum_a A_a^2=I; this identity alone is analysis bookkeeping and does not make arbitrary dot A physical",
        "direct_energy_law": "d||Q u||_2^2/dt+2nu||grad Q u||_2^2=<u,dot(Q^2)u>-2 Re<Q^2u,B(u,u)>",
        "resolved_repartition": "with V fixed and h=u-V, the direct law splits exactly into low-low work, q^2-weighted HH work and native interface <u,dot eta u>-2 Re<eta u,L_Vu>; selected support kills low-low",
        "outer_recombination": "native smooth-interface work equals the work of (dot Q+[L_V,Q])u minus the diagonal L_V work of Q u; the commutator must not be interpreted alone",
        "observer_gauge_quotient": "before assigning physical relink ownership require one common skew generator G with dot A_a+[G,A_a]=0 for every smooth role; this is the transported affine/Kelvin observer motion and its energy work cancels exactly",
        "physical_skew": "write the actual skew resolved operator as K=G+K_phys only after the common transport identity is verified; K_phys has antisymmetric pair flux, zero total work and is the only smooth skew contribution eligible for physical relink ownership",
        "symmetric": "the S work rows are symmetric synthesis-pair strain work and sum to the same full resolved strain/deformation; no new source or Xi is created",
        "forbidden_observer_motion": "a time-dependent square partition whose motion is not generated by the certified common PDE transport cannot be Hahn-routed to a physical owner, even though its channel-energy changes sum to zero",
        "forbidden_linear_complement": "I-Q is not the energy complement of a non-idempotent Q; its skew defect is 4 Re<Q(I-Q)u,K u>; use Q^2+R^2=I instead",
        "coefficient_reentry": "a large Duhamel/interface coefficient impulse only locates a stopping interval; actual Q^2 carrier energy and gauge-quotiented native interface work reenter the physical-energy gate, and the impulse magnitude is never used as work",
        "causal_route": "the energy gate returns inheritance, high strain, actual q^2-weighted HH generation, or positive gauge-quotiented native interface work; only the last is Hahn-routed to physical relink provenance or existing strain, with exact ties joint",
        "relation_to_donor_quotient": "the smooth Q^2 theorem and hard event-role donor theorem share resolved K/S provenance but remain separately typed measures; no same-measure claim is made without a physical-work pushforward",
        "scope": "this closes the local smooth-envelope/projector algebra, observer-gauge quotient and coefficient-obstruction energy-reentry boundary; it does not prove global recurrence termination, UV closure or Navier-Stokes regularity",
    }


@dataclass(frozen=True)
class SmoothQuadraticCarrierStress:
    samples: int
    worst_quadratic_partition_residual: float
    worst_partition_derivative_residual: float
    worst_gauge_transport_residual: float
    worst_gauge_work_cancellation_residual: float
    worst_native_outer_recombination_residual: float
    worst_relink_conservation_residual: float
    worst_strain_reconstruction_residual: float
    worst_pair_antisymmetry_residual: float
    worst_pair_row_sum_residual: float
    worst_direct_energy_identity_residual: float
    worst_resolved_repartition_residual: float
    worst_event_registration_residual: float
    linear_complement_counterexample_defect: float
    worst_quadratic_complement_skew_residual: float
    observer_motion_rejections: int
    minimum_energy_generation_margin: float
    branch_counts: dict[str, int]


def _random_transported_square_pair(
    rng: np.random.Generator,
    n: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    Z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    U, _ = np.linalg.qr(Z)
    theta = rng.uniform(0.0, 0.5 * math.pi, size=n)
    theta[0] = 0.0
    q = np.cos(theta)
    r = np.sin(theta)
    Q = (U * q) @ U.conj().T
    R = (U * r) @ U.conj().T
    H = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    G = 0.5 * (H - H.conj().T)
    dQ = -(G @ Q - Q @ G)
    dR = -(G @ R - R @ G)
    P = np.outer(U[:, 0], U[:, 0].conj())
    return Q, R, dQ, dR, P, G


def stress(samples: int = 50_000, seed: int = 20260810) -> SmoothQuadraticCarrierStress:
    rng = np.random.default_rng(seed)
    wp = wd = wg = wgc = wo = wr = ws = wanti = wrow = we = wrep = wreg = wq = 0.0
    min_generation = float("inf")
    counts: dict[str, int] = {}

    K0 = np.array([[0.0, -1.0], [1.0, 0.0]], dtype=complex)
    Q0 = np.diag([0.5, 0.0]).astype(complex)
    u0 = np.array([1.0, 1.0], dtype=complex)
    counter = hard_linear_complement_skew_defect(Q0, K0, u0)
    counter_defect = abs(float(counter["linear_complement_skew_defect"]))
    if counter_defect < 0.9 or abs(float(counter["overlap_formula_residual"])) > 2e-13:
        raise AssertionError("smooth-envelope linear-complement counterexample changed")

    # Permanent observer-motion guard: a user-chosen change of channel amplitudes
    # may conserve total partition energy but it is not generated by a common
    # skew transport and therefore cannot become a physical relink owner.
    obs_theta = np.array([0.35, 0.9])
    obs_rate = np.array([0.4, -0.2])
    Qobs = np.diag(np.cos(obs_theta)).astype(complex)
    Robs = np.diag(np.sin(obs_theta)).astype(complex)
    dQobs = np.diag(-np.sin(obs_theta) * obs_rate).astype(complex)
    dRobs = np.diag(np.cos(obs_theta) * obs_rate).astype(complex)
    observer_rejections = 0
    try:
        smooth_quadratic_interface_balance(
            (Qobs, Robs),
            (dQobs, dRobs),
            np.zeros((2, 2), dtype=complex),
            np.array([1.0, -0.7], dtype=complex),
            common_gauge_skew_operator=np.zeros((2, 2), dtype=complex),
        )
    except ValueError:
        observer_rejections = 1
    if observer_rejections != 1:
        raise AssertionError("arbitrary observer partition motion was admitted as physical relink")

    for _ in range(samples):
        n = int(rng.integers(2, 7))
        Q, R, dQ, dR, P, G = _random_transported_square_pair(rng, n)
        diag = quadratic_partition_diagnostics((Q, R), (dQ, dR))
        wp = max(wp, float(diag["quadratic_partition_residual"]))
        wd = max(wd, float(diag["quadratic_partition_derivative_residual"]))

        Hk = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        Kphys = 0.5 * (Hk - Hk.conj().T)
        Hs = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
        S = 0.5 * (Hs + Hs.conj().T)
        L = G + Kphys + S
        u = rng.normal(size=n) + 1j * rng.normal(size=n)
        bal = smooth_quadratic_interface_balance(
            (Q, R),
            (dQ, dR),
            L,
            u,
            common_gauge_skew_operator=G,
        )
        scale = float(bal["identity_scale"])
        wg = max(wg, float(bal["gauge_transport_operator_residual"]) / max(1.0, scale))
        wgc = max(wgc, abs(float(bal["gauge_work_cancellation_residual"])) / scale)
        wo = max(wo, abs(float(bal["native_outer_recombination_residual"])) / scale)
        wr = max(wr, abs(float(bal["total_physical_relink_work"])) / scale)
        ws = max(ws, abs(float(bal["global_strain_reconstruction_residual"])) / scale)
        wanti = max(
            wanti,
            float(bal["relink_pair_antisymmetry_residual"]) / scale,
            float(bal["strain_pair_symmetry_residual"]) / scale,
        )
        wrow = max(
            wrow,
            float(bal["relink_pair_row_sum_residual"]) / scale,
            float(bal["strain_pair_row_sum_residual"]) / scale,
        )

        T = rng.normal(size=(n, n, n)) + 1j * rng.normal(size=(n, n, n))
        V = rng.normal(size=n) + 1j * rng.normal(size=n)
        _eig, Uq = np.linalg.eigh(Q)
        D = (Uq * (-rng.uniform(0.0, 4.0, size=n))) @ Uq.conj().T
        nu = float(rng.uniform(0.0, 2.0))
        energy = smooth_carrier_energy_identity(
            tensor=T,
            state=u,
            resolved_state=V,
            smooth_envelope=Q,
            envelope_rate=dQ,
            viscosity_operator=D,
            viscosity=nu,
        )
        escale = max(1.0, abs(float(energy["carrier_energy_rate"])))
        we = max(we, abs(float(energy["direct_energy_identity_residual"])) / escale)
        wrep = max(wrep, abs(float(energy["resolved_repartition_residual"])) / escale)

        probe = rng.normal(size=n) + 1j * rng.normal(size=n)
        reg = abs(envelope_registration_residual(P, Q, u, probe))
        wreg = max(wreg, reg / max(1.0, float(np.linalg.norm(u)) * float(np.linalg.norm(probe))))
        if wreg > 5e-11:
            raise AssertionError("hard event role failed smooth-envelope plateau registration")

        Ktotal, _ = adjoint_split(L)
        comp = hard_linear_complement_skew_defect(Q, Ktotal, u)
        wq = max(wq, abs(float(comp["quadratic_complement_skew_residual"])) / scale)
        if wq > 5e-11:
            raise AssertionError("quadratic smooth complement lost skew conservation")

        work = bal["work_certificate"]
        Wpos = float(np.maximum(np.asarray(work.signed_native_interface_atoms), 0.0).sum())
        mode = int(rng.integers(0, 4))
        if mode == 2 and Wpos > 1e-12:
            E1 = max(1e-9, 3.0 * Wpos)
        else:
            E1 = max(float(rng.lognormal(mean=0.0, sigma=1.0)), 6.0 * Wpos + 1e-9)
        phase = float(rng.uniform(-math.pi, math.pi))
        z = math.sqrt(float(rng.uniform(0.05, 0.60)) * E1) * complex(math.cos(phase), math.sin(phase))
        Kaction = float(rng.uniform(0.0, 1.0 / 30.0))
        E0 = 0.05 * E1
        if mode == 0:
            Kaction = float(rng.uniform(1.0 / 30.0 + 1e-5, 0.05))
        elif mode == 1:
            E0 = 0.25 * E1
        out = coefficient_obstruction_energy_reentry(
            terminal_coefficient=z,
            terminal_probe_l2=1.0,
            terminal_carrier_energy=E1,
            initial_carrier_energy=E0,
            strain_action=Kaction,
            coefficient_obstruction_impulse=complex(rng.normal(), rng.normal()) * 1000.0,
            interface_work=work,
        )
        branch = str(out["branch"])
        counts[branch] = counts.get(branch, 0) + 1
        if bool(out["coefficient_impulse_used_as_physical_work"]):
            raise AssertionError("coefficient obstruction was promoted to physical work")
        if bool(out["observer_partition_motion_charged_as_physics"]):
            raise AssertionError("observer partition motion was promoted to physical work")
        if branch == "physical_high_high_transfer_generation":
            gate = dict(out["energy_gate"])
            margin = float(gate["physical_hh_work_lower"]) - float(gate["clean_threshold"])
            min_generation = min(min_generation, margin)

    if not math.isfinite(min_generation):
        min_generation = 0.0
    return SmoothQuadraticCarrierStress(
        samples,
        wp,
        wd,
        wg,
        wgc,
        wo,
        wr,
        ws,
        wanti,
        wrow,
        we,
        wrep,
        wreg,
        counter_defect,
        wq,
        observer_rejections,
        min_generation,
        counts,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-smooth-quadratic-carrier-interface"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    payload = {"certificate": cert, "stress": asdict(out)}
    (args.outdir / "smooth_quadratic_carrier_interface.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    md = f"""# Smooth quadratic-carrier energy with observer-gauge quotient

Status: **{cert['status']}**.

For a smooth self-adjoint PDE carrier `w=Q u`, the native carrier energy is

`E_Q=||Q u||_2^2=<u,Q^2u>`.

Smooth roles are therefore completed by `sum_a A_a^2=I`; `Q` and `I-Q` are not treated as hard projectors.  Direct Navier--Stokes differentiation gives

`d||Q u||_2^2/dt + 2 nu ||grad Q u||_2^2 = <u,dot(Q^2)u> - 2 Re<Q^2u,B(u,u)>`.

After resolved repartition the native interface is

`J_Q=<u,dot(Q^2)u>-2 Re<Q^2u,L_Vu>`,

and exactly equals outer commutator work minus the diagonal `L_V` work of `Q u`.  The commutator is never interpreted alone.

The crucial quotient is prior to ownership.  A moving square partition is accepted as transported PDE gauge only when one common skew generator `G` satisfies

`dot A_a + [G,A_a] = 0`

for every role.  Writing the actual skew resolved operator as `K=G+K_phys`, the `dot A` work cancels the common-`G` work exactly.  Only `K_phys` remains as antisymmetric conservative physical relink; `S` remains the existing symmetric strain/deformation work.  Arbitrary observer motion of the analysis windows fails this transport identity and cannot be Hahn-routed to a causal owner.

A coefficient obstruction remains only an interval locator.  Actual `Q^2` carrier energy and a **gauge-quotiented** interface-work certificate must pass through the physical-energy gate before inheritance, high strain, HH generation, physical relink provenance or strain ownership is named.  The coefficient impulse magnitude never enters a work threshold.

Stress: `{out.samples}` transported square-partition/interface/PDE/reentry states
- worst quadratic-partition residual: `{out.worst_quadratic_partition_residual:.3e}`
- worst differentiated-partition residual: `{out.worst_partition_derivative_residual:.3e}`
- worst common-gauge transport residual: `{out.worst_gauge_transport_residual:.3e}`
- worst gauge-work cancellation residual: `{out.worst_gauge_work_cancellation_residual:.3e}`
- worst native/outer recombination residual: `{out.worst_native_outer_recombination_residual:.3e}`
- worst physical-relink total residual: `{out.worst_relink_conservation_residual:.3e}`
- worst strain reconstruction residual: `{out.worst_strain_reconstruction_residual:.3e}`
- worst pair antisymmetry/symmetry residual: `{out.worst_pair_antisymmetry_residual:.3e}`
- worst pair row-sum residual: `{out.worst_pair_row_sum_residual:.3e}`
- worst direct carrier-energy residual: `{out.worst_direct_energy_identity_residual:.3e}`
- worst resolved-repartition residual: `{out.worst_resolved_repartition_residual:.3e}`
- worst hard-event plateau registration residual: `{out.worst_event_registration_residual:.3e}`
- linear-complement counterexample defect: `{out.linear_complement_counterexample_defect:.6f}`
- worst quadratic-complement skew residual: `{out.worst_quadratic_complement_skew_residual:.3e}`
- arbitrary observer-motion rejections: `{out.observer_motion_rejections}`
- minimum clean HH-generation margin: `{out.minimum_energy_generation_margin:.3e}`
- energy-reentry branches: `{out.branch_counts}`

This theorem closes the local Q-vs-Q^2 mismatch without charging observer-selected window motion as physics.  It remains complementary to the hard event-role donor/circulation quotient.  It does not prove global owner termination, UV closure, or Navier--Stokes regularity.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
