from __future__ import annotations

import math
from dataclasses import dataclass

STATUS = "DRAFT_NATIVE_SCALE_FREE_ACTION_SPEED__COMPACT_PRE_SINGULAR_NO_ZENO"

HIGH_STRAIN_THRESHOLD = 1.0 / 30.0


@dataclass(frozen=True)
class ActionSpeedLock:
    high_strain_rate_upper: float
    objective_rate_upper: float
    high_strain_time_floor: float
    objective_time_floor: float
    common_time_floor: float

    def __post_init__(self) -> None:
        vals = (
            self.high_strain_rate_upper,
            self.objective_rate_upper,
            self.high_strain_time_floor,
            self.objective_time_floor,
            self.common_time_floor,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in vals):
            raise ValueError("positive finite action-speed data required")
        if self.common_time_floor > min(self.high_strain_time_floor, self.objective_time_floor) + 1.0e-15:
            raise AssertionError("common action time floor exceeded one native face floor")


def objective_low_coeffs(
    global_energy: float,
    viscosity: float,
    scaled_lifetime: float,
    filter_kernel_l1: float,
    stress_reproducing_kernel_l32: float,
) -> tuple[float, float]:
    """Coefficients in dA_obj/dt <= a N^3 + b N^(5/2)."""
    E = float(global_energy)
    nu = float(viscosity)
    c = float(scaled_lifetime)
    g1 = float(filter_kernel_l1)
    crep = float(stress_reproducing_kernel_l32)
    if min(E, c, g1, crep) <= 0.0 or nu < 0.0:
        raise ValueError("positive energy/filter/lifetime and nonnegative viscosity required")
    if not all(math.isfinite(x) for x in (E, nu, c, g1, crep)):
        raise ValueError("finite objective-source data required")
    stress = crep * (g1 + 1.0)
    C3 = 1.0 / (1536.0 * math.pi**2) + 1.0 / 5700.0 + 2.0 * stress / 380.0
    a = c * C3 * E
    b = c * nu * math.sqrt(E) / 6000.0
    return a, b


def objective_global_rate_upper(
    global_energy: float,
    viscosity: float,
    scaled_lifetime: float,
    filter_kernel_l1: float,
    stress_reproducing_kernel_l32: float,
    smooth_objective_rate_coefficient: float,
) -> float:
    """Scale-free objective-action speed from simultaneous low/high scale bounds.

    Existing resolved-source calculus gives

        A'_obj,N <= a N^3 + b N^(5/2),

    while compact pre-singular smoothness of the same physical source gives

        A'_obj,N <= d N^-2.

    Splitting min(aN^3+bN^(5/2), dN^-2) into its two positive pieces yields

        sup_N A'_obj,N
        <= (2a)^(2/5)d^(3/5) + (2b)^(4/9)d^(5/9).

    No scale cutoff is inserted.
    """
    a, b = objective_low_coeffs(
        global_energy,
        viscosity,
        scaled_lifetime,
        filter_kernel_l1,
        stress_reproducing_kernel_l32,
    )
    c = float(scaled_lifetime)
    Ccirc = float(smooth_objective_rate_coefficient)
    if Ccirc <= 0.0 or not math.isfinite(Ccirc):
        raise ValueError("positive finite smooth objective-source coefficient required")
    d = c * Ccirc
    cubic_lock = (2.0 * a) ** (2.0 / 5.0) * d ** (3.0 / 5.0)
    viscous_lock = 0.0 if b == 0.0 else (2.0 * b) ** (4.0 / 9.0) * d ** (5.0 / 9.0)
    return cubic_lock + viscous_lock


def high_strain_global_rate_upper(filter_kernel_l1: float, grad_u_linf: float) -> float:
    """Scale-free dK_N/dt <= ||K_S||_1 ||grad u||_infty on a compact smooth interval."""
    g1 = float(filter_kernel_l1)
    G = float(grad_u_linf)
    if g1 <= 0.0 or G <= 0.0 or not all(math.isfinite(x) for x in (g1, G)):
        raise ValueError("positive finite filter L1 norm and gradient bound required")
    return g1 * G


def action_speed_lock(
    tau: float,
    *,
    global_energy: float,
    viscosity: float,
    scaled_lifetime: float,
    filter_kernel_l1: float,
    stress_reproducing_kernel_l32: float,
    smooth_objective_rate_coefficient: float,
    grad_u_linf: float,
) -> ActionSpeedLock:
    t = float(tau)
    if not (0.0 < t <= 0.1) or not math.isfinite(t):
        raise ValueError("0<tau<=1/10 required")
    Ck = high_strain_global_rate_upper(filter_kernel_l1, grad_u_linf)
    Co = objective_global_rate_upper(
        global_energy,
        viscosity,
        scaled_lifetime,
        filter_kernel_l1,
        stress_reproducing_kernel_l32,
        smooth_objective_rate_coefficient,
    )
    strain_floor = HIGH_STRAIN_THRESHOLD / Ck
    objective_floor = (t / 60.0) / Co
    return ActionSpeedLock(
        high_strain_rate_upper=Ck,
        objective_rate_upper=Co,
        high_strain_time_floor=strain_floor,
        objective_time_floor=objective_floor,
        common_time_floor=min(strain_floor, objective_floor),
    )


def action_time_floor(tau: float, **kwargs: float) -> float:
    return action_speed_lock(tau, **kwargs).common_time_floor


def maximum_action_vertices(physical_span: float, tau: float, **kwargs: float) -> int:
    T = float(physical_span)
    if T < 0.0 or not math.isfinite(T):
        raise ValueError("finite nonnegative physical span required")
    dt = action_time_floor(tau, **kwargs)
    return int(math.floor(T / dt + 2.0e-13)) + 1


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "high_strain": "dK_N/dt <= ||K_S||_1 ||grad u||_infty uniformly in N on every compact pre-singular interval; the native face is K=1/30",
        "objective_low_scale": "dA_obj,N/dt <= a N^3 + b N^(5/2) from the resolved physical-source calculus",
        "objective_high_scale": "the same objective source obeys dA_obj,N/dt <= d N^-2 by compact pre-singular smoothness",
        "scale_lock": "both bounds hold simultaneously, so sup_N dA_obj,N/dt <= (2a)^(2/5)d^(3/5)+(2b)^(4/9)d^(5/9)",
        "objective_face": "A_obj=tau/60 has a positive physical-time floor uniform in N",
        "no_zeno": "high-strain/objective-source ACTION first hits therefore occur only finitely many times on a compact pre-singular interval",
        "semantics": "this is a native physical-time speed lock, not an analyst scale partition, event gap assumption, or reset budget",
        "global_regularity_claimed": False,
    }
