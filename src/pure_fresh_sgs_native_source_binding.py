from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from src.descending_fresh_sgs_scale_epoch_telescope import FreshSGSScaleOwnerCertificate
from src.objective_source_routing_compiler import (
    objective_owner_weight_threshold,
    objective_sgs_integrated_square_service_lower,
)
from src.pure_fresh_sgs_pre_singular_exhaustion import (
    FixedSmoothSGSFilterEnvelope,
    PureFreshSGSFirstStopStep,
    PureFreshSGSPreSingularExhaustionCertificate,
    pure_fresh_sgs_pre_singular_exhaustion,
)

STATUS = (
    "DRAFT_NATIVE_OBJECTIVE_SGS_TO_PURE_FRESH_EXHAUSTION_BINDING__"
    "SOURCE_WEIGHT_AND_FORCED_SERVICE_ARE_ONE_CERTIFIED_LAW__FAIL_CLOSED"
)


def _positive(*xs: float) -> bool:
    return all(math.isfinite(float(x)) and float(x) > 0.0 for x in xs)


@dataclass(frozen=True)
class FreshSGSNativeSourceBinding:
    owner: FreshSGSScaleOwnerCertificate
    sgs_source_weight: float
    filter_l1: float
    lp_constant: float
    bernstein_constant: float

    def __post_init__(self) -> None:
        if not isinstance(self.owner, FreshSGSScaleOwnerCertificate):
            raise TypeError("typed fresh SGS scale owner required")
        if not _positive(
            self.sgs_source_weight,
            self.filter_l1,
            self.lp_constant,
            self.bernstein_constant,
        ):
            raise ValueError("positive finite objective-SGS binding data required")
        Y = objective_sgs_integrated_square_service_lower(
            self.sgs_source_weight,
            self.filter_l1,
            self.lp_constant,
            self.bernstein_constant,
        )
        tol = 8e-13 * max(1.0, Y, self.owner.forced_square_service_threshold)
        if abs(self.owner.forced_square_service_threshold - Y) > tol:
            raise ValueError(
                "fresh owner threshold is not the canonical image of its SGS source weight"
            )


def bind_fresh_owner_to_native_sgs_source(
    owner: FreshSGSScaleOwnerCertificate,
    sgs_source_weight: float,
    *,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> FreshSGSNativeSourceBinding:
    return FreshSGSNativeSourceBinding(
        owner,
        float(sgs_source_weight),
        float(filter_l1),
        float(lp_constant),
        float(bernstein_constant),
    )


def pure_fresh_source_floors_from_objective_action(
    objective_variation_action_floor: float,
    scaled_lifetime_upper: float,
    *,
    filter_l1: float,
    lp_constant: float,
    bernstein_constant: float,
) -> tuple[float, float]:
    A = float(objective_variation_action_floor)
    c = float(scaled_lifetime_upper)
    if not _positive(A, c, filter_l1, lp_constant, bernstein_constant):
        raise ValueError("positive finite objective-source floor data required")
    sigma = objective_owner_weight_threshold(A, c)
    Y = objective_sgs_integrated_square_service_lower(
        sigma, filter_l1, lp_constant, bernstein_constant
    )
    return sigma, Y


@dataclass(frozen=True)
class NativePureFreshSGSFirstStopStep:
    source_binding: FreshSGSNativeSourceBinding
    child_frequency: float
    child_critical_mass: float
    earlier_time: float
    later_time: float

    def __post_init__(self) -> None:
        if not isinstance(self.source_binding, FreshSGSNativeSourceBinding):
            raise TypeError("typed native objective-SGS source binding required")

    def core_step(self) -> PureFreshSGSFirstStopStep:
        return PureFreshSGSFirstStopStep(
            owner=self.source_binding.owner,
            child_frequency=float(self.child_frequency),
            child_critical_mass=float(self.child_critical_mass),
            sgs_source_weight=float(self.source_binding.sgs_source_weight),
            earlier_time=float(self.earlier_time),
            later_time=float(self.later_time),
        )


def native_pure_fresh_sgs_pre_singular_exhaustion(
    steps: Sequence[NativePureFreshSGSFirstStopStep],
    *,
    global_energy_upper: float,
    pre_singular_h1_seminorm_sq_upper: float,
    forced_square_service_floor: float,
    sgs_source_weight_floor: float,
    scaled_lifetime_upper: float,
    filter_envelope: FixedSmoothSGSFilterEnvelope,
) -> PureFreshSGSPreSingularExhaustionCertificate:
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty native fresh SGS word required")
    if any(not isinstance(x, NativePureFreshSGSFirstStopStep) for x in rows):
        raise TypeError("all fresh SGS steps must carry typed native source bindings")
    return pure_fresh_sgs_pre_singular_exhaustion(
        tuple(x.core_step() for x in rows),
        global_energy_upper=global_energy_upper,
        pre_singular_h1_seminorm_sq_upper=pre_singular_h1_seminorm_sq_upper,
        forced_square_service_floor=forced_square_service_floor,
        sgs_source_weight_floor=sgs_source_weight_floor,
        scaled_lifetime_upper=scaled_lifetime_upper,
        filter_envelope=filter_envelope,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "single_law": "the canonical compiler relation Y=C_Y Sigma_R is verified before a fresh source step can enter the exhaustion theorem",
        "objective_floor": "A_obj>=A_* and c<=c_* give Sigma_R>=A_*/(4c_*) on an SGS qualifying owner, then Y_*=C_Y Sigma_*",
        "fail_closed": "an unrelated positive scalar cannot be attached as source weight to a fresh service owner",
        "scope": "typed provenance adapter only; the PDE exhaustion inequalities remain in pure_fresh_sgs_pre_singular_exhaustion.py",
    }
