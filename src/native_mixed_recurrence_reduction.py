from __future__ import annotations

import math
from dataclasses import dataclass

from src.native_closed_triad_no_free_power_branching import (
    POWER_CONTINUATION_UPPER,
    POWER_LOG_COST_LOWER,
    STATUS as CLOSED_TRIAD_CURRENT_STATUS,
)
from src.native_scale_free_action_speed import (
    STATUS as ACTION_SPEED_STATUS,
    maximum_action_vertices as _maximum_action_vertices,
)
from src.native_stock_no_mint import STATUS as MODE_STOCK_STATUS

STATUS = (
    "DRAFT_NATIVE_INTRINSIC_LAW_SET_AND_MIXED_RECURRENCE_REDUCTION__"
    "CLOSED_TRIAD_CURRENT_STOCK_CONTINUITY_ACTION_SPEED_LOCK__"
    "MASTER_POWER_COMPOSITION_CONDITIONAL"
)


@dataclass(frozen=True)
class ConditionalMasterPowerBound:
    depth: int
    action_vertex_upper_bound: int
    guaranteed_free_power_vertices: int
    finite_summable_prefactor: float
    remainder_upper: float
    cyclic_log_cost_lower: float
    exact_action_power_ties_joint: bool = True
    power_cocharge_on_action_vertices: bool = False
    master_composition_certified: bool = False

    def __post_init__(self) -> None:
        if self.depth < 0 or self.action_vertex_upper_bound < 0 or self.guaranteed_free_power_vertices < 0:
            raise ValueError("nonnegative recurrence counts required")
        expected = max(0, self.depth - self.action_vertex_upper_bound)
        if self.guaranteed_free_power_vertices != expected:
            raise AssertionError("guaranteed free POWER count must use the ACTION upper bound")
        vals = (self.finite_summable_prefactor, self.remainder_upper, self.cyclic_log_cost_lower)
        if not all(math.isfinite(x) for x in vals) or self.finite_summable_prefactor <= 0.0:
            raise ValueError("positive finite conditional master data required")
        if not self.exact_action_power_ties_joint or self.power_cocharge_on_action_vertices:
            raise ValueError("ACTION/POWER exact ties stay joint and are not double charged")
        if self.master_composition_certified:
            raise ValueError("master-depth composition of same-event cyclic attenuation is still a theorem obligation")


def maximum_action_vertices(physical_span: float, tau: float, **speed_data: float) -> int:
    return _maximum_action_vertices(physical_span, tau, **speed_data)


def conditional_master_power_bound(
    depth: int,
    action_vertex_upper_bound: int,
    *,
    finite_summable_prefactor: float = 1.0,
) -> ConditionalMasterPowerBound:
    """Candidate master-facing algebra once the continuation-composition seam is proved.

    This function does *not* certify that the master continuation variable composes
    the local closed-triad ratios.  It records the exact algebraic consequence if
    that missing physical wiring theorem is established:

        R_L <= C_fin (10/13)^(L-N_ACTION),
        -log(R_L/C_fin) >= log(13/10)(L-N_ACTION).

    N_ACTION is an *upper bound* on ACTION-containing vertices.  Therefore at
    least max(0,L-N_ACTION) vertices of a depth-L survivor are free POWER.
    Exact ACTION/POWER ties remain joint causes; a tied vertex is covered by the
    ACTION count and receives no additional POWER co-charge.
    """
    L = int(depth)
    A = int(action_vertex_upper_bound)
    C = float(finite_summable_prefactor)
    if L < 0 or A < 0:
        raise ValueError("require L>=0 and N_ACTION>=0")
    if C <= 0.0 or not math.isfinite(C):
        raise ValueError("positive finite pre-existing prefactor required")
    free = max(0, L - A)
    upper = C * float(POWER_CONTINUATION_UPPER) ** free
    log_lower = POWER_LOG_COST_LOWER * free
    return ConditionalMasterPowerBound(
        depth=L,
        action_vertex_upper_bound=A,
        guaranteed_free_power_vertices=free,
        finite_summable_prefactor=C,
        remainder_upper=upper,
        cyclic_log_cost_lower=log_lower,
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "native_dependencies": {
            "closed_triad_current": CLOSED_TRIAD_CURRENT_STATUS,
            "mode_stock_continuity": MODE_STOCK_STATUS,
            "local_action_speed_lock": ACTION_SPEED_STATUS,
        },
        "primitive_law_set": {
            "I_closed_triad_current": "nonlinear modal energy routing on one closed helical triad is one scalar current R; energy and signed helicity are two null laws of that same current",
            "II_mode_stock_continuity": "E_A(t1)+D_A+Phi_out=E_A(t0)+Phi_in; nonlinearity transfers physical stock, viscosity sinks it, and owner relabelling cannot mint stock",
            "III_local_action_speed_lock": "high-strain and objective-source native thresholds have a scale-free positive physical-time price on every compact pre-singular interval",
        },
        "ontology_status": "STOCK | SIGNED POWER | LOCAL ACTION is the recurrence ontology generated above these smaller laws; material/source/strain/HH/contact/relink labels are manifestations rather than additional primitive mechanisms",
        "stock_consequence": "same-carrier inheritance/restriction and representation changes do not mint generation depth; fresh stock needs inherited energy or actual signed physical work",
        "action_consequence": "N_ACTION(I)<infinity on every compact pre-singular interval I",
        "power_bad": "geometry-bad/nonforward positive canonical work is already the existing TRANSFER_WORK_LOSS route",
        "power_good_local": "every free geometry-good canonical POWER junction has same-event continuation <=10/13 and nonforward side diversion >=3/13",
        "restriction_semantics": "the cyclic law survives any common q>=0 restriction on the lifted closed-triad occurrence, including allowed time/hard-cell/resolved-contact/donor restrictions, provided canonical dW+ is inherited and no re-Hahn is performed",
        "temporal_guard": "dW- remains same-time donor provenance; no FIFO/LIFO/proportional temporal matching, Markov ancestry, or old-stock-to-later-withdrawal assignment is used",
        "tie_policy": "exact ACTION/POWER ties remain joint causes; for a cost proof they may be counted among the finite ACTION-containing vertices and the POWER co-charge is omitted to prevent double charging",
        "candidate_master_bound": "for L>=N_ACTION(I), R_L <= C_fin (10/13)^(L-N_ACTION(I)); for all L use exponent max(0,L-N_ACTION(I)); equivalently -log(R_L/C_fin) >= log(13/10) max(0,L-N_ACTION(I))",
        "master_composition_status": "conditional: still requires a theorem identifying the physical master continuation variable on which consecutive same-event closed-triad attenuations compose; donor dW- must not be treated as between-time stock",
        "if_wired": "if that composition theorem passes review, the mixed genuine-owner pillar has no separate pure-POWER frontier: after finitely many ACTION vertices, bad POWER exits and good POWER pays the intrinsic cyclic branching cost",
        "global_regularity_claimed": False,
    }
