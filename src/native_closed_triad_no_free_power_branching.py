from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

from src.cyclic_helical_triad_donor_kernel import (
    CHILD_TO_DONOR_HI,
    SIDE_TO_DONOR_LO,
    SignedGoodSideRecipientCertificate,
)

STATUS = (
    "DRAFT_NATIVE_CLOSED_TRIAD_NO_FREE_POWER_BRANCHING__"
    "SAME_EVENT_GOOD_CONTINUATION_AT_MOST_10_OVER_13__"
    "MASTER_COMPOSITION_NOT_YET_CLAIMED"
)

POWER_CONTINUATION_UPPER = Fraction(10, 13)
SIDE_DIVERSION_LOWER = Fraction(3, 13)
POWER_LOG_COST_LOWER = math.log(13.0 / 10.0)

if POWER_CONTINUATION_UPPER != CHILD_TO_DONOR_HI:
    raise AssertionError("cyclic good-child continuation constant changed")
if SIDE_DIVERSION_LOWER != SIDE_TO_DONOR_LO:
    raise AssertionError("cyclic side-diversion constant changed")


@dataclass(frozen=True)
class ClosedTriadCurrentLaw:
    """One scalar nonlinear current with energy and signed-helicity null laws.

    Here ``a_i=s_i |k_i|``.  The three cyclic modal energy works are

        T0=(a1-a2)R,  T1=(a2-a0)R,  T2=(a0-a1)R.

    This is an Eulerian same-triad identity.  It is not a temporal ancestry law.
    """

    signed_frequencies: tuple[float, float, float]
    common_current: float
    works: tuple[float, float, float]
    energy_residual: float
    signed_helicity_residual: float
    temporal_matching_used: bool = False

    def __post_init__(self) -> None:
        vals = (*self.signed_frequencies, self.common_current, *self.works)
        if not all(math.isfinite(float(x)) for x in vals):
            raise ValueError("finite closed-triad current data required")
        if self.temporal_matching_used:
            raise ValueError("closed-triad current is same-time Eulerian physics, not temporal matching")
        scale = max(1.0, *(abs(float(x)) for x in self.works))
        if abs(self.energy_residual) > 1.0e-12 * scale:
            raise AssertionError("closed-triad scalar current lost energy conservation")
        hscale = max(
            1.0,
            *(abs(a * t) for a, t in zip(self.signed_frequencies, self.works)),
        )
        if abs(self.signed_helicity_residual) > 1.0e-12 * hscale:
            raise AssertionError("closed-triad scalar current lost signed-helicity conservation")


def closed_triad_current_law(
    frequencies: Sequence[float], helicities: Sequence[int], common_current: float
) -> ClosedTriadCurrentLaw:
    if len(frequencies) != 3 or len(helicities) != 3:
        raise ValueError("exactly three closed-triad frequencies/helicities required")
    freq = tuple(float(x) for x in frequencies)
    hel = tuple(int(s) for s in helicities)
    if not all(math.isfinite(x) and x > 0.0 for x in freq):
        raise ValueError("positive finite modal frequencies required")
    if any(s not in (-1, 1) for s in hel):
        raise ValueError("helicities must lie in {-1,+1}")
    R = float(common_current)
    if not math.isfinite(R):
        raise ValueError("finite common current required")
    a = tuple(float(s) * k for s, k in zip(hel, freq))
    T = (
        (a[1] - a[2]) * R,
        (a[2] - a[0]) * R,
        (a[0] - a[1]) * R,
    )
    energy = math.fsum(T)
    helicity = math.fsum(ai * ti for ai, ti in zip(a, T))
    return ClosedTriadCurrentLaw(
        signed_frequencies=(a[0], a[1], a[2]),
        common_current=R,
        works=(T[0], T[1], T[2]),
        energy_residual=energy,
        signed_helicity_residual=helicity,
    )


@dataclass(frozen=True)
class WeightedGoodPowerRestriction:
    """Same-weight restriction of signed-good closed-triad occurrences.

    ``good_work`` is the selected positive forward child work, ``donor_work`` is
    its unique same-triad energy-donor loss, and ``side_work`` is the simultaneous
    positive nonforward sibling.  All three are restricted by the *same*
    nonnegative physical weight before integration.
    """

    occurrences: int
    good_work: float
    donor_work: float
    side_work: float
    continuation_ratio: float
    side_diversion_ratio: float
    local_log_cost: float
    side_is_existing_transfer_work_loss: bool = True
    canonical_good_cause_replaced: bool = False
    side_separately_charged: bool = False
    between_time_stock_claimed: bool = False
    temporal_matching_used: bool = False
    later_hahn_used: bool = False
    master_composition_certified: bool = False

    def __post_init__(self) -> None:
        if self.occurrences <= 0:
            raise ValueError("nonempty weighted good-power family required")
        vals = (
            self.good_work,
            self.donor_work,
            self.side_work,
            self.continuation_ratio,
            self.side_diversion_ratio,
            self.local_log_cost,
        )
        if not all(math.isfinite(float(x)) for x in vals):
            raise ValueError("finite weighted good-power data required")
        if min(self.good_work, self.donor_work, self.side_work) <= 0.0:
            raise ValueError("positive restricted physical work required")
        scale = max(self.donor_work, self.good_work + self.side_work, 1.0e-300)
        if abs(self.donor_work - self.good_work - self.side_work) > 6.0e-10 * scale:
            raise AssertionError("same-weight cyclic restriction lost D=G+B")
        if self.continuation_ratio > float(POWER_CONTINUATION_UPPER) + 2.0e-12:
            raise AssertionError("good continuation exceeded 10/13 of same-event donor work")
        if self.side_diversion_ratio < float(SIDE_DIVERSION_LOWER) - 2.0e-12:
            raise AssertionError("nonforward side diversion fell below 3/13 of same-event donor work")
        if self.local_log_cost + 2.0e-12 < POWER_LOG_COST_LOWER:
            raise AssertionError("native cyclic continuation cost fell below log(13/10)")
        if not self.side_is_existing_transfer_work_loss:
            raise AssertionError("cyclic side sibling must remain the existing nonforward transfer-loss cause")
        if (
            self.canonical_good_cause_replaced
            or self.side_separately_charged
            or self.between_time_stock_claimed
            or self.temporal_matching_used
            or self.later_hahn_used
        ):
            raise ValueError("cyclic power cost may not change cause, double charge, mint stock, match time, or re-Hahn")
        if self.master_composition_certified:
            raise ValueError("local cyclic attenuation is certified; arbitrary recursive composition is still an open seam")


def weighted_signed_good_restriction(
    certificates: Sequence[SignedGoodSideRecipientCertificate],
    weights: Sequence[float],
) -> WeightedGoodPowerRestriction:
    """Integrate the atomwise cyclic branching law under one common positive restriction.

    At each signed-good occurrence the existing cyclic theorem gives

        D = G + B,
        G/D < 10/13,
        B/D > 3/13,

    where B is positive nonforward canonical work.  Applying the same q>=0 to
    all three members and integrating preserves these inequalities.  This is a
    same-event Radon statement only; it does not identify donor work with a
    between-time stock and does not by itself prove multiplicative master-depth
    composition.
    """
    rows = tuple(certificates)
    qs = tuple(float(q) for q in weights)
    if len(rows) != len(qs) or not rows:
        raise ValueError("matching nonempty certificate/weight families required")
    if any((not math.isfinite(q)) or q < 0.0 for q in qs):
        raise ValueError("finite nonnegative restriction weights required")

    good_terms: list[float] = []
    donor_terms: list[float] = []
    side_terms: list[float] = []
    active = 0
    for cert, q in zip(rows, qs):
        if q == 0.0:
            continue
        D = float(cert.donor_negative_work)
        G = float(cert.recipient_work)
        B = float(cert.side_positive_work)
        scale = max(D, G + B, 1.0e-300)
        if abs(D - G - B) > 6.0e-10 * scale:
            raise AssertionError("input signed-good certificate lost same-triad donor split")
        if not cert.side_is_positive_nonforward or not cert.side_terminal_transfer_loss_is_existing_router_consequence:
            raise AssertionError("input side sibling is not existing positive nonforward transfer-loss work")
        if G / D >= float(POWER_CONTINUATION_UPPER):
            raise AssertionError("input signed-good child/donor ratio lost strict 10/13 ceiling")
        if B / D <= float(SIDE_DIVERSION_LOWER):
            raise AssertionError("input signed-good side/donor ratio lost strict 3/13 floor")
        good_terms.append(q * G)
        donor_terms.append(q * D)
        side_terms.append(q * B)
        active += 1

    Gq = math.fsum(good_terms)
    Dq = math.fsum(donor_terms)
    Bq = math.fsum(side_terms)
    if active == 0 or Dq <= 0.0:
        raise ValueError("restriction must retain positive donor-row work")
    ratio = Gq / Dq
    side_ratio = Bq / Dq
    return WeightedGoodPowerRestriction(
        occurrences=active,
        good_work=Gq,
        donor_work=Dq,
        side_work=Bq,
        continuation_ratio=ratio,
        side_diversion_ratio=side_ratio,
        local_log_cost=-math.log(ratio),
    )


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "primitive_below_power": "canonical nonlinear POWER is a manifestation of a smaller closed-triad current law, not a primitive HH-event currency",
        "closed_triad_current": "T0=(a1-a2)R, T1=(a2-a0)R, T2=(a0-a1)R with a_i=s_i|k_i|; therefore sum T_i=0 and sum a_i T_i=0",
        "two_null_laws_one_current": "energy conservation and signed-helicity conservation are simultaneous null directions of the same one-dimensional triad current",
        "signed_good_branching": "for every positive forward signed-good root, the unique donor loss splits at the same event as D=G+B, with G/D<10/13 and B/D>3/13",
        "same_weight_radon_restriction": "any common measurable q>=0 on the lifted closed-triad family preserves D_q=G_q+B_q, G_q<=10 D_q/13, B_q>=3 D_q/13",
        "downstream_restriction_scope": "the same-event law survives allowed time, hard-cell, resolved-contact and donor restrictions exactly when the same positive physical weight is inherited before coarsening and canonical dW+ is not re-Hahn split",
        "hard_cell_young_status": "hard-cell compression and Young/Christ are downstream bookkeeping and are not load-bearing for the native POWER branching cost",
        "side_fate": "B_q is real positive canonical work on the nonforward J=0 sibling and is already TRANSFER_WORK_LOSS; it is not dissipation or a reset budget",
        "native_local_power_cost": f"C_cyc=-log(G_q/D_q) >= log(13/10)={POWER_LOG_COST_LOWER:.15g}",
        "single_charge": "if C_cyc is used as the POWER continuation cost, the same side sibling is not separately charged; canonical dW+ remains the cause and no later Hahn split is introduced",
        "forbidden_interpretation": "dW- is same-time donor provenance, not between-time stock; no FIFO/LIFO/proportional temporal matching or Markov ancestry is asserted",
        "helicity_guard": "signed helicity enters only as the second null law of the same closed-triad current; absolute-helicity magnitude is not a transfer cost or finite reset budget",
        "open_master_seam": "the local ratio is exact, but composition across arbitrary recursive depth still requires a theorem identifying the master's continuation variable on which these same-event attenuations compose",
        "global_regularity_claimed": False,
    }
