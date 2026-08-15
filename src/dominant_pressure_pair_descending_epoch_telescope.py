from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

from src.objective_pressure_pair_atomization import (
    DEFAULT_PAIR_DOMINANCE,
    clean_dominant_pair_shell_mass_lower,
)
from src.objective_source_routing_compiler import objective_owner_weight_threshold


STATUS = (
    "DRAFT_DOMINANT_PRESSURE_PAIR_DESCENDING_EPOCH_TELESCOPE__"
    "CERTIFIED_QUARTER_PAIR_CHILD_SCALE_AT_MOST_ONE_QUARTER__"
    "ACTUAL_CRITICAL_SHELL_MASS_PLUS_GLOBAL_ENERGY_FLOOR__"
    "NO_GENERIC_SHELL_PROGRESS_NO_RESET_COUNT_NO_PRESSURE_ENTROPY_COST"
)

PRESSURE_PAIR_SCALE_RATIO_UPPER = 0.25
CANONICAL_PAIR_DOMINANCE = float(DEFAULT_PAIR_DOMINANCE)


@dataclass(frozen=True)
class DominantPressurePairRenewalStep:
    """One consecutive resolved pressure-pair owner and its selected hard child.

    The step records the *actual* child hard-shell frequency selected from the
    certified resolved pair event.  The pressure theorem supplies two facts used
    here and nowhere else:

      child_frequency <= parent_frequency/4,
      child_critical_mass >= 80 pressure_source_weight / c

    on the canonical q_max>=1/4 dominant face.  ``child_critical_mass`` denotes
    the physical quantity M ||P_M u||_2^2 of the selected u-shell, not a generic
    carrier registration and not a downstream same-corridor service witness.
    """

    parent_frequency: float
    child_frequency: float
    pressure_source_weight: float
    normalized_pair_mass: float
    child_critical_mass: float
    scaled_lifetime: float

    def __post_init__(self) -> None:
        vals = (
            self.parent_frequency,
            self.child_frequency,
            self.pressure_source_weight,
            self.child_critical_mass,
            self.scaled_lifetime,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in vals):
            raise ValueError("positive finite dominant pressure-pair step data required")
        q = float(self.normalized_pair_mass)
        if not math.isfinite(q) or not (0.0 < q <= 1.0):
            raise ValueError("normalized pair mass must lie in (0,1]")
        if q + 2.0e-13 < CANONICAL_PAIR_DOMINANCE:
            raise ValueError("step is not on the certified quarter-dominant pressure-pair face")
        freq_tol = 8.0e-13 * max(1.0, self.parent_frequency, self.child_frequency)
        if self.child_frequency > PRESSURE_PAIR_SCALE_RATIO_UPPER * self.parent_frequency + freq_tol:
            raise ValueError("resolved pressure-pair child must satisfy N_next<=N/4")
        clean = clean_dominant_pair_shell_mass_lower(
            self.pressure_source_weight,
            self.scaled_lifetime,
        )
        mass_tol = 8.0e-13 * max(1.0, clean, self.child_critical_mass)
        if self.child_critical_mass + mass_tol < clean:
            raise ValueError("selected child shell lost the certified 80 Sigma_P/c critical-mass lower")


@dataclass(frozen=True)
class DominantPressurePairEpochCertificate:
    transition_count: int
    event_count: int
    root_frequency: float
    final_child_frequency: float
    global_energy_upper: float
    pressure_source_weight_floor: float
    scaled_lifetime: float
    child_critical_mass_floor: float
    physical_frequency_floor: float
    maximum_transition_count: int
    maximum_event_count: int
    maximum_observed_scale_ratio: float
    generic_shell_registration_used_as_progress: bool = False
    pressure_entropy_used_as_cost: bool = False
    critical_mass_used_as_additive_reset: bool = False
    global_time_clock_used: bool = False

    def __post_init__(self) -> None:
        if self.transition_count < 1:
            raise ValueError("nonempty dominant pressure-pair transition epoch required")
        if self.event_count != self.transition_count:
            # One row is one pressure-pair owner event followed by its child.
            raise ValueError("event count must equal the number of pressure-pair renewal rows")
        positive = (
            self.root_frequency,
            self.final_child_frequency,
            self.global_energy_upper,
            self.pressure_source_weight_floor,
            self.scaled_lifetime,
            self.child_critical_mass_floor,
            self.physical_frequency_floor,
        )
        if not all(math.isfinite(x) and x > 0.0 for x in positive):
            raise ValueError("positive finite epoch certificate values required")
        if self.maximum_transition_count < self.transition_count:
            raise ValueError("frequency/energy bound does not cover observed dominant pressure epoch")
        if self.maximum_event_count != self.maximum_transition_count:
            raise ValueError("event-count convention disagrees with transition rows")
        if self.maximum_observed_scale_ratio > PRESSURE_PAIR_SCALE_RATIO_UPPER + 8.0e-13:
            raise ValueError("epoch contains a non-pressure-pair scale transition")
        if (
            self.generic_shell_registration_used_as_progress
            or self.pressure_entropy_used_as_cost
            or self.critical_mass_used_as_additive_reset
            or self.global_time_clock_used
        ):
            raise ValueError("dominant pressure epoch used a forbidden non-native shortcut")


def dominant_pressure_pair_source_floor_from_objective_stop(
    objective_variation_action_floor: float,
    scaled_lifetime: float,
) -> float:
    """Uniform pressure-owner floor when pressure is a qualifying objective owner."""
    return objective_owner_weight_threshold(
        objective_variation_action_floor,
        scaled_lifetime,
    )


def dominant_pressure_pair_child_mass_floor(
    pressure_source_weight_floor: float,
    scaled_lifetime: float,
) -> float:
    """Canonical quarter-dominant physical u-shell lower, exactly 80 Sigma/c."""
    return clean_dominant_pair_shell_mass_lower(
        pressure_source_weight_floor,
        scaled_lifetime,
    )


def dominant_pressure_pair_physical_frequency_floor(
    global_energy_upper: float,
    pressure_source_weight_floor: float,
    scaled_lifetime: float,
) -> float:
    """Use M||P_Mu||^2>=mu_* and ||P_Mu||^2<=E_global.

    This is a lower bound on the actual selected child frequency,

        M >= mu_*/E_global,

    not an infrared cutoff inserted into Navier--Stokes.
    """
    E = float(global_energy_upper)
    if not math.isfinite(E) or E <= 0.0:
        raise ValueError("positive finite global kinetic-energy upper required")
    mu = dominant_pressure_pair_child_mass_floor(
        pressure_source_weight_floor,
        scaled_lifetime,
    )
    return mu / E


def _transition_count_upper(root_frequency: float, frequency_floor: float) -> int:
    N0 = float(root_frequency)
    Nmin = float(frequency_floor)
    if not all(math.isfinite(x) and x > 0.0 for x in (N0, Nmin)):
        raise ValueError("positive finite root/floor frequencies required")
    ratio = N0 / Nmin
    # L transitions require N_L<=N0*4^-L and N_L>=Nmin, hence 4^L<=ratio.
    if ratio < 4.0 * (1.0 - 8.0e-13):
        return 0
    raw = math.log(ratio) / math.log(4.0)
    k = max(0, int(math.floor(raw + 8.0e-13)))
    # Repair only floating-point boundary ambiguity; the mathematical bound is
    # floor(log_4(N0/Nmin)).
    while k > 0 and k * math.log(4.0) > math.log(ratio) + 2.0e-12:
        k -= 1
    while (k + 1) * math.log(4.0) <= math.log(ratio) + 2.0e-12:
        k += 1
    return k


def dominant_pressure_pair_epoch_telescope(
    steps: Sequence[DominantPressurePairRenewalStep],
    *,
    global_energy_upper: float,
    pressure_source_weight_floor: float,
    scaled_lifetime: float,
) -> DominantPressurePairEpochCertificate:
    """Close one maximal consecutive quarter-dominant pressure-pair epoch.

    Every row is one genuine pressure first-stop whose certified resolved pair
    owner is quarter-dominant and whose actual hard-shell child is selected as the
    next recursive state.  The proof uses only

      N_(j+1) <= N_j/4,
      N_(j+1)||P_(N_(j+1))u||_2^2 >= mu_*,
      ||P_Mu||_2^2 <= ||u||_2^2 <= E_global.

    Therefore N_(j+1)>=mu_*/E_global for every row, while geometrically
    N_L<=N_0 4^-L.  Hence

      4^L <= N_0 E_global / mu_*.

    No generic A=3M/4 shell registration, pressure H2, additive critical-mass
    reset, or artificial common time clock enters the count.
    """
    rows = tuple(steps)
    if not rows:
        raise ValueError("nonempty dominant pressure-pair epoch required")
    E = float(global_energy_upper)
    sigma_floor = float(pressure_source_weight_floor)
    c = float(scaled_lifetime)
    if not all(math.isfinite(x) and x > 0.0 for x in (E, sigma_floor, c)):
        raise ValueError("positive finite energy/source-floor/lifetime required")

    root = rows[0].parent_frequency
    max_ratio = 0.0
    for i, row in enumerate(rows):
        ctol = 8.0e-13 * max(1.0, c, row.scaled_lifetime)
        if abs(row.scaled_lifetime - c) > ctol:
            raise ValueError("one epoch must use the same registered scaled lifetime c")
        if row.pressure_source_weight + 8.0e-13 * max(1.0, sigma_floor) < sigma_floor:
            raise ValueError("pressure source owner fell below the uniform epoch floor")
        if i > 0:
            prev = rows[i - 1].child_frequency
            ftol = 8.0e-13 * max(1.0, prev, row.parent_frequency)
            if abs(row.parent_frequency - prev) > ftol:
                raise ValueError("consecutive pressure-pair epoch must recurse through the selected child shell")
        ratio = row.child_frequency / row.parent_frequency
        max_ratio = max(max_ratio, ratio)
        # The actual shell energy is bounded by the full kinetic energy.  Since
        # row.child_critical_mass is an actual lower registered on that shell,
        # it must be compatible with M E_global.
        physical_upper = row.child_frequency * E
        mtol = 8.0e-13 * max(1.0, physical_upper, row.child_critical_mass)
        if row.child_critical_mass > physical_upper + mtol:
            raise ValueError("pressure-pair child shell contradicts the supplied global energy upper")

    mu_floor = dominant_pressure_pair_child_mass_floor(sigma_floor, c)
    Nmin = mu_floor / E
    count_upper = _transition_count_upper(root, Nmin)
    if len(rows) > count_upper:
        raise ValueError("observed dominant pressure-pair epoch exceeds the physical scale/energy telescope")

    final_child = rows[-1].child_frequency
    if final_child + 8.0e-13 * max(1.0, final_child, Nmin) < Nmin:
        raise AssertionError("final dominant pressure child fell below its physical energy-imposed frequency floor")

    return DominantPressurePairEpochCertificate(
        transition_count=len(rows),
        event_count=len(rows),
        root_frequency=root,
        final_child_frequency=final_child,
        global_energy_upper=E,
        pressure_source_weight_floor=sigma_floor,
        scaled_lifetime=c,
        child_critical_mass_floor=mu_floor,
        physical_frequency_floor=Nmin,
        maximum_transition_count=count_upper,
        maximum_event_count=count_upper,
        maximum_observed_scale_ratio=max_ratio,
    )


def theorem_certificate() -> dict[str, object]:
    sigma = dominant_pressure_pair_source_floor_from_objective_stop(1.0, 1.0)
    mu = dominant_pressure_pair_child_mass_floor(sigma, 1.0)
    if abs(sigma - 0.25) > 1.0e-14:
        raise AssertionError("objective four-owner pressure floor lost A/(4c)")
    if abs(mu - 20.0) > 1.0e-12:
        raise AssertionError("objective-stop dominant pressure child lower lost 20 A/c^2")
    return {
        "status": STATUS,
        "source_floor": "a qualifying pressure owner at objective action floor A_* carries Sigma_P>=A_*/(4c)",
        "dominant_shell": "on the resolved pressure-pair owner with q_max>=1/4, the certified actual u-shell satisfies mu_child>=80 Sigma_P/c, hence >=20 A_*/c^2 under the objective owner floor",
        "scale_progress": "the same resolved pressure-pair supplier has an actual child hard-shell frequency N_next<=N/4; generic carrier A=3M/4 is not used as progress",
        "energy_floor": "mu_child=N_next||P_Nnext u||_2^2<=N_next E_global gives N_next>=mu_*/E_global",
        "epoch_bound": "for L consecutive dominant pressure-pair transitions, 4^L<=N_0 E_global/mu_*, so L<=floor(log_4(N_0 E_global/mu_*))",
        "pressure_entropy": "diffuse pressure H2 is not used as a cost or causal probability; this theorem closes only the quarter-dominant resolved-pair subepoch",
        "forbidden": "no critical-mass additive reset, no generic-shell scale progress, no theorem-depth event, no synthetic common clock",
        "scope": "draft eventually-pure dominant pressure-pair epoch telescope only; diffuse pressure, fresh SGS, cross-family source/strain/HH recurrence, initial/singular interfaces, and global regularity remain open",
    }
