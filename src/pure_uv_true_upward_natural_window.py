from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.critical_shell_service_reentry import (
    critical_shell_bounded_service_lower,
    critical_shell_integrated_service_lower,
)
from src.cyclic_helical_triad_donor_kernel import (
    cyclic_triad_measure_kernel,
    register_closed_helical_triad,
)
from src.hard_tail_true_upward_supply import (
    HardTailUpwardSupplySplit,
    UpwardSupplyAtom,
    hard_tail_upward_supply_split,
    mode_radius,
    upward_owner_support_alternative,
)
from src.high_tail_natural_window_reentry import natural_window_geometry, temporal_concentration_statistics
from src.high_tail_ultraviolet_locality import ultraviolet_hh_work_constant
from src.resolved_contact_native_binding import canonical_hh_edge_total_variation_young_bridge

STATUS = (
    "EXACT_PURE_UV_TRUE_UPWARD_FIRST_SHELL_NATURAL_WINDOW__"
    "CANONICAL_RECIPIENT_SUBMEASURE_SINGLE_CHARGED__H_EQUALS_U_ON_BOTH_PARENTS__"
    "P_SCALE_ONE__EDGE_TOTAL_VARIATION_BEFORE_YOUNG__COMMON_N_DW_UNIT"
)

PURE_UV_SHELL_INDEX = 1
PURE_UV_SHELL_RATIO = 2.0
PURE_UV_PARENT_UPPER_RATIO = 3.0 / 2.0


def _finite_positive(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"positive finite {name} required")
    return out


def _finite_nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"finite nonnegative {name} required")
    return out


@dataclass(frozen=True)
class PureUVRecipientSubmeasure:
    """One true-upward pure-UV recipient cause after donor sidecars are coalesced.

    The selected same-time donors remain as provenance, but their masses are summed
    before downstream capacity/owner logic.  This is a positive submeasure of the
    already-canonical recipient dW+, not a new Hahn law.
    """

    boundary: float
    recipient_closed_mode_index: int
    donor_closed_mode_indices: tuple[int, ...]
    canonical_positive_submeasure_mass: float
    recipient_radius: float
    recipient_shell_scale: float
    interaction_parent_radii: tuple[float, float]
    common_unit_scale: float
    recipient_shell_index: int = PURE_UV_SHELL_INDEX
    p_scale: float = 1.0
    h_equals_u_on_both_parents: bool = True
    donor_sidecars_coalesced: bool = True
    later_hahn_used: bool = False
    owner_mass_cloned: bool = False
    recipient_shell_reweighting_used: bool = False
    output_scale_locality_theorem_used: bool = False
    creates_new_event_depth: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        mass = _finite_positive(self.canonical_positive_submeasure_mass, "canonical pure-UV recipient mass")
        rr = _finite_positive(self.recipient_radius, "recipient radius")
        M = _finite_positive(self.recipient_shell_scale, "recipient shell scale")
        _finite_positive(self.common_unit_scale, "common-unit scale")
        if not self.donor_closed_mode_indices or len(set(self.donor_closed_mode_indices)) != len(self.donor_closed_mode_indices):
            raise ValueError("nonempty unique pure-UV donor sidecars required")
        if self.recipient_closed_mode_index in self.donor_closed_mode_indices:
            raise ValueError("recipient cannot be its own pure-UV donor")
        if len(self.interaction_parent_radii) != 2 or min(self.interaction_parent_radii) <= 0.0:
            raise ValueError("two positive pure-UV interaction-parent radii required")
        if self.recipient_shell_index != PURE_UV_SHELL_INDEX:
            raise AssertionError("pure-UV true-upward cause escaped the first dyadic shell")
        if abs(M - PURE_UV_SHELL_RATIO * N) > 5.0e-12 * max(M, N):
            raise AssertionError("pure-UV recipient shell is not exactly M=2N")
        if not (0.5 * M < rr <= M):
            raise AssertionError("pure-UV recipient left its actual first shell")
        if min(self.interaction_parent_radii) <= 0.25 * M:
            raise AssertionError("pure-UV recipient cause touched resolved parent support")
        if max(self.interaction_parent_radii) > PURE_UV_PARENT_UPPER_RATIO * M + 5.0e-12 * M:
            raise AssertionError("pure-UV parent escaped the automatic 3M/2 corridor")
        if self.common_unit_scale != N:
            raise ValueError("pure-UV cause must remain in the parent-tail N dW unit")
        if abs(self.p_scale - 1.0) > 5.0e-15:
            raise ValueError("pure-UV first-shell support has p_scale=1; no output-scale probability is permitted")
        if not self.h_equals_u_on_both_parents:
            raise AssertionError("resolved cutoff must vanish on both strict pure-UV parents, so h=u there")
        if not self.donor_sidecars_coalesced:
            raise ValueError("pure-UV capacity may be charged only after donor sidecars are coalesced by recipient")
        if (
            self.later_hahn_used
            or self.owner_mass_cloned
            or self.recipient_shell_reweighting_used
            or self.output_scale_locality_theorem_used
            or self.creates_new_event_depth
        ):
            raise ValueError(
                "pure-UV binding may not re-Hahn, clone a recipient owner, reweight by M, invoke output locality, or add event depth"
            )
        if mass <= 0.0:
            raise AssertionError("pure-UV recipient submeasure lost positive canonical work")

    @property
    def common_unit_work_mass(self) -> float:
        return self.common_unit_scale * self.canonical_positive_submeasure_mass


@dataclass(frozen=True)
class PureUVFirstShellLaw:
    boundary: float
    recipient_shell_scale: float
    recipient_submeasures: tuple[PureUVRecipientSubmeasure, ...]
    total_canonical_positive_mass: float
    total_common_unit_work: float
    native_work_mass_scale: float
    coalescing_native_residual: float
    p_scale: float = 1.0
    h_inf_output_scale: float = 0.0
    later_hahn_used: bool = False
    owner_mass_cloned: bool = False
    locality_theorem_used: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        M = _finite_positive(self.recipient_shell_scale, "pure-UV shell scale")
        total = _finite_positive(self.total_canonical_positive_mass, "pure-UV positive mass")
        common = _finite_positive(self.total_common_unit_work, "pure-UV common work")
        native = _finite_positive(self.native_work_mass_scale, "native pure-UV work scale")
        _finite_nonnegative(self.coalescing_native_residual, "pure-UV coalescing residual")
        if not self.recipient_submeasures:
            raise ValueError("nonempty pure-UV recipient law required")
        if abs(M - 2.0 * N) > 5.0e-12 * max(M, N):
            raise AssertionError("pure-UV law lost exact first-shell M=2N geometry")
        if any(c.boundary != N or c.recipient_shell_scale != M for c in self.recipient_submeasures):
            raise ValueError("pure-UV recipient submeasures do not share one N-anchored first shell")
        if abs(math.fsum(c.canonical_positive_submeasure_mass for c in self.recipient_submeasures) - total) > 5.0e-10 * native:
            raise AssertionError("recipient coalescing changed canonical pure-UV positive mass")
        if abs(common - N * total) > 5.0e-10 * N * native:
            raise AssertionError("pure-UV law left the parent-tail N dW unit")
        if self.coalescing_native_residual > 5.0e-10:
            raise AssertionError("pure-UV donor sidecars failed to coalesce on the native work scale")
        if self.p_scale != 1.0 or self.h_inf_output_scale != 0.0:
            raise ValueError("pure-UV first-shell law has no output-scale entropy or selection loss")
        if self.later_hahn_used or self.owner_mass_cloned or self.locality_theorem_used:
            raise ValueError("pure-UV shell law may not re-Hahn, clone owners, or invoke locality machinery")


def coalesce_pure_uv_recipient_submeasures(
    atoms: Sequence[UpwardSupplyAtom],
    *,
    native_work_mass_scale: float,
) -> tuple[PureUVRecipientSubmeasure, ...]:
    """Forget selected donor labels only after grouping them by physical recipient."""

    rows = tuple(atoms)
    native = _finite_positive(native_work_mass_scale, "native work-mass scale")
    if not rows:
        raise ValueError("nonempty pure-UV upward atoms required")
    if any(not atom.pure_uv_hh_by_support for atom in rows):
        raise ValueError("recipient coalescing accepts only pure-UV true-upward atoms")
    groups: dict[int, list[UpwardSupplyAtom]] = {}
    for atom in rows:
        groups.setdefault(atom.recipient_closed_mode_index, []).append(atom)
    out: list[PureUVRecipientSubmeasure] = []
    for recipient, group in sorted(groups.items()):
        ref = group[0]
        donors = tuple(sorted(atom.donor_closed_mode_index for atom in group))
        if len(set(donors)) != len(donors):
            raise ValueError("duplicate donor sidecar would double-charge one pure-UV recipient")
        for atom in group:
            scale = max(native, atom.physical_work_mass, ref.physical_work_mass)
            if (
                atom.boundary != ref.boundary
                or atom.recipient_mode != ref.recipient_mode
                or atom.recipient_shell_index != ref.recipient_shell_index
                or abs(atom.recipient_shell_scale - ref.recipient_shell_scale) > 5.0e-12 * ref.recipient_shell_scale
                or max(abs(a-b) for a,b in zip(sorted(atom.interaction_parent_radii), sorted(ref.interaction_parent_radii))) > 5.0e-12 * ref.recipient_shell_scale
            ):
                raise ValueError("same recipient donor sidecars disagree on physical edge geometry")
            if atom.later_hahn_used or atom.own_shell_causal_reweighting_used or atom.creates_new_event_depth:
                raise ValueError("upstream pure-UV atom already violated causal provenance")
            if scale <= 0.0:
                raise AssertionError("invalid native coalescing scale")
        out.append(
            PureUVRecipientSubmeasure(
                boundary=ref.boundary,
                recipient_closed_mode_index=recipient,
                donor_closed_mode_indices=donors,
                canonical_positive_submeasure_mass=math.fsum(atom.physical_work_mass for atom in group),
                recipient_radius=ref.recipient_radius,
                recipient_shell_scale=ref.recipient_shell_scale,
                interaction_parent_radii=ref.interaction_parent_radii,
                common_unit_scale=ref.common_unit_scale,
            )
        )
    return tuple(out)


def pure_uv_first_shell_law(split: HardTailUpwardSupplySplit) -> PureUVFirstShellLaw:
    pure = tuple(atom for atom in split.atoms if atom.pure_uv_hh_by_support)
    if not pure:
        raise ValueError("true-upward split carries no pure-UV positive submeasure")
    charges = coalesce_pure_uv_recipient_submeasures(pure, native_work_mass_scale=split.native_work_scale)
    total = math.fsum(charge.canonical_positive_submeasure_mass for charge in charges)
    expected = split.pure_uv_hh_physical_work
    native = max(split.native_work_scale, total, expected, 1.0e-300)
    residual = abs(total - expected) / native
    M = 2.0 * split.boundary
    return PureUVFirstShellLaw(
        boundary=split.boundary,
        recipient_shell_scale=M,
        recipient_submeasures=charges,
        total_canonical_positive_mass=total,
        total_common_unit_work=split.boundary * total,
        native_work_mass_scale=split.native_work_scale,
        coalescing_native_residual=residual,
    )


def pure_uv_natural_window_common_work_upper(
    window_peak_child_mass: float,
    parent_frequency: float,
    global_energy: float,
    scaled_lifetime: float,
) -> float:
    """Clean capacity for the canonical pure-UV recipient submeasure.

    Both physical parents are strictly above M/4, so every admissible resolved
    cutoff vanishes on them and h=u exactly.  Before Young, the already-canonical
    recipient submeasure is dominated by physical edge total variation.  The
    certified helicity-edge bridge has exact variation/clean-Young ratio
    sqrt(2) C_F < 1.  With automatic parent corridor R=3/2,

      N W_win <= (9/2) c sqrt(pi) N E_global sqrt(mu_win).
    """

    mu = _finite_nonnegative(window_peak_child_mass, "window peak child mass")
    N = _finite_positive(parent_frequency, "parent frequency")
    E = _finite_positive(global_energy, "global energy")
    c = _finite_positive(scaled_lifetime, "scaled lifetime")
    bridge = canonical_hh_edge_total_variation_young_bridge()
    if not bridge.canonical_positive_submeasure_dominated_by_edge_variation:
        raise AssertionError("canonical pure-UV positive submeasure lost edge-total-variation domination")
    if not bridge.edge_variation_dominated_by_capacity_measure:
        raise AssertionError("physical edge total variation escaped capacity measure")
    if not bridge.edge_variation_to_clean_young_ratio < 1.0:
        raise AssertionError("clean Young constant no longer dominates helicity-edge total variation")
    return PURE_UV_PARENT_UPPER_RATIO * ultraviolet_hh_work_constant() * c * N * E * math.sqrt(mu)


def pure_uv_direct_natural_window_reentry(
    law: PureUVFirstShellLaw,
    *,
    required_pure_common_work_lower: float,
    global_energy: float,
    scaled_lifetime: float,
    viscosity: float,
    maximum_window_common_work: float,
    window_length: float,
    window_peak_child_mass: float,
) -> dict[str, object]:
    """Bind the fixed first-shell pure-UV cause directly to sliding natural time."""

    lower = _finite_positive(required_pure_common_work_lower, "required pure-UV common-work lower")
    E = _finite_positive(global_energy, "global energy")
    c = _finite_positive(scaled_lifetime, "scaled lifetime")
    nu = _finite_positive(viscosity, "viscosity")
    Ww = _finite_positive(maximum_window_common_work, "maximum pure-UV window work")
    Tw = _finite_positive(window_length, "pure-UV natural-window length")
    mu = _finite_positive(window_peak_child_mass, "pure-UV window peak child mass")
    H = law.total_common_unit_work
    if H + 6.0e-13 * max(H, lower) < lower:
        raise ValueError("actual canonical pure-UV law is below its supplied owner lower")
    N = law.boundary
    geometry = natural_window_geometry(N, PURE_UV_SHELL_INDEX, c)
    M = float(geometry["selected_shell_frequency"])
    if abs(M - law.recipient_shell_scale) > 5.0e-12 * M:
        raise AssertionError("natural-window shell disagrees with pure-UV M=2N geometry")
    Tnatural = float(geometry["selected_natural_window"])
    if abs(Tw - Tnatural) > 8.0e-13 * max(Tw, Tnatural):
        raise ValueError("pure-UV temporal window is not the exact first-shell natural window c(2N)^-2")
    temporal = temporal_concentration_statistics(H, Ww, Tw)
    p_t = float(temporal["p_time"])
    capacity = pure_uv_natural_window_common_work_upper(mu, N, E, c)
    if Ww > capacity + 8.0e-13 * max(Ww, capacity):
        raise ValueError("pure-UV window work exceeds clean physical edge-variation Young capacity")
    weighted = math.sqrt(mu) / p_t
    clean = lower / (
        PURE_UV_PARENT_UPPER_RATIO * ultraviolet_hh_work_constant() * c * N * E
    )
    margin = weighted - clean
    if margin < -1.0e-11 * max(weighted, clean):
        raise AssertionError("pure-UV direct natural-window mass lower failed")
    return {
        "selected_shell_level": PURE_UV_SHELL_INDEX,
        "selected_shell_frequency": M,
        "forward_scale_ratio": PURE_UV_SHELL_RATIO,
        "p_scale": 1.0,
        "H_inf_output_scale": 0.0,
        "selected_natural_window": Tnatural,
        "natural_time_ratio": float(geometry["natural_time_ratio"]),
        "p_time": p_t,
        "H_inf_time": float(temporal["H_inf_time"]),
        "total_pure_uv_common_work": H,
        "maximum_window_common_work": Ww,
        "window_peak_child_mass": mu,
        "natural_window_common_work_capacity": capacity,
        "weighted_sqrt_child_mass": weighted,
        "clean_weighted_sqrt_child_mass_lower": clean,
        "scale_time_tradeoff_margin": margin,
        "parent_upper_comparability_ratio": PURE_UV_PARENT_UPPER_RATIO,
        "h_equals_u_on_both_parents": True,
        "recipient_charges_single": True,
        "later_hahn_used": False,
        "output_scale_locality_theorem_used": False,
        "recipient_shell_reweighting_used": False,
        "packet_persistence_used": False,
        "time_partition_used": False,
        "next_owner": "generic_critical_shell_first_stop",
        "full_survivor_own_scale_service_lower": critical_shell_bounded_service_lower(mu, c, nu),
        "full_survivor_integrated_service_lower": critical_shell_integrated_service_lower(mu, c, nu),
        "full_survivor_service_is_conditional": True,
        "status": STATUS,
    }


def hard_tail_pure_uv_natural_window_reentry(
    split: HardTailUpwardSupplySplit,
    *,
    physical_tail_dissipation: float,
    viscosity: float,
    global_energy: float,
    scaled_lifetime: float,
    maximum_window_common_work: float,
    window_length: float,
    window_peak_child_mass: float,
) -> dict[str, object]:
    """Compose the certified true-upward owner lower with the pure-UV branch."""

    D = _finite_positive(physical_tail_dissipation, "physical tail dissipation")
    nu = _finite_positive(viscosity, "viscosity")
    owner_threshold = nu * D
    alternative = upward_owner_support_alternative(split, owner_threshold=owner_threshold)
    if not alternative.pure_uv_owner:
        raise ValueError("the supplied true-upward split does not realize the pure-UV support owner")
    law = pure_uv_first_shell_law(split)
    out = pure_uv_direct_natural_window_reentry(
        law,
        required_pure_common_work_lower=alternative.threshold_half,
        global_energy=global_energy,
        scaled_lifetime=scaled_lifetime,
        viscosity=nu,
        maximum_window_common_work=maximum_window_common_work,
        window_length=window_length,
        window_peak_child_mass=window_peak_child_mass,
    )
    expected_clean = nu * D / (
        9.0 * scaled_lifetime * math.sqrt(math.pi) * split.boundary * global_energy
    )
    if abs(float(out["clean_weighted_sqrt_child_mass_lower"]) - expected_clean) > 2.0e-12 * max(expected_clean, 1.0e-300):
        raise AssertionError("pure-UV hard-tail corollary lost the exact 1/(9 c sqrt(pi) N E) coefficient")
    return {
        **out,
        "physical_tail_dissipation": D,
        "viscosity": nu,
        "true_upward_owner_threshold": owner_threshold,
        "pure_uv_owner_lower": alternative.threshold_half,
        "hard_tail_clean_weighted_sqrt_child_mass_lower": expected_clean,
    }


def theorem_certificate() -> dict[str, object]:
    bridge = canonical_hh_edge_total_variation_young_bridge()
    return {
        "status": STATUS,
        "support_geometry": "pure-UV true-upward means both recipient parents>M/4; donor<=N then forces j=1, M=2N, donor in (M/4,M/2], and both parents<=3M/2",
        "smooth_cutoff": "every admissible resolved cutoff is supported in B_(M/4), hence vanishes on both strict pure-UV parents and h=u exactly there",
        "canonical_cause": "selected same-time donor sidecars are coalesced by recipient before capacity; the result is a positive submeasure of the existing canonical dW+ and is never re-Hahn-split",
        "scale_law": "every pure-UV atom is already on the unique first shell M=2N, so p_scale=1 and H_inf_output_scale=0 exactly",
        "variation_bridge": "canonical positive submeasure <= physical edge total variation <= dA; exact helicity-edge variation/clean-Young ratio sqrt(2) C_F is strictly below one",
        "natural_window": "with R=3/2 and h=u on both parents, N W_win <= (9/2)c sqrt(pi) N E_global sqrt(mu_win)",
        "hard_tail_corollary": "if pure-UV wins the support split, H_pure>=nu D_tail/2 and sqrt(mu_win)/p_time>=nu D_tail/(9 c sqrt(pi) N E_global)",
        "causal_unit": "all work remains in parent-tail N dW; M=2N is geometry only and never a causal reweighting",
        "not_used": "no output-shell Hahn/locality selection, K/S interface ownership, packet persistence, temporal deposit matching, or event-depth clock",
        "edge_variation_to_clean_young_ratio": bridge.edge_variation_to_clean_young_ratio,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class PureUVNaturalWindowStress:
    status: str
    samples: int
    resolved_triads: int
    pure_uv_atoms_checked: int
    recipient_causes_checked: int
    maximum_parent_to_shell_ratio: float
    minimum_donor_to_shell_ratio: float
    maximum_donor_to_shell_ratio: float
    interior_donor_ratio_atoms: int
    maximum_coalescing_native_residual: float
    minimum_natural_window_margin: float
    maximum_scale_probability_residual: float
    edge_variation_to_clean_young_ratio: float


def _random_closed_triad(rng: np.random.Generator):
    while True:
        k0 = rng.normal(size=3)
        k1 = rng.normal(size=3)
        k2 = -(k0 + k1)
        if min(float(np.linalg.norm(v)) for v in (k0, k1, k2)) > 0.05:
            break
    helicities = tuple(int(v) for v in rng.choice((-1, 1), size=3))
    amplitudes = tuple(complex(v) for v in (rng.normal(size=3) + 1j * rng.normal(size=3)))
    return register_closed_helical_triad(
        wavevectors=(k0, k1, k2), helicities=helicities, amplitudes=amplitudes
    )


def stress(samples: int = 75_000, seed: int = 2026081303) -> PureUVNaturalWindowStress:
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = np.random.default_rng(int(seed))
    resolved = pure_checked = causes_checked = interior_donor_ratio_atoms = 0
    max_parent = max_donor = max_coalesce = max_p_scale = 0.0
    min_donor = math.inf
    min_margin = math.inf
    for _ in range(count):
        triad = _random_closed_triad(rng)
        kernel = cyclic_triad_measure_kernel(
            triad, quotient_measure_mass=float(10.0 ** rng.uniform(-2.0, 2.0))
        )
        if not kernel.numerically_resolved_transport:
            continue
        resolved += 1
        for flow in kernel.atoms:
            rd = mode_radius(flow.donor_child_mode)
            rr = mode_radius(flow.recipient_child_mode)
            if rr <= rd:
                continue
            recipient_slot = triad.slot_for_closed_mode_index(flow.recipient_closed_mode_index)
            parent_radii = tuple(
                mode_radius(triad.modes[index])
                for index in recipient_slot.parent_closed_indices
            )
            if len(parent_radii) != 2:
                raise AssertionError("closed-triad recipient lost its two interaction parents")
            # Choose a real tail boundary N from the full admissible pure-UV
            # corridor, not merely N=|k_d|.  For the selected donor->recipient
            # flow we need
            #   |k_d| <= N < |k_r| <= 2N
            # and both recipient parents > N/2.  These are exactly
            #   N >= max(|k_d|, |k_r|/2),
            #   N < min(|k_r|, 2 min(parent radii)).
            lower = max(rd, 0.5 * rr)
            upper = min(rr, 2.0 * min(parent_radii))
            if upper <= lower * (1.0 + 2.0e-12):
                continue
            frac = float(rng.uniform(0.02, 0.98))
            boundary = lower + frac * (upper - lower)
            split = hard_tail_upward_supply_split(triad, kernel, boundary=boundary)
            pure = tuple(atom for atom in split.atoms if atom.pure_uv_hh_by_support)
            if not pure:
                continue
            law = pure_uv_first_shell_law(split)
            pure_checked += len(pure)
            causes_checked += len(law.recipient_submeasures)
            max_coalesce = max(max_coalesce, law.coalescing_native_residual)
            max_p_scale = max(max_p_scale, abs(law.p_scale - 1.0))
            for atom in pure:
                M = atom.recipient_shell_scale
                max_parent = max(max_parent, atom.comparable_parent_upper_ratio)
                ratio = atom.donor_radius / M
                min_donor = min(min_donor, ratio)
                max_donor = max(max_donor, ratio)
                interior_donor_ratio_atoms += int(ratio < 0.49)

            H = law.total_common_unit_work
            threshold = 0.9 * min(split.upward_common_unit_work, 2.0 * H)
            lower = 0.5 * threshold
            E = float(10.0 ** rng.uniform(-1.0, 1.0))
            c = float(10.0 ** rng.uniform(-1.0, 0.7))
            nu = float(10.0 ** rng.uniform(-2.0, -0.2))
            p_t = float(rng.uniform(0.25, 1.0))
            Ww = p_t * H
            C = PURE_UV_PARENT_UPPER_RATIO * ultraviolet_hh_work_constant() * c * law.boundary * E
            mu = (1.0 + float(rng.uniform(0.0, 1.0))) ** 2 * (Ww / C) ** 2
            Tw = c / ((2.0 * law.boundary) ** 2)
            out = pure_uv_direct_natural_window_reentry(
                law,
                required_pure_common_work_lower=lower,
                global_energy=E,
                scaled_lifetime=c,
                viscosity=nu,
                maximum_window_common_work=Ww,
                window_length=Tw,
                window_peak_child_mass=mu,
            )
            min_margin = min(min_margin, float(out["scale_time_tradeoff_margin"]))
    if resolved == 0 or pure_checked == 0 or causes_checked == 0:
        raise AssertionError("random physical stress did not exercise pure-UV upward recipient causes")
    if interior_donor_ratio_atoms == 0 or min_donor >= 0.49:
        raise AssertionError("pure-UV stress failed to enter the physical donor corridor below the N=|k_d| endpoint")
    if min_donor <= 0.25 or max_donor > 0.5 + 5.0e-12:
        raise AssertionError("pure-UV stress donor/shell ratios escaped the exact (1/4,1/2] corridor")
    bridge = canonical_hh_edge_total_variation_young_bridge()
    return PureUVNaturalWindowStress(
        status=STATUS,
        samples=count,
        resolved_triads=resolved,
        pure_uv_atoms_checked=pure_checked,
        recipient_causes_checked=causes_checked,
        maximum_parent_to_shell_ratio=max_parent,
        minimum_donor_to_shell_ratio=min_donor,
        maximum_donor_to_shell_ratio=max_donor,
        interior_donor_ratio_atoms=interior_donor_ratio_atoms,
        maximum_coalescing_native_residual=max_coalesce,
        minimum_natural_window_margin=min_margin,
        maximum_scale_probability_residual=max_p_scale,
        edge_variation_to_clean_young_ratio=bridge.edge_variation_to_clean_young_ratio,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=75_000)
    ap.add_argument("--seed", type=int, default=2026081303)
    ap.add_argument("--outdir", type=Path, default=Path("results-pure-uv-natural-window"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples, args.seed)
    payload = {"certificate": theorem_certificate(), "stress": asdict(out)}
    (args.outdir / "pure_uv_true_upward_natural_window.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    md = f"""# Pure-UV true-upward natural-window binding

Status: **{STATUS}**.

The true-upward pure-UV branch is already the first physical dyadic shell: `M=2N`.  Both recipient interaction parents lie strictly above `M/4` and at or below `3M/2`, so every admissible resolved cutoff vanishes on both parents and `h=u` there.  Same-time donor sidecars are coalesced by recipient before capacity is read; no canonical `dW+` mass is cloned or re-Hahn-split.

Because there is only one possible output shell, `p_scale=1` exactly.  The canonical recipient submeasure passes first through physical edge total variation and then the clean Young envelope.  With `R=3/2`, one natural window obeys `N W_win <= (9/2)c sqrt(pi) N E_global sqrt(mu_win)`.  On the hard-tail pure-owner branch this gives `sqrt(mu_win)/p_time >= nu D_tail/(9 c sqrt(pi) N E_global)`.

Stress: `{out.samples}` random physical closed-triad draws
- resolved triads: `{out.resolved_triads}`
- pure-UV donor atoms / coalesced recipient causes: `{out.pure_uv_atoms_checked}` / `{out.recipient_causes_checked}`
- maximum parent/shell ratio: `{out.maximum_parent_to_shell_ratio:.12g}`
- donor/shell range: `{out.minimum_donor_to_shell_ratio:.12g}` / `{out.maximum_donor_to_shell_ratio:.12g}`
- interior donor-ratio atoms (`|k_d|/M<0.49`): `{out.interior_donor_ratio_atoms}`
- maximum coalescing residual: `{out.maximum_coalescing_native_residual:.3e}`
- minimum natural-window margin: `{out.minimum_natural_window_margin:.3e}`
- p_scale residual: `{out.maximum_scale_probability_residual:.3e}`
- edge-total-variation / clean-Young ratio: `{out.edge_variation_to_clean_young_ratio:.12g}`

No output-scale locality theorem, downstream Hahn split, packet persistence, temporal deposit matching, event-depth clock or Navier--Stokes global-regularity claim is used.
"""
    (args.outdir / "summary.md").write_text(md)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
