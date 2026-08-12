from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.continuum_helical_edge_measure_registration import HelicalModeIdentity
from src.cyclic_helical_triad_donor_kernel import (
    ClosedHelicalTriadRegistration,
    CyclicTriadMeasureKernel,
    cyclic_triad_measure_kernel,
    register_closed_helical_triad,
    signed_good_integer_triad,
)
from src.helical import stable_norm3
from src.helical_mode_set_energy_continuity import flow_atoms_from_cyclic_kernel
from src.radial_spectral_crossing_layer_cake import (
    equiradial_physical_transfer_triad,
    mode_radius,
    radial_exterior_balance,
)

STATUS = (
    "EXACT_HARD_TAIL_TRUE_UPWARD_SUPPLY__RADIAL_BOUNDARY_STOCK_FLOW_VISCOSITY__"
    "CANONICAL_LOW_TO_HIGH_SUPPLY_IN_COMMON_N_DW_UNIT__PURE_UV_HH_ONLY_FIRST_DYADIC_SHELL__"
    "DEEP_UPWARD_HAS_RESOLVED_SCALE_PARENT_CONTACT__NO_INTERNAL_CIRCULATION_REENTRY_OR_INTERFACE_OVERCLAIM"
)


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


def _native_residual(actual: float, expected: float, scale: float) -> float:
    a, b, s = float(actual), float(expected), float(scale)
    if not all(math.isfinite(v) for v in (a, b, s)) or s < 0.0:
        raise ValueError("finite nonnegative native comparison scale required")
    gap = abs(a - b)
    if s == 0.0:
        return 0.0 if gap == 0.0 else math.inf
    return gap / s


def _negated_mode(mode: HelicalModeIdentity) -> HelicalModeIdentity:
    return HelicalModeIdentity(tuple(-float(v) for v in mode.wavevector), mode.helicity)


def _recipient_shell(boundary: float, recipient_radius: float) -> tuple[int, float]:
    N = _finite_positive(boundary, "tail boundary")
    rho = _finite_positive(recipient_radius, "recipient radius")
    if not rho > N:
        raise ValueError("recipient shell is defined only for true upward tail recipients")
    j = 1
    M = 2.0 * N
    while rho > M:
        j += 1
        M *= 2.0
    return j, M


@dataclass(frozen=True)
class UpwardSupplyAtom:
    boundary: float
    donor_closed_mode_index: int
    recipient_closed_mode_index: int
    donor_mode: HelicalModeIdentity
    recipient_mode: HelicalModeIdentity
    physical_work_mass: float
    donor_radius: float
    recipient_radius: float
    recipient_shell_index: int
    recipient_shell_scale: float
    interaction_parent_radii: tuple[float, float]
    donor_is_interaction_parent: bool
    pure_uv_hh_by_support: bool
    resolved_scale_parent_contact: bool
    first_dyadic_shell: bool
    deep_upward_shell: bool
    comparable_parent_upper_ratio: float
    common_unit_scale: float
    later_hahn_used: bool = False
    own_shell_causal_reweighting_used: bool = False
    resolved_contact_declared_interface_owner: bool = False
    creates_new_event_depth: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        m = _finite_positive(self.physical_work_mass, "upward physical work mass")
        rd = _finite_positive(self.donor_radius, "donor radius")
        rr = _finite_positive(self.recipient_radius, "recipient radius")
        M = _finite_positive(self.recipient_shell_scale, "recipient shell scale")
        _finite_positive(self.common_unit_scale, "common work-unit scale")
        if self.donor_closed_mode_index == self.recipient_closed_mode_index:
            raise ValueError("upward supply donor and recipient must be distinct cyclic roots")
        if not (rd <= N < rr):
            raise AssertionError("upward supply atom does not actually cross the physical tail boundary")
        if self.recipient_shell_index < 1 or not (0.5 * M < rr <= M):
            raise AssertionError("recipient dyadic shell does not contain the actual recipient mode")
        if abs(M - (2.0 ** self.recipient_shell_index) * N) > 5.0e-13 * max(M, N):
            raise AssertionError("recipient shell scale changed from the boundary-anchored dyadic law")
        if len(self.interaction_parent_radii) != 2 or min(self.interaction_parent_radii) <= 0.0:
            raise ValueError("two positive interaction-parent radii required")
        if not self.donor_is_interaction_parent:
            raise AssertionError("cyclic energy donor must be one of the recipient edge's two interaction parents")
        support_pure = min(self.interaction_parent_radii) > 0.25 * M
        if self.pure_uv_hh_by_support != support_pure:
            raise AssertionError("pure-UV HH label changed from exact parent-frequency support")
        if self.resolved_scale_parent_contact == self.pure_uv_hh_by_support:
            raise AssertionError("upward supply must partition into pure-UV HH or resolved-scale parent contact")
        if self.first_dyadic_shell != (self.recipient_shell_index == 1):
            raise AssertionError("first-shell label changed")
        if self.deep_upward_shell != (self.recipient_shell_index >= 2):
            raise AssertionError("deep-shell label changed")
        if self.pure_uv_hh_by_support:
            if not self.first_dyadic_shell:
                raise AssertionError("true upward pure-UV HH supply cannot enter beyond the first dyadic shell")
            # M=2N, donor lies in (M/4,M/2].  The other parent obeys
            # |p_other| <= |child|+|donor| <= M+M/2 = 3M/2.
            if not (0.25 * M < rd <= 0.5 * M):
                raise AssertionError("pure upward donor left the first-shell parent corridor")
            if max(self.interaction_parent_radii) > 1.5 * M + 5.0e-12 * M:
                raise AssertionError("triad closure lost automatic comparable-parent support")
            if self.comparable_parent_upper_ratio > 1.5 + 5.0e-12:
                raise AssertionError("pure upward HH parents are not automatically comparable")
        if self.deep_upward_shell:
            if not self.resolved_scale_parent_contact:
                raise AssertionError("deep upward crossing must touch the recipient's resolved-scale parent region")
            if rd > 0.25 * M + 5.0e-12 * M:
                raise AssertionError("deep upward donor unexpectedly lies above the child-scale quarter boundary")
        if self.common_unit_scale != N:
            raise ValueError("high-tail upward supply must remain in the parent-tail common N dW unit")
        if (
            self.later_hahn_used
            or self.own_shell_causal_reweighting_used
            or self.resolved_contact_declared_interface_owner
            or self.creates_new_event_depth
        ):
            raise ValueError(
                "upward supply may not re-Hahn, reweight by recipient shell, overclaim resolved contact as interface ownership, or add event depth"
            )
        if m <= 0.0:
            raise AssertionError("upward supply atom lost positive canonical work")

    @property
    def common_unit_work_mass(self) -> float:
        return self.common_unit_scale * self.physical_work_mass


@dataclass(frozen=True)
class HardTailUpwardSupplySplit:
    boundary: float
    upward_physical_work: float
    upward_common_unit_work: float
    pure_uv_hh_physical_work: float
    resolved_contact_physical_work: float
    first_shell_physical_work: float
    deep_shell_physical_work: float
    atoms: tuple[UpwardSupplyAtom, ...]
    native_work_scale: float
    upward_partition_native_residual: float
    radial_upward_binding_native_residual: float
    common_unit_partition_native_residual: float
    pure_uv_nonfirst_shell_atoms: int
    deep_pure_uv_atoms: int
    resolved_contact_is_interface_owner: bool = False
    internal_high_high_included_as_supply: bool = False
    later_hahn_used: bool = False
    recipient_shell_reweighting_used: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        for name, value in (
            ("upward work", self.upward_physical_work),
            ("upward common-unit work", self.upward_common_unit_work),
            ("pure-UV HH work", self.pure_uv_hh_physical_work),
            ("resolved-contact work", self.resolved_contact_physical_work),
            ("first-shell work", self.first_shell_physical_work),
            ("deep-shell work", self.deep_shell_physical_work),
            ("native work scale", self.native_work_scale),
            ("upward partition residual", self.upward_partition_native_residual),
            ("radial upward binding residual", self.radial_upward_binding_native_residual),
            ("common-unit partition residual", self.common_unit_partition_native_residual),
        ):
            _finite_nonnegative(value, name)
        if self.native_work_scale <= 0.0:
            raise ValueError("positive native supply work scale required")
        if self.upward_physical_work <= 0.0 or not self.atoms:
            raise ValueError("nonempty true upward supply required")
        if max(
            self.upward_partition_native_residual,
            self.radial_upward_binding_native_residual,
            self.common_unit_partition_native_residual,
        ) > 5.0e-10:
            raise AssertionError("upward support split left the native physical work scale")
        if self.pure_uv_nonfirst_shell_atoms or self.deep_pure_uv_atoms:
            raise AssertionError("pure-UV HH upward supply escaped its exact first-shell support law")
        if abs(self.upward_common_unit_work - N * self.upward_physical_work) > 5.0e-10 * N * self.native_work_scale:
            raise AssertionError("upward supply changed from the common parent-tail N dW unit")
        if (
            self.resolved_contact_is_interface_owner
            or self.internal_high_high_included_as_supply
            or self.later_hahn_used
            or self.recipient_shell_reweighting_used
        ):
            raise ValueError("support disintegration may not manufacture an owner, include circulation, re-Hahn, or reweight causality")

    @property
    def pure_uv_common_unit_work(self) -> float:
        return self.boundary * self.pure_uv_hh_physical_work

    @property
    def resolved_contact_common_unit_work(self) -> float:
        return self.boundary * self.resolved_contact_physical_work


def classify_upward_supply_atom(
    *,
    boundary: float,
    donor_closed_mode_index: int,
    recipient_closed_mode_index: int,
    donor_mode: HelicalModeIdentity,
    recipient_mode: HelicalModeIdentity,
    physical_work_mass: float,
    recipient_interaction_parents: Sequence[HelicalModeIdentity],
) -> UpwardSupplyAtom:
    """Classify one already-canonical cyclic low-to-high donor atom by support.

    This function performs no Hahn split and does not construct the donor flow.
    It only reads the physical radii and the recipient edge's two already-bound
    interaction parents.
    """
    N=_finite_positive(boundary,"tail boundary")
    rd=mode_radius(donor_mode)
    rr=mode_radius(recipient_mode)
    if not (rd <= N < rr):
        raise ValueError("classify_upward_supply_atom requires a true low-to-high boundary crossing")
    parents=tuple(recipient_interaction_parents)
    if len(parents)!=2:
        raise ValueError("recipient edge requires exactly two interaction parents")
    donor_parent=_negated_mode(donor_mode)
    donor_is_parent=donor_parent in parents
    if not donor_is_parent:
        raise AssertionError("cyclic donor root is not one recipient interaction parent")
    parent_radii=tuple(mode_radius(mode) for mode in parents)
    j,M=_recipient_shell(N,rr)
    pure=min(parent_radii)>0.25*M
    return UpwardSupplyAtom(
        boundary=N,
        donor_closed_mode_index=int(donor_closed_mode_index),
        recipient_closed_mode_index=int(recipient_closed_mode_index),
        donor_mode=donor_mode,
        recipient_mode=recipient_mode,
        physical_work_mass=float(physical_work_mass),
        donor_radius=rd,
        recipient_radius=rr,
        recipient_shell_index=j,
        recipient_shell_scale=M,
        interaction_parent_radii=(parent_radii[0],parent_radii[1]),
        donor_is_interaction_parent=donor_is_parent,
        pure_uv_hh_by_support=pure,
        resolved_scale_parent_contact=not pure,
        first_dyadic_shell=(j==1),
        deep_upward_shell=(j>=2),
        comparable_parent_upper_ratio=max(parent_radii)/M,
        common_unit_scale=N,
    )


def upward_supply_atoms_from_closed_triad(
    triad: ClosedHelicalTriadRegistration,
    kernel: CyclicTriadMeasureKernel,
    *,
    boundary: float,
) -> tuple[UpwardSupplyAtom, ...]:
    N = _finite_positive(boundary, "tail boundary")
    if not kernel.numerically_resolved_transport:
        raise ValueError("upward support classification refuses numerically unresolved cyclic signs")
    slots = {slot.closed_mode_index: slot for slot in triad.slots}
    out: list[UpwardSupplyAtom] = []
    for atom in kernel.atoms:
        rd = mode_radius(atom.donor_child_mode)
        rr = mode_radius(atom.recipient_child_mode)
        if not (rd <= N < rr):
            continue
        rslot = slots[atom.recipient_closed_mode_index]
        out.append(
            classify_upward_supply_atom(
                boundary=N,
                donor_closed_mode_index=atom.donor_closed_mode_index,
                recipient_closed_mode_index=atom.recipient_closed_mode_index,
                donor_mode=atom.donor_child_mode,
                recipient_mode=atom.recipient_child_mode,
                physical_work_mass=atom.physical_work_mass,
                recipient_interaction_parents=rslot.edge_identity.parents,
            )
        )
    return tuple(out)


def hard_tail_upward_supply_split(
    triad: ClosedHelicalTriadRegistration,
    kernel: CyclicTriadMeasureKernel,
    *,
    boundary: float,
) -> HardTailUpwardSupplySplit:
    N = _finite_positive(boundary, "tail boundary")
    atoms = upward_supply_atoms_from_closed_triad(triad, kernel, boundary=N)
    if not atoms:
        raise ValueError("selected closed triad carries no true upward crossing at this boundary")
    flow_atoms = flow_atoms_from_cyclic_kernel(kernel)
    radial = radial_exterior_balance(flow_atoms, radius=N)
    native = kernel.native_work_mass_scale
    up = math.fsum(a.physical_work_mass for a in atoms)
    pure = math.fsum(a.physical_work_mass for a in atoms if a.pure_uv_hh_by_support)
    contact = math.fsum(a.physical_work_mass for a in atoms if a.resolved_scale_parent_contact)
    first = math.fsum(a.physical_work_mass for a in atoms if a.first_dyadic_shell)
    deep = math.fsum(a.physical_work_mass for a in atoms if a.deep_upward_shell)
    common = N * up
    return HardTailUpwardSupplySplit(
        boundary=N,
        upward_physical_work=up,
        upward_common_unit_work=common,
        pure_uv_hh_physical_work=pure,
        resolved_contact_physical_work=contact,
        first_shell_physical_work=first,
        deep_shell_physical_work=deep,
        atoms=atoms,
        native_work_scale=native,
        upward_partition_native_residual=_native_residual(pure + contact, up, native),
        radial_upward_binding_native_residual=_native_residual(up, radial.upward_crossing_flow, native),
        common_unit_partition_native_residual=_native_residual(N * (pure + contact), common, N * native),
        pure_uv_nonfirst_shell_atoms=sum(a.pure_uv_hh_by_support and not a.first_dyadic_shell for a in atoms),
        deep_pure_uv_atoms=sum(a.deep_upward_shell and a.pure_uv_hh_by_support for a in atoms),
    )


@dataclass(frozen=True)
class TailStockUpwardSupplyCertificate:
    boundary: float
    viscosity: float
    initial_tail_energy: float
    final_tail_energy: float
    integrated_upward_work: float
    integrated_downward_work: float
    normalized_tail_dissipation: float
    inherited_common_energy: float
    upward_common_work: float
    downward_common_work: float
    final_common_energy: float
    viscous_common_loss: float
    owner_threshold: float
    inherited_owner: bool
    true_upward_owner: bool
    continuity_native_residual: float
    native_common_energy_scale: float
    positive_tail_work_used_instead_of_upward_crossing: bool = False
    internal_high_high_counted_as_supply: bool = False
    fifo_matching_used: bool = False
    own_shell_causal_reweighting_used: bool = False

    def __post_init__(self) -> None:
        N = _finite_positive(self.boundary, "tail boundary")
        nu = _finite_positive(self.viscosity, "viscosity")
        for name, value in (
            ("initial tail energy", self.initial_tail_energy),
            ("final tail energy", self.final_tail_energy),
            ("integrated upward work", self.integrated_upward_work),
            ("integrated downward work", self.integrated_downward_work),
            ("normalized tail dissipation", self.normalized_tail_dissipation),
            ("inherited common energy", self.inherited_common_energy),
            ("upward common work", self.upward_common_work),
            ("downward common work", self.downward_common_work),
            ("final common energy", self.final_common_energy),
            ("viscous common loss", self.viscous_common_loss),
            ("owner threshold", self.owner_threshold),
            ("continuity residual", self.continuity_native_residual),
            ("native common-energy scale", self.native_common_energy_scale),
        ):
            _finite_nonnegative(value, name)
        if self.native_common_energy_scale <= 0.0:
            raise ValueError("positive native common-energy scale required")
        expected = (
            N * self.initial_tail_energy,
            N * self.integrated_upward_work,
            N * self.integrated_downward_work,
            N * self.final_tail_energy,
            2.0 * nu * self.normalized_tail_dissipation,
            nu * self.normalized_tail_dissipation,
        )
        actual = (
            self.inherited_common_energy,
            self.upward_common_work,
            self.downward_common_work,
            self.final_common_energy,
            self.viscous_common_loss,
            self.owner_threshold,
        )
        if max(abs(a-b) for a,b in zip(actual, expected)) > 5.0e-12 * self.native_common_energy_scale:
            raise AssertionError("tail supply certificate changed physical units")
        if not (self.inherited_owner or self.true_upward_owner):
            raise AssertionError("exact tail stock/supply law lost the inherited-or-upward owner cover")
        if self.inherited_owner != (self.inherited_common_energy + 5.0e-12*self.native_common_energy_scale >= self.owner_threshold):
            raise AssertionError("inherited tail owner flag changed")
        if self.true_upward_owner != (self.upward_common_work + 5.0e-12*self.native_common_energy_scale >= self.owner_threshold):
            raise AssertionError("true upward owner flag changed")
        if (
            self.positive_tail_work_used_instead_of_upward_crossing
            or self.internal_high_high_counted_as_supply
            or self.fifo_matching_used
            or self.own_shell_causal_reweighting_used
        ):
            raise ValueError("tail supply theorem may not substitute gross positive tail work, circulation, temporal matching, or own-shell causal units")


def tail_stock_upward_supply_certificate(
    *,
    boundary: float,
    viscosity: float,
    initial_tail_energy: float,
    final_tail_energy: float,
    integrated_upward_work: float,
    integrated_downward_work: float,
    normalized_tail_dissipation: float,
    native_common_energy_scale: float | None = None,
) -> TailStockUpwardSupplyCertificate:
    N = _finite_positive(boundary, "tail boundary")
    nu = _finite_positive(viscosity, "viscosity")
    E0 = _finite_nonnegative(initial_tail_energy, "initial tail energy")
    E1 = _finite_nonnegative(final_tail_energy, "final tail energy")
    up = _finite_nonnegative(integrated_upward_work, "integrated upward work")
    down = _finite_nonnegative(integrated_downward_work, "integrated downward work")
    D = _finite_nonnegative(normalized_tail_dissipation, "normalized tail dissipation")
    inherited = N * E0
    up_common = N * up
    down_common = N * down
    final_common = N * E1
    viscous = 2.0 * nu * D
    threshold = nu * D
    lhs = final_common + viscous + down_common
    rhs = inherited + up_common
    scale = max(lhs, rhs, inherited + up_common, final_common + viscous + down_common, 1.0e-300) if native_common_energy_scale is None else _finite_positive(native_common_energy_scale, "native common-energy scale")
    residual = _native_residual(lhs, rhs, scale)
    # Exact NS identity implies rhs >= 2 nu D, since final/down are nonnegative.
    # A finite-precision caller must remain inside its explicitly supplied native
    # physical envelope; no universal PDE quadrature tolerance is imposed here.
    if rhs + 5.0e-12 * scale < viscous:
        raise AssertionError("tail stock plus true upward supply did not cover physical viscous loss")
    inherited_owner = inherited + 5.0e-12 * scale >= threshold
    upward_owner = up_common + 5.0e-12 * scale >= threshold
    return TailStockUpwardSupplyCertificate(
        boundary=N,
        viscosity=nu,
        initial_tail_energy=E0,
        final_tail_energy=E1,
        integrated_upward_work=up,
        integrated_downward_work=down,
        normalized_tail_dissipation=D,
        inherited_common_energy=inherited,
        upward_common_work=up_common,
        downward_common_work=down_common,
        final_common_energy=final_common,
        viscous_common_loss=viscous,
        owner_threshold=threshold,
        inherited_owner=inherited_owner,
        true_upward_owner=upward_owner,
        continuity_native_residual=residual,
        native_common_energy_scale=scale,
    )


@dataclass(frozen=True)
class UpwardOwnerSupportAlternative:
    owner_threshold: float
    pure_uv_common_work: float
    resolved_contact_common_work: float
    pure_uv_owner: bool
    resolved_contact_owner: bool
    threshold_half: float
    support_partition_native_residual: float
    native_common_work_scale: float
    resolved_contact_declared_interface_owner: bool = False

    def __post_init__(self) -> None:
        for name, value in (
            ("owner threshold", self.owner_threshold),
            ("pure-UV common work", self.pure_uv_common_work),
            ("resolved-contact common work", self.resolved_contact_common_work),
            ("half threshold", self.threshold_half),
            ("support partition residual", self.support_partition_native_residual),
            ("native common-work scale", self.native_common_work_scale),
        ):
            _finite_nonnegative(value, name)
        if self.native_common_work_scale <= 0.0:
            raise ValueError("positive native common-work scale required")
        if not (self.pure_uv_owner or self.resolved_contact_owner):
            raise AssertionError("upward owner lost pure-UV/resolved-contact support alternative")
        if self.resolved_contact_declared_interface_owner:
            raise ValueError("resolved-scale parent contact is not yet a proved interface owner")


def upward_owner_support_alternative(
    split: HardTailUpwardSupplySplit,
    *,
    owner_threshold: float,
) -> UpwardOwnerSupportAlternative:
    threshold = _finite_nonnegative(owner_threshold, "upward owner threshold")
    native = max(split.upward_common_unit_work, split.boundary * split.native_work_scale, 1.0e-300)
    if split.upward_common_unit_work + 5.0e-12*native < threshold:
        raise ValueError("the supplied radial submeasure is not the true-upward owner at this threshold")
    half = 0.5 * threshold
    pure = split.pure_uv_common_unit_work
    contact = split.resolved_contact_common_unit_work
    pure_owner = pure + 5.0e-12*native >= half
    contact_owner = contact + 5.0e-12*native >= half
    return UpwardOwnerSupportAlternative(
        owner_threshold=threshold,
        pure_uv_common_work=pure,
        resolved_contact_common_work=contact,
        pure_uv_owner=pure_owner,
        resolved_contact_owner=contact_owner,
        threshold_half=half,
        support_partition_native_residual=_native_residual(pure + contact, split.upward_common_unit_work, native),
        native_common_work_scale=native,
    )


def deep_upward_resolved_contact_fixture() -> tuple[ClosedHelicalTriadRegistration, CyclicTriadMeasureKernel, HardTailUpwardSupplySplit]:
    """Deterministic physical triad with direct upward supply into M=4N.

    The geometry is fixed: interaction modes of radii 1 and sqrt(3) feed a
    recipient of radius sqrt(6).  At N=1 this recipient belongs to the second
    dyadic shell M=4, so the radius-one energy donor touches M/4 exactly.  We
    search only the finite helicity/phase gauge needed to orient actual NS work;
    the support statement itself is independent of that orientation.
    """
    wavevectors = (
        np.asarray((-2.0,-1.0,-1.0)),
        np.asarray((1.0,0.0,0.0)),
        np.asarray((1.0,1.0,1.0)),
    )
    phases = (
        1.0+0.0j, -1.0+0.0j, 0.0+1.0j, 0.0-1.0j,
        complex(2.0**-0.5,2.0**-0.5), complex(2.0**-0.5,-2.0**-0.5),
        complex(-2.0**-0.5,2.0**-0.5), complex(-2.0**-0.5,-2.0**-0.5),
    )
    for s0 in (-1,1):
        for s1 in (-1,1):
            for s2 in (-1,1):
                for phase in phases:
                    for position in range(3):
                        amps=[1.0+0.0j,1.0+0.0j,1.0+0.0j]
                        amps[position]=phase
                        triad=register_closed_helical_triad(
                            wavevectors=wavevectors,
                            helicities=(s0,s1,s2),
                            amplitudes=tuple(amps),
                        )
                        kernel=cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
                        if not kernel.numerically_resolved_transport:
                            continue
                        try:
                            split=hard_tail_upward_supply_split(triad,kernel,boundary=1.0)
                        except ValueError:
                            continue
                        if any(a.deep_upward_shell for a in split.atoms):
                            return triad,kernel,split
    raise AssertionError("finite physical helicity/phase search found no deep upward resolved-contact triad")


def theorem_certificate() -> dict[str, object]:
    return {
        "status": STATUS,
        "tail_identity": "N E_>(t1) + 2 nu D_tail + N Phi_down = N E_>(t0) + N Phi_up",
        "owner_cover": "N E_>(t0) >= nu D_tail OR N Phi_up >= nu D_tail; exact ties remain joint",
        "upward_measure": "Phi_up is the already-certified low-to-high restriction of cyclic dW-/dW+ donor flow; high-to-high circulation is excluded before the owner split",
        "common_unit": "all upward recipient shells remain in the parent-tail common N dW unit; no M_j/N causal reweighting",
        "support_split": "each true upward atom is either pure-UV HH-by-support (both recipient interaction parents > M/4) or has resolved-scale parent contact (at least one <= M/4)",
        "pure_uv_rigidity": "because the energy donor is itself one recipient interaction parent and |k_d|<=N, pure-UV upward supply forces M=2N; both parents then lie in (M/4,3M/2] by triad closure",
        "deep_upward": "every upward recipient shell M>=4N has an interaction parent at or below M/4 and therefore is not in the pure-UV HH support region",
        "interface_scope": "resolved-scale parent contact is a support fact only; without a separate positive binding through the smooth cutoff it is not called an interface owner",
        "anti_theorem": "positive high-to-high/internal tail circulation may be arbitrarily large while true upward supply is zero",
        "later_hahn_used": False,
        "fifo_lifo_used": False,
        "claims_global_regularity": False,
    }


@dataclass(frozen=True)
class HardTailTrueUpwardSupplyStress:
    samples: int
    resolved_triads: int
    unresolved_triads: int
    upward_oriented_atoms_checked: int
    pure_uv_atoms: int
    resolved_contact_atoms: int
    deep_upward_atoms: int
    first_shell_contact_atoms: int
    maximum_recipient_shell_index: int
    worst_upward_partition_native_residual: float
    worst_radial_upward_binding_native_residual: float
    worst_common_unit_partition_native_residual: float
    maximum_pure_uv_parent_ratio_to_shell: float
    minimum_pure_uv_donor_ratio_to_shell: float
    maximum_pure_uv_donor_ratio_to_shell: float
    pure_uv_nonfirst_shell_violations: int
    deep_pure_uv_violations: int
    owner_certificates_checked: int
    inherited_only_owner_cases: int
    upward_only_owner_cases: int
    joint_owner_cases: int
    worst_tail_continuity_native_residual: float
    internal_circulation_base: float
    internal_circulation_scaled: float
    internal_circulation_upward_base: float
    internal_circulation_upward_scaled: float


def _random_closed_triad(rng: np.random.Generator) -> ClosedHelicalTriadRegistration:
    while True:
        k0 = rng.normal(size=3)
        k1 = rng.normal(size=3)
        k2 = -(k0 + k1)
        if min(stable_norm3(v) for v in (k0,k1,k2)) > 0.05:
            break
    helicities = tuple(int(v) for v in rng.choice((-1,1), size=3))
    amplitudes = tuple(complex(v) for v in (rng.normal(size=3) + 1j*rng.normal(size=3)))
    return register_closed_helical_triad(wavevectors=(k0,k1,k2), helicities=helicities, amplitudes=amplitudes)


def _scaled_triad(triad: ClosedHelicalTriadRegistration, factor: float) -> ClosedHelicalTriadRegistration:
    lam = _finite_positive(factor, "amplitude scale")
    return register_closed_helical_triad(
        wavevectors=tuple(mode.wavevector for mode in triad.modes),
        helicities=tuple(mode.helicity for mode in triad.modes),
        amplitudes=tuple(lam*a for a in triad.amplitudes),
    )


def stress(samples: int = 75_000, seed: int = 2026081207) -> HardTailTrueUpwardSupplyStress:
    count = int(samples)
    if count <= 0:
        raise ValueError("positive stress sample count required")
    rng = np.random.default_rng(int(seed))
    resolved = unresolved = 0
    checked = pure_count = contact_count = deep_count = first_contact = 0
    max_shell = 0
    worst_part = worst_radial = worst_common = 0.0
    max_pure_parent = 0.0
    min_pure_donor = math.inf
    max_pure_donor = 0.0
    pure_bad = deep_bad = 0

    for _ in range(count):
        triad = _random_closed_triad(rng)
        kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=float(10.0 ** rng.uniform(-2.0,2.0)))
        if not kernel.numerically_resolved_transport:
            unresolved += 1
            continue
        resolved += 1
        for atom in kernel.atoms:
            rd = mode_radius(atom.donor_child_mode)
            rr = mode_radius(atom.recipient_child_mode)
            if rr <= rd:
                continue
            # Use the physical donor radius itself as boundary.  This maximizes
            # the visible recipient-shell depth without changing the atom.
            N = rd
            split = hard_tail_upward_supply_split(triad, kernel, boundary=N)
            matching = [a for a in split.atoms if a.donor_closed_mode_index == atom.donor_closed_mode_index and a.recipient_closed_mode_index == atom.recipient_closed_mode_index]
            if not matching:
                raise AssertionError("oriented upward donor atom disappeared from its own boundary restriction")
            checked += len(matching)
            worst_part = max(worst_part, split.upward_partition_native_residual)
            worst_radial = max(worst_radial, split.radial_upward_binding_native_residual)
            worst_common = max(worst_common, split.common_unit_partition_native_residual)
            for a in matching:
                max_shell = max(max_shell, a.recipient_shell_index)
                if a.pure_uv_hh_by_support:
                    pure_count += 1
                    max_pure_parent = max(max_pure_parent, a.comparable_parent_upper_ratio)
                    ratio = a.donor_radius / a.recipient_shell_scale
                    min_pure_donor = min(min_pure_donor, ratio)
                    max_pure_donor = max(max_pure_donor, ratio)
                    pure_bad += int(not a.first_dyadic_shell)
                else:
                    contact_count += 1
                    if a.first_dyadic_shell:
                        first_contact += 1
                if a.deep_upward_shell:
                    deep_count += 1
                    deep_bad += int(a.pure_uv_hh_by_support)

    if checked == 0 or pure_count == 0 or contact_count == 0 or deep_count == 0:
        raise AssertionError("random physical stress did not exercise all upward support regimes")

    # Exact owner algebra is stressed independently of triad geometry.  It is a
    # statement about the already-certified radial tail continuity law.
    owner_checked = inherited_only = upward_only = joint = 0
    worst_tail = 0.0
    for j in range(max(2_000, count // 10)):
        N = float(10.0 ** rng.uniform(-1.0,1.0))
        nu = float(10.0 ** rng.uniform(-2.0,-0.2))
        D = float(10.0 ** rng.uniform(-2.0,1.0))
        threshold = nu*D
        # Alternate exact physical owner regimes while leaving nonnegative final/down stock.
        mode = j % 3
        if mode == 0:  # inherited only
            inherited = 1.55*threshold
            upward = 0.55*threshold
        elif mode == 1:  # upward only
            inherited = 0.55*threshold
            upward = 1.55*threshold
        else:  # joint
            inherited = 1.1*threshold
            upward = 1.1*threshold
        down = float(rng.uniform(0.0,0.15)*threshold)
        final = inherited + upward - 2.0*threshold - down
        if final < 0.0:
            down = 0.0
            final = inherited + upward - 2.0*threshold
        cert = tail_stock_upward_supply_certificate(
            boundary=N,
            viscosity=nu,
            initial_tail_energy=inherited/N,
            final_tail_energy=final/N,
            integrated_upward_work=upward/N,
            integrated_downward_work=down/N,
            normalized_tail_dissipation=D,
        )
        owner_checked += 1
        worst_tail = max(worst_tail, cert.continuity_native_residual)
        inherited_only += int(cert.inherited_owner and not cert.true_upward_owner)
        upward_only += int(cert.true_upward_owner and not cert.inherited_owner)
        joint += int(cert.inherited_owner and cert.true_upward_owner)

    # Internal high-frequency traffic is not supply.  Put a physical equiradial
    # triad wholly above a lower radial boundary and scale actual amplitudes.
    base_triad = equiradial_physical_transfer_triad()
    scaled_triad = _scaled_triad(base_triad, 10.0)
    base_kernel = cyclic_triad_measure_kernel(base_triad, quotient_measure_mass=1.0)
    scaled_kernel = cyclic_triad_measure_kernel(scaled_triad, quotient_measure_mass=1.0)
    min_radius = min(mode_radius(mode) for mode in base_triad.modes)
    lower = 0.5*min_radius
    base_radial = radial_exterior_balance(flow_atoms_from_cyclic_kernel(base_kernel), radius=lower)
    scaled_radial = radial_exterior_balance(flow_atoms_from_cyclic_kernel(scaled_kernel), radius=lower)
    if base_radial.upward_crossing_flow != 0.0 or scaled_radial.upward_crossing_flow != 0.0:
        raise AssertionError("full-high equiradial circulation was incorrectly counted as true upward supply")
    if not (base_radial.high_internal_flow > 0.0 and scaled_radial.high_internal_flow > 100.0*base_radial.high_internal_flow):
        raise AssertionError("physical amplitude scaling did not expose arbitrarily larger internal tail circulation")

    return HardTailTrueUpwardSupplyStress(
        samples=count,
        resolved_triads=resolved,
        unresolved_triads=unresolved,
        upward_oriented_atoms_checked=checked,
        pure_uv_atoms=pure_count,
        resolved_contact_atoms=contact_count,
        deep_upward_atoms=deep_count,
        first_shell_contact_atoms=first_contact,
        maximum_recipient_shell_index=max_shell,
        worst_upward_partition_native_residual=worst_part,
        worst_radial_upward_binding_native_residual=worst_radial,
        worst_common_unit_partition_native_residual=worst_common,
        maximum_pure_uv_parent_ratio_to_shell=max_pure_parent,
        minimum_pure_uv_donor_ratio_to_shell=min_pure_donor,
        maximum_pure_uv_donor_ratio_to_shell=max_pure_donor,
        pure_uv_nonfirst_shell_violations=pure_bad,
        deep_pure_uv_violations=deep_bad,
        owner_certificates_checked=owner_checked,
        inherited_only_owner_cases=inherited_only,
        upward_only_owner_cases=upward_only,
        joint_owner_cases=joint,
        worst_tail_continuity_native_residual=worst_tail,
        internal_circulation_base=base_radial.high_internal_flow,
        internal_circulation_scaled=scaled_radial.high_internal_flow,
        internal_circulation_upward_base=base_radial.upward_crossing_flow,
        internal_circulation_upward_scaled=scaled_radial.upward_crossing_flow,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=STATUS)
    parser.add_argument("--samples", type=int, default=75_000)
    parser.add_argument("--seed", type=int, default=2026081207)
    parser.add_argument("--outdir", type=Path, default=Path("results-hard-tail-true-upward-supply"))
    args = parser.parse_args()
    result = stress(args.samples, args.seed)
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "certificate.json").write_text(json.dumps(theorem_certificate(), indent=2, sort_keys=True) + "\n")
    (args.outdir / "summary.json").write_text(json.dumps(asdict(result), indent=2, sort_keys=True) + "\n")
    lines = [
        "# Hard-tail true upward supply",
        "",
        f"Status: `{STATUS}`",
        f"- samples: `{result.samples}`",
        f"- resolved/unresolved triads: `{result.resolved_triads}` / `{result.unresolved_triads}`",
        f"- upward oriented atoms checked: `{result.upward_oriented_atoms_checked}`",
        f"- pure-UV / resolved-contact atoms: `{result.pure_uv_atoms}` / `{result.resolved_contact_atoms}`",
        f"- deep upward atoms: `{result.deep_upward_atoms}`",
        f"- first-shell resolved-contact atoms: `{result.first_shell_contact_atoms}`",
        f"- maximum recipient shell index: `{result.maximum_recipient_shell_index}`",
        f"- worst upward/radial-binding/common-unit residuals: `{result.worst_upward_partition_native_residual:.3e}` / `{result.worst_radial_upward_binding_native_residual:.3e}` / `{result.worst_common_unit_partition_native_residual:.3e}`",
        f"- maximum pure-UV parent/shell ratio: `{result.maximum_pure_uv_parent_ratio_to_shell:.12g}`",
        f"- pure-UV donor/shell ratio range: `{result.minimum_pure_uv_donor_ratio_to_shell:.12g}` .. `{result.maximum_pure_uv_donor_ratio_to_shell:.12g}`",
        f"- pure-UV nonfirst/deep violations: `{result.pure_uv_nonfirst_shell_violations}` / `{result.deep_pure_uv_violations}`",
        f"- owner certificates checked: `{result.owner_certificates_checked}`",
        f"- inherited-only / upward-only / joint owner cases: `{result.inherited_only_owner_cases}` / `{result.upward_only_owner_cases}` / `{result.joint_owner_cases}`",
        f"- worst exact tail-continuity native residual: `{result.worst_tail_continuity_native_residual:.3e}`",
        f"- internal high-tail circulation base/scaled: `{result.internal_circulation_base:.12g}` / `{result.internal_circulation_scaled:.12g}`",
        f"- corresponding true upward supply base/scaled: `{result.internal_circulation_upward_base:.12g}` / `{result.internal_circulation_upward_scaled:.12g}`",
        "",
        "True upward supply is the radial low-to-high restriction of the already-canonical cyclic donor law. High-to-high circulation is real work but is not counted as tail supply. Pure-UV HH-by-support upward atoms can only feed the first dyadic shell M=2N and are automatically comparable; deeper upward atoms have resolved-scale parent contact but are not thereby declared interface owners.",
    ]
    (args.outdir / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
