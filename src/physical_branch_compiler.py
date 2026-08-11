from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np


class PhysicalCurrency(str, Enum):
    """Fine physical destinations used by the branch compiler.

    The names are intentionally physical rather than theorem-module names.
    """

    MULTIPLICATIVE_TRANSFER = "multiplicative_transfer_loss"
    KELVIN_FLAT_EROSION = "kelvin_flat_barycentric_erosion"
    COHERENT_STICKY_ANCESTRY = "coherent_sticky_ancestry"
    RENYI_REUSE = "renyi_shannon_reuse_entropy_or_pair_mass"
    SIDEBAND_CURVATURE = "h1_h3_sideband_curvature_cost"
    RESOLVED_SOURCE_SGS = "resolved_source_or_sgs_service"
    CRITICAL_DISSIPATION = "critical_dissipation"
    UNIFORM_GLOBAL_RESET = "uniform_globally_bounded_resource_reset"
    INITIAL_BOUNDARY = "initial_boundary_coherent_fourier_tail"
    XI = "xi"


class PhysicalCause(str, Enum):
    """Causal roots, after quotienting duplicate theorem manifestations."""

    TRANSFER_WORK_LOSS = "transfer_work_loss"
    RESOLVED_SOURCE = "resolved_source"
    HIGH_STRAIN_DISSIPATION = "high_strain_dissipation"
    MATERIAL_RELINK = "material_relink"
    CAUSAL_REUSE = "causal_reuse"
    INTRINSIC_SIDEBAND = "intrinsic_sideband"
    NEW_COHERENT_ANCESTRY = "new_coherent_ancestry"
    UNIFORM_GLOBAL_RESOURCE = "uniform_global_resource"
    INITIAL_BOUNDARY = "initial_boundary"


CAUSE_TO_CURRENCY: Mapping[PhysicalCause, PhysicalCurrency] = {
    PhysicalCause.TRANSFER_WORK_LOSS: PhysicalCurrency.MULTIPLICATIVE_TRANSFER,
    PhysicalCause.RESOLVED_SOURCE: PhysicalCurrency.RESOLVED_SOURCE_SGS,
    PhysicalCause.HIGH_STRAIN_DISSIPATION: PhysicalCurrency.CRITICAL_DISSIPATION,
    PhysicalCause.MATERIAL_RELINK: PhysicalCurrency.COHERENT_STICKY_ANCESTRY,
    PhysicalCause.CAUSAL_REUSE: PhysicalCurrency.RENYI_REUSE,
    PhysicalCause.INTRINSIC_SIDEBAND: PhysicalCurrency.SIDEBAND_CURVATURE,
    PhysicalCause.NEW_COHERENT_ANCESTRY: PhysicalCurrency.COHERENT_STICKY_ANCESTRY,
    PhysicalCause.UNIFORM_GLOBAL_RESOURCE: PhysicalCurrency.UNIFORM_GLOBAL_RESET,
    PhysicalCause.INITIAL_BOUNDARY: PhysicalCurrency.INITIAL_BOUNDARY,
}


class MasterDisposition(str, Enum):
    TRANSFER_COST = "multiplicative_transfer_cost"
    FLAT = "flat_potential_erosion"
    RECURSE_CRITICAL = "sticky_or_scale_critical_recurse"
    ADDITIVE_RESET = "genuinely_uniform_additive_reset"
    BOUNDARY = "absorbing_initial_boundary"
    XI = "summable_xi"


CURRENCY_TO_MASTER: Mapping[PhysicalCurrency, MasterDisposition] = {
    PhysicalCurrency.MULTIPLICATIVE_TRANSFER: MasterDisposition.TRANSFER_COST,
    PhysicalCurrency.KELVIN_FLAT_EROSION: MasterDisposition.FLAT,
    PhysicalCurrency.COHERENT_STICKY_ANCESTRY: MasterDisposition.RECURSE_CRITICAL,
    PhysicalCurrency.RENYI_REUSE: MasterDisposition.TRANSFER_COST,
    PhysicalCurrency.SIDEBAND_CURVATURE: MasterDisposition.TRANSFER_COST,
    PhysicalCurrency.RESOLVED_SOURCE_SGS: MasterDisposition.RECURSE_CRITICAL,
    PhysicalCurrency.CRITICAL_DISSIPATION: MasterDisposition.RECURSE_CRITICAL,
    PhysicalCurrency.UNIFORM_GLOBAL_RESET: MasterDisposition.ADDITIVE_RESET,
    PhysicalCurrency.INITIAL_BOUNDARY: MasterDisposition.BOUNDARY,
    PhysicalCurrency.XI: MasterDisposition.XI,
}


class DoubleChargeRelation(str, Enum):
    MUTUALLY_EXCLUSIVE = "mutually_exclusive"
    DOWNSTREAM_NO_DOUBLE = "downstream_do_not_charge_separately"
    PRIMARY_DIAGNOSTIC = "coexist_choose_one_primary_other_diagnostic"
    INDEPENDENT_COCHARGE = "independent_resources_may_cocharge"


@dataclass(frozen=True)
class CauseHit:
    """One manifestation of a causal root at a physical time.

    Several theorem observables that arise from the same physical cause must be
    represented with the same ``cause``.  The compiler first quotients by cause,
    so adding another theorem name cannot create another charge.
    """

    time: float
    cause: PhysicalCause
    weight: float = 1.0
    witness: str = ""


@dataclass(frozen=True)
class BlockWitness:
    """Retained efficient block state seen by the compiler.

    ``fixed_transfer_loss`` is the transfer-selection stopping gate.  If it is
    true, no later evolution observable is needed to pay the retained block.
    Otherwise the first causal root is selected.  If no causal root is hit, a
    certified Kelvin-flat block is the only admissible zero-cost continuation.
    """

    fixed_transfer_loss: bool
    kelvin_flat_certified: bool
    hits: tuple[CauseHit, ...] = ()


@dataclass(frozen=True)
class TransferPartition:
    total_mass: float
    xi_mass: float
    retained_mass: float
    currency_mass: dict[str, float]
    first_time: float | None
    tied_causes: tuple[str, ...]


@dataclass(frozen=True)
class UniformResourceCertificate:
    threshold: float
    total_budget: float
    scale_independent_threshold: bool
    globally_bounded_resource: bool

    def valid(self) -> bool:
        return (
            self.threshold > 0
            and self.total_budget >= 0
            and self.scale_independent_threshold
            and self.globally_bounded_resource
        )


@dataclass(frozen=True)
class DuhamelTransferKernelCertificate:
    """The deliberately unresolved continuum bridge.

    The adjoint positive measure dGamma is an amplitude-generation law.  It is
    not automatically the same as the physical positive child-transfer law.
    """

    normalized_positive_kernel: bool
    same_material_labels: bool
    identified_with_physical_transfer_law: bool

    def master_ready(self) -> bool:
        return (
            self.normalized_positive_kernel
            and self.same_material_labels
            and self.identified_with_physical_transfer_law
        )


@dataclass(frozen=True)
class PhysicalEnergyCausalBridgeCertificate:
    """Preferred master-facing bridge after the raw dGamma=dT countermodel.

    Duhamel supplies the same-time quadratic parent-pair support.  The selected
    child energy balance supplies the positive physical HH work weights.
    """

    same_material_labels: bool
    exact_selected_coefficient_equation: bool
    positive_physical_hh_work_law: bool
    measure_agnostic_synchronization: bool
    residual_work_delegates_to_existing_ledger: bool

    def master_ready(self) -> bool:
        return (
            self.same_material_labels
            and self.exact_selected_coefficient_equation
            and self.positive_physical_hh_work_law
            and self.measure_agnostic_synchronization
            and self.residual_work_delegates_to_existing_ledger
        )


class UnresolvedCompilerBridge(RuntimeError):
    pass


def _validate_hits(hits: Iterable[CauseHit]) -> tuple[CauseHit, ...]:
    out = tuple(hits)
    for hit in out:
        if not math.isfinite(hit.time):
            raise ValueError("causal hit time must be finite")
        if not math.isfinite(hit.weight) or hit.weight <= 0:
            raise ValueError("causal hit weights must be finite and positive")
    return out


def _first_causal_split(witness: BlockWitness) -> tuple[dict[PhysicalCurrency, float], float | None, tuple[str, ...]]:
    """Partition one unit of retained transfer mass at the first causal defect.

    Exact ties between *independent causal roots* are not broken by a theorem-name
    priority.  Their positive stopping weights are normalized, giving a symmetric
    Radon--Nikodym split whose total is one.  Duplicate manifestations of the same
    causal root are combined before normalization.
    """

    if witness.fixed_transfer_loss:
        return {PhysicalCurrency.MULTIPLICATIVE_TRANSFER: 1.0}, None, (
            PhysicalCause.TRANSFER_WORK_LOSS.value,
        )

    hits = _validate_hits(witness.hits)
    if not hits:
        if witness.kelvin_flat_certified:
            return {PhysicalCurrency.KELVIN_FLAT_EROSION: 1.0}, None, ()
        raise UnresolvedCompilerBridge(
            "retained low-transfer-cost block has neither a first causal defect nor a certified Kelvin-flat continuation"
        )

    first_time = min(hit.time for hit in hits)
    tol = 64.0 * math.ulp(max(1.0, abs(first_time)))
    first_hits = [hit for hit in hits if abs(hit.time - first_time) <= tol]

    # t=0 is an absorbing physical boundary.  It is not an interior fresh event.
    if first_time <= tol and any(hit.cause is PhysicalCause.INITIAL_BOUNDARY for hit in first_hits):
        return {PhysicalCurrency.INITIAL_BOUNDARY: 1.0}, first_time, (
            PhysicalCause.INITIAL_BOUNDARY.value,
        )

    by_cause: dict[PhysicalCause, float] = {}
    for hit in first_hits:
        by_cause[hit.cause] = by_cause.get(hit.cause, 0.0) + hit.weight

    total_weight = sum(by_cause.values())
    if total_weight <= 0:
        raise AssertionError("positive causal quotient has zero total weight")

    split: dict[PhysicalCurrency, float] = {}
    for cause, weight in by_cause.items():
        currency = CAUSE_TO_CURRENCY[cause]
        split[currency] = split.get(currency, 0.0) + weight / total_weight

    if not math.isclose(sum(split.values()), 1.0, rel_tol=2e-15, abs_tol=2e-15):
        raise AssertionError("first-causal split failed to conserve unit transfer mass")

    return split, first_time, tuple(sorted(cause.value for cause in by_cause))


def compile_transfer_measure(
    *,
    total_mass: float,
    xi_mass: float,
    witness: BlockWitness,
) -> TransferPartition:
    """Single-charge partition of positive physical child-transfer mass.

    The selected cross-cell/interface measure is excised first and charged once
    to Xi.  The retained measure is then routed by the transfer gate / first
    causal stopping rule / Kelvin-flat continuation.
    """

    if not math.isfinite(total_mass) or total_mass <= 0:
        raise ValueError("total physical transfer mass must be positive and finite")
    if not math.isfinite(xi_mass) or xi_mass < 0 or xi_mass > total_mass:
        raise ValueError("Xi mass must lie in [0,total_mass]")

    retained = total_mass - xi_mass
    currency_mass: dict[str, float] = {}
    if xi_mass:
        currency_mass[PhysicalCurrency.XI.value] = xi_mass

    if retained == 0:
        return TransferPartition(total_mass, xi_mass, retained, currency_mass, None, ())

    split, first_time, tied_causes = _first_causal_split(witness)
    for currency, fraction in split.items():
        currency_mass[currency.value] = currency_mass.get(currency.value, 0.0) + retained * fraction

    if not math.isclose(sum(currency_mass.values()), total_mass, rel_tol=2e-14, abs_tol=2e-14 * total_mass):
        raise AssertionError("compiler double-counted or lost physical transfer mass")

    return TransferPartition(total_mass, xi_mass, retained, currency_mass, first_time, tied_causes)


def master_disposition(
    currency: PhysicalCurrency,
    *,
    uniform_certificate: UniformResourceCertificate | None = None,
) -> MasterDisposition:
    """Project one physical currency into the existing master telescope.

    A uniform reset is rejected unless its threshold is scale-independent in a
    genuinely globally bounded resource.  Critical N E and D_V therefore cannot
    enter the additive reset count through this function.
    """

    disposition = CURRENCY_TO_MASTER[currency]
    if disposition is MasterDisposition.ADDITIVE_RESET:
        if uniform_certificate is None or not uniform_certificate.valid():
            raise UnresolvedCompilerBridge(
                "uniform reset requested without a scale-independent threshold in a genuinely globally bounded resource"
            )
    return disposition


def require_duhamel_transfer_kernel(cert: DuhamelTransferKernelCertificate) -> None:
    """Legacy stronger bridge.

    Raw dGamma and physical energy transfer are generally different measures, so
    this equality certificate is no longer required by the preferred compiler
    route.  The function is retained only to make the stronger assertion explicit
    when a caller actually has such a special-case identification.
    """
    if not cert.master_ready():
        raise UnresolvedCompilerBridge(
            "raw adjoint dGamma has not been identified with physical transfer; use the preferred physical-energy causal bridge instead of assuming equality"
        )


def require_physical_energy_causal_bridge(cert: PhysicalEnergyCausalBridgeCertificate) -> None:
    if not cert.master_ready():
        raise UnresolvedCompilerBridge(
            "physical-energy causal recursion is not ready: require the exact selected coefficient equation, physical HH work weights, canonical labels, measure-agnostic synchronization, and residual delegation"
        )


def forbidden_double_charge_matrix() -> dict[str, str]:
    """Canonical pair decisions for the physically important overlap cases."""

    return {
        "H1/H3 dephasing | causing pressure/SGS/viscous source": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "high resolved strain | forced D_V lower bound": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "physical covariance deformation | strain/source that causes it": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "representative covariance change | physical covariance deformation": DoubleChargeRelation.PRIMARY_DIAGNOSTIC.value,
        "coherent relinking in retained graph | omitted cross-cell Xi": DoubleChargeRelation.MUTUALLY_EXCLUSIVE.value,
        "new coherent Moyal mass | fresh affine-radius certificate": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "pair rescue | reuse-cycle/Renyi currency": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "transfer deficit | backscatter/cancellation manifestation": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "physical cross-cell transfer | symbol-freezing approximation": DoubleChargeRelation.INDEPENDENT_COCHARGE.value,
        "Duhamel classified residual | underlying source/interface term": DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value,
        "initial-boundary termination | fresh interior packet": DoubleChargeRelation.MUTUALLY_EXCLUSIVE.value,
    }


@dataclass(frozen=True)
class CompilerStress:
    samples: int
    minimum_mass_margin: float
    maximum_mass_residual: float
    tie_samples: int
    flat_samples: int
    transfer_stop_samples: int
    boundary_samples: int


def stress(samples: int = 50_000, seed: int = 20260808) -> CompilerStress:
    rng = np.random.default_rng(seed)
    min_margin = float("inf")
    max_resid = 0.0
    tie_samples = flat_samples = transfer_samples = boundary_samples = 0
    causes = list(PhysicalCause)

    for _ in range(samples):
        total = float(10 ** rng.uniform(-8, 2))
        xi = total * float(rng.uniform(0.0, 0.2))
        fixed = bool(rng.random() < 0.2)
        flat = False
        hits: list[CauseHit] = []

        if not fixed:
            if rng.random() < 0.18:
                flat = True
                flat_samples += 1
            else:
                n = int(rng.integers(1, 6))
                times = rng.integers(1, 5, size=n).astype(float)
                if rng.random() < 0.08:
                    times[0] = 0.0
                    hits.append(CauseHit(0.0, PhysicalCause.INITIAL_BOUNDARY, float(rng.uniform(.1, 2.0)), "boundary"))
                    boundary_samples += 1
                for j in range(n):
                    cause = causes[int(rng.integers(0, len(causes) - 1))]
                    hits.append(CauseHit(float(times[j]), cause, float(rng.uniform(.01, 3.0)), f"w{j}"))
                first = min(h.time for h in hits)
                if sum(abs(h.time - first) < 1e-15 for h in hits) > 1:
                    tie_samples += 1
        else:
            transfer_samples += 1

        out = compile_transfer_measure(
            total_mass=total,
            xi_mass=xi,
            witness=BlockWitness(fixed, flat, tuple(hits)),
        )
        subtotal = sum(out.currency_mass.values())
        resid = abs(subtotal - total)
        max_resid = max(max_resid, resid)
        min_margin = min(min_margin, total + 2e-12 * max(1.0, total) - subtotal)
        if subtotal > total + 2e-12 * max(1.0, total):
            raise AssertionError("compiler created physical transfer mass")
        if any(v < -1e-15 for v in out.currency_mass.values()):
            raise AssertionError("compiler created negative branch mass")
        if fixed:
            non_xi = {k for k, v in out.currency_mass.items() if k != PhysicalCurrency.XI.value and v > 0}
            if non_xi != {PhysicalCurrency.MULTIPLICATIVE_TRANSFER.value}:
                raise AssertionError("transfer stopping gate did not absorb later branches")

    return CompilerStress(samples, min_margin, max_resid, tie_samples, flat_samples, transfer_samples, boundary_samples)


def theorem_certificate() -> dict[str, object]:
    return {
        "status": "EXACT_CAUSAL_QUOTIENT__JOINT_STOP_MASTER_PREFERRED__FINE_RN_PARTITION_OPTIONAL",
        "transfer_partition": "dT=dT_Xi+dT_ret exactly; the legacy fine compiler can further partition dT_ret by physical RN stopping weights when such common-unit weights are independently available",
        "stopping_rule": "fixed transfer loss is a stage-zero terminal gate; otherwise use the first physical cause set; no hit => registered generated continuation or certified Kelvin-flat erosion",
        "duplicate_rule": "theorem manifestations sharing one physical causal root are quotiented before any ledger projection",
        "tie_rule": "preferred master keeps an exact simultaneous first hit as one unsplit joint physical stop; no common-unit tie weights are required; RN splitting is optional fine-subledger metadata only",
        "boundary_rule": "t=0 is absorbing and cannot be renamed fresh interior ancestry",
        "master_rule": "a joint stop is terminal if it already contains a fixed transfer certificate, otherwise a valid genuinely global reset may terminate it, otherwise source/dissipation/relink/HH-regeneration remain one recursive-critical state",
        "duhamel_warning": "raw dGamma is neither the child-energy law nor the preferred parent-pair law; physical pair productivity is derived directly under dT by sharp Young plus KL positivity",
        "energy_bridge": "on K<=1/30, E0<E1/5 and classified residual positive work<E1/5 force actual positive HH work>=8E1/15; all causal entropy weights remain physical child-energy work",
        "survival_rule": "if registered continuation retains fraction q of generated dT, its physical logarithmic productivity is exactly q times the unconditioned Lambda; q>=1/2 costs at most log2 at that layer",
        "continuum_status": "outer moving roles, event-role registration, interface provenance, physical pair weighting, survival algebra and exact-tie master projection are supplied; the remaining PDE bridge is measurable first-hit cause-set extraction on every recursive smooth-SGS block",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-physical-branch-compiler"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    payload = {
        "certificate": cert,
        "double_charge_matrix": forbidden_double_charge_matrix(),
        "stress": asdict(out),
    }
    (args.outdir / "physical_branch_compiler.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    md = f"""# Single-charge physical branch compiler\n\nStatus: **{cert['status']}**.\n\nThe compiler first removes the selected physical cross-cell/interface transfer measure into `Xi`.  On the retained positive child-transfer measure, a definite transfer loss is an absorbing gate.  Otherwise theorem manifestations are quotiented by their physical causal root and the first causal root is charged.  If independent roots hit at the same first time, their positive stopping weights partition the transfer measure; no lexicographic theorem priority is used.  If there is no hit, the block is admissible only with a certified Kelvin-flat continuation.\n\nRaw `dGamma` from the adjoint Kelvin--Duhamel gate is kept as an amplitude-generation diagnostic, not forced to equal energy transfer.  The preferred master-facing bridge uses Duhamel only for same-time quadratic parent-pair support and uses the actual positive high--high child-energy work as the causal weights; the half-slab/parabolic synchronization geometry is measure agnostic.\n\nCritical `N E`, resolved-source service and `D_V` remain recursive/scale-critical currencies.  They are never promoted to the finite additive reset count without an independent scale-independent threshold in a genuinely globally bounded resource.\n\nStress: `{out.samples}` synthetic causal block states\n- maximum transfer-mass residual: `{out.maximum_mass_residual:.3e}`\n- exact/quantized first-time tie samples: `{out.tie_samples}`\n- Kelvin-flat samples: `{out.flat_samples}`\n- transfer-stop samples: `{out.transfer_stop_samples}`\n- initial-boundary samples: `{out.boundary_samples}`\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
