import math

import pytest

from src.physical_branch_compiler import (
    BlockWitness,
    CauseHit,
    DoubleChargeRelation,
    DuhamelTransferKernelCertificate,
    MasterDisposition,
    PhysicalCause,
    PhysicalCurrency,
    UniformResourceCertificate,
    UnresolvedCompilerBridge,
    compile_transfer_measure,
    forbidden_double_charge_matrix,
    master_disposition,
    require_duhamel_transfer_kernel,
)


def test_transfer_gate_absorbs_later_manifestations_and_xi_is_disjoint():
    out = compile_transfer_measure(
        total_mass=10.0,
        xi_mass=1.0,
        witness=BlockWitness(
            fixed_transfer_loss=True,
            kelvin_flat_certified=False,
            hits=(CauseHit(1.0, PhysicalCause.RESOLVED_SOURCE),),
        ),
    )
    assert out.currency_mass[PhysicalCurrency.XI.value] == 1.0
    assert out.currency_mass[PhysicalCurrency.MULTIPLICATIVE_TRANSFER.value] == 9.0
    assert PhysicalCurrency.RESOLVED_SOURCE_SGS.value not in out.currency_mass
    assert math.isclose(sum(out.currency_mass.values()), 10.0)


def test_duplicate_theorem_manifestations_of_one_causal_root_do_not_double_charge():
    out = compile_transfer_measure(
        total_mass=1.0,
        xi_mass=0.0,
        witness=BlockWitness(
            fixed_transfer_loss=False,
            kelvin_flat_certified=False,
            hits=(
                CauseHit(2.0, PhysicalCause.RESOLVED_SOURCE, 1.0, "H1 dephasing source"),
                CauseHit(2.0, PhysicalCause.RESOLVED_SOURCE, 3.0, "pressure/SGS manifestation"),
            ),
        ),
    )
    assert out.currency_mass == {PhysicalCurrency.RESOLVED_SOURCE_SGS.value: 1.0}


def test_independent_exact_tie_partitions_mass_without_lexicographic_double_charge():
    out = compile_transfer_measure(
        total_mass=8.0,
        xi_mass=0.0,
        witness=BlockWitness(
            fixed_transfer_loss=False,
            kelvin_flat_certified=False,
            hits=(
                CauseHit(1.0, PhysicalCause.RESOLVED_SOURCE, 1.0),
                CauseHit(1.0, PhysicalCause.CAUSAL_REUSE, 3.0),
            ),
        ),
    )
    assert math.isclose(out.currency_mass[PhysicalCurrency.RESOLVED_SOURCE_SGS.value], 2.0)
    assert math.isclose(out.currency_mass[PhysicalCurrency.RENYI_REUSE.value], 6.0)
    assert math.isclose(sum(out.currency_mass.values()), 8.0)


def test_initial_boundary_is_absorbing_not_fresh_interior_mass():
    out = compile_transfer_measure(
        total_mass=3.0,
        xi_mass=0.0,
        witness=BlockWitness(
            fixed_transfer_loss=False,
            kelvin_flat_certified=False,
            hits=(
                CauseHit(0.0, PhysicalCause.INITIAL_BOUNDARY),
                CauseHit(0.0, PhysicalCause.NEW_COHERENT_ANCESTRY),
            ),
        ),
    )
    assert out.currency_mass == {PhysicalCurrency.INITIAL_BOUNDARY.value: 3.0}


def test_no_hit_requires_certified_flatness():
    flat = compile_transfer_measure(
        total_mass=2.0,
        xi_mass=0.0,
        witness=BlockWitness(False, True, ()),
    )
    assert flat.currency_mass == {PhysicalCurrency.KELVIN_FLAT_EROSION.value: 2.0}
    with pytest.raises(UnresolvedCompilerBridge):
        compile_transfer_measure(
            total_mass=2.0,
            xi_mass=0.0,
            witness=BlockWitness(False, False, ()),
        )


def test_critical_dissipation_is_recursive_not_additive_reset():
    assert master_disposition(PhysicalCurrency.CRITICAL_DISSIPATION) is MasterDisposition.RECURSE_CRITICAL
    assert master_disposition(PhysicalCurrency.RESOLVED_SOURCE_SGS) is MasterDisposition.RECURSE_CRITICAL


def test_uniform_reset_requires_real_global_resource_certificate():
    bad = UniformResourceCertificate(1.0, 10.0, False, True)
    with pytest.raises(UnresolvedCompilerBridge):
        master_disposition(PhysicalCurrency.UNIFORM_GLOBAL_RESET, uniform_certificate=bad)
    good = UniformResourceCertificate(1.0, 10.0, True, True)
    assert master_disposition(PhysicalCurrency.UNIFORM_GLOBAL_RESET, uniform_certificate=good) is MasterDisposition.ADDITIVE_RESET


def test_duhamel_kernel_is_explicitly_not_auto_identified_with_physical_transfer():
    with pytest.raises(UnresolvedCompilerBridge):
        require_duhamel_transfer_kernel(DuhamelTransferKernelCertificate(True, True, False))
    require_duhamel_transfer_kernel(DuhamelTransferKernelCertificate(True, True, True))


def test_forbidden_double_charge_matrix_contains_required_physical_pairs():
    matrix = forbidden_double_charge_matrix()
    assert matrix["H1/H3 dephasing | causing pressure/SGS/viscous source"] == DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value
    assert matrix["high resolved strain | forced D_V lower bound"] == DoubleChargeRelation.DOWNSTREAM_NO_DOUBLE.value
    assert matrix["coherent relinking in retained graph | omitted cross-cell Xi"] == DoubleChargeRelation.MUTUALLY_EXCLUSIVE.value
    assert matrix["physical cross-cell transfer | symbol-freezing approximation"] == DoubleChargeRelation.INDEPENDENT_COCHARGE.value
    assert matrix["initial-boundary termination | fresh interior packet"] == DoubleChargeRelation.MUTUALLY_EXCLUSIVE.value
