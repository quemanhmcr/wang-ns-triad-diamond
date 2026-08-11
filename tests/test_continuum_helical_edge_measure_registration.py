import math

import numpy as np
import pytest

from src.continuum_helical_edge_measure_registration import (
    CLEAN_CHANGE_OF_MEASURE,
    LOW_COST_DEFICIT_CEILING,
    continuum_edge_measure_ledger,
    edge_measure_to_service_or_flat,
    helical_coefficients,
    helical_reconstruction,
    register_continuum_triad_fiber,
    signed_good_core_physical_law,
    unitary_fourier_convolution_factor,
    unitary_sharp_young_physical_work_upper,
)
from src.helical import coupling_g, helical_basis
from src.helical_physical_edge_registration import leray_project
from src.physical_pair_weighted_productivity import physical_work_capacity_constant
from src.triad_extremizer import symmetric_gamma, symmetric_rstar


def _divfree(k: np.ndarray, raw: np.ndarray) -> np.ndarray:
    return leray_project(k, np.asarray(raw, complex))


def _pure_helical_fiber(
    x: np.ndarray,
    y: np.ndarray,
    signs: tuple[int, int, int] = (1, -1, 1),
    *,
    work_sign: float = 1.0,
    quotient_mass: float = 1.0,
):
    z = x + y
    sx, sy, sz = signs
    g = coupling_g(x, y, -z, sx, sy, sz)
    signed_frequency = sx * np.linalg.norm(x) - sy * np.linalg.norm(y)
    target = work_sign * (1.0 if signed_frequency >= 0.0 else -1.0)
    az = target * np.exp(-1j * np.angle(g))
    return register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=helical_basis(x, sx),
        uy=helical_basis(y, sy),
        uz=az * helical_basis(z, sz),
        quotient_measure_mass=quotient_mass,
    )


def _extremal_fiber(*, quotient_mass: float = 1.0):
    r = symmetric_rstar()
    gamma = symmetric_gamma(r)
    n = math.exp(-gamma)
    xx = 0.5
    yy = math.sqrt(n * n - xx * xx)
    x = np.array([xx, yy, 0.0])
    y = np.array([1.0 - xx, -yy, 0.0])
    assert math.isclose(np.linalg.norm(x), r, rel_tol=2e-12, abs_tol=2e-12)
    assert math.isclose(np.linalg.norm(y), r, rel_tol=2e-12, abs_tol=2e-12)
    return _pure_helical_fiber(x, y, quotient_mass=quotient_mass)


def test_unitary_fourier_factor_is_native_and_clean_productivity_upper_deliberately_dominates():
    cf = unitary_fourier_convolution_factor()
    assert math.isclose(cf, (2.0 * math.pi) ** (-1.5), rel_tol=0.0, abs_tol=0.0)
    assert 0.0 < cf < 1.0
    assert unitary_sharp_young_physical_work_upper(1.0) < physical_work_capacity_constant(1.0)


def test_arbitrary_divergence_free_vector_fiber_reconstructs_all_eight_helicity_channels():
    x = np.array([1.2, -0.3, 0.4])
    y = np.array([-0.2, 1.1, 0.7])
    z = x + y
    ux = _divfree(x, np.array([1.0 + 0.2j, -0.4 + 0.7j, 0.3 - 0.5j]))
    uy = _divfree(y, np.array([-0.3 + 0.9j, 0.8 - 0.1j, 0.2 + 0.4j]))
    uz = _divfree(z, np.array([0.5 - 0.6j, -0.7 + 0.2j, 0.9 + 0.3j]))

    coeff = helical_coefficients(x, ux)
    rebuilt, residual = helical_reconstruction(x, ux)
    assert set(coeff) == {-1, 1}
    assert residual < 3e-10
    assert np.linalg.norm(rebuilt - ux) < 3e-10

    fiber = register_continuum_triad_fiber(
        x=x,
        y=y,
        z=z,
        ux=ux,
        uy=uy,
        uz=uz,
        quotient_measure_mass=0.7,
    )
    assert len(fiber.modal_atoms) == 8
    assert fiber.ordered_quotient_source_residual < 3e-11
    assert fiber.parent_swap_residual < 3e-11
    assert abs(fiber.signed_work_reconstruction_residual) < 5e-10
    assert abs(fiber.signed_progress_reconstruction_residual) < 5e-10

    ledger = continuum_edge_measure_ledger((fiber,))
    assert abs(ledger.signed_direct_work - ledger.signed_modal_work) < 1e-9
    assert abs(ledger.signed_direct_progress - ledger.signed_registered_progress) < 1e-9
    assert abs(ledger.polarization_residual) < 1e-12
    assert ledger.capacity_is_causal_law is False
    assert ledger.parent_orientation_chosen is False


def test_parent_swap_is_a_quotient_not_a_second_physical_edge():
    x = np.array([0.8, 0.3, 0.1])
    y = np.array([0.1, 0.7, -0.2])
    z = x + y
    ux = _divfree(x, np.array([0.4 + 0.3j, -0.2 + 0.8j, 0.7 - 0.1j]))
    uy = _divfree(y, np.array([-0.6 + 0.2j, 0.5 - 0.4j, 0.3 + 0.9j]))
    uz = _divfree(z, np.array([0.2 - 0.1j, 0.7 + 0.5j, -0.4 + 0.3j]))
    a = register_continuum_triad_fiber(
        x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, quotient_measure_mass=1.3
    )
    b = register_continuum_triad_fiber(
        x=y, y=x, z=z, ux=uy, uy=ux, uz=uz, quotient_measure_mass=1.3
    )
    la = continuum_edge_measure_ledger((a,))
    lb = continuum_edge_measure_ledger((b,))
    assert math.isclose(la.signed_direct_work, lb.signed_direct_work, rel_tol=2e-10, abs_tol=2e-10)
    assert math.isclose(la.capacity_mass, lb.capacity_mass, rel_tol=2e-10, abs_tol=2e-10)
    assert math.isclose(la.signed_registered_progress, lb.signed_registered_progress, rel_tol=2e-10, abs_tol=2e-10)
    assert math.isclose(la.positive_edge_work, lb.positive_edge_work, rel_tol=2e-10, abs_tol=2e-10)


def test_signed_edge_measure_is_reconstructed_before_hahn_and_positive_atoms_only_dominate_aggregate():
    r = symmetric_rstar()
    theta = math.acos(1.0 / (2.0 * r * r) - 1.0)
    x = r * np.array([math.cos(theta / 2), math.sin(theta / 2), 0.0])
    y = r * np.array([math.cos(theta / 2), -math.sin(theta / 2), 0.0])
    positive = _pure_helical_fiber(x, y, work_sign=1.0, quotient_mass=1.0)
    negative = _pure_helical_fiber(x, y, work_sign=-1.0, quotient_mass=0.8)
    ledger = continuum_edge_measure_ledger((positive, negative))
    assert math.isclose(
        ledger.positive_edge_work - ledger.negative_edge_work,
        ledger.signed_direct_work,
        rel_tol=5e-10,
        abs_tol=5e-10,
    )
    assert ledger.positive_edge_work >= ledger.fiber_positive_work - 1e-10
    assert ledger.fiber_positive_work >= ledger.aggregate_positive_work - 1e-10
    assert ledger.positive_dominance_over_aggregate >= -1e-10
    assert ledger.positive_dominance_over_fibers >= -1e-10


def test_positive_nonforward_work_remains_physical_and_routes_to_native_transfer_cost():
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([-0.8, 0.6, 0.0])
    fiber = _pure_helical_fiber(x, y)
    ledger = continuum_edge_measure_ledger((fiber,))
    assert np.linalg.norm(x + y) < max(np.linalg.norm(x), np.linalg.norm(y))
    assert ledger.positive_nonforward_work > 0.0
    assert ledger.positive_forward_work == pytest.approx(0.0, abs=1e-12)
    assert ledger.signed_registered_progress == pytest.approx(0.0, abs=1e-12)
    assert ledger.block_transfer_deficit == pytest.approx(1.0, abs=2e-10)

    out = edge_measure_to_service_or_flat(
        ledger,
        tau=0.01,
        objective_variation_action=0.0,
        total_strain_action=0.0,
        coherent_deformation_action=0.0,
        aspect=1.0,
        scale_radius=1.0,
        has_predecessor=True,
        scaled_lifetime=1.0,
    )
    causes = tuple(out["triggered_causes"])
    assert any(c["cause"] == "physical_transfer_cost" for c in causes)
    assert out["transfer_deficit_source"] == "actual_continuum_signed_A_J_c_edge_measure"
    assert out["capacity_used_as_causal_law"] is False


def test_low_deficit_extremal_physical_edge_supplies_actual_good_core_change_of_measure():
    ledger = continuum_edge_measure_ledger((_extremal_fiber(),))
    assert 0.0 <= ledger.block_transfer_deficit < LOW_COST_DEFICIT_CEILING
    core = signed_good_core_physical_law(ledger)
    assert core.certified_capacity_fraction_lower > 0.5
    assert core.realized_capacity_fraction >= core.certified_capacity_fraction_lower - 5e-10
    assert core.density_condition_number < CLEAN_CHANGE_OF_MEASURE
    assert core.clean_normalized_rn_lower == pytest.approx(50.0 / 53.0)
    assert core.clean_normalized_rn_upper == pytest.approx(53.0 / 50.0)
    assert core.realized_normalized_rn_min >= 50.0 / 53.0 - 5e-9
    assert core.realized_normalized_rn_max <= 53.0 / 50.0 + 5e-9
    assert core.physical_hodge_defect_mean_upper >= 0.0


def test_continuum_layer_rejects_nonphysical_measure_and_nondivergencefree_inputs():
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.2, 1.0, 0.0])
    z = x + y
    ux = helical_basis(x, 1)
    uy = helical_basis(y, -1)
    uz = helical_basis(z, 1)
    with pytest.raises(ValueError, match="nonnegative finite quotient-measure"):
        register_continuum_triad_fiber(
            x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, quotient_measure_mass=-1.0
        )
    with pytest.raises(ValueError, match="divergence free"):
        register_continuum_triad_fiber(
            x=x,
            y=y,
            z=z,
            ux=np.array([1.0, 0.0, 0.0], complex),
            uy=uy,
            uz=uz,
            quotient_measure_mass=1.0,
        )
