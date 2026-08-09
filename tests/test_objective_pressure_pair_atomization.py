import math
from fractions import Fraction

import numpy as np
import pytest

from src.objective_pressure_pair_atomization import (
    STATUS,
    canonical_all_pair_absolute_capacity_upper,
    clean_dominant_pair_shell_mass_lower,
    dominant_pair_peak_shell_mass_lower,
    frobenius_dual,
    objective_pressure_pair_route,
    ordered_pair_sharp_clean_coefficient,
    pair_collision_entropy,
    scalarized_pressure_split,
    scalarized_unordered_pair_atoms,
    theorem_certificate,
    u_shell_mass_lower_from_resolved_shell_mass,
    unordered_pair_capacity_upper,
    unordered_pair_matrices,
)


def test_frobenius_dual_and_exact_pressure_split():
    HV = np.diag([0.4, -0.1, 0.2])
    HR = np.array([[0.1, 0.02, 0.0], [0.02, 0.0, 0.0], [0.0, 0.0, -0.03]])
    H = HV + HR
    Z = frobenius_dual(H)
    assert np.linalg.norm(Z, "fro") == pytest.approx(1.0)
    out = scalarized_pressure_split(H, HV, HR)
    assert out["exact_scalar_residual"] == pytest.approx(0.0, abs=1e-14)
    assert out["positive_cover_margin"] >= -1e-14


def test_unordered_pair_grouping_is_exact_and_orientation_free():
    rng = np.random.default_rng(7)
    A = rng.normal(size=(4, 4, 3, 3))
    atoms = unordered_pair_matrices(A)
    rec = sum((M for _, M in atoms), np.zeros((3, 3)))
    assert np.linalg.norm(rec - A.sum(axis=(0, 1)), "fro") < 1e-12
    atomsT = unordered_pair_matrices(np.swapaxes(A, 0, 1))
    assert [k for k, _ in atoms] == [k for k, _ in atomsT]
    for (_, M), (_, MT) in zip(atoms, atomsT):
        assert np.linalg.norm(M - MT, "fro") < 1e-12
    Z = frobenius_dual(rec)
    scalars = scalarized_unordered_pair_atoms(A, Z)
    assert sum(x for _, x in scalars) == pytest.approx(float(np.sum(Z * rec)), abs=1e-12)


def test_clean_pair_constant_comes_only_from_existing_certificates():
    c = ordered_pair_sharp_clean_coefficient()
    assert c == Fraction(256, 1425)
    assert c < Fraction(1, 5)
    assert 2 * c < Fraction(2, 5)


def test_pair_capacity_has_derivative_suppression_and_no_coherent_scale_claim():
    N = 8.0
    M = 2.0
    cap_diag = unordered_pair_capacity_upper(3.0, 5.0, M, M, N, diagonal=True)
    cap_cross = unordered_pair_capacity_upper(3.0, 5.0, M, M / 2, N, diagonal=False)
    assert cap_diag == pytest.approx((1 / 5) * (M / N) ** 4 * math.sqrt(15.0))
    assert cap_cross == pytest.approx((2 / 5) * (M / N) ** 4 * math.sqrt(15.0))


def test_countable_canonical_pair_law_has_finite_absolute_capacity():
    assert canonical_all_pair_absolute_capacity_upper(3.0, 8.0) == pytest.approx(24.0 / 2560.0)


def test_source_unit_scale_cannot_change_dimensionless_pair_dominance():
    # This is a regression against using a source-unit tolerance in q_ab comparisons.
    total = 1.0e12
    weights = [total / 5.0] * 5
    out = objective_pressure_pair_route(
        total,
        1.0,
        4.0,
        sgs_positive_source_weight=0.0,
        pair_positive_weights=weights,
        pair_shell_indices=[(i, i) for i in range(5)],
        pair_frequencies=[(1.0, 1.0)] * 5,
    )
    assert out["joint_pair_routes"] == ("diffuse_pair_source_entropy",)


def test_resolved_shell_lower_passes_to_u_only_by_contraction():
    assert u_shell_mass_lower_from_resolved_shell_mass(2.75) == pytest.approx(2.75)


def test_actual_pair_weight_exposes_peak_hard_shell_mass():
    N = 16.0
    M = 4.0
    R = 0.05
    c = 0.8
    lower = dominant_pair_peak_shell_mass_lower(R, c, M, M / 2, N, diagonal=False)
    assert lower == pytest.approx((5 / 2) * (N / M) ** 4 * R / c)


def test_clean_quarter_dominant_pair_gives_80_sigma_over_c():
    assert clean_dominant_pair_shell_mass_lower(1.0, 1.0) == pytest.approx(80.0)
    assert clean_dominant_pair_shell_mass_lower(0.7, 0.4) == pytest.approx(80.0 * 0.7 / 0.4)


def test_duplicate_pair_records_cannot_manufacture_source_entropy():
    with pytest.raises(ValueError, match="integrated once"):
        objective_pressure_pair_route(
            1.0,
            1.0,
            4.0,
            sgs_positive_source_weight=0.0,
            pair_positive_weights=[0.25, 0.25],
            pair_shell_indices=[(0, 1), (0, 1)],
            pair_frequencies=[(1.0, 0.5), (1.0, 0.5)],
        )


def test_one_shell_label_has_one_physical_frequency():
    with pytest.raises(ValueError, match="multiple observer-assigned frequencies"):
        objective_pressure_pair_route(
            1.0,
            1.0,
            4.0,
            sgs_positive_source_weight=0.0,
            pair_positive_weights=[0.25, 0.25],
            pair_shell_indices=[(0, 1), (0, 2)],
            pair_frequencies=[(1.0, 0.5), (0.75, 0.25)],
        )


def test_pressure_primary_owner_tie_is_joint():
    sigma = 1.0
    out = objective_pressure_pair_route(
        sigma,
        1.0,
        4.0,
        sgs_positive_source_weight=0.5,
        pair_positive_weights=[0.5],
        pair_shell_indices=[(0, 0)],
        pair_frequencies=[(1.0, 1.0)],
    )
    assert set(out["joint_primary_owners"]) == {"sgs_pressure_source", "resolved_pressure_pair_law"}
    assert out["sgs_stress_l32_lower_if_owner"] == pytest.approx(190.0)
    assert out["aggregate_muV_is_canonical_route"] is False


def test_exact_quarter_pair_boundary_keeps_dominance_and_entropy_joint():
    sigma = 1.0
    weights = [0.125] * 4  # pair total=sigma/2 and each normalized pair mass=1/4
    out = objective_pressure_pair_route(
        sigma,
        1.0,
        4.0,
        sgs_positive_source_weight=0.5,
        pair_positive_weights=weights,
        pair_shell_indices=[(0, 1), (2, 3), (4, 5), (6, 7)],
        pair_frequencies=[(1.0, 0.5)] * 4,
    )
    assert set(out["joint_pair_routes"]) == {
        "dominant_hard_pair_to_critical_shell",
        "diffuse_pair_source_entropy",
    }
    assert out["pair_source_entropy_lower"] == pytest.approx(math.log(4.0))
    for wit in out["dominant_pair_witnesses"]:
        assert wit["critical_shell_mass_lower"] == pytest.approx(80.0)
        assert wit["parent_to_child_natural_lifetime_ratio_at_least"] >= 16.0


def test_five_way_diffuse_pair_law_has_more_than_log4_entropy():
    ent = pair_collision_entropy([1, 1, 1, 1, 1])
    assert ent["H2_pair_source"] == pytest.approx(math.log(5.0))
    out = objective_pressure_pair_route(
        1.0,
        1.0,
        4.0,
        sgs_positive_source_weight=0.0,
        pair_positive_weights=[0.2] * 5,
        pair_shell_indices=[(i, i) for i in range(5)],
        pair_frequencies=[(1.0, 1.0)] * 5,
    )
    assert out["joint_pair_routes"] == ("diffuse_pair_source_entropy",)
    assert out["pair_source_entropy_is_causal_probability"] is False


def test_certificate_removes_aggregate_muV_from_canonical_pressure_route():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "80 Sigma_P/c" in cert["dominant_pair"]
    assert "not the canonical pressure renewal route" in cert["coarse_muV"]
    assert "not causal" in cert["positive_source_cover"]
    assert "|S_(N/4)|<=1" in cert["dominant_pair"]
