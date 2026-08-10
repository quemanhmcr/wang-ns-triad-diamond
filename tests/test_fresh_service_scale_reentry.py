import math

import pytest

from src.fresh_service_scale_reentry import (
    STATUS,
    canonical_annular_frame_registration,
    fresh_service_scale_route,
    pushforward_fresh_edges_to_bands,
    scale_law_statistics,
    selected_band_hard_shell_mass_lower,
    smooth_band_to_hard_shell_mass_lower,
    theorem_certificate,
)


def test_canonical_frame_registers_two_hard_shell_support():
    c = canonical_annular_frame_registration()
    assert c["support_lower_ratio"] == pytest.approx(0.5)
    assert c["support_upper_ratio"] == pytest.approx(2.0)
    assert "|phi_j|<=1" in c["pointwise_multiplier_bound"]
    assert "sum_j phi_j" in c["square_partition"]


def test_fresh_band_pushforward_quotients_cell_refinement():
    coarse = pushforward_fresh_edges_to_bands(
        [0.7, 0.3, 0.4],
        [-2, -2, -1],
        [False, False, True],
        [False, False, False],
    )
    fine = pushforward_fresh_edges_to_bands(
        [0.2, 0.5, 0.1, 0.2, 0.4],
        [-2, -2, -2, -2, -1],
        [False, False, False, False, True],
        [False, False, False, False, False],
    )
    assert coarse == {-2: pytest.approx(1.0)}
    assert fine == {-2: pytest.approx(1.0)}


def test_interface_or_old_edges_do_not_enter_fresh_scale_law():
    out = pushforward_fresh_edges_to_bands(
        [1.0, 2.0, 3.0],
        [0, -1, -2],
        [True, True, False],
        [True, False, False],
    )
    assert out == {-2: pytest.approx(3.0)}


def test_two_hard_shell_cover_has_clean_two_thirds_factor():
    assert smooth_band_to_hard_shell_mass_lower(3.0) == pytest.approx(2.0)
    assert selected_band_hard_shell_mass_lower(24.0, 2.0) == pytest.approx(2.0)


def test_scale_statistics_use_native_max_atom_and_h2_only_as_corollary():
    out = scale_law_statistics({0: 0.5, -1: 0.3, -2: 0.2})
    assert out["selected_band"] == 0
    assert out["p_max"] == pytest.approx(0.5)
    assert math.exp(-out["H_inf_scale"]) == pytest.approx(0.5)
    assert out["H2_scale"] >= out["H_inf_scale"]


def test_fresh_service_threshold_gives_one_over_24_hinf_shell_tradeoff():
    Y = 1.2
    c = 0.5
    # Exact fresh threshold Y/4, with pmax=1/2.
    law = {0: 0.15, -1: 0.09, -2: 0.06}
    out = fresh_service_scale_route(Y, c, 8.0, law, viscosity=1.0)
    clean = out["p_max"] * Y / (24.0 * c)
    base = Y / (24.0 * c)
    assert out["hard_shell_mass_lower"] >= clean - 1e-14
    assert out["H_inf_weighted_hard_shell_mass_lower"] >= base - 1e-14
    assert out["H2_weighted_hard_shell_mass_lower"] >= base - 1e-14
    assert out["next_owner"] == "generic_critical_shell_first_stop"
    assert out["maximum_candidate_shell_over_block_scale"] <= 2.0


def test_full_survivor_service_preserves_scale_conjugacy():
    out = fresh_service_scale_route(
        1.0,
        1.0,
        4.0,
        {0: 0.1, -1: 0.075, -2: 0.075},
        viscosity=1.0,
    )
    assert out["H_inf_weighted_full_survivor_service_lower"] >= out["clean_H_inf_weighted_full_survivor_service_lower"] - 1e-13
    assert out["H_inf_weighted_full_survivor_integrated_service_lower"] >= out["clean_H_inf_weighted_full_survivor_integrated_service_lower"] - 1e-13


def test_fresh_provenance_is_not_promoted_to_whole_shell_materiality():
    out = fresh_service_scale_route(1.0, 1.0, 4.0, {0: 0.25}, viscosity=1.0)
    assert "whole hard u-shell is not declared fresh material" in out["material_semantics"]
    assert "not a child-energy causal probability" in out["probability_semantics"]


def test_certificate_removes_cell_dominance_only_from_renewal_entrance():
    cert = theorem_certificate()
    assert cert["status"] == STATUS
    assert "invariant under coherent-cell refinement" in cert["material_quotient"]
    assert "not required for fresh-service renewal entrance" in cert["no_cell_dominance"]
    assert "not promoted to whole-shell freshness" in cert["material_scope"]
    assert "signed-good scale progress is not asserted" in cert["scale_scope"]
