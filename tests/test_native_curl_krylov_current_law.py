import math
import random

import pytest

from src.helical import coupling_magnitude_closed

from src.native_curl_krylov_current_law import (
    critical_boost_logistic_bound,
    critical_hilbert_square_balance,
    fixed_curl_cocycle_rhs,
    sharp_helicity_flip_boost_geometry,
    radial_fitness_selection_balance,
    parabolic_energy_clock_from_endpoints,
    self_return_operator_geometry,
    curl_spectral_bundle_base_velocity,
    canonical_spectral_triple_source,
    critical_beltrami_split,
    critical_determinant_log_rate,
    critical_determinant_live_balance,
    critical_tangent_correlation_geometry,
    critical_escape_balance,
    curl_nijenhuis_torsion_eigenfactor,
    nijenhuis_curvature_from_root_work,
    curl_krylov_state,
    hankel_vandermonde_determinant,
    heterochiral_frontier_progress_side_bound,
    intrinsic_spectral_reynolds_barrier,
    isolated_three_point_euler_law,
    martingale_observable_increment,
    modal_euler_action_decomposition,
    observable_tangent_gram,
    pairwise_curl_shear_capacity,
    rms_curl_scale_log_rate,
    spectral_source_action,
    spectral_fitness_replicator_law,
    symmetric_heterochiral_stationarity_residual,
    symmetric_heterochiral_upward_efficiency,
    theorem_certificate,
    three_point_spectral_speed_law,
    three_point_observable_volume_law,
    observable_source_speed_bound,
    three_point_current_law,
    three_point_martingale_spread,
)


def test_delta_is_exact_pairwise_signed_curl_shear_capacity():
    a = (-3.0, -0.5, 1.25, 4.0)
    e = (0.7, 1.1, 2.3, 0.4)
    st = curl_krylov_state(a, e)
    pair = sum(e[i] * e[j] * (a[i] - a[j]) ** 2 for i in range(4) for j in range(i + 1, 4))
    assert pairwise_curl_shear_capacity(a, e) == pytest.approx(pair)
    assert st.critical_determinant == pytest.approx(st.energy * st.enstrophy - st.helicity**2)
    assert st.critical_determinant == pytest.approx(pair)


def test_second_krylov_residual_is_exact_d2_over_d1():
    a = (-2.7, -0.4, 1.1, 2.2, 4.3)
    e = (0.8, 1.2, 0.5, 1.7, 0.9)
    st = curl_krylov_state(a, e)
    assert st.third_hankel_determinant == pytest.approx(
        st.critical_determinant * st.second_residual_energy
    )
    assert st.third_hankel_determinant == pytest.approx(hankel_vandermonde_determinant(a, e, 2))
    assert st.curl_defect_energy == pytest.approx(
        st.defect_energy**2 / st.energy
        + st.defect_curl_moment**2 / st.defect_energy
        + st.second_residual_energy
    )


def test_two_signed_curl_nodes_have_zero_second_krylov_residual():
    st = curl_krylov_state((-3.0, 2.0), (1.3, 0.7))
    assert st.defect_energy > 0.0
    assert st.second_residual_energy == pytest.approx(0.0, abs=1.0e-12)
    assert st.third_hankel_determinant == pytest.approx(0.0, abs=1.0e-12)
    assert st.beta2 == pytest.approx(0.0, abs=1.0e-12)


def test_log_delta_law_has_one_strain_term_and_two_native_rayleigh_terms():
    st = curl_krylov_state((-2.0, 0.75, 3.5), (0.9, 1.4, 0.6))
    out = critical_determinant_log_rate(st, defect_strain_inner=0.31, viscosity=0.08)
    direct = 2.0 * 0.31 / st.defect_energy - 0.16 * (
        st.enstrophy / st.energy + st.curl_defect_energy / st.defect_energy
    )
    assert out["log_delta_rate"] == pytest.approx(direct)
    assert out["viscous_rayleigh"] == pytest.approx(out["lanczos_rayleigh"])









def test_intrinsic_spectral_reynolds_two_is_a_necessary_delta_growth_barrier():
    st = curl_krylov_state((-2.2, -0.5, 0.9, 2.7), (0.8, 1.3, 0.7, 0.6))
    nu = 0.4
    N2 = st.enstrophy / st.energy
    # Put the spectral action strictly below the universal Re=2 threshold.
    A = 3.5 * nu * nu * N2 * st.defect_energy
    out = intrinsic_spectral_reynolds_barrier(st, A, nu)
    assert out["spectral_reynolds"] == pytest.approx(math.sqrt(3.5))
    assert out["spectral_reynolds"] < 2.0
    assert out["relaxed_log_delta_upper"] < 0.0
    assert out["fixed_state_log_delta_upper"] <= out["relaxed_log_delta_upper"] + 1.0e-12
    assert out["intrinsic_critical_mass"] < out["critical_mass_threshold"]


def test_interaction_volume_recovers_same_scale_invariant_critical_mass():
    st = curl_krylov_state((-1.4, 0.6, 2.0), (0.4, 1.0, 0.8))
    A = 0.73
    nu = 0.2
    out = intrinsic_spectral_reynolds_barrier(st, A, nu)
    N2 = st.enstrophy / st.energy
    assert out["interaction_volume"] == pytest.approx(st.critical_determinant / A)
    assert out["intrinsic_critical_mass"] == pytest.approx(
        st.energy / (N2 * out["interaction_volume"])
    )
    assert out["spectral_reynolds_squared"] == pytest.approx(
        out["intrinsic_critical_mass"] / (nu * nu)
    )


def test_isolated_three_wave_euler_speed_is_delta_circle_and_cubic_quadrature_is_invariant():
    a = (-1.9, 0.45, 2.3)
    g = 0.21 * complex(math.cos(0.37), math.sin(0.37))
    base_radii = (0.7, 1.4, 0.9)
    totals = []
    for phase in (0.0, 0.3, 0.8, 1.2, math.pi / 2):
        A = (
            math.sqrt(base_radii[0]) * complex(math.cos(0.2), math.sin(0.2)),
            math.sqrt(base_radii[1]) * complex(math.cos(-0.4), math.sin(-0.4)),
            math.sqrt(base_radii[2]) * complex(
                math.cos(phase + 0.37 + 0.2), math.sin(phase + 0.37 + 0.2)
            ),
        )
        out = isolated_three_point_euler_law(a, A, g)
        totals.append(out["total_source_speed_squared"] )
        assert out["radial_energy_speed_squared"] + out["phase_shape_speed_squared"] == pytest.approx(
            out["total_source_speed_squared"]
        )
        assert out["cubic_product_rate"].imag == pytest.approx(0.0, abs=1.0e-12)
        assert out["energy_rate_residual"] == pytest.approx(0.0, abs=1.0e-12)
        assert out["helicity_rate_residual"] == pytest.approx(0.0, abs=1.0e-12)
    assert max(totals) - min(totals) == pytest.approx(0.0, abs=1.0e-12)


def test_isolated_phase_only_rotates_fixed_nonlinear_speed_between_radial_and_phase_motion():
    a = (-1.2, 0.5, 1.8)
    g = 0.3 + 0.0j
    mags = (math.sqrt(0.8), math.sqrt(1.1), math.sqrt(0.6))
    radial = isolated_three_point_euler_law(a, mags, g)
    quadrature = isolated_three_point_euler_law(a, (mags[0], mags[1], 1j * mags[2]), g)
    assert radial["phase_cosine"] == pytest.approx(1.0)
    assert radial["phase_shape_speed_squared"] == pytest.approx(0.0, abs=1.0e-12)
    assert quadrature["phase_cosine"] == pytest.approx(0.0, abs=1.0e-12)
    assert quadrature["radial_energy_speed_squared"] == pytest.approx(0.0, abs=1.0e-12)
    assert quadrature["phase_shape_speed_squared"] == pytest.approx(
        radial["total_source_speed_squared"]
    )


def test_closed_triad_spread_is_exact_signed_curl_martingale_transport():
    a, m, b = -1.4, 0.35, 2.2
    q = 3.7
    out = three_point_martingale_spread((a, m, b), q)
    assert out["left_recipient_fraction"] + out["right_recipient_fraction"] == pytest.approx(1.0)
    assert out["conditional_mean"] == pytest.approx(m)
    assert out["conditional_variance"] == pytest.approx((m - a) * (b - m))
    assert out["quadratic_variation_rate"] == pytest.approx(q * (m - a) * (b - m))

    quadratic = martingale_observable_increment((a, m, b), q, (a * a, m * m, b * b))
    assert quadratic["second_divided_difference"] == pytest.approx(1.0)
    assert quadratic["observable_increment"] == pytest.approx(out["quadratic_variation_rate"] )


def test_critical_absolute_curl_increment_is_tanaka_reading_of_same_spread():
    a, m, b = -1.7, 0.4, 2.1
    q = 1.3
    out = martingale_observable_increment((a, m, b), q, (abs(a), abs(m), abs(b)))
    assert out["observable_increment"] > 0.0
    # If all three nodes have one sign, |a| is affine and the same generator vanishes.
    same = martingale_observable_increment((0.2, 0.7, 1.8), q, (0.2, 0.7, 1.8))
    assert same["observable_increment"] == pytest.approx(0.0, abs=1.0e-12)


def test_true_uv_log_progress_is_continuously_paid_by_same_time_side_work():
    rng = random.Random(314159)
    for _ in range(5000):
        D = rng.uniform(0.01, 0.999)
        Slo = max(1.0 - D + 1.0e-9, 0.01)
        if Slo >= 0.999:
            continue
        S = rng.uniform(Slo, 0.999)
        out = heterochiral_frontier_progress_side_bound(D, S)
        assert out["child_log_progress_per_common_current"] <= (
            out["progress_upper_from_side"] + 1.0e-12
        )
    near = heterochiral_frontier_progress_side_bound(0.999999, 0.999999)
    ratio = near["child_log_progress_per_common_current"] / near["side_work_per_common_current"]
    assert ratio == pytest.approx(2.0, rel=2.0e-6)


def test_three_point_spectral_velocity_has_exact_waleffe_delta_speed():
    a = (-2.1, 0.55, 2.7)
    e = (0.31, 2.4, 0.87)
    out = three_point_spectral_speed_law(
        a, e, waleffe_magnitude=0.37, phase_cosine=-0.61
    )
    assert out["spectral_velocity_norm_squared"] == pytest.approx(
        4.0 * 0.37**2 * 0.61**2 * out["critical_determinant"]
    )
    assert out["speed_efficiency"] == pytest.approx(2.0 * 0.37 * 0.61)
    assert out["speed_efficiency"] <= 1.0 + 1.0e-12


def test_one_spectral_velocity_controls_every_observable_speed():
    a = (-2.0, -0.4, 1.1, 2.6)
    e = (0.8, 1.2, 0.7, 1.1)
    # Build an arbitrary tangent source by projecting raw coefficients off the
    # two affine constraints in the weighted sqrt-energy representation.
    raw = [0.7, -1.1, 0.3, 0.9]
    import numpy as np
    sqrt_e = np.sqrt(np.asarray(e))
    av = np.asarray(a) * sqrt_e
    x = np.asarray(raw, dtype=float)
    q0 = sqrt_e / np.linalg.norm(sqrt_e)
    x = x - np.dot(x, q0) * q0
    av = av - np.dot(av, q0) * q0
    q1 = av / np.linalg.norm(av)
    x = x - np.dot(x, q1) * q1
    rates = tuple(2.0 * sqrt_e * x)
    act = spectral_source_action(a, e, rates)
    assert act["spectral_velocity_norm_squared"] == pytest.approx(float(np.dot(x, x)))
    for phi in (tuple(v * v for v in a), tuple(abs(v) for v in a), tuple(math.exp(0.2 * v) for v in a)):
        out = observable_source_speed_bound(a, e, rates, phi)
        assert abs(out["observable_rate"]) <= out["observable_rate_upper"] + 1.0e-11


def test_three_node_tangent_space_makes_observable_speed_cauchy_exact():
    a = (-1.7, 0.6, 2.2)
    e = (0.45, 1.8, 0.72)
    speed = three_point_spectral_speed_law(a, e, waleffe_magnitude=0.28, phase_cosine=0.77)
    rates = speed["modal_energy_rates"]
    for phi in (tuple(v * v for v in a), tuple(abs(v) for v in a), tuple(v**4 for v in a)):
        out = observable_source_speed_bound(a, e, rates, phi)
        assert abs(out["observable_rate"]) == pytest.approx(
            out["observable_rate_upper"], rel=2.0e-11, abs=2.0e-11
        )


def test_universal_observable_volume_law_contains_invariants_curvature_and_critical_stock():
    a = (-1.8, 0.65, 2.4)
    e = (0.7, 1.1, 0.9)
    kwargs = dict(waleffe_magnitude=0.29, phase_cosine_abs=0.73)

    constant = three_point_observable_volume_law(a, e, (1.0, 1.0, 1.0), **kwargs)
    affine = three_point_observable_volume_law(a, e, a, **kwargs)
    assert constant["gram_determinant"] == pytest.approx(0.0, abs=1.0e-11)
    assert affine["gram_determinant"] == pytest.approx(0.0, abs=1.0e-11)
    assert constant["observable_current_magnitude"] == pytest.approx(0.0, abs=1.0e-12)
    assert affine["observable_current_magnitude"] == pytest.approx(0.0, abs=1.0e-12)

    quadratic = three_point_observable_volume_law(a, e, tuple(x * x for x in a), **kwargs)
    st = curl_krylov_state(a, e)
    assert quadratic["gram_determinant"] == pytest.approx(st.third_hankel_determinant)

    critical = three_point_observable_volume_law(a, e, tuple(abs(x) for x in a), **kwargs)
    assert critical["gram_determinant"] > 0.0
    assert critical["observable_current_magnitude"] == pytest.approx(
        critical["represented_current_magnitude"]
    )


def test_critical_tangent_volume_vanishes_on_each_intrinsic_boundary_face():
    homochiral_a = (0.6, 1.3, 2.1)
    homochiral = observable_tangent_gram(
        homochiral_a, (0.8, 1.1, 0.7), tuple(abs(x) for x in homochiral_a)
    )
    assert homochiral["gram_determinant"] == pytest.approx(0.0, abs=1.0e-11)

    equiradial_a = (-2.0, 2.0)
    equiradial = observable_tangent_gram(
        equiradial_a, (0.9, 1.4), tuple(abs(x) for x in equiradial_a)
    )
    assert equiradial["gram_determinant"] == pytest.approx(0.0, abs=1.0e-11)



def test_critical_tangent_volume_is_radial_variance_times_sign_radius_decorrelation():
    rng = random.Random(2026081412)
    for _ in range(1000):
        n = rng.randint(3, 8)
        a = tuple((1 if rng.random() < 0.5 else -1) * 10 ** rng.uniform(-2, 2) for _ in range(n))
        e = tuple(10 ** rng.uniform(-4, 4) for _ in range(n))
        out = critical_tangent_correlation_geometry(a, e)
        E = sum(e)
        represented = (
            out["critical_determinant"]
            * out["radial_variance_component"]
            * out["decorrelation_factor"]
            / E
        )
        assert out["critical_tangent_gram"] == pytest.approx(
            represented, rel=3.0e-9, abs=3.0e-9 * max(1.0, represented)
        )

    hom = critical_tangent_correlation_geometry((0.4, 1.1, 2.8), (0.7, 1.0, 0.5))
    assert abs(hom["signed_radius_correlation"]) == pytest.approx(1.0, abs=2.0e-12)
    assert hom["critical_tangent_gram"] == pytest.approx(0.0, abs=1.0e-12)


def test_critical_beltrami_determinant_has_radial_plus_helicity_faces():
    a = (-3.0, -1.0, 0.75, 2.5)
    e = (0.4, 1.1, 0.8, 1.3)
    out = critical_beltrami_split(a, e)
    assert out["critical_determinant"] == pytest.approx(
        out["radial_variance_component"] + out["helicity_coexistence_component"]
    )
    assert out["helicity_coexistence_component"] == pytest.approx(
        4.0 * out["positive_helicity_stock"] * out["negative_helicity_stock"]
    )

    homochiral = critical_beltrami_split((1.0, 2.0, 4.0), (0.7, 1.2, 0.5))
    assert homochiral["helicity_coexistence_component"] == pytest.approx(0.0, abs=1.0e-12)

    equiradial = critical_beltrami_split((-2.0, 2.0), (0.7, 1.3))
    assert equiradial["radial_variance_component"] == pytest.approx(0.0, abs=1.0e-12)


def test_rms_curl_scale_can_only_move_up_through_euler_defect_stretching():
    st = curl_krylov_state((-2.4, -0.7, 1.2, 3.1), (0.6, 1.4, 0.9, 0.5))
    stretch = 0.37
    out = rms_curl_scale_log_rate(st, stretch, viscosity=0.11)
    theta = st.defect_energy / st.enstrophy
    sigma = stretch / st.defect_energy
    assert out["nonlinear_log_rate"] == pytest.approx(theta * sigma)
    assert out["viscous_log_rate"] <= 1.0e-14
    # Euler conversion from critical determinant growth to RMS scale growth is exact.
    delta = critical_determinant_log_rate(st, stretch, viscosity=0.0)
    assert out["nonlinear_log_rate"] == pytest.approx(
        0.5 * theta * delta["nonlinear_log_rate"]
    )


def test_three_point_barycentric_efficiency_is_sharp_and_median_gets_half_energy():
    a = (-1.6, 0.4, 2.3)
    # Ask once for the theorem's exact equality weights, then instantiate them.
    trial = three_point_current_law(a, (1.0, 1.0, 1.0), waleffe_magnitude=0.23)
    p = trial.equality_weights
    out = three_point_current_law(a, p, waleffe_magnitude=0.23)
    assert p[out.median_index] == pytest.approx(0.5)
    assert out.barycentric_efficiency == pytest.approx(1.0, abs=2.0e-12)
    assert out.gross_energy_current_magnitude == pytest.approx(out.gross_current_upper)
    assert out.gross_energy_current_magnitude <= out.global_gross_current_upper + 1.0e-12


def test_three_point_curvature_current_is_waleffe_phase_times_sqrt_d2():
    out = three_point_current_law(
        (-1.7, 0.8, 2.6),
        (0.4, 1.3, 0.9),
        waleffe_magnitude=0.31,
        phase_cosine_abs=0.6,
    )
    expected = 4.0 * 0.31 * 0.6 * math.sqrt(out.third_hankel_determinant)
    assert out.curvature_current_magnitude == pytest.approx(expected)
    assert out.curvature_current_magnitude <= 2.0 * math.sqrt(out.third_hankel_determinant) + 1.0e-12


def test_random_three_point_amplitudes_never_exceed_barycentric_efficiency_one():
    rng = random.Random(20260814)
    for _ in range(2000):
        a = sorted(rng.uniform(-4.0, 4.0) for _ in range(3))
        if min(a[1] - a[0], a[2] - a[1]) < 1.0e-5:
            continue
        e = tuple(math.exp(rng.uniform(-5.0, 5.0)) for _ in range(3))
        out = three_point_current_law(a, e, waleffe_magnitude=0.5)
        assert out.barycentric_efficiency <= 1.0 + 2.0e-12
        assert out.gross_energy_current_magnitude <= out.global_gross_current_upper + 2.0e-11


def test_symmetric_heterochiral_native_slice_recovers_the_new_interior_ratio():
    r = 0.5981296065824204
    assert symmetric_heterochiral_stationarity_residual(r) == pytest.approx(0.0, abs=3.0e-7)
    val = symmetric_heterochiral_upward_efficiency(r)
    assert val == pytest.approx(0.124802578529693, rel=2.0e-12)
    assert val > symmetric_heterochiral_upward_efficiency(0.61090410159)


def test_certificate_refuses_case_taxonomy_and_global_overclaim():
    cert = theorem_certificate()
    assert "pairwise signed-curl shear" in cert["critical_determinant"]
    assert "[C,J_u]" in cert["commutator"]
    assert "spatial concentration" in cert["global_extension_guard"]
    assert "0.5981296" in cert["symmetric_transport_note"]
    assert cert["case_taxonomy_used"] is False
    assert cert["temporal_matching_used"] is False
    assert cert["global_regularity_claimed"] is False



def test_nijenhuis_torsion_factor_is_exact_root_curvature_multiplier():
    a = (0.7, -0.9, 1.2)
    R = 1.3
    T = ((a[1] - a[2]) * R, (a[2] - a[0]) * R, (a[0] - a[1]) * R)
    Q = sum(ai * ai * ti for ai, ti in zip(a, T))
    for i in range(3):
        j, k = [x for x in range(3) if x != i]
        factor = curl_nijenhuis_torsion_eigenfactor(a[j], a[k], a[i])
        assert nijenhuis_curvature_from_root_work(a[j], a[k], a[i], T[i]) == pytest.approx(Q)
        assert factor * T[i] == pytest.approx(Q)


def test_nijenhuis_torsion_vanishes_on_parent_curl_level():
    assert curl_nijenhuis_torsion_eigenfactor(-2.0, 3.0, -2.0) == pytest.approx(0.0)
    assert curl_nijenhuis_torsion_eigenfactor(-2.0, 3.0, 3.0) == pytest.approx(0.0)


def test_critical_live_balance_is_exact_and_radial_impedance_implies_two_barrier():
    a = (-3.2, -0.8, 0.5, 2.7, 4.1)
    e = (0.6, 1.1, 0.9, 0.7, 0.4)
    # Build one nontrivial Euler source and project off the two affine null directions.
    raw = [0.8, -0.3, 0.6, -0.7, 0.2]
    n = len(a)
    S0 = sum(raw)
    S1 = sum(ai * si for ai, si in zip(a, raw))
    A1 = sum(a)
    A2 = sum(ai * ai for ai in a)
    det = n * A2 - A1 * A1
    x0 = (S0 * A2 - S1 * A1) / det
    x1 = (S1 * n - S0 * A1) / det
    rates = tuple(si - x0 - x1 * ai for ai, si in zip(a, raw))
    out = critical_escape_balance(a, e, rates, viscosity=0.17)
    assert out["critical_alignment"] <= 1.0 + 1.0e-12
    assert out["critical_alignment"] >= -1.0 - 1.0e-12
    assert out["cubic_viscous_factor"] + 1.0e-12 >= out["cubic_viscous_lower_from_radial_variance"]
    if math.isfinite(out["radial_impedance_threshold"]):
        assert out["radial_impedance_threshold"] >= 2.0 - 1.0e-12
    represented = out["live_escape_number"] - out["cubic_viscous_factor"]
    assert out["normalized_critical_rate"] == pytest.approx(represented)


def test_state_sharp_delta_reynolds_threshold_is_stronger_than_two():
    st = curl_krylov_state((-2.5, -0.4, 1.3, 3.2), (0.7, 1.0, 0.9, 0.5))
    out = intrinsic_spectral_reynolds_barrier(st, spectral_velocity_norm_squared=0.8, viscosity=0.2)
    assert out["state_sharp_reynolds_threshold"] >= 2.0
    assert out["state_sharp_reynolds_threshold"] == pytest.approx(
        2.0 * math.sqrt(1.0 + out["defect_fraction"] + out["defect_skew_fraction"])
    )



def test_uv_log_progress_is_bounded_by_same_martingale_quadratic_variation():
    rng = random.Random(2026081415)
    for _ in range(10000):
        D = rng.uniform(1.0e-7, 1.0 - 1.0e-7)
        S = rng.uniform(max(1.0e-7, 1.0 + 1.0e-7 - D), 1.0 - 1.0e-7)
        out = heterochiral_frontier_progress_side_bound(D, S)
        assert out["child_log_progress_per_common_current"] <= out["log_progress_upper_from_curvature"] + 1.0e-12
        expected = (1.0 - D) * (1.0 + S) * (D + S)
        assert out["signed_curl_curvature_per_common_current_at_unit_child_radius"] == pytest.approx(expected)


def test_modal_euler_action_has_exact_spectral_within_phase_birth_pythagorean_split():
    a = (1.0, 1.0, -2.0, -2.0, 3.0)
    z = (1.0 + 2.0j, -0.7 + 0.4j, 0.9 - 1.1j, 0.0 + 0.0j, 1.2 + 0.3j)
    f = (0.3 - 0.8j, 1.1 + 0.2j, -0.6 + 0.7j, 0.4 - 0.9j, -0.2 + 0.5j)
    out = modal_euler_action_decomposition(a, z, f)
    assert out["total_euler_action"] == pytest.approx(out["represented_total_action"])
    assert out["within_eigenspace_radial_action"] >= 0.0
    assert out["phase_rotation_action"] >= 0.0
    assert out["new_amplitude_birth_action"] == pytest.approx(abs(f[3]) ** 2)
    assert 0.0 <= out["curl_spectral_fraction"] <= 1.0



def test_spectral_fitness_is_exact_replicator_score_with_affine_euler_nulls():
    a = (-3.0, -0.7, 0.9, 2.8)
    e = (0.8, 1.1, 0.7, 0.5)
    raw = [0.6, -0.4, 0.2, -0.1]
    n = len(a); S0 = sum(raw); S1 = sum(ai * si for ai, si in zip(a, raw)); A1=sum(a); A2=sum(ai*ai for ai in a); det=n*A2-A1*A1
    x0=(S0*A2-S1*A1)/det; x1=(S1*n-S0*A1)/det
    rates=tuple(si-x0-x1*ai for ai,si in zip(a,raw))
    out=spectral_fitness_replicator_law(a,e,rates,viscosity=0.13)
    assert out["fitness_energy_mean"] == pytest.approx(0.0, abs=1.0e-12)
    assert out["fitness_helicity_mean"] == pytest.approx(0.0, abs=1.0e-12)
    assert out["fitness_action"] == pytest.approx(out["spectral_action"])
    assert out["normalized_mass_residual"] == pytest.approx(0.0, abs=1.0e-12)


def test_nonaffine_spectral_observable_has_no_phase_independent_euler_sign():
    a = (-1.7, 0.4, 2.2)
    e = (0.8, 1.1, 0.6)
    phi = tuple(abs(x) for x in a)
    plus = three_point_observable_volume_law(a,e,phi,waleffe_magnitude=0.31,phase_cosine_abs=0.73)
    # The helper returns magnitude.  Exact cyclic work is linear in cos(Phi), so a pi phase shift
    # preserves a,e,g and reverses every non-affine observable response.
    assert plus["observable_current_magnitude"] > 0.0
    det = plus["oriented_observable_determinant_abs"]
    current = 4.0 * 0.31 * 0.73 * math.sqrt(math.prod(e)) * det
    assert current > 0.0
    assert -current < 0.0



def test_krylov_impedance_is_exact_productive_coordinate_balance():
    st = curl_krylov_state((-3.1, -0.6, 1.2, 2.9, 4.0), (0.5, 1.2, 0.9, 0.7, 0.3))
    A = 0.83
    # choose a legal q2 alignment strictly inside [-1,1]
    gamma = 0.37
    c2 = gamma * math.sqrt(A)
    strain = math.sqrt(st.defect_energy) * st.beta2 * c2
    out = critical_determinant_live_balance(st, strain, A, viscosity=0.14)
    assert out["q2_alignment"] == pytest.approx(gamma)
    assert out["spectral_reconfiguration_action"] == pytest.approx(A * (1.0 - gamma * gamma))
    assert out["krylov_impedance"] + 1.0e-12 >= out["minimum_krylov_impedance"]
    direct = critical_determinant_log_rate(st, strain, 0.14)["log_delta_rate"]
    N2 = st.enstrophy / st.energy
    assert out["normalized_log_delta_rate"] == pytest.approx(direct / (2.0 * 0.14 * N2))



def test_canonical_full_pde_triple_current_is_barycentric_and_affine_null():
    a, m, b, tau = -2.3, 0.4, 3.1, 0.73
    out = canonical_spectral_triple_source(a, m, b, tau)
    assert out["energy_residual"] == pytest.approx(0.0, abs=1.0e-13)
    assert out["helicity_residual"] == pytest.approx(0.0, abs=1.0e-13)
    donor = tau * (b - a)
    assert out["median_source"] == pytest.approx(-donor)
    assert out["left_source"] == pytest.approx(donor * (b - m) / (b - a))
    assert out["right_source"] == pytest.approx(donor * (m - a) / (b - a))


def test_curl_spectral_bundle_base_law_preserves_euler_energy_and_helicity_levels():
    a = (-2.1, -0.3, 1.2, 3.4)
    raw_q = (0.5, 0.7, 0.9, 0.4)
    norm = math.sqrt(sum(x * x for x in raw_q))
    q = tuple(x / norm for x in raw_q)
    triples = (
        (0, 1, 2, 0.31),
        (0, 1, 3, -0.17),
        (0, 2, 3, 0.23),
        (1, 2, 3, -0.41),
    )
    out = curl_spectral_bundle_base_velocity(a, q, triples, viscosity=0.19)
    assert out["euler_energy_tangent_residual"] == pytest.approx(0.0, abs=1.0e-13)
    assert out["euler_helicity_tangent_residual"] == pytest.approx(0.0, abs=1.0e-13)
    assert out["normalized_mass_residual"] == pytest.approx(0.0, abs=1.0e-13)
    N2 = sum(ai * ai * qi * qi for ai, qi in zip(a, q))
    assert out["rms_curl_squared"] == pytest.approx(N2)
    assert out["viscous_velocity"] == pytest.approx(
        tuple(-0.19 * (ai * ai - N2) * qi for ai, qi in zip(a, q))
    )


def test_one_bundle_triple_matches_canonical_energy_source_exactly():
    a = (-1.8, 0.5, 2.7)
    raw_q = (0.6, 0.9, 0.7)
    norm = math.sqrt(sum(x * x for x in raw_q))
    q = tuple(x / norm for x in raw_q)
    E = 2.4
    chi = -0.37
    bundle = curl_spectral_bundle_base_velocity(a, q, ((0, 1, 2, chi),))
    rates = tuple(2.0 * E * qi * vi for qi, vi in zip(q, bundle["euler_velocity"]))
    tau = 2.0 * E * q[0] * q[1] * q[2] * chi
    src = canonical_spectral_triple_source(a[0], a[1], a[2], tau)
    assert rates == pytest.approx((src["left_source"], src["median_source"], src["right_source"]))


def test_self_return_operator_form_is_exactly_the_existing_interaction_volume_and_reynolds():
    st = curl_krylov_state((-2.4, -0.6, 1.0, 3.1), (0.8, 1.2, 0.7, 0.4))
    A = 0.91
    nu = 0.16
    old = intrinsic_spectral_reynolds_barrier(st, A, nu)
    g2 = A / st.critical_determinant
    N = math.sqrt(st.enstrophy / st.energy)
    new = self_return_operator_geometry(st.energy, N, g2, nu)
    assert new["interaction_volume"] == pytest.approx(old["interaction_volume"])
    assert new["spectral_reynolds"] == pytest.approx(old["spectral_reynolds"])
    assert new["intrinsic_critical_mass"] == pytest.approx(old["intrinsic_critical_mass"])



def test_energy_law_supplies_exact_parabolic_clock_without_event_counting():
    E0, E1, nu = 3.7, 1.4, 0.23
    clock = parabolic_energy_clock_from_endpoints(E0, E1, nu)
    assert 2.0 * nu * clock == pytest.approx(math.log(E0 / E1))


def test_critical_log_rate_is_exact_live_factor_on_energy_parabolic_clock():
    a = (-3.0, -0.7, 0.8, 2.6, 4.2)
    e = (0.7, 1.0, 0.8, 0.6, 0.3)
    raw = [0.7, -0.2, 0.5, -0.8, 0.1]
    n=len(a); S0=sum(raw); S1=sum(ai*si for ai,si in zip(a,raw)); A1=sum(a); A2=sum(ai*ai for ai in a); det=n*A2-A1*A1
    x0=(S0*A2-S1*A1)/det; x1=(S1*n-S0*A1)/det
    rates=tuple(si-x0-x1*ai for ai,si in zip(a,raw))
    nu=0.11
    out=critical_escape_balance(a,e,rates,nu)
    K=out["critical_stock"]
    direct=out["total_critical_rate"] / K
    E=sum(e); Z=sum(ai*ai*ei for ai,ei in zip(a,e)); N2=Z/E
    assert out["normalized_log_critical_rate_per_parabolic_clock"] == pytest.approx(
        direct/(2.0*nu*N2)
    )



def test_radial_geometry_optimizes_out_to_one_productive_l2_score():
    rng = random.Random(2026081416)
    for _ in range(400):
        a = tuple(sorted(rng.uniform(-4.0, 4.0) for _ in range(5)))
        if min(a[i+1]-a[i] for i in range(4)) < 1.0e-4:
            continue
        e = tuple(math.exp(rng.uniform(-2.0, 2.0)) for _ in a)
        raw = [rng.uniform(-1.0, 1.0) for _ in a]
        n=len(a); S0=sum(raw); S1=sum(ai*si for ai,si in zip(a,raw)); A1=sum(a); A2=sum(ai*ai for ai in a); det=n*A2-A1*A1
        if abs(det)<1.0e-10: continue
        x0=(S0*A2-S1*A1)/det; x1=(S1*n-S0*A1)/det
        rates=tuple(si-x0-x1*ai for ai,si in zip(a,raw))
        out=critical_escape_balance(a,e,rates,viscosity=0.17)
        lhs=out["normalized_log_critical_rate_per_parabolic_clock"]
        if math.isfinite(lhs):
            assert lhs <= out["radial_optimized_log_critical_upper"] + 2.0e-9
            assert out["energy_clock_mean_abs_curl_log_upper"] == pytest.approx(
                out["radial_optimized_log_critical_upper"] + 1.0
            )



def test_radial_fitness_selection_is_exact_forced_gradient_balance():
    rng = random.Random(2026081417)
    for _ in range(500):
        a = tuple(sorted(rng.uniform(-5.0, 5.0) for _ in range(6)))
        if min(a[i+1]-a[i] for i in range(5)) < 1.0e-4:
            continue
        e = tuple(math.exp(rng.uniform(-2.0, 2.0)) for _ in a)
        raw = [rng.uniform(-1.0, 1.0) for _ in a]
        n=len(a); S0=sum(raw); S1=sum(ai*si for ai,si in zip(a,raw)); A1=sum(a); A2=sum(ai*ai for ai in a); det=n*A2-A1*A1
        if abs(det)<1.0e-10: continue
        x0=(S0*A2-S1*A1)/det; x1=(S1*n-S0*A1)/det
        rates=tuple(si-x0-x1*ai for ai,si in zip(a,raw))
        out=radial_fitness_selection_balance(a,e,rates,viscosity=0.13)
        assert out["viscous_radial_selection"] + 1.0e-11 >= out["viscous_selection_lower"]
        assert out["normalized_log_mean_curl_rate"] <= out["moment_upper"] + 1.0e-10
        assert out["moment_upper"] <= out["quadratic_upper"] + 1.0e-10
        assert out["direct_normalized_log_mean_curl_rate"] == pytest.approx(
            out["normalized_log_mean_curl_rate"], rel=5.0e-9, abs=5.0e-9
        )


def test_homochiral_radial_fitness_is_exact_affine_null_boundary():
    a=(-4.0,-2.5,-1.1,-0.4);e=(0.7,1.0,0.5,0.9)
    # Any energy/helicity-null source also annihilates |a|=-a.
    rates=(0.3,-0.6,0.45,-0.15)
    # project exactly off 1,a for the fixture
    n=len(a);S0=sum(rates);S1=sum(ai*si for ai,si in zip(a,rates));A1=sum(a);A2=sum(ai*ai for ai in a);det=n*A2-A1*A1
    x0=(S0*A2-S1*A1)/det;x1=(S1*n-S0*A1)/det
    rates=tuple(si-x0-x1*ai for ai,si in zip(a,rates))
    out=radial_fitness_selection_balance(a,e,rates,viscosity=0.21)
    assert out["direct_normalized_log_mean_curl_rate"] == pytest.approx(
        -out["viscous_radial_selection"] + 0.0, abs=2.0e-12
    )



def test_normalized_curvature_height_is_the_same_single_productive_score():
    a=(-3.0,-0.8,0.5,2.7,4.1);e=(0.6,1.1,0.9,0.7,0.4)
    raw=[0.8,-0.3,0.6,-0.7,0.2];n=len(a);S0=sum(raw);S1=sum(ai*si for ai,si in zip(a,raw));A1=sum(a);A2=sum(ai*ai for ai in a);det=n*A2-A1*A1
    x0=(S0*A2-S1*A1)/det;x1=(S1*n-S0*A1)/det;rates=tuple(si-x0-x1*ai for ai,si in zip(a,raw))
    out=radial_fitness_selection_balance(a,e,rates,viscosity=0.17)
    assert out["normalized_curvature_height_score"] == pytest.approx(out["productive_fitness_score"])
    assert 2.0*out["curvature_height"] == pytest.approx(out["nonlinear_critical_rate"])


def test_sharp_helicity_flip_boost_has_uniform_galilean_derivative_bound():
    rng=random.Random(2026081419)
    c=3.0*math.sqrt(6.0)/16.0
    for _ in range(20000):
        q=10**rng.uniform(-5.0,3.0);l=10**rng.uniform(-5.0,3.0)
        lo=abs(l-q)+1.0e-10*max(l,q);hi=l+q-1.0e-10*max(l,q)
        if lo>=hi: continue
        k=rng.uniform(lo,hi)
        out=sharp_helicity_flip_boost_geometry(q,l,k)
        assert out["max_helicity_flip_boost"] <= c*q + 2.0e-11*max(1.0,q)
        assert out["sharp_constant"] == pytest.approx(c)


def test_sharp_helicity_flip_constant_is_approached_by_low_high_difference_ratio_half():
    q=1.0
    c=3.0*math.sqrt(6.0)/16.0
    # k-l=q/2 and k+l >> q approach the exact supremum.
    for m in (100.0,1000.0,10000.0):
        k=0.5*(m+0.5);l=0.5*(m-0.5)
        out=sharp_helicity_flip_boost_geometry(q,l,k)
        assert out["normalized_boost"] < c + 1.0e-12
    out=sharp_helicity_flip_boost_geometry(q,0.5*(10000.0-0.5),0.5*(10000.0+0.5))
    assert out["normalized_boost"] == pytest.approx(c,rel=2.0e-8)



def test_sharp_galilean_boost_bound_majorizes_actual_opposite_helicity_waleffe_matrix_element():
    rng=random.Random(2026081420)
    for _ in range(8000):
        q=10**rng.uniform(-4.0,2.0);l=10**rng.uniform(-4.0,2.0)
        lo=abs(l-q)+1.0e-9*max(l,q);hi=l+q-1.0e-9*max(l,q)
        if lo>=hi: continue
        k=rng.uniform(lo,hi)
        bound=sharp_helicity_flip_boost_geometry(q,l,k)
        s=rng.choice((-1,1));out_s=-s
        for low_s in (-1,1):
            g=coupling_magnitude_closed(q,l,k,low_s,s,out_s)
            actual=2.0*math.sqrt(k*l)*g
            assert actual <= bound["max_helicity_flip_boost"] + 2.0e-11*max(1.0,q)
            assert actual <= bound["sharp_upper"] + 2.0e-11*max(1.0,q)



def test_fixed_cartan_tensor_rhs_is_energy_helicity_null_and_heat_contracting():
    lam=(-3.1,-1.0,0.4,2.2,4.0);z=(0.7,-1.1,0.5,0.9,-0.4)
    triples=((0,1,2,0.31),(0,1,4,-0.27),(0,3,4,0.18),(1,2,3,-0.42),(2,3,4,0.23))
    nu=0.17
    out=fixed_curl_cocycle_rhs(lam,z,triples,viscosity=nu)
    assert out["euler_energy_rate"] == pytest.approx(0.0,abs=1.0e-13)
    assert out["euler_helicity_rate"] == pytest.approx(0.0,abs=1.0e-13)
    assert out["euler_phase_space_divergence"] == 0.0
    assert out["phase_space_divergence"] == pytest.approx(-nu*sum(x*x for x in lam))


def test_one_fixed_cartan_triple_reconstructs_same_canonical_spectral_energy_source():
    lam=(-2.0,0.6,3.1);z=(0.8,-0.5,1.2);f=0.37
    out=fixed_curl_cocycle_rhs(lam,z,((0,1,2,f),))
    rates=tuple(2.0*zi*fi for zi,fi in zip(z,out["euler_rhs"]))
    tau=-2.0*f*z[0]*z[1]*z[2]
    src=canonical_spectral_triple_source(lam[0],lam[1],lam[2],tau)
    assert rates == pytest.approx((src["left_source"],src["median_source"],src["right_source"]))



def test_critical_hilbert_square_completion_is_exact_and_norm_threshold_only_necessary():
    nu=0.17;M3=4.2;Q=3.7;kap=-0.8
    out=critical_hilbert_square_balance(kap,M3,Q,nu)
    assert out["critical_rate"] == pytest.approx(2*kap-2*nu*M3)
    assert out["represented_rate"] == pytest.approx(out["critical_rate"])
    assert out["critical_rate"] <= out["radial_companion_upper"]
    # A large companion norm can coexist with negative rate because orientation matters.
    assert out["companion_to_viscous_norm_ratio"] > 1.0
    assert out["critical_rate"] < 0.0


def test_critical_boost_must_beat_quadratic_mean_curl_heat_when_stock_grows():
    E=2.0;K=1.4;nu=0.1;lower=K**3/E**2
    # Pick M3 just above Jensen and kappa large enough for positive growth.
    M3=1.05*lower;kap=nu*M3+0.2
    out=critical_boost_logistic_bound(E,K,kap,M3,nu)
    assert out["critical_rate"] > 0.0
    assert out["boost_rayleigh_rate"] > out["quadratic_heat_rate"]
    assert out["critical_rate"] <= out["logistic_upper"] + 1.0e-12



def test_productive_curvature_action_is_exact_radial_fisher_regression_and_below_total_action():
    rng=random.Random(2026081421)
    for _ in range(400):
        a=tuple(sorted(rng.uniform(-5.0,5.0) for _ in range(6)))
        if min(a[i+1]-a[i] for i in range(5))<1e-4:continue
        e=tuple(math.exp(rng.uniform(-2.0,2.0)) for _ in a)
        raw=[rng.uniform(-1.0,1.0) for _ in a];n=len(a);S0=sum(raw);S1=sum(ai*si for ai,si in zip(a,raw));A1=sum(a);A2=sum(ai*ai for ai in a);det=n*A2-A1*A1
        if abs(det)<1e-10:continue
        x0=(S0*A2-S1*A1)/det;x1=(S1*n-S0*A1)/det;rates=tuple(si-x0-x1*ai for ai,si in zip(a,raw))
        nu=0.19;out=radial_fitness_selection_balance(a,e,rates,nu)
        assert out["productive_fisher_action"] <= out["total_euler_fitness_variance"] + 1e-10
        assert 0.0 <= out["productive_fisher_fraction"] <= 1.0 + 1e-10
        N2=out["rms_curl_squared"]
        assert out["physical_log_mean_curl_upper"] == pytest.approx(
            out["productive_fisher_action"]/(2.0*nu*N2)
        )
        # The exact rate must lie below the physical-time productive-action upper.
        exact_physical=out["normalized_log_mean_curl_rate"]*(2.0*nu*N2)
        assert exact_physical <= out["physical_log_mean_curl_upper"] + 2e-9
