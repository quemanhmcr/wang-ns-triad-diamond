import math
import random

import pytest

from src.native_curl_krylov_current_law import (
    critical_beltrami_split,
    critical_determinant_log_rate,
    critical_tangent_correlation_geometry,
    curl_krylov_state,
    hankel_vandermonde_determinant,
    heterochiral_frontier_progress_side_bound,
    intrinsic_spectral_reynolds_barrier,
    isolated_three_point_euler_law,
    martingale_observable_increment,
    observable_tangent_gram,
    pairwise_curl_shear_capacity,
    rms_curl_scale_log_rate,
    spectral_source_action,
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
