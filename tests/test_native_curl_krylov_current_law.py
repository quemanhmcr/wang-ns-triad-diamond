import math
import random

import numpy as np
import pytest

from src.helical import coupling_g, coupling_magnitude_closed, helical_basis

from src.native_curl_krylov_current_law import (
    continuum_midpoint_operator_sobolev_dictionary,
    critical_energy_sphere_escape_geometry,
    spacetime_critical_hom_coefficients,
    endogenous_euler_pair_projection,
    continuum_primitive_critical_channel_constants,
    continuum_critical_carre_du_champ_constants,
    continuum_critical_gauss_bianchi_constants,
    continuum_critical_operator_isometry_constant,
    closed_triad_critical_action_bound,
    sobolev_spectral_hilbert_square,
    critical_spectral_hminus_half_square,
    poisson_critical_scale_measure_moments,
    poisson_energy_boundary_jet_geometry,
    poisson_scalar_reversible_counterprofile,
    critical_loggap_collective_bound,
    critical_logscale_strain_kernel,
    sobolev_strain_transfer_multiplier,
    radial_mean_resolvent_balance,
    normalized_triple_orientation_heat_law,
    cartan_hminus1_gradient_split,
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
    assert "nu_E=kappa/M3" in cert["primitive_regeneration_persistence"]
    assert "2/5<=alpha<1/2" in cert["primitive_affine_core_budget_falsifier"]
    assert "||A_u||_pair^2~tau^(-4+6alpha)" in cert["primitive_affine_core_budget_falsifier"]
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

    ratios=[]
    for eps in (1e-2,1e-3,1e-4):
        out=heterochiral_frontier_progress_side_bound(1-eps*eps,eps)
        assert out["high_child_retained_fraction"] > 10.0/13.0
        assert out["child_log_progress_per_common_current"] < 1.02*eps*eps
        ratios.append(out["child_log_progress_per_common_current"]/out["log_progress_upper_from_curvature"])
    assert ratios[-1] > 0.9998 and ratios == sorted(ratios)


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



def test_dirichlet_profile_does_not_quotient_euler_orientation_or_pay_escape_action():
    import numpy as np

    # One conditioned physical closed-triad witness.  The common Waleffe scalar is normalized
    # to one because both compared actions are homogeneous of degree two in that scalar.
    a=np.array((-0.12348794255990123,-0.1198258560874539,0.12269735720091647))
    e=np.array((0.10384664767284751,0.04215218685969257,0.03798502944543756))
    r=np.abs(a)
    assert r[0]+r[1]>r[2] and r[1]+r[2]>r[0] and r[2]+r[0]>r[1]
    signs=tuple(1 if x>0 else -1 for x in a)
    assert coupling_magnitude_closed(*r,*signs) > 0.15
    T=np.array((a[1]-a[2],a[2]-a[0],a[0]-a[1]))
    assert float(np.sum(T)) == pytest.approx(0.0,abs=2e-16)
    assert float(np.sum(a*T)) == pytest.approx(0.0,abs=2e-16)

    # A pi common-phase reversal fixes every radial q datum but sends j_E,kappa to their negatives.
    for y in (0.0,0.7/r.min(),3.0/r.min()):
        decay=np.exp(-2.0*r*y)
        q=2.0*float(np.sum(r*r*e*decay))
        jp=float(np.sum(r*T*decay)); jm=float(np.sum(r*(-T)*decay))
        assert q > 0.0
        assert jm == pytest.approx(-jp,rel=2e-15,abs=2e-15)
    kappa=0.5*float(np.sum(r*T))
    assert kappa != 0.0
    assert -kappa == pytest.approx(0.5*float(np.sum(r*(-T))),rel=2e-15,abs=2e-15)

    # The tempting lower comparison A_escape <= int j_E^2/q is false with a wide margin.
    E=float(np.sum(e)); K=float(np.sum(r*e)); Z=float(np.sum(r*r*e)); det=E*Z-K*K
    relative_defect=det/(E*Z)
    Aescape=kappa*kappa/((Z/E)*det)
    assert relative_defect == pytest.approx(1.4642530609247162e-4,rel=2e-12)
    nodes,weights=np.polynomial.legendre.leggauss(256); r0=float(r.min())
    def fx(x):
        decay=np.exp(-2.0*r*x/r0)
        q=2.0*float(np.sum(r*r*e*decay)); j=float(np.sum(r*T*decay))
        return 0.0 if q==0.0 else (j*j/q)/r0
    transport=0.0
    for lo,hi in ((0,1),(1,4),(4,12),(12,30),(30,80),(80,160),(160,320),(320,400)):
        xx=0.5*(hi-lo)*nodes+0.5*(hi+lo)
        transport += 0.5*(hi-lo)*float(weights@np.array([fx(float(x)) for x in xx]))
    ratio=transport/Aescape
    assert transport == pytest.approx(6.13902345443830e-4,rel=2e-12,abs=2e-15)
    assert ratio == pytest.approx(3.39325211604616e-6,rel=3e-12)
    assert ratio < 1.0e-5


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


def test_euler_acceleration_normals_and_critical_leaf_turning_decomposition_are_exact():
    lam=np.array((-3.1,-1.0,0.4,2.2,4.0));z=np.array((0.7,-1.1,0.5,0.9,-0.4))
    triples=((0,1,2,0.31),(0,1,4,-0.27),(0,3,4,0.18),(1,2,3,-0.42),(2,3,4,0.23))
    def Q(x): return np.array(fixed_curl_cocycle_rhs(lam,x,triples)["euler_rhs"])
    F=Q(z);Fp=Q(z+F)-Q(z)-Q(F)  # exact polarization of the quadratic Euler field
    assert float(z@Fp) == pytest.approx(-float(F@F),abs=3e-14)
    assert float((lam*z)@Fp) == pytest.approx(-float((lam*F)@F),abs=3e-14)

    abs_lam=np.abs(lam);E=float(z@z);H=float(z@(lam*z));Z=float((lam*z)@(lam*z))
    K=float(z@(abs_lam*z));X=float((lam*z)@(abs_lam*z));G=np.array(((E,H),(H,Z)))
    a,b=np.linalg.solve(G,np.array((K,X)));g=abs_lam*z-a*z-b*lam*z
    assert float(g@z) == pytest.approx(0.0,abs=2e-14)
    assert float(g@(lam*z)) == pytest.approx(0.0,abs=2e-14)
    c,d=np.linalg.solve(G,np.array((float(z@Fp),float((lam*z)@Fp))))
    PT=Fp-c*z-d*lam*z
    direct=float(F@(abs_lam*F)+(abs_lam*z)@Fp)
    represented=float(F@((abs_lam-a-b*lam)*F)+g@PT)
    assert direct == pytest.approx(represented,abs=3e-14)


def test_normalized_projector_lax_krein_radial_flatness_and_covariant_orientation_are_exact():
    def J(x):
        x=np.asarray(x); return np.array(((0.,-x[2],x[1]),(x[2],0.,-x[0]),(-x[1],x[0],0.)))

    C=np.diag((0.8,-1.7,3.1));L=np.diag((0.8,1.7,3.1));L2=L@L;nu=0.19
    u=np.array((0.7,-1.1,0.6));E=float(u@u);e=u/math.sqrt(E);P=np.outer(e,e)
    F=J(u)@(C@u);f=F/math.sqrt(E);N2=float((C@u)@(C@u))/E
    v=f-nu*(L2@e-N2*e);Pt=np.outer(v,e)+np.outer(e,v)
    Ae=math.sqrt(E)*(J(e)@C+C@J(e))
    assert np.linalg.norm(Ae+Ae.T)<2e-14 and np.linalg.norm(Ae@e-f)<2e-14
    Aeminus=math.sqrt(E)*(J(-e)@C+C@J(-e))
    assert np.linalg.norm(Aeminus+Ae)<2e-14  # same P, opposite Euler orientation generator
    comm=Ae@P-P@Ae;PL=P@L2-L2@P;double=P@PL-PL@P
    assert np.linalg.norm(Pt-comm+nu*double)<3e-14
    heat=-2*nu*float(((L2-N2*np.eye(3))@e)@((L2-N2*np.eye(3))@e))
    assert 2*float((-nu*(L2-N2*np.eye(3))@e)@(L2@e)) == pytest.approx(heat,abs=3e-14)

    # Euler co-rotating frame: at U=I the connection cancels Euler exactly and only heat remains.
    du=Ae@u-nu*L2@u
    assert np.linalg.norm(-Ae@u+du+nu*L2@u)<3e-14
    Bframe=L@Ae-Ae@L
    assert np.linalg.norm((-Ae@L+L@Ae)-Bframe)<3e-14

    # Jacobi/Lax plus operator heat intertwining in the finite cross-product Lie algebra.
    omega=np.array((0.4,-0.8,1.2));adv=np.array((-0.5,0.3,0.9))
    assert np.linalg.norm(J(np.cross(omega,adv))-(J(omega)@J(adv)-J(adv)@J(omega)))<3e-14
    basis=np.eye(3);Ds=[-1j*J(q) for q in basis];B=J(omega).astype(complex)
    dop=sum((D@(D@B-B@D)-(D@B-B@D)@D) for D in Ds)
    assert np.linalg.norm(dop-2*B)<4e-14

    # B-invariance implies the positive-metric Krein adjoint formula.
    Js=np.diag((1.,1.,-1.,-1.));G=np.diag((0.7,1.3,2.1,3.4));Bmet=G@Js
    K=np.array(((0.,1.,-2.,.5),(-1.,0.,.3,.8),(2.,-.3,0.,-1.1),(-.5,-.8,1.1,0.)))
    U=np.linalg.solve(Bmet,K);Udag=np.linalg.solve(G,U.T@G)
    assert np.linalg.norm(U.T@Bmet+Bmet@U)<2e-14
    assert np.linalg.norm(Udag+Js@U@Js)<3e-14

    # Canonical radial exponential-tilt generators commute as vector fields.
    r=np.array((0.6,1.1,2.3,4.0));p0=np.array((.1,.2,.3,.4))
    def V(g): return 2*(g-float(p0@g))*p0
    def DV(g,q): return 2*((g-float(p0@g))*q-p0*float(g@q))
    Vr,Vq=V(r),V(r*r)
    assert np.linalg.norm(DV(r*r,Vr)-DV(r,Vq))<2e-14

    # Full normalized NS uses the same absolute-curl Krylov frame and closes D_t w exactly.
    m=float(e@(L@e));delta=math.sqrt(float(e@(L2@e))-m*m);n=(L@e-m*e)/delta
    a=float(n@f);w=f-a*n;alpha=float(n@(L@n));h=L@n-delta*e-alpha*n

    # One self-adjoint frame commutator carries scale growth, spread and turning.
    mE=float(e@Bframe@e);dE=float(n@Bframe@e)
    nE=Ae@n+(Bframe@e-mE*e-dE*n)/delta
    assert mE == pytest.approx(2*a*delta,abs=3e-14)
    assert np.linalg.norm(Bframe@e-mE*e-dE*n-delta*(nE-Ae@n))<4e-14
    heat_m=-2*nu*(float(e@(L@L@L@e))-m*float(e@(L2@e)))
    assert heat_m <= 2e-14
    A=a-nu*delta*(m+alpha);W=w-nu*delta*h;Pp=np.eye(3)-np.outer(e,e)-np.outer(n,n)
    assert np.linalg.norm((L2-N2*np.eye(3))@e-delta*(m+alpha)*n-delta*h)<3e-14
    assert np.linalg.norm(v-A*n-W)<3e-14
    md=2*float((L@e)@v);dd=(float((L2@e)@v)-m*md)/delta
    nd=((L-m*np.eye(3))@v-md*e)/delta-(dd/delta)*n
    nperp=(A*h+Pp@((L-m*np.eye(3))@W))/delta
    assert np.linalg.norm(Pp@nd-nperp)<5e-14
    ft=-nu*N2*f+math.sqrt(E)*(J(v)@(C@e)+J(e)@(C@v))
    Dtw=Pp@ft-a*nperp
    wrhs=math.sqrt(E)*Pp@(J(v)@(C@e)+J(e)@(C@v))-nu*N2*w-a*nperp
    assert np.linalg.norm(Dtw-wrhs)<5e-14


def test_absolute_curl_productive_frame_turning_helicity_toda_and_two_radius_rigidity():
    r=np.array((0.7,1.1,1.9,2.4,3.3,4.2));sgn=np.array((1.,-1.,1.,-1.,1.,-1.));C=sgn*r
    e=np.array((0.31,-0.27,0.44,0.18,-0.53,0.39));e=e/np.linalg.norm(e)
    raw=np.array((0.8,-0.2,0.6,-0.7,0.15,0.4));A=np.stack((e,C*e),axis=1)
    f=raw-A@np.linalg.solve(A.T@A,A.T@raw)  # Euler-like E/H tangent velocity
    m=float(e@(r*e));delta=math.sqrt(float(e@(r*r*e))-m*m);n=(r-m)*e/delta
    a=float(n@f);w=f-a*n;alpha=float(n@(r*n));Ln=r*n
    h=Ln-e*float(e@Ln)-n*float(n@Ln)
    md=2*float((r*e)@f);zd=2*float((r*r*e)@f);dd=(zd-2*m*md)/(2*delta)
    P=lambda x:x-e*float(e@x)-n*float(n@x)
    nd=((r-m)*f-md*e)/delta-(dd/delta)*n
    represented=-a*e+(a*h+P((r-m)*w))/delta
    assert md == pytest.approx(2*a*delta,abs=3e-14)
    assert dd == pytest.approx(a*(alpha-m)+float(h@w),abs=3e-14)
    assert np.linalg.norm(nd-represented) < 8e-14

    cH=float((C*e)@n);vH=P(C*e)
    assert a*cH+float(vH@w) == pytest.approx(0.0,abs=5e-14)
    assert float(w@w)+2e-14 >= a*a*cH*cH/float(vH@vH)
    theta=P(nd+a*e);Rw=P((r-m)*w)
    assert np.linalg.norm(delta*theta-a*h-Rw) < 8e-14
    assert a*a*float(h@h) <= 2*delta*delta*float(theta@theta)+2*float(Rw@Rw)+2e-13

    # Perfect radial productivity is the first Toda spectral-measure vector field after ds/dt=a/delta.
    rr=np.array((0.8,1.7,3.6,4.4));p0=np.array((0.17,0.29,0.31,0.23));mt=float(p0@rr)
    dt=math.sqrt(float(p0@((rr-mt)**2)));at=0.73
    pdot=2*(at/dt)*(rr-mt)*p0;alphat=float((p0*((rr-mt)**2))@rr)/(dt*dt)
    assert float(pdot.sum()) == pytest.approx(0.0,abs=2e-14)
    assert float(pdot@rr) == pytest.approx(2*at*dt,abs=3e-14)
    vard=float(pdot@((rr-mt)**2)); assert vard/(2*dt) == pytest.approx(at*(alphat-mt),abs=3e-14)
    assert np.max(np.abs((dt/at)*pdot-2*(rr-mt)*p0)) < 2e-14

    # Two occupied absolute-curl radii make span{e,n} Lambda-invariant; a third generic radius opens h.
    rr=np.array((1.,1.,3.,3.));ee=np.array((.3,.4,.5,math.sqrt(.5)));mt=float(ee@(rr*ee));dt=math.sqrt(float(ee@(rr*rr*ee))-mt*mt);nn=(rr-mt)*ee/dt
    hh=rr*nn-ee*float(ee@(rr*nn))-nn*float(nn@(rr*nn)); assert np.linalg.norm(hh)<5e-14
    rr=np.array((1.,2.,4.));ee=np.sqrt(np.array((.2,.3,.5)));mt=float(ee@(rr*ee));dt=math.sqrt(float(ee@(rr*rr*ee))-mt*mt);nn=(rr-mt)*ee/dt
    hh=rr*nn-ee*float(ee@(rr*nn))-nn*float(nn@(rr*nn)); assert np.linalg.norm(hh)>0.6


def test_one_fixed_cartan_triple_reconstructs_same_canonical_spectral_energy_source():
    lam=(-2.0,0.6,3.1);z=(0.8,-0.5,1.2);f=0.37
    out=fixed_curl_cocycle_rhs(lam,z,((0,1,2,f),))
    rates=tuple(2.0*zi*fi for zi,fi in zip(z,out["euler_rhs"]))
    tau=-2.0*f*z[0]*z[1]*z[2]
    src=canonical_spectral_triple_source(lam[0],lam[1],lam[2],tau)
    assert rates == pytest.approx((src["left_source"],src["median_source"],src["right_source"]))


def test_one_fixed_physical_triad_has_no_one_step_productive_turning_sign_law():
    k=np.array((1.0,0.0,0.0));p=np.array((0.0,1.0,0.0));q=-(k+p);signs=(1,-1,1)
    g=coupling_g(k,p,q,*signs);r=tuple(float(np.linalg.norm(x)) for x in (k,p,q));a=tuple(s*ri for s,ri in zip(signs,r))
    witnesses=(
        ((-0.8279616875058876+0.3946891248953539j),(1.5792101934448373+0.9842738001575556j),(-0.8004144025010349+0.6368837919875677j)),
        ((-0.8443800889130898+1.3788468817563762j),(1.2550724935537503+0.9000320968569572j),(-1.383873821313951+0.7176558245071815j)),
        ((-0.5271256718618944-0.28732080739415755j),(-0.8791133630368287-2.193513141320069j),(0.23100288463833551-0.5310286435217219j)),
        ((0.11754085029501977-0.6626806610393868j),(-0.44037096703025486+0.2573129155878022j),(0.17955721141905312-0.30347351520429744j)),
    )
    quadrants=set()
    for A in witnesses:
        out=isolated_three_point_euler_law(a,A,g);F=out["rhs"];d=(a[1]-a[2],a[2]-a[0],a[0]-a[1])
        Fdot=(2*d[0]*g*(F[1]*A[2]+A[1]*F[2]).conjugate(),2*d[1]*g*(F[2]*A[0]+A[2]*F[0]).conjugate(),2*d[2]*g*(F[0]*A[1]+A[0]*F[1]).conjugate())
        kap=sum(ri*(Ai.conjugate()*Fi).real for ri,Ai,Fi in zip(r,A,F))
        kdot=sum(ri*(abs(Fi)**2+(Ai.conjugate()*Fdi).real) for ri,Ai,Fi,Fdi in zip(r,A,F,Fdot))
        h=2e-7
        def kap_at(B):
            G=isolated_three_point_euler_law(a,B,g)["rhs"]
            return sum(ri*(Bi.conjugate()*Gi).real for ri,Bi,Gi in zip(r,B,G))
        fd=(kap_at(tuple(Ai+h*Fi for Ai,Fi in zip(A,F)))-kap_at(tuple(Ai-h*Fi for Ai,Fi in zip(A,F))))/(2*h)
        assert kdot == pytest.approx(fd,rel=2e-7,abs=2e-8)
        quadrants.add((1 if kap>0 else -1,1 if kdot>0 else -1))
    assert quadrants == {(-1,-1),(-1,1),(1,-1),(1,1)}

    A=((-0.386654644603378-0.40492849999157815j),(-0.10954380694909455-0.33416886975520926j),(-0.11456682205321557-0.4793283181315438j))
    out=isolated_three_point_euler_law(a,A,g);F=out["rhs"];d=(a[1]-a[2],a[2]-a[0],a[0]-a[1])
    Fdot=(2*d[0]*g*(F[1]*A[2]+A[1]*F[2]).conjugate(),2*d[1]*g*(F[2]*A[0]+A[2]*F[0]).conjugate(),2*d[2]*g*(F[0]*A[1]+A[0]*F[1]).conjugate())
    E=sum(abs(x)**2 for x in A);K=sum(ri*abs(x)**2 for ri,x in zip(r,A));Z=sum(ri*ri*abs(x)**2 for ri,x in zip(r,A));D=E*Z-K*K
    kap=sum(ri*(x.conjugate()*f).real for ri,x,f in zip(r,A,F));kdot=sum(ri*(abs(f)**2+(x.conjugate()*fd).real) for ri,x,f,fd in zip(r,A,F,Fdot))
    Zdot=2*sum(ri*ri*(x.conjugate()*f).real for ri,x,f in zip(r,A,F));Ddot=E*Zdot-4*K*kap
    eta=kap*math.sqrt(E/(Z*D));etad=math.sqrt(E/(Z*D))*(kdot-0.5*kap*(Zdot/Z+Ddot/D))
    assert eta>0.04 and etad>0.02


def test_weighted_jacobi_is_bracket_level_and_full_euler_companions_can_vanish():
    def mode(k,s): return (np.asarray(k,float),int(s))
    def lam(m): return m[1]*float(np.linalg.norm(m[0]))
    def ccoef(out,a,b):
        ko,so=out;ka,sa=a;kb,sb=b
        if np.linalg.norm(ko-ka-kb)>1e-12:return 0j
        ho,ha,hb=helical_basis(ko,so),helical_basis(ka,sa),helical_basis(kb,sb)
        return np.vdot(ho,1j*(np.dot(ha,kb)*hb-np.dot(hb,ka)*ha))
    def fcoef(out,a,b): return ccoef(out,a,b)/lam(out)
    def weighted(final,a,b,c):
        km=a[0]+b[0]
        return sum(lam(M:=mode(km,sm))*fcoef(M,a,b)*fcoef(final,c,M) for sm in (-1,1))
    def full(final,a,b,c):
        km=a[0]+b[0]
        return sum((lam(b)-lam(a))*(lam(M:=mode(km,sm))-lam(c))*fcoef(M,a,b)*fcoef(final,c,M) for sm in (-1,1))

    I,J,K,L=mode((1,0,0),1),mode((0,1,0),-1),mode((0,0,1),1),mode((1,1,1),-1)
    A,B,C=weighted(L,J,K,I),weighted(L,K,I,J),weighted(L,I,J,K)
    scale=max(abs(A),abs(B),abs(C))
    assert abs(A+B+C)<=8e-15*scale
    assert abs(A)**2<=2*(abs(B)**2+abs(C)**2)+1e-14
    assert max(abs(B),abs(C))+1e-14>=abs(A)/2
    outer=(0.7-0.2j)*(-0.4+0.6j)*(1.1+0.3j)*(0.5-0.8j)
    assert abs(outer*(A+B+C))<=8e-15*abs(outer)*scale

    # Restoring the actual Euler signed-curl gaps destroys any universal companion lower bound.
    I,J,K,L=mode((2,2,-1),1),mode((0,0,3),1),mode((-2,-2,2),-1),mode((0,0,4),-1)
    Aw,Bw,Cw=weighted(L,J,K,I),weighted(L,K,I,J),weighted(L,I,J,K)
    Af,Bf,Cf=full(L,J,K,I),full(L,K,I,J),full(L,I,J,K)
    assert abs(Aw+Bw+Cw)<=3e-15*max(abs(Aw),abs(Bw),abs(Cw))
    assert abs(Af)>5.0
    assert abs(Bf)<2e-14 and abs(Cf)<2e-14
    assert abs(Cw)>0.25  # the bracket companion exists; a curl gap kills its full Euler path

    m1=J[0]+K[0];m2=K[0]+I[0];m3=I[0]+J[0];ell=L[0]
    assert np.linalg.norm(m1+m2+m3-2*ell)<1e-14
    assert max(np.linalg.norm(m1),np.linalg.norm(m2),np.linalg.norm(m3))+1e-14 >= 2*np.linalg.norm(ell)/3
    assert np.linalg.norm(m3)>=2*np.linalg.norm(ell)/3 and abs(Cf)<2e-14



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


def test_hminus1_cartan_motion_is_orthogonal_to_energy_gradient_heat():
    a=(-3.0,-1.2,0.7,2.4)
    z=(0.8+0.2j,-0.6+0.5j,0.9-0.3j,-0.4-0.7j)
    raw=(0.3-0.4j,-0.2+0.7j,0.5+0.1j,-0.8+0.2j)
    # Remove the single real energy-radial component from the arbitrary source.
    E=sum(abs(x)**2 for x in z)
    c=sum((zi.conjugate()*fi).real for zi,fi in zip(z,raw))/E
    f=tuple(fi-c*zi for zi,fi in zip(z,raw))
    nu=0.19
    out=cartan_hminus1_gradient_split(a,z,f,nu)
    assert out["cartan_gradient_metric_cross"] == pytest.approx(0.0,abs=1e-12)
    assert out["ns_hminus1_action"] == pytest.approx(
        out["euler_hminus1_action"]+out["viscous_hminus1_action"]
    )


def test_normalized_cartan_triple_orientation_is_viscosity_invariant():
    out=normalized_triple_orientation_heat_law((0.7,2.1,4.3),(0.8,1.4,0.5),-0.37,0.23)
    assert out["raw_current_heat_rate"] != 0.0
    assert out["amplitude_root_heat_rate"] != 0.0
    assert out["normalized_orientation_heat_rate"] == pytest.approx(0.0,abs=1e-14)


def test_radial_mean_resolvent_square_is_exact_and_sharper_than_full_weighted_action():
    rng=random.Random(2026081423)
    for _ in range(800):
        a=tuple(rng.choice((-1,1))*10**rng.uniform(-2,2) for _ in range(6))
        e=tuple(10**rng.uniform(-3,3) for _ in a)
        p=[x/sum(e) for x in e]
        raw=[rng.uniform(-1,1) for _ in a]
        A01=sum(pi*ai for pi,ai in zip(p,a));A11=sum(pi*ai*ai for pi,ai in zip(p,a))
        b0=sum(pi*x for pi,x in zip(p,raw));b1=sum(pi*ai*x for pi,ai,x in zip(p,a,raw));det=A11-A01*A01
        if det<1e-14:continue
        c1=(b1-A01*b0)/det;c0=b0-c1*A01
        # Convert a weighted fitness to physical energy rates S=2 E f.
        rates=tuple(2.0*ei*(fi-c0-c1*ai) for ai,ei,fi in zip(a,e,raw))
        out=radial_mean_resolvent_balance(a,e,rates,0.17)
        assert out["mean_curl_rate"] == pytest.approx(out["represented_mean_curl_rate"],rel=2e-9,abs=2e-9)
        assert out["mean_curl_rate"] <= out["mean_curl_optimal_upper"] + 2e-9
        assert out["mean_curl_optimal_upper"] <= out["mean_curl_resolvent_upper"] + 2e-9
        assert 0.0 <= out["weighted_cauchy_fraction"] <= 1.0 + 1e-10


def test_master_sobolev_strain_transfer_has_energy_enstrophy_endpoints_and_critical_midpoint():
    rng=random.Random(2026081426)
    for _ in range(3000):
        r=10**rng.uniform(-6,6);s=10**rng.uniform(-6,6)
        for same in (False,True):
            e0=sobolev_strain_transfer_multiplier(0.0,r,s,same_helicity=same)
            em=sobolev_strain_transfer_multiplier(0.5,r,s,same_helicity=same)
            e1=sobolev_strain_transfer_multiplier(1.0,r,s,same_helicity=same)
            assert e0["strain_multiplier"] == pytest.approx(-1.0,abs=2e-12)
            assert e1["strain_multiplier"] == pytest.approx(1.0 if same else -1.0,abs=2e-12)
            if same:
                assert em["strain_multiplier"] == pytest.approx(0.0,abs=2e-12)
            else:
                ker=critical_logscale_strain_kernel(r,s)
                assert em["strain_multiplier"] == pytest.approx(-ker["strain_to_critical_multiplier"],rel=2e-12,abs=2e-12)


def test_critical_sech_kernel_is_sylvester_multiplier_and_loggap_collective_bound_decays():
    rng=random.Random(2026081427)
    for _ in range(5000):
        r=10**rng.uniform(-8,8);s=10**rng.uniform(-8,8)
        out=critical_logscale_strain_kernel(r,s)
        assert out["strain_to_critical_multiplier"] == pytest.approx(
            2.0*math.sqrt(r*s)/(r+s),rel=3e-12,abs=3e-12
        )
        assert 0.0 < out["strain_to_critical_multiplier"] <= 1.0 + 1e-12
    b2=critical_loggap_collective_bound(2.0)
    b5=critical_loggap_collective_bound(5.0)
    assert 0.0 < b5["collective_multiplier_bound"] < b2["collective_multiplier_bound"] <= 1.0
    assert b5["collective_multiplier_bound"] == pytest.approx(1.0/math.sinh(2.5))


def test_poisson_critical_scale_measure_is_probability_with_mean_inverse_mean_curl():
    rng=random.Random(2026081428)
    for _ in range(1000):
        a=tuple(rng.choice((-1,1))*10**rng.uniform(-3,3) for _ in range(7))
        e=tuple(10**rng.uniform(-4,4) for _ in a)
        out=poisson_critical_scale_measure_moments(a,e)
        assert out["probability_mass"] == pytest.approx(1.0,abs=2e-12)
        assert out["mean_poisson_scale"] == pytest.approx(
            1.0/(2.0*out["mean_absolute_curl"]),rel=2e-12,abs=2e-12
        )
        assert out["poisson_scale_variance"] >= -1e-12


def test_resolvent_mean_curl_action_is_never_weaker_than_previous_productive_fisher_upper():
    rng=random.Random(2026081429)
    for _ in range(1200):
        a=tuple(rng.choice((-1,1))*10**rng.uniform(-2,2) for _ in range(6))
        e=tuple(10**rng.uniform(-3,3) for _ in a)
        p=[x/sum(e) for x in e]
        raw=[rng.uniform(-1,1) for _ in a]
        A01=sum(pi*ai for pi,ai in zip(p,a));A11=sum(pi*ai*ai for pi,ai in zip(p,a))
        b0=sum(pi*x for pi,x in zip(p,raw));b1=sum(pi*ai*x for pi,ai,x in zip(p,a,raw));det=A11-A01*A01
        if det<1e-14:continue
        c1=(b1-A01*b0)/det;c0=b0-c1*A01
        rates=tuple(2.0*ei*(fi-c0-c1*ai) for ai,ei,fi in zip(a,e,raw))
        nu=0.23
        newer=radial_mean_resolvent_balance(a,e,rates,nu)
        older=radial_fitness_selection_balance(a,e,rates,nu)
        assert newer["log_mean_curl_optimal_upper"] <= older["physical_log_mean_curl_upper"] + 2e-8


def test_critical_midpoint_uniquely_minimizes_blockwise_symmetric_euler_growth_on_zero_to_one():
    rng=random.Random(2026081430)
    for _ in range(10000):
        r=10**rng.uniform(-7,7);q=10**rng.uniform(-7,7);s=rng.random()
        for same in (False,True):
            cur=sobolev_strain_transfer_multiplier(s,r,q,same_helicity=same)["strain_multiplier"]
            mid=sobolev_strain_transfer_multiplier(0.5,r,q,same_helicity=same)["strain_multiplier"]
            assert abs(mid) <= abs(cur) + 2e-12
            if same:
                assert mid == pytest.approx(0.0,abs=1e-13)


def test_helicity_polar_krein_generator_has_reciprocal_pseudounitary_singular_values():
    import numpy as np
    import scipy.linalg as la
    rng=np.random.default_rng(2026081433)
    for n in (4,8,12):
        signs=np.array(([1,-1]*(n//2)),float)
        J=np.diag(signs)
        # Build the general J-skew form: skew J-even compact part plus selfadjoint J-odd boost.
        R=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        A=(R-R.conj().T)/2
        A=(A+J@A@J)/2
        Q=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        S=(Q+Q.conj().T)/2
        S=(S-J@S@J)/2
        L=A+S
        assert np.linalg.norm(L.conj().T@J+J@L) <= 2e-12*max(1.0,np.linalg.norm(L))
        U=la.expm(0.17*L)
        assert np.linalg.norm(U.conj().T@J@U-J) <= 2e-11*max(1.0,np.linalg.norm(U)**2)
        sv=np.linalg.svd(U,compute_uv=False)
        assert np.max(np.abs(sv*sv[::-1]-1.0)) <= 2e-10


def test_critical_spectral_hminus_half_square_is_exact_and_critical_probability_action():
    rng=random.Random(2026081434)
    for _ in range(1200):
        a=tuple(rng.choice((-1,1))*10**rng.uniform(-2,2) for _ in range(6))
        e=tuple(10**rng.uniform(-3,3) for _ in a)
        # Produce a legitimate tangent fitness by removing its weighted affine part.
        E=sum(e);p=[x/E for x in e];raw=[rng.uniform(-1,1) for _ in a]
        A01=sum(pi*ai for pi,ai in zip(p,a));A11=sum(pi*ai*ai for pi,ai in zip(p,a))
        b0=sum(pi*x for pi,x in zip(p,raw));b1=sum(pi*ai*x for pi,ai,x in zip(p,a,raw));det=A11-A01*A01
        if det<1e-14:continue
        c1=(b1-A01*b0)/det;c0=b0-c1*A01
        f=tuple(fi-c0-c1*ai for ai,fi in zip(a,raw))
        rates=tuple(2.0*ei*fi for ei,fi in zip(e,f))
        nu=0.17
        out=critical_spectral_hminus_half_square(a,e,rates,nu)
        assert out["critical_rate"] == pytest.approx(out["represented_critical_rate"],rel=2e-9,abs=2e-9)
        assert out["critical_rate"] <= out["critical_scalar_optimal_upper"] + 2e-9
        assert out["critical_scalar_optimal_upper"] <= out["critical_rate_upper"] + 2e-9
        assert 0.0 <= out["critical_cauchy_fraction"] <= 1.0 + 1e-10
        K=sum(abs(ai)*ei for ai,ei in zip(a,e))
        pi=[abs(ai)*ei/K for ai,ei in zip(a,e)]
        expected=sum(pii*(fi/abs(ai))**2 for pii,fi,ai in zip(pi,f,a))
        assert out["critical_probability_fitness_action"] == pytest.approx(expected,rel=2e-10,abs=2e-10)


def test_master_sobolev_hilbert_square_uses_one_local_fitness_currency_for_all_s():
    rng=random.Random(2026081435)
    for _ in range(500):
        a=tuple(rng.choice((-1,1))*10**rng.uniform(-2,2) for _ in range(6))
        e=tuple(10**rng.uniform(-2,2) for _ in a)
        E=sum(e);p=[x/E for x in e];raw=[rng.uniform(-1,1) for _ in a]
        A01=sum(pi*ai for pi,ai in zip(p,a));A11=sum(pi*ai*ai for pi,ai in zip(p,a))
        b0=sum(pi*x for pi,x in zip(p,raw));b1=sum(pi*ai*x for pi,ai,x in zip(p,a,raw));det=A11-A01*A01
        if det<1e-14:continue
        c1=(b1-A01*b0)/det;c0=b0-c1*A01
        rates=tuple(2.0*ei*(fi-c0-c1*ai) for ai,ei,fi in zip(a,e,raw))
        for s in (0.0,0.25,0.5,0.73,1.0,1.4):
            out=sobolev_spectral_hilbert_square(a,e,rates,s,0.19)
            assert out["full_rate"] == pytest.approx(out["represented_rate"],rel=2e-9,abs=2e-9)
            assert out["full_rate"] <= out["rate_upper"] + 2e-9
            assert out["sobolev_probability_local_fitness_action"] >= -1e-12
        # s=0 nonlinear cross vanishes exactly by energy conservation.
        e0=sobolev_spectral_hilbert_square(a,e,rates,0.0,0.19)
        assert e0["nonlinear_half_rate"] == pytest.approx(0.0,abs=2e-10)
        # s=1/2 reproduces the dedicated critical spectral square.
        mid=sobolev_spectral_hilbert_square(a,e,rates,0.5,0.19)
        crit=critical_spectral_hminus_half_square(a,e,rates,0.19)
        assert mid["sobolev_stock"] == pytest.approx(crit["critical_stock"])
        assert mid["spectral_shifted_action"] == pytest.approx(crit["spectral_hminus_half_action"])


def test_closed_triad_critical_hminus_half_action_is_sharply_scale_free_by_energy():
    rng=random.Random(2026081440)
    for _ in range(5000):
        r0=10**rng.uniform(-2,2);r1=10**rng.uniform(-2,2);theta=rng.uniform(0.01,math.pi-0.01)
        r2=math.sqrt(r0*r0+r1*r1+2*r0*r1*math.cos(theta))
        signs=tuple(rng.choice((-1,1)) for _ in range(3))
        a=tuple(si*ri for si,ri in zip(signs,(r0,r1,r2)))
        e=tuple(10**rng.uniform(-4,4) for _ in range(3))
        out=closed_triad_critical_action_bound(a,e,phase_cosine_abs=rng.random())
        assert out["action_to_energy_critical_ratio"] <= 0.5 + 2e-10
        assert max(out["root_geometric_ratios"]) <= 1.0 + 2e-10
        assert out["log_critical_rate_upper_at_viscosity_one"] <= out["scale_free_log_rate_upper_at_viscosity_one"] + 2e-8


def test_closed_triad_one_half_action_constant_is_approached_at_low_high_high_boundary():
    # Homochiral radii (1,eps,1), with eps << delta << 1 in the modal energies,
    # approach the sharp root geometry and make the positive diagonal part of E*K negligible.
    vals=[]
    for eps,delta in ((1e-3,3e-2),(1e-5,3e-3),(1e-7,3e-4)):
        out=closed_triad_critical_action_bound((1.0,eps,1.0),(1e-12,1.0,delta))
        vals.append(out["action_to_energy_critical_ratio"])
    assert vals[0] < vals[1] < vals[2] < 0.5 + 1e-10
    assert vals[-1] > 0.49


def test_productive_and_viscous_volumes_give_exact_critical_reynolds_race_and_ns_scaling():
    # Start from one exact tangent spectral source.
    a=(-2.4,-0.7,1.1,3.2)
    e=(0.8,1.3,0.9,0.5)
    raw=(0.7,-0.4,0.5,-0.8)
    E=sum(e);p=[x/E for x in e];A01=sum(pi*ai for pi,ai in zip(p,a));A11=sum(pi*ai*ai for pi,ai in zip(p,a));b0=sum(pi*x for pi,x in zip(p,raw));b1=sum(pi*ai*x for pi,ai,x in zip(p,a,raw));det=A11-A01*A01;c1=(b1-A01*b0)/det;c0=b0-c1*A01
    rates=tuple(2.0*ei*(fi-c0-c1*ai) for ai,ei,fi in zip(a,e,raw))
    nu=0.23
    base=critical_spectral_hminus_half_square(a,e,rates,nu)
    assert base["productive_reynolds"]**2 == pytest.approx(base["viscous_to_productive_volume_ratio"])
    exact_normalized=base["critical_rate"]/(2.0*nu*base["third_absolute_curl_moment"])
    represented=base["productive_growth_sign"]*base["productive_reynolds"]-1.0
    assert exact_normalized == pytest.approx(represented)
    assert base["productive_action_volume"] + 1e-12 >= base["critical_action_volume"]

    # Exact Navier--Stokes dilation at one time: a->lambda a, e->lambda^-1 e,
    # Euler energy rate S->lambda S.  Both intrinsic volumes scale lambda^-3.
    for lam in (0.17,0.8,2.5,11.0):
        aa=tuple(lam*x for x in a);ee=tuple(x/lam for x in e);ss=tuple(lam*x for x in rates)
        out=critical_spectral_hminus_half_square(aa,ee,ss,nu)
        assert out["critical_action_volume"] == pytest.approx(base["critical_action_volume"]/lam**3,rel=3e-10,abs=3e-10)
        assert out["productive_action_volume"] == pytest.approx(base["productive_action_volume"]/lam**3,rel=3e-10,abs=3e-10)
        assert out["critical_viscous_volume"] == pytest.approx(base["critical_viscous_volume"]/lam**3,rel=3e-10,abs=3e-10)
        assert out["productive_reynolds"] == pytest.approx(base["productive_reynolds"],rel=3e-10,abs=3e-10)


def test_dimensionless_critical_reynolds_operator_is_selfadjoint_odd_and_exactly_competes_with_identity_heat():
    import numpy as np
    rng=np.random.default_rng(2026081444)
    for n in (4,8,12):
        r=np.exp(rng.uniform(-3,3,n));sgn=np.array(([1,-1]*(n//2)),float);J=np.diag(sgn);Lam=np.diag(r)
        M=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));X=M-M.conj().T
        Sigma=.5*np.diag(np.sqrt(r))@(X@J-J@X)@np.diag(np.sqrt(r))
        nu=.37
        R=(np.diag(1/r)@Sigma@np.diag(1/r))/nu
        assert np.linalg.norm(R-R.conj().T) <= 2e-12*max(1.0,np.linalg.norm(R))
        assert np.linalg.norm(R@J+J@R) <= 2e-12*max(1.0,np.linalg.norm(R))
        y=rng.normal(size=n)+1j*rng.normal(size=n);z=Lam@y
        kappa=np.vdot(y,Sigma@y).real;M3=np.vdot(z,z).real
        direct=2*kappa-2*nu*M3
        represented=2*nu*np.vdot(z,(R-np.eye(n))@z).real
        assert direct == pytest.approx(represented,rel=2e-11,abs=2e-11)
        ray=np.vdot(z,R@z).real/M3
        assert ray == pytest.approx(kappa/(nu*M3),rel=2e-11,abs=2e-11)
        ev=np.linalg.eigvalsh(R)
        assert np.max(np.abs(ev+ev[::-1])) <= 3e-10*max(1.0,np.max(np.abs(ev)))


def test_helicity_odd_reynolds_eigenvectors_are_neutral_and_top_heat_is_dirichlet_nonnegative():
    import numpy as np
    rng=np.random.default_rng(2026081446)
    for n in (6,10,14):
        J=np.diag(np.array(([1,-1]*(n//2)),float))
        Q=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));R=(Q+Q.conj().T)/2;R=(R-J@R@J)/2
        ev,V=np.linalg.eigh(R);lam=ev[-1];v=V[:,-1]
        if abs(lam)>1e-10:
            assert abs(np.vdot(v,J@v)) <= 3e-10
        # Three arbitrary commuting momentum generators (diagonal in a Fourier basis).
        Dvals=rng.integers(-4,5,size=(n,3)).astype(float)
        lap=np.zeros_like(R)
        rhs=0.0
        for j in range(3):
            D=np.diag(Dvals[:,j]);lap += D@(D@R-R@D)-(D@R-R@D)@D
            Dv=D@v;rhs += 2*np.vdot(Dv,(lam*np.eye(n)-R)@Dv).real
        lhs=np.vdot(v,lap@v).real
        assert lhs == pytest.approx(rhs,rel=3e-10,abs=3e-10)
        assert lhs >= -3e-10


def test_continuum_critical_operator_hilbert_schmidt_constant_is_exact_one_over_64():
    out=continuum_critical_operator_isometry_constant()
    assert out["raw_fixed_wavevector_integral"] == pytest.approx(math.pi**3/8.0,rel=1e-15)
    assert out["unitary_fourier_factor_squared"] == pytest.approx((2*math.pi)**-3,rel=1e-15)
    assert out["hilbert_schmidt_norm_squared_coefficient"] == pytest.approx(1.0/64.0,rel=2e-15)
    assert out["critical_norm_to_hs_isometry_factor"] == pytest.approx(8.0)


def test_continuum_hs_isometry_implies_paired_reynolds_capacity_and_small_data_radius():
    # Pure algebraic consequence of a self-adjoint J-odd Hilbert--Schmidt operator.
    import numpy as np
    rng=np.random.default_rng(2026081448);nu=.41
    for n in (6,10,16):
        J=np.diag(np.array(([1,-1]*(n//2)),float))
        Q=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n));Qc=(Q+Q.conj().T)/2;Qc=(Qc-J@Qc@J)/2
        K=64.0*np.vdot(Qc,Qc).real
        R=Qc/nu;ev=np.linalg.eigvalsh(R)
        pos=ev[ev>1e-10]
        assert K == pytest.approx(128.0*nu*nu*np.sum(pos*pos),rel=3e-10,abs=3e-10)
        assert np.linalg.norm(R,2) <= math.sqrt(K)/(8.0*math.sqrt(2.0)*nu)+3e-10
        assert int(np.sum(pos>1.0)) <= K/(128.0*nu*nu)+1e-10


def test_continuum_midpoint_operator_sobolev_scale_places_critical_stock_at_plain_hs_energy():
    expected={0.0:-0.25,0.5:0.0,1.0:0.25,1.5:0.5,2.0:0.75}
    for s,pow_ in expected.items():
        out=continuum_midpoint_operator_sobolev_dictionary(s)
        assert out["operator_laplacian_power"] == pytest.approx(pow_)
        assert out["operator_norm_squared_multiplier"] == pytest.approx(64.0)
    # Reflection around s=1/2 becomes opposite operator Sobolev powers.
    for delta in (0.1,0.35,0.8):
        lo=continuum_midpoint_operator_sobolev_dictionary(0.5-delta)["operator_laplacian_power"]
        hi=continuum_midpoint_operator_sobolev_dictionary(0.5+delta)["operator_laplacian_power"]
        assert hi == pytest.approx(-lo)



def test_certificate_connects_critical_reynolds_operator_to_the_graded_current_parent_without_overclaim():
    cert=theorem_certificate()
    assert "degree imbalance" in cert["graded_current_strain_parent"]
    assert "helicity-odd Gram imbalance" in cert["critical_gram_self_frustration"]
    assert "curvature floor" in cert["critical_curvature_floor"]
    assert "same critical metric" in cert["critical_projected_gauss_tax"]
    assert "two-by-two Gram projection" in cert["critical_two_null_gauss_tax"]
    assert "does not prove large-data regularity" in cert["graded_current_persistence_guard"]
    assert cert["global_regularity_claimed"] is False


def test_certificate_records_actual_Q_as_primitive_and_midpoint_as_reading():
    cert=theorem_certificate()
    assert "Q=nu delta+i_u" in cert["primitive_actual_current_operator"]
    assert "Q^2=nu(beta wedge)^*" in cert["primitive_current_curvature_square"]
    assert "chord/tangent contraction squares to zero" in cert["primitive_current_nilpotent_chords"]
    assert "2nu S" in cert["primitive_full_hodge_gradient"]
    assert "[Q*,Q^2]" in cert["primitive_intertwining_lamb"]
    assert "alpha -> nu beta" in cert["primitive_finite_current_chains"]
    assert "reading rather than the fundamental state operator" in cert["primitive_midpoint_status"]
    assert "persistent critical near-kernel" in cert["primitive_turning_frontier"]
    assert cert["global_regularity_claimed"] is False


def test_primitive_Q_adjacent_degree_intertwining_recovers_lamb_interaction():
    # Smooth low-bandwidth periodic referee for
    # H1(Q eta)-Q(H2 eta)=nu^2 d(omega.b)+nu b x (u x omega).
    # Frequencies are <=1 on a 16^3 grid, so all products in the identity stay
    # below Nyquist and the FFT calculation is an exact trigonometric referee.
    import numpy as np
    N=16
    rng=np.random.default_rng(202608151036)
    k=np.fft.fftfreq(N)*N
    kx,ky,kz=np.meshgrid(k,k,k,indexing="ij")
    K=np.stack([kx,ky,kz],axis=-1)
    k2=(K*K).sum(axis=-1)
    mask=(np.max(np.abs(K),axis=-1)<=1)&(k2>0)

    def fftv(v): return np.fft.fftn(v,axes=(0,1,2))
    def ifftv(vh): return np.fft.ifftn(vh,axes=(0,1,2)).real
    def divfree():
        a=rng.normal(size=(N,N,N,3))
        ah=fftv(a)*mask[...,None]
        dot=(ah*K).sum(-1)
        ah=ah-K*dot[...,None]/np.where(k2==0,1,k2)[...,None]
        ah[k2==0]=0
        return ifftv(ah)
    def dvec(v,j): return ifftv((1j*K[...,j])[...,None]*fftv(v))
    def grad(f):
        fh=np.fft.fftn(f)
        return np.stack([np.fft.ifftn(1j*K[...,j]*fh).real for j in range(3)],axis=-1)
    def lap(v): return ifftv(-k2[...,None]*fftv(v))
    def curl(v):
        vh=fftv(v); out=np.empty_like(vh)
        out[...,0]=1j*(ky*vh[...,2]-kz*vh[...,1])
        out[...,1]=1j*(kz*vh[...,0]-kx*vh[...,2])
        out[...,2]=1j*(kx*vh[...,1]-ky*vh[...,0])
        return ifftv(out)
    def strain(v):
        g=np.empty(v.shape[:-1]+(3,3))
        for j in range(3): g[..., :, j]=dvec(v,j)
        return .5*(g+np.swapaxes(g,-1,-2))

    nu=.37
    u=divfree(); b=divfree(); omega=curl(u); S=strain(u)
    Qeta=nu*curl(b)-np.cross(u,b)
    H1Q=-nu**2*lap(Qeta)+(u*u).sum(-1)[...,None]*Qeta+2*nu*np.einsum("...ij,...j->...i",S,Qeta)
    H2b=-nu**2*lap(b)+(u*u).sum(-1)[...,None]*b-2*nu*np.einsum("...ij,...j->...i",S,b)
    QH2=nu*curl(H2b)-np.cross(u,H2b)
    lhs=H1Q-QH2
    rhs=nu**2*grad((omega*b).sum(-1))+nu*np.cross(b,np.cross(u,omega))
    rel=np.linalg.norm(lhs-rhs)/max(1.0,np.linalg.norm(lhs),np.linalg.norm(rhs))
    assert rel < 2e-11


def test_primitive_Q_state_chords_are_square_zero_and_tangents_anticommute():
    # Exterior-algebra matrix referee for i_v^2=0 and {i_v,i_w}=0 on Lambda^* R^3.
    import numpy as np
    basis=[(),(0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)]
    index={b:i for i,b in enumerate(basis)}
    def contraction(v):
        M=np.zeros((8,8))
        for col,I in enumerate(basis):
            for pos,j in enumerate(I):
                J=I[:pos]+I[pos+1:]
                M[index[J],col]+=((-1.0)**pos)*v[j]
        return M
    v=np.array([.7,-1.2,.4]); w=np.array([-.3,.5,1.1])
    Iv=contraction(v); Iw=contraction(w)
    assert np.linalg.norm(Iv@Iv) < 1e-14
    assert np.linalg.norm(Iw@Iw) < 1e-14
    assert np.linalg.norm(Iv@Iw+Iw@Iv) < 1e-14



def test_primitive_critical_channel_exact_continuum_constants_and_factor_two_heat_identity():
    out=continuum_primitive_critical_channel_constants()
    assert out["raw_translation_integral"] == pytest.approx(4.0*math.pi,rel=0,abs=1e-14)
    assert out["one_form_creation_graded_multiplicity"] == pytest.approx(4.0)
    assert out["two_form_creation_graded_multiplicity"] == pytest.approx(2.0)
    assert out["critical_channel_hs_coefficient"] == pytest.approx(2.0/math.pi**2,rel=1e-14)
    assert out["curvature_channel_hs_coefficient"] == pytest.approx(1.0/math.pi**2,rel=1e-14)
    assert out["translation_dirichlet_hs_coefficient"] == pytest.approx(2.0/math.pi**2,rel=1e-14)
    assert out["dirichlet_to_curvature_hs_ratio"] == pytest.approx(2.0)


def test_primitive_critical_channel_anticommutator_and_operator_derivative_are_exact():
    import numpy as np
    rng=np.random.default_rng(2026081551)
    # The universal algebra: A=[L,Q*] implies {A,Q*}=[L,(Q*)^2].
    for n in (4,7,11):
        lam=np.exp(rng.uniform(-2,2,n))
        L=np.diag(1.0/lam)
        Qs=rng.normal(size=(n,n))+1j*rng.normal(size=(n,n))
        A=L@Qs-Qs@L
        lhs=A@Qs+Qs@A
        rhs=L@(Qs@Qs)-(Qs@Qs)@L
        assert np.linalg.norm(lhs-rhs) <= 3e-12*max(1.0,np.linalg.norm(rhs))

    # One translation block on the 3D exterior algebra checks
    # B=i sum_j (dx^j wedge)[D_j,A] with beta_q=i q wedge alpha_q.
    basis=[tuple(i for i in range(3) if mask>>i&1) for mask in range(8)]
    index={b:i for i,b in enumerate(basis)}
    def wedge1(a):
        W=np.zeros((8,8),complex)
        for col,b in enumerate(basis):
            for i,ai in enumerate(a):
                if i in b: continue
                sign=(-1)**sum(x<i for x in b)
                W[index[tuple(sorted((i,)+b))],col]+=sign*ai
        return W
    e=[np.eye(3)[j] for j in range(3)]
    E=[wedge1(x) for x in e]
    for _ in range(40):
        q=rng.normal(size=3)
        alpha=rng.normal(size=3)+1j*rng.normal(size=3)
        # incompressible velocity coefficient
        alpha=alpha-q*np.vdot(q,alpha)/(np.dot(q,q))
        ell=rng.normal(size=3)
        if np.linalg.norm(ell)<0.3: ell[0]+=1.0
        k=ell+q
        if np.linalg.norm(k)<0.3: k[1]+=1.0
        diff=1.0/np.linalg.norm(k)-1.0/np.linalg.norm(ell)
        Ablock=diff*wedge1(alpha)
        beta=1j*sum(q[j]*(E[j]@wedge1(alpha)) for j in range(3))
        Bblock=diff*beta
        derived=1j*sum(E[j]@(q[j]*Ablock) for j in range(3))
        assert np.linalg.norm(Bblock-derived) <= 4e-12*max(1.0,np.linalg.norm(Bblock))


def test_certificate_records_primitive_critical_channel_as_self_dissipative_derivative_without_closure_claim():
    cert=theorem_certificate()
    assert "Q*1=alpha" in cert["primitive_vacuum_state_chain"]
    assert "current/metric commutator" in cert["primitive_metric_current_law"]
    assert "A=[Lambda^-1,Q*]" in cert["primitive_critical_two_way_channel"]
    assert "2/pi^2" in cert["primitive_critical_channel_isometry"]
    assert "{A,Q*}" in cert["primitive_critical_channel_derivative"]
    assert "not a uniform instantaneous gap" in cert["primitive_critical_channel_rigidity"]
    assert "A_t=A(F_E)-nu Delta_op A" in cert["primitive_critical_channel_dynamics"]
    assert "infinite productive regeneration" in cert["primitive_turning_frontier"]
    assert cert["global_regularity_claimed"] is False


def test_primitive_critical_carre_du_champ_constants_and_cauchy_binet_are_exact():
    import itertools
    import numpy as np

    out=continuum_critical_carre_du_champ_constants()
    assert out["fractional_laplacian_kernel_coefficient"] == pytest.approx(1.0/math.pi**2,rel=1e-15)
    assert out["gauss_hs_to_vorticity_metric_factor"] == pytest.approx(1.0/(4.0*math.pi**2),rel=1e-15)
    assert out["second_elementary_symmetric_prefactor"] == pytest.approx(1.0/(2.0*math.pi**4),rel=1e-15)
    assert out["determinant_prefactor"] == pytest.approx(1.0/(6.0*math.pi**6),rel=1e-15)

    v=[np.array([1.0,0.2,-0.1]),np.array([0.3,-0.7,0.4]),np.array([-0.2,0.5,0.9]),np.array([0.6,0.1,0.2])]
    w=[0.4,0.7,0.3,0.6]
    G=sum((wi/math.pi**2)*np.outer(vi,vi) for wi,vi in zip(w,v))
    assert np.linalg.eigvalsh(G).min() >= -2e-15
    e2=0.5*(np.trace(G)**2-np.trace(G@G))
    e2_cb=sum(w[i]*w[j]*np.linalg.norm(np.cross(v[i],v[j]))**2
              for i,j in itertools.product(range(len(v)),repeat=2))/(2.0*math.pi**4)
    det_cb=sum(w[i]*w[j]*w[k]*np.linalg.det(np.stack([v[i],v[j],v[k]],axis=1))**2
               for i,j,k in itertools.product(range(len(v)),repeat=3))/(6.0*math.pi**6)
    assert e2 == pytest.approx(e2_cb,rel=2e-14,abs=2e-14)
    assert np.linalg.det(G) == pytest.approx(det_cb,rel=2e-14,abs=2e-14)


def test_critical_carre_du_champ_heat_product_law_has_negative_positive_sink():
    import numpy as np

    n=64
    x=2.0*math.pi*np.arange(n)/n
    k=np.fft.fftfreq(n,1.0/n)
    def op(a,mult):
        return np.fft.ifft(mult*np.fft.fft(a)).real
    def lam(a): return op(a,np.abs(k))
    def lap(a): return op(a,k*k)
    def der(a): return op(a,1j*k)
    def gam(a,b): return a*lam(b)+b*lam(a)-lam(a*b)
    f=np.sin(x)+0.23*np.cos(2*x)
    g=0.7*np.cos(x)-0.19*np.sin(3*x)
    lhs=lap(gam(f,g))-gam(lap(f),g)-gam(f,lap(g))
    rhs=-2.0*gam(der(f),der(g))
    assert np.max(np.abs(lhs-rhs)) <= 2e-11


def test_certificate_collapses_gauss_source_to_intrinsic_positive_critical_metric_without_claiming_bridge():
    cert=theorem_certificate()
    assert "int tr Gamma_u=2K" in cert["primitive_critical_carre_du_champ"]
    assert "omega^T Gamma_u omega" in cert["primitive_critical_vorticity_metric"]
    assert "Loewner-positive" in cert["primitive_critical_metric_heat_law"]
    assert "finite-energy R3 has only vacuum" in cert["primitive_critical_rank_persistence"]
    assert "proposed factor 1/2 is numerically false" in cert["primitive_critical_stretching_bridge_guard"]
    assert "factor 1 from below without a proof" in cert["primitive_critical_stretching_bridge_guard"]
    assert cert["global_regularity_claimed"] is False



def test_primitive_two_particle_roads_and_inversion_geometry_are_exact():
    import numpy as np

    # Affine divergence-free state: the local pair law reduces to exact algebra,
    # including the fact that pressure enters only through the center road.
    A=np.array([[0.2,0.7,-0.1],[-0.3,-0.4,0.5],[0.6,-0.2,0.2]],float)
    assert abs(np.trace(A)) < 1e-15
    P=np.array([[0.4,-0.2,0.1],[-0.2,0.3,0.05],[0.1,0.05,-0.7]],float)
    r=np.array([0.8,-0.6,0.9])
    du=A@r
    q=0.5*float(du@du)
    dut=-(A@A)@r-P@r
    q_t=float(du@dut)
    center_pressure=float((P@r)@du)
    relative_flux=float(du@(A.T@du))+q*float(np.trace(A))
    nu=0.37
    viscous=nu*(2.0*float(np.trace(A.T@A))-2.0*float(np.sum(A*A)))
    assert q_t+center_pressure+relative_flux == pytest.approx(viscous,abs=2e-14)

    # Inversion is a scaled reflection; its Hessian cubic is exactly the
    # critical inward-compression integrand.
    rr=np.array([1.1,-0.7,0.9])
    v=np.array([0.4,0.8,-0.3])
    rho=float(np.linalg.norm(rr)); e=rr/rho
    DI=(np.eye(3)-2.0*np.outer(e,e))/rho**2
    a=float(v@e); w=v-a*e
    D2=(-4.0*a*w+2.0*(a*a-float(w@w))*e)/rho**3
    assert float(np.linalg.norm(DI@v)) == pytest.approx(float(np.linalg.norm(v))/rho**2,rel=2e-15)
    assert float(np.linalg.norm(D2)) == pytest.approx(2.0*float(v@v)/rho**3,rel=2e-15)
    assert float((DI@v)@D2) == pytest.approx(-2.0*float(v@v)*float(v@rr)/rho**6,rel=2e-15)

    # Common/relative endpoint algebra has no third energy reservoir.
    up=np.array([1.2,-0.4,0.7]); um=np.array([-0.3,0.9,0.2])
    U=0.5*(up+um); vv=up-um
    assert float(U@vv) == pytest.approx(0.5*(float(up@up)-float(um@um)),rel=2e-15)
    assert float(U@U)+0.25*float(vv@vv) == pytest.approx(0.5*(float(up@up)+float(um@um)),rel=2e-15)

    # Endpoint chain rule: A=(Ax+Ay)/2, C=Ax-Ay and the square/pressure source
    # differences are homogeneous in the affine defect C.
    Ax=np.array([[0.4,0.2,-0.1],[-0.3,-0.1,0.5],[0.2,-0.4,-0.3]])
    Ay=np.array([[-0.2,0.1,0.3],[0.4,0.2,-0.1],[-0.5,0.2,0.0]])
    assert abs(np.trace(Ax)) < 1e-15 and abs(np.trace(Ay)) < 1e-15
    AA=0.5*(Ax+Ay); C=Ax-Ay
    assert np.max(np.abs((Ax@Ax-Ay@Ay)-(AA@C+C@AA))) <= 3e-16
    assert np.trace(Ax@Ax)-np.trace(Ay@Ay) == pytest.approx(2.0*np.trace(AA@C),abs=3e-16)


def test_primitive_pair_radial_transverse_split_and_pressure_converter_constants_are_exact():
    # For k=e3 and a divergence-free polarization e1, radial integration of
    # (1-cos(k.r))/|r|^4 leaves |cos(theta)| on S^2.  The projected angular
    # fraction is exactly 1/4; the complementary transverse fraction is 3/4.
    total_ang=2.0*math.pi  # 2pi * int_-1^1 |mu| dmu
    parallel_ang=0.5*math.pi  # pi * int_-1^1 (1-mu^2)|mu| dmu
    assert parallel_ang/total_ang == pytest.approx(0.25,rel=1e-15)
    assert 1.0-parallel_ang/total_ang == pytest.approx(0.75,rel=1e-15)

    # Differentiate a^2 |r|^-4 along dot r=a n+b.  The non-pressure Euler
    # contribution to K_parallel' is pi^-2(B3-2A3).  Since K_parallel=K/4
    # and K'_E=2 kappa, pressure is forced to be the exact converter below.
    A3,B3=1.7,-0.43
    kappa=-(A3+B3)/math.pi**2
    adv_parallel=(B3-2.0*A3)/math.pi**2
    pressure_parallel=0.5*kappa-adv_parallel
    expected=3.0*(A3-B3)/(2.0*math.pi**2)
    assert pressure_parallel == pytest.approx(expected,rel=1e-15)
    assert pressure_parallel + (-pressure_parallel) == 0.0

    # div_r(delta u)=1/2 div u(x)+1/2 div u(y).
    tr_x=(0.7-0.2-0.5)
    tr_y=(-0.4+0.9-0.5)
    assert 0.5*(tr_x+tr_y) == pytest.approx(0.0,abs=1e-15)


def test_poisson_energy_boundary_jets_and_hankel_escape_action_are_exact():
    E,K,Z,M3,kappa,nu=3.2,2.1,2.7,4.4,-0.63,0.18
    assert E*Z-K*K > 0.0
    out=poisson_energy_boundary_jet_geometry(E,K,Z,M3,kappa,nu)
    assert out["energy_profile_value"] == pytest.approx(E)
    assert out["energy_profile_first_derivative"] == pytest.approx(-2.0*K)
    assert out["energy_profile_second_derivative"] == pytest.approx(4.0*Z)
    assert out["energy_profile_third_derivative"] == pytest.approx(-8.0*M3)
    assert out["euler_profile_value"] == 0.0
    assert out["euler_profile_first_derivative"] == pytest.approx(-4.0*kappa)
    assert out["hankel_determinant"] == pytest.approx(4.0*(E*Z-K*K),rel=2e-15)
    assert out["energy_rate"] == pytest.approx(-2.0*nu*Z,rel=2e-15)
    assert out["critical_rate"] == pytest.approx(2.0*kappa-2.0*nu*M3,rel=2e-15)
    expected=kappa*kappa/((Z/E)*(E*Z-K*K))
    assert out["escape_action_from_profile"] == pytest.approx(expected,rel=2e-15)
    assert out["escape_action_physical"] == pytest.approx(expected,rel=2e-15)


def test_poisson_scalar_counterprofile_is_zero_net_but_not_physically_realizable():
    out=poisson_scalar_reversible_counterprofile(0.037,0.31,2.4,0.12)
    assert out["energy_loss_integrability_power"] > 0.0
    assert abs(out["hankel_determinant"]) <= 2e-12*max(1.0,out["critical_stock"]**2)
    assert abs(out["reversible_source_area"]) <= 2e-15*max(1.0,abs(out["energy_rate"]))
    assert out["fake_reversible_boundary_current"] == pytest.approx(
        2.0*out["fake_curvature_height"],rel=2e-15
    )
    assert out["fake_reversible_boundary_current"] > 0.0
    assert out["physical_equiradial_curvature_height"] == 0.0
    assert out["fake_curvature_height"] > out["physical_equiradial_curvature_height"]


def test_poisson_product_defect_realizes_depth_euler_work_and_boundary_derivative():
    import numpy as np

    n=12
    rng=np.random.default_rng(20260815)
    ks=np.fft.fftfreq(n,1.0/n)
    kx,ky,kz=np.meshgrid(ks,ks,ks,indexing="ij")
    k2=kx*kx+ky*ky+kz*kz
    radius=np.sqrt(k2)
    mask=(np.maximum.reduce([np.abs(kx),np.abs(ky),np.abs(kz)])<=1)&(k2>0)

    u0=rng.normal(size=(3,n,n,n))
    uh=np.fft.fftn(u0,axes=(1,2,3))*mask[None]
    den=np.where(k2>0,k2,1.0)
    dot=kx*uh[0]+ky*uh[1]+kz*uh[2]
    for j,kj in enumerate((kx,ky,kz)):
        uh[j]-=np.where(k2>0,kj*dot/den,0.0)
    u=np.fft.ifftn(uh,axes=(1,2,3)).real

    def multiplier(v,m):
        vh=np.fft.fftn(v,axes=(1,2,3))
        return np.fft.ifftn(vh*m,axes=(1,2,3)).real

    def py(v,y): return multiplier(v,np.exp(-y*radius))
    def lam(v): return multiplier(v,radius)
    def curl(v):
        vh=np.fft.fftn(v,axes=(1,2,3))
        wh=np.empty_like(vh)
        wh[0]=1j*(ky*vh[2]-kz*vh[1])
        wh[1]=1j*(kz*vh[0]-kx*vh[2])
        wh[2]=1j*(kx*vh[1]-ky*vh[0])
        return np.fft.ifftn(wh,axes=(1,2,3)).real
    def leray(v):
        vh=np.fft.fftn(v,axes=(1,2,3))
        dd=kx*vh[0]+ky*vh[1]+kz*vh[2]
        out=vh.copy()
        for j,kj in enumerate((kx,ky,kz)):
            out[j]-=np.where(k2>0,kj*dd/den,0.0)
        out[:,0,0,0]=0.0
        return np.fft.ifftn(out,axes=(1,2,3)).real
    def cross(a,b): return np.cross(a,b,axisa=0,axisb=0,axisc=0)
    def inner(a,b): return float(np.mean(np.sum(a*b,axis=0)))

    omega=curl(u)
    lamb=cross(u,omega)
    FE=leray(lamb)
    kappa=inner(lam(u),FE)
    assert abs(kappa) > 1e-8

    y=0.29
    uy=py(u,y); wy=py(omega,y)
    Dy=py(lamb,y)-cross(uy,wy)
    RE_direct=2.0*inner(uy,py(FE,y))
    RE_defect=2.0*inner(uy,leray(Dy))
    same_depth=inner(uy,cross(uy,wy))
    assert abs(same_depth) <= 2e-13*max(1.0,abs(RE_direct))
    assert RE_defect == pytest.approx(RE_direct,rel=3e-12,abs=3e-12)

    # Simpson audit of the exact integral D_y=int P_(y-s) Gamma_cross(P_su,P_somega) ds.
    count=81
    ss=np.linspace(0.0,y,count)
    vals=[]
    for sv in ss:
        us=py(u,sv); ws=py(omega,sv)
        gam=cross(us,lam(ws))+cross(lam(us),ws)-lam(cross(us,ws))
        vals.append(py(gam,y-sv))
    vals=np.stack(vals)
    weights=np.ones(count); weights[1:-1:2]=4.0; weights[2:-1:2]=2.0
    Dint=(y/(count-1))/3.0*np.tensordot(weights,vals,axes=(0,0))
    assert np.linalg.norm(Dy-Dint) <= 2e-9*max(1.0,np.linalg.norm(Dy))

    gamma0=cross(u,lam(omega))+cross(lam(u),omega)-lam(lamb)
    Rprime_gamma=2.0*inner(u,gamma0)
    Rprime_semigroup=-2.0*inner(lam(u),FE)-2.0*inner(u,lam(FE))
    assert Rprime_gamma == pytest.approx(-4.0*kappa,rel=3e-12,abs=3e-12)
    assert Rprime_semigroup == pytest.approx(-4.0*kappa,rel=3e-12,abs=3e-12)



def test_poisson_covariance_stress_two_reservoir_heat_and_fisher_laws_are_exact():
    import numpy as np

    n=12
    rng=np.random.default_rng(20260815)
    ks=np.fft.fftfreq(n,1.0/n)
    kx,ky,kz=np.meshgrid(ks,ks,ks,indexing="ij")
    k2=kx*kx+ky*ky+kz*kz
    radius=np.sqrt(k2)
    den=np.where(k2>0.0,k2,1.0)

    # Random real low-band state; products remain below Nyquist, so the tensor
    # and Lamb identities are not testing an aliasing artifact.
    raw=rng.normal(size=(3,n,n,n))
    uh=np.fft.fftn(raw,axes=(1,2,3))*(radius<=2.2)[None]
    dot=kx*uh[0]+ky*uh[1]+kz*uh[2]
    for j,kj in enumerate((kx,ky,kz)):
        uh[j]-=np.where(k2>0.0,kj*dot/den,0.0)
    uh[:,0,0,0]=0.0
    u=np.fft.ifftn(uh,axes=(1,2,3)).real
    u/=math.sqrt(float(np.mean(np.sum(u*u,axis=0))))

    def spatial_axes(a): return tuple(range(a.ndim-3,a.ndim))
    def fft(a): return np.fft.fftn(a,axes=spatial_axes(a))
    def ifft(h): return np.fft.ifftn(h,axes=spatial_axes(h)).real
    def py(a,y): return ifft(fft(a)*np.exp(-y*radius))
    def lam(a): return ifft(fft(a)*radius)
    def lap(a): return ifft(fft(a)*k2)
    def tensor(a,b): return np.einsum("i...,j...->ij...",a,b)
    def cross(a,b): return np.cross(a,b,axisa=0,axisb=0,axisc=0)
    def inner(a,b): return float(np.mean(np.sum(a*b,axis=0)))
    def tensor_inner(a,b): return float(np.mean(np.sum(a*b,axis=(0,1))))
    def curl(a):
        ah=fft(a)
        wh=np.empty_like(ah)
        wh[0]=1j*(ky*ah[2]-kz*ah[1])
        wh[1]=1j*(kz*ah[0]-kx*ah[2])
        wh[2]=1j*(kx*ah[1]-ky*ah[0])
        return ifft(wh)
    def leray(a):
        ah=fft(a); dd=kx*ah[0]+ky*ah[1]+kz*ah[2]
        out=ah.copy()
        for j,kj in enumerate((kx,ky,kz)):
            out[j]-=np.where(k2>0.0,kj*dd/den,0.0)
        out[:,0,0,0]=0.0
        return ifft(out)
    def grad_scalar(a):
        ah=fft(a)
        return np.stack([ifft(1j*kx*ah),ifft(1j*ky*ah),ifft(1j*kz*ah)])
    def grad_vector(a):
        ah=fft(a)
        return np.stack([np.stack([ifft(1j*kj*ah[i]) for i in range(3)]) for kj in (kx,ky,kz)])
    def div_vector(a):
        ah=fft(a)
        return ifft(1j*(kx*ah[0]+ky*ah[1]+kz*ah[2]))
    def div_tensor(a):
        ah=fft(a)
        return ifft(1j*(kx*ah[:,0]+ky*ah[:,1]+kz*ah[:,2]))
    def strain(a):
        ah=fft(a); G=np.empty((3,3,n,n,n),float)
        for i in range(3):
            for j,kj in enumerate((kx,ky,kz)):
                G[i,j]=ifft(1j*kj*ah[i])
        return 0.5*(G+np.swapaxes(G,0,1))

    omega=curl(u)
    y=0.31; z=0.23; nu=0.41
    v=py(u,y); wy=py(omega,y)
    tau=py(tensor(u,u),y)-tensor(v,v)
    trtau=np.trace(tau,axis1=0,axis2=1)

    # The canonical Poisson lift is a local harmonic half-space system.
    vy=-lam(v); vyy=lap(v)
    assert np.linalg.norm(vyy-lap(v)) <= 2e-13*max(1.0,np.linalg.norm(vyy))
    tauyy=py(lap(tensor(u,u)),y)-tensor(vyy,v)-2.0*tensor(vy,vy)-tensor(v,vyy)
    minus_delta4_tau=lap(tau)-tauyy
    gv=grad_vector(v)
    dirichlet_source=tensor(vy,vy)
    for j in range(3):
        dirichlet_source += tensor(gv[j],gv[j])
    assert np.linalg.norm(minus_delta4_tau-2.0*dirichlet_source) <= 8e-12*max(1.0,np.linalg.norm(dirichlet_source))
    source_eigs=np.linalg.eigvalsh(np.moveaxis(minus_delta4_tau,(0,1),(-2,-1)))
    assert float(source_eigs.min()) >= -3e-12

    source_C=np.eye(3)[:,:,None,None,None]*np.trace(minus_delta4_tau,axis1=0,axis2=1)-minus_delta4_tau
    direct_C=np.zeros_like(source_C)
    for g in (vy,gv[0],gv[1],gv[2]):
        g2=np.sum(g*g,axis=0)
        direct_C += 2.0*(np.eye(3)[:,:,None,None,None]*g2-tensor(g,g))
    assert np.linalg.norm(source_C-direct_C) <= 8e-12*max(1.0,np.linalg.norm(direct_C))
    axis=np.array([0.4,-0.7,1.1]); axis/=np.linalg.norm(axis)
    contracted=np.einsum("i,ij...,j->...",axis,source_C,axis)
    cross_source=np.zeros_like(contracted)
    for g in (vy,gv[0],gv[1],gv[2]):
        cross_source += 2.0*np.sum(np.cross(axis[:,None,None,None],g,axisa=0,axisb=0,axisc=0)**2,axis=0)
    assert np.max(np.abs(contracted-cross_source)) <= 8e-12*max(1.0,float(np.max(np.abs(cross_source))))

    # Positive Markov covariance and exact two-depth Germano cocycle.
    eig=np.linalg.eigvalsh(np.moveaxis(tau,(0,1),(-2,-1)))
    assert float(eig.min()) >= -2e-12
    tau_yz=py(tensor(u,u),y+z)-tensor(py(u,y+z),py(u,y+z))
    tau_z_v=py(tensor(v,v),z)-tensor(py(v,z),py(v,z))
    cocycle=py(tau,z)+tau_z_v
    assert np.linalg.norm(tau_yz-cocycle) <= 3e-12*max(1.0,np.linalg.norm(tau_yz))

    # (partial_y+Lambda)tau=Gamma_v.
    partial_y_tau=-lam(py(tensor(u,u),y))+tensor(lam(v),v)+tensor(v,lam(v))
    gamma_v=tensor(v,lam(v))+tensor(lam(v),v)-lam(tensor(v,v))
    assert np.linalg.norm(partial_y_tau+lam(tau)-gamma_v) <= 4e-12*max(1.0,np.linalg.norm(gamma_v))

    # Product defect is exactly the covariance stress divergence.
    lamb=cross(u,omega)
    D=py(lamb,y)-cross(v,wy)
    stress_D=0.5*grad_scalar(trtau)-div_tensor(tau)
    assert np.linalg.norm(D-stress_D) <= 5e-12*max(1.0,np.linalg.norm(D))
    filtered=py(leray(lamb),y)
    stress_rhs=leray(cross(v,curl(v)))-leray(div_tensor(tau))
    assert np.linalg.norm(filtered-stress_rhs) <= 5e-12*max(1.0,np.linalg.norm(filtered))

    # The pressure version removes Leray from the bulk equation.
    raw0=-div_tensor(tensor(u,u))
    div_raw0_h=fft(div_vector(raw0))
    p0_h=np.where(k2>0.0,-div_raw0_h/den,0.0)
    p0=ifft(p0_h); pi=py(p0,y)
    m=tau+tensor(v,v)
    local_euler=-div_tensor(m)-grad_scalar(pi)
    assert np.linalg.norm(local_euler-filtered) <= 7e-12*max(1.0,np.linalg.norm(filtered))
    assert lap(pi) == pytest.approx(div_vector(div_tensor(m)),rel=8e-12,abs=8e-12)

    # Exactly two scalar reservoirs and a genuinely nonzero reversible exchange.
    E=inner(u,u); U=inner(v,v); V=float(np.mean(trtau))
    assert U+V == pytest.approx(E,rel=3e-14,abs=3e-14)
    Ut_euler=2.0*inner(v,py(leray(lamb),y))
    stress_work=2.0*tensor_inner(tau,strain(v))
    assert abs(Ut_euler) > 1e-5
    assert Ut_euler == pytest.approx(stress_work,rel=3e-12,abs=3e-12)

    Z=inner(omega,omega); Zy=inner(wy,wy)
    assert -2.0*nu*Zy < 0.0
    assert -2.0*nu*(Z-Zy) < 0.0

    # Pure heat Loewner law and exact unresolved gradient-variance bill.
    ut=-nu*lap(u); vt=py(ut,y)
    tau_t=py(tensor(ut,u)+tensor(u,ut),y)-tensor(vt,v)-tensor(v,vt)
    heat_lhs=tau_t+nu*lap(tau)
    heat_rhs=np.zeros_like(tau)
    uh2=fft(u)
    for kj in (kx,ky,kz):
        du=ifft(1j*kj*uh2); dv=py(du,y)
        heat_rhs += -2.0*nu*(py(tensor(du,du),y)-tensor(dv,dv))
    assert np.linalg.norm(heat_lhs-heat_rhs) <= 8e-12*max(1.0,np.linalg.norm(heat_rhs))
    grad_cov=-float(np.mean(np.trace(heat_rhs,axis1=0,axis2=1)))/(2.0*nu)
    assert grad_cov == pytest.approx(Z-Zy,rel=3e-13,abs=3e-13)

    # Constant-one covariance Fisher estimate.
    cov_omega=py(np.sum(omega*omega,axis=0),y)-np.sum(wy*wy,axis=0)
    d2=np.sum(D*D,axis=0)
    assert np.max(d2-trtau*cov_omega) <= 2e-11*max(1.0,float(np.max(trtau*cov_omega)))
    ratio=np.zeros_like(trtau); live=trtau>1e-13
    ratio[live]=d2[live]/trtau[live]
    fisher=float(np.mean(ratio))
    assert fisher <= (Z-Zy)+2e-12

    # The product defect has a canonical covariance least-squares axis.
    C=np.eye(3)[:,:,None,None,None]*trtau-tau
    activity=float(np.mean(np.einsum("i...,ij...,j...->...",v,C,v)))
    RE=2.0*inner(v,D)
    Cm=np.moveaxis(C,(0,1),(-2,-1)); Dm=np.moveaxis(D,0,-1)
    ceig=np.linalg.eigvalsh(Cm)
    assert float(ceig.min()) > 1e-5  # regular-rank referee for the axis-gradient identity below
    am=np.einsum("...ij,...j->...i",np.linalg.inv(Cm),Dm); a=np.moveaxis(am,-1,0)
    assert np.linalg.norm(D-np.einsum("ij...,j...->i...",C,a)) <= 4e-12*max(1.0,np.linalg.norm(D))

    axis_density=np.einsum("i...,ij...,j...->...",a,C,a)
    reg_density=cov_omega-axis_density
    assert float(reg_density.min()) >= -2e-11
    axis=float(np.mean(axis_density)); Rreg=float(np.mean(reg_density))
    beta=RE/(2.0*activity)
    align=float(np.mean(np.einsum("i...,ij...,j...->...",a-beta*v,C,a-beta*v)))
    productive=RE*RE/(4.0*activity)
    assert axis == pytest.approx(productive+align,rel=4e-12,abs=4e-12)
    assert Z-Zy == pytest.approx(productive+align+Rreg,rel=4e-12,abs=4e-12)
    assert RE*RE <= 4.0*activity*(Z-Zy)+2e-12

    # Exact relative-reservoir completed square: one positive term, three nonpositive terms.
    Vt=-RE-2.0*nu*(Z-Zy)
    Vt_square=(activity/(2.0*nu)
               -nu*(RE+activity/nu)**2/(2.0*activity)
               -2.0*nu*align-2.0*nu*Rreg)
    assert Vt_square == pytest.approx(Vt,rel=4e-12,abs=4e-12)

    # Pointwise rank-one e2 identity and Maclaurin depth-area integrand bound.
    Mfull=tau+tensor(v,v)
    def e2(T):
        tr=np.trace(T,axis1=0,axis2=1)
        return 0.5*(tr*tr-np.einsum("ij...,ji...->...",T,T))
    point_activity=np.einsum("i...,ij...,j...->...",v,C,v)
    assert np.max(np.abs(point_activity-(e2(Mfull)-e2(tau)))) <= 5e-12
    assert np.max(point_activity-(np.trace(Mfull,axis1=0,axis2=1)**2)/3.0) <= 5e-12

    # On this positive-rank periodic referee the conditional axis integration-by-parts sign is exact.
    ah=fft(a); grad_a=np.empty((3,3,n,n,n),float)
    for i in range(3):
        for j,kj in enumerate((kx,ky,kz)):
            grad_a[i,j]=ifft(1j*kj*ah[i])
    stress=tau-0.5*np.eye(3)[:,:,None,None,None]*trtau
    axis_ibp=float(np.mean(np.sum(stress*grad_a,axis=(0,1))))
    assert axis_ibp == pytest.approx(axis,rel=6e-12,abs=6e-12)

    # Germano also splits transverse activity into inherited/new nonnegative pieces.
    w=py(v,z)
    inherited_tau=py(tau,z)
    inherited_M=np.eye(3)[:,:,None,None,None]*np.trace(inherited_tau,axis1=0,axis2=1)-inherited_tau
    inherited_activity=float(np.mean(np.einsum("i...,ij...,j...->...",w,inherited_M,w)))
    new_M=np.eye(3)[:,:,None,None,None]*np.trace(tau_z_v,axis1=0,axis2=1)-tau_z_v
    new_activity=float(np.mean(np.einsum("i...,ij...,j...->...",w,new_M,w)))
    tau_total=tau_yz
    total_M=np.eye(3)[:,:,None,None,None]*np.trace(tau_total,axis1=0,axis2=1)-tau_total
    total_activity=float(np.mean(np.einsum("i...,ij...,j...->...",w,total_M,w)))
    assert inherited_activity >= -2e-13
    assert new_activity >= -2e-13
    assert total_activity == pytest.approx(inherited_activity+new_activity,rel=3e-12,abs=3e-12)

    # Boundary jets: tau/y->Gamma_u, D/y->Gamma_cross, and the Fisher bill is 2 M3.
    gamma_u=tensor(u,lam(u))+tensor(lam(u),u)-lam(tensor(u,u))
    gamma_x=cross(u,lam(omega))+cross(lam(u),omega)-lam(lamb)
    trgamma=np.trace(gamma_u,axis1=0,axis2=1)
    assert np.linalg.norm(gamma_x-(0.5*grad_scalar(trgamma)-div_tensor(gamma_u))) <= 6e-12*max(1.0,np.linalg.norm(gamma_x))
    M3=inner(omega,lam(omega))

    # Full half-space Dirichlet and Hessian energies have constants K and 2 M3.
    modal=np.sum(np.abs(fft(u))**2,axis=0)/(n**6)
    K_fourier=float(np.sum(radius*modal))
    M3_fourier=float(np.sum((radius**3)*modal))
    assert K_fourier == pytest.approx(inner(u,lam(u)),rel=3e-13,abs=3e-13)
    assert 2.0*M3_fourier == pytest.approx(2.0*M3,rel=3e-13,abs=3e-13)

    boundary_ratio=np.zeros_like(trgamma); live=trgamma>1e-13
    boundary_ratio[live]=np.sum(gamma_x*gamma_x,axis=0)[live]/trgamma[live]
    boundary_fisher=float(np.mean(boundary_ratio))
    assert boundary_fisher <= 2.0*M3+2e-12

    FE=leray(lamb); kappa=inner(lam(u),FE)

    # The whole critical hierarchy collapses to one scalar Dirichlet-depth continuity law.
    Fy=py(FE,y)
    Lv=lam(v); L2v=lam(Lv)
    qdir=2.0*inner(Lv,Lv)
    qdir_y=-4.0*inner(Lv,L2v)
    qdir_yy=8.0*inner(L2v,L2v)
    jE=2.0*inner(Lv,Fy)
    jE_y=-4.0*inner(L2v,Fy)
    jnu=0.5*nu*qdir_y
    vt_full=Fy-nu*lap(v)
    qdir_t=4.0*inner(Lv,lam(vt_full))
    assert qdir_t+jE_y+0.5*nu*qdir_yy == pytest.approx(0.0,abs=6e-12)
    assert qdir_y <= 1e-13
    assert 2.0*inner(lam(u),FE) == pytest.approx(2.0*kappa,rel=3e-14,abs=3e-14)
    assert -2.0*nu*M3 == pytest.approx(0.5*nu*(-4.0*M3),rel=3e-14,abs=3e-14)
    assert inner(u,FE) == pytest.approx(0.0,abs=5e-13)  # int_0^infinity j_E dy

    # One completely monotone Dirichlet density encodes the full Sobolev jet/moment ladder.
    live_r=radius>0.0; rp=radius[live_r]; mp=modal[live_r]
    for order in range(5):
        jet=float(np.sum(2.0*((-2.0*rp)**order)*(rp**2)*mp))
        sob=(((-1.0)**order)*(2.0**(order+1))
             *float(np.sum((rp**(order+2))*mp)))
        assert jet == pytest.approx(sob,rel=3e-14,abs=3e-14)
        moment=float(np.sum(2.0*(rp**2)*math.factorial(order)/((2.0*rp)**(order+1))*mp))
        hnorm=(math.factorial(order)/(2.0**order))*float(np.sum((rp**(1-order))*mp))
        assert moment == pytest.approx(hnorm,rel=3e-14,abs=3e-14)
    assert float(np.sum(rp*mp)) == pytest.approx(inner(u,lam(u)),rel=3e-13,abs=3e-13)
    assert 0.5*float(np.sum(mp)) == pytest.approx(0.5*E,rel=3e-13,abs=3e-13)

    # The tail is exactly the same critical stock at Poisson depth y.
    Ktail=inner(v,Lv)
    Ktail_t=2.0*inner(vt_full,Lv)
    assert Ktail_t == pytest.approx(jE+0.5*nu*qdir_y,rel=4e-12,abs=4e-12)

    # Signed helicity is the difference channel of the same positive depth density.
    def jop(a):
        ah=fft(curl(a)); out=np.zeros_like(ah)
        np.divide(ah,radius,out=out,where=radius>0.0)
        return ifft(out)
    Jv=jop(v); JFy=jop(Fy)
    vp=0.5*(v+Jv); vm=0.5*(v-Jv)
    Fp=0.5*(Fy+JFy); Fm=0.5*(Fy-JFy)
    qH=2.0*inner(Lv,curl(v))
    qp=0.5*(qdir+qH); qm=0.5*(qdir-qH)
    assert qp >= -3e-13 and qm >= -3e-13
    assert qp == pytest.approx(2.0*inner(lam(vp),lam(vp)),rel=4e-12,abs=4e-12)
    assert qm == pytest.approx(2.0*inner(lam(vm),lam(vm)),rel=4e-12,abs=4e-12)
    jH=2.0*inner(curl(v),Fy)
    jp=0.5*(jE+jH); jm=0.5*(jE-jH)
    assert jp == pytest.approx(2.0*inner(lam(vp),Fp),rel=4e-12,abs=4e-12)
    assert jm == pytest.approx(2.0*inner(lam(vm),Fm),rel=4e-12,abs=4e-12)
    qHt_euler=4.0*inner(lam(Fy),curl(v))
    jH_y=-4.0*inner(lam(curl(v)),Fy)
    assert qHt_euler+jH_y == pytest.approx(0.0,abs=6e-12)
    jH0=2.0*inner(omega,FE)
    assert jH0 == pytest.approx(0.0,abs=6e-13)
    assert 0.5*(2.0*kappa+jH0) == pytest.approx(kappa,rel=4e-14,abs=4e-14)
    assert 0.5*(2.0*kappa-jH0) == pytest.approx(kappa,rel=4e-14,abs=4e-14)

    M0dir=inner(u,lam(u)); M1dir=0.5*E; q0dir=2.0*Z
    defect=M1dir*q0dir-M0dir*M0dir
    assert M0dir == pytest.approx(K_fourier,rel=3e-13,abs=3e-13)
    assert q0dir == pytest.approx(2.0*Z,rel=3e-14,abs=3e-14)
    assert defect == pytest.approx(E*Z-M0dir*M0dir,rel=3e-14,abs=3e-14)
    N2=Z/E
    escape=kappa*kappa/(N2*(E*Z-M0dir*M0dir))
    flux_escape=(2.0*kappa)**2*M1dir/(q0dir*defect)
    assert flux_escape == pytest.approx(escape,rel=3e-13,abs=3e-13)

    # Local parent rho_t+div_4 J=0; its integrated normal current is j_E+j_nu.
    vh=fft(v); vth=fft(vt_full)
    grad_v=[ifft(1j*kx*vh),ifft(1j*ky*vh),ifft(1j*kz*vh),-lam(v)]
    grad_vt=[ifft(1j*kx*vth),ifft(1j*ky*vth),ifft(1j*kz*vth),-lam(vt_full)]
    rho_t=2.0*sum(np.sum(a0*b0,axis=0) for a0,b0 in zip(grad_v,grad_vt))
    lap4v=ifft(-k2*vh)+lap(v)
    divJ=(-2.0*sum(np.sum(a0*b0,axis=0) for a0,b0 in zip(grad_vt,grad_v))
          -2.0*np.sum(vt_full*lap4v,axis=0))
    assert np.linalg.norm(rho_t+divJ) <= 7e-12*max(1.0,np.linalg.norm(rho_t))
    Jy=float(np.mean(np.sum(-2.0*vt_full*(-lam(v)),axis=0)))
    assert Jy == pytest.approx(jE+jnu,rel=4e-12,abs=4e-12)

    # Natural Euler transport action returns to the critical H^{-1/2} action.
    point_action=jE*jE/qdir if qdir>1e-30 else 0.0
    assert point_action <= 2.0*inner(Fy,Fy)+2e-12
    Fh=fft(FE); inv_radius=np.zeros_like(radius); live_radius=radius>1e-13; inv_radius[live_radius]=1.0/radius[live_radius]
    hminus_half=float(np.sum(np.sum(np.abs(Fh)**2,axis=0)*inv_radius)/(n**6))
    assert hminus_half >= point_action-2e-12

    # Modewise depth integration gives the exact moments and zero-integral Euler current.
    modal_u=np.sum(np.abs(fft(u))**2,axis=0)/(n**6)
    M0_modes=float(np.sum(radius*modal_u))
    M1_modes=0.5*float(np.sum(modal_u[live_radius]))
    assert M0_modes == pytest.approx(M0dir,rel=3e-13,abs=3e-13)
    assert M1_modes == pytest.approx(M1dir,rel=3e-13,abs=3e-13)
    current_depth_exact=float(np.mean(np.sum(u*FE,axis=0)))
    assert current_depth_exact == pytest.approx(0.0,abs=5e-13)

    boundary_work=inner(u,gamma_x)
    gamma_strain=tensor_inner(gamma_u,strain(u))
    assert boundary_work == pytest.approx(gamma_strain,rel=5e-12,abs=5e-12)
    assert boundary_work == pytest.approx(-2.0*kappa,rel=5e-12,abs=5e-12)

    # The activity slope is the exact abstract boundary pair-area quantity.
    M0=np.eye(3)[:,:,None,None,None]*trgamma-gamma_u
    activity_prime=float(np.mean(np.einsum("i...,ij...,j...->...",u,M0,u)))
    assert kappa*kappa <= 0.5*activity_prime*M3+2e-12
    eps=1.0e-4
    ve=py(u,eps); te=py(tensor(u,u),eps)-tensor(ve,ve)
    tre=np.trace(te,axis1=0,axis2=1)
    Me=np.eye(3)[:,:,None,None,None]*tre-te
    Ae=float(np.mean(np.einsum("i...,ij...,j...->...",ve,Me,ve)))
    assert Ae/eps == pytest.approx(activity_prime,rel=8e-4,abs=8e-4)

def test_endogenous_euler_coefficient_pair_projection_and_static_cancellation_are_exact():
    # Scalar audit of A_u=nu_E Domega+C_perp in the canonical pair Hilbert metric.
    m3=3.7; nu_e=0.41; cperp2=1.9; nu=0.23
    kappa=nu_e*m3
    a2=nu_e*nu_e*m3+cperp2
    out=endogenous_euler_pair_projection(a2,m3,kappa,nu)
    assert out["endogenous_euler_coefficient"] == pytest.approx(nu_e,rel=1e-15)
    assert out["signed_productive_reynolds"] == pytest.approx(nu_e/nu,rel=1e-15)
    assert out["orthogonal_reconfiguration_squared"] == pytest.approx(cperp2,rel=2e-15)
    assert out["critical_rate"] == pytest.approx(2.0*m3*(nu_e-nu),rel=2e-15)
    # Changing only the orthogonal 3D road cannot change the instantaneous K rate.
    out2=endogenous_euler_pair_projection(a2+8.0,m3,kappa,nu)
    assert out2["critical_rate"] == pytest.approx(out["critical_rate"],rel=1e-15)
    assert out2["orthogonal_reconfiguration_squared"] == pytest.approx(cperp2+8.0,rel=2e-15)

    # Triangle coboundary identity behind equality rigidity.
    import numpy as np
    ux=np.array([0.7,-0.4,1.2]); uy=np.array([-0.2,0.9,0.3]); uz=np.array([0.5,0.1,-0.6])
    lhs=np.cross(ux,uy)+np.cross(uy,uz)+np.cross(uz,ux)
    rhs=np.cross(ux-uz,uy-uz)
    assert np.max(np.abs(lhs-rhs)) <= 4e-16


def test_lossless_energy_sphere_projection_collapses_productive_action_exactly():
    import numpy as np
    rng=np.random.default_rng(2026081556)
    for n in (4,7,11):
        ell=np.exp(rng.uniform(-1.2,1.4,n))
        X=rng.normal(size=n)
        Y0=rng.normal(size=n)
        Y=Y0-X*float(X@Y0)/float(X@X)  # exact Euler tangent direction
        LX=ell*X
        E=64.0*float(X@X)
        K=64.0*float(X@LX)
        Z=64.0*float(LX@LX)
        m=K/E
        g=(ell-m)*X
        kappa=64.0*float(g@Y)
        out=critical_energy_sphere_escape_geometry(E,K,Z,kappa)
        assert abs(float(X@Y)) <= 3e-15*max(1.0,float(np.linalg.norm(X)*np.linalg.norm(Y)))
        assert out["x_norm_squared"] == pytest.approx(float(X@X),rel=3e-14)
        assert out["x_Lx_pairing"] == pytest.approx(float(X@LX),rel=3e-14)
        assert out["Lx_norm_squared"] == pytest.approx(float(LX@LX),rel=3e-14)
        assert out["uphill_tangent_norm_squared"] == pytest.approx(float(g@g),rel=3e-13,abs=3e-13)
        assert out["uphill_euler_pairing"] == pytest.approx(float(g@Y),rel=3e-13,abs=3e-13)
        direct=(float((g/np.linalg.norm(g))@Y)**2)/float(LX@LX)
        assert out["escape_action"] == pytest.approx(direct,rel=5e-13,abs=5e-13)


def test_physical_scalar_triple_product_is_exact_global_L2_escape_angle():
    import numpy as np
    rng=np.random.default_rng(2026081557)
    u=rng.normal(size=(19,3))
    Lambda_u=rng.normal(size=(19,3))
    omega=rng.normal(size=(19,3))
    m=0.73
    s=Lambda_u-m*u
    kappa=float(np.sum(s*np.cross(u,omega)))
    cyclic=float(np.sum(omega*np.cross(s,u)))
    assert kappa == pytest.approx(cyclic,rel=2e-15,abs=2e-15)
    assert np.max(np.abs(np.cross(s,u)-np.cross(Lambda_u,u))) <= 8e-16
    shat=s/np.linalg.norm(s)
    omegahat=omega/np.linalg.norm(omega)
    action=float(np.sum(omegahat*np.cross(shat,u)))**2
    expected=kappa*kappa/(float(np.sum(omega*omega))*float(np.sum(s*s)))
    assert action == pytest.approx(expected,rel=3e-15,abs=3e-15)


def test_spacetime_critical_hom_bianchi_electric_split_and_transgression_are_exact():
    import numpy as np

    def wedge_basis(j):
        W=np.zeros((16,16),complex)
        bit=1<<j
        for mask in range(16):
            if mask & bit:
                continue
            swaps=sum(1 for k in range(j) if mask & (1<<k))
            W[mask|bit,mask]=(-1.0)**swaps
        return W

    W=[wedge_basis(j) for j in range(4)]  # dt, dx, dy, dz
    n=7
    x=2*math.pi*np.arange(n)/n
    kk=np.fft.fftfreq(n,1.0/n)
    Fm=np.exp(2j*math.pi*np.outer(np.arange(n),np.arange(n))/n)/math.sqrt(n)
    D=Fm@np.diag(1j*kk)@Fm.conj().T
    rr=np.zeros(n,float); nz=np.abs(kk)>1e-12; rr[nz]=1.0/np.abs(kk[nz])
    Rsp=Fm@np.diag(rr)@Fm.conj().T
    I16=np.eye(16); In=np.eye(n)
    d=np.kron(D,W[1])
    R=np.kron(Rsp,I16)
    Wt=np.kron(In,W[0])

    alpha=np.stack([np.sin(x)+0.11*np.cos(2*x),0.23*np.sin(2*x),0.19*np.cos(x)],axis=1)
    Ealpha=np.zeros((16*n,16*n),complex)
    for i,a in enumerate(alpha):
        Ealpha[16*i:16*(i+1),16*i:16*(i+1)]=a[0]*W[1]+a[1]*W[2]+a[2]*W[3]
    Bvals=0.37*np.cos(x)+0.13*np.sin(2*x)
    MB=np.kron(np.diag(Bvals),I16)
    E4=Ealpha-MB@Wt  # mathbb A wedge = alpha wedge - B dt wedge
    F4=d@E4+E4@d
    A4=R@E4-E4@R
    nu=0.31
    coeff=spacetime_critical_hom_coefficients(nu)
    G4=coeff["curvature_commutator_coefficient"]*(R@F4-F4@R)+E4@A4
    DL=nu*d+0.5*E4; DR=nu*d
    nabla_A4=DL@A4+A4@DR
    nabla_G4=DL@G4-G4@DR
    scale=max(1.0,float(np.linalg.norm(G4)),float(np.linalg.norm(nu*F4@A4)))
    assert np.linalg.norm(2.0*nabla_A4-G4) <= 5e-11*scale
    assert np.linalg.norm(nabla_G4-nu*F4@A4) <= 5e-11*scale

    # Freeze alpha in time.  Then e=-(d B), enough to referee every electric sign.
    Fsp=d@Ealpha+Ealpha@d
    Asp=R@Ealpha-Ealpha@R
    Gsp=2.0*nu*(R@Fsp-Fsp@R)+Ealpha@Asp
    Gtemp=G4-Gsp
    CRB=R@MB-MB@R
    eop=-(d@MB-MB@d)
    predicted=Wt@(
        -2.0*nu*(R@eop-eop@R)
        +Ealpha@CRB
        -MB@Asp
    )
    assert np.linalg.norm(Gtemp-predicted) <= 6e-11*max(1.0,float(np.linalg.norm(Gtemp)))

    # Graded transgression: nabla(A4 G4)=(nabla A4)G4-A4(nabla G4).
    trans=nabla_A4@G4-A4@nabla_G4
    rhs=0.5*(G4@G4)-nu*A4@F4@A4
    assert np.linalg.norm(trans-rhs) <= 8e-11*max(1.0,float(np.linalg.norm(rhs)))
    assert coeff["electric_e_commutator_coefficient"] == pytest.approx(-2.0*nu)
    assert coeff["transgression_square_coefficient"] == pytest.approx(0.5)


def test_affine_core_budget_falsifier_exponents_leave_only_dynamic_compatibility():
    # This is a scaling countergeometry, explicitly not an NS solution.
    for alpha in (0.4,0.43,0.49):
        eE=-2.0+5.0*alpha
        eK=-2.0+4.0*alpha
        eZ=-2.0+3.0*alpha
        eM3=-2.0+2.0*alpha
        eKap=-3.0+4.0*alpha
        ePairArea=-4.0+6.0*alpha
        eRatio=eKap-eM3
        assert eE >= -2e-15          # bounded core energy
        assert eZ > -1.0             # int Z dt finite
        assert eK < 0.0              # K diverges
        assert eM3 <= -1.0           # int M3 dt diverges
        assert ePairArea <= -1.0      # int ||A_u||_pair^2 dt diverges
        assert eRatio < 0.0           # kappa/(nu M3) diverges


def test_primitive_pair_critical_scaling_is_neutral_only_at_K():
    # Exact 3D NS dilation u_lambda=lambda u(lambda x,lambda^2 t).
    E,Z,K,M3,kappa,A2,N2,dt=2.3,1.7,0.91,4.2,-1.4,3.6,1.7/2.3,0.08
    det=E*Z-K*K
    for lam in (0.17,0.8,2.5,11.0):
        El=E/lam; Zl=lam*Z; Kl=K; M3l=lam*lam*M3; kapl=lam*lam*kappa; A2l=lam*lam*A2
        N2l=lam*lam*N2; dtl=dt/(lam*lam)
        assert Kl == K
        assert El*Zl-Kl*Kl == pytest.approx(det,rel=2e-15)
        assert M3l*dtl == pytest.approx(M3*dt,rel=2e-15)
        assert A2l*dtl == pytest.approx(A2*dt,rel=2e-15)
        action_l=kapl*kapl/(N2l*(El*Zl-Kl*Kl))*dtl
        action=kappa*kappa/(N2*det)*dt
        assert action_l == pytest.approx(action,rel=3e-14,abs=3e-14)


def test_certificate_records_two_particle_critical_history_without_closure_claim():
    cert=theorem_certificate()
    assert "pressure has only the common-coordinate road" in cert["primitive_two_particle_transport"]
    assert "Euler only exchanging their energies" in cert["primitive_common_relative_pair_law"]
    assert "Delta_c v=4 Delta_r v" in cert["primitive_pair_endpoint_compatibility"]
    assert "center-Dirichlet cost of inverted-pair velocity" in cert["primitive_pair_critical_field"]
    assert "nu_E=kappa/M3" in cert["endogenous_euler_coefficient"]
    assert "X_t=Y-nu L^2 X" in cert["primitive_energy_sphere_dynamics"]
    assert "tangent Rayleigh-gradient direction" in cert["primitive_critical_uphill_projection"]
    assert "A_escape=kappa^2/[N^2(EZ-K^2)]" in cert["primitive_scalar_triple_escape"]
    assert "only the first term has fixed sign" in cert["primitive_euler_kappa_derivative"]
    assert "remaining freedom is tangent turning" in cert["primitive_euler_acceleration_leaf"]
    assert "P_t=[A_e,P]-nu[P,[P,Lambda^2]]" in cert["primitive_normalized_projector_two_road"]
    assert "K'_E=2(omega,Sym_G U omega)_G=2kappa" in cert["primitive_krein_lax_boost"]
    assert "[V_r,V_(r^2)]=0" in cert["primitive_radial_tilt_commutation"]
    assert "D_t w=" in cert["primitive_orientation_covariant_law"]
    assert "self-rotating isospectral curl frame" in cert["primitive_corotating_heat_frame"]
    assert "scale growth, spread and turning" in cert["primitive_frame_commutator_current"]
    assert "A_escape=a^2/N^2" in cert["primitive_absolute_curl_productive_frame"]
    assert "self-turning" in cert["primitive_absolute_curl_self_turning"]
    assert "no 0/0 quotient" in cert["primitive_helicity_side_motion"]
    assert "first Toda spectral-measure flow" in cert["primitive_toda_saturation_flow"]
    assert "at most two radii" in cert["primitive_two_radius_krylov_rigidity"]
    assert "extra |C| weight" in cert["primitive_turning_reconfiguration_balance"]
    assert "all four sign quadrants" in cert["primitive_one_step_turning_guard"]
    assert "(lambda_M-lambda_I)" in cert["primitive_full_two_step_gap"]
    assert "bracket-level A+B+C=0" in cert["primitive_weighted_jacobi_continuation"]
    assert "both cyclic full companions vanish" in cert["primitive_full_continuation_companion_guard"]
    assert "max |m_r|>=2|ell|/3" in cert["primitive_fourier_diamond_carrier_guard"]
    assert "representation-free positive composition bridge is unproved" in cert["primitive_jacobi_biot_savart_seam"]
    assert "retained child/donor work tends 1" in cert["uv_curvature"]
    assert "10/13" in cert["primitive_fixed_loss_window_guard"]
    assert "gross-traffic budget" in cert["primitive_positive_curvature_composition_guard"]
    assert "whole critical hierarchy is one boundary jet" in cert["primitive_poisson_energy_profile"]
    assert "P_y(fg)-P_yf P_yg" in cert["primitive_poisson_product_defect"]
    assert "R'_E(0)=-4kappa" in cert["primitive_poisson_boundary_critical_law"]
    assert "zero-depth-integral return current" in cert["primitive_poisson_depth_two_road"]
    assert "M1 q(0)-M0^2=EZ-K^2" in cert["primitive_poisson_hankel_escape"]
    assert "fake scalar profile" in cert["primitive_poisson_realizability_guard"]
    assert "positive covariance-stress lift" in cert["primitive_poisson_covariance_stress"]
    assert "Delta4 v=0" in cert["primitive_harmonic_halfspace_normal_form"]
    assert "-Delta4 tau=2 sum_A" in cert["primitive_harmonic_covariance_source"]
    assert "int_R4+ |grad_4 v|^2=K" in cert["primitive_harmonic_critical_energy"]
    assert "dynamical no-concentration theorem remains unproved" in cert["primitive_harmonic_concentration_guard"]
    assert "tau_(y+z)(u)=P_z tau_y(u)+tau_z(P_yu)" in cert["primitive_poisson_germano_cocycle"]
    assert "(1/2) partial_y V|_0=K" in cert["primitive_poisson_two_reservoir_exchange"]
    assert "stress-energy tensor" in cert["primitive_poisson_stress_divergence"]
    assert "constant 1" in cert["primitive_poisson_covariance_heat_fisher"]
    assert "<=2M3" in cert["primitive_poisson_covariance_boundary"]
    assert "A_(y+z)=A_z(P_yu)" in cert["primitive_poisson_transverse_activity"]
    assert "rho_t+div_4 J=0" in cert["primitive_poisson_dirichlet_local_current"]
    assert "M0=int mathfrak q dy=K" in cert["primitive_poisson_dirichlet_moments"]
    assert "(-1)^n partial_y^n q(0)" in cert["primitive_poisson_dirichlet_sobolev_ladder"]
    assert "K_y=int_y^infinity mathfrak q" in cert["primitive_poisson_dirichlet_tail_law"]
    assert "j_+^E(0)=j_-^E(0)=kappa" in cert["primitive_poisson_helicity_depth_split"]
    assert "int j_E(y)^2/mathfrak q(y) dy" in cert["primitive_poisson_transport_action_guard"]
    assert "pi common-phase reversal" in cert["primitive_poisson_orientation_incompleteness"]
    assert "transport/action ratio 3.393252e-6" in cert["primitive_poisson_transport_lower_guard"]
    assert "R_E^2/(4A_y)+R_align+R_reg" in cert["primitive_poisson_covariance_pythagoras"]
    assert "V_t=A_y/(2nu)" in cert["primitive_poisson_relative_completed_square"]
    assert "rank changes of the pseudoinverse" in cert["primitive_poisson_axis_gradient_guard"]
    assert "int_0^infinity A_y dy" in cert["primitive_poisson_activity_depth_area_guard"]
    assert "G4=2 nabla_4 A4" in cert["primitive_spacetime_critical_bianchi"]
    assert "Ec=-2nu[R,e wedge]" in cert["primitive_spacetime_critical_electric"]
    assert "not supply int||Ec||^2 dt" in cert["primitive_spacetime_transgression_guard"]
    assert "argmin_lambda" in cert["primitive_productive_pair_projection"]
    assert "affine line" in cert["primitive_productive_alignment_rigidity"]
    assert "J_lambda=V+lambda B" in cert["primitive_hom_current_lambda_family"]
    assert "cancels exactly" in cert["primitive_static_reconfiguration_cancellation"]
    assert "K_parallel=K/4" in cert["primitive_pair_radial_transverse_split"]
    assert "K_perp=3K/4" in cert["primitive_pair_radial_transverse_split"]
    assert "P_parallel=3(A3-B3)/(2pi^2)" in cert["primitive_pressure_pair_converter"]
    assert "mints no total critical energy" in cert["primitive_pressure_pair_converter"]
    assert "Q_t+div_r J=2nu Delta_r Q-2nu Z" in cert["primitive_relative_pressure_free_law"]
    assert "weighted inward relative transport toward r=0" in cert["primitive_relative_pressure_free_law"]
    assert "K=(2pi^2)^-1" in cert["primitive_inverted_pair_kinetic"]
    assert "finite Hilbert path length" in cert["primitive_inverted_pair_history"]
    assert "not endpoint-velocity control or regularity" in cert["primitive_inverted_pair_history"]
    assert "homogeneous in C" in cert["primitive_pair_affine_defect"]
    assert "8E/(pi R)" in cert["primitive_pair_collision_concentration"]
    assert "exactly the enstrophy boundary law" in cert["primitive_pair_collision_boundary"]
    assert "K->K" in cert["primitive_pair_scale_neutrality"]
    assert "regeneration/accumulation of critical mass" in cert["primitive_pair_scale_neutrality"]
    assert "pair critical stock, physical energy-loss speed and critical heat" in cert["primitive_pair_material_scale_lock"]
    assert "int kappa^2/[N^2(EZ-K^2)] dt=infinity" in cert["primitive_regeneration_persistence"]
    assert "A(F_E)=-2 d_op^* V" in cert["primitive_regeneration_shortcut_guard"]
    assert "successfully regenerated/accumulated infinitely often" in cert["primitive_two_road_frontier"]
    assert "unproved" in cert["primitive_two_road_frontier"]
    assert cert["global_regularity_claimed"] is False

def test_primitive_critical_gauss_bianchi_continuum_constants_are_exact():
    out=continuum_critical_gauss_bianchi_constants()
    assert out["riesz_kernel_coefficient"] == pytest.approx(1.0/(2.0*math.pi**2),rel=1e-15)
    assert out["pair_area_hs_raw_coefficient"] == pytest.approx(1.0/(2.0*math.pi**4),rel=1e-15)
    assert out["gauss_source_hs_raw_coefficient"] == pytest.approx(1.0/(4.0*math.pi**4),rel=1e-15)
    assert out["operator_to_pair_norm_squared_factor"] == pytest.approx(1.0/math.pi**2,rel=1e-15)
    assert out["critical_square_hs_prefactor_times_viscosity"] == pytest.approx(math.pi**2/2.0,rel=1e-15)
    assert out["top_channel_normal_potential_coefficient"] == pytest.approx(0.25,rel=1e-15)


def test_primitive_critical_gauss_bianchi_is_exact_hom_connection_algebra():
    import numpy as np

    # Full exterior algebra on R^3, represented by bit masks 0,...,7.
    def wedge_basis(j):
        W=np.zeros((8,8),complex)
        bit=1<<j
        for mask in range(8):
            if mask & bit:
                continue
            swaps=sum(1 for k in range(j) if mask & (1<<k))
            W[mask|bit,mask]=(-1.0)**swaps
        return W

    W=[wedge_basis(j) for j in range(3)]
    n=7
    x=2*math.pi*np.arange(n)/n
    kk=np.fft.fftfreq(n,1.0/n)
    Fm=np.exp(2j*math.pi*np.outer(np.arange(n),np.arange(n))/n)/math.sqrt(n)
    D=Fm@np.diag(1j*kk)@Fm.conj().T
    rr=np.zeros(n,float)
    nz=np.abs(kk)>1e-12
    rr[nz]=1.0/np.abs(kk[nz])
    Rsp=Fm@np.diag(rr)@Fm.conj().T
    I8=np.eye(8)
    d=np.kron(D,W[0])
    R=np.kron(Rsp,I8)

    alpha=np.stack([np.zeros(n),np.sin(x)+0.17*np.cos(2*x),0.31*np.sin(2*x)],axis=1)
    E=np.zeros((8*n,8*n),complex)
    for i,a in enumerate(alpha):
        block=sum(a[j]*W[j] for j in range(3))
        E[8*i:8*(i+1),8*i:8*(i+1)]=block
    F=d@E+E@d
    A=R@E-E@R
    B=R@F-F@R
    V=E@R@E
    nu=0.37
    D0=nu*d
    Dh=nu*d+0.5*E
    G=V+2.0*nu*B

    scale=max(1.0,np.linalg.norm(G),np.linalg.norm(nu*F@A))
    assert np.linalg.norm(d@A+A@d-B) <= 3e-11*scale
    assert np.linalg.norm(A@E+E@A) <= 3e-11*scale
    assert np.linalg.norm(E@A-V) <= 3e-11*scale
    nabla_A=Dh@A+A@D0  # degree(A)=1
    assert np.linalg.norm(2.0*nabla_A-G) <= 4e-11*scale
    nabla_G=Dh@G-G@D0  # degree(G)=2
    assert np.linalg.norm(nabla_G-nu*F@A) <= 5e-11*scale


def test_critical_gauss_null_source_contains_2d_and_shear_but_not_generic_3d_geometry():
    import numpy as np

    # Embedded 2D example: omega is vertical while all velocity differences are horizontal.
    pts=[(0.2,0.4,0.7),(1.1,0.6,2.0),(2.2,1.4,0.3)]
    def u2(p):
        x,y,z=p
        return np.array([math.sin(y),math.sin(x),0.0])
    def w2(p):
        x,y,z=p
        return np.array([0.0,0.0,math.cos(x)-math.cos(y)])
    for p0 in pts:
        for p1 in pts:
            assert abs(float(w2(p0)@(u2(p1)-u2(p0)))) <= 2e-14

    # One-direction shear: u=phi(y)e1, (u.grad)u=0 and the same source vanishes.
    def us(p):
        x,y,z=p
        return np.array([math.sin(y),0.0,0.0])
    def ws(p):
        x,y,z=p
        return np.array([0.0,0.0,-math.cos(y)])
    for p0 in pts:
        for p1 in pts:
            assert abs(float(ws(p0)@(us(p1)-us(p0)))) <= 2e-14

    # ABC flow is genuinely 3D and Beltrami; one generic pair has nonzero source.
    def u3(p):
        x,y,z=p
        return np.array([math.sin(z)+math.cos(y),math.sin(x)+math.cos(z),math.sin(y)+math.cos(x)])
    p0,p1=pts[0],pts[1]
    source=float(u3(p0)@(u3(p1)-u3(p0)))  # omega=u for A=B=C=1 ABC
    assert abs(source) > 1e-3


def test_certificate_records_critical_gauss_bianchi_residual_without_stability_overclaim():
    cert=theorem_certificate()
    assert "kappa(0)=-pi^2" in cert["primitive_critical_pair_area_loop"]
    assert "Gc=V+2nu B" in cert["primitive_critical_gauss_residual"]
    assert "nabla Gc=nu(beta wedge)A" in cert["primitive_critical_gauss_bianchi"]
    assert "embedded 2D" in cert["primitive_critical_gauss_null_geometry"]
    assert "L6=-nu^2(Delta_x+Delta_y)+|u(x)|^2/4" in cert["primitive_critical_six_dimensional_gauss_floor"]
    assert "separate stability/compactness theorem" in cert["primitive_critical_gauss_stability_guard"]
    assert "quantitative stability/history" in cert["primitive_turning_frontier"]
    assert cert["global_regularity_claimed"] is False



def test_critical_six_dimensional_gauss_normal_form_has_no_hidden_drift_cross_term():
    import numpy as np

    # Smooth periodic divergence-free ABC velocity.  A separable low-mode top-form
    # test phi(x,y)=f(x)g(y) makes the 6D quadratic form inexpensive to referee.
    n=10
    grid=2*math.pi*np.arange(n)/n
    X=np.stack(np.meshgrid(grid,grid,grid,indexing="ij"),axis=-1)
    x,y,z=X[...,0],X[...,1],X[...,2]
    u=np.stack([np.sin(z)+np.cos(y),np.sin(x)+np.cos(z),np.sin(y)+np.cos(x)],axis=-1)
    f=1.0+0.21*np.cos(x)+0.17*np.sin(y)+0.13*np.cos(z)
    gradf=np.stack([-0.21*np.sin(x),0.17*np.cos(y),-0.13*np.sin(z)],axis=-1)
    g=1.0+0.19*np.sin(x)-0.11*np.cos(y)+0.07*np.sin(z)
    gradg=np.stack([0.19*np.cos(x),0.11*np.sin(y),0.07*np.cos(z)],axis=-1)
    nu=0.43

    # Product averages factor in y.  The output top->2 channel is Hodge-equivalent
    # to nu grad_x(phi)+(1/2)u phi; the input channel is nu grad_y(phi), and the
    # two Hom channels are orthogonal.
    mean_g2=float(np.mean(g*g))
    mean_gradg2=float(np.mean(np.sum(gradg*gradg,axis=-1)))
    out_vec=nu*gradf+0.5*u*f[...,None]
    lhs=mean_g2*float(np.mean(np.sum(out_vec*out_vec,axis=-1))) \
        +nu*nu*float(np.mean(f*f))*mean_gradg2
    rhs=nu*nu*(mean_g2*float(np.mean(np.sum(gradf*gradf,axis=-1)))
               +float(np.mean(f*f))*mean_gradg2) \
        +0.25*mean_g2*float(np.mean(np.sum(u*u,axis=-1)*(f*f)))
    cross=float(np.mean(np.sum(u*gradf,axis=-1)*f))
    assert abs(cross) <= 2e-15
    assert lhs == pytest.approx(rhs,rel=3e-14,abs=3e-14)


def test_critical_increment_stress_work_angle_and_nonlocal_transport_referee():
    import numpy as np

    n=16
    rng=np.random.default_rng(7341)
    kk=np.fft.fftfreq(n,1.0/n)
    kx,ky,kz=np.meshgrid(kk,kk,kk,indexing="ij")
    ks=(kx,ky,kz)
    k2=kx*kx+ky*ky+kz*kz
    kmag=np.sqrt(k2)
    kvec=np.stack(ks)

    def fft(a):
        return np.fft.fftn(a,axes=tuple(range(a.ndim-3,a.ndim)))

    def ifft(a):
        return np.fft.ifftn(a,axes=tuple(range(a.ndim-3,a.ndim))).real

    def leray(vh):
        dot=np.sum(kvec*vh,axis=0)
        return vh-kvec*dot[None]/np.where(k2==0.0,1.0,k2)[None]

    def curlh(vh):
        return 1j*np.stack([
            ky*vh[2]-kz*vh[1],
            kz*vh[0]-kx*vh[2],
            kx*vh[1]-ky*vh[0],
        ])

    def lambda_h(h):
        return kmag*h if h.ndim==3 else kmag[None]*h

    def cross(a,b):
        return np.cross(a,b,axisa=0,axisb=0,axisc=0)

    # Low support <=2 leaves room for every product used below, so this is an
    # exact finite Fourier referee rather than a wrap-around aliasing artifact.
    uh=fft(rng.normal(size=(3,n,n,n)))
    uh*=((kmag<=2.0)&(kmag>0.0))
    uh=leray(uh)
    uh[:,k2==0.0]=0.0
    u=ifft(uh)
    uh=leray(fft(u))
    omega_h=curlh(uh)
    omega=ifft(omega_h)
    lambda_u_h=lambda_h(uh)
    lambda_u=ifft(lambda_u_h)
    lambda_omega=ifft(lambda_h(omega_h))
    fe_h=leray(fft(cross(u,omega)))
    fe=ifft(fe_h)
    lambda_fe_h=lambda_h(fe_h)

    def gamma_tensor(a,b):
        la=ifft(lambda_h(fft(a)))
        lb=ifft(lambda_h(fft(b)))
        out=np.empty((3,3,n,n,n))
        for i in range(3):
            for j in range(3):
                out[i,j]=(a[i]*lb[j]+b[j]*la[i]
                          -ifft(lambda_h(fft(a[i]*b[j]))))
        return out

    gamma=gamma_tensor(u,u)
    div_gamma_h=np.zeros_like(uh)
    for i in range(3):
        for j in range(3):
            div_gamma_h[i]+=1j*ks[j]*fft(gamma[i,j])
    p_div_gamma_h=leray(div_gamma_h)
    u_cross_lambda_omega_h=leray(fft(cross(u,lambda_omega)))
    lambda_u_cross_omega_h=leray(fft(cross(lambda_u,omega)))
    corrected=(lambda_fe_h-u_cross_lambda_omega_h-lambda_u_cross_omega_h)

    def rel(a,b):
        return float(np.linalg.norm(a-b)/max(1.0e-30,np.linalg.norm(a),np.linalg.norm(b)))

    assert rel(p_div_gamma_h,corrected) < 2.0e-12
    # Two tempting simplifications are genuinely false on this actual state.
    assert rel(p_div_gamma_h,lambda_fe_h) > 1.0e-2
    assert rel(lambda_fe_h,u_cross_lambda_omega_h-lambda_u_cross_omega_h) > 1.0e-2

    grad_u=np.empty((3,3,n,n,n))
    for i in range(3):
        for j in range(3):
            grad_u[i,j]=ifft(1j*ks[j]*uh[i])
    strain=0.5*(grad_u+np.swapaxes(grad_u,0,1))
    div_gamma=ifft(div_gamma_h)
    kappa=float(np.mean(np.sum(lambda_u*fe,axis=0)))
    u_div_gamma=float(np.mean(np.sum(u*div_gamma,axis=0)))
    gamma_strain=float(np.mean(np.sum(gamma*strain,axis=(0,1))))
    scale=max(1.0,abs(kappa),abs(u_div_gamma),abs(gamma_strain))
    assert abs(u_div_gamma-2.0*kappa) < 3.0e-12*scale
    assert abs(gamma_strain+2.0*kappa) < 3.0e-12*scale

    # Exact Euler material derivative of Gamma from the operator law.  The
    # local upper-convected candidate has an order-one residual.
    gamma_t=gamma_tensor(fe,u)+gamma_tensor(u,fe)
    adv_gamma=np.empty_like(gamma)
    for i in range(3):
        for j in range(3):
            gh=fft(gamma[i,j])
            adv_gamma[i,j]=sum(u[a]*ifft(1j*ks[a]*gh) for a in range(3))
    material_gamma=gamma_t+adv_gamma
    gamma_grad=np.einsum('ik...,kj...->ij...',gamma,grad_u)
    grad_t_gamma=np.einsum('ki...,kj...->ij...',grad_u,gamma)
    upper_residual=material_gamma+gamma_grad+grad_t_gamma
    rms=lambda a: float(np.sqrt(np.mean(np.abs(a)**2)))
    assert rms(upper_residual) > 0.5*rms(material_gamma)

    # The single angle cos(theta)=K/sqrt(EZ) also does not telescope: its
    # exact Euler derivative contains both kappa and vortex-stretch work P.
    E=float(np.mean(np.sum(u*u,axis=0)))
    K=float(np.mean(np.sum(u*lambda_u,axis=0)))
    Z=float(np.mean(np.sum(omega*omega,axis=0)))
    stretch=float(np.mean(np.sum(omega*np.einsum('ij...,j...->i...',strain,omega),axis=0)))
    cosine=K/math.sqrt(E*Z)
    assert 0.0 < cosine < 1.0
    theta=math.acos(cosine)
    rhs=2.0*kappa/K-stretch/Z
    theta_dot=-rhs/math.tan(theta)
    assert theta_dot/(-1.0/math.tan(theta)) == pytest.approx(rhs,rel=2e-14,abs=2e-14)


def test_certificate_records_increment_stress_and_rejects_fake_persistence_shortcuts():
    cert=theorem_certificate()
    assert "P div Gamma_u=Lambda F_E-P(u cross Lambda omega)-P((Lambda u) cross omega)" in cert["primitive_critical_increment_stress_work"]
    assert "kappa=-(1/2) int Gamma_u:S" in cert["primitive_critical_increment_stress_work"]
    assert "not a local conformation tensor" in cert["primitive_critical_increment_transport_guard"]
    assert "det Gamma is not a material invariant" in cert["primitive_critical_increment_transport_guard"]
    assert "-tan(theta) theta'_E=2kappa/K-P/Z" in cert["primitive_single_angle_guard"]
    assert "one-angle monotone or telescope" in cert["primitive_single_angle_guard"]
    assert cert["global_regularity_claimed"] is False
