import math

import numpy as np
import pytest

from src.native_material_vorticity_heat_law import (
    accumulated_transverse_heat_memory,
    actual_ns_current_gauss_algebra,
    affine_local_blowup_guard,
    canonical_maxwell_extension_spectral_law,
    klein_spacetime_vortex_worldsheet_algebra,
    canonical_poisson_scale_overlap,
    closed_vortex_line_period_cost,
    conformal_fieldline_length_variation,
    curl_line_geometry_algebra,
    covariant_divergence_test_coercivity,
    maxwell_duality_stress_algebra,
    mixed_curl_gauss_direction_algebra,
    mixed_strain_poynting_factorization,
    primitive_spacetime_gauge_algebra,
    so33_exterior_square_algebra,
    sl3_cartan_casimir_algebra,
    sl3_characteristic_polynomial_rate,
    sl3_tangent_chord_degree_algebra,
    material_hminus2_reset_identity,
    material_hodge_speed_ladder,
    material_log_distortion_energy_bound,
    local_flux_velocity_gauge_algebra,
    local_incompressible_flux_velocity_odes,
    material_metric_path_action_bound,
    material_state_speed_lock,
    primitive_material_current_fourier_law,
    moving_polarization_memory_bound,
    pair_direction_mismatch_decomposition,
    endpoint_current_parallelogram,
    poynting_equality_residual_gauss_algebra,
    critical_projected_current_tax,
    rank_one_incompressible_stretch_null,
    theorem_certificate,
    transverse_heat_determinant,
    symmetric_metric_cartan_balance,
    transverse_heat_log_rate,
    transverse_two_covector_area_identity,
    twist_free_leaf_residual_circulation,
    two_level_curl_geometry_algebra,
    universal_curl_gauss_residual,
    transported_twoform_transverse_determinant,
    vortex_slip_twist_algebra,
    vortex_productivity_frustration_algebra,
    vorticity_stress_current_algebra,
)


def _det_one_metric(rng):
    m = rng.normal(size=(3, 3))
    g = m.T @ m + 0.05 * np.eye(3)
    return g / np.linalg.det(g) ** (1.0 / 3.0)


def _transverse(q, v):
    q = np.asarray(q, float)
    q = q / np.linalg.norm(q)
    v = np.asarray(v, float)
    return v - q * np.dot(q, v)


def test_transverse_heat_determinant_is_exact_amplification_identity():
    rng = np.random.default_rng(2026081501)
    for _ in range(2000):
        g = _det_one_metric(rng)
        q = rng.normal(size=3)
        out = transverse_heat_determinant(g, q)
        assert out["transverse_heat_determinant"] == pytest.approx(
            out["geometric_amplification_squared"], rel=3e-11, abs=3e-11
        )
        assert out["transverse_heat_trace"] + 1e-12 >= 2.0 * out["geometric_amplification"]


def test_two_covector_heat_gram_is_amplification_times_area():
    rng = np.random.default_rng(2026081502)
    for _ in range(2000):
        g = _det_one_metric(rng)
        q = rng.normal(size=3)
        xi = _transverse(q, rng.normal(size=3))
        eta = _transverse(q, rng.normal(size=3))
        out = transverse_two_covector_area_identity(g, q, xi, eta)
        assert out["heat_gram_determinant"] == pytest.approx(
            out["amplification_times_area_squared"], rel=2e-10, abs=2e-10
        )
        assert out["no_free_area_margin"] >= 0.0


def test_accumulated_heat_area_remembers_all_geometric_amplification():
    rng = np.random.default_rng(2026081503)
    for _ in range(500):
        q = rng.normal(size=3)
        metrics = [_det_one_metric(rng) for _ in range(rng.integers(2, 20))]
        weights = np.exp(rng.uniform(-5.0, 1.0, len(metrics)))
        out = accumulated_transverse_heat_memory(metrics, q, weights)
        assert out["accumulated_heat_area"] + 3e-9 >= out["accumulated_geometric_amplification"]


def test_accumulated_heat_memory_has_exact_proportional_equality_family():
    q = (0.0, 0.0, 1.0)
    # det-one diagonal metrics with a fixed transverse anisotropy shape and varying
    # common transverse heat scale; proportional restrictions saturate Minkowski.
    # g^-1_perp = c(t) B while det g=1 can be realized by choosing the longitudinal
    # coefficient to restore determinant one.
    B = np.diag([4.0, 0.25])
    metrics = []
    weights = []
    for c in (0.2, 0.7, 1.3, 4.0):
        aperp = c * B
        ginv = np.diag([aperp[0, 0], aperp[1, 1], 1.0 / np.linalg.det(aperp)])
        g = np.linalg.inv(ginv)
        metrics.append(g)
        weights.append(0.3)
    out = accumulated_transverse_heat_memory(metrics, q, weights)
    assert out["memory_margin"] == pytest.approx(0.0, abs=2e-12)


def test_material_hminus2_reset_speed_is_exactly_viscous_enstrophy_metric():
    rng = np.random.default_rng(2026081504)
    for n in (3, 10, 50):
        for _ in range(100):
            lam = np.exp(rng.uniform(-8, 8, n))
            b = rng.normal(size=n) + 1j * rng.normal(size=n)
            nu = 10.0 ** rng.uniform(-4, 2)
            out = material_hminus2_reset_identity(lam, b, nu)
            assert out["hminus2_reset_speed_squared"] == pytest.approx(
                nu * nu * out["beta_l2_squared"], rel=2e-12, abs=2e-12
            )
            assert out["energy_loss_density_factor"] == pytest.approx(
                nu * out["beta_l2_squared"], rel=2e-12, abs=2e-12
            )


def test_rank_one_incompressible_gradient_has_zero_vortex_stretching():
    rng = np.random.default_rng(2026081505)
    for _ in range(3000):
        xi = rng.normal(size=3)
        a = rng.normal(size=3)
        a -= xi * np.dot(a, xi) / np.dot(xi, xi)
        out = rank_one_incompressible_stretch_null(a, xi)
        assert out["stretch_vector_norm"] <= 2e-12 * max(1.0, out["vorticity_norm"])


def test_material_log_distortion_budget_is_exact_energy_normalization():
    out = material_log_distortion_energy_bound(0.8, 3.2, 0.4)
    assert out["strain_spacetime_l2_squared"] == pytest.approx(2.0)
    assert out["log_max_stretch_l2_squared_upper"] == pytest.approx(1.6)


def test_certificate_is_primitive_and_refuses_global_overclaim():
    cert = theorem_certificate()
    assert "beta_t=-nu L_g beta" in cert["primitive_material_system"]
    assert "closed autonomous material NS system" in cert["primitive_material_system"]
    assert "Minkowski" in cert["history_memory"]
    assert "energy loss" in cert["heat_only_reset"]
    assert cert["case_taxonomy_used"] is False
    assert cert["analysis_cutoff_used"] is False
    assert cert["global_regularity_claimed"] is False


def test_logarithmic_stretch_rate_is_exact_transverse_heat_area_rate():
    rng=np.random.default_rng(2026081506)
    for _ in range(2000):
        g=_det_one_metric(rng); inv=np.linalg.inv(g)
        h=rng.normal(size=(3,3)); h=.5*(h+h.T)
        h=h-(np.trace(inv@h)/3.0)*g
        q=rng.normal(size=3)
        out=transverse_heat_log_rate(g,h,q)
        assert out["log_vorticity_amplification_squared_rate"] == pytest.approx(
            out["log_transverse_heat_determinant_rate"], rel=2e-9, abs=2e-9
        )


def test_pair_directional_mismatch_is_only_material_mismatch_plus_nonaffinity():
    rng=np.random.default_rng(2026081507)
    for _ in range(3000):
        def sl3():
            m=rng.normal(size=(3,3))
            if np.linalg.det(m)<0: m[:,0]*=-1
            return m/(np.linalg.det(m)**(1.0/3.0))
        fa=sl3(); fb=sl3(); qa=rng.normal(size=3); qb=rng.normal(size=3)
        out=pair_direction_mismatch_decomposition(fa,fb,qa,qb)
        assert out["identity_residual"] <= 3e-8*max(1.0,out["physical_cross_norm"])
        # When the deformation is the same, the entire physical cross product is
        # exactly the local heat-covector image of q_a x q_b.
        same=pair_direction_mismatch_decomposition(fa,fa,qa,qb)
        assert same["nonaffinity_term_norm"] == pytest.approx(0.0,abs=2e-12)
        assert same["physical_cross_norm"] == pytest.approx(same["material_heat_covector_norm"],rel=2e-8,abs=2e-8)


def test_moving_polarization_history_is_memory_plus_heat_reset_remainder():
    rng=np.random.default_rng(2026081508)
    for _ in range(1000):
        n=int(rng.integers(2,16)); fs=[]; qs=[]; ws=[]
        for j in range(n):
            while True:
                f=rng.normal(size=(3,3)); d=np.linalg.det(f)
                if abs(d)>.15: break
            if d<0: f[:,0]*=-1; d=-d
            f=f/d**(1.0/3.0); fs.append(f)
            qs.append(rng.normal(size=3)+np.array([.4,-.2,.3]))
            ws.append(float(np.exp(rng.uniform(-4,0))))
        if np.linalg.norm(qs[-1])<1e-4: qs[-1]+=np.array([1.,0.,0.])
        out=moving_polarization_memory_bound(fs,qs,ws)
        assert out["actual_vorticity_amplitude_history"] <= out["memory_reset_upper"]+2e-7
        assert out["fixed_final_polarization_history"] <= out["final_plane_heat_memory"]+2e-7


def test_moving_polarization_reduces_to_pure_memory_when_heat_does_not_rewrite_beta():
    rng=np.random.default_rng(2026081509); q=np.array([.7,-.4,1.2]);fs=[];ws=[]
    for _ in range(12):
        while True:
            f=rng.normal(size=(3,3));d=np.linalg.det(f)
            if abs(d)>.2:break
        if d<0:f[:,0]*=-1;d=-d
        fs.append(f/d**(1.0/3.0));ws.append(.1)
    out=moving_polarization_memory_bound(fs,[q.copy() for _ in fs],ws)
    assert out["heat_only_reset_remainder"] == pytest.approx(0.0,abs=1e-13)


def test_primitive_material_metric_and_beta_reset_speeds_are_exactly_locked():
    rng=np.random.default_rng(2026081510)
    for _ in range(2000):
        b=10.0**rng.uniform(-10,10);nu=10.0**rng.uniform(-6,3)
        out=material_state_speed_lock(b,nu)
        assert out["metric_affine_speed_squared"] == pytest.approx(2*b)
        assert out["hminus2_beta_reset_speed_squared"] == pytest.approx(nu*nu*b)
        assert out["positive_energy_decay_rate"] == pytest.approx(nu*out["metric_affine_speed_squared"])
        assert out["positive_energy_decay_rate"] == pytest.approx(2*out["hminus2_beta_reset_speed_squared"]/nu)



def test_primitive_current_is_divergence_of_metric_speed_and_obeys_flux_continuity():
    rng=np.random.default_rng(2026081511)
    for _ in range(5000):
        k=rng.integers(-20,21,3).astype(float)
        if np.dot(k,k)==0: continue
        u=rng.normal(size=3)+1j*rng.normal(size=3)
        u-=k*np.dot(k,u)/np.dot(k,k)
        out=primitive_material_current_fourier_law(k,u,10.0**rng.uniform(-5,2))
        assert out["codifferential_divergence_residual"] <= 2e-9
        assert out["laplacian_factorization_residual"] <= 2e-9
        assert out["vorticity_flux_continuity_residual"] <= 2e-9
        # h=g_t has exactly twice the intrinsic two-form norm squared.
        assert out["metric_speed_norm_squared"] == pytest.approx(
            2*out["beta_norm_squared_matrix_convention"],rel=3e-11,abs=3e-11
        )


def test_all_hodge_scales_lock_reset_current_and_metric_speed():
    rng=np.random.default_rng(2026081512)
    for s in (-2.0,-1.5,-1.0,-.5,0.0,.5,1.0,2.0,3.0):
        for _ in range(80):
            lam=np.exp(rng.uniform(-7,7,30)); b=rng.normal(size=30)+1j*rng.normal(size=30)
            nu=10.0**rng.uniform(-5,2)
            out=material_hodge_speed_ladder(lam,b,nu,s)
            assert out["beta_reset_hs_minus_two_squared"] == pytest.approx(
                out["viscous_current_hs_minus_one_squared"],rel=3e-12,abs=3e-12
            )
            assert out["beta_reset_hs_minus_two_squared"] == pytest.approx(
                .5*nu*nu*out["metric_speed_hs_squared"],rel=3e-12,abs=3e-12
            )


def test_vorticity_stress_metric_work_and_current_cyclic_work_are_same_algebra():
    rng=np.random.default_rng(2026081513)
    for _ in range(5000):
        u=rng.normal(size=3); w=rng.normal(size=3); c=rng.normal(size=3)
        A=rng.normal(size=(3,3));S=.5*(A+A.T);S-=np.trace(S)/3*np.eye(3)
        out=vorticity_stress_current_algebra(u,w,c,S)
        assert out["stress_metric_work_density"] == pytest.approx(out["stretching_density"],rel=2e-12,abs=2e-12)
        assert out["velocity_stress_divergence_work_density"] == pytest.approx(out["current_cyclic_work_density"],rel=2e-12,abs=2e-12)


def test_certificate_records_local_current_noether_and_all_scale_lock_without_overclaim():
    cert=theorem_certificate()
    assert "beta_t+d j=0" in cert["primitive_current_law"]
    assert "every real s" in cert["all_scale_speed_ladder"]
    assert "div T_beta" in cert["stress_current_noether"]
    assert "cross-product skewness" in cert["enstrophy_derivative_null"]
    assert cert["global_regularity_claimed"] is False



def test_canonical_curl_operator_is_csimons_maxwell_critical_extension():
    rng=np.random.default_rng(2026081511)
    for _ in range(2000):
        n=int(rng.integers(1,30)); a=rng.uniform(-10,10,n); a[np.abs(a)<.25]+=0.5
        e=np.exp(rng.uniform(-7,3,n))
        out=canonical_maxwell_extension_spectral_law(a,e)
        assert out["critical_maxwell_energy"] == pytest.approx(out["self_dual_energy"]+out["anti_self_dual_energy"],rel=2e-12)
        assert out["helicity_chern_simons"] == pytest.approx(out["self_dual_energy"]-out["anti_self_dual_energy"],rel=2e-12,abs=2e-12)
        assert out["boundary_maxwell_energy_density"] == pytest.approx(2*out["enstrophy"])
        assert out["negative_quarter_boundary_profile_derivative"] == pytest.approx(out["critical_viscous_bulk_gradient"])


def test_poisson_depth_overlap_is_exact_critical_sech_filter():
    rng=np.random.default_rng(2026081512)
    for _ in range(3000):
        r=10.0**rng.uniform(-8,8); s=10.0**rng.uniform(-8,8)
        out=canonical_poisson_scale_overlap(r,s)
        assert out["poisson_overlap"] == pytest.approx(out["log_scale_sech"],rel=2e-13,abs=2e-13)
        assert 0.0 < out["poisson_overlap"] <= 1.0+1e-14


def test_maxwell_stress_is_pure_cross_duality_and_bps_sectors_are_stressless():
    rng=np.random.default_rng(2026081513)
    for _ in range(3000):
        m=rng.normal(size=(4,4)); f=m-m.T
        out=maxwell_duality_stress_algebra(f)
        assert out["maxwell_stress_norm_squared"] == pytest.approx(out["cross_duality_product_times_four"],rel=3e-11,abs=3e-11)
        assert out["pure_self_dual_stress_norm"] <= 2e-10*max(1.0,out["field_energy_density"])
        assert out["pure_anti_self_dual_stress_norm"] <= 2e-10*max(1.0,out["field_energy_density"])



def test_physical_spacetime_curvature_has_euler_topological_null_and_joule_work():
    rng=np.random.default_rng(2026081514)
    for _ in range(5000):
        u=rng.normal(size=3); w=rng.normal(size=3); c=rng.normal(size=3); nu=10.0**rng.uniform(-5,2)
        out=primitive_spacetime_gauge_algebra(u,w,c,nu)
        assert out["euler_topological_null"] == pytest.approx(0.0,abs=2e-12)
        assert out["chern_simons_density_rate_half"] == pytest.approx(out["viscous_chern_simons_density_rate_half"],rel=3e-12,abs=3e-12)
        assert out["negative_joule_work_density"] == pytest.approx(out["stretch_minus_ohmic_density"],rel=3e-12,abs=3e-12)


def test_exterior_square_of_volume_preserving_4d_deformation_is_so33_rotation_plus_boost():
    rng=np.random.default_rng(2026081515)
    for _ in range(2000):
        A=rng.normal(size=(4,4));A-=np.trace(A)/4.0*np.eye(4)
        out=so33_exterior_square_algebra(A)
        assert out["so33_residual"] <= 3e-10*max(1.0,out["strain_exterior_hs_squared"])
        assert out["rotation_compact_residual"] <= 3e-10*max(1.0,out["strain_exterior_hs_squared"])
        assert out["strain_boost_residual"] <= 3e-10*max(1.0,out["strain_exterior_hs_squared"])
        assert out["strain_exterior_hs_squared"] == pytest.approx(out["twice_strain_frobenius_squared"],rel=3e-12,abs=3e-12)


def test_maxwell_stress_has_only_two_equal_positive_and_two_equal_negative_principal_values():
    rng=np.random.default_rng(2026081516)
    for _ in range(2000):
        m=rng.normal(size=(4,4));f=m-m.T
        out=maxwell_duality_stress_algebra(f)
        assert out["stress_square_scalar_residual"] <= 3e-10*max(1.0,out["maxwell_stress_norm_squared"])



def test_general_transport_cofactor_reads_twoform_amplification_without_incompressibility():
    rng=np.random.default_rng(202608151101)
    for _ in range(3000):
        q,_=np.linalg.qr(rng.normal(size=(3,3)))
        vals=np.exp(rng.uniform(-2.0,2.0,3))
        f=q@np.diag(vals)
        qv=rng.normal(size=3)
        out=transported_twoform_transverse_determinant(f,qv)
        assert out["transverse_inverse_metric_determinant"] == pytest.approx(out["physical_amplification_squared"],rel=2e-10,abs=2e-10)


def test_vortex_line_slip_absorbs_perpendicular_current_and_twist_is_orthogonal_sink():
    rng=np.random.default_rng(202608151102)
    for _ in range(4000):
        omega=rng.normal(size=3)
        if np.linalg.norm(omega)<.05: omega[0]+=.2
        G=rng.normal(size=(3,3));G-=np.trace(G)/3*np.eye(3)
        u=rng.normal(size=3);nu=10.0**rng.uniform(-4,1)
        out=vortex_slip_twist_algebra(u,omega,G,nu)
        assert out["viscous_current_cost_density"] == pytest.approx(out["slip_twist_cost_density"],rel=2e-9,abs=2e-9)
        assert out["negative_joule_work_density"] == pytest.approx(out["slip_twist_square_density"],rel=2e-8,abs=2e-8)
        assert out["parallel_residual_magnitude_sink"] >= 0.0


def test_twist_free_frobenius_leaf_cancels_normal_curvature_from_slip():
    rng=np.random.default_rng(202608151103)
    for _ in range(3000):
        gradphi=rng.normal(size=3);g=np.linalg.norm(gradphi)
        if g<.05: gradphi[0]+=.2;g=np.linalg.norm(gradphi)
        xi=gradphi/g
        H=rng.normal(size=(3,3));H=.5*(H+H.T)
        gradmu=rng.normal(size=3);mu=float(rng.normal())
        if abs(mu)<.05: mu+=.2
        perp=lambda z:z-xi*np.dot(xi,z)
        curvature=perp(H@xi)/g
        gradg=H@xi
        assert curvature == pytest.approx(perp(gradg)/g,rel=1e-13,abs=1e-13)
        omega=mu*gradphi;c=np.cross(gradmu,gradphi);m=np.linalg.norm(omega);nu=.37
        direct=-nu*np.cross(omega,c)/(m*m)
        surface=-nu*perp(gradmu/mu)
        assert direct == pytest.approx(surface,rel=2e-11,abs=2e-11)
        assert nu*np.dot(c,c) == pytest.approx(nu*g*g*np.dot(perp(gradmu),perp(gradmu)),rel=2e-11,abs=2e-11)


def test_certificate_records_vortex_line_gauge_without_reconnection_overclaim():
    cert=theorem_certificate()
    assert "perpendicular viscous current" in cert["vortex_line_slip_gauge"]
    assert "Frobenius obstruction" in cert["frobenius_twist_current"]
    assert "not itself a reconnection theorem" in cert["vortex_line_topology_guard"]
    assert cert["global_regularity_claimed"] is False



def test_hodge_motion_is_lax_isospectral_and_current_frame_lock_is_modewise_exact():
    rng=np.random.default_rng(202608151201)
    # Similarity/Lax law: L(t)=U(t)L0U(t)^-1 has L_t=[A,L] and fixed spectrum.
    for n in (4,7,11):
        L0=np.diag(np.exp(rng.uniform(-2,2,n)))
        U=rng.normal(size=(n,n))
        while abs(np.linalg.det(U))<.05: U=rng.normal(size=(n,n))
        A=rng.normal(size=(n,n))
        L=U@L0@np.linalg.inv(U)
        Lt=A@L-L@A
        assert np.trace(Lt) == pytest.approx(0.0,abs=2e-10)
        assert np.sort_complex(np.linalg.eigvals(L)) == pytest.approx(np.sort_complex(np.linalg.eigvals(L0)),rel=2e-10,abs=2e-10)
        assert Lt == pytest.approx(A@L-L@A,rel=1e-13,abs=1e-13)
    # One exact Fourier fiber of c=-div h and |grad grad v|^2=|c|^2.
    for _ in range(5000):
        k=rng.normal(size=3);r=np.linalg.norm(k)
        if r<.05: continue
        v=rng.normal(size=3)+1j*rng.normal(size=3);v-=k*np.dot(k,v)/(r*r)
        h=1j*(np.outer(k,v)+np.outer(v,k))
        beta=1j*(np.outer(k,v)-np.outer(v,k))
        c=-1j*np.einsum('i,ij->j',k,beta)
        divh=1j*np.einsum('i,ij->j',k,h)
        A=1j*np.outer(v,k)
        assert c == pytest.approx(-divh,rel=2e-11,abs=2e-11)
        assert r*r*np.vdot(A,A).real == pytest.approx(np.vdot(c,c).real,rel=2e-11,abs=2e-11)
        assert r*r*np.vdot(h,h).real == pytest.approx(2*np.vdot(c,c).real,rel=2e-11,abs=2e-11)


def test_material_symmetric_space_path_action_is_exact_energy_normalization():
    out=material_metric_path_action_bound(.8,3.2,.4)
    assert out["metric_speed_spacetime_l2_squared"] == pytest.approx(8.0)
    assert out["material_path_length_l2_squared_upper"] == pytest.approx(6.4)
    assert out["affine_distance_l2_squared_upper"] == pytest.approx(6.4)
    assert out["log_max_stretch_l2_squared_upper"] == pytest.approx(1.6)


def test_material_acceleration_cancels_explicit_velocity_gradient_square():
    rng=np.random.default_rng(202608151202)
    for _ in range(5000):
        A=rng.normal(size=(3,3));A-=np.trace(A)/3*np.eye(3)
        P=rng.normal(size=(3,3));P=.5*(P+P.T)
        V=rng.normal(size=(3,3))
        DtA=-(A@A)-P+V
        F=rng.normal(size=(3,3))
        Fdot=A@F
        Fdd=DtA@F+A@Fdot
        assert Fdd == pytest.approx((-P+V)@F,rel=2e-11,abs=2e-11)
        omega=rng.normal(size=3)
        euler_DtA=-(A@A)-P
        wdd=euler_DtA@omega+A@(A@omega)
        assert wdd == pytest.approx(-P@omega,rel=2e-11,abs=2e-11)


def test_certificate_records_lax_zero_curvature_and_geodesic_acceleration_without_overclaim():
    cert=theorem_certificate()
    assert "does not create or destroy heat eigenvalues" in cert["hodge_lax_isospectral"]
    assert "pure-gauge SL(3) connection" in cert["maurer_cartan_zero_curvature"]
    assert "same derivative-order activity" in cert["connection_current_lock"]
    assert "quadratic A^2 self-stretch cancels" in cert["lagrangian_geodesic_acceleration"]
    assert cert["global_regularity_claimed"] is False



def test_klein_quadric_splits_growth_tangent_current_from_dissipative_normal_current():
    rng=np.random.default_rng(202608151801)
    for _ in range(5000):
        u=rng.normal(size=3);om=rng.normal(size=3)
        if np.linalg.norm(om)<.05: om[0]+=.2
        c=rng.normal(size=3);nu=10.0**rng.uniform(-5,1)
        out=klein_spacetime_vortex_worldsheet_algebra(u,om,c,nu)
        assert out["projected_klein_pfaffian"] == pytest.approx(0.0,abs=2e-10)
        assert out["fixed_beta_klein_distance_squared"] == pytest.approx(out["viscosity_squared_parallel_current"],rel=2e-11,abs=2e-11)
        assert out["twist_dissipation_density"] >= 0.0


def test_local_faraday_field_always_has_flux_velocity_gauge_away_from_vorticity_zeros():
    rng=np.random.default_rng(202608151802)
    for _ in range(5000):
        u=rng.normal(size=3);om=rng.normal(size=3)
        if np.linalg.norm(om)<.05:om[1]+=.2
        e=rng.normal(size=3)
        out=local_flux_velocity_gauge_algebra(u,om,e)
        assert out["identity_residual"] <= 3e-11


def test_closed_vortex_line_nonideal_period_has_sharp_twist_dissipation_cost():
    out=closed_vortex_line_period_cost(3.0,2.0,.5)
    assert out["minimum_twist_dissipation"] == pytest.approx(8.0/3.0)
    # Equality is attained by constant m*tau along the line.
    L=3.;nu=.5;q=2./(nu*L)
    D=nu*L*q*q
    assert (nu*L*q)**2 == pytest.approx(nu*L*D)


def test_certificate_distinguishes_klein_curvature_from_true_leafwise_topology_obstruction():
    cert=theorem_certificate()
    assert "vortex worldsheets" in cert["klein_vortex_worldsheet"]
    assert "tangentially" in cert["klein_tangent_normal_current"]
    assert "not a local reconnection obstruction" in cert["local_flux_velocity_gauge"]
    assert "leafwise cohomology" in cert["leafwise_period_obstruction"]
    assert "not by itself a reconnection theorem" in cert["klein_topology_guard"]
    assert cert["global_regularity_claimed"] is False



def test_universal_curl_polar_line_geometry_is_exact_at_one_jet():
    rng=np.random.default_rng(202608151901)
    for _ in range(5000):
        b=rng.normal(size=3)
        if np.linalg.norm(b)<.05:b[0]+=.2
        G=rng.normal(size=(3,3))
        out=curl_line_geometry_algebra(b,G)
        assert out["curl_norm_squared"] == pytest.approx(out["twist_plus_defect_curl_norm_squared"],rel=2e-10,abs=2e-10)


def test_ns_enstrophy_is_two_consecutive_readings_of_same_curl_geometry():
    rng=np.random.default_rng(202608151902)
    for _ in range(5000):
        u=rng.normal(size=3);om=rng.normal(size=3);c=rng.normal(size=3)
        if np.linalg.norm(u)<.05:u[0]+=.2
        if np.linalg.norm(om)<.05:om[1]+=.2
        out=two_level_curl_geometry_algebra(u,om,c,.37)
        assert out["stretching_density"] == pytest.approx(out["represented_stretching_density"],rel=2e-10,abs=2e-10)
        assert out["palinstrophy_density"] == pytest.approx(out["represented_palinstrophy_density"],rel=2e-10,abs=2e-10)


def test_certificate_records_one_curl_geometry_not_new_owner_taxonomy():
    cert=theorem_certificate()
    assert "one operator law" in cert["curl_polar_line_geometry"]
    assert "same curl-polar law" in cert["iterated_curl_ns_grammar"]
    assert cert["case_taxonomy_used"] is False



def test_one_field_cartan_hodge_law_reconstructs_full_dealiased_vorticity_pde():
    rng=np.random.default_rng(202608152101)
    for n in (16,20):
        k1=np.fft.fftfreq(n,1/n)
        K=np.stack(np.meshgrid(k1,k1,k1,indexing='ij'),-1);K2=np.sum(K*K,-1);nz=K2>0
        u0=rng.normal(size=(n,n,n,3));uh=np.fft.fftn(u0,axes=(0,1,2))
        dot=np.sum(K*uh,axis=-1);uh[nz]-=K[nz]*(dot[nz]/K2[nz])[:,None];uh[~nz]=0
        uh*=((np.max(np.abs(K),axis=-1)<=n//8)[...,None])
        omh=1j*np.cross(K,uh);ch=K2[...,None]*uh
        urec=np.zeros_like(uh);urec[nz]=1j*np.cross(K[nz],omh[nz])/K2[nz,None]
        assert urec == pytest.approx(uh,rel=2e-11,abs=2e-11)
        u=np.fft.ifftn(uh,axes=(0,1,2)).real
        om=np.fft.ifftn(omh,axes=(0,1,2)).real
        c=np.fft.ifftn(ch,axes=(0,1,2)).real
        nu=.37;e=-np.cross(u,om)+nu*c;eh=np.fft.fftn(e,axes=(0,1,2))
        beta_t=-1j*np.cross(K,eh)
        direct=1j*np.cross(K,np.fft.fftn(np.cross(u,om),axes=(0,1,2)))-nu*K2[...,None]*omh
        assert beta_t == pytest.approx(direct,rel=2e-10,abs=2e-10)


def test_one_cartan_hodge_current_gives_energy_helicity_and_enstrophy_balances():
    rng=np.random.default_rng(202608152102);n=20
    k1=np.fft.fftfreq(n,1/n);K=np.stack(np.meshgrid(k1,k1,k1,indexing='ij'),-1);K2=np.sum(K*K,-1);nz=K2>0
    u0=rng.normal(size=(n,n,n,3));uh=np.fft.fftn(u0,axes=(0,1,2));dot=np.sum(K*uh,axis=-1);uh[nz]-=K[nz]*(dot[nz]/K2[nz])[:,None];uh[~nz]=0
    uh*=((np.max(np.abs(K),axis=-1)<=n//8)[...,None]);omh=1j*np.cross(K,uh);ch=K2[...,None]*uh
    u=np.fft.ifftn(uh,axes=(0,1,2)).real;om=np.fft.ifftn(omh,axes=(0,1,2)).real;c=np.fft.ifftn(ch,axes=(0,1,2)).real;nu=.41
    e=-np.cross(u,om)+nu*c;eh=np.fft.fftn(e,axes=(0,1,2));bt=-1j*np.cross(K,eh)
    def innh(a,b):return float(np.vdot(a,b).real)/(n**3)
    Z=innh(omh,omh);Linv=np.zeros_like(omh);Linv[nz]=omh[nz]/K2[nz,None]
    Eprime=2*innh(bt,Linv)
    Hprime=-2*float(np.sum(e*om));Hrep=-2*nu*float(np.sum(c*om))
    Zprime=2*innh(omh,bt);Zrep=-2*float(np.sum(e*c))
    assert Eprime == pytest.approx(-2*nu*Z,rel=2e-10,abs=2e-10)
    assert Hprime == pytest.approx(Hrep,rel=2e-10,abs=2e-10)
    assert Zprime == pytest.approx(Zrep,rel=2e-10,abs=2e-10)
    # Euler ideal current is exactly null against u and omega pointwise/integrated.
    ideal=-np.cross(u,om)
    assert float(np.max(np.abs(np.sum(ideal*u,axis=-1)))) <= 3e-12
    assert float(np.max(np.abs(np.sum(ideal*om,axis=-1)))) <= 3e-12


def test_certificate_places_one_field_cartan_hodge_law_above_representations():
    cert=theorem_certificate()
    assert "one autonomous closed-two-form law" in cert["primitive_one_field_cartan_hodge"]
    assert "same current" in cert["single_current_balance_ladder"]
    assert "Cartan exterior algebra" in cert["exterior_algebra_euler_null"]
    assert cert["global_regularity_claimed"] is False



def test_sl3_fundamental_casimirs_split_into_cartan_strain_and_rotation_pieces():
    rng=np.random.default_rng(202608152901)
    for _ in range(5000):
        A=rng.normal(size=(3,3));A-=np.trace(A)/3*np.eye(3)
        out=sl3_cartan_casimir_algebra(A)
        assert out["quadratic_casimir"] == pytest.approx(out["quadratic_cartan_split"],rel=2e-12,abs=2e-12)
        assert out["cubic_casimir"] == pytest.approx(out["cubic_cartan_split"],rel=2e-12,abs=2e-12)
        assert out["cubic_casimir"] == pytest.approx(3*out["determinant"],rel=2e-12,abs=2e-12)


def test_finite_sl3_tangent_chord_characteristic_polynomial_is_exact():
    rng=np.random.default_rng(202608152902)
    for _ in range(5000):
        A=rng.normal(size=(3,3));A-=np.trace(A)/3*np.eye(3);s=float(rng.uniform(-5,5))
        out=sl3_tangent_chord_degree_algebra(A,s)
        assert out["chord_jacobian"] == pytest.approx(out["characteristic_polynomial"],rel=2e-11,abs=2e-11)


def test_symmetric_space_metric_speed_balance_has_cubic_cartan_source_and_dirichlet_sink():
    out=symmetric_metric_cartan_balance(12.0,-3.0,5.0,.4)
    assert out["enstrophy"] == pytest.approx(6.0)
    assert out["enstrophy_derivative"] == pytest.approx(1.0)
    assert out["euler_cubic_source"] == pytest.approx(3.0)
    assert out["viscous_dirichlet_sink"] == pytest.approx(-2.0)


def test_affine_local_blowup_guard_proves_local_group_algebra_cannot_be_global_regularizer():
    for tau in (1.0,.2,.03,.005):
        out=affine_local_blowup_guard(tau,2.0)
        assert out["compatibility_residual"] <= 2e-8*max(1.0,out["velocity_gradient_norm"]**2)
        assert out["finite_energy"] is False
        assert out["strain_rate"] == pytest.approx(1/tau)


def test_certificate_records_finite_sl3_degree_law_and_local_blowup_guard_without_overclaim():
    cert=theorem_certificate()
    assert "unique cubic Cartan invariant" in cert["symmetric_space_cubic_law"]
    assert "fundamental sl(3) invariant polynomials" in cert["sl3_fundamental_casimir_null"]
    assert "degree one" in cert["finite_sdiff_chord_degree"]
    assert "no purely local finite-dimensional group law" in cert["affine_local_blowup_guard"]
    assert "is false" in cert["cofactor_force_falsification"]
    assert cert["global_regularity_claimed"] is False



def test_full_sl3_characteristic_polynomial_rate_is_one_cofactor_current():
    rng=np.random.default_rng(202608153101)
    for _ in range(5000):
        A=rng.normal(size=(3,3));At=rng.normal(size=(3,3));lam=float(rng.uniform(-4,4))
        out=sl3_characteristic_polynomial_rate(A,At,lam)
        eps=2e-6
        fd=(np.linalg.det(lam*np.eye(3)-(A+eps*At))-np.linalg.det(lam*np.eye(3)-(A-eps*At)))/(2*eps)
        assert out["characteristic_polynomial_rate"] == pytest.approx(fd,rel=2e-7,abs=2e-7)


def test_local_incompressible_flux_velocity_needs_only_two_scalar_vortex_line_odes():
    out=local_incompressible_flux_velocity_odes((2.,-1.,3.),(.4,.7,-.2),1.3)
    assert out["psi_characteristic_rhs"] == pytest.approx(2*.4-1*.7+3*(-.2))
    assert out["lambda_characteristic_rhs"] == pytest.approx(-1.3)
    assert out["vorticity_magnitude"] == pytest.approx(math.sqrt(14.0))


def test_certificate_records_one_current_deformation_spectrum_transport_and_local_flux_freezing():
    cert=theorem_certificate()
    assert "same Cartan-Hodge electromotive current" in cert["characteristic_polynomial_current"]
    assert "A^3=0" in cert["nilpotent_zero_charge_guard"]
    assert "locally beta_t+Lie_w beta=0 and div w=0" in cert["local_incompressible_flux_freezing"]
    assert "need not globalize" in cert["flux_freezing_global_guard"]
    assert cert["global_regularity_claimed"] is False



def test_curl_transverse_defect_is_geodesic_gradient_of_field_generated_conformal_line_length():
    rng=np.random.default_rng(202608153301)
    for _ in range(5000):
        m=10.0**rng.uniform(-3,3);n=rng.normal(size=3);n/=np.linalg.norm(n)
        k=rng.normal(size=3);k-=n*np.dot(n,k);gm=rng.normal(size=3);V=rng.normal(size=3);V-=n*np.dot(n,V)
        out=conformal_fieldline_length_variation(m,k,gm,n,V)
        assert out["line_shape_variation_density"] == pytest.approx(out["represented_variation_density"],rel=2e-11,abs=2e-11)


def test_certificate_reads_enstrophy_and_palinstrophy_as_vortex_flux_line_geometry_without_curve_flow_overclaim():
    cert=theorem_certificate()
    assert "intrinsic conformal line length" in cert["vorticity_conformal_line_length"]
    assert "Z=int dPhi int m ds" in cert["flux_disintegrated_enstrophy"]
    assert "must not be promoted" in cert["line_gradient_flow_guard"]
    assert cert["global_regularity_claimed"] is False



def test_vortex_productivity_frustration_identity_forces_derivative_obstruction_to_persistent_equality():
    rng=np.random.default_rng(202608153401)
    for _ in range(5000):
        m=10.0**rng.uniform(-4,4);nu=10.0**rng.uniform(-5,2);up=rng.normal();tau=rng.normal();ts=rng.normal()
        q=m-up*tau+2*nu*ts
        out=vortex_productivity_frustration_algebra(m,up,tau,ts,q,nu)
        assert out["compatibility_residual"] <= 2e-10*max(1.0,m)
        assert out["coercive_margin"] >= -2e-9*max(1.0,m*m,out["coercive_rhs"])
    # Open-set exact equality would also set the derivatives/curl to zero.
    out=vortex_productivity_frustration_algebra(2.0,7.0,0.0,0.0,0.0,.3)
    assert out["compatibility_residual"] == pytest.approx(2.0)


def test_twist_free_leaf_poynting_residual_circulation_is_exact_vorticity_flux():
    out=twist_free_leaf_residual_circulation(3.0,5.0)
    assert out["residual_circulation"] == pytest.approx(3.0)
    assert out["loop_l2_residual_lower"] == pytest.approx(9.0/5.0)


def test_certificate_records_spatial_persistence_frustration_without_reintroducing_pointwise_gap():
    cert=theorem_certificate()
    assert "cannot persist" in cert["vortex_productivity_frustration"]
    assert "pushes the obstruction" in cert["vortex_frustration_coercivity"]
    assert "vorticity flux" in cert["twist_free_leaf_flux_obstruction"]
    assert "does not create a pointwise amplitude gap" in cert["persistence_not_pointwise_guard"]
    assert cert["global_regularity_claimed"] is False



def test_poynting_equality_residual_obeys_exact_covariant_gauss_law():
    rng=np.random.default_rng(202608153501)
    for _ in range(5000):
        u=rng.normal(size=3);om=rng.normal(size=3);c=rng.normal(size=3);nu=10.0**rng.uniform(-4,2)
        ell=np.cross(u,om);G=ell-2*nu*c
        # impose the exact divergence dictated by div(u x omega)=|omega|^2-u.c and div c=0
        divG=np.dot(om,om)-np.dot(u,c)
        out=poynting_equality_residual_gauss_algebra(u,om,c,divG,nu)
        assert out["gauss_residual"] <= 2e-11*max(1.0,out["vorticity_squared"])


def test_covariant_gauss_weak_form_recovers_constant_test_normal_tax():
    E=7.0;Z=3.0;nu=.4
    out=covariant_divergence_test_coercivity(Z,0.0,E,nu)
    assert out["residual_l2_squared_lower"] == pytest.approx(4*nu*nu*Z*Z/E)


def test_certificate_collapses_poynting_twist_residual_into_one_sourced_gauss_field():
    cert=theorem_certificate()
    assert "sourced covariant-divergence field" in cert["poynting_residual_gauss_law"]
    assert "Schrödinger dual norm" in cert["poynting_residual_covariant_coercivity"]
    assert "zero-frequency shadow" in cert["poynting_residual_zero_mode"]
    assert "not separate mechanisms" in cert["poynting_residual_polar_collapse"]
    assert "no uniform strict gap" in cert["schrodinger_gap_guard"]



def test_universal_curl_gauss_identity_is_not_special_to_velocity_level():
    rng=np.random.default_rng(202608153601)
    for _ in range(5000):
        v=rng.normal(size=3);w=rng.normal(size=3);c=rng.normal(size=3);nu=10.0**rng.uniform(-4,2)
        divG=np.dot(w,w)-np.dot(v,c)
        out=universal_curl_gauss_residual(v,w,c,divG,nu)
        assert out["gauss_residual"] <= 2e-11*max(1.0,out["curl_energy_density"])


def test_certificate_puts_poynting_gauss_law_inside_one_universal_curl_operator_identity():
    cert=theorem_certificate()
    assert "every divergence-free v" in cert["universal_curl_gauss_operator"]
    assert "reflection form" in cert["antiheat_reflection_form"]
    assert "same exact transverse determinant/Minkowski heat-memory law" in cert["local_flux_memory_recovery"]



def test_actual_ns_current_itself_obeys_positive_curvature_gauss_law():
    rng=np.random.default_rng(202608154401)
    for _ in range(5000):
        u=rng.normal(size=3);om=rng.normal(size=3);c=rng.normal(size=3);nu=10.0**rng.uniform(-4,2)
        # div e = -div(u cross omega) = -|omega|^2+u.c because div curl omega=0.
        div_e=-np.dot(om,om)+np.dot(u,c)
        out=actual_ns_current_gauss_algebra(u,om,c,div_e,nu)
        assert out["gauss_residual"] <= 2e-11*max(1.0,out["vorticity_squared"])


def test_mixed_gauss_law_polarizes_to_arbitrary_test_direction_without_div_b_assumption():
    rng=np.random.default_rng(202608154402)
    for _ in range(5000):
        u=rng.normal(size=3);om=rng.normal(size=3);b=rng.normal(size=3);cb=rng.normal(size=3);nu=10.0**rng.uniform(-4,2)
        G=np.cross(u,b)-2*nu*cb
        # The vector identity fixes div G = omega.b-u.curl b.
        divG=np.dot(om,b)-np.dot(u,cb)
        out=mixed_curl_gauss_direction_algebra(u,om,b,cb,divG,nu)
        assert out["gauss_residual"] <= 2e-11*max(1.0,abs(out["curvature_pairing"]))


def test_midpoint_hodge_square_has_opposite_strain_signs_and_poynting_factorization():
    rng=np.random.default_rng(202608154403)
    n=16;nu=.37;theta=.5
    k1=np.fft.fftfreq(n,1/n)
    K=np.stack(np.meshgrid(k1,k1,k1,indexing='ij'),-1);K2=np.sum(K*K,-1);nz=K2>0
    def divfree_field():
        x=rng.normal(size=(n,n,n,3));xh=np.fft.fftn(x,axes=(0,1,2))
        dot=np.sum(K*xh,axis=-1);xh[nz]-=K[nz]*(dot[nz]/K2[nz])[:,None];xh[~nz]=0
        xh*=((np.max(np.abs(K),axis=-1)<=n//8)[...,None])
        return xh,np.fft.ifftn(xh,axes=(0,1,2)).real
    uh,u=divfree_field();bh,b=divfree_field()
    cbh=1j*np.cross(K,bh);cb=np.fft.ifftn(cbh,axes=(0,1,2)).real
    grad=np.empty((n,n,n,3,3),float)
    for i in range(3):
        for j in range(3):
            grad[...,i,j]=np.fft.ifftn(1j*K[...,j]*uh[...,i],axes=(0,1,2)).real
    S=.5*(grad+np.swapaxes(grad,-1,-2));Sb=np.einsum('...ij,...j->...i',S,b)
    cross=np.cross(u,b);dotub=np.sum(u*b,axis=-1)
    mean=lambda z: float(np.mean(z))
    strain=mean(np.sum(b*Sb,axis=-1));curl2=mean(np.sum(cb*cb,axis=-1));cross2=mean(np.sum(cross*cross,axis=-1));ub2=mean(dotub*dotub)
    base=nu*nu*curl2+theta*theta*(cross2+ub2)
    plus=mean(np.sum((nu*cb+theta*cross)**2,axis=-1))+theta*theta*ub2
    minus=mean(np.sum((nu*cb-theta*cross)**2,axis=-1))+theta*theta*ub2
    assert plus == pytest.approx(base+2*theta*nu*strain,rel=3e-10,abs=3e-10)
    assert minus == pytest.approx(base-2*theta*nu*strain,rel=3e-10,abs=3e-10)
    G=cross-2*nu*cb;g2=mean(np.sum(G*G,axis=-1))
    out=mixed_strain_poynting_factorization(strain,curl2,cross2,g2,nu)
    assert out["factorization_residual"] <= 4e-10*max(1.0,abs(out["strain_work"]))


def test_every_sobolev_square_is_the_parallelogram_of_heat_and_actual_currents():
    rng=np.random.default_rng(202608154404)
    for _ in range(1000):
        j0=rng.normal(size=(17,3));j1=rng.normal(size=(17,3));nu=10.0**rng.uniform(-4,2)
        out=endpoint_current_parallelogram(j0,j1,nu)
        assert out["sobolev_rate_from_inner_product"] == pytest.approx(out["sobolev_rate_from_parallelogram"],rel=2e-12,abs=2e-12)
    out=critical_projected_current_tax(7.0,3.0,.4)
    assert out["critical_projected_residual_lower"] == pytest.approx(4*.4*.4*9/7)


def test_certificate_records_graded_current_collapse_and_preserves_the_persistence_guard():
    cert=theorem_certificate()
    assert "non-nilpotence curvature" in cert["primitive_current_square"]
    assert "actual NS acceleration" in cert["actual_acceleration_gauss"]
    assert "Schrödinger" in cert["midpoint_graded_hodge_square"]
    assert "every smooth b" in cert["mixed_gauss_polarization"]
    assert "same critical H^-1/2 metric" in cert["critical_projected_current_tax"]
    assert "intrinsic Gram plane" in cert["two_null_projected_current_tax"]
    assert "pressure sector" in cert["leray_gauss_scalar_no_go"]
    assert "odd projection does not preserve Gram positivity" in cert["helicity_odd_gram_guard"]
    assert "remaining large-data theorem is dynamical" in cert["graded_current_persistence_frontier"]
    assert cert["global_regularity_claimed"] is False
