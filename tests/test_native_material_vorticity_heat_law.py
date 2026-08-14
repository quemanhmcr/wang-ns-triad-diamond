import math

import numpy as np
import pytest

from src.native_material_vorticity_heat_law import (
    accumulated_transverse_heat_memory,
    canonical_maxwell_extension_spectral_law,
    canonical_poisson_scale_overlap,
    maxwell_duality_stress_algebra,
    material_hminus2_reset_identity,
    material_hodge_speed_ladder,
    material_log_distortion_energy_bound,
    material_state_speed_lock,
    primitive_material_current_fourier_law,
    moving_polarization_memory_bound,
    pair_direction_mismatch_decomposition,
    rank_one_incompressible_stretch_null,
    theorem_certificate,
    transverse_heat_determinant,
    transverse_heat_log_rate,
    transverse_two_covector_area_identity,
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
