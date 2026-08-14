import math

import numpy as np
import pytest

from src.native_material_vorticity_heat_law import (
    accumulated_transverse_heat_memory,
    material_hminus2_reset_identity,
    material_log_distortion_energy_bound,
    material_state_speed_lock,
    moving_polarization_memory_bound,
    pair_direction_mismatch_decomposition,
    rank_one_incompressible_stretch_null,
    theorem_certificate,
    transverse_heat_determinant,
    transverse_heat_log_rate,
    transverse_two_covector_area_identity,
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
