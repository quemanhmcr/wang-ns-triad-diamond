from dataclasses import replace
import math

import pytest

from src.cyclic_helical_triad_donor_kernel import cyclic_triad_measure_kernel, signed_good_integer_triad
from src.hard_tail_true_upward_supply import deep_upward_resolved_contact_fixture, hard_tail_upward_supply_split
from src.high_tail_natural_window_reentry import comparable_natural_window_common_work_upper
from src.pure_uv_true_upward_natural_window import (
    PURE_UV_PARENT_UPPER_RATIO,
    PURE_UV_SHELL_INDEX,
    coalesce_pure_uv_recipient_submeasures,
    hard_tail_pure_uv_natural_window_reentry,
    pure_uv_direct_natural_window_reentry,
    pure_uv_first_shell_law,
    pure_uv_natural_window_common_work_upper,
    theorem_certificate,
)


def _signed_good_pure_split():
    triad, _ = signed_good_integer_triad()
    kernel = cyclic_triad_measure_kernel(triad, quotient_measure_mass=1.0)
    return hard_tail_upward_supply_split(triad, kernel, boundary=8.0)


def test_certificate_is_first_shell_canonical_and_has_no_scale_locality_or_rehahn():
    cert = theorem_certificate()
    assert "M=2N" in cert["support_geometry"]
    assert "h=u" in cert["smooth_cutoff"]
    assert "coalesced by recipient" in cert["canonical_cause"]
    assert "p_scale=1" in cert["scale_law"]
    assert "1/2" not in cert["scale_law"]
    assert "9 c sqrt(pi)" in cert["hard_tail_corollary"]
    assert "no output-shell Hahn/locality" in cert["not_used"]
    assert cert["edge_variation_to_clean_young_ratio"] < 1.0
    assert cert["claims_global_regularity"] is False


def test_actual_signed_good_pure_atoms_coalesce_once_on_the_unique_first_shell():
    split = _signed_good_pure_split()
    law = pure_uv_first_shell_law(split)
    assert law.recipient_shell_scale == pytest.approx(16.0)
    assert law.p_scale == 1.0
    assert law.h_inf_output_scale == 0.0
    assert law.total_canonical_positive_mass == pytest.approx(split.pure_uv_hh_physical_work)
    assert law.total_common_unit_work == pytest.approx(8.0 * split.pure_uv_hh_physical_work)
    assert law.coalescing_native_residual <= 1e-15
    assert all(c.recipient_shell_index == PURE_UV_SHELL_INDEX for c in law.recipient_submeasures)
    assert all(c.h_equals_u_on_both_parents for c in law.recipient_submeasures)
    assert all(c.p_scale == 1.0 for c in law.recipient_submeasures)
    assert all(max(c.interaction_parent_radii) <= PURE_UV_PARENT_UPPER_RATIO*c.recipient_shell_scale + 5e-12*c.recipient_shell_scale for c in law.recipient_submeasures)
    assert sum(c.common_unit_work_mass for c in law.recipient_submeasures) == pytest.approx(law.total_common_unit_work)


def test_recipient_cause_cannot_be_cloned_or_declared_uncoalesced():
    law = pure_uv_first_shell_law(_signed_good_pure_split())
    charge = law.recipient_submeasures[0]
    with pytest.raises(ValueError):
        replace(charge, donor_sidecars_coalesced=False)
    with pytest.raises(ValueError):
        replace(charge, owner_mass_cloned=True)
    with pytest.raises(ValueError):
        replace(charge, p_scale=0.5)
    with pytest.raises(ValueError):
        replace(law, p_scale=0.5)
    with pytest.raises(ValueError):
        replace(law, locality_theorem_used=True)


def test_resolved_contact_atom_is_rejected_from_pure_uv_coalescing():
    _, _, split = deep_upward_resolved_contact_fixture()
    contact = next(a for a in split.atoms if a.resolved_scale_parent_contact)
    with pytest.raises(ValueError):
        coalesce_pure_uv_recipient_submeasures((contact,), native_work_mass_scale=split.native_work_scale)


def test_pure_uv_capacity_uses_gross_edge_variation_and_exactly_removes_old_factor_four():
    args = dict(
        window_peak_child_mass=1.7,
        parent_frequency=2.0,
        global_energy=3.0,
        scaled_lifetime=0.8,
    )
    new = pure_uv_natural_window_common_work_upper(**args)
    old = comparable_natural_window_common_work_upper(
        **args, locality_radius=PURE_UV_PARENT_UPPER_RATIO
    )
    assert 4.0 * new == pytest.approx(old)


def test_direct_pure_uv_window_has_no_output_scale_loss_and_exact_first_shell_time():
    law = pure_uv_first_shell_law(_signed_good_pure_split())
    H = law.total_common_unit_work
    N = law.boundary
    E = 2.0
    c = 1.0
    p_time = 0.4
    Ww = p_time * H
    C = PURE_UV_PARENT_UPPER_RATIO * 3.0 * math.sqrt(math.pi) * c * N * E
    mu = (Ww / C) ** 2
    out = pure_uv_direct_natural_window_reentry(
        law,
        required_pure_common_work_lower=H,
        global_energy=E,
        scaled_lifetime=c,
        viscosity=0.1,
        maximum_window_common_work=Ww,
        window_length=c / ((2.0*N)**2),
        window_peak_child_mass=mu,
    )
    assert out["selected_shell_level"] == 1
    assert out["selected_shell_frequency"] == pytest.approx(2.0*N)
    assert out["p_scale"] == 1.0
    assert out["H_inf_output_scale"] == 0.0
    assert out["p_time"] == pytest.approx(p_time)
    assert out["selected_natural_window"] == pytest.approx(c/(4.0*N*N))
    assert out["scale_time_tradeoff_margin"] == pytest.approx(0.0, abs=5e-12*max(1.0,out["weighted_sqrt_child_mass"]))
    assert out["output_scale_locality_theorem_used"] is False
    assert out["recipient_shell_reweighting_used"] is False


def test_hard_tail_pure_owner_corollary_has_exact_one_over_9_coefficient():
    split = _signed_good_pure_split()
    law = pure_uv_first_shell_law(split)
    H = law.total_common_unit_work
    nu = 1.0
    D = H
    N = law.boundary
    E = 2.0
    c = 0.75
    Ww = 0.5 * H
    cap_coeff = PURE_UV_PARENT_UPPER_RATIO * 3.0 * math.sqrt(math.pi) * c * N * E
    mu = (2.0 * Ww / cap_coeff) ** 2
    out = hard_tail_pure_uv_natural_window_reentry(
        split,
        physical_tail_dissipation=D,
        viscosity=nu,
        global_energy=E,
        scaled_lifetime=c,
        maximum_window_common_work=Ww,
        window_length=c / ((2.0*N)**2),
        window_peak_child_mass=mu,
    )
    expected = nu * D / (9.0 * c * math.sqrt(math.pi) * N * E)
    assert out["pure_uv_owner_lower"] == pytest.approx(0.5*nu*D)
    assert out["hard_tail_clean_weighted_sqrt_child_mass_lower"] == pytest.approx(expected)
    assert out["p_scale"] == 1.0


def test_wrong_window_or_later_hahn_fails_closed():
    law = pure_uv_first_shell_law(_signed_good_pure_split())
    H = law.total_common_unit_work
    with pytest.raises(ValueError):
        replace(law, later_hahn_used=True)
    with pytest.raises(ValueError):
        pure_uv_direct_natural_window_reentry(
            law,
            required_pure_common_work_lower=0.5*H,
            global_energy=2.0,
            scaled_lifetime=1.0,
            viscosity=0.1,
            maximum_window_common_work=0.25*H,
            window_length=0.5/(law.boundary**2),
            window_peak_child_mass=1.0,
        )
