import math

import numpy as np

from src.high_strain_resolved_ancestor import high_strain_ancestor_mass_threshold
from src.nn_critical_heat_carrier_seed import (
    NNCriticalHeatAtom,
    normalized_shell_probe_coefficient,
    persistent_seed_low_low_gap,
    pushforward_nn_critical_heat_law,
    renewal_carrier_critical_mass_lower,
    renewal_natural_lifetime_ratio,
    renewal_scale,
    shell_relative_support,
    theorem_certificate,
)


def test_renewal_scale_places_hard_shell_inside_outer_role_geometry():
    M = 20.0
    A = renewal_scale(M)
    assert math.isclose(A, 15.0)
    lo, hi = shell_relative_support()
    assert math.isclose(lo, 2 / 3)
    assert math.isclose(hi, 4 / 3)
    assert lo > 3 / 5
    assert hi < 3 / 2


def test_critical_shell_becomes_clean_renewed_critical_carrier():
    c = 1.7
    mu = high_strain_ancestor_mass_threshold(c)
    M = 12.0
    E = mu / M
    A = renewal_scale(M)
    assert A * E >= renewal_carrier_critical_mass_lower(c) - 1e-14
    assert math.isclose(renewal_carrier_critical_mass_lower(c), 8 * math.pi**2 / (25 * c**2))


def test_renewed_scale_has_longer_natural_lifetime_than_child():
    N = 80.0
    M = N / 4.0
    assert math.isclose(renewal_natural_lifetime_ratio(N, M), 256 / 9)


def test_smooth_seed_keeps_strict_low_low_moat():
    assert persistent_seed_low_low_gap(1 / 30) > 0


def test_shell_own_direction_registers_exactly_without_packet_choice():
    f = np.array([1 + 2j, -3 + 0.5j, 0.7 - 1.1j])
    out = normalized_shell_probe_coefficient(f, np.ones(3))
    assert math.isclose(out["coefficient_energy"], out["shell_energy"])
    assert out["registration_residual"] < 1e-13


def test_nn_critical_heat_law_pushes_forward_without_atom_floor():
    c = 1.0
    N = 128.0
    mu = high_strain_ancestor_mass_threshold(c)
    atoms = (
        NNCriticalHeatAtom(1e-9, N, N / 4, 1.2 * mu / (N / 4), 0.001),
        NNCriticalHeatAtom(3.0, N, N / 8, 2.0 * mu / (N / 8), 0.002),
        NNCriticalHeatAtom(0.4, N, N / 16, 1.1 * mu / (N / 16), 0.003),
    )
    seeds = pushforward_nn_critical_heat_law(atoms, scaled_lifetime=c)
    assert math.isclose(sum(s.probability for s in seeds), 1.0)
    assert math.isclose(sum(s.heat_mass for s in seeds), sum(a.mass for a in atoms))
    assert min(s.renewal_critical_mass for s in seeds) >= renewal_carrier_critical_mass_lower(c) - 1e-12
    assert seeds[0].probability > 0


def test_certificate_keeps_v_heat_and_u_shell_marks_distinct():
    cert = theorem_certificate()
    assert "distinct exact marks" in cert["material_provenance"]
    assert "does not yet extend" in cert["scope"]
    assert "full A-natural slab" in cert["scope"]
