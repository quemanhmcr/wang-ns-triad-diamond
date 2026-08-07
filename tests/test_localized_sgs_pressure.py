from src.localized_sgs_pressure import (
    ResolvedWindowStep,
    critical_boundary_charge_lower_bound,
    infer_combined_work,
    positive_work_depletion_bound,
    pressure_cancellation_trichotomy,
    weighted_chain_identity,
)


def test_localized_resolved_energy_identity():
    w = infer_combined_work(3.0, 2.1, 0.4, -0.2)
    s = ResolvedWindowStep(3.0, 2.1, w, 0.4, -0.2)
    assert abs(s.residual) < 1e-15


def test_pressure_cancellation_trichotomy():
    a = pressure_cancellation_trichotomy(2.0, 1.2)
    assert a["branch"] == "combined_work"
    b = pressure_cancellation_trichotomy(2.0, 0.3)
    assert b["branch"] == "pressure_boundary"
    assert b["pressure_cancellation"] >= 1.0
    assert abs(critical_boundary_charge_lower_bound(2.0, 4.0) - 0.25) < 1e-15


def test_weighted_chain_telescopes_and_depletes():
    energies = [4.0, 3.0, 2.5, 1.7]
    weights = [1.0, 0.8, 0.5]
    steps = []
    for j, (d, l) in enumerate([(0.2, 0.1), (0.3, -0.05), (0.1, 0.2)]):
        w = infer_combined_work(energies[j], energies[j+1], d, l)
        steps.append(ResolvedWindowStep(energies[j], energies[j+1], w, d, l))
    c = weighted_chain_identity(steps, weights)
    assert abs(c["residual"]) < 1e-14
    dep = positive_work_depletion_bound(steps, weights)
    assert dep["margin"] >= -1e-14
