from src.nested_grains import TriadEdge
from src.reservoir_pool_erosion import (
    exact_pool_certificate,
    first_forced_relink_generation,
    old_pool_service_capacity_upper,
    relinking_incidence_route,
    total_old_pool_service_upper,
)


def test_old_pool_half_life():
    c0 = old_pool_service_capacity_upper(0, 0.25, 4.0, 2.0, 3.0)
    for q in range(12):
        assert old_pool_service_capacity_upper(q, 0.25, 4.0, 2.0, 3.0) <= c0 * 0.5**q + 1e-14
    assert total_old_pool_service_upper(0.25, 4.0, 2.0, 3.0) < 2 * c0


def test_uniform_service_forces_relink_eventually():
    q = first_forced_relink_generation(0.1, 0.25, 4.0, 2.0, 3.0)
    assert q >= 0
    assert old_pool_service_capacity_upper(q, 0.25, 4.0, 2.0, 3.0) < 0.05


def test_relink_incidence_fresh_or_cycle():
    edges = [
        TriadEdge(("a", "b", "c"), 0.0, 1.0, 1.0),
        TriadEdge(("c", "d", "e"), 0.0, 1.0, 1.0),
        TriadEdge(("e", "a", "f"), 0.0, 1.0, 1.0),
    ]
    row = relinking_incidence_route(edges)
    assert max(row["fresh_units"], row["cycle_rank"]) >= row["triads"]


def test_certificate():
    c = exact_pool_certificate()
    assert c["clean_old_pool_ratio"] == "<1/2"
