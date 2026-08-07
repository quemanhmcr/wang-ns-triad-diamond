import math

import numpy as np

from src.nested_grains import (
    GrainNode,
    Packet,
    TriadEdge,
    build_root,
    extract_tree,
    gaussian_triad_edge,
    make_synthetic_branches,
    scalar_width_ratio,
    split_node,
    verify_nested,
    walk,
)


def test_exact_equal_width_resonant_overlap():
    p = Packet("X", "X", np.array([0.5, 0.4, 0.0]), np.zeros(3), 0.02, 1.0)
    q = Packet("Y", "Y", np.array([0.5, -0.4, 0.0]), np.zeros(3), 0.02, 1.0)
    z = Packet("Z", "Z", p.kappa + q.kappa, np.zeros(3), 0.02, 1.0)
    edge = gaussian_triad_edge(p, q, z)
    assert abs(edge.defect_sq) < 1e-14
    assert abs(edge.weight - 1.0) < 1e-14
    assert abs(scalar_width_ratio([1, 1, 1]) - 1.0) < 1e-14


def test_frequency_and_space_costs_add():
    sigma = 0.1
    p = Packet("X", "X", np.zeros(3), np.zeros(3), sigma, 1.0)
    q = Packet("Y", "Y", np.zeros(3), np.array([1.0, 0.0, 0.0]), sigma, 1.0)
    z = Packet("Z", "Z", np.array([0.2, 0.0, 0.0]), np.zeros(3), sigma, 1.0)
    edge = gaussian_triad_edge(p, q, z)
    expected_freq = 0.2**2 / (12 * sigma**2)
    expected_space = sigma**2 * 2.0 / 3.0
    assert abs(edge.defect_sq - expected_freq - expected_space) < 1e-12
    assert abs(edge.weight - math.exp(-edge.defect_sq)) < 1e-12


def test_annular_certificate_on_adversarial_edges():
    vertices = {"X", "Y", "Z", "X2", "Y2", "Z2"}
    edges = []
    # Populate every annulus and a long tail; weights obey base*exp(-defect_sq).
    for i, defect in enumerate(np.linspace(0.0, 8.0, 80)):
        verts = ("X", "Y", "Z") if i % 2 == 0 else ("X2", "Y2", "Z2")
        base = 1.0 + (i % 5)
        edges.append(TriadEdge(verts, float(defect**2), base, base * math.exp(-defect**2)))
    # Add cross edges so the certificate is exercised.
    for defect in np.linspace(2.0, 9.0, 40):
        base = 0.7
        edges.append(TriadEdge(("X", "Y2", "Z2"), float(defect**2), base, base * math.exp(-defect**2)))
    node = GrainNode("root", 0, vertices, edges, sum(e.weight for e in edges), sum(e.base_weight for e in edges))
    split_node(node, depth=0, bins=12, tail_tolerance=1e-4, bin_width=0.3)
    assert node.cross_loss <= node.certified_relative_bound * node.transfer + 1e-10


def test_nested_tree_and_no_rejoining():
    xs, ys, zs = make_synthetic_branches(6, sigma=0.02, separation=10.0, seed=4)
    root = build_root(xs, ys, zs)
    extract_tree(root, max_depth=4, bin_schedule=[9, 16, 25, 36], tail_schedule=[1e-2, 1e-3, 1e-4, 1e-5], bin_width=0.35)
    assert verify_nested(root)
    for node in walk(root):
        child_sets = [child.vertices for child in node.children]
        for i in range(len(child_sets)):
            for j in range(i + 1, len(child_sets)):
                assert child_sets[i].isdisjoint(child_sets[j])


def test_single_coherent_branch_is_not_artificially_split():
    xs, ys, zs = make_synthetic_branches(1, sigma=0.02, separation=12.0, seed=1)
    root = build_root(xs, ys, zs)
    extract_tree(root, max_depth=5, bin_schedule=[9, 16, 25], tail_schedule=[1e-3, 1e-4, 1e-5], bin_width=0.35)
    assert len(root.children) <= 1
    assert verify_nested(root)


def test_unequal_widths_are_rejected_in_exact_edge_module():
    p = Packet("Xw", "X", np.zeros(3), np.zeros(3), 0.02, 1.0)
    q = Packet("Yw", "Y", np.zeros(3), np.zeros(3), 0.03, 1.0)
    z = Packet("Zw", "Z", np.zeros(3), np.zeros(3), 0.02, 1.0)
    import pytest
    with pytest.raises(ValueError):
        gaussian_triad_edge(p, q, z)
