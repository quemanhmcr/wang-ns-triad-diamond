from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class Packet:
    """One Gaussian atom in the scale-critical atomic packet model."""

    name: str
    side: str
    kappa: np.ndarray
    center: np.ndarray
    sigma: float
    coefficient: float


@dataclass(frozen=True)
class TriadEdge:
    """A positive scalar Gaussian triad contribution.

    `base_weight` is the coefficient product times the exact width factor.
    `weight = base_weight * exp(-defect_sq)` is the exact scalar envelope
    contribution in the stated Gaussian model.
    """

    vertices: tuple[str, str, str]
    defect_sq: float
    base_weight: float
    weight: float


@dataclass
class GrainNode:
    node_id: str
    depth: int
    vertices: set[str]
    edges: list[TriadEdge]
    transfer: float
    base_budget: float
    children: list["GrainNode"] = field(default_factory=list)
    selected_bin: int | None = None
    lower_cut: float | None = None
    upper_cut: float | None = None
    moat_loss: float = 0.0
    tail_loss: float = 0.0
    cross_loss: float = 0.0
    certified_relative_bound: float = 0.0


class UnionFind:
    def __init__(self, items: Iterable[str]):
        self.parent = {x: x for x in items}
        self.rank = {x: 0 for x in items}

    def find(self, x: str) -> str:
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
        return self.parent[x]

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def scalar_width_ratio(sigmas: Sequence[float], dimension: int = 3) -> float:
    s = np.asarray(sigmas, dtype=float)
    if s.shape != (3,) or np.min(s) <= 0:
        raise ValueError("three positive widths required")
    a = s * s / float(np.sum(s * s))
    return float((27.0 * np.prod(a)) ** (dimension / 4.0))


def gaussian_triad_edge(a: Packet, b: Packet, c: Packet) -> TriadEdge:
    if (a.side, b.side, c.side) != ("X", "Y", "Z"):
        raise ValueError("edge packets must be ordered X,Y,Z")
    sigma_ref = math.sqrt((a.sigma * a.sigma + b.sigma * b.sigma + c.sigma * c.sigma) / 3.0)
    mismatch = a.kappa + b.kappa - c.kappa
    freq_cost = float(mismatch @ mismatch) / (12.0 * sigma_ref * sigma_ref)
    xs = (a.center, b.center, c.center)
    pair = sum(float((xs[i] - xs[j]) @ (xs[i] - xs[j])) for i in range(3) for j in range(i + 1, 3))
    spatial_cost = (sigma_ref * sigma_ref / 3.0) * pair
    width = scalar_width_ratio([a.sigma, b.sigma, c.sigma])
    width_cost = -math.log(max(width, 1e-300))
    defect_sq = max(0.0, freq_cost + spatial_cost + width_cost)
    base = abs(a.coefficient * b.coefficient * c.coefficient) * width
    return TriadEdge(
        vertices=(a.name, b.name, c.name),
        defect_sq=defect_sq,
        base_weight=base,
        weight=base * math.exp(-defect_sq),
    )


def all_gaussian_edges(x_packets: Sequence[Packet], y_packets: Sequence[Packet], z_packets: Sequence[Packet]) -> list[TriadEdge]:
    return [gaussian_triad_edge(a, b, c) for a in x_packets for b in y_packets for c in z_packets]


def _components(vertices: set[str], edges: Sequence[TriadEdge], cutoff_sq: float) -> dict[str, set[str]]:
    uf = UnionFind(vertices)
    for edge in edges:
        if edge.defect_sq < cutoff_sq:
            u, v, w = edge.vertices
            uf.union(u, v)
            uf.union(v, w)
    result: dict[str, set[str]] = {}
    for vertex in vertices:
        result.setdefault(uf.find(vertex), set()).add(vertex)
    return result


def split_node(
    node: GrainNode,
    *,
    depth: int,
    bins: int,
    tail_tolerance: float,
    bin_width: float = 0.5,
    min_child_transfer: float = 1e-15,
) -> GrainNode:
    """Split a node by a transfer-selected annulus in Gaussian defect space.

    The start radius is adaptive: exp(-L^2) * base_budget / transfer <= tail_tolerance.
    Among `bins` adjacent annuli one has transfer at most transfer/bins.
    Short edges generate connected components. Any cross-component edge is then
    either in the selected moat or beyond its upper edge, where the Gaussian
    envelope gives the tail estimate.
    """
    if node.transfer <= 0 or not node.edges or bins < 1:
        return node
    ratio = max(1.0, node.base_budget / max(node.transfer, 1e-300))
    start_sq = max(0.0, math.log(ratio / max(tail_tolerance, 1e-300)))
    start = math.sqrt(start_sq)
    annulus_masses: list[float] = []
    for m in range(bins):
        lo = start + m * bin_width
        hi = lo + bin_width
        annulus_masses.append(sum(e.weight for e in node.edges if lo * lo <= e.defect_sq < hi * hi))
    selected = int(np.argmin(np.asarray(annulus_masses)))
    lower = start + selected * bin_width
    upper = lower + bin_width
    lower_sq, upper_sq = lower * lower, upper * upper

    components = _components(node.vertices, node.edges, lower_sq)
    vertex_to_component = {v: root for root, members in components.items() for v in members}

    internal: dict[str, list[TriadEdge]] = {root: [] for root in components}
    moat_loss = 0.0
    tail_loss = 0.0
    cross_loss = 0.0
    for edge in node.edges:
        roots = {vertex_to_component[v] for v in edge.vertices}
        if len(roots) == 1:
            internal[next(iter(roots))].append(edge)
            continue
        cross_loss += edge.weight
        if lower_sq <= edge.defect_sq < upper_sq:
            moat_loss += edge.weight
        elif edge.defect_sq >= upper_sq:
            tail_loss += edge.weight
        else:
            raise AssertionError("short cross edge contradicts component construction")

    children: list[GrainNode] = []
    for idx, (root, members) in enumerate(sorted(components.items(), key=lambda kv: sorted(kv[1]))):
        child_edges = internal[root]
        transfer = sum(e.weight for e in child_edges)
        if transfer <= min_child_transfer:
            continue
        children.append(
            GrainNode(
                node_id=f"{node.node_id}.{idx}",
                depth=depth + 1,
                vertices=set(members),
                edges=child_edges,
                transfer=transfer,
                base_budget=sum(e.base_weight for e in child_edges),
            )
        )

    moat_bound = node.transfer / bins
    tail_bound = node.base_budget * math.exp(-upper_sq)
    node.children = children
    node.selected_bin = selected
    node.lower_cut = lower
    node.upper_cut = upper
    node.moat_loss = moat_loss
    node.tail_loss = tail_loss
    node.cross_loss = cross_loss
    node.certified_relative_bound = (moat_bound + tail_bound) / node.transfer

    # Exact finite-model certificates.
    if moat_loss > moat_bound + 1e-11 * max(1.0, node.transfer):
        raise AssertionError("pigeonhole moat bound failed")
    if tail_loss > tail_bound + 1e-11 * max(1.0, node.base_budget):
        raise AssertionError("Gaussian tail bound failed")
    if cross_loss > moat_bound + tail_bound + 1e-11 * max(1.0, node.transfer):
        raise AssertionError("cross-transfer certificate failed")
    return node


def extract_tree(
    root: GrainNode,
    *,
    max_depth: int,
    bin_schedule: Sequence[int],
    tail_schedule: Sequence[float],
    bin_width: float = 0.5,
) -> GrainNode:
    frontier = [root]
    for depth in range(max_depth):
        new_frontier: list[GrainNode] = []
        bins = bin_schedule[min(depth, len(bin_schedule) - 1)]
        tail_tol = tail_schedule[min(depth, len(tail_schedule) - 1)]
        for node in frontier:
            split_node(node, depth=depth, bins=bins, tail_tolerance=tail_tol, bin_width=bin_width)
            # A one-child node is already sticky at this resolution; recurse only
            # when an actual decomposition was produced.
            if len(node.children) > 1:
                new_frontier.extend(node.children)
        frontier = new_frontier
        if not frontier:
            break
    return root


def walk(root: GrainNode) -> list[GrainNode]:
    out: list[GrainNode] = []
    stack = [root]
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(reversed(node.children))
    return out


def level_loss_bounds(root: GrainNode) -> dict[int, dict[str, float]]:
    stats: dict[int, dict[str, float]] = {}
    for node in walk(root):
        if node.lower_cut is None:
            continue
        row = stats.setdefault(node.depth, {"parent_transfer": 0.0, "cross_loss": 0.0, "certified_loss": 0.0, "nodes": 0.0})
        row["parent_transfer"] += node.transfer
        row["cross_loss"] += node.cross_loss
        row["certified_loss"] += node.certified_relative_bound * node.transfer
        row["nodes"] += 1.0
    return stats


def verify_nested(root: GrainNode) -> bool:
    for node in walk(root):
        seen: set[str] = set()
        for child in node.children:
            if not child.vertices <= node.vertices:
                return False
            if seen & child.vertices:
                return False
            seen |= child.vertices
        if any(not set(e.vertices) <= node.vertices for e in node.edges):
            return False
    return True


def make_synthetic_branches(branches: int, sigma: float, separation: float, seed: int = 0) -> tuple[list[Packet], list[Packet], list[Packet]]:
    """Create coherent Gaussian triads whose branches are separated in phase space."""
    rng = np.random.default_rng(seed)
    rstar = 0.6109041018281888
    h = math.sqrt(rstar * rstar - 0.25)
    p0 = np.array([0.5, h, 0.0])
    q0 = np.array([0.5, -h, 0.0])
    z0 = p0 + q0
    xs: list[Packet] = []
    ys: list[Packet] = []
    zs: list[Packet] = []
    for i in range(branches):
        # Generic frequency and spatial displacements avoid accidental shared children.
        direction = np.array([math.cos(1.3 * i), math.sin(1.3 * i), 0.35 * i])
        direction /= max(np.linalg.norm(direction), 1.0)
        dk = separation * sigma * direction
        dx = (separation / sigma) * np.roll(direction, 1)
        jitter = 0.02 * sigma * rng.normal(size=3)
        coeff = 1.0 / (branches ** (2.0 / 3.0))
        xs.append(Packet(f"X{i}", "X", p0 + dk + jitter, dx, sigma, coeff))
        ys.append(Packet(f"Y{i}", "Y", q0 + dk - jitter, dx, sigma, coeff))
        zs.append(Packet(f"Z{i}", "Z", z0 + 2.0 * dk, dx, sigma, coeff))
    return xs, ys, zs


def build_root(xs: Sequence[Packet], ys: Sequence[Packet], zs: Sequence[Packet]) -> GrainNode:
    edges = all_gaussian_edges(xs, ys, zs)
    vertices = {p.name for p in [*xs, *ys, *zs]}
    return GrainNode(
        node_id="root",
        depth=0,
        vertices=vertices,
        edges=edges,
        transfer=sum(e.weight for e in edges),
        base_budget=sum(e.base_weight for e in edges),
    )


def experiment(outdir: Path) -> dict[str, object]:
    outdir.mkdir(parents=True, exist_ok=True)
    scenarios = []
    for branches in (1, 2, 4, 8):
        for separation in (4.0, 8.0, 12.0):
            xs, ys, zs = make_synthetic_branches(branches, sigma=0.02, separation=separation, seed=branches)
            root = build_root(xs, ys, zs)
            extract_tree(
                root,
                max_depth=6,
                bin_schedule=[9, 16, 25, 36, 49, 64],
                tail_schedule=[1 / 16, 1 / 64, 1 / 256, 1 / 1024, 1 / 4096, 1 / 16384],
                bin_width=0.35,
            )
            levels = level_loss_bounds(root)
            scenarios.append(
                {
                    "branches": branches,
                    "separation": separation,
                    "root_transfer": root.transfer,
                    "root_base_budget": root.base_budget,
                    "root_children": len(root.children),
                    "tree_nodes": len(walk(root)),
                    "nested": verify_nested(root),
                    "total_cross_loss": sum(v["cross_loss"] for v in levels.values()),
                    "total_certified_loss": sum(v["certified_loss"] for v in levels.values()),
                    "levels": levels,
                }
            )
    result = {
        "model": "finite L^(3/2)-atomic Gaussian triad model",
        "scenarios": scenarios,
        "schedules": {
            "bins": [9, 16, 25, 36, 49, 64],
            "tail": [1 / 16, 1 / 64, 1 / 256, 1 / 1024, 1 / 4096, 1 / 16384],
        },
    }
    (outdir / "nested_grains.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    lines = [
        "# Nested Gaussian grain extraction",
        "",
        "The annular moat and Gaussian tail certificates are exact in the stated finite atomic model.",
        "",
        "| branches | separation | root children | tree nodes | cross loss | certified bound | nested |",
        "|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in scenarios:
        lines.append(
            f"| {row['branches']} | {row['separation']:.1f} | {row['root_children']} | {row['tree_nodes']} | "
            f"{row['total_cross_loss']:.3e} | {row['total_certified_loss']:.3e} | {row['nested']} |"
        )
    lines += [
        "",
        "This is not a Navier--Stokes regularity theorem. The missing bridge is an atomic Gaussian extraction from an arbitrary near-extremal block with controlled synthesis constants.",
    ]
    (outdir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("results-nested-grains"))
    args = parser.parse_args()
    experiment(args.outdir)


if __name__ == "__main__":
    main()
