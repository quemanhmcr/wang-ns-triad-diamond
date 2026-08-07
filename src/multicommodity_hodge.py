from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


@dataclass(frozen=True)
class DirectedEdge:
    tail: int
    head: int
    conductance: float
    target: float


def incidence(nv: int, edges: list[DirectedEdge]) -> np.ndarray:
    D = np.zeros((nv, len(edges)), dtype=float)
    for j, e in enumerate(edges):
        D[e.tail, j] = -1.0
        D[e.head, j] = 1.0
    return D


def hodge_energy(nv: int, edges: list[DirectedEdge]) -> dict[str, object]:
    D = incidence(nv, edges)
    c = np.array([e.conductance for e in edges], dtype=float)
    a = np.array([e.target for e in edges], dtype=float)
    W = np.diag(c)
    L = D @ W @ D.T
    rhs = D @ W @ a
    phi = np.linalg.pinv(L, rcond=1e-13) @ rhs
    residual = a - D.T @ phi
    energy = float(np.sum(c * residual * residual))
    div = D @ (c * residual)
    return {
        "energy": energy,
        "potential": phi,
        "residual": residual,
        "weighted_divergence_norm": float(np.linalg.norm(div)),
        "D": D,
        "conductance": c,
        "target": a,
    }


def electrical_flow(nv: int, edges: list[DirectedEdge], s: int, t: int) -> tuple[np.ndarray, float]:
    D = incidence(nv, edges)
    c = np.array([e.conductance for e in edges], dtype=float)
    W = np.diag(c)
    L = D @ W @ D.T
    b = np.zeros(nv, dtype=float)
    b[s] = -1.0
    b[t] = 1.0
    psi = np.linalg.pinv(L, rcond=1e-13) @ b
    f = W @ D.T @ psi
    energy = float(np.sum((f * f) / c))
    if np.linalg.norm(D @ f - b) > 1e-8:
        raise RuntimeError("electrical flow divergence mismatch")
    return f, energy


def multicommodity_rayleigh(
    nv: int,
    edges: list[DirectedEdge],
    cycles: list[np.ndarray],
    masses: np.ndarray,
) -> dict[str, float]:
    data = hodge_energy(nv, edges)
    c = np.asarray(data["conductance"], dtype=float)
    a = np.asarray(data["target"], dtype=float)
    D = np.asarray(data["D"], dtype=float)
    masses = np.asarray(masses, dtype=float)
    num = 0.0
    den = 0.0
    for mu, z in zip(masses, cycles):
        z = np.asarray(z, dtype=float)
        if np.linalg.norm(D @ z) > 1e-8:
            raise ValueError("commodity is not a cycle flow")
        circ = float(z @ a)
        resistance_energy = float(np.sum((z * z) / c))
        num += float(mu) * circ * circ
        den += float(mu) * resistance_energy
    bound = 0.0 if den == 0.0 else num / den
    return {
        "hodge_energy": float(data["energy"]),
        "weighted_circulation_sq": num,
        "weighted_cycle_resistance": den,
        "rayleigh_bound": bound,
        "margin": float(data["energy"]) - bound,
    }


def colored_union(
    nv: int,
    old_edges: list[DirectedEdge],
    new_edges: list[DirectedEdge],
) -> list[DirectedEdge]:
    return list(old_edges) + list(new_edges)


def pair_cycle_electrical(
    nv: int,
    old_edges: list[DirectedEdge],
    new_edges: list[DirectedEdge],
    s: int,
    t: int,
) -> tuple[np.ndarray, float, float, float]:
    f_old, r_old = electrical_flow(nv, old_edges, s, t)
    f_new, r_new = electrical_flow(nv, new_edges, s, t)
    z = np.concatenate([-f_old, f_new])
    return z, r_old + r_new, r_old, r_new


def pair_weighted_union_certificate(
    nv: int,
    old_edges: list[DirectedEdge],
    new_edges: list[DirectedEdge],
    terminals: list[int],
    probs: np.ndarray,
) -> dict[str, float]:
    probs = np.asarray(probs, dtype=float)
    probs = probs / probs.sum()
    union = colored_union(nv, old_edges, new_edges)
    cycles: list[np.ndarray] = []
    masses: list[float] = []
    rbar = 0.0
    for ai, s in enumerate(terminals):
        for bi, t in enumerate(terminals):
            mu = float(probs[ai] * probs[bi])
            if s == t or mu == 0.0:
                continue
            z, r, _, _ = pair_cycle_electrical(nv, old_edges, new_edges, s, t)
            cycles.append(z)
            masses.append(mu)
            rbar += mu * r
    cert = multicommodity_rayleigh(nv, union, cycles, np.array(masses))
    cert["pair_resistance_budget"] = rbar
    cert["pair_mass"] = float(sum(masses))
    return cert


def exact_gradient_edges(nv: int, undirected_pairs: list[tuple[int, int]], heights: np.ndarray, conductance: float = 1.0) -> list[DirectedEdge]:
    edges: list[DirectedEdge] = []
    for u, v in undirected_pairs:
        if heights[u] <= heights[v]:
            tail, head = u, v
        else:
            tail, head = v, u
        target = float(heights[head] - heights[tail])
        edges.append(DirectedEdge(tail, head, conductance, target))
    return edges


def synchronization_certificate(
    nv: int,
    old_pairs: list[tuple[int, int]],
    new_pairs: list[tuple[int, int]],
    h_old: np.ndarray,
    h_new: np.ndarray,
    terminals: list[int],
    probs: np.ndarray,
    conductance: float = 1.0,
) -> dict[str, float]:
    old_edges = exact_gradient_edges(nv, old_pairs, h_old, conductance)
    new_edges = exact_gradient_edges(nv, new_pairs, h_new, conductance)
    cert = pair_weighted_union_certificate(nv, old_edges, new_edges, terminals, probs)
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    d = np.array([h_new[v] - h_old[v] for v in terminals], float)
    mean = float(np.sum(probs * d))
    var = float(np.sum(probs * (d - mean) ** 2))
    numerator_exact = 2.0 * var
    cert["shift_variance"] = var
    cert["pair_shift_sq"] = numerator_exact
    cert["synchronization_bound"] = 0.5 * cert["hodge_energy"] * cert["pair_resistance_budget"]
    cert["sync_margin"] = cert["synchronization_bound"] - var
    return cert


def integer_gauge_mismatch_bound(variance: float, gamma: float) -> float:
    # If shifts belong to gamma*Z then P[d_I != d_J] <= E[(dI-dJ)^2]/gamma^2 = 2 Var/gamma^2.
    return min(1.0, max(0.0, 2.0 * variance / (gamma * gamma)))


def tree_pair_resistance_identity(
    n: int,
    tree_edges: list[tuple[int, int, float]],
    probs: np.ndarray,
) -> dict[str, float]:
    probs = np.asarray(probs, float)
    probs = probs / probs.sum()
    adj: list[list[tuple[int, int]]] = [[] for _ in range(n)]
    for idx, (u, v, resistance) in enumerate(tree_edges):
        adj[u].append((v, idx))
        adj[v].append((u, idx))
    if len(tree_edges) != n - 1:
        raise ValueError("tree must have n-1 edges")
    contribs = []
    for idx, (u, v, r) in enumerate(tree_edges):
        stack = [u]
        seen = {v}
        side = set()
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            side.add(x)
            for y, _ in adj[x]:
                if y not in seen:
                    stack.append(y)
        s = float(sum(probs[i] for i in side))
        contribs.append(2.0 * r * s * (1.0 - s))
    formula = float(sum(contribs))
    # direct all-pairs path resistance via Floyd-Warshall on a tree-sized matrix
    dist = np.full((n, n), np.inf)
    np.fill_diagonal(dist, 0.0)
    for u, v, r in tree_edges:
        dist[u, v] = dist[v, u] = r
    for k in range(n):
        dist = np.minimum(dist, dist[:, [k]] + dist[[k], :])
    direct = float(np.sum(probs[:, None] * probs[None, :] * dist))
    return {"direct": direct, "formula": formula, "margin": direct - formula, "max_edge_contribution": max(contribs, default=0.0)}


def _toy_graphs(gamma: float) -> dict[str, object]:
    # terminals 0,1,2. Old and new trees encode either a synchronized shift or a mismatch.
    nv = 5
    old_pairs = [(0, 3), (1, 3), (3, 4), (2, 4)]
    new_pairs = [(0, 3), (2, 3), (3, 4), (1, 4)]
    h_old = gamma * np.array([0, 0, 0, 1, 2], float)
    h_new_sync = gamma * np.array([1, 1, 1, 2, 3], float)
    h_new_mismatch = gamma * np.array([1, 2, 1, 2, 3], float)
    probs = np.array([0.45, 0.35, 0.20])
    sync = synchronization_certificate(nv, old_pairs, new_pairs, h_old, h_new_sync, [0,1,2], probs)
    curved = synchronization_certificate(nv, old_pairs, new_pairs, h_old, h_new_mismatch, [0,1,2], probs)
    sync["integer_mismatch_probability_bound"] = integer_gauge_mismatch_bound(sync["shift_variance"], gamma)
    curved["integer_mismatch_probability_bound"] = integer_gauge_mismatch_bound(curved["shift_variance"], gamma)
    return {"synchronized": sync, "mismatched": curved}


def random_rayleigh_checks(samples: int, seed: int = 7) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    worst = 1e9
    checked = 0
    for _ in range(samples):
        nv = int(rng.integers(4, 9))
        # start from a cycle graph so cycle space is nontrivial
        edges: list[DirectedEdge] = []
        for v in range(nv):
            edges.append(DirectedEdge(v, (v+1)%nv, float(rng.uniform(0.2,2.0)), float(rng.normal())))
        extra = int(rng.integers(0, nv+1))
        for _k in range(extra):
            u, v = rng.choice(nv, 2, replace=False)
            edges.append(DirectedEdge(int(u), int(v), float(rng.uniform(0.2,2.0)), float(rng.normal())))
        D = incidence(nv, edges)
        # project random vectors to kernel(D)
        P = np.eye(len(edges)) - D.T @ np.linalg.pinv(D @ D.T, rcond=1e-12) @ D
        cycles=[]; masses=[]
        for _j in range(4):
            z=P @ rng.normal(size=len(edges))
            if np.linalg.norm(z)<1e-10:
                continue
            cycles.append(z)
            masses.append(float(rng.uniform(0.1,1.0)))
        if not cycles:
            continue
        cert=multicommodity_rayleigh(nv,edges,cycles,np.array(masses))
        worst=min(worst,cert["margin"])
        checked+=1
    return {"checked":checked,"worst_margin":float(worst)}


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument("--samples",type=int,default=20000)
    p.add_argument("--outdir",type=Path,default=Path("results-multicommodity-hodge"))
    args=p.parse_args(); args.outdir.mkdir(parents=True,exist_ok=True)
    rstar=0.6109041018281888
    gamma=math.log(1.0/rstar)
    random_checks=random_rayleigh_checks(args.samples)
    toys=_toy_graphs(gamma)
    tree=tree_pair_resistance_identity(5,[(0,1,1.0),(1,2,0.7),(1,3,2.0),(3,4,0.4)],np.array([0.3,0.25,0.2,0.15,0.1]))
    result={"gamma":gamma,"random_checks":random_checks,"toy":toys,"tree_identity":tree}
    (args.outdir/"multicommodity_hodge.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    lines=[
        "# Multicommodity Hodge routing",
        "",
        f"gamma*: `{gamma:.12f}`",
        f"Random Rayleigh checks: `{random_checks['checked']}`; worst energy-bound margin `{random_checks['worst_margin']:.3e}`.",
        "",
        "## Gauge synchronization toys",
    ]
    for name,cert in toys.items():
        lines.append(f"- {name}: Hodge `{cert['hodge_energy']:.9f}`, pair resistance `{cert['pair_resistance_budget']:.9f}`, shift variance `{cert['shift_variance']:.9f}`, sync margin `{cert['sync_margin']:.3e}`, integer mismatch bound `{cert['integer_mismatch_probability_bound']:.6f}`")
    lines += ["", "## Tree pair-resistance identity", f"Direct `{tree['direct']:.12f}`, edge-cut formula `{tree['formula']:.12f}`, difference `{tree['margin']:.3e}`.", "", "The multicommodity Rayleigh inequality, exact-gradient synchronization identity, integer gauge mismatch consequence, and tree resistance identity are exact finite-dimensional statements. Random tests only check implementation."]
    (args.outdir/"summary.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines))

if __name__ == "__main__":
    main()
