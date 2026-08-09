from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.high_strain_resolved_ancestor import high_strain_ancestor_mass_threshold
from src.old_incident_heat_erosion import canonical_nn_critical_thresholds

RENEWAL_SCALE_FACTOR = 3.0 / 4.0
HARD_SHELL_LOWER_RELATIVE_TO_RENEWAL = 2.0 / 3.0
HARD_SHELL_UPPER_RELATIVE_TO_RENEWAL = 4.0 / 3.0
SMOOTH_ENVELOPE_LOWER = 3.0 / 5.0
SMOOTH_ENVELOPE_UPPER = 3.0 / 2.0
TRANSPORT_CUT = 1.0 / 4.0
LOW_STRAIN_ACTION = 1.0 / 30.0


def renewal_scale(shell_upper_frequency: float) -> float:
    M = float(shell_upper_frequency)
    if M <= 0 or not math.isfinite(M):
        raise ValueError("positive finite shell frequency required")
    return RENEWAL_SCALE_FACTOR * M


def shell_relative_support() -> tuple[float, float]:
    """A_j={M/2<|xi|<=M} becomes {2A/3<|xi|<=4A/3} for A=3M/4."""
    return HARD_SHELL_LOWER_RELATIVE_TO_RENEWAL, HARD_SHELL_UPPER_RELATIVE_TO_RENEWAL


def renewal_carrier_critical_mass_lower(scaled_lifetime: float) -> float:
    """A||P_j u||_2^2 >= (3/4)mu_* on every critical shell-time atom."""
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    return RENEWAL_SCALE_FACTOR * high_strain_ancestor_mass_threshold(c)


def renewal_natural_lifetime_ratio(child_frequency: float, shell_upper_frequency: float) -> float:
    """T_A/T_N=(N/A)^2; M<=N/4 gives the clean lower 256/9."""
    N = float(child_frequency)
    M = float(shell_upper_frequency)
    if N <= 0 or M <= 0 or M > N / 4.0 * (1.0 + 1e-13) or not math.isfinite(N + M):
        raise ValueError("require positive resolved shell M<=N/4")
    A = renewal_scale(M)
    return (N / A) ** 2


def persistent_seed_low_low_gap(
    strain_action: float = LOW_STRAIN_ACTION,
    envelope_lower: float = SMOOTH_ENVELOPE_LOWER,
    transport_cut: float = TRANSPORT_CUT,
) -> float:
    K = float(strain_action)
    lo = float(envelope_lower)
    cut = float(transport_cut)
    if K < 0 or lo <= 0 or cut < 0 or not all(math.isfinite(x) for x in (K, lo, cut)):
        raise ValueError("valid finite support data required")
    return lo * math.exp(-K) - 2.0 * cut


def normalized_shell_probe_coefficient(shell_state: np.ndarray, smooth_envelope: np.ndarray) -> dict[str, float]:
    """Canonical shell direction registers exactly into any envelope equal to one on it.

    In a finite spectral model ``shell_state`` already lies in the hard shell and
    ``smooth_envelope`` is the diagonal scalar value of Q on those coordinates.
    The continuum statement is: with f=P_j u, psi=f/||f|| and QP_j=P_j,

        <Q u,psi>=<P_j u,psi>=||P_j u||_2.

    This chooses no spatial/coherent packet: the dual direction is the shell's
    own normalized physical state.
    """
    f = np.asarray(shell_state, complex)
    q = np.asarray(smooth_envelope, float)
    if f.ndim != 1 or q.shape != f.shape or np.any(~np.isfinite(f)) or np.any(~np.isfinite(q)):
        raise ValueError("matching finite one-dimensional shell/envelope data required")
    if np.any(np.abs(q - 1.0) > 2e-12):
        raise ValueError("envelope must equal one on the represented hard shell")
    norm = float(np.linalg.norm(f))
    if norm <= 0:
        raise ValueError("nonzero shell state required")
    psi = f / norm
    z = complex(np.vdot(psi, q * f))
    return {
        "shell_energy": norm * norm,
        "coefficient_abs": abs(z),
        "coefficient_energy": abs(z) ** 2,
        "registration_residual": abs(z - norm),
    }


@dataclass(frozen=True)
class NNCriticalHeatAtom:
    mass: float
    child_frequency: float
    shell_upper_frequency: float
    shell_energy_u: float
    time: float
    nn_endpoint_mark: str = "new_new_heat_edge"


@dataclass(frozen=True)
class SmoothCarrierSeed:
    probability: float
    heat_mass: float
    time: float
    child_frequency: float
    shell_upper_frequency: float
    renewal_frequency: float
    shell_critical_mass: float
    renewal_critical_mass: float
    natural_lifetime_ratio: float
    nn_endpoint_mark: str


def pushforward_nn_critical_heat_law(
    atoms: Sequence[NNCriticalHeatAtom],
    *,
    scaled_lifetime: float,
) -> tuple[SmoothCarrierSeed, ...]:
    """Push the positive NN-critical heat law to lower-scale smooth carrier seeds.

    Input atoms are already restricted to the canonical `NN intersect G` heat
    sublaw.  Their heat masses remain the law weights.  The map forgets the
    coherent edge coordinates only after retaining their NN provenance mark and
    sends `(j,t)` to the whole hard shell plus a smooth Fourier plateau envelope.

    No lower bound on one heat atom or one coherent cell is required.  Critical
    mass comes from membership in G: M||P_j u(t)||^2>=mu_*.
    """
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    rows = tuple(atoms)
    if not rows:
        raise ValueError("nonempty NN-critical positive heat law required")
    mu_star = high_strain_ancestor_mass_threshold(c)
    total = 0.0
    for a in rows:
        vals = (a.mass, a.child_frequency, a.shell_upper_frequency, a.shell_energy_u, a.time)
        if a.mass <= 0 or a.child_frequency <= 0 or a.shell_upper_frequency <= 0 or a.shell_energy_u < 0:
            raise ValueError("positive heat mass/frequencies and nonnegative shell energy required")
        if not all(math.isfinite(x) for x in vals):
            raise ValueError("finite NN-critical heat atom data required")
        if a.shell_upper_frequency > a.child_frequency / 4.0 * (1.0 + 1e-13):
            raise ValueError("heat atom shell must lie in the resolved ball")
        if a.shell_upper_frequency * a.shell_energy_u < mu_star - 2e-13 * max(1.0, mu_star):
            raise ValueError("atom is not on the critical shell-time set G")
        total += a.mass
    seeds: list[SmoothCarrierSeed] = []
    for a in rows:
        M = a.shell_upper_frequency
        A = renewal_scale(M)
        mu = M * a.shell_energy_u
        seeds.append(
            SmoothCarrierSeed(
                probability=a.mass / total,
                heat_mass=a.mass,
                time=a.time,
                child_frequency=a.child_frequency,
                shell_upper_frequency=M,
                renewal_frequency=A,
                shell_critical_mass=mu,
                renewal_critical_mass=A * a.shell_energy_u,
                natural_lifetime_ratio=renewal_natural_lifetime_ratio(a.child_frequency, M),
                nn_endpoint_mark=a.nn_endpoint_mark,
            )
        )
    if not math.isclose(sum(s.probability for s in seeds), 1.0, rel_tol=2e-14, abs_tol=2e-14):
        raise AssertionError("NN-critical heat-law pushforward failed to normalize")
    lower = renewal_carrier_critical_mass_lower(c)
    for s in seeds:
        if s.renewal_critical_mass < lower - 3e-13 * max(1.0, lower):
            raise AssertionError("critical heat atom lost its lower-scale carrier mass")
    return tuple(seeds)


def theorem_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    lower = renewal_carrier_critical_mass_lower(c)
    frac = canonical_nn_critical_thresholds()["nn_critical_intersection_fraction_lower"]
    lo, hi = shell_relative_support()
    gap = persistent_seed_low_low_gap()
    if abs(lower - 8.0 * math.pi * math.pi / (25.0 * c * c)) > 2e-13 * max(1.0, lower):
        raise AssertionError("clean renewed critical mass identity changed")
    if gap <= 0:
        raise AssertionError("smooth shell carrier seed lost low-low support separation")
    return {
        "status": "EXACT_NN_CRITICAL_HEAT_LAW_TO_SMOOTH_SHELL_CARRIER_SEEDS__NO_CELL_MASS_FLOOR__TEMPORAL_RENEWAL_REMAINS",
        "input_law": f"on every sufficiently old supplied epoch, the existing positive NN-intersect-critical heat sublaw has mass at least (1/4)e^(-1/32)={frac:.12g} of S_heat",
        "selection": "normalize that physical heat sublaw itself; no shell argmax, coherent-cell maximizer or uniform atom mass is introduced",
        "scale": "for A=3M/4, the hard shell {M/2<|xi|<=M} becomes {2A/3<|xi|<=4A/3}",
        "critical_coefficient": f"G gives M||P_j u||^2>=mu_*; hence A||P_j u||^2>=(3/4)mu_*=8pi^2/(25c^2)={lower:.12g} at c={c:.12g}",
        "smooth_envelope": f"choose scalar Q_A=1 on the hard shell with support lower >=3A/5 and upper <=3A/2; the hard shell [{lo:.12g},{hi:.12g}]A lies strictly inside this plateau envelope",
        "dual_probe": "with psi=P_j u/||P_j u||, <Q_A u,psi>=||P_j u|| exactly; the smooth carrier seed inherits the whole-shell critical coefficient without choosing a packet",
        "low_low": f"if the renewed A-scale strain action stays <=1/30, the envelope lower edge remains (3/5)e^(-1/30)A>A/2; dimensionless gap={gap:.12g}",
        "lifetime": "M<=N/4 implies A<=3N/16 and T_A/T_N=(N/A)^2>=256/9",
        "material_provenance": "each seed retains the NN coherent heat-edge mark that selected its shell-time atom, but the theorem does not claim the whole u-shell is new material; V-shell heat provenance and u-shell critical energy remain distinct exact marks",
        "scope": "this supplies a positive law of lower-scale smooth Fourier carrier seeds with critical coefficients and support moat; it does not yet extend a seed through a full A-natural slab, identify its whole energy as NN material, or prove universal slab renewal",
    }


@dataclass(frozen=True)
class CarrierSeedStress:
    samples: int
    minimum_renewal_critical_mass_margin: float
    minimum_natural_lifetime_ratio: float
    minimum_low_low_gap: float
    worst_shell_probe_registration_residual: float
    worst_probability_residual: float
    minimum_heat_mass_conservation_margin: float


def stress(samples: int = 50_000, seed: int = 20260809) -> CarrierSeedStress:
    rng = np.random.default_rng(seed)
    mm = ml = mg = float("inf")
    wr = wp = 0.0
    mh = float("inf")
    for _ in range(samples):
        c = float(math.exp(rng.uniform(-2.0, 2.0)))
        N = float(math.exp(rng.uniform(0.0, 8.0)))
        mu_star = high_strain_ancestor_mass_threshold(c)
        n = int(rng.integers(1, 20))
        atoms: list[NNCriticalHeatAtom] = []
        total = 0.0
        for j in range(n):
            shell_index = int(rng.integers(0, 15))
            M = (N / 4.0) * 2.0 ** (-shell_index)
            E = float(rng.uniform(1.0, 5.0)) * mu_star / M
            mass = float(math.exp(rng.uniform(-8.0, 4.0)))
            total += mass
            atoms.append(NNCriticalHeatAtom(mass, N, M, E, float(rng.uniform(0.0, c / (N*N)))))
        seeds = pushforward_nn_critical_heat_law(atoms, scaled_lifetime=c)
        lower = renewal_carrier_critical_mass_lower(c)
        mm = min(mm, min(s.renewal_critical_mass - lower for s in seeds))
        ml = min(ml, min(s.natural_lifetime_ratio for s in seeds))
        pres = abs(sum(s.probability for s in seeds) - 1.0)
        wp = max(wp, pres)
        pushed_mass = sum(s.heat_mass for s in seeds)
        mh = min(mh, total - pushed_mass)
        if abs(pushed_mass - total) > 3e-13 * max(1.0, total):
            raise AssertionError("carrier-seed pushforward changed NN-critical heat mass")
        if min(s.renewal_critical_mass for s in seeds) < lower - 3e-12 * max(1.0, lower):
            raise AssertionError("renewed shell carrier lost critical mass")
        if min(s.natural_lifetime_ratio for s in seeds) < 256.0 / 9.0 - 2e-12:
            raise AssertionError("renewed carrier lifetime separation fell below 256/9")

        K = float(rng.uniform(0.0, LOW_STRAIN_ACTION))
        gap = persistent_seed_low_low_gap(K)
        mg = min(mg, gap)
        if gap <= 0:
            raise AssertionError("renewed smooth seed met low-low output")

        dim = int(rng.integers(2, 40))
        f = rng.normal(size=dim) + 1j * rng.normal(size=dim)
        q = np.ones(dim, float)
        reg = normalized_shell_probe_coefficient(f, q)
        rr = float(reg["registration_residual"]) / max(1.0, math.sqrt(float(reg["shell_energy"])))
        wr = max(wr, rr)
        if rr > 3e-12:
            raise AssertionError("canonical shell direction did not register into smooth envelope")
    return CarrierSeedStress(samples, mm, ml, mg, wr, wp, mh)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-nn-critical-heat-carrier-seed"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "nn_critical_heat_carrier_seed.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# NN-critical heat law becomes a law of smooth lower-scale carrier seeds\n\nStatus: **{cert['status']}**.\n\nThe preceding theorem leaves a positive high-strain heat sublaw which is simultaneously `NN` in material endpoint provenance and on the critical shell-time set `G`.  Do **not** select a largest coherent cell.  Normalize this positive heat sublaw itself and use it as the physical selector.  Every atom already carries a deterministic shell `A_j={{M/2<|xi|<=M}}` and a physical time `t`.\n\nOn `G`,\n\n`M ||P_j u(t)||_2^2 >= mu_* = 32 pi^2/(75 c^2)`.\n\nDefine the renewed frequency\n\n`A=3M/4`.\n\nThen the hard shell becomes\n\n`2A/3 < |xi| <= 4A/3`,\n\nand its critical mass at the renewed scale is\n\n`A ||P_j u(t)||_2^2 >= (3/4)mu_* = 8 pi^2/(25 c^2)`.\n\nChoose a smooth scalar Fourier envelope `Q_A` equal to one on this whole hard shell, with lower support `>=3A/5` and upper support `<=3A/2`.  There is no packet synthesis.  The canonical terminal dual direction is simply\n\n`psi=P_j u(t)/||P_j u(t)||_2`,\n\nso exactly\n\n`<Q_A u(t),psi>=||P_j u(t)||_2`.\n\nThus every atom of the NN-critical heat law gives a **smooth whole-shell carrier seed with a scale-independent critical coefficient**, without asserting any positive mass for one coherent cell.  The atom's NN coherent endpoints are retained as material provenance, but are not used to falsely declare the entire `u` shell new material; the `V`-heat mark and `u`-shell critical mass stay as simultaneous, distinct physical marks.\n\nThe support geometry is already compatible with the outer moving-role PDE.  If the renewed `A`-scale strain action remains `<=1/30`,\n\n`(3/5)e^(-1/30) A > A/2`,\n\nwhile `S_(A/4)u tensor S_(A/4)u` lies below `A/2`.  Hence the seed has a strict low--low moat.  Also `M<=N/4` implies `A<=3N/16`, so\n\n`T_A/T_N=(N/A)^2 >= 256/9`.\n\nPushforward of the normalized NN-critical heat law preserves its full positive mass; coalescing many edge atoms onto the same `(j,t)` seed does not create or destroy physical weight.  The existing theorem supplies at least `(1/4)e^(-1/32)` of high-strain heat service to this input law.\n\nStress: `{out.samples}` critical-shell/law/support/probe states\n- minimum renewed critical-mass margin: `{out.minimum_renewal_critical_mass_margin:.3e}`\n- minimum renewed/child natural-lifetime ratio: `{out.minimum_natural_lifetime_ratio:.9f}`\n- minimum low-low support gap: `{out.minimum_low_low_gap:.6e}`\n- worst canonical shell-probe registration residual: `{out.worst_shell_probe_registration_residual:.3e}`\n- worst seed-law probability residual: `{out.worst_probability_residual:.3e}`\n- minimum heat-mass conservation margin: `{out.minimum_heat_mass_conservation_margin:.3e}`\n\nThis closes the **eventwise carrier-seed selection** of the high-strain route.  What remains is temporal/material renewal: extend almost every such smooth seed through a full `A`-natural slab while preserving the appropriate material provenance, or hit an already named strain/source/relink/HH-generation/boundary first stop.  No full-shell `NN` claim, packet mass floor, or Navier--Stokes global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
