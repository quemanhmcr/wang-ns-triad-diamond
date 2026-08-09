from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.first_hit_heat_reservoir_erosion import (
    CLEAN_HEAT_POOL_RATIO,
    PHYSICAL_HEAT_POOL_RATIO_UPPER,
    old_pool_heat_capacity_upper,
)
from src.high_strain_heat_increment_service import (
    critical_ancestor_heat_service_fraction_lower,
    high_strain_heat_service_lower,
)


def band_addressed_material_partition(
    edge_weights: Sequence[float],
    old_shell: Sequence[bool],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
) -> dict[str, float]:
    """Positive OO/ON/NN partition with the exact eventwise shell mark retained.

    Heat service is first disintegrated by deterministic Fourier shell j because
    delta_r commutes exactly with P_j.  Moyal then disintegrates each shell into
    coherent spatial edges.  The material old pool at that event is represented
    shellwise: an endpoint may be old only on a shell belonging to the current
    transported old-frequency envelope.

    The shell mark is *not* promoted to a persistent material label.  It is an
    exact event-fiber provenance mark, while intrinsic zeta remains the material
    identity transported between events.
    """
    w = np.asarray(edge_weights, float)
    os = np.asarray(old_shell, bool)
    a = np.asarray(old_here, bool)
    b = np.asarray(old_neighbor, bool)
    if w.ndim != 1 or os.shape != w.shape or a.shape != w.shape or b.shape != w.shape:
        raise ValueError("matching one-dimensional edge/shell/ownership data required")
    if np.any(~np.isfinite(w)) or np.any(w < 0):
        raise ValueError("finite nonnegative edge weights required")
    if np.any((a | b) & (~os)):
        raise ValueError("an old material endpoint must lie on a currently old-addressed shell")

    oo = a & b
    on = np.logical_xor(a, b)
    nn = (~a) & (~b)
    old_incident = a | b

    total = float(w.sum())
    old_shell_service = float(w[os].sum())
    old_old = float(w[oo].sum())
    interface = float(w[on].sum())
    new_new = float(w[nn].sum())
    incident = float(w[old_incident].sum())
    return {
        "total": total,
        "old_shell_service": old_shell_service,
        "old_old": old_old,
        "old_new_interface": interface,
        "new_new": new_new,
        "old_incident": incident,
        "ownership_partition_residual": old_old + interface + new_new - total,
        "incident_identity_residual": old_old + interface - incident,
        "old_shell_capacity_margin": old_shell_service - incident,
    }


def old_incident_heat_capacity_upper(
    *,
    generation: int,
    initial_low_cut_ratio: float,
    initial_block_frequency: float,
    frame_energy_bound: float,
    global_energy: float,
    scaled_lifetime: float = 1.0,
) -> float:
    """Whole old-incident heat capacity on a supplied signed-good epoch.

    Every OO or ON edge touches an old endpoint, hence its deterministic shell
    belongs to the transported old shell envelope.  Its positive mass is thus
    bounded by the *whole* heat service on those old shells.  Orthogonal shell
    energy and the global frame budget give the same c M_old^2 E/N capacity as
    the preceding OO theorem, now for OO+ON together.
    """
    return old_pool_heat_capacity_upper(
        generation=generation,
        initial_low_cut_ratio=initial_low_cut_ratio,
        initial_block_frequency=initial_block_frequency,
        frame_energy_bound=frame_energy_bound,
        global_energy=global_energy,
        scaled_lifetime=scaled_lifetime,
    )


def forced_nn_service_lower(total_heat_service: float, old_incident_capacity: float) -> float:
    """If OO+ON <= C_inc, exact ownership forces NN >= S_heat-C_inc."""
    S = float(total_heat_service)
    C = float(old_incident_capacity)
    if S < 0 or C < 0 or not math.isfinite(S + C):
        raise ValueError("finite nonnegative heat service/capacity required")
    return max(0.0, S - C)


def positive_measure_intersection_lower(total: float, first_mass: float, second_mass: float) -> float:
    """Sharp inclusion-exclusion lower bound for two submeasures of one positive law."""
    T = float(total)
    A = float(first_mass)
    B = float(second_mass)
    if T < 0 or A < 0 or B < 0 or A > T or B > T or not math.isfinite(T + A + B):
        raise ValueError("finite submeasure masses in [0,total] required")
    return max(0.0, A + B - T)


def canonical_nn_critical_thresholds() -> dict[str, float]:
    """Choose epsilon=g/2 so old incidence smallness forces a critical NN overlap.

    The existing high-strain heat theorem gives S_G >= g S_total on the
    critical-shell-time set G, where g=e^(-1/32)/2.  If C_inc<=epsilon S_* and
    S_total>=S_*, then S_NN >= (1-epsilon)S_total.  Inclusion-exclusion yields
    S_(NN intersect G) >= (g-epsilon)S_total.  Taking epsilon=g/2 leaves g/2.
    """
    g = float(critical_ancestor_heat_service_fraction_lower())
    eps = 0.5 * g
    return {
        "critical_heat_fraction": g,
        "old_incident_fraction_target": eps,
        "new_new_fraction_lower": 1.0 - eps,
        "nn_critical_intersection_fraction_lower": g - eps,
    }


def first_nn_critical_generation(
    *,
    scaled_lifetime: float,
    initial_old_capacity: float,
) -> int:
    """First q with C_inc(q)<=epsilon S_*(c), epsilon=g/2."""
    c = float(scaled_lifetime)
    C0 = float(initial_old_capacity)
    if c <= 0 or C0 < 0 or not math.isfinite(c + C0):
        raise ValueError("positive finite scaled lifetime and nonnegative capacity required")
    eps = canonical_nn_critical_thresholds()["old_incident_fraction_target"]
    target = eps * high_strain_heat_service_lower(c)
    if C0 <= target:
        return 0
    r = float(CLEAN_HEAT_POOL_RATIO)
    q = max(0, int(math.ceil(math.log(target / C0) / math.log(r))))
    while C0 * r**q > target:
        q += 1
    while q > 0 and C0 * r ** (q - 1) <= target:
        q -= 1
    return q


def theorem_certificate(scaled_lifetime: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    th = canonical_nn_critical_thresholds()
    g = th["critical_heat_fraction"]
    eps = th["old_incident_fraction_target"]
    if not (0 < eps < g < 0.5):
        raise AssertionError("canonical NN/critical overlap thresholds lost their ordering")
    if not (PHYSICAL_HEAT_POOL_RATIO_UPPER < float(CLEAN_HEAT_POOL_RATIO) < 0.7):
        raise AssertionError("old-shell heat contraction changed")
    return {
        "status": "EXACT_OLD_INCIDENT_HEAT_EROSION__OO_PLUS_ON_DECAY__POSITIVE_NN_CRITICAL_HEAT_SEED_ON_SUPPLIED_SIGNED_GOOD_EPOCH",
        "shell_fiber": "delta_r P_j=P_j delta_r exactly; deterministic shell j is retained as an eventwise heat-law mark before Moyal and is not promoted to a persistent material identity",
        "material_fiber": "capacity routing specializes the arbitrary-O ownership theorem to the canonical orthogonal-band-addressed old reservoir pool already used by Moyal/reservoir erosion; within each event shell an old endpoint can occur only inside the current old-frequency envelope",
        "old_incident": "pointwise old_incident=OO+ON and old_incident is a positive submeasure of the whole old-shell heat law",
        "capacity": "S_OO+S_ON <= S_old-shell <= c M_old^2 P E_global/N",
        "erosion": f"on a supplied signed-good first-hit epoch rho_phys<={PHYSICAL_HEAT_POOL_RATIO_UPPER:.12g}<441/640<7/10",
        "critical_heat": f"the existing high-strain law puts at least g=e^(-1/32)/2={g:.12g} of total heat service on critical shell-time ancestors",
        "canonical_stop": f"once C_old-incident<=epsilon S_* with epsilon=g/2={eps:.12g}, NN carries at least {th['new_new_fraction_lower']:.12g} of total heat service and NN intersect critical carries at least g/2={th['nn_critical_intersection_fraction_lower']:.12g}",
        "interpretation": "same-time ON heat service need not be invented into a temporal relink charge: because the edge keeps its exact shell, any edge touching old material remains old-shell incident and erodes with the old frequency envelope",
        "scope": "for the canonical band-addressed old reservoir pool, this produces a positive epoch-new NN heat sublaw simultaneously carrying critical resolved-shell ancestry on every sufficiently old supplied signed-good epoch; it does not give a per-cell mass floor, a selected transfer parent, or universal epoch/slab renewal",
    }


@dataclass(frozen=True)
class OldIncidentHeatStress:
    samples: int
    worst_shell_increment_commutator_residual: float
    worst_ownership_partition_residual: float
    worst_old_incident_identity_residual: float
    minimum_old_shell_capacity_margin: float
    maximum_physical_epoch_ratio: float
    minimum_clean_ratio_margin: float
    minimum_nn_fraction_margin: float
    minimum_nn_critical_intersection_margin: float
    minimum_stopping_generation_margin: float


def _periodic_shell_project(f: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return np.fft.ifft(np.fft.fft(f) * mask)


def stress(samples: int = 50_000, seed: int = 20260809) -> OldIncidentHeatStress:
    rng = np.random.default_rng(seed)
    wc = wp = wi = 0.0
    mcap = mratio = mnn = mint = mstop = float("inf")
    maxepoch = 0.0
    rclean = float(CLEAN_HEAT_POOL_RATIO)
    th = canonical_nn_critical_thresholds()
    g = th["critical_heat_fraction"]
    eps = th["old_incident_fraction_target"]

    for i in range(samples):
        # Exact shell preservation by a translation increment, checked on a
        # representative subset because the analytic identity is multiplier
        # commutation.
        if i < min(samples, 3000):
            n = int(rng.integers(8, 80))
            f = rng.normal(size=n) + 1j * rng.normal(size=n)
            mask = rng.random(n) < 0.45
            shift = int(rng.integers(-n // 3, n // 3 + 1))
            Pf = _periodic_shell_project(f, mask)
            lhs = np.roll(Pf, shift) - Pf
            df = np.roll(f, shift) - f
            rhs = _periodic_shell_project(df, mask)
            res = float(np.linalg.norm(lhs - rhs) / max(1.0, np.linalg.norm(f)))
            wc = max(wc, res)
            if res > 3e-11:
                raise AssertionError("translation increment changed deterministic Fourier shell provenance")

        # Positive event-fiber ownership: old endpoints are allowed only on old shells.
        nedge = int(rng.integers(2, 150))
        w = rng.lognormal(mean=-1.0, sigma=1.3, size=nedge)
        old_shell = rng.random(nedge) < 0.55
        a = old_shell & (rng.random(nedge) < 0.5)
        b = old_shell & (rng.random(nedge) < 0.5)
        part = band_addressed_material_partition(w, old_shell, a, b)
        scale = max(1.0, float(part["total"]))
        wp = max(wp, abs(float(part["ownership_partition_residual"])) / scale)
        wi = max(wi, abs(float(part["incident_identity_residual"])) / scale)
        mcap = min(mcap, float(part["old_shell_capacity_margin"]))
        if abs(float(part["ownership_partition_residual"])) > 3e-12 * scale:
            raise AssertionError("OO/ON/NN positive partition lost heat mass")
        if abs(float(part["incident_identity_residual"])) > 3e-12 * scale:
            raise AssertionError("old incident heat service is not OO+ON")
        if float(part["old_shell_capacity_margin"]) < -3e-12 * scale:
            raise AssertionError("old-incident service escaped whole old-shell heat service")

        # First-hit/signed-good physical coefficient is strictly below the clean envelope.
        mratio_step = math.exp(float(rng.uniform(0.0, 1.0 / 30.0)))
        nratio_step = float(rng.uniform(8.0 / 5.0 * (1.0 + 1e-12), 4.0))
        eratio = mratio_step * mratio_step / nratio_step
        maxepoch = max(maxepoch, eratio)
        mratio = min(mratio, rclean - eratio)
        if eratio >= rclean + 2e-13:
            raise AssertionError("old-incident heat coefficient escaped clean 441/640 envelope")

        # Once old incidence <=epsilon*S*, NN occupies >=1-epsilon of total.
        c = float(math.exp(rng.uniform(-2.0, 2.0)))
        Sstar = high_strain_heat_service_lower(c)
        total = float(rng.uniform(1.0, 4.0)) * Sstar
        oldinc = float(rng.uniform(0.0, 1.0)) * eps * Sstar
        nn = forced_nn_service_lower(total, oldinc)
        nn_margin = nn / total - (1.0 - eps)
        mnn = min(mnn, nn_margin)
        if nn_margin < -3e-13:
            raise AssertionError("forced NN fraction fell below 1-epsilon")

        good = float(rng.uniform(g, 1.0)) * total
        overlap_lower = positive_measure_intersection_lower(total, nn, good)
        target_overlap = (g - eps) * total
        imargin = overlap_lower - target_overlap
        mint = min(mint, imargin)
        if imargin < -4e-13 * max(1.0, total):
            raise AssertionError("NN/critical positive-measure intersection lower failed")

        # Finite stopping generation with clean geometric old-incident capacity.
        C0 = float(math.exp(rng.uniform(-5.0, 6.0)))
        q = first_nn_critical_generation(scaled_lifetime=c, initial_old_capacity=C0)
        target = eps * Sstar
        Cq = C0 * rclean**q
        smargin = target - Cq
        mstop = min(mstop, smargin)
        if smargin < -4e-13 * max(1.0, target):
            raise AssertionError("NN-critical stopping generation still has too much old incidence")
        if q > 0 and C0 * rclean ** (q - 1) <= target - 4e-13 * max(1.0, target):
            raise AssertionError("NN-critical stopping generation was not minimal")

    return OldIncidentHeatStress(
        samples,
        wc,
        wp,
        wi,
        mcap,
        maxepoch,
        mratio,
        mnn,
        mint,
        mstop,
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-old-incident-heat-erosion"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate()
    out = stress(args.samples)
    (args.outdir / "old_incident_heat_erosion.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    th = canonical_nn_critical_thresholds()
    md = f"""# Old-incident heat erosion: mixed ON service cannot form a neutral third regime\n\nStatus: **{cert['status']}**.\n\nThe high-strain heat law already carries an exact deterministic Fourier shell mark before coherent localization.  Because translations and shell multipliers commute,\n\n`delta_r P_j = P_j delta_r`.\n\nMoyal is then applied **inside each shell**.  Keep this shell index only as an event-fiber provenance mark; intrinsic `zeta` remains the material identity transported between events.  Represent the current old reservoir pool shellwise by its transported intrinsic-zeta cells.\n\nAn OO or ON edge has at least one old endpoint.  By definition of the current old pool that endpoint lies on a shell inside the transported old-frequency envelope.  Since the heat edge remains in the same exact shell, pointwise\n\n`old-incident = OO + ON`,\n\nand as positive measures\n\n`S_OO+S_ON <= S_old-shell`.\n\nThe preceding first-hit heat theorem bounds the whole old-shell law on `T(N)=cN^-2` by\n\n`S_old-shell <= c M_old^2 P E_global/N`.\n\nHence **OO and ON together** inherit the same physical contraction\n\n`rho_phys <= (5/8)exp(1/15) = {PHYSICAL_HEAT_POOL_RATIO_UPPER:.12g} < 441/640 < 7/10`.\n\nThere is therefore no need to manufacture a temporal relink charge merely because a same-time heat edge is mixed.  ON service remains physically real, but any heat service touching the old material pool is still old-frequency incident and its total capacity decays geometrically.\n\nThe high-strain heat theorem also supplies a second restriction of the same positive law: at least\n\n`g=e^(-1/32)/2 = {th['critical_heat_fraction']:.12g}`\n\nof total heat service lies on shell-time atoms with critical resolved mass.  Choose\n\n`epsilon=g/2 = {th['old_incident_fraction_target']:.12g}`.\n\nAfter the first material age for which `C_old-incident<=epsilon S_*`, every high-strain event in that supplied epoch obeys\n\n`S_NN >= (1-epsilon) S_heat >= {th['new_new_fraction_lower']:.12g} S_heat`.\n\nSince the critical-ancestor set carries at least `g S_heat`, sharp inclusion--exclusion of two submeasures of the **same** heat law yields\n\n`S_(NN intersect critical) >= (g-epsilon)S_heat = {th['nn_critical_intersection_fraction_lower']:.12g} S_heat`.\n\nThus sufficiently old high-strain events contain a fixed positive **epoch-new NN coherent heat sublaw simultaneously marked by a critical lower-frequency resolved ancestor**.  No independence assumption, packet argmax or per-cell mass floor appears.\n\nStress: `{out.samples}` shell/ownership/epoch/intersection states\n- worst shell-increment commutator residual: `{out.worst_shell_increment_commutator_residual:.3e}`\n- worst OO/ON/NN partition residual: `{out.worst_ownership_partition_residual:.3e}`\n- worst old-incident identity residual: `{out.worst_old_incident_identity_residual:.3e}`\n- minimum old-shell capacity margin: `{out.minimum_old_shell_capacity_margin:.3e}`\n- maximum sampled physical epoch ratio: `{out.maximum_physical_epoch_ratio:.9f}`\n- minimum clean `441/640` margin: `{out.minimum_clean_ratio_margin:.3e}`\n- minimum forced NN-fraction margin: `{out.minimum_nn_fraction_margin:.3e}`\n- minimum NN/critical intersection margin: `{out.minimum_nn_critical_intersection_margin:.3e}`\n- minimum stopping-generation margin: `{out.minimum_stopping_generation_margin:.3e}`\n\nThis removes ON as a potentially neutral long-lived high-strain regime on a supplied signed-good material epoch.  The next seam is sharper: turn the positive NN-critical heat sublaw into a renewed smooth material carrier/epoch, or hit an already named source/relink/interface cause, without asserting a fixed cell mass.  Universal slab renewal is still open, and no global-regularity claim is made.\n"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
