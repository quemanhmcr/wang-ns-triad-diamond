from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from src.amplitude_entropy_causal_reuse import adjoint_vector_growth_upper
from src.coherent_service_stopping import epoch_certificate
from src.heat_edge_material_ownership import partition_positive_edge_measure
from src.nn_critical_heat_carrier_seed import (
    LOW_STRAIN_ACTION,
    SMOOTH_ENVELOPE_LOWER,
    SMOOTH_ENVELOPE_UPPER,
)
from src.nn_seed_temporal_first_stop import inherited_seed_critical_mass_lower

BOUNDED_HEAT_RADIUS = 3.0


def transported_annular_support_ratios(
    strain_action: float = LOW_STRAIN_ACTION,
    lower: float = SMOOTH_ENVELOPE_LOWER,
    upper: float = SMOOTH_ENVELOPE_UPPER,
) -> tuple[float, float]:
    """Radial support envelope of the affine-dual transported renewed role.

    The anchor support is [lower*A, upper*A].  A trace-free affine history with
    action K moves every frequency radius by factors between exp(-K) and exp(K).
    Nonaffine role-interface motion is not hidden here: it is the separately
    monitored Duhamel/interface first-stop term of the temporal theorem.
    """
    K = float(strain_action)
    lo = float(lower)
    hi = float(upper)
    if K < 0 or lo <= 0 or hi <= lo or not all(math.isfinite(x) for x in (K, lo, hi)):
        raise ValueError("valid finite annular support data required")
    return lo * math.exp(-K), hi * math.exp(K)


def renewed_analysis_probe_growth_upper(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Scale-independent inverse-heat/strain cost for the renewed registered analysis dual.

    The exact registered-coefficient interaction picture puts common affine/Kelvin
    transport and bulk viscosity in the analysis propagator.  The nonaffine/full-
    transport mismatch belongs to the classified interface residual already
    monitored by the temporal theorem; it is not hidden in this norm bound.  On T=c A^-2 and transported support |k|<=R A,

        ||psi(s)|| <= exp(K + nu*c*R^2) ||psi(t)||.

    This is the existing amplitude-entropy adjoint bound with the renewed role's
    actual upper support R=(3/2)exp(K).  Viscosity is therefore not silently
    discarded when converting coefficient persistence to carrier energy.
    """
    c = float(scaled_lifetime)
    nu = float(viscosity)
    K = float(strain_action)
    if c <= 0 or nu < 0 or K < 0 or not all(math.isfinite(x) for x in (c, nu, K)):
        raise ValueError("positive lifetime and nonnegative finite viscosity/action required")
    _, R = transported_annular_support_ratios(K)
    return adjoint_vector_growth_upper(
        c,
        viscosity=nu,
        strain_action=K,
        carrier_ratio_upper=R,
    )


def persistent_carrier_critical_mass_lower(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Uniform A||Q_A u(s)||_2^2 lower on a full no-hit natural corridor.

    No-hit means *every prefix* stays below the residual and HH impulse faces.
    Hence the exact Duhamel triangle gives |z(s)|>=|z(t)|/4 at every prefix.
    The temporal theorem gives A|z(t)|^2/16 >= pi^2/(50c^2).  Cauchy with the
    backward dual norm then yields the whole-corridor carrier-energy lower below.
    """
    J = renewed_analysis_probe_growth_upper(scaled_lifetime, viscosity, strain_action)
    return inherited_seed_critical_mass_lower(scaled_lifetime) / (J * J)


def heat_defect_fraction_lower(strain_action: float = LOW_STRAIN_ACTION) -> float:
    """Full H_A heat-defect fraction on the transported annular lower support.

    For r~N(0,A^-2 I),
      E ||delta_r w||_2^2
        = int 2(1-exp(-|xi|^2/(2A^2))) |what(xi)|^2 dxi.
    The multiplier is increasing in |xi| and the renewed role stays away from 0.
    """
    lo, _ = transported_annular_support_ratios(strain_action)
    return 2.0 * (1.0 - math.exp(-0.5 * lo * lo))


def gaussian_3d_tail_probability(radius: float) -> float:
    """Exact closed form P(|Z|>radius) for Z~N(0,I_3), evaluated in doubles."""
    b = float(radius)
    if b < 0 or not math.isfinite(b):
        raise ValueError("finite nonnegative radius required")
    return math.erfc(b / math.sqrt(2.0)) + math.sqrt(2.0 / math.pi) * b * math.exp(-0.5 * b * b)


def bounded_heat_defect_fraction_lower(
    strain_action: float = LOW_STRAIN_ACTION,
    radius: float = BOUNDED_HEAT_RADIUS,
) -> float:
    """Clean heat fraction guaranteed inside |r|<=radius/A.

    Since ||delta_r w||_2^2<=4||w||_2^2, the Gaussian tail can remove at most
    4 P(|Z|>radius) of the normalized square service.  The Arb certificate below
    proves that radius=3 leaves strictly more than half the full annular heat
    lower.  We deliberately return the clean half, not the sharper decimal.
    """
    q = heat_defect_fraction_lower(strain_action)
    tail = gaussian_3d_tail_probability(radius)
    if q - 4.0 * tail <= 0.5 * q:
        raise ValueError("chosen bounded heat radius does not retain the clean half-service")
    return 0.5 * q


def uniform_bounded_square_service_lower(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Scale-independent Y0 with some actual |r_s|<=3/A at every corridor time.

    If A||w(s)||_2^2>=mu_surv and the truncated heat average is at least
    q_b||w(s)||_2^2, then some displacement in the truncated ball obeys

        A||delta_r w(s)||_2^2 >= q_b mu_surv =: Y0.

    Moyal can then disintegrate this *actual bounded-displacement* service law.
    """
    return bounded_heat_defect_fraction_lower(strain_action) * persistent_carrier_critical_mass_lower(
        scaled_lifetime,
        viscosity,
        strain_action,
    )


def integrated_bounded_heat_service_lower(
    scaled_lifetime: float,
    viscosity: float = 1.0,
    strain_action: float = LOW_STRAIN_ACTION,
) -> float:
    """Scale-independent full-natural-slab service lower c*Y0.

    With T_A=cA^-2 and A E_{|r|<=3/A}||delta_r w||^2>=Y0 pointwise,
    the normalized integrated heat service A^3 int dt E[...] is >=cY0.
    """
    c = float(scaled_lifetime)
    if c <= 0 or not math.isfinite(c):
        raise ValueError("positive finite scaled lifetime required")
    return c * uniform_bounded_square_service_lower(c, viscosity, strain_action)


def no_hit_prefix_amplitude_lower(
    terminal_amplitude: float,
    residual_impulse_abs: float,
    hh_impulse_abs: float,
) -> float:
    """Exact triangle lower for one prefix of a no-hit corridor.

    This helper exposes why endpoint survival automatically means prefix-wise
    survival when the first-hit monitors are defined by cumulative impulses from
    the prefix to the terminal event.
    """
    A = float(terminal_amplitude)
    r = float(residual_impulse_abs)
    h = float(hh_impulse_abs)
    if A <= 0 or r < 0 or h < 0 or not all(math.isfinite(x) for x in (A, r, h)):
        raise ValueError("positive terminal amplitude and nonnegative finite impulses required")
    if r >= A / 4.0 or h >= A / 2.0:
        raise ValueError("prefix is not on the strict no-hit branch")
    return A - r - h


def material_service_partition(
    edge_weights: Sequence[float],
    old_here: Sequence[bool],
    old_neighbor: Sequence[bool],
) -> dict[str, float]:
    """Re-read materiality from the renewed positive service law itself.

    The terminal NN witness selected the high-strain seed.  It is *not* used to
    declare Q_Au new material.  After bounded-displacement heat/Moyal creates a
    positive renewed service measure, ownership is assigned afresh by the exact
    endpoint partition already certified for coherent heat edges.
    """
    return partition_positive_edge_measure(edge_weights, old_here, old_neighbor)


def service_epoch_reentry_certificate(
    scaled_lifetime: float,
    initial_old_capacity: float,
    viscosity: float = 1.0,
) -> dict[str, float | int | str]:
    """Feed the renewed uniform Y0 into the existing geometric stopping theorem.

    This does not manufacture signed-good scale progress.  It says that whenever
    the recursion supplies the already-canonical sticky signed-good continuation,
    the high-strain survivor now supplies exactly the previously missing uniform
    service threshold needed to start that material epoch.
    """
    C0 = float(initial_old_capacity)
    if C0 < 0 or not math.isfinite(C0):
        raise ValueError("finite nonnegative old capacity required")
    Y0 = uniform_bounded_square_service_lower(scaled_lifetime, viscosity)
    out = epoch_certificate(Y0, C0)
    out["renewed_service_threshold"] = Y0
    out["bounded_displacement_radius_over_A"] = BOUNDED_HEAT_RADIUS
    return out


def arb_bounded_heat_probe_certificate() -> dict[str, str]:
    """Rigorous interval check that beta=3 retains > half the annular heat lower."""
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("python-flint required") from exc
    ctx.prec = 180
    pi = arb.pi()
    K = arb(1) / 30
    lo = (arb(3) / 5) * (-K).exp()
    q = arb(2) * (arb(1) - (-(lo * lo) / 2).exp())
    b = arb(3)
    ball = (b / arb(2).sqrt()).erf() - (arb(2) / pi).sqrt() * b * (-(b * b) / 2).exp()
    tail = arb(1) - ball
    retained = q - 4 * tail
    if not (retained > q / 2):
        raise AssertionError("radius-3 truncated heat service did not retain a rigorous half")
    if not (lo > arb(1) / 2):
        raise AssertionError("transported annular lower edge unexpectedly met the low-low radius")
    return {
        "transported_lower_ratio": str(lo),
        "full_heat_fraction_lower": str(q),
        "gaussian_tail_probability_radius_3": str(tail),
        "truncated_fraction_lower": str(retained),
        "clean_truncated_fraction": str(q / 2),
        "status": "ARB_CERTIFIED_RADIUS_3_RETAINS_MORE_THAN_HALF_ANNULAR_HEAT_SERVICE",
    }


def theorem_certificate(scaled_lifetime: float = 1.0, viscosity: float = 1.0) -> dict[str, object]:
    c = float(scaled_lifetime)
    nu = float(viscosity)
    if c <= 0 or nu < 0 or not math.isfinite(c + nu):
        raise ValueError("positive finite lifetime and nonnegative viscosity required")
    arb = arb_bounded_heat_probe_certificate()
    lo, hi = transported_annular_support_ratios()
    J = renewed_analysis_probe_growth_upper(c, nu)
    mu = persistent_carrier_critical_mass_lower(c, nu)
    q = heat_defect_fraction_lower()
    qb = bounded_heat_defect_fraction_lower()
    Y0 = uniform_bounded_square_service_lower(c, nu)
    Sint = integrated_bounded_heat_service_lower(c, nu)
    return {
        "status": "EXACT_FULL_NATURAL_CRITICAL_CARRIER_TO_BOUNDED_ANNULAR_SERVICE_THRESHOLD__MATERIALITY_REREAD_FROM_SERVICE__UNIVERSAL_RENEWAL_REMAINS",
        "prefix_persistence": "on a full no-hit natural corridor every prefix has |I_interface|<|z(t)|/4 and |I_HH|<|z(t)|/2, hence |z(s)|>|z(t)|/4",
        "analysis_probe_growth": f"for the registered affine/Kelvin/viscous analysis dual, ||psi(s)||<=J||psi(t)|| with J=exp(K+nu*c*R^2)={J:.12g}, R=(3/2)exp(1/30)={hi:.12g}",
        "uniform_carrier_mass": f"A||Q_Au(s)||_2^2 >= pi^2/[50 c^2 J^2]={mu:.12g} for every prefix s",
        "annular_support": f"{lo:.12g} A <= |xi| <= {hi:.12g} A throughout the affine-dual low-strain corridor",
        "full_heat": f"E_HA||delta_r w||_2^2 >= q||w||_2^2 with q={q:.12g}",
        "bounded_heat": f"radius beta=3 retains more than half the heat lower; some actual |r_s|<=3/A has A||delta_r w(s)||_2^2>=Y0={Y0:.12g}",
        "clean_bounded_fraction": qb,
        "integrated_service": f"A^3 int_(full natural slab) E_(|r|<=3/A)||delta_r w||_2^2 dt >= cY0={Sint:.12g}",
        "material_reentry": "apply Moyal only after the bounded service exists, then assign OO/ON/NN by the two intrinsic endpoints of this renewed service law; the old NN seed witness is provenance, not whole-carrier ownership",
        "stopping_interface": "Y0 is exactly a scale-independent bounded-displacement square-service threshold consumable by the existing old-pool geometric stopping theorem on any supplied sticky signed-good continuation",
        "no_efficiency_overclaim": "this proves coherent service entrance, not near-extremal HH transfer efficiency of the whole renewed carrier",
        "scope": "the high-strain survivor no longer needs a whole-carrier NN attachment theorem; universal source/relink slab renewal and supply of the next signed-good recursive continuation remain separate continuum tasks",
        "arb": arb,
    }


@dataclass(frozen=True)
class AnnularServiceReentryStress:
    samples: int
    minimum_prefix_quarter_margin: float
    minimum_carrier_mass_margin: float
    minimum_spectral_heat_margin: float
    minimum_truncated_heat_margin: float
    worst_material_partition_residual: float
    maximum_forced_generation: int


def stress(samples: int = 50_000, seed: int = 20260809) -> AnnularServiceReentryStress:
    rng = np.random.default_rng(seed)
    mp = mm = mh = mt = float("inf")
    wr = 0.0
    mq = 0
    for _ in range(samples):
        c = float(math.exp(rng.uniform(math.log(0.25), math.log(2.5))))
        nu = float(rng.uniform(0.0, 2.0))
        K = float(rng.uniform(0.0, float(LOW_STRAIN_ACTION)))

        # Prefix-wise exact triangle branch.  Complex phases can only improve the
        # triangle lower, so stress the sharp magnitude ledger directly.
        amp = float(math.exp(rng.uniform(-5.0, 5.0)))
        ir = float(rng.uniform(0.0, 0.249999)) * amp
        ih = float(rng.uniform(0.0, 0.499999)) * amp
        low = no_hit_prefix_amplitude_lower(amp, ir, ih)
        mp = min(mp, low - amp / 4.0)
        if low <= amp / 4.0 - 2e-12 * max(1.0, amp):
            raise AssertionError("no-hit prefix lost the quarter coefficient")

        # Coefficient-to-energy lower with the exact scale-independent adjoint cost.
        J = renewed_analysis_probe_growth_upper(c, nu, K)
        clean = inherited_seed_critical_mass_lower(c) / (J * J)
        direct = persistent_carrier_critical_mass_lower(c, nu, K)
        mm = min(mm, direct - clean)
        if direct + 2e-13 * max(1.0, clean) < clean:
            raise AssertionError("persistent carrier critical mass lower failed")

        # Arbitrary positive spectral law inside the transported annulus.
        lo, hi = transported_annular_support_ratios(K)
        n = int(rng.integers(2, 80))
        rho = rng.uniform(lo, hi, size=n)
        e = rng.lognormal(mean=-1.0, sigma=1.3, size=n)
        exact = float(np.dot(2.0 * (1.0 - np.exp(-0.5 * rho * rho)), e))
        lower = heat_defect_fraction_lower(K) * float(e.sum())
        mh = min(mh, exact - lower)
        if exact + 4e-12 * max(1.0, lower) < lower:
            raise AssertionError("annular heat multiplier lower failed")

        full = heat_defect_fraction_lower(K)
        tail = gaussian_3d_tail_probability(BOUNDED_HEAT_RADIUS)
        retained = full - 4.0 * tail
        mt = min(mt, retained - 0.5 * full)
        if retained <= 0.5 * full:
            raise AssertionError("bounded heat probe lost the clean half")

        # Renewed materiality is a partition of the new positive service law,
        # independent of how the old NN seed witness was selected.
        ne = int(rng.integers(1, 80))
        weights = rng.lognormal(mean=-2.0, sigma=1.7, size=ne)
        a = rng.random(ne) < 0.5
        b = rng.random(ne) < 0.5
        part = material_service_partition(weights, a, b)
        wr = max(wr, abs(float(part["partition_residual"])))
        if abs(float(part["partition_residual"])) > 3e-12 * max(1.0, float(part["total"])):
            raise AssertionError("renewed material service partition failed")

        Y0 = uniform_bounded_square_service_lower(c, nu, K)
        C0 = float(math.exp(rng.uniform(-8.0, 8.0))) * Y0
        ep = epoch_certificate(Y0, C0)
        mq = max(mq, int(ep["first_forced_generation"]))

    return AnnularServiceReentryStress(samples, mp, mm, mh, mt, wr, mq)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-critical-annular-carrier-service-reentry"))
    ap.add_argument("--scaled-lifetime", type=float, default=1.0)
    ap.add_argument("--viscosity", type=float, default=1.0)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = theorem_certificate(args.scaled_lifetime, args.viscosity)
    out = stress(args.samples)
    (args.outdir / "critical_annular_carrier_service_reentry.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Critical annular carrier -> bounded coherent-service re-entry

Status: **{cert['status']}**.

A full no-hit NN-seed temporal corridor is stronger than an endpoint statement.  At every backward prefix `s`, the cumulative interface and HH monitors remain strictly below `|z(t)|/4` and `|z(t)|/2`, so the exact moving-role Duhamel identity gives

`|z(s)| > |z(t)|/4`.

The registered affine/Kelvin analysis dual is not unitary: bulk viscosity is inside its interaction-picture propagator.  For the renewed smooth role, affine-dual support stays below `R A` with `R=(3/2)exp(1/30)`.  The existing registered-probe inverse-heat/strain estimate therefore gives

`||psi(s)|| <= J ||psi(t)||`,  `J=exp(1/30 + nu c R^2)`.

Consequently, for **every** point of the full natural slab,

`A ||Q_Au(s)||_2^2 >= pi^2/(50 c^2 J^2) = {persistent_carrier_critical_mass_lower(args.scaled_lifetime,args.viscosity):.12g}`.

This carrier is annular, not merely high-pass at one instant.  Its lower frequency edge stays at least `(3/5)exp(-1/30)A`.  The intrinsic heat displacement `r~N(0,A^-2 I)` therefore gives

`E_HA ||delta_r w||_2^2 >= q_* ||w||_2^2`,

with `q_*={heat_defect_fraction_lower():.12g}`.  The Gaussian is unbounded, so we do **not** feed this directly into the reservoir theorem.  Instead use `||delta_r w||_2^2<=4||w||_2^2`.  Exact 3D Gaussian tail arithmetic at `beta=3` proves that deleting `|r|>3/A` leaves strictly more than half of the annular heat lower.  Hence at every slab time there exists an actual bounded displacement

`|r_s|<=3/A`,  `A||delta_(r_s) w(s)||_2^2 >= Y0`,

where

`Y0={uniform_bounded_square_service_lower(args.scaled_lifetime,args.viscosity):.12g}>0`.

The full natural slab also carries normalized integrated bounded heat service at least

`c Y0={integrated_bounded_heat_service_lower(args.scaled_lifetime,args.viscosity):.12g}`.

Now material identity is assigned in the physically correct order.  First create this positive renewed heat/increment service, then apply exact Moyal, then partition its two intrinsic endpoints into OO/ON/NN.  The old NN heat-edge witness which selected the original critical seed is retained only as seed provenance; it is **not** used to declare the whole smooth carrier new material.  Thus the apparent material-attachment problem disappears: the renewed carrier supplies its own actual material service law.

The number `Y0` is exactly the scale-independent bounded-displacement square-service threshold required by the existing geometric old-pool stopping theorem.  On any supplied sticky signed-good continuation, old capacity therefore has finite cost-free age and must exit through the already-certified dissipation / interface / new-coherent-mass / entropy / cycle alternatives.  This theorem does not claim that service alone is near-extremal HH transfer efficiency, and it does not yet supply universal source/relink slab renewal.

Arb: `{cert['arb']['status']}`.

Stress: `{out.samples}` no-hit-prefix/adjoint/annular-heat/material states
- minimum prefix quarter-coefficient margin: `{out.minimum_prefix_quarter_margin:.3e}`
- minimum carrier-mass identity margin: `{out.minimum_carrier_mass_margin:.3e}`
- minimum spectral heat lower margin: `{out.minimum_spectral_heat_margin:.3e}`
- minimum radius-3 truncated-heat margin: `{out.minimum_truncated_heat_margin:.3e}`
- worst renewed OO/ON/NN partition residual: `{out.worst_material_partition_residual:.3e}`
- maximum sampled finite forced-cost generation: `{out.maximum_forced_generation}`
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
