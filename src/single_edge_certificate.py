from __future__ import annotations

import argparse
import json
import math
from collections import deque
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any

# The theorem constants are deliberately rational and conservative.
U0 = Fraction(2, 25)          # |log(parent ratio)| <= 0.08
V0 = Fraction(2, 25)          # |mean scale residual| <= 0.08
A_CUSP = Fraction(1, 50)      # linear cusp coefficient
B_TANGENT = Fraction(1, 1)    # quadratic common-scale coefficient
GLOBAL_GAP = Fraction(1, 100) # deficit outside local box
Y_CUTOFF = Fraction(9, 10)    # y >= .9 handled analytically
SMOOTH_SHELL = Fraction(2, 25)   # log-shell halfwidth for the packet block
SMOOTH_DELTA = Fraction(1, 20)   # smooth log-filter transition halfwidth
SMOOTH_MOAT = Fraction(9, 250)   # certified residual common moat (>0.036)
PHYSICAL_GOOD_ETA = Fraction(1, 10_000)
PHYSICAL_GAP_RADIUS = Fraction(1, 80)   # sqrt(eta)+25 eta exactly at eta=1e-4
PHYSICAL_WEIGHT_COND = Fraction(53, 50) # child-transfer/capacity condition number

# A rational bracket for the unique symmetric critical point.
RSTAR_LO = Fraction(61090410158, 100_000_000_000)
RSTAR_HI = Fraction(61090410160, 100_000_000_000)


def float_envelope(x: float, y: float) -> float:
    """Maximum single-edge J over helicity signs for 0 < x <= y < 1.

    Child magnitude is normalized to one.  This is the exact sign-reduced
    envelope proved in docs/single_edge_stability_certificate.md.
    """
    if not (0.0 < x <= y < 1.0 and x + y > 1.0):
        return 0.0
    s = x + y
    d = y - x
    j2 = (
        math.log(1.0 / y) ** 2
        * s * s
        * (s * s - 1.0)
        * (1.0 + d) ** 3
        * (1.0 - d)
        / (8.0 * (s * s - d * d) ** 2)
    )
    return math.sqrt(max(0.0, j2))


def float_rstar() -> float:
    lo, hi = float(RSTAR_LO), float(RSTAR_HI)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f = -math.log(mid) - 4.0 * mid * mid + 1.0
        if f > 0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def float_jstar() -> float:
    r = float_rstar()
    return math.sqrt(4.0 * r * r - 1.0) * math.log(1.0 / r) / (4.0 * math.sqrt(2.0) * r)


def residual_coordinates(x: float, y: float) -> tuple[float, float]:
    """Return (u,v) with u=log(y/x)>=0 and v=child-mean-parent-gamma*."""
    r = float_rstar()
    gamma = math.log(1.0 / r)
    u = math.log(y / x)
    v = -0.5 * (math.log(x) + math.log(y)) - gamma
    return u, v


def local_mixed_rhs(x: float, y: float) -> float:
    u, v = residual_coordinates(x, y)
    return float(A_CUSP) * u + float(B_TANGENT) * v * v


def hodge_residual_energy(x: float, y: float) -> float:
    u, v = residual_coordinates(x, y)
    return 0.5 * u * u + 2.0 * v * v


@dataclass(frozen=True)
class Box:
    y0: Fraction
    y1: Fraction
    l0: Fraction
    l1: Fraction
    depth: int = 0

    def split(self) -> tuple["Box", "Box"]:
        # Split the relatively wider normalized coordinate.
        wy = (self.y1 - self.y0) / (Y_CUTOFF - Fraction(1, 2))
        wl = self.l1 - self.l0
        if wy >= wl:
            m = (self.y0 + self.y1) / 2
            return (
                Box(self.y0, m, self.l0, self.l1, self.depth + 1),
                Box(m, self.y1, self.l0, self.l1, self.depth + 1),
            )
        m = (self.l0 + self.l1) / 2
        return (
            Box(self.y0, self.y1, self.l0, m, self.depth + 1),
            Box(self.y0, self.y1, m, self.l1, self.depth + 1),
        )


def _qstr(q: Fraction) -> str:
    return f"{q.numerator}/{q.denominator}"


def _run_arb_certificate(max_depth: int = 28) -> dict[str, Any]:
    try:
        from flint import arb, ctx
    except ImportError as exc:  # pragma: no cover - exercised in Actions
        raise RuntimeError("python-flint is required for the rigorous certificate") from exc

    ctx.prec = 160

    def aq(q: Fraction):
        return arb(_qstr(q))

    def interval(a: Fraction, b: Fraction):
        return aq(a).union(aq(b))

    zero = arb(0)
    one = arb(1)
    two = arb(2)
    four = arb(4)
    sqrt2 = arb(2).sqrt()

    rlo, rhi = aq(RSTAR_LO), aq(RSTAR_HI)

    def critical_eq(r):
        return -r.log() - 4 * r * r + 1

    flo = critical_eq(rlo)
    fhi = critical_eq(rhi)
    if not (flo > zero and fhi < zero):
        raise AssertionError(f"r* bracket failed: f(lo)={flo}, f(hi)={fhi}")
    # f'(r)=-1/r-8r < 0 for r>0, so the bracket contains the unique root.
    rstar = rlo.union(rhi)
    gamma = -rstar.log()
    jstar = ((4 * rstar * rstar - 1).sqrt() * gamma) / (4 * sqrt2 * rstar)
    jstar_lo = jstar.lower()

    # A common smooth midgap moat for log-shell halfwidth 2/25 and
    # filter halfwidth 1/20.  If parent and child shells each have this
    # halfwidth, a transfer-weighted midgap can shift by at most one shell
    # halfwidth, leaving gamma/2 - 2*sigma - delta.
    smooth_moat_margin = gamma / 2 - 2 * aq(SMOOTH_SHELL) - aq(SMOOTH_DELTA)
    if not (smooth_moat_margin > aq(SMOOTH_MOAT)):
        raise AssertionError(f"smooth common-midgap moat failed: {smooth_moat_margin}")

    # On signed-good edges r=m*c >= 1-eta, one has Def=1-m<=eta.
    # The mixed theorem gives u<=50 eta=1/200 and |v|<=sqrt(eta)=1/100,
    # hence |log(q/p)-gamma*|<=1/80.  Since F=A J* r=T log(q/p),
    # the child-transfer/capacity density condition number is the expression below.
    gap_radius = aq(PHYSICAL_GAP_RADIUS)
    eta_good = aq(PHYSICAL_GOOD_ETA)
    physical_weight_cond = (gamma + gap_radius) / ((1 - eta_good) * (gamma - gap_radius))
    if not (physical_weight_cond < aq(PHYSICAL_WEIGHT_COND)):
        raise AssertionError(f"physical transfer/capacity weight comparison failed: {physical_weight_cond}")

    # --- Local certificate in exact log coordinates. ---
    # Direct interval evaluation on the whole rectangle has severe dependency
    # wrapping, so certify the derivative by adaptive dyadic subdivision.
    @dataclass(frozen=True)
    class LocalBox:
        u0: Fraction
        u1: Fraction
        v0: Fraction
        v1: Fraction
        depth: int = 0

        def split(self):
            wu = (self.u1 - self.u0) / U0
            wv = (self.v1 - self.v0) / (2 * V0)
            if wu >= wv:
                m = (self.u0 + self.u1) / 2
                return (LocalBox(self.u0, m, self.v0, self.v1, self.depth + 1),
                        LocalBox(m, self.u1, self.v0, self.v1, self.depth + 1))
            m = (self.v0 + self.v1) / 2
            return (LocalBox(self.u0, self.u1, self.v0, m, self.depth + 1),
                    LocalBox(self.u0, self.u1, m, self.v1, self.depth + 1))

    def transverse_derivative(box: LocalBox):
        u = interval(box.u0, box.u1)
        v = interval(box.v0, box.v1)
        R = rstar * (-v).exp()
        halfu = u / 2
        c = halfu.cosh()
        sh = halfu.sinh()
        S = 2 * R * c
        delta = 2 * R * sh
        L = gamma + v - u / 2
        if not (L.lower() > zero and (S * S - 1).lower() > zero and (1 - delta).lower() > zero):
            raise AssertionError(f"local subbox left forward-triad domain: {box}")
        Hsqrt = (((1 + delta) ** 3) * (1 - delta)).sqrt()
        J = L * c * (4 * R * R * c * c - 1).sqrt() * Hsqrt / (4 * sqrt2 * R)
        g = (
            -1 / (2 * L)
            + halfu.tanh() / 2
            + S * delta / (2 * (S * S - 1))
            + S * (1 - 2 * delta) / (2 * (1 - delta * delta))
        )
        return -(J / jstar) * g

    local_queue: deque[LocalBox] = deque([LocalBox(Fraction(0), U0, -V0, V0, 0)])
    local_boxes = 0
    local_max_depth = 0
    worst_du_lower = None
    while local_queue:
        box = local_queue.pop()
        du = transverse_derivative(box)
        if du.lower() > aq(A_CUSP):
            local_boxes += 1
            local_max_depth = max(local_max_depth, box.depth)
            if worst_du_lower is None or du.lower() < worst_du_lower:
                worst_du_lower = du.lower()
            continue
        if box.depth >= 16:
            raise AssertionError(f"transverse derivative subdivision failed: box={box}, du={du}")
        b0, b1 = box.split()
        local_queue.append(b0)
        local_queue.append(b1)

    assert worst_du_lower is not None

    # Along u=0, D(v)=1-J(0,v)/J*, D(0)=D'(0)=0 by the defining root equation.
    # Here q=sqrt(4-R^{-2}); differentiating twice gives the exact expression below.
    # As for the transverse direction, subdivide to eliminate interval wrapping.
    @dataclass(frozen=True)
    class TangentBox:
        v0: Fraction
        v1: Fraction
        depth: int = 0

        def split(self):
            m = (self.v0 + self.v1) / 2
            return TangentBox(self.v0, m, self.depth + 1), TangentBox(m, self.v1, self.depth + 1)

    def tangent_second_derivative(box: TangentBox):
        v = interval(box.v0, box.v1)
        Rt = rstar * (-v).exp()
        aa = 1 / (Rt * Rt)
        q = (4 - aa).sqrt()
        Lv = gamma + v
        return (
            2 * aa / q + Lv * (2 * aa / q + aa * aa / (q ** 3))
        ) / (4 * sqrt2 * jstar)

    tangent_queue: deque[TangentBox] = deque([TangentBox(-V0, V0, 0)])
    tangent_boxes = 0
    tangent_max_depth = 0
    worst_d2_lower = None
    while tangent_queue:
        box = tangent_queue.pop()
        d2 = tangent_second_derivative(box)
        if d2.lower() > 2 * aq(B_TANGENT):
            tangent_boxes += 1
            tangent_max_depth = max(tangent_max_depth, box.depth)
            if worst_d2_lower is None or d2.lower() < worst_d2_lower:
                worst_d2_lower = d2.lower()
            continue
        if box.depth >= 16:
            raise AssertionError(f"tangent curvature subdivision failed: box={box}, d2={d2}")
        b0, b1 = box.split()
        tangent_queue.append(b0)
        tangent_queue.append(b1)

    assert worst_d2_lower is not None

    # This proves D >= A u + B v^2 on the whole local rectangle by integration.
    # Since u<=2/25, A u >= u^2/4, hence D >= 1/2*(u^2/2+2v^2).
    c_hodge = min(
        aq(B_TANGENT) / 2,
        2 * aq(A_CUSP) / aq(U0),
    )
    if not (c_hodge >= arb("1/2")):
        raise AssertionError(f"unexpected Hodge coefficient: {c_hodge}")

    # --- Exact sharp-cutoff Mellin-flux bridge on the same local box. ---
    # For the adverse maximizing orbit (+,-,-), after factoring out the common
    # triad phase/amplitude factor, the upper forward segment is
    #   (x+y) log(1/y),
    # while the lower cutoff segment is -(1-y) log(y/x).
    # Monotone endpoint bounds give a rigorous uniform leakage ratio < 1/10.
    exp_minus_v0 = (-aq(V0)).exp()
    y_lower = rlo * exp_minus_v0
    sum_lower = 2 * rlo * exp_minus_v0
    progress_lower = -rhi.log() - aq(V0) - aq(U0) / 2
    mellin_adverse_ratio_upper = aq(U0) * (1 - y_lower) / (sum_lower * progress_lower)
    if not (mellin_adverse_ratio_upper < arb("1/10")):
        raise AssertionError(f"Mellin lower-segment leakage too large: {mellin_adverse_ratio_upper}")

    # --- Global exclusion outside the local rectangle. ---
    # Exact sign reduction leaves the envelope.  For y>=0.9,
    # J <= log(1/y)/sqrt(2) <= log(10/9)/sqrt(2).
    corner_upper = (arb("10/9").log()) / sqrt2
    target = (one - aq(GLOBAL_GAP)) * jstar_lo
    if not (corner_upper < target):
        raise AssertionError(f"analytic y>=0.9 exclusion failed: {corner_upper} vs {target}")
    target2 = target * target

    initial = Box(Fraction(1, 2), Y_CUTOFF, Fraction(0), Fraction(1), 0)
    queue: deque[Box] = deque([initial])
    certified_gap_boxes = 0
    certified_local_boxes = 0
    max_seen_depth = 0
    max_queue = 1
    worst_gap_ratio = 0.0

    while queue:
        box = queue.pop()
        max_seen_depth = max(max_seen_depth, box.depth)
        y = interval(box.y0, box.y1)
        lam = interval(box.l0, box.l1)
        a = 2 * y - 1
        s = 1 + lam * a
        d = (1 - lam) * a
        x = (s - d) / 2

        # The square parameterization enforces x+y>=1 and x<=y exactly.
        ell = (1 / y).log()
        numerator = ell * ell * s * s * (s * s - 1) * ((1 + d) ** 3) * (1 - d)
        denominator = 8 * (s * s - d * d) ** 2
        j2 = numerator / denominator

        ubox = (y / x).log()
        vbox = -((x.log() + y.log()) / 2) - gamma
        inside_local = (
            ubox.upper() <= aq(U0)
            and vbox.lower() >= -aq(V0)
            and vbox.upper() <= aq(V0)
        )
        if inside_local:
            certified_local_boxes += 1
            continue

        if j2.upper() < target2.lower():
            certified_gap_boxes += 1
            # Reporting only; theorem logic above uses Arb comparisons.
            try:
                ratio = float(j2.upper()) / float(target2.lower())
                worst_gap_ratio = max(worst_gap_ratio, ratio)
            except Exception:
                pass
            continue

        if box.depth >= max_depth:
            raise AssertionError(
                "global branch-and-bound reached max depth without certification: "
                f"box={box}, j2={j2}, target2={target2}, u={ubox}, v={vbox}"
            )
        b0, b1 = box.split()
        queue.append(b0)
        queue.append(b1)
        max_queue = max(max_queue, len(queue))

    return {
        "precision_bits": int(ctx.prec),
        "rstar_bracket": [_qstr(RSTAR_LO), _qstr(RSTAR_HI)],
        "rstar_ball": str(rstar),
        "gamma_ball": str(gamma),
        "jstar_ball": str(jstar),
        "local_box": {"u_max": _qstr(U0), "v_abs_max": _qstr(V0)},
        "mixed_stability": {"A": _qstr(A_CUSP), "B": _qstr(B_TANGENT)},
        "transverse_derivative_lower_bound": str(worst_du_lower),
        "local_derivative_boxes": local_boxes,
        "local_derivative_max_depth": local_max_depth,
        "tangent_second_derivative_lower_bound": str(worst_d2_lower),
        "tangent_derivative_boxes": tangent_boxes,
        "tangent_derivative_max_depth": tangent_max_depth,
        "hodge_coefficient_lower_bound": "1/2",
        "mellin_adverse_ratio_upper_bound": str(mellin_adverse_ratio_upper),
        "mellin_flux_retention_lower_bound": "9/10",
        "smooth_midgap_shell_halfwidth": _qstr(SMOOTH_SHELL),
        "smooth_midgap_filter_halfwidth": _qstr(SMOOTH_DELTA),
        "smooth_midgap_moat_lower_bound": _qstr(SMOOTH_MOAT),
        "smooth_midgap_moat_ball": str(smooth_moat_margin),
        "physical_good_eta": _qstr(PHYSICAL_GOOD_ETA),
        "physical_gap_radius": _qstr(PHYSICAL_GAP_RADIUS),
        "physical_weight_condition_upper": _qstr(PHYSICAL_WEIGHT_COND),
        "physical_weight_condition_ball": str(physical_weight_cond),
        "global_gap": _qstr(GLOBAL_GAP),
        "y_cutoff": _qstr(Y_CUTOFF),
        "corner_upper_ball": str(corner_upper),
        "global_boxes_gap": certified_gap_boxes,
        "global_boxes_local": certified_local_boxes,
        "global_max_depth": max_seen_depth,
        "global_max_queue": max_queue,
        "worst_reported_gap_ratio": worst_gap_ratio,
        "status": "CERTIFIED",
    }


def numerical_stress(samples: int = 100_000, seed: int = 20260807) -> dict[str, float]:
    """Adversarial/random check only; never used by the theorem certificate."""
    import numpy as np

    rng = np.random.default_rng(seed)
    jstar = float_jstar()
    worst_local = float("inf")
    worst_global = float("inf")
    max_ratio = 0.0

    # Local samples in the exact (u,v) coordinates.
    r = float_rstar()
    gamma = math.log(1.0 / r)
    us = rng.uniform(0.0, float(U0), samples)
    vs = rng.uniform(-float(V0), float(V0), samples)
    for u, v in zip(us, vs):
        R = r * math.exp(-v)
        x = R * math.exp(-u / 2)
        y = R * math.exp(u / 2)
        J = float_envelope(x, y)
        deficit = 1.0 - J / jstar
        margin = deficit - float(A_CUSP) * u - v * v
        worst_local = min(worst_local, margin)

    # Global square parameterization; only test points outside the local box.
    ys = rng.uniform(0.5, 0.999999, samples)
    ls = rng.uniform(0.0, 1.0, samples)
    for y, lam in zip(ys, ls):
        a = 2.0 * y - 1.0
        s = 1.0 + lam * a
        d = (1.0 - lam) * a
        x = 0.5 * (s - d)
        J = float_envelope(x, y)
        ratio = J / jstar
        max_ratio = max(max_ratio, ratio)
        u = math.log(y / x)
        v = -0.5 * (math.log(x) + math.log(y)) - gamma
        if u > float(U0) or abs(v) > float(V0):
            worst_global = min(worst_global, (1.0 - ratio) - float(GLOBAL_GAP))

    return {
        "samples_local": samples,
        "samples_global": samples,
        "worst_local_mixed_margin": float(worst_local),
        "worst_global_gap_margin": float(worst_global),
        "largest_global_ratio_seen": float(max_ratio),
    }


def render_summary(cert: dict[str, Any], stress: dict[str, float] | None = None) -> str:
    lines = [
        "# Certified single-edge stability",
        "",
        f"Status: **{cert['status']}** (Arb / python-flint, {cert['precision_bits']} bits).",
        "",
        "## Theorem constants",
        "",
        f"- local imbalance radius: `{cert['local_box']['u_max']}`",
        f"- local mean-scale radius: `{cert['local_box']['v_abs_max']}`",
        f"- mixed bound: `Def >= ({cert['mixed_stability']['A']}) |u| + ({cert['mixed_stability']['B']}) v^2`",
        f"- local Hodge conversion: `Def >= {cert['hodge_coefficient_lower_bound']} (r_p^2+r_q^2)`",
        f"- adverse sharp-cutoff Mellin retention: `>= {cert['mellin_flux_retention_lower_bound']}` of the upper progress segment",
        f"- smooth common-midgap moat: shell `{cert['smooth_midgap_shell_halfwidth']}`, filter `{cert['smooth_midgap_filter_halfwidth']}`, residual `>= {cert['smooth_midgap_moat_lower_bound']}`",
        f"- physical good-core threshold: `eta={cert['physical_good_eta']}`, gap radius `<= {cert['physical_gap_radius']}`, transfer/capacity condition `< {cert['physical_weight_condition_upper']}`",
        f"- global exclusion outside the local box: `Def >= {cert['global_gap']}`",
        "",
        "## Certified enclosures",
        "",
        f"- r*: `{cert['rstar_ball']}`",
        f"- gamma*: `{cert['gamma_ball']}`",
        f"- J*: `{cert['jstar_ball']}`",
        f"- certified lower bound for local transverse derivative: `{cert['transverse_derivative_lower_bound']}`",
        f"- local derivative leaf boxes: `{cert['local_derivative_boxes']}` (max depth `{cert['local_derivative_max_depth']}`)",
        f"- certified lower bound for symmetric second derivative: `{cert['tangent_second_derivative_lower_bound']}`",
        f"- tangent derivative leaf boxes: `{cert['tangent_derivative_boxes']}` (max depth `{cert['tangent_derivative_max_depth']}`)",
        f"- adverse lower/upper Mellin segment ratio upper bound: `{cert['mellin_adverse_ratio_upper_bound']}`",
        f"- smooth common-midgap moat enclosure: `{cert['smooth_midgap_moat_ball']}`",
        f"- physical transfer/capacity condition enclosure: `{cert['physical_weight_condition_ball']}`",
        f"- y>=0.9 analytic upper bound: `{cert['corner_upper_ball']}`",
        "",
        "## Global branch-and-bound",
        "",
        f"- gap-certified boxes: `{cert['global_boxes_gap']}`",
        f"- boxes absorbed by the local theorem: `{cert['global_boxes_local']}`",
        f"- maximum subdivision depth: `{cert['global_max_depth']}`",
        "",
        "The random stress test below is adversarial evidence only; it is not used in the proof.",
    ]
    if stress is not None:
        lines += [
            "",
            "## Numerical stress test",
            "",
            f"- local samples: `{stress['samples_local']}`",
            f"- global samples: `{stress['samples_global']}`",
            f"- worst mixed-bound margin: `{stress['worst_local_mixed_margin']:.6e}`",
            f"- worst global-gap margin: `{stress['worst_global_gap_margin']:.6e}`",
            f"- largest J/J* seen globally: `{stress['largest_global_ratio_seen']:.12f}`",
        ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results-single-edge-certificate")
    ap.add_argument("--max-depth", type=int, default=28)
    ap.add_argument("--stress-samples", type=int, default=100_000)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    cert = _run_arb_certificate(max_depth=args.max_depth)
    stress = numerical_stress(args.stress_samples) if args.stress_samples > 0 else None
    payload = {"certificate": cert, "stress": stress}
    (out / "single_edge_certificate.json").write_text(json.dumps(payload, indent=2))
    summary = render_summary(cert, stress)
    (out / "summary.md").write_text(summary)
    print(summary)


if __name__ == "__main__":
    main()
