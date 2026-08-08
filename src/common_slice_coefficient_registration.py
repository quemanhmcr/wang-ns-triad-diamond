from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

INHERIT_FRACTION = 1.0 / 4.0
RESIDUAL_FRACTION = 1.0 / 4.0
GENERATED_FRACTION = 1.0 / 2.0
CONTINUING_COEFFICIENT_FRACTION = INHERIT_FRACTION
CONTINUING_PRODUCT_FRACTION = INHERIT_FRACTION**2
ASYNC_CONE = 10.0 / 39.0
COMMON_SLICE_OFFSET = 2.0 / 5.0


def common_slice_natural_window_margin(
    support_span_ratio: float = ASYNC_CONE,
    slice_offset: float = COMMON_SLICE_OFFSET,
) -> float:
    """Fraction of the shortest parent natural lifetime left beyond every event."""
    a = float(support_span_ratio)
    b = float(slice_offset)
    if a < 0 or b < 0:
        raise ValueError("nonnegative synchronization geometry required")
    return 1.0 - a - b


def exact_adjoint_residual(z_event: complex, z_slice: complex, i_hh: complex, i_r: complex) -> complex:
    """Residual of z_event=z_slice+I_HH+I_R on one parent registration interval."""
    return complex(z_event) - complex(z_slice) - complex(i_hh) - complex(i_r)


def registration_first_stop(
    z_event: complex,
    z_slice: complex,
    i_hh: complex,
    i_r: complex,
    material_relink: bool = False,
) -> dict[str, object]:
    """Classify an event-to-common-slice registration by the exact adjoint gate.

    The exact triangle gate guarantees at least one of inherited coefficient,
    classified residual, or HH generation.  A genuine material relink is another
    physical obstruction.  All threshold-crossing obstructions are returned as a
    set; this registration helper never chooses a primary cause by theorem-name
    order.  The single-charge compiler owns simultaneous-stop partitioning.  When
    no obstruction fires, the coefficient is registered on the common slice with
    the clean inherited fraction 1/4.
    """
    A = abs(complex(z_event))
    if A <= 0:
        raise ValueError("positive event coefficient required")
    res = abs(exact_adjoint_residual(z_event, z_slice, i_hh, i_r))
    if res > 2e-12 * max(1.0, A, abs(z_slice), abs(i_hh), abs(i_r)):
        raise ValueError("adjoint decomposition is not exact within tolerance")
    hits: list[str] = []
    if material_relink:
        hits.append("material_relink")
    if abs(i_r) >= RESIDUAL_FRACTION * A:
        hits.append("classified_residual")
    if abs(i_hh) >= GENERATED_FRACTION * A:
        hits.append("hh_generation")
    if hits:
        branch = {
            "material_relink": "material_relink_stop",
            "classified_residual": "classified_residual_stop",
            "hh_generation": "hh_generation_stop",
        }[hits[0]] if len(hits) == 1 else "multiple_causal_stops_before_common_slice"
        return {
            "branch": branch,
            "stop_causes": tuple(hits),
            "continuing": False,
            "event_amplitude": A,
            "residual_impulse": abs(i_r),
            "residual_threshold": RESIDUAL_FRACTION * A,
            "hh_impulse": abs(i_hh),
            "hh_threshold": GENERATED_FRACTION * A,
            "registered_amplitude_lower": 0.0,
            "primary_selected": False,
        }
    # If neither source branch fires, triangle inequality forces inheritance.
    inherited = abs(z_slice)
    clean = INHERIT_FRACTION * A
    if inherited < clean - 3e-12 * max(1.0, A):
        raise AssertionError("exact adjoint gate left no inherit/generate/residual branch")
    return {
        "branch": "registered_material_inheritance",
        "continuing": True,
        "event_amplitude": A,
        "slice_amplitude": inherited,
        "registered_amplitude_lower": clean,
        "registration_fraction": INHERIT_FRACTION,
    }


def continuing_pair_product_lower(event_product: float) -> float:
    if event_product < 0:
        raise ValueError("nonnegative event product required")
    return CONTINUING_PRODUCT_FRACTION * event_product


def register_pair(
    first: tuple[complex, complex, complex, complex],
    second: tuple[complex, complex, complex, complex],
    first_relink: bool = False,
    second_relink: bool = False,
) -> dict[str, object]:
    a = registration_first_stop(*first, material_relink=first_relink)
    b = registration_first_stop(*second, material_relink=second_relink)
    if not bool(a["continuing"]) or not bool(b["continuing"]):
        stops = [str(x["branch"]) for x in (a, b) if not bool(x["continuing"])]
        return {"branch": "pair_stops_before_common_slice", "stops": stops, "first": a, "second": b}
    event_product = float(a["event_amplitude"]) * float(b["event_amplitude"])
    actual_slice_product = float(a["slice_amplitude"]) * float(b["slice_amplitude"])
    clean = continuing_pair_product_lower(event_product)
    if actual_slice_product < clean - 4e-12 * max(1.0, event_product):
        raise AssertionError("continuing pair lost more than the clean product registration factor")
    return {
        "branch": "pair_registered_on_common_slice",
        "continuing": True,
        "event_product": event_product,
        "slice_product": actual_slice_product,
        "registered_product_lower": clean,
        "product_fraction": CONTINUING_PRODUCT_FRACTION,
        "first": a,
        "second": b,
    }


def theorem_certificate() -> dict[str, object]:
    margin = common_slice_natural_window_margin()
    return {
        "status": "EXACT_FIRST_STOP_COMMON_SLICE_COEFFICIENT_REGISTRATION_ON_SELECTED_ROLE_MODEL",
        "geometry": f"alpha<=10/39 and offset=2/5 leave natural-window margin {margin:.12g}=67/195",
        "identity": "z_event=z_slice+I_HH+I_R under the same adjoint Kelvin interaction picture",
        "continuing": "if |I_R|<A/4 and |I_HH|<A/2 and no material relink occurs, then |z_slice|>=A/4",
        "pair": "two continuing parents retain at least 1/16 of their event coefficient product on the common slice",
        "first_stop": "failure to register returns the complete set of earlier HH-generation, classified-residual/source, and material-relink obstructions; no primary is selected here",
        "single_charge": "simultaneous registration obstructions are delegated as a set to the physical branch compiler, which alone owns tie partitioning",
        "weights": "this registration theorem does not replace physical energy weights; expanded HH nodes are still weighted by actual positive child-energy work",
        "continuum_status": "exact after the selected moving parent-role coefficient equation exists; constructing that equation for every recursive continuum SGS block remains the outer-role extraction bridge",
    }


@dataclass(frozen=True)
class RegistrationStress:
    samples: int
    minimum_natural_window_margin: float
    minimum_continuing_fraction_margin: float
    minimum_pair_product_margin: float
    worst_adjoint_identity_residual: float
    branch_counts: dict[str, int]


def _random_exact_decomposition(rng: np.random.Generator, force_mode: int) -> tuple[complex, complex, complex, complex]:
    # Normalize event amplitude to a random positive size and choose phases freely.
    A = float(math.exp(rng.uniform(-5.0, 5.0)))
    theta = float(rng.uniform(-math.pi, math.pi))
    z_event = A * complex(math.cos(theta), math.sin(theta))
    if force_mode == 0:
        # Continuing case: keep both impulses strictly below their thresholds and
        # define z_slice by the exact identity.  Rejection sampling enforces the
        # inherited lower bound automatically rather than building it in.
        for _ in range(1000):
            rmag = float(rng.uniform(0.0, 0.245)) * A
            hmag = float(rng.uniform(0.0, 0.49)) * A
            rp = float(rng.uniform(-math.pi, math.pi))
            hp = float(rng.uniform(-math.pi, math.pi))
            ir = rmag * complex(math.cos(rp), math.sin(rp))
            ih = hmag * complex(math.cos(hp), math.sin(hp))
            zs = z_event - ir - ih
            if abs(zs) >= INHERIT_FRACTION * A - 1e-13 * max(1.0, A):
                return z_event, zs, ih, ir
        raise RuntimeError("could not sample continuing exact gate")
    if force_mode == 1:
        rmag = float(rng.uniform(0.25, 1.2)) * A
        rp = float(rng.uniform(-math.pi, math.pi))
        ir = rmag * complex(math.cos(rp), math.sin(rp))
        ih = 0j
        return z_event, z_event - ir, ih, ir
    hmag = float(rng.uniform(0.5, 1.4)) * A
    hp = float(rng.uniform(-math.pi, math.pi))
    ih = hmag * complex(math.cos(hp), math.sin(hp))
    ir = 0j
    return z_event, z_event - ih, ih, ir


def stress(samples: int = 50_000, seed: int = 20260809) -> RegistrationStress:
    rng = np.random.default_rng(seed)
    mn = common_slice_natural_window_margin()
    if mn <= 0:
        raise AssertionError("common slice escaped parent natural window")
    mf = mp = float("inf")
    wr = 0.0
    branches: dict[str, int] = {}
    for _ in range(samples):
        mode = int(rng.integers(0, 3))
        data = _random_exact_decomposition(rng, mode)
        out = registration_first_stop(*data)
        b = str(out["branch"])
        branches[b] = branches.get(b, 0) + 1
        wr = max(wr, abs(exact_adjoint_residual(*data)))
        if bool(out["continuing"]):
            margin = float(out["slice_amplitude"]) - INHERIT_FRACTION * float(out["event_amplitude"])
            mf = min(mf, margin)
            if margin < -3e-12 * max(1.0, float(out["event_amplitude"])):
                raise AssertionError("continuing registration fraction failed")

        # Independently stress two-parent product registration on continuing pairs.
        d1 = _random_exact_decomposition(rng, 0)
        d2 = _random_exact_decomposition(rng, 0)
        pair = register_pair(d1, d2)
        margin = float(pair["slice_product"]) - float(pair["registered_product_lower"])
        mp = min(mp, margin)
        if margin < -5e-12 * max(1.0, float(pair["event_product"])):
            raise AssertionError("pair product registration failed")

    return RegistrationStress(samples, mn, mf, mp, wr, branches)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--samples", type=int, default=50_000)
    ap.add_argument("--outdir", type=Path, default=Path("results-common-slice-coefficient-registration"))
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    out = stress(args.samples)
    cert = theorem_certificate()
    (args.outdir / "common_slice_coefficient_registration.json").write_text(
        json.dumps({"certificate": cert, "stress": asdict(out)}, indent=2), encoding="utf-8"
    )
    md = f"""# Common-slice coefficient registration by first causal stop

Status: **{cert['status']}**.

The asynchronous cone already gives `alpha<=10/39`.  With the common reference slice `s=a-(2/5)T_min`, every event lies inside its parent natural adjoint window with clean leftover margin

`1-2/5-10/39 = {out.minimum_natural_window_margin:.12g} = 67/195`.

For one eventwise parent coefficient, integrate the **same** adjoint Kelvin interaction picture from the common slice `s` to the event time `t`:

`z(t)=z(s)+I_HH[s,t]+I_R[s,t]`.

Let `A=|z(t)|`.  The exact triangle gate says

- `|z(s)|>=A/4`, or
- `|I_R|>=A/4`, or
- `|I_HH|>=A/2`.

Therefore coefficient persistence is not an assumption.  If no classified residual/source, HH-generation, or genuine material-relink stop occurs before the common slice, the event coefficient is **registered** with

`|z(s)|>=A/4`.

Two continuing parents consequently retain at least

`(1/4)^2 = 1/16`

of their event coefficient product on the synchronized slice.  This is exactly the registration factor used by the amplitude--entropy productivity theorem.

If the factor fails, the parent is not allowed to continue as an uncharged root.  Its first obstruction is itself the earlier causal event: HH generation is re-entered through the physical-energy causal gate, `R_class` delegates to its existing source/interface owner, and a genuine material-cell switch is relink/fresh ancestry.  Small frequency/covariance representative changes remain the already-summable representation `Xi` and are not the physical cause of coefficient loss.

Stress: `{out.samples}` exact complex adjoint decompositions
- common-window margin: `{out.minimum_natural_window_margin:.12g}`
- minimum continuing coefficient margin: `{out.minimum_continuing_fraction_margin:.3e}`
- minimum continuing pair-product margin: `{out.minimum_pair_product_margin:.3e}`
- worst exact adjoint identity residual: `{out.worst_adjoint_identity_residual:.3e}`
- branches: `{out.branch_counts}`

This closes common-slice **coefficient registration at the selected-role model level**.  It does not yet construct the selected moving role for every continuum smooth-SGS block; that outer-role extraction remains the PDE bridge.  Physical causal probabilities remain the actual positive child-energy work law, not raw Duhamel amplitude weights.  No global-regularity claim is made.
"""
    (args.outdir / "summary.md").write_text(md, encoding="utf-8")
    print(md)


if __name__ == "__main__":
    main()
