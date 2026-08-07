# Certified single-edge stability

Status: **CERTIFIED** (Arb / python-flint, 160 bits).

## Theorem constants

- local imbalance radius: `2/25`
- local mean-scale radius: `2/25`
- mixed bound: `Def >= (1/50) |u| + (1/1) v^2`
- local Hodge conversion: `Def >= 1/2 (r_p^2+r_q^2)`
- global exclusion outside the local box: `Def >= 1/100`

## Certified enclosures

- r*: `[0.6109041016 +/- 2.01e-11]`
- gamma*: `[0.4928152853 +/- 5.33e-11]`
- J*: `[0.1001101759 +/- 5.38e-11]`
- certified lower bound for local transverse derivative: `[0.046811393544502673280151283652372740107239224017 +/- 3.34e-49]`
- local derivative leaf boxes: `6` (max depth `3`)
- certified lower bound for symmetric second derivative: `6.8056815564632415771484375000000000000000000000`
- tangent derivative leaf boxes: `1` (max depth `0`)
- y>=0.9 analytic upper bound: `[0.07450113509096039668367846364159696261312372323 +/- 3.40e-48]`

## Global branch-and-bound

- gap-certified boxes: `637`
- boxes absorbed by the local theorem: `93`
- maximum subdivision depth: `15`

The random stress test below is adversarial evidence only; it is not used in the proof.

## Numerical stress test

- local samples: `100000`
- global samples: `100000`
- worst mixed-bound margin: `9.024931e-06`
- worst global-gap margin: `2.398768e-02`
- largest J/J* seen globally: `0.999860906800`
