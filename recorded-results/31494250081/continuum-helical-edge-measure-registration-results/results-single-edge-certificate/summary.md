# Certified single-edge stability

Status: **CERTIFIED** (Arb / python-flint, 160 bits).

## Theorem constants

- local imbalance radius: `2/25`
- local mean-scale radius: `2/25`
- mixed bound: `Def >= (1/50) |u| + (1/1) v^2`
- local Hodge conversion: `Def >= 1/2 (r_p^2+r_q^2)`
- adverse sharp-cutoff Mellin retention: `>= 9/10` of the upper progress segment
- smooth common-midgap moat: shell `2/25`, filter `1/20`, residual `>= 9/250`
- physical good-core threshold: `eta=1/10000`, gap radius `<= 1/80`, transfer/capacity condition `< 53/50`
- full all-scale Mellin countermodel: `x=13/40`, `y=17/20`, signs `(-,+,-)`, coefficient `> 3/2 J*`
- global exclusion outside the local box: `Def >= 1/100`

## Certified enclosures

- r*: `[0.6109041016 +/- 2.01e-11]`
- gamma*: `[0.4928152853 +/- 5.33e-11]`
- J*: `[0.1001101759 +/- 5.38e-11]`
- certified lower bound for local transverse derivative: `[0.046811393544502673280151283652372740107239224017 +/- 3.34e-49]`
- local derivative leaf boxes: `6` (max depth `3`)
- certified lower bound for symmetric second derivative: `6.8056815564632415771484375000000000000000000000`
- tangent derivative leaf boxes: `1` (max depth `0`)
- adverse lower/upper Mellin segment ratio upper bound: `[0.08296357711552840288318294177429262864534165217 +/- 2.33e-48]`
- smooth common-midgap moat enclosure: `[0.0364076427 +/- 3.98e-11]`
- physical transfer/capacity condition enclosure: `[1.0521543606 +/- 9.28e-11]`
- full Mellin countermodel enclosure: `[0.15717672262614943884620576785724102344384092081 +/- 8.74e-48]`
- y>=0.9 analytic upper bound: `[0.07450113509096039668367846364159696261312372323 +/- 3.40e-48]`

## Global branch-and-bound

- gap-certified boxes: `637`
- boxes absorbed by the local theorem: `93`
- maximum subdivision depth: `15`

The random stress test below is adversarial evidence only; it is not used in the proof.

## Numerical stress test

- local samples: `50000`
- global samples: `50000`
- worst mixed-bound margin: `6.690045e-05`
- worst global-gap margin: `2.392458e-02`
- largest J/J* seen globally: `0.999881415798`
