# Strain coherence / objective-gradient dichotomy

Status: **CERTIFIED**.

- coherent variation threshold: `<= 1/20` of the initial non-conformal strain
- dimensionless strain-time: `dT <= 1/30`
- pointwise Hodge coefficient: `>= 1/4`
- certified bracket enclosure: `[0.30018184166647135303892042974521300461699362944 +/- 5.14e-48]`
- local-coordinate radius enclosure: `[0.070000000000000000000000000000000000000000000000 +/- 4.04e-49]`
- time-averaged edge deficit: `>= 1/24 (dT)^2`
- coherence failure requires objective-strain variation `>= d/20`
- variable-strain traces: `20000`
- worst numerical Hodge ratio `H/(dT)^2`: `0.737209699`
- worst local shape coordinate: `0.065424692`
- worst material-strain identity residual: `8.677e-16`
- worst corotational identity residual: `5.439e-16`

The low-cost alternative is therefore no longer "the strain was not frozen".
Either coherent non-conformal strain pays a multiplier cost, or the co-rotating
strain changes by a definite fraction.  Navier--Stokes identifies the source of
that change through `D_t grad u + (grad u)^2 = -Hess p + nu Delta grad u`.
