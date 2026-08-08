# Affine coherent-state Bessel budget

Status: **CERTIFIED_EQUAL_COVARIANCE_COHERENT_BESSEL_RIESZ**.

For equal physical covariance factor `L`, normalized affine Gaussian packets satisfy exactly

`|<g_(X,k,L),g_(Y,l,L)>| = exp[-|L^-1(X-Y)|^2/8-|L^T(k-l)|^2/2]`.

With intrinsic phase point `zeta=(L^-1 X/2,L^T k)`, this is simply `exp(-|zeta_a-zeta_b|^2/2)` and is invariant under a common physical affine change.

If the intrinsic phase points are `4`-separated in R^6, disjoint-ball packing gives at most `(2n+3)^6` centers in shell `[4n,4(n+1))`. Arb certifies the infinite absolute Gram row sum is `<25/4`. Schur therefore gives

`sum_a |<f,g_a>|^2 <= (25/4)||f||_2^2`.

At separation `5`, Arb further certifies the off-diagonal Gram row `<3/50`, hence the Gram spectrum lies in `[47/50,53/50]`. Therefore arbitrary synthesis coefficients in that separated equal-covariance family satisfy `sum|c_a|^2 <= (50/47)||sum c_a g_a||_2^2`. This directly supplies the Bessel/frame coefficient budget needed by old-pool erosion inside one covariance cell. Changes of covariance cell and transfer extraction before the separated coherent synthesis remain iterative-interface issues.

Stress: `50000` affine/pair checks plus 1500 finite Gram probes
- worst affine-coordinate invariance residual: `6.361e-12`
- worst overlap-coordinate residual: `8.327e-16`
- maximum sampled Gram operator norm: `1.000094490`
- maximum sampled absolute row sum: `1.000094490`
- minimum separation margin: `3.051e-01`
- minimum sampled 5-separated Gram eigenvalue: `0.999999172`
- maximum sampled 5-separated Gram eigenvalue: `1.000000828`
