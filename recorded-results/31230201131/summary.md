# Affine coherent Moyal energy ledger

Status: **EXACT_BY_PLANCHEREL_FUBINI**.

For any normalized affine Gaussian window `g_L`, the coherent/STFT transform obeys exactly

`integral |V_L f(X,k)|^2 dX dk/(2pi)^3 = ||f||_2^2`.

Any measurable phase-space partition therefore defines positive reservoir energies `E_C` with `sum_C E_C=||f||_2^2`.  This supplies an exact analysis-level reservoir budget `P=1`, independent of affine aspect and without synthesis-coefficient cancellation.

With an orthogonal dyadic band partition the whole-old-pool erosion theorem applies with `P=1`; a smooth LP partition pays only its fixed square-function overlap constant.  The 5-separated coherent Riesz theorem remains useful for a discrete synthesis realization.

Stress: `5000` discrete periodic Moyal regressions
- worst relative Moyal residual: `6.061e-16`
- worst cell-partition residual: `5.787e-16`
- minimum old-pool half-life margin: `0.000e+00`
