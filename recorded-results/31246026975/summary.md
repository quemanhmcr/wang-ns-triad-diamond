# Coherent covariance interface: changing Gaussian analysis cells has a quantified Xi cost

Status: **EXACT_GAUSSIAN_FIDELITY_AND_MOYAL_WINDOW_STABILITY**.

For centered normalized Gaussian windows with physical covariances `Sigma,Theta`,

`<g_Sigma,g_Theta> = 2^(3/2)(det Sigma det Theta)^(1/4)/det(Sigma+Theta)^(1/2)`.

If `exp(a_i)` are the generalized covariance eigenvalues, this is `prod_i cosh(a_i/2)^(-1/2)`.  Since `log cosh x <= x^2/2`,

`||g_Sigma-g_Theta||_2 <= ||a||_2/(2 sqrt(2))`.

Polarized Moyal in the **window slot** is exact:

`||V_(g_Sigma) f - V_(g_Theta) f||_(L2 phase)^2 = ||f||_2^2 ||g_Sigma-g_Theta||_2^2`.

Hence the total variation of the positive coherent energy density obeys

`int ||V_Sigma f|^2-|V_Theta f|^2| dmu <= d_log(Sigma,Theta) ||f||_2^2/sqrt(2)`.

The same bound holds after any common phase-space partition.  Thus changing the representative Gaussian covariance inside a small covariance cell creates a deterministic transfer-interface error `Xi_cov`, while common affine center/carrier transport remains exact gauge.  Large covariance jumps are not hidden in this estimate and remain a genuine relink/strain/source branch.

Stress: `50000`
- maximum exact window distance / clean log bound: `0.999320817`
- worst overlap-formula residual: `1.438e-14`
- worst window-Moyal relative residual: `7.042e-16`
- maximum coherent energy-TV / bound: `0.473502584`
