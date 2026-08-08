# Canonical material coherent labels and summable representation Xi

Status: **EXACT_CANONICAL_MATERIAL_LABELS_AND_SUMMABLE_REPRESENTATION_XI**.

Use the intrinsic material coordinate `zeta=(L^-1 X/2,L^T k)` and a nested dyadic grid. Common affine/Kelvin transport leaves `zeta` and hence every cell address exactly unchanged. Refining a selected cell into its dyadic children is a partition identity, so Moyal energy is additive and refinement itself carries **zero** switch/relink charge.

At causal depth `j`, choose normalized frequency-cell diameter `h_j=h_0 2^-j`. Smooth SGS symbol freezing then has total representation error

`Xi_sym <= 2 A_3 L_* B_* h_0`.

Choose covariance representative mesh `delta_j=delta_0 2^-j`. The exact coherent covariance-interface theorem gives

`Xi_cov <= sqrt(2) delta_0 E_global`.

Thus the causal Duhamel pushforward, Shannon/Renyi reuse and Hodge/resistance/holonomy graph may use the **same material dyadic cell address by construction**. Frequency/covariance representatives are auxiliary and have a summable, tunably small error. A genuine physical change of selected material cell is not hidden in this representation theorem; it remains the existing Moyal switch/fresh/relink/backflow currency.

Stress: `50000` affine/nested-grid/schedule checks
- worst affine zeta residual: `1.877e-15`
- maximum nested-address failure: `0.000e+00`
- worst refinement-energy residual: `0.000e+00`
- minimum frequency schedule margin: `1.789e-16`
- minimum covariance schedule margin: `2.954e-16`
