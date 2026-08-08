# Band-limited source sampling: replication cannot reuse one field value for free

Status: **EXACT_SCALING_GIVEN_STANDARD_PLANCHEREL_POLYA_BERNSTEIN_RIESZ**.

For `N^-1`-separated source centers, the standard Plancherel--Polya sampling theorem gives `sum_a |f(x_a)|^p <= C_PP N^3 ||f||_p^p` for a fixed band limit.  Combined with the affine source factor `||L^-1||||L||^2<=kappa^2 r_g`, every power of `N` cancels exactly.

Differentiated SGS source:

`sum_a rho_R,a^(3/2) <= C_samp (kappa^2 s0 C_D3)^(3/2) ||R||_(3/2)^(3/2)`.

The exact Germano increment estimate then routes the right-hand side to the global cubic velocity-increment charge at that filter scale.

Viscous-fourth source:

`sum_a rho_nu,a^2 <= C_samp,2 (kappa^2 s0 nu C_41)^2 d_V`,  `d_V=N^-1||grad V||_2^2`.

Thus many separated viscous source grains pay additive resolved dissipation.

For the strict filtered pressure,

`-Delta P=partial_i partial_j(V_i V_j+R_ij)`

and `P` remains band limited.  Riesz + Bernstein gives

`||P||_(3/2)^(3/2) <= sqrt(2) C_R^(3/2)[C_B^3 mu_V^(3/2)+||R||_(3/2)^(3/2)]`.

Hence separated pressure-third source grains route to resolved low-pass critical mass or the same cubic increment charge; pressure-third near-field is not an independent unpriced source.

Stress: `50000`
- worst scale-invariance residual: `6.262e-16`
- minimum SGS routing margin: `0.000e+00`
- minimum pressure routing margin: `0.000e+00`
- minimum viscous routing margin: `0.000e+00`
