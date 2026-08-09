# Coherent averaged resolved-strain source calculus

Status: **EXACT_COHERENT_AVERAGED_RESOLVED_STRAIN_SOURCE_IDENTITY__REYNOLDS_CORRECTIONS_ROUTE_TO_CRITICAL_DISSIPATION**.

Let the Gaussian analysis eddy move by its own affine regression,

`Xdot=bar V`, `Ldot=bar A L`, `bar V=<V>_gamma`, `bar A=<grad V>_gamma`.

For any field `f`, fixed intrinsic `z` gives the exact moving-average identity

`d<f>/dt = <D_t^V f-r.grad f>`,

where `r=V-bar V-bar A(x-X)`.  Applying this to the resolved gradient equation yields

`dot bar A = -bar A^2 - <a^2> - <Hess P> - <grad div R_sgs> + nu <Delta A> - <r.grad A>`,

with `a=A-bar A`.  The corresponding corotational symmetric-strain identity is exact.

The two terms created only by coherent averaging are not new currencies.  In intrinsic variables `R=L^-1 r`, the affine projection gives `E[R tensor z]=0`, while both `V` and the affine regression are divergence free.  Gaussian integration by parts therefore gives

`<r.grad A> = <(z.R)(A-bar A)>`.

With `K_coh^2=E||L^-1(A-bar A)L||^2`,

`||<a^2>|| <= kappa^2 K_coh^2`,
`||<r.grad A>|| <= sqrt(7) kappa K_coh^2`.

Since `int N^-2 K_coh^2 dt <= C_var D_V`, the integrated normalized source weight of both Reynolds corrections is at most

`Sigma_Reynolds <= 1.18115356379 D_V`

on `kappa<=567/500`.  The averaged quadratic `bar A` source is also linearly bounded by `D_V`; together the local quadratic+averaging contribution is at most

`Sigma_local <= 2.0383176489 D_V`.

Thus changing from the point-sampled affine jet to the coherent whole-eddy jet does not open a new source ledger.  Its only new terms are critical dissipation.

The remaining filtered sources do not worsen under averaging: `||<Hess P>||<=||Hess P||_infty`, and similarly for differentiated SGS and viscosity.  Therefore the existing resolved pressure/mass-or-SGS, SGS/coherent-service and viscous-`D_V` collision thresholds apply unchanged.  These sources are not charged once before averaging and again after averaging; the average is the source used by this transporter.

Stress: `50000` source/reynolds/collision states
- worst corotational identity residual: `9.172e-16`
- minimum Reynolds bound margin: `1.789e-07`
- minimum scaled-weight margin: `2.261e-07`
- minimum quadratic-weight margin: `2.633e-07`
- clean Reynolds source coefficient: `1.18115356379`
- clean total local coefficient: `2.0383176489`

This closes the averaged-jet source calculus at the exact Gaussian-analysis level.  The remaining continuum step is no longer a missing source formula: it is the **assembly theorem** showing that every recursively selected efficient smooth-SGS block may use this coherent averaged transporter, one-shot near-Gaussian profile, physical-energy causal gate and exact coherent binary work measure with the same selected transfer normalization and only the already summable `Xi`.  No global-regularity claim is made.
