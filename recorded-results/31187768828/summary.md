# Affine SGS boundary ledger

Status: **CERTIFIED**.

For an affine window transported by the strict low-pass affine jet, write
`W=U-V_aff`.  The nonviscous combined-work leakage is
`int grad chi . (e W + R U)`.  Pointwise Young gives

`.5 |U|^2 |W| + |R||U| <= (2/3)|U|^3 +(1/6)|W|^3 +(2/3)|R|^(3/2)`.

Hence large differential-advection/SGS-transport leakage forces the scale-critical annular charge
`int_A(|U|^3+|W|^3+|R|^(3/2)) >= |L_cubic|/||grad chi||`.
For the affine shell window this is at least `2M |L_cubic|/(3 N Cchi)`.
Pressure cancellation separately forces `S M/(3 N Cchi)` in `|U|^3+|P|^(3/2)`.

The resolved viscous boundary term is not a new source: on a parabolic lifetime and the smooth spectral support,
`|L_nu| <= (33/20)c nu Cchi M^-1 sup ||U||_2^2` for `delta<=1/20`.
Thus it renormalizes the existing `1/M` localization coefficient.

A quadratic spatial partition `sum eta_alpha=1` has exact total boundary cancellation:
`sum_alpha(e partial_t eta_alpha + grad eta_alpha.F)=0`.  Overlap can matter after selecting a lineage, but there is no global packet-count loss from the partition itself.

Stress checks: `50000`
- worst cubic Young ratio: `0.999999481`
- minimum cubic margin: `4.865e-06`
- worst partition cancellation residual: `5.997e-15`
- minimum clean viscous-bound margin: `9.481e-08`
- worst SGS-stress / increment-cubic bound ratio: `0.638619634`
- worst exact increment-identity residual: `1.349e-13`
