# High-frequency service dissipation -> inherited critical shell or physical regeneration

Status: **EXACT_SQUARE_LP_HIGH_SERVICE_TO_PHYSICAL_TAIL_ENERGY__INHERITED_CRITICAL_SHELL_OR_PHYSICAL_REGENERATION__NO_RESOLVED_DV_RELABEL**.

The physical theorem is formulated first in the PDE's own currency

`D_tail=N int ||grad P_>N u||_2^2 dt`.

This orthogonal hard-tail dissipation is **not** resolved low-pass `D_V`, and it is not silently identified with the smooth LP `d_high` used by coherent increment service.  For high LP analysis multipliers supported above `a M_j` with square-Bessel upper `B`, Plancherel gives

`D_tail >= c_LP D_high`,  `c_LP=a^2/B`.

The canonical smooth square-normalized analysis--synthesis frame registered upstream has `a=1/2`, `B=1`, hence `c_LP=1/4`.  Independently, for the auxiliary hard annuli `M_j=2^jN`, `M_j/2<|xi|<=M_j`, the same one-quarter spectral lower holds.  Writing `mu_j=M_j||P_j u||_2^2`, the hard annular currency is `D_>^hard=int sum_j2^j mu_j d tau`.

Indeed `D_high` alone cannot force a critical shell: placing all currency at level `j` with `mu_j=2^-j D_high` keeps `2^j mu_j=D_high` while `mu_j->0`.

The hard-annulus comparison gives

`D_>^hard/4 <= D_tail <= D_>^hard`.

The missing physics is then viscosity plus the exact hard-tail energy law.  With `w=P_>N u` and actual positive nonlinear tail work

`W_>^+=int 2[Re<w,-P_>N P div(u tensor u)>]_+ dt`,

the Navier--Stokes energy identity gives

`N||w(s)||_2^2 + N W_>^+ >= 2 nu D_tail`.

Therefore at least one native owner carries `nu D_tail`.  With an LP supplier this is at least `nu c_LP D_high`; for the hard-annulus lower it is `nu D_>^hard/4` (exact ties remain joint):

1. **Inherited tail energy.**  Since

   `N||w(s)||_2^2=sum_j 2^-j mu_j(s)`

   and `sum_j2^-j=1`, some actual high shell `M_j>=2N` obeys

   `M_j||P_j u(s)||_2^2 >= nu D_tail`,

   hence at least `nu c_LP D_high` for an LP supplier.

   This is a genuine critical-shell seed and enters the existing generic shell first-stop/service theorem.  That theorem's observed-history guard is unchanged: the present theorem does not manufacture a full natural survivor.

2. **Positive nonlinear regeneration.**  Orthogonality disintegrates tail work into hard shell positive works `W_j^+` with `sum W_j^+>=W_>^+`.  Because every `M_j/N>=2`,

   `sum_j M_j W_j^+ >= 2N W_>^+`.

   At shell `M_j`, choose `V=S_(M_j/4)u`.  The hard shell is strictly above `M_j/2`, whereas `B(V,V)` is supported at or below `M_j/2`; low--low work is absent.  Expanding `u=V+h`, the signed shell work is exactly HH work plus resolved mixed/cross-interface work.  Hence positive parts satisfy

   `W_shell^+ <= W_HH^+ + W_interface^+`.

   On the clean regeneration branch the aggregate own-scale shell work is at least `2 nu D_tail`, so HH or resolved-interface positive work carries at least `nu D_tail` (at least `nu c_LP D_high` for the LP supplier).

The last HH statement is **not** the generated-energy gate.  Actual positive HH work still must pass its own energy comparison before the physical KL productivity theorem may be invoked.  Likewise interface work remains interface/strain provenance rather than being declared free.

Thus the anonymous high-frequency enstrophy exit has been replaced by native physical owners:

`smooth-LP D_high -> physical D_tail -> inherited critical shell OR actual positive HH/interface regeneration`.

No `D_high -> D_V` relabel, no packet persistence, no additive reset.

Stress: `50000` high-tail/shell/work states
- minimum hard-tail gradient lower margin: `7.821e-08`
- minimum hard-tail gradient upper margin: `7.351e-09`
- minimum two-owner energy margin: `2.946e-07`
- minimum inherited-shell margin: `2.020e-08`
- minimum shell-work disintegration margin: `5.521e-08`
- minimum HH/interface owner margin: `1.806e-07`
- maximum sampled joint energy-owner count: `2`
- maximum sampled joint regeneration-owner count: `2`

This theorem closes the unit mismatch in the coherent-service `D_high` branch at the physical-energy level.  Supplier-specific continuation of a positive regeneration owner and the low-frequency pressure-reservoir lineage remain separate master-facing questions.  No Navier--Stokes global-regularity conclusion is asserted.
