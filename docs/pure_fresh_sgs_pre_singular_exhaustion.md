# Pure fresh-SGS recurrence is exhausted on every compact pre-singular interval

## Status

**DRAFT_PURE_FRESH_SGS_PRE_SINGULAR_EXHAUSTION__WEIGHTED_SCALE_MOMENT_FORBIDS_ARBITRARILY_HIGH_DOUBLING__FIXED_SMOOTH_FILTER_ENERGY_ENVELOPE_GIVES_POSITIVE_FIRST_STOP_TIME__NO_ENTROPY_COST_NO_FINITE_RESET_NO_GLOBAL_H1_ASSUMPTION**

This note attacks the hard survivor left after `fresh_service_scale_reentry.md` and the strictly-descending companion theorem.  It uses only the already selected positive fresh service law, the fixed smooth resolved cutoff, the Navier--Stokes energy bound, and regularity on a **fixed compact interval strictly before a hypothetical singular time**.  It does not turn scale entropy, shell mass, source weight, or elapsed analysis horizon into a synthetic reset.

## 1. The fresh service law has a weighted scale moment

Write the canonical square-normalized annular bands as

\[
M_j=2^jN,\qquad j\le0,
\]

and let `F_j` be the actual integrated fresh NN service pushed to band `j`.  The certified increment-energy bound is

\[
F_j\le 4\int M_j\|u_j\|_2^2\,d\tau.
\]

Multiply by `2^{-j}` and sum.  Square normalization gives

\[
\boxed{
\sum_{j\le0}2^{-j}F_j
\le 4N\int\sum_j\|u_j\|_2^2d\tau
\le 4cNE_{\rm global}.
}
\]

On a fresh owner `F=sum F_j>=Y/4`.  Therefore its deterministic normalized scale law `p_j=F_j/F` obeys

\[
\boxed{
\sum_{j\le0}2^{-j}p_j
\le K_N:=\frac{16cNE_{\rm global}}{Y}.
}
\]

This normalized law is still **not** a child-energy causal probability and its `H_inf` is not a cost.

Put `n=-j>=0` and

\[
L_N=\lceil\log_2(2K_N)\rceil.
\]

Markov gives

\[
\Pr(n\ge L_N)\le K_N2^{-L_N}\le\frac12.
\]

Hence at least half the actual fresh service lies on the first `L_N` bands, and one actual band satisfies

\[
\boxed{p_{\max}\ge\frac1{2L_N}.}
\]

The content is a service-moment consequence of finite energy.  No entropy has been scalarized.

## 2. A fresh doubling cannot survive to arbitrarily high scale

The only way the fresh hard-shell route can increase scale is the selected `j=0` upper shell

\[
N_{\rm next}=2N.
\]

Because the selected band is the max atom of the actual scale law, the certified two-shell theorem gives

\[
\mu_{2N}
\ge p_{\max}\frac{Y}{24c}
\ge \frac{Y}{48cL_N}.
\]

Fix a starting event time `t_bar<T_*` in a hypothetical maximal classical solution.  On `[0,t_bar]`, ordinary classical regularity gives

\[
H_*:=\sup_{0\le t\le t_{bar}}\|\nabla u(t)\|_2^2<\infty.
\]

For the exact doubling shell `{N<|xi|<=2N}`,

\[
\mu_{2N}
=2N\|P_{N<|D|\le2N}u\|_2^2
\le \frac{2H_*}{N}.
\]

Thus a fresh doubling requires simultaneously

\[
\frac{Y_*}{48c_*L_N}
\le \frac{2H_*}{N}.
\]

But `L_N=O(log N)`.  Therefore this fails at a finite scale.  **Fresh doubling has a finite parent-frequency ceiling on every compact pre-singular interval.**  This does not assume an `H^1` bound uniform as `t_bar upward T_*`.

Since every fresh child has frequency at most `2N`, once high-scale doubling is excluded an arbitrary pure-fresh word — even one alternating downward, same-scale and upward moves — remains in a bounded frequency range.

## 3. The fixed filtered NS source gives a derived positive event duration

Let

\[
R=S_{N/4}(u\otimes u)-V\otimes V,
\qquad V=S_{N/4}u.
\]

The chosen strict cutoff is smooth and fixed.  Write `C_S=||K^S||_1` for its unit-scale convolution-kernel norm.  Since `|S_{N/4}(xi)|<=1`,

\[
\|R\|_1
\le (C_S+1)E_{\rm global}.
\]

Also `supp Rhat subset B_(N/2)`.  Choose once and for all a smooth reproducing multiplier equal to one on `B_(1/2)` and let `C_rep` be the `L^(3/2)` norm of its unit-scale kernel.  Scaling and Young give

\[
\|R\|_{3/2}
\le C_{rep}N\|R\|_1.
\]

The already-certified order-two Fourier estimate on `B_(N/2)` gives

\[
\rho_R:=N^{-4}\|\operatorname{sym}\nabla\nabla\!\cdot R\|_\infty
\le C_2\|R\|_{3/2},
\]

hence

\[
\boxed{
\rho_R
\le C_2C_{rep}(C_S+1)NE_{\rm global}.
}
\]

These are fixed analysis constants of the smooth resolved representation, not new PDE currencies.

For one actual SGS-source first-stop slab,

\[
\Sigma_R=\int\rho_R\,d\tau,
\qquad d\tau=N^2dt.
\]

If the native source face supplies the uniform positive owner floor `Sigma_R>=sigma_*`, then

\[
\boxed{
\Delta t
\ge
\frac{\sigma_*}
{C_2C_{rep}(C_S+1)N^3E_{\rm global}}.
}
\]

This is **not an assumed event clock**.  It is a consequence of the pointwise filtered Navier--Stokes source bound.

## 4. Exhaustion of an eventually-pure fresh-SGS word

Take a recursively ordered fresh-SGS word backwards from a fixed pre-singular event.  Consecutive first-stop slabs meet at their physical endpoints and otherwise do not overlap.

The doubling argument gives a finite scale ceiling `N_max`.  Therefore every event in the word obeys the common lower

\[
\Delta t\ge
\delta_*:=
\frac{\sigma_*}
{C_2C_{rep}(C_S+1)N_{max}^3E_{\rm global}}>0.
\]

If the entire word lies in a finite physical span `T_span`, then

\[
\boxed{
\#\{\text{fresh SGS first stops}\}
\le \left\lfloor\frac{T_{span}}{\delta_*}\right\rfloor.
}
\]

Thus an eventually-pure fresh-SGS source tail cannot contain infinitely many genuine recursive first stops on a fixed compact pre-singular interval.

## 5. What this removes from the mixed-owner frontier

Together with the existing source compiler and pressure-pair telescope, this removes another free source-only escape:

- local coherent / averaged-Reynolds source already delegates to critical `D_V`;
- viscous source delegates to critical `D_V`;
- resolved pressure-pair epochs have their own descending telescope;
- high-frequency SGS service delegates to physical hard-tail ownership;
- selected-interface `Xi` is a transfer/interface ledger, not a new recursive source currency;
- **fresh SGS itself now has no infinite eventually-pure first-stop tail on a compact pre-singular interval.**

What remains is the genuinely mixed cross-family problem: source descendants can still alternate with strain/dissipation and actual HH/high-tail work, and any old-reservoir/capacity route must be kept in its correct inherited/reuse ledger rather than promoted to a causal reset.  This note makes no global-regularity claim.

## Native source binding

The executable native entry is `pure_fresh_sgs_native_source_binding.py`.  It refuses an arbitrary positive source scalar: the `FreshSGSScaleOwnerCertificate` must satisfy the certified compiler identity

`Y = objective_sgs_integrated_square_service_lower(Sigma_R) = C_Y Sigma_R`

for the same typed SGS owner before the pre-singular exhaustion theorem is called.  A fixed objective first-stop action floor `A_*` and slab bound `c<=c_*` give the uniform pair of floors

`Sigma_* = A_*/(4c_*)`,  `Y_* = C_Y Sigma_*`.

Thus the positive time lower is attached to the same local NS source law that generated the fresh service; it cannot be manufactured by attaching an unrelated scalar after the fact.
