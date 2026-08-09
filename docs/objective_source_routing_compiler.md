# Coherent objective-source owner compiler

## Why this is a compiler, not a new packet theorem

The coherent service-or-flat gate already isolates a large **objective strain variation action**.  The exact coherent averaged transporter equation already identifies the physical terms that can produce that variation.  Nothing in the PDE asks us to construct another persistent packet merely to synchronize those terms.

The correct operation is therefore to group source terms by **physical owner** and hand each owner to the theorem that already controls it.

For the coherent averaged strain equation, use four owner classes:

1. `local_DV`: the averaged quadratic term together with the two coherent Reynolds/transport corrections;
2. `pressure`: the averaged filtered pressure Hessian;
3. `SGS`: the averaged differentiated filtered SGS stress;
4. `viscosity`: the averaged viscous term.

If `A_obj` is the objective variation action on a scaled lifetime `c`, positivity and the triangle inequality give

\[
\Sigma_{local}+\Sigma_P+\Sigma_R+\Sigma_\nu\ge \frac{A_{obj}}{c}.
\]

Hence at least one owner has

\[
\Sigma_*\ge \frac{A_{obj}}{4c}.
\]

All owners meeting the threshold are retained as a joint set.  There is no lexicographic theorem priority.

## Local coherent source is already resolved dissipation

The coherent averaged source theorem proves

\[
\Sigma_{local}\le C_{local}D_V,
\]

where `C_local` is the already-certified sum of the averaged quadratic and coherent Reynolds coefficients.  Therefore

\[
D_V\ge \frac{\Sigma_{local}}{C_{local}}.
\]

This immediately enters the generic critical-shell theorem with

\[
\mu_0=\frac{D_0}{c}.
\]

The resulting shell is a recursive scale-critical supplier.  It is not a globally bounded additive reset.

## Viscosity is even more direct

For the resolved filtered transporter,

\[
\rho_\nu\le \frac{\nu}{1500}\sqrt{d_V}.
\]

Thus if `int rho_nu >= Sigma_nu` on a scaled interval of length `c`, Cauchy gives

\[
D_V\ge \frac1c\left(\frac{1500\Sigma_\nu}{\nu}\right)^2.
\]

This enters the same generic critical-shell reentry.  Temporal concentration only increases the dissipation price.

## The clean objective-SGS cancellation

The resolved objective source collision gives

\[
\|R\|_{3/2}\ge 380\rho_R.
\]

For a normalized filter with `g1=||G||_1`, the exact Germano estimate is

\[
\|R\|_{3/2}^{3/2}
\le (1+g_1)^{3/2}g_1^{1/2}Q,
\]

where `Q` is the cubic velocity-increment charge.  The coherent square-service theorem then gives

\[
Y=\frac{(Q/g_1)^{2/3}}{(C_{LP}C_B)^2}.
\]

The powers `3/2` and `2/3` cancel **exactly**:

\[
\boxed{
Y\ge C_Y\rho_R,
\qquad
C_Y=\frac{380}{g_1(1+g_1)(C_{LP}C_B)^2}.
}
\]

This route contains no affine radius and no temporal-persistence assumption.  It belongs directly to the filtered objective source, rather than to an auxiliary packet geometry.

If the integrated SGS owner weight is `Sigma_R`, then

\[
Y_{tot}\ge C_Y\Sigma_R.
\]

The existing high/low estimate is linear under time integration:

\[
S_{low}\ge Y_{tot}-2D_{high}.
\]

Therefore either

\[
D_{high}\ge \frac14Y_{tot},
\]

or

\[
S_{low}\ge \frac12Y_{tot}.
\]

On the low branch, once old integrated capacity is at most `Y_tot/8`, either selected-interface service is at least `Y_tot/8`, or fresh new-new service is at least `Y_tot/4`.  With the existing quarter-dominance threshold, a dominant fresh edge gives pair critical-mass occupation at least `Y_tot/32`, whole-shell occupation at least `Y_tot/64`, and therefore a pointwise whole-shell event at least

\[
\frac{Y_{tot}}{64c}.
\]

That event enters the generic critical-shell theorem.  If no fresh edge dominates, the existing collision chain gives atomic entropy `>=log 4`, ancestry entropy `>=log 2`, or same-ancestry pair/cycle mass `>=1/4`.

Crucially, the high-frequency dissipation branch is kept under its own owner.  It is **not** silently renamed resolved `D_V`.

## Pressure keeps its real two-way structure

For the strict filtered pressure,

\[
\rho_P\le \frac{\mu_V}{5700}+\frac{\|R\|_{3/2}}{380}.
\]

If an integrated pressure owner has weight `Sigma_P`, then positivity gives the exact alternative

\[
\boxed{
\int \mu_V\,d\tau\ge 2850\Sigma_P
\quad\text{or}\quad
\int\|R\|_{3/2}\,d\tau\ge190\Sigma_P.
}
\]

The stress alternative is precisely an effective objective-SGS source weight `Sigma_P/2`, so it enters the same linear coherent-service route above.

The resolved-mass alternative is intentionally **not** fed to the generic critical-shell theorem.  `mu_V=N||V||_2^2` is low-pass mass occupation and may sit arbitrarily far below the current block scale.  Calling it a critical shell would be a false geometric assertion.

Instead it remains the pressure-reservoir owner.  On a supplied signed-good low-strain lineage, every fixed materially reused low-low pair has pressure-service coefficient contraction

\[
\left(\frac{21}{20}\right)^4\left(\frac58\right)^3
=\frac{194481}{655360}<\frac13,
\]

so its total future capacity is less than `3/2` times its generation-zero coefficient.  Persistent pressure service must therefore relink to new reservoir pairs, fragment into atomic/component entropy or cycle structure, leave the low-strain reuse regime, or use the SGS part.

The material-label sidecar quotient is compatible with this statement: a selected reservoir family may change and pay its Moyal/ancestry charge without automatically destroying an unchanged smooth carrier.  A genuine role/probe change still delegates to the physical interface theorem.

## Master-facing semantics

The compiler creates no new currency:

- local coherent source -> resolved `D_V` -> generic critical-shell recursion;
- viscosity -> resolved `D_V` -> generic critical-shell recursion;
- SGS -> coherent service -> high-frequency owner / old-pool erosion / `Xi` / fresh shell / entropy-cycle;
- pressure -> SGS service or pressure-reservoir occupation;
- pure material sidecars remain sidecars; genuine role change remains interface/relink.

No branch above is promoted to an additive finite reset.  `D_V`, critical shell mass, high-frequency dissipation, source service and pressure reservoir occupation are all scale-sensitive recursive currencies unless a separate globally bounded scale-independent theorem says otherwise.

## Scope

This theorem is an assembly theorem for source ownership.  It deliberately does not prove supplier-specific signed-good scale progress, does not convert low-pass pressure mass into a shell, and does not identify high-frequency enstrophy with resolved low-pass dissipation.

The remaining continuum work is the final master assembly/topology that composes these native owners with physical positive HH causality, common-slice propagation, service/flat erosion and the existing entropy/Hodge machinery without double charging.

No 3D Navier--Stokes global-regularity conclusion is asserted.

## High-frequency SGS service is handed to physical tail energy, not to `D_V`

The `D_high` exit of coherent increment service remains a standard smooth-LP normalized-enstrophy observable.  It is **not** identified with an orthogonal hard tail.  The companion theorem requires the chosen LP partition to supply a fixed certified comparison

\[
D_{tail}:=N\int\|\nabla P_{>N}u\|_2^2dt\ge c_{LP}D_{high}.
\]

For hard orthogonal dyadic annuli the spectral lower is exactly `c_LP=1/4`; a smooth LP implementation carries its own fixed constant.  The hard-tail Navier--Stokes energy law then gives

\[
N\|P_{>N}u(s)\|_2^2+NW_>^+\ge2\nu D_{tail}.
\]

Hence either inherited high-tail energy or actual positive nonlinear regeneration work carries at least `nu D_tail`, and therefore at least `nu c_LP D_high`.  The inherited owner exposes a genuine hard shell with the same critical-mass lower.  The regeneration owner disintegrates into actual positive high-shell work and, because low--low is support-excluded on those shells, into HH or resolved cross/interface work.

Thus the compiler's high-frequency branch now has the explicit physical owner interface

`SGS source -> coherent service -> smooth-LP D_high -> physical D_tail -> inherited critical shell OR positive HH/interface regeneration`.

This still creates no finite reset and does not make a large HH work atom satisfy the generated-energy productivity gate automatically.
