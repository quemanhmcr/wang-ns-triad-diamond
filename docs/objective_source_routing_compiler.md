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

## Pressure is pair-native, not reservoir-native

The coarse resolved estimate

\[
\rho_P\le \frac{\mu_V}{5700}+\frac{\|R\|_{3/2}}{380}
\]

remains correct and is retained as a **diagnostic API**.  It is no longer the canonical compiler state, because `mu_V=N||V||_2^2` erases the elliptic derivative suppression carried by the actual pressure interaction.

The canonical route scalarizes the actual coherent-averaged pressure Hessian by its event Frobenius dual and uses the exact positive source cover

\[
\boxed{
\rho_P\le[r_{SGS}]_+ + \sum_{a\le b}[p_{ab}]_+.
}
\]

Thus either positive SGS pressure source or the resolved unordered hard-pair law carries at least `Sigma_P/2`; exact ties remain joint.

If the SGS source law carries actual integrated weight `R_SGS`, the direct order-two source bound gives

\[
\int\|R\|_{3/2}\,d\tau\ge380R_{SGS},
\]

so **that actual weight**, not merely the clean half-threshold, enters the same linear coherent-service map used by the objective SGS owner.

If the resolved pair law is an owner, the pressure-pair theorem already gives

\[
\boxed{
\mu_{child}e^{H_2^P}\ge320\frac{\Sigma_P}{c}.
}
\]

The strict low-pass contraction transfers the resolved `V`-shell lower to the same hard shell of `u`, so every positive resolved pair law supplies an input to the generic critical-shell first-stop theorem.  There is no separate diffuse-pressure fate: large `H_2^P` weakens the shell threshold but does not change its owner.

On a full no-hit natural shell corridor, the generic shell theorem therefore gives a conditional own-scale service lower.  Because that theorem is linear in the shell critical mass, the pressure conjugacy survives composition exactly:

\[
\boxed{
e^{H_2^P}Y_{shell}
\ge Y_{shell}^{generic}\!\left(320\frac{\Sigma_P}{c}
ight),
}
\]

and for the full-natural integrated service,

\[
\boxed{
e^{H_2^P}S_{shell}
\ge S_{shell}^{generic}\!\left(320\frac{\Sigma_P}{c}
ight).
}
\]

The compiler records these as **conditional full-survivor** bounds; it does **not** claim service before the first-stop corridor is resolved.

The quarter split remains only a readable corollary:

\[
q_{max}\ge1/4\Rightarrow\mu_{child}\ge80\Sigma_P/c,
\qquad
q_{max}\le1/4\Rightarrow H_2^P\ge\log4.
\]

Neither side creates a new stop class.  Pressure-source `H_2^P` is not a Shannon/Renyi child-energy probability.

The derivative-correct fixed-material-pair `<1/5` lifetime remains useful **after** material sidecars are attached as a reuse/capacity refinement, but it is no longer the pressure renewal entrance.  The old `mu_V` split likewise remains available for diagnostics and compatibility only.

## Master-facing semantics

The compiler creates no new currency:

- local coherent source -> resolved `D_V` -> generic critical-shell recursion;
- viscosity -> resolved `D_V` -> generic critical-shell recursion;
- SGS -> coherent service -> high-frequency owner / old-pool erosion / `Xi` / fresh shell / entropy-cycle;
- pressure -> actual positive SGS coherent service or entropy-weighted critical shell;
- pure material sidecars remain sidecars; genuine role change remains interface/relink.

No branch above is promoted to an additive finite reset.  `D_V`, critical shell mass, high-frequency dissipation and source service are scale-sensitive recursive quantities unless a separate globally bounded scale-independent theorem says otherwise.  Aggregate pressure `mu_V` is not a recursive state at all in the canonical route.

## Scope

This theorem is an assembly theorem for source ownership.  It deliberately does not prove supplier-specific signed-good scale progress and does not identify high-frequency enstrophy with resolved low-pass dissipation.  Pressure shell renewal is now supplied directly by the pair theorem rather than by converting aggregate low-pass mass.

The remaining continuum work is the final master assembly/topology that composes these native owners with physical positive HH causality, common-slice propagation, service/flat erosion and the existing entropy/Hodge machinery without double charging.

No 3D Navier--Stokes global-regularity conclusion is asserted.

## High-frequency SGS service is handed to physical tail energy, not to `D_V`

The `D_high` exit of coherent increment service remains a smooth-LP normalized-enstrophy observable.  It is **not** identified with an orthogonal hard tail.  We fix the same square-normalized smooth LP analysis--synthesis frame used by coherent service; its lower high-band support is `M_j/2` and its `L^2` square-Bessel constant is one.  Hence the companion theorem has the canonical certified comparison

\[
D_{tail}:=N\int\|\nabla P_{>N}u\|_2^2dt\ge c_{LP}D_{high}.
\]

Thus `c_LP=1/4` already for the canonical smooth LP frame; hard orthogonal dyadic annuli independently exhibit the same one-quarter spectral lower.  Alternative LP frames must carry their own fixed comparison constant.  The hard-tail Navier--Stokes energy law then gives

\[
N\|P_{>N}u(s)\|_2^2+NW_>^+\ge2\nu D_{tail}.
\]

Hence either inherited high-tail energy or actual positive nonlinear regeneration work carries at least `nu D_tail`, and therefore at least `nu c_LP D_high`.  The inherited owner exposes a genuine hard shell with the same critical-mass lower.  The regeneration owner disintegrates into actual positive high-shell work and, because low--low is support-excluded on those shells, into HH or resolved cross/interface work.

Thus the compiler's high-frequency branch now has the explicit physical owner interface

`SGS source -> coherent service -> smooth-LP D_high -> physical D_tail -> inherited critical shell OR positive HH/interface regeneration`.

This still creates no finite reset and does not make a large HH work atom satisfy the generated-energy productivity gate automatically.
