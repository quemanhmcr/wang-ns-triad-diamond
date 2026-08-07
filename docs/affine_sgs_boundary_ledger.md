# Affine SGS boundary residuals: differential transport, `RU`, viscosity and partition cancellation

Use the pressure-free **combined-work** localized energy identity

\[
\int\!\!\int\chi G
+\nu\int\!\!\int\chi|\nabla U|^2
=E_\chi(t_0)-E_\chi(t_1)+\widetilde L_\chi,
\]

\[
\widetilde L_\chi
=\int\!\!\int e\,\partial_t\chi
+\int\!\!\int\nabla\chi\cdot(eU+RU-\nu\nabla e).
\]

Let the affine window be transported by the strict low-frequency jet
`V_aff`.  Then `partial_t chi+V_aff.grad chi=0` and

\[
\boxed{
\widetilde L_\chi
=\int\!\!\int\nabla\chi\cdot(eW+RU-\nu\nabla e),
\qquad W=U-V_{aff}.
}
\]

This displays the remaining physical leakage channels without double counting.

## 1. Differential resolved advection and `RU` force a critical annular charge

Pointwise,

\[
\frac12|U|^2|W|+|R|_F|U|
\le
\frac23|U|^3+rac16|W|^3+rac23|R|_F^{3/2}.
\]

Therefore, for `A=supp grad chi`, if

\[
L_{cub}=\int_A\nabla\chi\cdot(eW+RU),
\]

then

\[
\boxed{
\int_A\left(|U|^3+|W|^3+|R|_F^{3/2}\right)
\ge \frac{|L_{cub}|}{\|\nabla\chi\|_\infty}.
}
\]

Every term is Navier--Stokes scale critical.  Thus a large `RU` or differential
transport leakage is not an uncontrolled error: it creates a concrete annular
critical charge.

For the affine grain window,

\[
\|\nabla\chi\|_\infty\le\frac{3NC_\chi}{2M},
\]

hence

\[
\boxed{
\int_A\left(|U|^3+|W|^3+|R|^{3/2}\right)
\ge\frac{2M}{3NC_\chi}|L_{cub}|.
}
\]

The existing raw-SGS pressure-cancellation branch similarly becomes

\[
\boxed{
\int_A(|U|^3+|P|^{3/2})
\ge \frac{SM}{3NC_\chi}.
}
\]

## 2. Resolved viscous boundary flux renormalizes the `1/M` ledger

The combined-work leakage contains

\[
L_\nu=-\nu\int\!\!\int\nabla\chi\cdot\nabla e.
\]

For the smooth low-pass field `U`, the filter support is contained in
`|xi|<=e^delta N`; therefore

\[
\|\nabla U\|_2\le e^\delta N\|U\|_2.
\]

On a lifetime `T=cN^{-2}`,

\[
|L_\nu|
\le
\frac32e^\delta\frac{c\nu C_\chi}{M}
\sup_t\|U(t)\|_2^2.
\]

For the certified smooth transition `delta<=1/20`, Arb gives
`e^(1/20)<11/10`, hence

\[
\boxed{
|L_\nu|
\le
\frac{33}{20}\frac{c\nu C_\chi}{M}
\sup_t\|U(t)\|_2^2.
}
\]

Thus resolved viscous boundary transport introduces no new scale: it adds to the
same `a/M` coefficient already used by the filter/window commutator.

This is distinct from **bulk viscosity of an individual Gaussian role**, which is
already tangent to its Gaussian manifold.
\n\n## SGS stress is a cubic velocity-increment charge\n\nThe `R` term can be returned to a primitive velocity observable.  For every\nnormalized convolution filter `G_l` (`int G_l=1`), let\n\n\[\n\delta_r u(x)=u(x-r)-u(x).\n\]\n\nThen the SGS stress has the exact increment representation\n\n\[\n\boxed{\nR(x)=\int G_l(r)\,\delta_ru\otimes\delta_ru\,dr\n-\left(\int G_l(r)\delta_ru\,dr\right)^{\!\otimes2}.\n}\n\]\n\nThis identity is algebraic and does not require a positive filter.  Put\n`g_1=||G||_1`.  With\n\n\[\nA_2=\int |G_l(r)|\,|\delta_ru|^2dr,\n\]\n\nCauchy gives `|int G_l delta u|^2 <= g_1 A_2`, hence\n\n\[\n|R|_F\le(1+g_1)A_2.\n\]\n\nA second Holder/Jensen step gives\n\n\[\nA_2^{3/2}\le g_1^{1/2}\int |G_l(r)|\,|\delta_ru|^3dr.\n\]\n\nTherefore\n\n\[\n\boxed{\n|R(x)|_F^{3/2}\n\le (1+g_1)^{3/2}g_1^{1/2}\n\int |G_l(r)|\,|\delta_ru(x)|^3dr.\n}\n\]\n\nThus a large `RU` boundary branch forces an **Onsager-type cubic increment\ncharge at the actual SGS filter scale**.  This is a physical velocity observable,\nnot an independent stress-tensor mystery.  The remaining theorem is to convert\na persistent increment charge into fresh/reused affine grains without assuming\na global `L^3` bound.\n
## 3. A quadratic partition has no global boundary-count loss

Let nonnegative weights `eta_alpha` satisfy

\[
\sum_\alpha\eta_\alpha(x,t)=1.
\]

Then

\[
\sum_\alpha\partial_t\eta_\alpha=0,
\qquad
\sum_\alpha\nabla\eta_\alpha=0.
\]

For any energy density `e` and flux `F`, therefore

\[
\boxed{
\sum_\alpha
\left(e\,\partial_t\eta_\alpha+\nabla\eta_\alpha\cdot F\right)=0.
}
\]

Taking `eta_alpha=chi_alpha^2` gives the natural energy partition.  Spatial
overlap is not, by itself, a global Bellman/packet-count penalty.  It becomes a
real issue only after one selects a lineage/subfamily or compares different
frequency-role partitions, where the uncancelled interfaces must be entered into
the cross-error ledger.

## 4. What this closes and what it does not

This note closes the **algebraic and scale-critical classification** of the
resolved window residuals:

- large differential-advection/`RU` leakage -> critical annular charge;
- pressure cancellation -> critical pressure/velocity annular charge;
- viscous boundary -> explicit `O(1/M)` energy term;
- complete quadratic partition -> exact global cancellation.

It does not yet prove that the `R` charge or a selected-subfamily interface charge
always produces a fresh affine Gaussian grain with a uniform constant.  That is
now an ancestry/extraction problem rather than an unidentified PDE term.
