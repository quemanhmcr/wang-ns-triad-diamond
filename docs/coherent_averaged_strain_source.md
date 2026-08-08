# Coherent averaged resolved-strain source: whole-eddy Kelvin transport has no new source currency

The coherent affine projection identifies the natural common transporter of a
material Gaussian eddy:

\[
\bar V=\langle V\rangle_\gamma,
\qquad
\bar A=\langle\nabla V\rangle_\gamma.
\]

A remaining concern is whether replacing the point-sampled center jet by this
whole-eddy affine regression creates an uncontrolled source term.  It does not.
The moving-average identity is exact, and the two new Reynolds terms are
quadratic in the already identified coherent deformation variance.

## 1. Move the Gaussian measure with its own affine regression

Let

\[
x=X+Lz,
\qquad z\sim N(0,I_3),
\]

and choose

\[
\boxed{
\dot X=\bar V,
\qquad
\dot L=\bar A L.
}
\]

Since `V` is incompressible,

\[
\operatorname{tr}\bar A
=\langle\operatorname{div}V\rangle=0.
\]

For any smooth tensor field `f`, differentiation at fixed intrinsic `z` gives
exactly

\[
\boxed{
\frac d{dt}\langle f\rangle_\gamma
=
\left\langle
\partial_t f+ar U\cdot\nabla f
\right\rangle_\gamma,
}
\]

where

\[
\bar U(x)=\bar V+\bar A(x-X).
\]

Put

\[
r=V-\bar U.
\]

Then both `V` and `bar U` are divergence free, so `r` is divergence free, and

\[
\boxed{
\frac d{dt}\langle f\rangle_\gamma
=
\langle D_t^V f-r\cdot\nabla f\rangle_\gamma.
}
\]

No frozen packet or time-independent window is used.

## 2. Average the exact resolved gradient equation

For the strict transporter `V=S_(N/4)u`,

\[
D_t^V A
=-A^2-\nabla^2P-\nabla\nabla\cdot R_{SGS}
+\nu\Delta A,
\qquad A=\nabla V.
\]

Write

\[
\bar A=\langle A\rangle,
\qquad
a=A-\bar A.
\]

Since `E a=0`,

\[
\langle A^2\rangle
=\bar A^2+\langle a^2\rangle.
\]

Therefore the coherent averaged affine jet obeys the exact equation

\[
\boxed{
\dot{\bar A}
=-\bar A^2
-\langle a^2\rangle
-\langle\nabla^2P\rangle
-\langle\nabla\nabla\cdot R_{SGS}\rangle
+\nu\langle\Delta A\rangle
-\langle r\cdot\nabla A\rangle.
}
\]

This is the resolved Navier--Stokes source calculus for the whole coherent eddy.

## 3. Objective strain identity

Let

\[
\bar A=\bar S+\bar\Omega.
\]

Using the same corotating convention as the center resolved-strain theorem,

\[
\boxed{
\begin{aligned}
\mathring{\bar S}
={}&-\bar S^2-\bar\Omega^2+[\bar S,\bar\Omega]\\
&-\operatorname{sym}\langle a^2\rangle
-\langle\nabla^2P\rangle\\
&-\operatorname{sym}\langle\nabla\nabla\cdot R_{SGS}\rangle
+\nu\langle\Delta S\rangle\\
&-\operatorname{sym}\langle r\cdot\nabla A\rangle.
\end{aligned}}
\]

The filtered pressure, SGS and viscous sources are the same physical sources as
before, only averaged over the selected coherent eddy.  They must not be charged
again under new names.

## 4. The transport Reynolds term has an exact Gaussian form

In intrinsic coordinates

\[
W=L^{-1}V(X+Lz),
\qquad
R=L^{-1}r.
\]

The coherent affine projection gives

\[
\mathbb E[R\otimes z]=0.
\]

Moreover

\[
\operatorname{div}_zR=0.
\]

For the physical Gaussian density,

\[
\nabla_x\log\rho
=-L^{-T}z.
\]

Integration by parts therefore gives

\[
\begin{aligned}
\langle r\cdot\nabla A\rangle
&=\langle(z\cdot R)A\rangle\\
&=\boxed{
\langle(z\cdot R)(A-\bar A)\rangle,
}
\end{aligned}
\]

because `E[z dot R]=0`.

So the moving-average correction is not a derivative loss.  It is a covariance
of two already measured deformation fluctuations.

## 5. Both new terms are quadratic in coherent deformation variance

Recall

\[
\mathcal K_C^2
=
\mathbb E
\|L^{-1}(A-\bar A)L\|_F^2.
\]

If `kappa=cond(L)`, then

\[
\mathbb E\|a\|_F^2
\le\kappa^2\mathcal K_C^2.
\]

Hence

\[
\boxed{
\|\langle a^2\rangle\|_F
\le\kappa^2\mathcal K_C^2.
}
\]

The coherent affine theorem also gives

\[
\mathbb E|z|^2|R|^2\le7\mathcal K_C^2.
\]

Cauchy then yields

\[
\boxed{
\|\langle r\cdot\nabla A\rangle\|_F
\le\sqrt7\,\kappa\,\mathcal K_C^2.
}
\]

Thus the total source created solely by coherent averaging satisfies

\[
\boxed{
\|S_{Reynolds}\|
\le
(\kappa^2+\sqrt7\,\kappa)\mathcal K_C^2.
}
\]

There is no new independent source norm.

## 6. The entire new source weight is critical dissipation

The coherent deformation theorem gives

\[
\int N^{-2}\mathcal K_C^2dt
\le
C_{var}D_V,
\]

with

\[
C_{var}
=\kappa^2(2\pi)^{-3/2}(3/2)^3
\]

on the shell lower-axis bound.  Consequently

\[
\boxed{
\Sigma_{Reynolds}
\le
(\kappa^2+\sqrt7\kappa)C_{var}D_V.
}
\]

At `kappa<=567/500`, this is one fixed scale-independent coefficient multiplying
the **scale-critical** normalized dissipation `D_V`.

The averaged quadratic base term is also harmless.  Jensen and the Gaussian
density peak give

\[
\int N^{-2}|\bar A|^2dt
\le
(2\pi)^{-3/2}(3/2)^3D_V.
\]

Since the quadratic objective-strain source is bounded by `4|bar A|^2`, it too
is linear in `D_V` after source-weight integration.

These are critical currencies.  They are not uniform finite-count resets.

## 7. Averaged pressure, SGS and viscosity inherit the old collisions

Probability averaging cannot increase a global supremum:

\[
\|\langle\nabla^2P\rangle\|
\le\|\nabla^2P\|_\infty,
\]

and likewise for the filtered differentiated-SGS and viscous terms.  Therefore
the existing resolved-source bounds remain valid without an extra averaging
factor:

\[
\rho_P
\le\frac{\mu_V}{5700}
+\frac{\|R_{SGS}\|_{3/2}}{380},
\]

\[
\rho_{R,2}\le\frac{\|R_{SGS}\|_{3/2}}{380},
\]

and

\[
\rho_\nu\le\frac\nu{1500}\sqrt{\mathfrak d_V}.
\]

Thus averaged pressure routes to coherent resolved mass or SGS service, averaged
SGS routes to Germano/Onsager coherent ancestry, and averaged viscosity routes to
critical dissipation exactly as before.

## 8. Single-charge interpretation

Changing transporter from the center jet to the coherent averaged jet does not
create five new source names.  The causal roots are:

- filtered pressure: the existing pressure/SGS/mass root;
- differentiated SGS: the existing coherent-service root;
- viscosity: the existing critical-dissipation root;
- quadratic averaged strain: critical `D_V`;
- averaging Reynolds covariance: the **same coherent deformation event** already
  measured by `K_coh`, hence critical `D_V` when large.

The last two must not be separately counted as both “coherent deformation” and
“source variation.”  The first physical cause owns the charge.

## 9. What remains

The averaged affine transporter and its source calculus are now structurally
closed at the Gaussian analysis level.  The remaining continuum task is an
**assembly/extraction theorem**:

> every recursively selected efficient smooth-SGS block must enter the same
> coherent averaged transporter, one-shot near-Gaussian profile, physical-energy
> causal gate and exact coherent binary-work construction with the original
> selected physical transfer normalization and only the already summable `Xi`.

This note does not prove that final recursive assembly, and it makes no global
regularity claim.
