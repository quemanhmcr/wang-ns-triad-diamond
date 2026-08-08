# Coherent affine projection: all higher spatial deformation collapses to one variance

The resolved-role Egorov identity shows that the moving outer role is controlled
by spatial variation of the resolved deformation, while the center-jet shear
countermodel shows that one pointwise Hessian cannot see all of that variation.
The natural replacement is not an infinite list of `H4,H5,...` currencies.  It
is the orthogonal affine projection of the resolved flow on the **whole coherent
eddy**.

## 1. The coherent affine jet is a Gaussian regression

In intrinsic coordinates write

\[
W(z)=L^{-1}V(X+Lz),
\qquad z\sim\gamma=N(0,I_3).
\]

Define

\[
\boxed{
\bar v=\mathbb E_\gamma W,
\qquad
\bar A=\mathbb E_\gamma\nabla_z W.
}
\]

Gaussian integration by parts gives

\[
\boxed{
\bar A_{ij}=\mathbb E_\gamma[W_i z_j].
}
\]

Thus `vbar+Abar z` is exactly the orthogonal projection of `W` onto constant and
linear Gaussian Hermite modes.  Put

\[
R(z)=W(z)-\bar v-\bar A z.
\]

Then `R` has no Hermite degree zero or one.  This is a physical choice of gauge:
the common affine transporter is the best affine deformation seen by the entire
coherent eddy.

For incompressible `V`, `tr Abar=0`, so this averaged affine flow is still volume
preserving.

## 2. Ornstein--Uhlenbeck spectral gap removes the high-Hermite tower

Define coherent deformation variance

\[
\boxed{
\mathcal K_C^2
=\mathbb E_\gamma
\|\nabla_zW-\bar A\|_F^2.
}
\]

Since every Hermite component of `R` has degree at least two,

\[
\boxed{
\mathbb E_\gamma|R|^2
\le\frac12\mathcal K_C^2.
}
\]

There is also a clean weighted estimate.  In one coordinate the Gaussian
creation/annihilation identity gives

\[
\|z_jf\|_2^2
\le4\|\partial_jf\|_2^2+2\|f\|_2^2.
\]

Summing in three dimensions and using the spectral gap,

\[
\boxed{
\mathbb E_\gamma|z|^2|R|^2
\le7\mathcal K_C^2.
}
\]

These statements involve the **entire Hermite spectrum**.  No truncation at
`H3` is required.

## 3. Full Gaussian-core non-affine forcing

For the unchirped intrinsic Gaussian carrier

\[
\psi(z)=g(z)e^{iq\cdot z},
\qquad
\nabla_z\psi=(iq-z/2)\psi,
\]

the scalar non-affine advection is `R dot grad psi`.  Hence

\[
\frac{\|R\cdot\nabla\psi\|_2}{\|\psi\|_2}
\le
\left(\frac{|q|}{\sqrt2}+\frac{\sqrt7}{2}\right)\mathcal K_C.
\]

The vector rapid-distortion amplitude mismatch is
`(grad W-Abar) psi`, whose relative norm is at most `K_C`.  Therefore

\[
\boxed{
\frac{\|F_{nonaff}\|_2}{\|\psi\|_2}
\le
\left(1+\frac{|q|}{\sqrt2}+\frac{\sqrt7}{2}\right)\mathcal K_C.
}
\]

This is the desired spectrum collapse: all spatial non-affine Gaussian-core
forcing is controlled by one affine-invariant physical observable.

The earlier `H1/H3` analysis remains valuable because it gives sharper
transfer-facing structure for the lowest non-affine mode.  It is no longer
necessary to invent master currencies for each higher Hermite degree merely to
control the outer moving role.

## 4. Intrinsic carrier is bounded on the scale-matched branch

Let the physical Gaussian axes be `ell_i`,

\[
r_g=(\ell_1\ell_2\ell_3)^{1/3},
\qquad
\kappa=\ell_{max}/\ell_{min}.
\]

Elementary geometry gives

\[
\boxed{
\ell_{max}\le\kappa^{2/3}r_g.
}
\]

Thus for a selected carrier with `|k|/N<=R_k`,

\[
\boxed{
|q|=|L^Tk|
\le\kappa^{2/3}(Nr_g)R_k.
}
\]

On the transition-aspect, scale-matched branch `kappa<=567/500` and
`Nr_g<=s_0`, this is a fixed number.  Large radius is already an affine critical
mass/ancestry branch rather than a hidden failure of the forcing estimate.

## 5. Large coherent deformation is critical dissipation

The same observable also has a direct global collision.  Pointwise

\[
\mathcal K_C^2
\le
\mathbb E_\gamma
\|L^{-1}\nabla V(X+Lz)L\|_F^2
\le
\kappa^2\mathbb E_\gamma|\nabla V(X+Lz)|_F^2.
\]

The normalized physical Gaussian density has maximum

\[
(2\pi)^{-3/2}r_g^{-3}.
\]

The shell theorem gives every physical axis `>2/(3N)`, hence

\[
r_g>\frac2{3N}.
\]

Consequently

\[
\mathcal K_C^2
\le
\kappa^2(2\pi)^{-3/2}\frac{27}{8}N^3\|\nabla V\|_2^2.
\]

On `T=cN^{-2}`, Cauchy in time gives for

\[
I_K=\int_0^T\mathcal K_C(t)dt
\]

the scale-free collision

\[
\boxed{
I_K^2
\le C_{coh}\,cD_V,
\qquad
C_{coh}
=\kappa^2(2\pi)^{-3/2}\frac{27}{8}.
}
\]

On the full transition strip `kappa<=567/500`, this is one fixed universal
coefficient.  Therefore

\[
\boxed{
D_V\ge\frac{I_K^2}{C_{coh}c}.
}
\]

This is a **critical dissipation** currency.  As with the high-strain theorem,
it must not be promoted to a scale-independent finite reset count.

## 6. Physical dichotomy

The outer-role non-affinity now has a natural two-way interpretation:

\[
\boxed{
\text{small }I_K
\Longrightarrow
\text{the full Gaussian-core low--high residual is perturbative},
}
\]

while

\[
\boxed{
\text{large }I_K
\Longrightarrow
\text{definite critical }D_V.
}
\]

Thus the center-flat cubic shear is not an escape.  It simply has nonzero
coherent deformation variance; if accumulated strongly it pays dissipation, and
if accumulated weakly it cannot strongly move the Gaussian core away from its
best coherent affine transport.

## 7. New remaining bridge

Using `Abar` as the common affine gauge is more natural than returning to the
center jet, but it changes the source question.  The existing resolved objective
strain identity was derived for

\[
A(X)=\nabla V(X).
\]

The next theorem should derive the **resolved Navier--Stokes evolution of the
coherent averaged affine jet** `Abar(t)`, including the motion of the Gaussian
analysis measure, and show that its corotational variation routes once into the
same filtered pressure, SGS, viscous, coherent service and critical-dissipation
currencies.

That averaged-transporter source calculus is not asserted here.  No global
regularity claim is made.
