# Resolved-role Egorov identity: the missing defect is coherent deformation variance

The recursive coherent-work theorem shows that binary causality is already
inside the quadratic Navier--Stokes work once the outer selected roles exist.
The remaining structural question is therefore the moving **outer role** under
the strict resolved low--high flow.

A fixed Fourier cell is not physical: resolved advection may move energy across
that cell without generating a new high-frequency packet.  Common material
transport must be quotiented.  The affine case is already exact.  This note
identifies the full non-affine remainder and, importantly, gives a countermodel
to a tempting but false center-Hessian closure.

## 1. Exact affine-subtracted commutator

Let `Q_t` be a smooth scalar Fourier multiplier with convolution kernel `K_N`.
For a divergence-free resolved transporter `V`,

\[
[V\cdot\nabla,Q]f(x)
=
\int K_N(y)[V(x)-V(x-y)]\cdot\nabla f(x-y)\,dy.
\]

Fix a material center `X` and affine jet

\[
A=\nabla V(X).
\]

Transport the multiplier by the dual affine flow, so its instantaneous motion
cancels the affine commutator.  Then exactly

\[
\boxed{
(\partial_tQ+[V\cdot\nabla,Q])f(x)
=
\int K_N(y)
[V(x)-V(x-y)-Ay]\cdot\nabla f(x-y)\,dy.
}
\]

This is an operator identity.  It is not a semiclassical asymptotic expansion.
All common translation and affine strain/rotation have disappeared before any
estimate is made.

## 2. Only velocity curvature remains

If the velocity Hessian is bounded by `H` on the relevant line segments, the
integral form of Taylor's theorem gives

\[
\begin{aligned}
V(x)-V(x-y)-Ay
&=\int_0^1[\nabla V(x-sy)-A]y\,ds,
\end{aligned}
\]

hence

\[
\boxed{
|V(x)-V(x-y)-Ay|
\le
H\left(|x-X||y|+\frac12|y|^2\right).
}
\]

The vector linearized equation has one additional affine-amplitude remainder,

\[
Q[(\nabla V-A)f],
\]

and this is also controlled by the same Hessian:

\[
\|Q[(\nabla V-A)f]\|_2
\le \|K_N\|_1 H\,\||x-X|f\|_2.
\]

Suppose the dimensionless kernel moments satisfy

\[
\|K_N\|_1\le k_0,
\qquad
\int|y||K_N(y)|dy\le\frac{m_1}{N},
\qquad
\int|y|^2|K_N(y)|dy\le\frac{m_2}{N^2}.
\]

For a role/probe define

\[
M_x=N\frac{\||x-X|f\|_2}{\|f\|_2},
\quad
M_g=\frac{\|\nabla f\|_2}{N\|f\|_2},
\quad
M_{xg}=\frac{\||x-X|\nabla f\|_2}{\|f\|_2}.
\]

Minkowski gives the clean full scalar-cell/vector-amplitude estimate

\[
\boxed{
\frac{\|R_{Eg}f\|_2}{\|f\|_2}
\le
\frac{H}{N}
\left(k_0M_x+m_1M_{xg}+\frac32m_2M_g\right).
}
\]

Thus on a parabolic lifetime the natural dimensionless quantity is precisely
`int H/N dt`, the affine-curvature scaling.  There is no independent tower of
third, fourth, ... Taylor remainders in the moving-multiplier equation.

For the isotropic normalized coherent Gaussian with
`z=N(x-X)` standard Gaussian and `q=|k|/N`, the moments are exact:

\[
M_x=\sqrt3,
\qquad
M_g=\sqrt{q^2+3/4},
\qquad
M_{xg}=\sqrt{15/4+3q^2}.
\]

Smooth cell kernels have fixed finite `k_0,m_1,m_2`; Leray/helical variation is
still handled by the already summable symbol-freezing ledger.

## 3. Center curvature alone is not enough

A tempting conclusion would be to replace `H` above by the single center tensor
`H(X)`, hence by the existing `B(X)=L^{-1}H(X)[L,L]`.  This is false.

Take the smooth divergence-free shear

\[
\boxed{
V(x)=
\left(
0,
 a\left[\sin(rx_1)-\frac12\sin(2rx_1)\right],
0
\right),
\qquad r=N/8.
}
\]

Both Fourier modes lie in the strict transporter ball `B_(N/4)`.  At `X=0`,

\[
V(0)=0,
\qquad
\nabla V(0)=0,
\qquad
\nabla^2V(0)=0,
\]

but

\[
\boxed{
\partial_1^3V_2(0)=3ar^3\ne0.
}
\]

Indeed

\[
\sin z-\frac12\sin2z
=\frac12z^3-\frac18z^5+\cdots.
\]

So a coherent role centered at the origin sees genuine non-affine transport even
though the **point-sampled affine and quadratic jets vanish**.  Any theorem that
claims the full moving-role residual is bounded solely by `B(X)` is therefore
wrong under the current hypotheses.

This is a countermodel to an intermediate bridge, not to Navier--Stokes
regularity.  The field is smooth, divergence free and strict-lowpass, and may
occur as the resolved part of smooth divergence-free initial data.

## 4. The physically natural repair: coherent deformation variance

The countermodel says the affine jet should be associated with the **whole
coherent eddy**, not one spatial point.

In intrinsic Gaussian coordinates put

\[
F(z)=L^{-1}\nabla V(X+Lz)L.
\]

Define its coherent affine jet and deformation variance by

\[
\boxed{
\bar A_C=\mathbb E_\gamma F,
\qquad
\mathcal K_C^2
=\mathbb E_\gamma\|F-\bar A_C\|_F^2.
}
\]

This quantity has an immediate physical meaning:

> `K_C` measures how much the resolved deformation gradient varies across the
> actual coherent eddy after common affine deformation has been quotiented.

For a genuinely affine resolved flow, `F` is constant and `K_C=0` exactly.
Under a common invertible affine change of physical coordinates,
`L->ML` and `grad V->M(grad V)M^-1`, so `F` and `K_C` are unchanged.

## 5. Gaussian Poincare identifies the curvature controlling it

Differentiate in intrinsic coordinates:

\[
\nabla_zF
=
L^{-1}(\nabla^2V)(X+Lz)[L,L]
=:B(z).
\]

The standard Gaussian Poincare inequality applies componentwise:

\[
\boxed{
\mathcal K_C^2
\le
\mathbb E_\gamma\|B(z)\|_F^2.
}
\]

In Hermite coordinates this is transparent.  If

\[
F-\bar A_C=\sum_{|\alpha|\ge1}c_\alpha H_\alpha,
\]

then

\[
\operatorname{Var}_\gamma(F)=\sum|c_\alpha|^2,
\qquad
\mathbb E|\nabla_zF|^2
=
\sum|\alpha||c_\alpha|^2.
\]

Hence every nonconstant spatial deformation mode is seen by the coherent
curvature energy, including the center-flat cubic shear above.

This suggests that the correct no-escape observable for the outer moving role is
not just `B(X)` but a **coherent/cellwise curvature measure**.

## 6. Relation to the current H1/H3 architecture

Pointwise, every differentiated-incompressible `B(z)` still has the exact
irreducible split

\[
15=(7\oplus3)_{H3}\oplus5_{H1},
\]

with

\[
\|\operatorname{Sym}B(z)\|^2+\|C_H(z)\|^2
\ge\frac16\|B(z)\|^2.
\]

After Gaussian averaging the same inequality remains true for the **curvature
energy**.  What is not yet proved is that the existing time-dependent H1/H3
sideband no-escape theorems, which were formulated around an osculating packet
jet, consume this full coherent/cellwise curvature energy with no missing higher
Hermite channel.

That is now the precise next mathematical bridge.

There are two honest possibilities:

1. extend H1/H3 no-escape from center curvature to the coherent curvature
   measure above; or
2. prove that the higher coherent Hermite part is uniformly small/transfer-costly
   on the strict-lowpass near-extremal branch and can be absorbed into flat
   erosion or multiplicative loss.

## 7. Scope

This note proves the exact affine-subtracted Egorov identity and gives a precise
Navier--Stokes-compatible countermodel to center-Hessian-only closure.  It also
identifies an affine-invariant coherent deformation observable controlled by
coherent curvature through Gaussian Poincare.

It does **not** yet prove the final collision of that coherent curvature with the
existing master currencies, and it makes no global-regularity claim.
