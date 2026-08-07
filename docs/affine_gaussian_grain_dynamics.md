# Affine Gaussian grain dynamics and strain rigidity

This module is the first spacetime step after the frozen SGS/profile bridge.  It
uses the exact kinematics of an incompressible affine flow rather than imposing a
packet-motion ansatz.  The statement is still a local affine/Kelvin-Gaussian
model theorem, not a global Navier--Stokes regularity result.

## 1. Why absolute strain is the wrong defect

A sharp-Young Gaussian extremizer has affine symmetries, while a Fourier triad is
transported by the velocity gradient.  Therefore a common rigid rotation or a
common scalar dilation of the triad plane must not be charged as a defect.  The
physical quantity that changes the Euclidean shape of the helical triad is the
**non-conformal planar strain**.

For an affine resolved velocity

\[
U(x,t)=U_0(t)+A(t)(x-X(t)),\qquad \operatorname{tr}A=0,
\]

a Kelvin carrier satisfies

\[
\dot k=-A^T k.
\]

Only the symmetric part changes `|k|`; the antisymmetric part is rigid rotation.
If `S=(A+A^T)/2`, then

\[
\frac{d}{dt}\log|k|=-\widehat k^T S\widehat k.
\]

This is the classical Kelvin/shearing-wave mechanism.  Exact Navier--Stokes
shearing-wave solutions in linear backgrounds are a useful physical sanity
check for using this kinematics rather than an arbitrary packet ODE.

## 2. Extremal triad measures planar conformal strain

At the symmetric single-edge optimizer, let `phi=theta_*/2`; then

\[
\cos\phi=\frac1{2r_*},\qquad
n_a=(\cos\phi,\sin\phi),\quad
n_b=(\cos\phi,-\sin\phi),\quad
n_c=(1,0).
\]

Use the signed Hodge coordinates

\[
u=\ell_b-\ell_a,
\qquad
v=\ell_c-\frac{\ell_a+\ell_b}{2}-\gamma_*.
\]

Let the trace-free restriction of the symmetric carrier driver to the triad
plane be

\[
D=\begin{pmatrix}\delta&\beta\\ \beta&-\delta\end{pmatrix}.
\]

A direct differentiation gives

\[
\dot u=4\beta\cos\phi\sin\phi,
\qquad
\dot v=-2\delta\sin^2\phi.
\]

Hence

\[
\frac12\dot u^2+2\dot v^2
=8\cos^2\phi\sin^2\phi\,\beta^2
+8\sin^4\phi\,\delta^2.
\]

Since `||D||_F^2=2(delta^2+beta^2)` and `sin phi < cos phi`, one obtains the
exact coercivity coefficient

\[
\boxed{
\frac12\dot u^2+2\dot v^2
\ge
4\sin^4\phi\,\|D\|_F^2.
}
\]

The Arb certificate proves

\[
\boxed{
4\sin^4\phi>\frac{43}{100}.
}
\]

Thus a maximally efficient triad can ride a common rotation and a planar scalar
strain for free, but any trace-free in-plane stretching/shear immediately moves
it away from the optimal multiplier geometry.  The kernel is exactly the
planar-conformal direction at first order.

## 3. A finite packet-lifetime theorem for frozen principal strain

Take a frozen symmetric trace-free planar strain with eigenvalues `+/- d` and
start from the exact extremal triad.  For any unit carrier direction under
`k(t)=exp(-Dt)k(0)`, if `alpha=d log|k|/dt`, then

\[
|\dot\alpha|\le2d^2.
\]

In the Hodge vector

\[
z=(u/\sqrt2,\sqrt2 v),
\]

the preceding coercivity gives

\[
|\dot z(0)|^2>0.86d^2,
\]

while

\[
|\dot z(t)-\dot z(0)|\le\sqrt{40}\,d^2t.
\]

Therefore

\[
|z(t)|
\ge dt\left(\sqrt{0.86}-\sqrt{10}\,dt\right).
\]

If

\[
\boxed{dT\le\frac1{25}},
\]

then Arb certifies that the squared bracket is `>3/5`.  Moreover
`|u|,|v|<=2dT<=2/25`, so the whole episode stays inside the already-certified
single-edge local box.  Thus

\[
\mathcal H(t):=\frac{u(t)^2}{2}+2v(t)^2
\ge\frac35d^2t^2
\]

and the previous edge theorem `Def>=H/2` yields

\[
\boxed{
\frac1T\int_0^T\operatorname{Def}(t)\,dt
\ge\frac1{10}(dT)^2.
}
\]

For a Navier--Stokes packet lifetime `T=c N^{-2}`, this is a dimensionless
quadratic cost in the non-conformal strain number `d/N^2` whenever the local
affine gradient is coherent over the packet lifetime.  If the gradient is not
coherent, its time/space variation becomes the next explicit error branch rather
than being hidden inside the packet ansatz.

## 4. Viscosity preserves the common Gaussian role structure

For the scalar Gaussian envelope of

\[
f_t+(Ax)\cdot\nabla f=\nu\Delta f,
\]

write the Fourier Gaussian with precision matrix `P` and spectral peak `kappa`.
The exact coefficient equations are

\[
\boxed{
\dot P=AP+PA^T+2\nu I,
\qquad
\dot\kappa=-A^T\kappa-2\nu P^{-1}\kappa.
}
\]

The dual center

\[
b=P\kappa
\]

has the striking exact law

\[
\boxed{\dot b=A b.}
\]

Therefore, if three roles have the same Gaussian precision and satisfy
`kappa_a+kappa_b=kappa_c`, resonance is preserved exactly by the common affine
Gaussian evolution.  Equal widths also remain equal because all three copies
solve the same matrix ODE.

Incompressibility gives

\[
\boxed{
\frac{d}{dt}\log\det P
=2\nu\operatorname{tr}(P^{-1})\ge0.
}
\]

Strain does not change this determinant directly; viscosity narrows the Fourier
packet (equivalently broadens its physical envelope).

If `P=pI` at the beginning of the packet lifetime, the viscous contribution to
the carrier-length driver is the scalar matrix `2 nu p^{-1} I`.  Hence its
trace-free planar part is zero: **viscosity is first-order neutral for triad
shape at an isotropic Gaussian instant**.  The first-order shape cost above is
therefore genuinely a strain effect, not an artifact of dissipation.

## 5. Vector Kelvin polarization and physical work

For a transverse Kelvin mode of the linearized incompressible Navier--Stokes
equation around `U=Ax`, pressure enforces transversality through

\[
\dot a=-Aa+2k\frac{k\cdot Aa}{|k|^2}-\nu|k|^2a.
\]

Together with `dot k=-A^T k`, this preserves `k dot a=0`.  Since the pressure
term is parallel to `k`, it does no direct work on the transverse amplitude:

\[
\boxed{
\frac d{dt}|a|^2
=-2a\cdot Sa-2\nu|k|^2|a|^2.
}
\]

Thus the same strain tensor that deforms triad geometry also performs physical
Reynolds/Kelvin work on the packet polarization.  This is the correct dynamical
object to couple to the SGS transfer ledger.

## 6. Affine-deforming moving windows cancel advection exactly

The previous localized SGS ledger used a moving packet window.  The natural
window is not merely translated; it should deform with the resolved velocity
gradient.  Let

\[
\dot X=U(X,t),\qquad
\dot F=(\nabla U)(X,t)F,
\]

and define

\[
\chi(x,t)=\chi_0\!\left(F^{-1}(x-X)/R\right).
\]

For the affine Taylor velocity

\[
U_{aff}(x)=U(X)+(\nabla U)(X)(x-X),
\]

one has the exact material cancellation

\[
\boxed{
(\partial_t+U_{aff}\cdot\nabla)\chi=0.
}
\]

If `||nabla^2 U||<=H` on the window, Taylor's theorem gives

\[
\boxed{
|(\partial_t+U\cdot\nabla)\chi|
\le
\frac H2 R\,\|F\|^2\|F^{-1}\|\,\|\nabla\chi_0\|_\infty.
}
\]

So the dominant translational **and linear-strain** leakage is removed exactly.
The remaining advective leakage is a curvature quantity.  This makes the next
PDE target precise: control the dimensionless Hessian/strain-variation error on
an `N^{-2}` lifetime, or show that its failure produces fresh/Bellman critical
mass.

## 7. What remains

The current theorem does not assert that a genuine Navier--Stokes grain sees a
frozen affine gradient throughout its lifetime.  It identifies the exact neutral
symmetries and the quantity that must be controlled:

- persistent trace-free strain in the triad plane pays an explicit average
  multiplier cost;
- common affine transport can be absorbed into the moving Gaussian/window frame;
- viscosity preserves the common Gaussian role structure and is first-order
  shape-neutral at isotropic width;
- spatial curvature and temporal rotation/variation of the strain frame are the
  true new perturbation variables.

The next closure should therefore be a **strain-coherence or curvature-cost
dichotomy** on the `N^{-2}` lifetime, not an attempt to penalize all strain.
