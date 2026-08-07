# Affine-covariant Gaussian forcing: only the third Hermite chaos is transverse

Let a physical Gaussian grain be written in its own affine coordinates

\[
x=X+Lz,\qquad \Sigma_x=LL^T,
\]

so that its unchirped envelope is proportional to `exp(-|z|^2/4)` and
`|psi|^2` makes `z` a standard three-dimensional Gaussian.  If

\[
H_{ijk}=\partial_j\partial_k U_i(X),
\]

define the **grain-normalized velocity curvature**

\[
\boxed{
B_{abc}=(L^{-1})_{ai}H_{ijk}L_{jb}L_{kc},
\qquad q=L^Tk.
}
\]

These variables are covariant under an arbitrary common affine change of
physical coordinates.  They are also independent of the orthogonal choice of
frame factor in `Sigma_x=LL^T` up to the corresponding tensor rotation.

For the quadratic Taylor remainder

\[
R_2(x)=\frac12 H[x-X,x-X],
\]

one has

\[
R_2\cdot\nabla\psi
=
\frac{i}{2}(q_aB_{abc})z_bz_c\psi
-
\frac14 B_{abc}z_az_bz_c\psi.
\]

The first term is a quadratic phase: it is tangent to the Gaussian manifold and
changes the wavefront chirp/covariance.  Only the full symmetrization
`T=Sym B` enters the cubic polynomial.  Wick's formula gives the exact
unprojected norm

\[
\frac{\|R_2\cdot\nabla\psi\|_2^2}{\|\psi\|_2^2}
=
\frac14\left(2\|q\cdot B\|_F^2+\operatorname{tr}(q\cdot B)^2\right)
+
\frac1{16}\left(6\|T\|_F^2+9\|\operatorname{tr}T\|^2\right).
\]

But

\[
T_{abc}z_az_bz_c
=
T:H_3(z)+3(\operatorname{tr}T)\cdot z.
\]

The linear part is also tangent to the Gaussian manifold: it changes the
center/carrier parameters.  After projecting onto the complement of the full
center--carrier--covariance--chirp tangent space, the first genuine shape
forcing is therefore exactly

\[
\boxed{
\frac{\|F_\perp\|_2^2}{\|\psi\|_2^2}
=
\frac38\|\operatorname{Sym}B\|_F^2,
}
\]

hence

\[
\boxed{
\frac{\|F_\perp\|_2}{\|\psi\|_2}
\le \frac{\sqrt6}{4}\|B\|_F.
}
\]

This is the appropriate Gaussian-beam interpretation of non-affine Navier--
Stokes advection.  A quadratic wavefront phase is not a defect; it is an
osculating packet parameter.  The transverse error starts in the third Hermite
chaos.

A large Euclidean condition number of `L` is likewise not itself a cost.  If
physical coordinates, velocity and the grain are transformed together, `B`
and `q` are unchanged.  The dynamic question is whether the physical velocity
Hessian is small in the grain's own metric.
