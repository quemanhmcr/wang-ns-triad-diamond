# The quadratic swirl kernel of affine Gaussian forcing

The third-Hermite forcing theorem has a nontrivial kernel.  That kernel is not a
failure of coercivity; incompressibility gives it a precise physical meaning.

Let

\[
V_a(z)=\frac12B_{abc}z_bz_c,
\qquad B_{abc}=B_{acb},
\]

and impose the differentiated incompressibility condition

\[
B_{aac}=0.
\]

If the fully symmetric tensor vanishes,

\[
\operatorname{Sym}_{abc}B=0,
\]

then there is a unique symmetric trace-free matrix `M` such that

\[
\boxed{
B_{abc}=\varepsilon_{abd}M_{dc}+\varepsilon_{acd}M_{db}.
}
\]

It is recovered by

\[
\boxed{
M_{dc}=\frac13\varepsilon_{abd}B_{abc}.
}
\]

Consequently

\[
\boxed{V(z)=z\times(Mz).}
\]

This field is divergence free and tangent to every normalized Gaussian sphere:

\[
V(z)\cdot z=0.
\]

Therefore

\[
V\cdot\nabla e^{-|z|^2/4}=0.
\]

For a carrier `q`, the only scalar action is

\[
q\cdot V(z),
\]

which is a quadratic form in `z` and hence belongs to the chirp/covariance
tangent space of the Gaussian manifold.

Thus the five-dimensional symmetric-tracefree `M` family is an exact
**ellipsoid-normalized quadratic swirl**.  It must not be charged as scalar
non-affine envelope forcing.  In the vector Navier--Stokes packet it can still
act through spatially varying strain/vorticity and helical polarization; that is
where this mode must be tested next.
