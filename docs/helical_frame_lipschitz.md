# Uniform helical-frame regularity on the signed-good extremal core

The helical line bundle over the full frequency sphere has nonzero Chern number,
so no global smooth phase convention exists.  A local packet theorem does not
need a global gauge.  The signed-good extremal core is uniformly separated from
triad collinearity, and the triad-normal gauge is uniformly regular there.

At the repository good threshold `eta_0=10^-4`, single-edge stability gives

\[
0\le u\le\frac1{200},
\qquad |v|\le\frac1{100}.
\]

With child length normalized to one,

\[
x=r_*e^{-v-u/2},
\qquad
y=r_*e^{-v+u/2}.
\]

If `theta` is the angle between the two parent vectors, the law of cosines
simplifies exactly to

\[
\boxed{
\cos\theta
=\frac{e^{2v}}{2r_*^2}-\cosh u.
}
\]

Monotonicity in `u,v,r_*`, together with the certified `r_*` bracket, yields by
Arb

\[
\boxed{
\frac14<\cos\theta<\frac25,
\qquad
\sin\theta>\frac9{10}.
}
\]

Let `a,b` be the unit parent directions and

\[
n=\frac{a\times b}{|a\times b|}.
\]

For any differentiable variation through the good core,

\[
\dot n
=\frac{(I-nn^T)(\dot a\times b+a\times\dot b)}{|a\times b|},
\]

hence

\[
\boxed{
\|\dot n\|
\le\frac{10}{9}(\|\dot a\|+\|\dot b\|).
}
\]

For the triad-normal helical gauge

\[
h_s(a;n)=\frac{n\times a+i s n}{\sqrt2},
\]

\[
\dot h_s
=\frac{\dot n\times a+n\times\dot a+i s\dot n}{\sqrt2}.
\]

Therefore

\[
\boxed{
\|\dot h_s(a;n)\|
\le\frac52(\|\dot a\|+\|\dot b\|).
}
\]

For the child role, using the same triad normal, the identical `5/2` constant
controls the sum of the two parent-direction rates and the child-direction
rate.

This does not close the continuum moving-packet theorem by itself.  It removes
a different possible obstruction: **helical topology does not create a local
chart singularity on a near-extremal packet block**.  Once carrier-direction
variation is controlled by frequency-cell freezing, affine transport and the
curvature-balanced spatial ledger, helical-frame variation follows linearly
with a scale-free constant.
