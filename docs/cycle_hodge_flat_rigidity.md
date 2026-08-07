# Cycle Hodge energy and flat-cell rigidity

This is a finite-dimensional structural module, not a Navier--Stokes regularity proof.

## Hodge projection replaces cycle packing

Each active triad `p+q -> c` contributes two directed arcs `p->c` and `q->c`.
Let `B` be the oriented incidence matrix and `ell_v=log |k_v|`. Every arc near
the symmetric triad extremizer prefers increment

\[
\ell_{\rm child}-\ell_{\rm parent}=\gamma_*,\qquad
\gamma_* = \log(1/r_*).
\]

For positive arc weights `W`,

\[
\mathcal E_G(\ell)=\|W^{1/2}(B\ell-\gamma_*\mathbf 1)\|_2^2.
\]

The exact minimum is

\[
\inf_\ell \mathcal E_G(\ell)
=\operatorname{dist}_W^2(\gamma_*\mathbf1,\operatorname{ran}B)
=\sup_{z\in\ker B^T\setminus\{0\}}
\frac{\gamma_*^2(z\cdot\mathbf1)^2}{z^TW^{-1}z}.
\]

This uses the complete cycle space and loses no constant through overlapping
cycle bases. The connection is flat exactly when there exists a level function
`h` satisfying `h_child=h_parent+1` on every active arc.

## Single-triad residual identity

Put

\[
u=\ell_p-\ell_q,
\qquad
v=\ell_c-\frac{\ell_p+\ell_q}{2}-\gamma_*.
\]

The two arc residuals satisfy

\[
\boxed{r_p^2+r_q^2=2v^2+u^2/2.}
\]

Thus local multiplier stability in the imbalance and mean-shift variables
controls the global Hodge energy.

## Curved reuse motif

For

\[
a+b\to m,\qquad m+c\to d,\qquad b+c\to n,
\]

the cycle rank is one, but no level function exists. With unit arc weights,

\[
\inf\mathcal E_G=\gamma_*^2/5.
\]

## Flat butterfly countermodel

For

\[
a+b\to m,\qquad a+c\to n,\qquad m+n\to d,
\]

the levels `a,b,c=0`, `m,n=1`, `d=2` make every preferred increment exact.
It has cycle rank one but zero scale Hodge energy.

The flat graph has an exact extremal vector realization. Let unit vectors obey

\[
a\cdot b=a\cdot c=\chi=\cos\theta,
\]

put `R=sqrt(2+2 chi)`, and define

\[
m=(a+b)/R,\qquad n=(a+c)/R.
\]

If the third triad is also extremal, `m dot n=chi`, then

\[
b\cdot c=2\chi^2-1,
\qquad b+c=2\chi a,
\qquad \boxed{d=(m+n)/R=a}.
\]

The internal level-one directions lie at angle `theta/2` from `d`; a new
extremal companion for `d` must lie at angle `theta`. Hence an isolated flat
butterfly cannot continue using only its internal packets. It requires fresh
boundary data at the new scale.

The corrected dichotomy is therefore

\[
\boxed{
\text{curved cycle}\Rightarrow\text{Hodge cost},
\qquad
\text{flat extremal cell}\Rightarrow\text{fresh boundary cost}.
}
\]

A general flat graded network may contain larger equal-length diamonds. The
remaining task is to decompose such networks into flat cells and prove that a
positive density of their boundary companions is fresh, unless an even more
rigid geometric class occurs.
