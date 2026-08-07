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


The single-edge certificate now makes this quantitative.  On
`|u|<=2/25`, `|v|<=2/25`,

\[
\operatorname{Def}_e\ge\frac12(r_p^2+r_q^2).
\]

Hence, if a triad of transfer weight `w_e` gives weight `w_e` to each of its
two parent-to-child Hodge arcs,

\[
\sum_e w_e\operatorname{Def}_e
\ge\frac12\mathcal E_H.
\]

For the normalized block ratio `R_block=sum_e w_e J_e/J_*`,

\[
-\log R_{\rm block}\ge1-R_{\rm block}
\ge\frac12\mathcal E_H.
\]

Thus a Hodge threshold `E_H>=h_H` supplies the explicit theorem-level block
cost `c_{0,H}=h_H/2`.  Edges outside the certified local rectangle have
`Def_e>=1/100` and are charged before the Hodge network is formed.


For the three-triad nonflat reuse motif below, unit arc weights give exactly
`E_H=gamma_*^2/5`.  If the three triads have equal normalized transfer weights
`w_e=1/3` and **all three lie in the certified local rectangle**, each of their
two Hodge arcs inherits weight `1/3`, so linear scaling gives

\[
\mathcal E_H\ge\frac{\gamma_*^2}{15},
\qquad
-\log R_{\rm block}\ge\frac{\gamma_*^2}{30}
\approx0.00809556352.
\]

If at least one of the three triads is outside the local rectangle, its weight
`1/3` and the global edge gap `Def_e>=1/100` instead give

\[
-\log R_{\rm block}\ge\frac1{300}.
\]

Therefore the equal-transfer nonflat motif has the unconditional certified cost

\[
\boxed{
c_{0,\rm motif}
\ge\min\!\left\{\frac{\gamma_*^2}{30},\frac1{300}\right\}
=\frac1{300}>0.
}
\]

This is a concrete theorem-level positive cost for that normalized motif; it is
not asserted to be the universal master `c_0` across all costly branches.

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

## Planar flat erosion theorem

The flat butterfly suggests a stronger ledger when a component is contained in
a lifted planar angular chart. Write directions by real angles. An exact
extremal triad takes two parent angles `x<y` with

\[
y-x=\theta_*
\]

and produces the midpoint angle

\[
z=(x+y)/2.
\]

For a finite parent set `S`, every eligible midpoint lies in

\[
[\min S+\theta_*/2,\;\max S-\theta_*/2].
\]

Hence

\[
\boxed{\operatorname{diam} S_{j+1}\le
\operatorname{diam} S_j-\theta_*}
\]

whenever no fresh directions are inserted. A flat exact cascade in a fixed
angular chart therefore dies after at most `diam(S_0)/theta_*` generations.

More generally, let `F_j` be fresh boundary directions and put

\[
E_j=\operatorname{diam}(S_j\cup F_j)-\operatorname{diam}S_j.
\]

Then

\[
\operatorname{diam}S_{j+1}
\le \operatorname{diam}S_j+E_j-\theta_*.
\]

After `L` generations,

\[
\boxed{
\sum_{j=0}^{L-1}E_j
\ge L\theta_*+\operatorname{diam}S_L-\operatorname{diam}S_0.
}
\]

Thus a planar flat cascade can avoid Hodge curvature only by importing angular
boundary span at a linear rate.

For near-extremal pairs with angular separation at least
`theta_*-delta_j` and child-angle error at most `epsilon_j`, the same proof gives

\[
\operatorname{diam}S_{j+1}
\le \operatorname{diam}S_j+E_j
-(\theta_*-\delta_j)+2\epsilon_j.
\]

This is a second additive cocycle: curved networks pay Hodge energy, while
planar flat networks pay fresh angular span. Converting fresh span into
critical-mass entropy requires the packet inverse theorem's comparable-mass
conclusion.
