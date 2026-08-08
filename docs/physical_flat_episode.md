# Physical flat episode: replacing the abstract master perturbation

The master theorem permits a near-flat step

\[
P_{j+1}\le P_j-\kappa_0+\zeta_j,
\]

but until now `zeta_j` was an abstract perturbation.  Signed-good triad geometry makes it explicit.

Let the parent magnitudes be

\[
x=r_*e^{-v-u/2},\qquad y=r_*e^{-v+u/2},
\]

with child magnitude normalized to one.  If `theta_e` is the angle between the parent **directions** and `c_e=cos(theta_e/2)`, the law of cosines gives exactly

\[
\boxed{c_e^2=c_*^2e^{2v}-\sinh^2(u/2),\qquad c_*=(2r_*)^{-1}.}
\]

On the certified local box `|u|,|v|<=2/25`, Arb gives

\[
\boxed{|c_e-c_*|\le\frac65|v|+\frac15u^2.}
\]

There is a second exact effect which is absent from the ideal equal-parent spherical model.  Let `p,q` be the unit parent directions, `m=(p+q)/|p+q|`, and `n_c` the direction of `xp+yq`.  In the midpoint/tangent basis,

\[
xp+yq=c_e(x+y)m+\sin(\theta_e/2)(x-y)d.
\]

Hence

\[
\tan\angle(n_c,m)
=\tan(\theta_e/2)|\tanh(u/2)|.
\]

The local half-angle is above `3/4`, so

\[
\boxed{\|n_c-m\|\le |u|/2.}
\]

## Transfer-weighted barycenter perturbation

Let a normalized positive-transfer coupling have old-parent barycenter `b_1`, companion barycenter `b_2`, and actual child barycenter `b_c`.  Put

\[
H=\mathbb E(2v^2+u^2/2)=E_H^{phys},
\qquad
\Delta_b=|b_2-b_1|.
\]

Combining child/midpoint error, angle error, `c_*>4/5`, and Cauchy gives

\[
\boxed{
\left|b_c-\frac{b_1}{c_*}\right|
\le e,
\qquad
e\le2\sqrt H+\frac12H+\frac58\Delta_b.
}
\]

This is a physical near-version of the exact flat identity `b_child=b_parent/c_*`.

On the concentrated branch `|b_1|>=c_*`, with `P=-log|b|`,

\[
P_c
\le
P_1-\kappa_*
-\log\left(1-\frac{c_*e}{|b_1|}\right)
\le
P_1-\kappa_* -\log(1-e).
\]

If `e<=1/2`, `-log(1-e)<=2e`, so the master perturbation can be chosen as

\[
\boxed{
\zeta
\le4\sqrt{E_H^{phys}}+E_H^{phys}+\frac54\Delta_b.
}
\]

## Uniform physical erosion on a 1%-flat block

The service-or-flat theorem at `tau=1/100` gives

\[
\sqrt{E_H^{phys}}\le\tau/3.
\]

If the two parent marginals are synchronized to

\[
\Delta_b\le\tau,
\]

then

\[
\boxed{
\zeta\le\frac{31}{12}\tau+\frac{\tau^2}{9}.
}
\]

Arb certifies

\[
\boxed{
\kappa_0
:=\kappa_*-\zeta
>\frac{17}{100}.
}
\]

Thus a physical `1%`-Kelvin-flat synchronized block is already a quantitative flat step in the master barycentric episode, with a fixed erosion rate.  The remaining interface task is sharply reduced: either parent marginal mismatch exceeds the chosen synchronization threshold and must be charged to selected-interface/fresh/entropy bookkeeping, or the flat step consumes at least `0.17` units of barycentric potential.
