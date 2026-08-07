# Multicommodity Hodge routing and gauge synchronization

This module is finite-dimensional. It is not a Navier--Stokes regularity theorem.

Let a directed graph have incidence matrix `D`, positive conductances `W`, and a target edge cochain `a`. Define

\[
E_H=\inf_\phi \|W^{1/2}(a-D^T\phi)\|_2^2.
\]

For every cycle flow `z` with `Dz=0`,

\[
(z\cdot a)^2\le E_H\, z^TW^{-1}z.
\]

Hence for any family of cycle commodities `(z_r,mu_r)`,

\[
\boxed{
E_H\ge
\frac{\sum_r\mu_r(z_r\cdot a)^2}
{\sum_r\mu_r z_r^TW^{-1}z_r}.
}
\]

The denominator is the total electrical congestion of the chosen cycle routing. This is basis-free: overlapping cycles are handled by one Rayleigh quotient rather than by summing individual cycle gaps.

## Old/new ancestry routing

Let `G_old` and `G_new` be two coloured networks on the same reused terminals. Route one unit electrically from terminal `i` to `j` in each network. The difference

\[
z_{ij}=f^{new}_{ij}-f^{old}_{ij}
\]

is divergence free in the coloured union, and

\[
\|z_{ij}\|_{W^{-1}}^2
=R^{new}_{eff}(i,j)+R^{old}_{eff}(i,j)
\]

when the two edge colours are disjoint.

For terminal probabilities `p_i`, define

\[
\mathfrak R=
\sum_{i,j}p_ip_j
\left(R^{new}_{eff}(i,j)+R^{old}_{eff}(i,j)\right).
\]

If a positive fraction of pair mass has scale circulation at least `delta`, then

\[
E_H\ge \frac{\delta^2\,\text{(that pair fraction)}}{\mathfrak R}.
\]

## Gauge synchronization

Suppose each coloured network is exactly graded:

\[
a_{old}=D_{old}^Th_{old},\qquad
a_{new}=D_{new}^Th_{new}.
\]

On reused terminals set

\[
d_i=h_{new}(i)-h_{old}(i).
\]

For every electrical pair cycle,

\[
z_{ij}\cdot a=d_j-d_i.
\]

Therefore the multicommodity numerator is

\[
\sum_{i,j}p_ip_j(d_i-d_j)^2=2\operatorname{Var}_p(d),
\]

and

\[
\boxed{
\operatorname{Var}_p(d)
\le \frac12 E_H\mathfrak R.
}
\]

If each `d_i` belongs to `gamma_* Z`, then

\[
\mathbb P(d_I\ne d_J)
\le \frac{2\operatorname{Var}_p(d)}{\gamma_*^2}
\le \frac{E_H\mathfrak R}{\gamma_*^2}.
\]

Thus low Hodge energy plus low pair-resistance budget forces a single generation gauge on most reused pair mass. This is the synchronization needed to concatenate flat spherical cells across blocks.

## The new obstruction: resistance

The only quantity in the theorem not already controlled by the previous ledgers is the pair resistance budget `mathfrak R`. This is not a technical nuisance: a tree-like network can make `mathfrak R` large even when every individual edge is locally good.

For a weighted tree with edge resistances `r_e`, deleting edge `e` splits terminal probability into masses `s_e` and `1-s_e`. Exactly

\[
\boxed{
\sum_{i,j}p_ip_jR_T(i,j)
=2\sum_e r_e s_e(1-s_e).
}
\]

So large routing resistance is itself a weighted cut ledger: it must come from edges that simultaneously have large resistance and separate nontrivial terminal mass. The next module should convert this cut-resistance quantity into either a Bellman split or a summable cross-edge moat.
