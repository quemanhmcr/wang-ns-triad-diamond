# Reservoir pool erosion: relinking cannot be postponed forever

The single-reservoir spectral half-life still leaves a possible loophole: perhaps an efficient cascade owns a very large **pool** of old low-frequency atoms and simply chooses a different old atom at every generation.  This note shows that the entire old pool erodes geometrically once its packet family has the standard Bessel/bounded-overlap energy budget.

Let generation zero have block scale `N_0` and old reservoir atoms with

\[
M_{a,0}\le\alpha N_0.
\]

Assume their transfer-adapted packet frame satisfies at every service time

\[
\boxed{
\sum_aE_{a,q}\le P\|u(t_q)\|_2^2
\le P E_{global}.
}
\]

The constant `P` is kept explicit; orthogonality is not assumed.

On a signed-good low-strain lineage,

\[
M_{a,q}\le(21/20)^qM_{a,0},
\qquad
N_q\ge(8/5)^qN_0.
\]

A low band contributes increment-square service proportional to

\[
{M_{a,q}^3\over N_q^2}E_{a,q}.
\]

Summing over the whole old pool,

\[
\boxed{
\mathsf C_{old}(q)
\le
\alpha^3N_0P E_{global}
\left({231525\over512000}\right)^q
<2^{-q}\alpha^3N_0P E_{global}.
}
\]

Therefore

\[
\boxed{
\sum_{q\ge0}\mathsf C_{old}(q)
<2\alpha^3N_0P E_{global}.
}
\]

This bound is adversarial: the distribution of energy among old atoms may change completely from one generation to the next.  Only the global packet-frame energy budget and material frequency-growth law are used.

If every efficient block requires low/base service at least `sigma_*>0`, then after finitely many generations the old pool cannot provide even half of `sigma_*`.  A quantitative amount of service must come from **newly relinked spectral capacity**.

## Relinking enters the existing fresh-or-cycle graph

In the finite atomic quadratic model, represent each non-negligible relinking event by an active triad after the existing tiny-edge pruning.  For one connected 3-uniform incidence component with `m` triads and `n` packet vertices,

\[
\beta=2m-n+1,
\qquad
\boxed{(n-1)+\beta=2m.}
\]

Hence

\[
\boxed{
n-1\ge m
\quad\text{or}\quad
\beta\ge m.
}
\]

So a long sequence of relinking events cannot form a third neutral regime: it is either fresh-rich or cycle-rich.  Combining this Euler theorem with old-pool erosion is the finite-atomic version of reservoir relinking synchronization.

## Remaining continuum point

The theorem assumes the iterative packet realization supplies a Bessel/bounded-overlap coefficient budget `P` and that newly generated service capacity is represented by the active quadratic interaction graph.  Establishing those properties with summable selected-interface `Xi` for the genuine PDE packet extraction remains the continuum bridge.
