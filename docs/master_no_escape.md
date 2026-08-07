# Master finite-dimensional no-escape theorem

This theorem glues the previously proved finite-dimensional modules. It is an
abstract closure theorem for the grain-cascade model, not a Navier--Stokes
regularity proof.

## 1. Cost-or-episode architecture

For each block `j=0,...,L-1` let `R_j` be its normalized transfer efficiency.
After transfer-adapted extraction suppose

\[
R_j\le e^{-C_j}+\eta_j,
\qquad C_j\ge0,
\]

and define the exact cross penalty

\[
\xi_j=\log(1+\eta_j e^{C_j}).
\]

Then

\[
-\log R_j\ge C_j-\xi_j.
\]

The Bellman, ancestry, resistance, Hodge, and balanced-spherical modules
supply a fixed threshold `c0>0` so that every **costly** block satisfies

\[
C_j\ge c_0.
\]

The only low-cost blocks are synchronized, flat, hemispherical blocks on a
dominant directional core.

## 2. Barycentric entropy-or-reset lemma

Let

\[
\mu=\sum_iw_i\delta_{x_i},
\qquad b=\int x\,d\mu(x),
\]

and fix `beta in (0,1)`.  The exact barycenter--collision inequality from the
spherical module is

\[
H_2(\mu)\ge\log\frac{2}{1+|b|}.
\]

Hence there are two cases.

If `|b|<=beta`, the block pays the fixed atomic entropy cost

\[
\boxed{H_2(\mu)\ge h_\beta:=\log\frac{2}{1+\beta}>0.}
\]

If `|b|>beta`, define the barycentric potential

\[
P(\mu):=-\log|b(\mu)|.
\]

Then it has the uniform reset bound

\[
\boxed{0\le P(\mu)<P_{max}:=-\log\beta.}
\]

No mass is discarded and no new cross-error is created by this reset.

For exact equal-marginal flat propagation the spherical midpoint identity gives

\[
b_{j+1}=b_j/c_*,
\qquad c_*=\cos(\theta_*/2),
\]

and therefore

\[
P_{j+1}=P_j-\kappa_*,
\qquad \kappa_*=-\log c_*.
\]

The nonsymmetric companion identity

\[
b(\nu)=2c_*b(\lambda)-b(\mu)
\]

is used as the gate into this episode: if old/companion marginals are not
synchronized well enough for the barycentric erosion inequality, then either
the companion barycenter is small (hence it pays collision entropy) or it is
concentrated and is sent to the fresh/reuse ancestry classification.  Thus the
only low-cost transition is precisely the synchronized flat transition for
which the potential erosion below is valid.

A particularly natural threshold is

\[
\beta=c_*.
\]

Then

\[
P_{max}=\kappa_*,
\qquad
h_\beta=\log\frac{2}{1+c_*}\approx0.09516.
\]

In the exact model this means a zero-cost synchronized flat episode contains
at most one transition before another entropy/fresh/reuse event must occur.

## 3. Erosion inside one flat episode

For an exact optimal midpoint step with no fresh direction,

\[
P_{j+1}\le P_j-\kappa_*,
\qquad
\kappa_*=-\log\cos(\theta_*/2)>0.
\]

For a near-extremal block allow a nonnegative perturbation `zeta_j` and a
fixed `kappa0>0`:

\[
\boxed{P_{j+1}\le P_j-\kappa_0+\zeta_j.}
\]

If a flat episode contains `m` low-cost blocks, summing while `P>=0` gives

\[
m\kappa_0\le P_{max}+\sum_{episode}\zeta_j.
\]

## 4. Master no-escape theorem with potential resets

Let `N_K` be the number of costly blocks and `N_F` the number of low-cost flat
blocks.  Costly blocks separate the flat cascade into at most `N_K+1`
episodes. Hence

\[
N_F\kappa_0
\le
(N_K+1)P_{max}+Z,
\qquad Z:=\sum_j\zeta_j.
\]

Since `L=N_K+N_F`,

\[
\boxed{
N_K\ge
\frac{\kappa_0L-P_{max}-Z}{\kappa_0+P_{max}}.
}
\]

Therefore

\[
\sum_jC_j
\ge
\frac{c_0}{\kappa_0+P_{max}}
(\kappa_0L-P_{max}-Z).
\]

Combining with the logarithmic cross penalties gives

\[
\boxed{
-\log\prod_{j<L}R_j
\ge
c_{eff}L
-\frac{c_0(P_{max}+Z)}{\kappa_0+P_{max}}
-\Xi,
}
\]

where

\[
\boxed{
c_{eff}=\frac{c_0\kappa_0}{\kappa_0+P_{max}}>0,
\qquad
\Xi=\sum_j\xi_j.
}
\]

Thus, if the perturbation sum `Z` and cross-error penalty `Xi` remain bounded,

\[
\boxed{
\prod_{j<L}R_j\lesssim e^{-c_{eff}L}.
}
\]

This is the corrected closure theorem.  The potential may reset after every
costly block, but every reset creates a new episode and the number of episodes
is itself controlled by the number of costly blocks.

## 5. Where each costly block comes from

The previously proved finite-dimensional modules supply the alternatives:

1. transfer-weighted component splitting -> Bellman collision entropy;
2. hidden atomic entropy -> component entropy or same-ancestry pair cycles;
3. high resistance -> Poisson Bellman cut or conductance collision entropy;
4. low resistance plus curved reuse -> Hodge cost;
5. low resistance plus low Hodge -> integer gauge synchronization;
6. synchronized flat but barycentrically balanced -> atomic collision entropy;
7. synchronized flat and concentrated -> a dominant cap; nested extraction
   either splits it (Bellman cost) or produces a trackable fresh/reused grain.

A fresh grain gives a fresh ancestry split; a reused grain returns to 2--5.
The only branch not immediately paying a cost is therefore the cap-confined,
synchronized-flat, no-fresh transition treated by the episode erosion above.

The single-edge certificate removes the previous placeholder in branch 4.  If
the good-edge Hodge threshold is `E_H>=h_H`, then

\[
-\log R_{block}\ge\frac12E_H\ge\frac{h_H}{2}.
\]

If instead a set of outside-neighborhood edges has transfer weight at least
`beta_bad`, it pays at least `beta_bad/100`.  Thus, after fixing the positive
classification thresholds used by the other branches, one may take

\[
\boxed{
c_0=\min\{c_{Bell},c_{anc},c_{res},h_H/2,\beta_{bad}/100,
h_\beta,c_{fresh},\ldots\}>0.
}
\]

This formula does not claim that the threshold choices have been optimized; it
records that the Hodge/single-edge entry is now a theorem-level positive number
rather than a numerical candidate.

## 6. Remaining analytic bridge

The finite-dimensional single-edge stability input is now certified.  Moreover,
the sharp-cutoff single-triad Mellin identity explains the `log(q/p)` progress
exactly and retains at least `9/10` of the upper forward segment on the certified
local box.  For the actual Navier--Stokes PDE one still has to globalize this
identity: prove a transfer-adapted Gaussian atomic extraction with uniform
synthesis constants, control cancellation among many triads, produce the
positive normalized `w_e` used above, compare sharp and packet-compatible
localized fluxes, and keep the errors summable.  Pressure, time synchronization,
and the critical local-energy charge remain part of that PDE bridge.
