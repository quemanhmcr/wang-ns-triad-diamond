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

## 2. Cap-or-entropy reset lemma

Let

\[
\mu=\sum_iw_i\delta_{x_i},
\qquad b=\int x\,d\mu(x).
\]

Fix an angular radius `alpha in (0,pi/2)` and a tolerated leakage
`epsilon_cap in (0,1)`, and put

\[
\eta_{cap}=\epsilon_{cap}(1-\cos\alpha).
\]

There are two cases.

If `|b|<=1-eta_cap`, the exact barycenter--collision inequality gives

\[
\boxed{
H_2(\mu)\ge
h_{cap}:=\log\frac{2}{2-\eta_{cap}}>0.
}
\]

This is a costly entropy block.

If `|b|>1-eta_cap`, set `n=b/|b|`. Since

\[
\int(1-n\cdot x)\,d\mu=1-|b|<\eta_{cap},
\]

Markov gives

\[
\boxed{
\mu\{d_{S^2}(x,n)\le\alpha\}
\ge1-\epsilon_{cap}.
}
\]

Thus every low-entropy synchronized flat episode has a dominant core of mass
at least `1-epsilon_cap` inside a fixed cap.  On this core the spherical
potential

\[
P=\sup_x[-\log(n\cdot x)]
\]

starts below the universal reset value

\[
\boxed{P_{max}=-\log\cos\alpha.}
\]

The leaked mass is assigned to the already-summable cross/remainder budget.
This reset lemma is what prevents a costly block from recharging the spherical
potential without bound.

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

## 6. Remaining analytic bridge

The finite-dimensional atomic model is now closed modulo the quantitative
single-edge stability input.  For the actual Navier--Stokes PDE one still has
to prove a Gaussian atomic extraction with uniform synthesis constants,
summable near-extremal errors, and a certified stability constant converting
Hodge residual energy into an actual transfer deficit.
