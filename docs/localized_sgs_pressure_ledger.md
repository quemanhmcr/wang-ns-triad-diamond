# Localized SGS / pressure-work depletion ledger

This module identifies what pressure can and cannot do once the smooth SGS
scale-transfer observable has been localized to a packet tube. It is an exact
resolved-energy calculation; the remaining analytic issue is to control the
resulting boundary/annular charges in a genuine Navier--Stokes extraction.

## 1. Resolved local energy equation

For a smooth convolution filter at scale `ell`, write

\[
U=\bar u^\ell,\qquad P=\bar p^\ell,
\qquad R=\overline{u\otimes u}-U\otimes U,
\]

and

\[
\Pi=-R:\nabla U,
\qquad e=\frac12|U|^2.
\]

Filtering Navier--Stokes and using incompressibility gives the exact resolved
energy balance

\[
\partial_t e+\nabla\cdot J=-\Pi-\nu|\nabla U|^2,
\]

where

\[
J=(e+P)U+RU-\nu\nabla e.
\]

Let `chi(x,t)>=0` be a smooth moving packet window. Multiplying by `chi` and
integrating over space-time gives

\[
\boxed{
\int\!\!\int \chi\Pi
+\nu\int\!\!\int\chi|\nabla U|^2
=E_\chi(t_0)-E_\chi(t_1)+L_\chi,
}
\]

with

\[
E_\chi(t)=\int\chi e,
\]

\[
L_\chi=\int\!\!\int e\,\partial_t\chi
+\int\!\!\int\nabla\chi\cdot
\left[(e+P)U+RU-\nu\nabla e\right].
\]

This formula is exact but pressure still appears in the spatial leakage.

## 2. Combined work removes pressure exactly

Define the signed combined pressure--flux work density

\[
\boxed{
G=\Pi+\nabla\cdot(PU).
}
\]

Since

\[
\int\chi\,\nabla\cdot(PU)=-\int PU\cdot\nabla\chi,
\]

the pressure term cancels from the localized leakage identity:

\[
\boxed{
\int\!\!\int\chi G
+\nu\int\!\!\int\chi|\nabla U|^2
=E_\chi(t_0)-E_\chi(t_1)+\widetilde L_\chi,
}
\]

where

\[
\widetilde L_\chi
=\int\!\!\int e\,\partial_t\chi
+\int\!\!\int\nabla\chi\cdot
\left[eU+RU-\nu\nabla e\right].
\]

Thus pressure is not an interscale source. It is a spatial transport channel,
and using `G` places it on the transport side of the ledger exactly rather than
estimating it prematurely.

## 3. Weighted finite-chain depletion

For consecutive packet-time intervals let

\[
W_j+D_j=E_j-E_{j+1}+L_j
\]

be the integrated combined-work identity. For arbitrary nonnegative weights
`a_j`, summation by parts gives

\[
\boxed{
\sum_{j=0}^{N-1}a_j(W_j+D_j-L_j)
=a_0E_0-a_{N-1}E_N
+\sum_{j=1}^{N-1}(a_j-a_{j-1})E_j.
}
\]

If the weights are nonincreasing, all energies are nonnegative, and
`W=W_+-W_-`, then

\[
\boxed{
\sum_ja_jW_{j,+}+\sum_ja_jD_j
\le a_0E_0+\sum_ja_j(L_j)_++\sum_ja_jW_{j,-}.
}
\]

So persistent forward combined work is paid by initial localized energy,
positive window leakage, or negative combined work/backscatter. Pressure cannot
reset this budget because it has already been included in `G`.

This identity is especially compatible with the transfer-weighted Bellman
architecture: decreasing lineage weights do not create a packet-count factor.

## 4. Pressure cancellation of raw SGS flux creates a critical annular charge

The Hodge/midgap module naturally first sees the raw SGS transfer. At a fixed
time/window, write

\[
S=\int\chi\Pi,
\qquad
W=\int\chi G
=S-\int PU\cdot\nabla\chi.
\]

Assume `S>=0`. Then there is an exact dichotomy:

\[
\boxed{
W\ge\frac S2
\quad\text{or}\quad
\left|\int PU\cdot\nabla\chi\right|\ge\frac S2.
}
\]

In the second branch, if `A=supp grad chi`, Holder and Young with exponents
`3/2` and `3` give

\[
\frac S2
\le \|\nabla\chi\|_\infty
\|P\|_{L^{3/2}(A)}\|U\|_{L^3(A)}
\]

and

\[
ab\le\frac23a^{3/2}+\frac13b^3.
\]

Therefore

\[
\boxed{
\int_A\left(|U|^3+|P|^{3/2}\right)
\ge
\frac{S}{2\|\nabla\chi\|_\infty}.
}
\]

This is precisely a scale-critical CKN-type boundary charge. Hence pressure can
cancel a forward SGS event only by pushing a definite critical quantity into
the spatial moat. That is a physical alternative which can be handed to the
fresh/leakage/profile ledger instead of being treated as an uncontrolled sign
error.

## 5. Multipole collision upgrade

The annular charge alternative can be sharpened before any spacetime counting.
Because `div U=0`,

\[
\int U\cdot\nabla\chi=0,
\]

so a constant pressure performs no window-boundary work.  For the Newtonian
pressure kernel this permits the far-field subtraction
`K_ij(x-y)-K_ij(x_0-y)`.  The explicit kernel satisfies

\[
\sum_{ij}|\nabla K_{ij}(z)|\,|u_i u_j|
\le 10|z|^{-4}|u|^2.
\]

Thus remote packet shells carry a `2^{-4n}` pressure weight while their
three-dimensional packing grows only as `2^{3n}`.  The spare `2^{-n}` is
summable.  In particular, under a packet packing bound, absence of a fresh
packet of critical mass `mu_*` gives `W_far<=C_far mu_*`.

For the local pressure source, Calderon--Zygmund plus a frequency-localized
Bernstein estimate gives

\[
r^{-1}\|V\|_2^2
\ge C_B^{-2}
\left(\frac{\rho}{(r\|\nabla\chi\|_\infty)C_R}\right)^{2/3}
\]

when normalized local pressure work is at least `rho`.  Consequently, in the
finite packet model,

\[
W_{cancel}\le C_{near}\mu_*^{3/2}+C_{far}\mu_*
\]

and pressure cancellation at least `rho` forces a positive fresh critical-mass
threshold.  The continuum kernel decay is exact; the remaining PDE task is to
construct a transfer-adapted packet frame satisfying the stated packing and
near/far decomposition with summable errors.

## 6. Moving windows and transport geometry

The term

\[
\int e\,\partial_t\chi+\int eU\cdot\nabla\chi
\]

suggests choosing the packet window approximately transported by the resolved
velocity, so the largest resolved advective leakage cancels. This is not yet
used as a theorem hypothesis: the window may stretch under strain, and one must
control that distortion over the packet lifetime. But it identifies the correct
transport geometry for a future spacetime grain extraction.

The remaining leakage terms are then the SGS transport `RU`, viscous boundary
flux, window-transport mismatch, and any pressure-cancellation annular charge.
The smooth-filter commutator bound from the preceding module is summable for an
expanding moat, but this does **not** by itself control the moving-window Taylor
error.  The curvature-balanced moat theorem gives the combined law
`a/M+b kappa M` and requires an adaptive spatial transition width.

## 7. Remaining bridge

To turn the exact ledger into a PDE closure one still needs to prove that a
near-extremal packet lineage admits moving windows for which the positive
leakage and pressure-cancellation annular charges are summable or force a
Bellman/fresh critical-mass event. The identity itself has no pressure tail or
packet-count loss.


## Affine-window compatibility

For the affine Gaussian grain window

\[
\chi_{L,M}(x)=\chi_0(L^{-1}(x-X)/M),
\]

the shell/aspect theorem gives the clean physical gradient bound

\[
\boxed{
N^{-1}\|\nabla\chi_{L,M}\|_\infty
\le 3C_\chi/(2M).
}
\]

Therefore every pressure-boundary estimate in this note that depends linearly
on `||grad chi||_inf` can be transplanted to the affine window without an aspect
penalty.  This observation does not by itself close the local pressure-Hessian,
`RU`, viscous-boundary or packet-overlap coefficients; it only removes the
geometric obstruction from the window shape.
