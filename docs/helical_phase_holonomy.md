# Relative helical phase holonomy on a reuse diamond

For a forward transfer edge `x+y=z`, let `theta_v` be the phase of the physical
helical amplitude at mode `v`.  With the repository reality convention, the
phase entering the child-energy transfer is

\[
\delta_e
=\arg g_e-\theta_x-\theta_y+\theta_z-\tau_e,
\]

where `tau_e` is the forward-transfer target (`0` or `pi` in the present sign
normalization).  This is invariant under any phase change of the helical basis.

Consider the reuse diamond

\[
a+b=m,\qquad m+c=d,
\]
\[
b+c=n,\qquad a+n=d.
\]

Then all six modal phases cancel **exactly**:

\[
\boxed{
\delta_{abm}+\delta_{mcd}-\delta_{bcn}-\delta_{and}
\equiv H_\phi\pmod{2\pi},
}
\]

with

\[
H_\phi=
\arg g_{abm}+\arg g_{mcd}-\arg g_{bcn}-\arg g_{and}
-(\tau_{abm}+\tau_{mcd}-\tau_{bcn}-\tau_{and}).
\]

`H_phi` is a genuine observable of the reuse motif.  A different helical gauge
at any physical Fourier mode changes two incident coupling phases with opposite
signs and cancels from `H_phi`.

The triad-normal gauge identifies the geometry behind this invariant.  Each
triad has a constant quadrature coupling phase, while a mode reused by two
incident triads acquires the spin-one transition factor

\[
e^{-is\psi},
\]

where `psi` is their signed dihedral angle about that mode.  Thus `H_phi` is a
relative **incidence spin holonomy**: a signed combination of dihedral spin
rotations plus the constant sign/target phases.  A rigid rotation of the whole
diamond changes none of it.

There is a sharp phase-lock cost.  For principal `|H|<=pi`, among all four
residuals satisfying

\[
\delta_1+\delta_2-\delta_3-\delta_4\equiv H\pmod{2\pi},
\]

the minimum of the unweighted polarization loss is

\[
\boxed{
\sum_{i=1}^4(1-\cos\delta_i)
\ge
4\left(1-\cos\frac{|H|}{4}\right).
}
\]

The bound is sharp, attained by
`(H/4,H/4,-H/4,-H/4)`.  A direct proof classifies the stationary points of
`sum cos(delta_i)` on the compact constraint torus: each sine is equal, so each
angle is either `alpha` or `pi-alpha`; for `|H|<=pi` the all-equal lift has value
`4 cos(H/4)>2`, strictly above every mixed branch.

If every diamond edge has capacity weight at least `beta` and multiplier
`m_e>=1-eta`, the existing polarization identity therefore gives

\[
\boxed{
D_{phase}
\ge
4\beta(1-\eta)
\left(1-\cos\frac{|H_\phi|}{4}\right).
}
\]

At the repository good threshold `eta=10^-4`, Arb certifies that the concrete
classification `|H_phi|>=1/5` pays

\[
\boxed{D_{phase}\ge \beta/250.}
\]

This is a new finite-packet no-escape branch.  It does not yet assert that every
continuum packet lineage contains a diamond with a uniform lower edge weight;
that remains part of the spatial packet/reuse extraction.
