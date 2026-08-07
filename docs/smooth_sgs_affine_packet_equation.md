# Smooth-SGS affine packet equation: microscopic roles versus macroscopic energy transport

The smooth-SGS transfer block uses Fourier/helical roles of the **full velocity**
`u`, while the exact localized SGS energy ledger uses the coarse-grained field
`U=G_N*u`.  These are two different levels of description and must not be
collapsed into one forcing term.

## 1. Exact microscopic role equation

Write Navier--Stokes in Leray form

\[
\partial_tu+\mathbb P\nabla\!\cdot(u\otimes u)=\nu\Delta u,
\qquad \nabla\cdot u=0.
\]

Let `Q_i(t,D)` be a smooth shell/cell/helicity multiplier whose range is
solenoidal and which commutes with spatial derivatives.  Put

\[
w_i=Q_i u.
\]

Then exactly

\[
\boxed{
\partial_t w_i
=-Q_i\mathbb P\nabla\!\cdot(u\otimes u)
+\nu\Delta w_i
+(\partial_tQ_i)u.
}
\]

There is **no pressure forcing** in this equation: the Leray projector removes it
before packetization.  Reintroducing a separate pressure force into `F_i` would
double count the pressure already encoded by the solenoidal nonlinear operator.

## 2. A strict low-pass produces the physical Kelvin transporter

On the signed-good extremal core, with child frequency normalized to `N`, the
smaller parent satisfies

\[
\frac{|k_{parent}|}{N}
\ge r_*e^{-1/80} > \frac35.
\]

The final inequality is Arb-certified.  Choose a smooth transport field

\[
V=S_{N/4}u,
\qquad \operatorname{supp}\widehat V\subset\{|\xi|\le N/4\},
\]

and write `u=V+h`.  Since

\[
\operatorname{supp}\widehat{V\otimes V}
\subset\{|\xi|\le N/2\},
\]

one has for every selected role

\[
\boxed{Q_i\mathbb P\nabla\cdot(V\otimes V)=0.}
\]

Hence the selected-role nonlinearity separates **exactly** into

\[
\boxed{
\text{low--high transport}
+\text{high--high transfer}.
}
\]

This is the physically correct origin of the Kelvin/rapid-distortion background:
it is the strict low-frequency part of the same Navier--Stokes velocity, not an
externally imposed affine field.

If `w_j,w_k` are the selected companion roles, define the designated high--high
triad source

\[
\mathcal T_i^{jk}
=-Q_i\mathbb P\nabla\cdot
(w_j\otimes w_k+w_k\otimes w_j).
\]

All other high--high interactions remain a genuine packet-forcing/cross-error
term.  Near-extremal transfer extraction must make their **transfer contribution**
small; this note does not assume an unproved `L^2` forcing bound for them.

## 3. Moving frequency cells: the correct residual is a Heisenberg commutator

For scalar affine transport `L_A=(Ax)\cdot\nabla`, `tr A=0`, a Fourier multiplier
`m(t,D)` obeys

\[
\widehat{[L_A,m(D)]f}
=-(A^T\xi)\cdot\nabla_\xi m\,\widehat f.
\]

Therefore the packet-cell motion enters through

\[
\boxed{
\partial_tm-(A^T\xi)\cdot\nabla_\xi m,
}
\]

not through `partial_t m` alone.  A multiplier transported by the dual affine
flow solves

\[
\partial_tm-(A^T\xi)\cdot\nabla_\xi m=0,
\]

so common affine carrier/covariance motion creates **zero** spectral-cell forcing.
The vector/helical part of the same statement is the objective polarization and
common-`SL(2)` gauge already proved in the helical modules.

## 4. Two-level no-double-counting rule

The microscopic full-velocity role equation contains

- the designated high--high triad source;
- other high--high interactions;
- cross-cell low--high interactions;
- moving-projector Heisenberg residual;
- non-affine third-Hermite forcing;
- relative-polarization curvature;
- spatial multiplier/Leray commutators if the role is strongly windowed;
- viscous window terms if one evolves the strongly localized field.

By contrast, the **resolved SGS energy equation**

\[
\partial_te+\nabla\cdot[(e+P)U+RU-\nu\nabla e]
=-\Pi-\nu|\nabla U|^2
\]

contains `RU`, pressure boundary work, resolved window transport and the resolved
viscous boundary flux.  These macroscopic terms must not be inserted again into
`F_i` of the microscopic helical roles.

The coupling of the two levels is the already-certified trilinear SGS transfer
measure: smooth-symbol freezing identifies the selected full-velocity triads
with the physical resolved flux up to the existing summable `O(h)` multiplier
error.

This separation is exact bookkeeping, not a turbulence closure.
