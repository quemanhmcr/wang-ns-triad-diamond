# Curvature-balanced moving moat

This note corrects an overly optimistic localization heuristic in the previous
ledger.  A wider spatial moat reduces the smooth-filter commutator, but it makes
the affine Taylor approximation of the velocity field worse.  The two effects
must be balanced.

## 1. Two localization errors pull in opposite directions

At frequency `N`, use a material transition radius

\[
R=\frac{M}{N}.
\]

The smooth-filter commutator from the previous SGS module has normalized form

\[
E_{comm}\le\frac{a}{M}.
\]

For an affine-deforming window, the exact cancellation leaves the Taylor
remainder

\[
|(\partial_t+u\cdot\nabla)\chi|
\le\frac12\|\nabla^2u\|_\infty R
\|F\|^2\|F^{-1}\|\|\nabla\chi_0\|_\infty.
\]

On a parabolic packet lifetime `T=cN^{-2}`, introduce the scale-invariant
curvature

\[
\kappa=N^{-3}\|\nabla^2u\|_\infty.
\]

Then the integrated Taylor factor has the form

\[
E_{curv}\le b\kappa M,
\qquad
b=\frac c2\|F\|^2\|F^{-1}\|\|\nabla\chi_0\|_\infty,
\]

up to the already-separated local energy/transfer normalization.

Thus the geometric localization ledger is

\[
\boxed{
E(M)=\frac aM+b\kappa M.
}
\]

## 2. Exact optimal balance

AM--GM, or one derivative, gives

\[
\boxed{
M_*=\sqrt{\frac{a}{b\kappa}},
\qquad
E_*=2\sqrt{ab\kappa}.
}
\]

This has a clear physical interpretation: the packet transition region is set by
a local balance between filter nonlocality and curvature of the advecting flow,
not by generation number alone.

## 3. Countermodel to the old fixed quadratic moat schedule

The previous commutator discussion used

\[
M_j=(j+3)^2,
\]

which indeed gives `sum 1/M_j<infinity`.  But take the perfectly decaying
curvature sequence

\[
\kappa_j=(j+3)^{-3}.
\]

Then

\[
\sum_j\frac1{M_j}<\infty,
\qquad
\sum_jM_j\kappa_j
=\sum_j\frac1{j+3}=\infty.
\]

So the statement "the same expanding moat schedule makes all spatial
localization errors summable" is false.  This is a structural countermodel and
is now recorded in the ledger.

The balanced choice instead gives

\[
M_j=\kappa_j^{-1/2}=(j+3)^{3/2},
\]

and both errors are

\[
(j+3)^{-3/2},
\]

which are summable.

## 4. Curvature failure forces critical mass in a band-limited packet model

In the unitary Fourier convention, for a scalar component supported in
`|xi|<=Lambda N`, Cauchy--Schwarz gives

\[
\|\partial_{ij}f\|_\infty
\le
(2\pi)^{-3/2}
\left(\frac{4\pi}{7}\right)^{1/2}
(\Lambda N)^{7/2}\|f\|_2.
\]

Writing

\[
\mu=N\|f\|_2^2,
\]

this is

\[
\boxed{
N^{-3}\|\partial_{ij}f\|_\infty
\le C_B\sqrt\mu,
\qquad
C_B=(2\pi)^{-3/2}\sqrt{4\pi/7}\,\Lambda^{7/2}.
}
\]

Hence a curvature level `kappa` forces

\[
\mu\ge(\kappa/C_B)^2.
\]

If even the **optimally balanced** moat costs at least `eta`, then

\[
2\sqrt{ab\kappa}\ge\eta
\quad\Longrightarrow\quad
\kappa\ge\frac{\eta^2}{4ab},
\]

so

\[
\boxed{
\mu
\ge
\frac{\eta^4}{16a^2b^2C_B^2}.
}
\]

Thus failure of the curvature-balanced localization branch is not neutral: in a
frequency-localized packet model it forces a definite scale-critical energy
charge, which can be handed to the fresh/reuse ancestry ledger.

## 5. Consequence for the program

Spatial moat widths must now be chosen **adaptively from the local curvature**.
The remaining continuum work is to build this adaptive choice into the nested
packet partition without breaking transfer weights, and to distinguish whether
the forced critical mass is fresh or already belongs to the reused ancestry.
\n\n## Helical localization renormalizes the same curvature coefficient\n\nThe localized relative-polarization theorem shows that, on the low-strain
lifetime branch `c sigma_0<=1/30`, carrier-direction and triad-normal frame
variation contribute

\[
E_{pol}\le3h+\frac{15}{2}c\,\kappa M.
\]

Therefore the physical spatial optimizer remains of the same form.  If the
previous non-helical coefficient is `b`, use

\[
E_{loc}(M)\le \frac aM+(b+15c/2)\kappa M+3h.
\]

The optimizer is still square-root in curvature.  No separate polarization moat
schedule is introduced.\n

## Affine ellipsoidal successor

The isotropic calculation has an affine version.  For a grain frame `L` and
`z=L^{-1}(x-X)`, transporting `Xdot=U(X)`, `Ldot=A(X)L` gives

\[
D_t z=L^{-1}[U(X+Lz)-U(X)-A(X)Lz].
\]

Thus the relevant curvature is the grain-normalized tensor

\[
\kappa_{aff}=\sup\|L^{-1}(\nabla^2U)[L,L]\|.
\]

The shell certificate controls the smallest physical grain axis, so an
ellipsoidal window `chi_0(z/M)` still has normalized physical gradient `O(1/M)`.
For a convolution filter the commutator bound is exactly

\[
\|[\chi,G_N*]f\|_2
\le \frac{3m_1(G)C_\chi}{2M}\|f\|_2.
\]

The material Taylor remainder is `O(kappa_aff M)`.  Hence the same square-root
balance survives in affine coordinates.  See `docs/affine_window_balance.md`.
