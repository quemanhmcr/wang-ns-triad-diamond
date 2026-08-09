# High-frequency service dissipation: energy forces inheritance or regeneration

## The false shortcut

The high-frequency exit in coherent increment service is a normalized enstrophy law, not resolved low-pass dissipation.  On hard dyadic annuli

\[
M_j=2^jN,\qquad M_j/2<|\xi|\le M_j,\qquad j\ge1,
\]

set

\[
\mu_j=M_j\|P_ju\|_2^2,
\qquad
D_>^{LP}=\int\sum_{j\ge1}2^j\mu_j\,d\tau,
\qquad \tau=N^2t.
\]

There is **no** direct implication `D_high -> critical shell mass`.  At one arbitrarily high level `j`, choose

\[
\mu_j=2^{-j}D.
\]

Then `2^j mu_j=D` while `mu_j->0`.  Any proof which silently calls `D_high` resolved `D_V`, or extracts a uniform shell mass directly from it, is imposing a false geometry.

## The physical identity that fixes the problem

Let

\[
w=P_{>N}u.
\]

On the `j`th hard annulus,

\[
(M_j/2)^2\|P_ju\|_2^2
\le \|\nabla P_ju\|_2^2
\le M_j^2\|P_ju\|_2^2.
\]

Since `M_j/N=2^j` and `d tau=N^2dt`, summing gives

\[
\boxed{
\frac14D_>^{LP}
\le N\int\|\nabla w\|_2^2dt
\le D_>^{LP}.
}
\]

The hard Fourier projection is orthogonal and commutes with viscosity and Leray projection.  Its exact energy identity is

\[
\|w(t)\|_2^2-\|w(s)\|_2^2
+2\nu\int_s^t\|\nabla w\|_2^2dt
=
\int_s^t2\Re\langle w,F_>\rangle dt,
\]

with

\[
F_>=-P_{>N}\mathbb P\nabla\cdot(u\otimes u).
\]

Define actual positive nonlinear work

\[
W_>^+=\int_s^t2[\Re\langle w,F_>\rangle]_+dt.
\]

Dropping the nonnegative terminal energy gives

\[
\boxed{
N\|w(s)\|_2^2+NW_>^+
\ge \frac\nu2D_>^{LP}.
}
\]

Thus at least one of two **physical owners** carries

\[
\boxed{\frac\nu4D_>^{LP}}.
\]

Exact equality keeps both owners jointly.

## Owner 1: inherited tail energy is a real critical shell

At the initial slice,

\[
N\|w(s)\|_2^2
=
\sum_{j\ge1}\frac{N}{M_j}\mu_j(s)
=
\sum_{j\ge1}2^{-j}\mu_j(s).
\]

The dyadic weights sum to one.  Therefore

\[
\sup_j\mu_j(s)\ge N\|w(s)\|_2^2.
\]

On the inherited owner branch there is an actual shell with

\[
\boxed{
M_j\|P_ju(s)\|_2^2\ge\frac\nu4D_>^{LP}.
}
\]

This is precisely the deterministic input of the generic critical-shell theorem.  No causal probability and no packet label is introduced.  The generic theorem retains its observed-history guard and may return a named strain/interface/HH stop, `t=0`, or own-scale service.

## Owner 2: actual nonlinear regeneration work

Write `w=sum_j w_j` with `w_j=P_ju`.  Orthogonality gives the signed work as the sum of shell signed works, hence

\[
W_>^+\le\sum_jW_j^+,
\qquad
W_j^+=\int2[\Re\langle w_j,F_j\rangle]_+dt.
\]

Because every high shell has `M_j/N=2^j>=2`,

\[
\boxed{
\sum_jM_jW_j^+
\ge2N\sum_jW_j^+
\ge2NW_>^+.
}
\]

So the regeneration owner supplies an actual positive **own-scale shell-work law** of total mass at least

\[
\frac\nu2D_>^{LP}.
\]

This remains work, not amplitude.

For each shell choose the resolved field `V=S_(M_j/4)u` and `h=u-V`.  The low-low source is support-excluded:

\[
P_j\mathbb P\nabla\cdot(V\otimes V)=0,
\]

because `V tensor V` lies at frequencies at most `M_j/2` while the hard annulus is strictly above `M_j/2`.  Therefore the exact shell work is the sum of

- high-high work from `h tensor h`;
- mixed resolved/high cross work, equivalently the physical transport/interface/strain owner in the outer-role formulation.

For signed densities `r=r_HH+r_I`,

\[
[r]_+\le[r_{HH}]_++[r_I]_+.
\]

After summing shell-time atoms, at least one of HH or resolved-interface positive work carries

\[
\boxed{\frac\nu4D_>^{LP}}.
\]

This is a physical work statement.  A large interface branch is not made free.  A large HH branch is not automatically declared to satisfy the generated-energy condition `W_HH>=8E_1/15`; it must pass the existing physical energy gate before KL productivity is used.

## Master-facing route

The high-frequency coherent-service exit now reads

\[
\boxed{
D_{high}
\longrightarrow
\text{inherited critical shell}
\quad\text{or}\quad
\text{actual positive HH/interface regeneration}.
}
\]

No branch is converted to a globally bounded reset.  No high enstrophy is relabeled resolved `D_V`.  No frozen packet is introduced.

The remaining issue is not the unit mismatch of `D_high`; it is continuation of the actual regeneration owner through its own energy/interface first-stop logic, together with the separate low-frequency pressure-reservoir lineage.

No 3D Navier--Stokes global-regularity conclusion is asserted.
