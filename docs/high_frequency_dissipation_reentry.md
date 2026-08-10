# High-frequency service dissipation: physical tail energy forces inheritance or regeneration

## Keep the LP observable and the PDE currency distinct

The high-frequency exit in coherent increment service is a **standard Littlewood--Paley normalized-enstrophy observable**.  Write it schematically as

\[
D_{high}=\int\sum_{j\ge1}2^j\mu_j\,d\tau,
\qquad \mu_j=M_j\|u_j\|_2^2,
\qquad M_j=2^jN.
\]

This is not resolved low-pass `D_V`.  It is also not literally the orthogonal hard-tail gradient unless a particular LP partition has been fixed and its comparison constants have been certified.

That distinction matters.  There is no direct implication `D_high -> critical shell mass`: if all high currency sits at level `j` and

\[
\mu_j=2^{-j}D,
\]

then `2^j mu_j=D` while `mu_j->0` as `j->infinity`.  High normalized enstrophy by itself contains frequency leverage, not a scale-independent mass floor.

The theorem therefore uses the actual PDE currency

\[
\boxed{
D_{tail}:=N\int_s^t\|\nabla P_{>N}u\|_2^2\,dt.
}
\]

A smooth LP high observable may enter only through an explicit, fixed comparison

\[
\boxed{
D_{tail}\ge c_{LP}D_{high}.
}
\]

For high multipliers with `supp phi_j subset {|xi|>=a M_j}` and `sum_j|phi_j|^2<=B`, Plancherel gives the exact formula

\[
\boxed{c_{LP}=a^2/B.}
\]

The canonical smooth analysis--synthesis frame registered in coherent increment service has `a=1/2` and `B=1`, hence `c_LP=1/4`.  No hidden equality between the LP observable and the hard tail is assumed.

For the auxiliary **hard orthogonal annuli**

\[
M_j/2<|\xi|\le M_j,
\qquad M_j=2^jN,
\]

put `mu_j=M_j||P_j u||_2^2` and

\[
D_>^{hard}=\int\sum_{j\ge1}2^j\mu_j\,d\tau.
\]

Plancherel gives the exact spectral comparison

\[
\boxed{
\frac14D_>^{hard}\le D_{tail}\le D_>^{hard}.
}
\]

Thus the auxiliary hard-annulus supplier has the same clean constant `1/4`.  The canonical smooth square-LP frame already supplies `c_LP=1/4`; only an alternative LP frame needs a different certified constant.

## The physical identity

Let

\[
w=P_{>N}u.
\]

The hard Fourier projection is orthogonal and commutes with viscosity and Leray projection.  Its exact energy identity is

\[
\|w(t)\|_2^2-\|w(s)\|_2^2
+2\nu\int_s^t\|\nabla w\|_2^2dt
=
\int_s^t2\Re\langle w,F_>\rangle dt,
\]

where

\[
F_>=-P_{>N}\mathbb P\nabla\cdot(u\otimes u).
\]

Define the actual positive nonlinear tail work

\[
W_>^+=\int_s^t2[\Re\langle w,F_>\rangle]_+dt.
\]

Discarding the nonnegative terminal energy gives

\[
\boxed{
N\|w(s)\|_2^2+NW_>^+
\ge2\nu D_{tail}.
}
\]

Therefore at least one of two physical owners carries

\[
\boxed{
\nu D_{tail}.
}
\]

If the entrance came from a smooth LP high-frequency exit, this is at least

\[
\boxed{
\nu c_{LP}D_{high}.
}
\]

Exact equality keeps both owners jointly.

## Owner 1: inherited tail energy exposes a real critical shell

At the initial slice, using the hard orthogonal annuli,

\[
N\|w(s)\|_2^2
=
\sum_{j\ge1}\frac{N}{M_j}\mu_j(s)
=
\sum_{j\ge1}2^{-j}\mu_j(s).
\]

The weights sum to one, so

\[
\sup_j\mu_j(s)\ge N\|w(s)\|_2^2.
\]

On the inherited owner branch there is therefore an actual high shell `M_j>=2N` with

\[
\boxed{
M_j\|P_ju(s)\|_2^2\ge\nu D_{tail}
\ge\nu c_{LP}D_{high}.
}
\]

This is exactly the deterministic input of the generic critical-shell theorem.  No causal probability, packet persistence or synthetic material label is introduced.  The generic shell theorem still enforces its observed-history guard and may return strain/interface/HH first stop, `t=0`, or own-scale service.

## Owner 2: actual nonlinear regeneration work

Decompose the hard tail orthogonally as `w=sum_j w_j`, `w_j=P_j u`.  Signed tail work is the sum of signed shell works, hence positivity gives

\[
W_>^+\le\sum_jW_j^+,
\qquad
W_j^+=\int2[\Re\langle w_j,F_j\rangle]_+dt.
\]

Every hard high shell has `M_j/N=2^j>=2`.  Therefore

\[
\boxed{
\sum_jM_jW_j^+
\ge2N\sum_jW_j^+
\ge2NW_>^+.
}
\]

On the regeneration owner branch this produces an actual positive **own-scale shell-work law** with total mass at least

\[
\boxed{2\nu D_{tail}}.
\]

For each shell choose the resolved field `V=S_(M_j/4)u` and `h=u-V`.  The low--low source is excluded by support:

\[
P_j\mathbb P\nabla\cdot(V\otimes V)=0,
\]

because `V tensor V` lies at frequencies at most `M_j/2`, while the hard annulus is strictly above `M_j/2`.

Thus the signed shell nonlinear work is exactly the sum of

- high--high work from `h tensor h`;
- mixed resolved/high cross work, i.e. the resolved transport/interface/strain owner in the outer-role formulation.

For signed densities `r=r_HH+r_I`,

\[
[r]_+\le[r_{HH}]_++[r_I]_+.
\]

After summing shell-time atoms, at least one of the two physical work owners carries

\[
\boxed{
u D_{tail}}
\]

and therefore at least `nu c_LP D_high` for an LP-supplied exit.

This is a work statement, not an amplitude statement.  A large interface branch remains physical interface/strain provenance.  A large HH branch is **not** automatically the generated-energy condition `W_HH>=8E_1/15`; its child-energy gate must still be supplied before the current KL productivity theorem is invoked.

The continuation of that branch is now sharper.  Decompose the actual shell resolved operator as `L_V=K+S` before taking positive parts.  The `K` component is conservative role-to-role flux and is quotiented by same-event donor tracing; the `S` component is the already existing strain/deformation owner.  Thus resolved interface is no longer an independent recursive-generation destination.

## Master-facing route

The coherent-service high-frequency exit is now factored without a unit fiction:

\[
\boxed{
D_{high}^{LP}
\xrightarrow{\;D_{tail}\ge c_{LP}D_{high}\;}
D_{tail}^{physical}
\longrightarrow
\text{inherited critical shell}
\quad\text{or}\quad
\text{actual positive HH/interface regeneration}.
}
\]

No branch becomes an additive reset.  No high-frequency enstrophy is relabeled resolved `D_V`.  No smooth LP partition is silently identified with a hard orthogonal projector.  No frozen packet is introduced.

The remaining high-frequency question is continuation of the actual regeneration owner through its own work/energy/interface first-stop logic.  The low-frequency pressure-reservoir lineage is separate.

No 3D Navier--Stokes global-regularity conclusion is asserted.
