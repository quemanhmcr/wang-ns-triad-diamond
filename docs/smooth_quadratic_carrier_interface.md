# Smooth quadratic-carrier interface: read the energy law at `Q^2`

## Status

Candidate theorem awaiting exact-SHA certification on the integration branch.
The draft PR's earlier green runs are historical input only; they do not certify
the observer-gauge correction introduced here.

The candidate status string is

`EXACT_SMOOTH_QUADRATIC_CARRIER_INTERFACE__Q2_ENERGY_LAW__COMMON_GAUGE_QUOTIENT_BEFORE_PHYSICAL_RELINK__SYMMETRIC_WORK_EXISTING_STRAIN__COEFFICIENT_OBSTRUCTION_ENERGY_REENTRY`.

This theorem repairs two representation mismatches in the continuum spine.
First, the propagated PDE carrier is a smooth non-idempotent envelope `Q`, so
its energy effect is `Q^2`, not `Q`.  Second, a moving analysis partition is not
automatically a physical relink law.  Observer motion is quotiented before any
causal owner is assigned.

The event reader remains the hard orthogonal projector `P` with `QP=P`.  The
hard-projector donor theorem remains an event theorem; the smooth theorem below
uses its own `Q^2` energy law and never borrows hard idempotence.

## 1. Hard events and smooth carriers keep different jobs

At an actual nonlinear work event, use the exact hard Fourier/helical role
projector `P`.  It disintegrates physical energy and work without synthesis
error.  Between events, choose a real scalar smooth Fourier multiplier `Q`
such that

\[
0\le Q\le I,
\qquad Q^*=Q,
\qquad QP=P.
\]

The hard boundary is never differentiated.  The smooth envelope is never
promoted to a second transfer probability law.

To account for the complete energy, place `Q` in a smooth quadratic analysis
partition

\[
\boxed{
\sum_a A_a^2=I,
\qquad A_a^*=A_a,
\qquad 0\le A_a\le I.
}
\]

For one selected envelope the clean two-role construction is made before the
proof: choose a smooth angle `theta` with `theta=0` on the hard-role plateau and
`theta=pi/2` outside the outer support, then put

\[
Q=\cos\theta,
\qquad R=\sin\theta.
\]

Both multipliers are smooth and

\[
Q^2+R^2=I.
\]

Transport both by the same coherent affine dual flow.  Since they solve the
same scalar transport equation, the quadratic partition remains exact.  In
differential form,

\[
\boxed{
\sum_a\partial_t(A_a^2)
=\sum_a(\dot A_aA_a+A_a\dot A_a)=0.
}
\]

This is the same analysis--synthesis geometry already used by the canonical
smooth Littlewood--Paley frame.  Overlap is part of the analysis; it is not an
error to be paid once per carrier.

## 2. Direct smooth-carrier Navier--Stokes energy law

On a smooth pre-singular interval write Leray Navier--Stokes as

\[
\partial_tu+\mathcal B(u,u)=\nu\Delta u,
\qquad \nabla\cdot u=0.
\]

For one analysis operator `A=A(t,D)`, put

\[
w=Au,
\qquad \eta=A^2.
\]

Because `A` is self-adjoint,

\[
2\operatorname{Re}\langle Au,\dot A u\rangle
=\langle u,\dot\eta u\rangle.
\]

Because scalar Fourier `A` commutes with `Delta`, direct differentiation gives

\[
\boxed{
\frac d{dt}\|Au\|_2^2
+2\nu\|\nabla Au\|_2^2
=\langle u,\dot\eta u\rangle
-2\operatorname{Re}\langle\eta u,\mathcal B(u,u)\rangle.
}
\]

This identity is upstream of every resolved cutoff, moving-role commutator and
causal label.  It is the native conservation law which those representations
must reproduce.

Summing over the complete quadratic partition recovers the usual physical
energy identity exactly:

\[
\sum_a\|A_au\|_2^2=\|u\|_2^2,
\qquad
\sum_a\langle u,\partial_t(A_a^2)u\rangle=0.
\]

## 3. Exact resolved/HH/interface repartition

Fix the strict resolved transporter first:

\[
V=S_{N/4}u,
\qquad h=u-V,
\qquad
\mathcal L_Vf=\mathcal B(V,f)+\mathcal B(f,V).
\]

The exact quadratic identity is

\[
-\mathcal B(u,u)
=\mathcal B(V,V)-\mathcal B(h,h)-\mathcal L_Vu.
\]

Therefore the direct carrier energy law becomes

\[
\frac d{dt}\|A_au\|_2^2
+2\nu\|\nabla A_au\|_2^2
=W_{LL,a}+W_{HH,a}+J_a,
\]

where

\[
W_{LL,a}
=2\operatorname{Re}\langle A_a^2u,\mathcal B(V,V)\rangle,
\]

\[
\boxed{
W_{HH,a}
=-2\operatorname{Re}\langle A_a^2u,\mathcal B(h,h)\rangle,
}
\]

and the native moving-interface work is

\[
\boxed{
J_a
=\langle u,\partial_t(A_a^2)u\rangle
-2\operatorname{Re}\langle A_a^2u,\mathcal L_Vu\rangle.
}
\]

For the selected outer carrier, the already certified support moat gives
`A_a B(V,V)=0`, hence also `A_a^2 B(V,V)=0`.  Thus `W_LL,a=0`.

The HH term is weighted by the actual energy effect `a_a^2`.  This is exactly
the `q^2`-weighted physical HH work used by the smooth material-carrier relay.

## 4. What the outer-role commutator means at energy level

The exact outer-role equation is

\[
(\partial_t+\mathcal L_V-\nu\Delta)(Au)
=A\mathcal B(V,V)-A\mathcal B(h,h)
+(\dot A+[\mathcal L_V,A])u.
\]

Pairing this equation with `Au` gives the work of the Heisenberg term

\[
2\operatorname{Re}
\langle Au,(\dot A+[\mathcal L_V,A])u\rangle.
\]

But the role equation also has the diagonal resolved work

\[
2\operatorname{Re}\langle Au,\mathcal L_VAu\rangle
\]

on its left-hand side.  The elementary identity

\[
\boxed{
\begin{aligned}
&2\operatorname{Re}
\langle Au,(\dot A+[\mathcal L_V,A])u\rangle
-2\operatorname{Re}\langle Au,\mathcal L_VAu\rangle
\\
&\qquad
=\langle u,\partial_t(A^2)u\rangle
-2\operatorname{Re}\langle A^2u,\mathcal L_Vu\rangle
=J_A
\end{aligned}
}
\]

shows the correct handoff.

The commutator work alone is not the native smooth-interface energy work.  It
must be recombined with the diagonal role work already present in the same
outer equation.  This recombination is exact and returns the `A^2`-weighted
Navier--Stokes law.

## 5. Quotient common observer transport before physical relink

Take the adjoint split of the actual resolved linearized operator,

\[
\mathcal L_V=K+S,
\qquad K^*=-K,
\qquad S^*=S.
\]

The identity `sum_a A_a^2=I` by itself does **not** make `dot A_a` physical.
A user may vary an overlapping partition while the Navier--Stokes state is
unchanged; the channel energies then exchange mass even though nature has done
nothing.  Such motion must not become a causal owner merely because the total
exchange sums to zero.

The canonical smooth-role motion is instead tied to the actual common
transport.  Let `G` be the skew generator of the common affine/Kelvin transport
used to propagate every smooth role.  Before charging any relink work require

\[
\boxed{
\dot A_a+[G,A_a]=0
\qquad\text{for every }a.
}
\]

Equivalently, for `eta_a=A_a^2`,

\[
\partial_t\eta_a+[G,\eta_a]=0.
\]

This is the exact observer-gauge quotient.  Its channel energy work cancels:

\[
\boxed{
\langle u,\partial_t\eta_a u\rangle
-2\operatorname{Re}\langle\eta_a u,Gu\rangle=0.
}
\]

Only after this identity has been verified do we decompose the actual skew
resolved operator as

\[
K=G+K_{phys}.
\]

Then the native interface row becomes

\[
\boxed{
J_a
=-2\operatorname{Re}\langle\eta_a u,K_{phys}u\rangle
-2\operatorname{Re}\langle\eta_a u,Su\rangle.
}
\]

The first term is now eligible for physical relink ownership because the common
observer transport has already been removed.  With `v_a=eta_a u`, define

\[
T_{ab}^{phys}
=-2\operatorname{Re}\langle v_a,K_{phys}v_b\rangle.
\]

Then

\[
T_{ab}^{phys}=-T_{ba}^{phys},
\qquad
\sum_bT_{ab}^{phys}
=-2\operatorname{Re}\langle\eta_a u,K_{phys}u\rangle,
\]

and therefore the total physical relink work is zero.  Internal circulation is
redistribution, not generation.

The symmetric synthesis-pair work

\[
D_{ab}=-2\operatorname{Re}\langle v_a,Sv_b\rangle
\]

satisfies

\[
D_{ab}=D_{ba},
\qquad
\sum_bD_{ab}
=-2\operatorname{Re}\langle\eta_a u,Su\rangle,
\]

and reconstructs the same global resolved strain/deformation work.

A useful anti-example is an arbitrary time-dependent angle partition
`Q=cos(theta(t))`, `R=sin(theta(t))` with `L_V=0`.  Its channel energies can move
while total energy remains fixed.  Unless that motion is generated by the
certified common `G`, the theorem refuses to form a physical relink certificate.
This is the permanent guard against charging nature for a change of observer.

## 6. Why `I-Q` is the wrong smooth complement

For a hard projector, `Q^2=Q` and `I-Q` is its orthogonal energy complement.
For a smooth envelope this is false:

\[
Q^2+(I-Q)^2
=I-2Q(I-Q)\ne I.
\]

For `K^*=-K`, define the commutator work

\[
\mathcal I_A(K)
=2\operatorname{Re}\langle Au,[K,A]u\rangle.
\]

Direct algebra gives

\[
\boxed{
\mathcal I_Q(K)+\mathcal I_{I-Q}(K)
=4\operatorname{Re}\langle Q(I-Q)u,Ku\rangle.
}
\]

It need not vanish.  For

\[
K=\begin{pmatrix}0&-1\\1&0\end{pmatrix},
\qquad
Q=\operatorname{diag}(1/2,0),
\qquad
u=(1,1),
\]

the defect equals `-1`.

If instead `R` is chosen with `Q^2+R^2=I`, then

\[
\boxed{
\mathcal I_Q(K)+\mathcal I_R(K)=0.
}
\]

This counterexample is a permanent guard against silently importing hard-role
orthogonality into a smooth PDE carrier.

## 7. Hard-event registration remains exact

The selected hard event role lies in the plateau where `Q=1`.  Hence

\[
QP=P
\]

and for every terminal probe `phi`,

\[
\boxed{
\langle Pu,\phi\rangle
=\langle Qu,P\phi\rangle.
}
\]

The quadratic complement vanishes on that plateau.  No coefficient is lost and
no second event measure is created.  At the next actual nonlinear interaction,
hard Fourier/helical roles are read again from physical work; they are not
propagated through the slab.

## 8. Coefficient obstruction triggers energy reentry, not ownership

Suppose the backward adjoint identity encounters a large coefficient impulse

\[
I_{interface}
=\int\langle\psi,R_Q\rangle dt.
\]

Its size is an amplitude fact.  It is not automatically physical work.  The
correct semantics is:

> the coefficient impulse locates the first interval on which clean coefficient
> continuation fails; it does not name the physical causal owner.

At that first interval, use the actual terminal carrier energy

\[
E_1=\|Q(t)u(t)\|_2^2
\ge\frac{|z(t)|^2}{\|\psi(t)\|_2^2},
\]

the actual initial carrier energy `E_0`, and the positive part of the native
interface work `J_Q`.  Reenter the existing physical-energy gate.

On the low-strain branch it gives exactly one of the already certified physical
alternatives:

1. `E_0>=E_1/5`: material/carrier energy inheritance;
2. `W_interface^+>=E_1/5`: actual interface physical work;
3. otherwise
   \[
   W_{HH}^+\ge\frac8{15}E_1.
   \]

Only in alternative 2 is the interface Hahn law consulted.  Since pointwise

\[
J=J^{rel}+J^{str},
\qquad
[J]_+\le[J^{rel}]_++[J^{str}]_+,
\]

positive native interface work is covered by conservative relink and existing
strain work.  One carries at least half of the positive interface law; exact
ties remain joint.

The magnitude of `I_interface` appears in no energy threshold and in no causal
probability.  This is the same rule already imposed on HH Duhamel amplitude:

\[
\boxed{
\text{Duhamel locates support/obstruction; physical energy work supplies causality.}
}
\]

## 9. Relation to the resolved donor/circulation quotient

The certified resolved-interface donor/circulation quotient begins from actual
resolved low--high work on a complete hard event-role partition.  It then traces
positive skew gain to simultaneous donor roles and quotients role circulation
from recursive depth.

The present theorem is complementary and logically upstream at a different
representation layer:

- this theorem supplies the exact energy law of the smooth propagated PDE
  carrier and prevents a non-idempotent `Q` from borrowing projector algebra;
- the donor quotient supplies same-event hard-role provenance once actual
  resolved work is disintegrated at an event.

The two physical measures are not identified by name.  Smooth commutator work,
native `Q^2` carrier-interface work and hard event-role resolved work retain
their exact definitions.  Any future composition must use an explicit physical
work pushforward, not a theorem-name substitution.

## 10. Analytic scope

The identities are first stated on a smooth pre-singular Navier--Stokes interval
with smooth compact scalar multipliers, where all differentiated pairings are
classical.  Equivalently one may prove them on spectral Galerkin truncations;
the formulas are exact at every truncation and pass by the regularity available
on the supplied smooth interval.

No discontinuous boundary is differentiated.  No limiting hard approximation
to `Q` is used.  The smooth angle partition may be chosen once and transported
by the same affine dual flow as the selected carrier.

## 11. What this closes and what remains

This candidate closes the following local seam if exact-SHA CI confirms it:

- the propagated envelope energy is read at `Q^2` rather than pretending `Q`
  is a hard projector;
- smooth overlap is completed by a square partition;
- common affine/Kelvin role motion is proved and quotiented as observer gauge
  before physical ownership;
- only the residual skew operator `K_phys` is eligible for conservative physical
  relink;
- symmetric work remains existing strain/deformation;
- a coefficient obstruction is only a first-stop locator and cannot enter the
  canonical master until actual energy/work reentry resolves it.

It does not prove that repeated strain, relink, service or donor owners terminate
globally.  It does not close the UV-unbounded full-survivor alternative, and it
does not prove Navier--Stokes regularity.

## 12. Certification

The corrected observer-gauge theorem has not yet been promoted.  Its dedicated
workflow and the full physical-energy causal integration workflow are wired to
run on the exact implementation SHA.  Certification data and stored artifacts
will be recorded here only after those exact-SHA runs succeed.

The earlier draft-PR runs remain useful regression history for the underlying
`Q^2` algebra, but they are **not** used as certification evidence for the new
`dot A+[G,A]=0` gauge quotient or the canonical master obstruction barrier.
