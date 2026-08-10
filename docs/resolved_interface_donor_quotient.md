# Resolved interface donor/circulation quotient

## Status

Candidate theorem.  Certification is pending dedicated GitHub Actions and the full physical-energy causal integration sweep.

The intended status string is

`EXACT_RESOLVED_INTERFACE_DONOR_QUOTIENT__POSITIVE_INTERFACE_TO_CONSERVATIVE_SKEW_DONOR_OR_EXISTING_STRAIN__FINITE_SAME_EVENT_DONOR_EXHAUSTION__CIRCULATION_ZERO_RECURSION_DEPTH`.

This theorem addresses the interface-owner seam in the continuum master.  Its purpose is not to create a new stopping currency.  It shows that the apparently separate resolved role-interface branch is already made of two physical objects present in Navier--Stokes: conservative redistribution between simultaneous roles and resolved strain/deformation work.

## 1. Start from the actual resolved low--high energy operator

At one shell-time event the strict resolved transporter is fixed first.  The mixed low--high Navier--Stokes term is linear in the selected high field, so write its actual energy-side operator as

\[
L_V=K+S,
\qquad K^*=-K,
\qquad S^*=S.
\]

Let the exact event roles be a complete orthogonal partition `P_a`, and put `w_a=P_a h` for the selected high field `h`.  The signed resolved work delivered to role `a` is

\[
R_a=2\operatorname{Re}\langle w_a,L_V h\rangle.
\]

Disintegrate it before taking positive parts:

\[
R_a=R_a^K+R_a^S,
\]

where

\[
T_{ab}=2\operatorname{Re}\langle w_a,K w_b\rangle,
\qquad
D_{ab}=2\operatorname{Re}\langle w_a,S w_b\rangle,
\]

and

\[
R_a^K=\sum_bT_{ab},
\qquad
R_a^S=\sum_bD_{ab}.
\]

This is the master-facing form of the theorem.  It applies directly to the high-tail mixed resolved/high work because `V_j=S_{M_j/4}u` is fixed at each shell and `B(V_j,h)+B(h,V_j)` is linear in `h`.

The already certified moving-projector commutator identity is a companion specialization.  For its interface work one obtains the same adjoint split, with the skew piece equal to conservative off-diagonal role flux and the symmetric piece equal to off-diagonal strain work.  We do **not** identify the full mixed-work observable with the commutator observable; they share the same physical operator decomposition.

At every physical work atom,

\[
R=R_K+R_S,
\]

hence

\[
[R]_+\le[R_K]_++[R_S]_+.
\]

After integrating any shell/time disintegration in the same physical work unit,

\[
W_{resolved}^+\le W_{skew}^+ + W_{strain}^+.
\]

Thus a positive resolved cross/interface owner has only two native continuations: conservative `K` redistribution or existing `S` strain/deformation work.  One carries at least half of its positive mass; exact ties remain joint.  There is no third interface currency.

## 2. The symmetric branch is already strain/deformation

Because `S*=S`,

\[
D_{ab}=D_{ba}.
\]

For the full resolved low--high work this includes both diagonal and off-diagonal strain contributions.  In the moving-projector commutator specialization, the diagonal pieces cancel from the commutator and the remaining symmetric interface is exactly the off-diagonal strain work already certified in `nonaffine_role_interface_work.py`.

Either way, `S` is the same physical symmetric deformation operator.  Therefore a large `S` branch delegates once to the existing coherent `K_coh`, high-strain, objective-source, and critical-`D_V` machinery.  It must not be charged again as an independent source merely because it appeared after a role or shell decomposition.

## 3. The skew branch is a divergence of role flux

Because `K*=-K`,

\[
T_{ab}=-T_{ba},
\qquad
T_{aa}=0.
\]

The net skew work on role `a` is

\[
R_a^K=\sum_bT_{ab}.
\]

Orient the positive physical flux by

\[
F_{b\to a}=[T_{ab}]_+.
\]

Then exactly

\[
\boxed{R_a^K=\sum_bF_{b\to a}-\sum_bF_{a\to b}.}
\]

So positive skew work on one role is energy received from other simultaneous roles.  It is not generation.  Summing over all roles gives

\[
\sum_aR_a^K=0.
\]

More generally, for every role set `C`, all internal skew transfers cancel pairwise:

\[
\boxed{
\sum_{a\in C}R_a^K
=F_{C^c\to C}-F_{C\to C^c}.
}
\]

This is the exact finite-role divergence theorem for resolved skew work.  An overall sign convention for the PDE work merely reverses every directed edge; antisymmetry and the donor theorem are unchanged.

## 4. Finite donor exhaustion

Take any role `a_0` with

\[
R_{a_0}^K>0.
\]

Build its backward donor closure `C`: whenever `a` is already in `C`, include every `b` with

\[
F_{b\to a}>0.
\]

Because the event role partition is finite, this closure stabilizes after finitely many additions.  By construction **no positive skew flux enters `C` from outside**.

Suppose, for contradiction, that every role in `C` has nonnegative net skew work.  Since `a_0` has strictly positive gain,

\[
\sum_{a\in C}R_a^K>0.
\]

But the subset divergence theorem and zero external inflow give

\[
\sum_{a\in C}R_a^K
=-F_{C\to C^c}\le0,
\]

which is impossible.

Therefore `C` contains at least one role `b` with

\[
R_b^K<0.
\]

That role is an actual energy donor for the simultaneous conservative redistribution.  Any directed donor walk reaching it can have cycles deleted, leaving a simple donor path of length at most

\[
\#\text{roles}-1.
\]

No role index is selected as canonical.  The theorem retains the entire set of reachable negative-net donors.

## 5. Why cycles are free in recursion depth

A skew role cycle can carry nonzero circulating flux while every internal contribution cancels in the divergence.  This circulation is physically real redistribution, but it does not create energy, advance physical time, or create another Navier--Stokes generation.

Therefore a chain of same-event skew donor relays is **zero recursion depth** after the physical role-flux quotient.  It may change which event role is read as the donor, but it does not create a new causal measure and it does not supply scale progress.

This is different from declaring interface “free”.  The interface work is physical.  The statement is sharper:

- its symmetric part is existing strain/deformation work;
- its skew part is conservative role-to-role energy flux whose source is another simultaneous role.

The master should charge the physical owner reached after this quotient, not the number of role labels traversed.

## 6. High-tail consequence in the unchanged causal unit

The certified high-tail regeneration theorem uses the common physical unit `N dW`, where `N` is the parent block scale.  If resolved interface is a clean high-tail owner, then

\[
W_{interface}^+\ge \frac{\nu D_{tail}}2.
\]

The exact interface split gives

\[
W_{skew}^++W_{sym}^+\ge W_{interface}^+.
\]

Hence

\[
\boxed{
W_{skew}^+\ge\frac{\nu D_{tail}}4
\quad\text{or}\quad
W_{sym}^+\ge\frac{\nu D_{tail}}4,
}
\]

with exact ties joint.

The first alternative is same-event donor tracing.  The second is existing strain/deformation provenance.  No shell-dependent factor `M/N` is introduced, so the high-tail causal law is unchanged.

## 7. What this removes from the master

Before this theorem, `positive_resolved_cross_interface` could appear as a named recursive owner whose physical continuation was only described qualitatively.  After the quotient it cannot form an independent recursive generation chain:

- symmetric interface is a witness relay to existing strain/deformation ownership;
- skew interface is a same-time conservative donor relay;
- skew circulation is internal divergence-free flux and cannot accumulate recursion depth;
- neither branch creates a new entropy, reset, clock, or scale coordinate.

This closes the **interface-owner completion** seam without introducing another master currency.

## 8. Scope

This theorem does not prove that the donor role reached by the skew quotient must itself terminate.  It also does not prove that every repeated strain/service/reuse owner terminates, nor does it close the UV-unbounded full-survivor alternative.

Its contribution is narrower and structural: resolved interface work is no longer an independent escape mechanism.  The remaining global problem is pushed back to the genuine physical owners already present in Navier--Stokes.

There is no claim here of a proof of global regularity for 3D Navier--Stokes.
