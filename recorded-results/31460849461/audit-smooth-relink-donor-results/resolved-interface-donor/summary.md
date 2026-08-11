# Resolved interface donor/circulation quotient

Status: **EXACT_RESOLVED_INTERFACE_DONOR_QUOTIENT__POSITIVE_INTERFACE_TO_CONSERVATIVE_SKEW_DONOR_OR_EXISTING_STRAIN__FINITE_SAME_EVENT_DONOR_EXHAUSTION__CIRCULATION_ZERO_RECURSION_DEPTH**.

The resolved low--high cross/interface work is not an independent energy source.  Write its actual physical linearized operator as `L=K+S`, with `K*=-K` and `S*=S`.  Event by event,

`R = R_K + R_S`,

so positive work obeys

`[R]_+ <= [R_K]_+ + [R_S]_+`.

Therefore positive resolved cross/interface work has only two native continuations: conservative skew role redistribution, or the already existing symmetric strain/deformation work.  One carries at least half of the positive interface law; exact ties remain joint.  The symmetric branch is not charged again: it delegates to the existing coherent-deformation / high-strain / objective-source / critical-`D_V` owners.

For a complete orthogonal event-role partition, the skew pair work is

`T_ab=2 Re <w_a,K w_b>=-T_ba`,

and the net skew work of role `a` is `R_a^K=sum_b T_ab`.  With directed physical flux `F[b->a]=[T_ab]_+`,

`R_a^K = incoming_a - outgoing_a`,

and `sum_a R_a^K=0`.  More generally, for every role set `C`, all internal transfers cancel and

`sum_(a in C) R_a^K = boundary inflow(C) - boundary outflow(C)`.

This immediately gives a finite donor theorem.  Start from any role with positive skew gain and close the set backward under every positive incoming donor edge.  If every role in that closure had nonnegative net skew work, the closure would have strictly positive total gain.  But by construction no positive flux enters it from outside, so its exact boundary balance is nonpositive.  Contradiction.  Hence the closure contains an actual negative-net donor role.  Removing cycles gives a donor path of at most `#roles-1` edges.

The whole donor trace occurs at the **same physical event time**.  A role cycle is circulation, not a new Navier--Stokes generation.  It creates no extra causal charge and no master recursion depth.  No donor is selected by role index; all reachable negative-net donors remain a set-valued physical provenance mark.

For the high-tail interface owner this yields, in the unchanged common `N dW` unit,

`W_interface^+ >= nu D_tail/2`

implies

`W_skew^+ >= nu D_tail/4`

or

`W_symmetric^+ >= nu D_tail/4`,

with ties joint.  The first is same-event donor tracing; the second is existing strain/deformation provenance.  No `M/N` reweighting is introduced.

Stress: `100000` split/flux/donor/high-tail states
- worst signed `R=R_K+R_S` residual: `8.882e-16`
- minimum positive-cover margin: `-1.776e-15`
- worst role-divergence residual: `2.138e-15`
- worst total skew-work residual: `3.284e-15`
- worst donor-closure balance residual: `3.284e-15`
- minimum recipient incoming-flux margin: `0.000e+00`
- donor-existence failures: `0`
- maximum sampled shortest donor path: `3`
- high-tail component failures: `0`

This theorem closes **resolved interface work as an independent recursive-generation loophole**.  It does not claim that the eventual donor or strain owner globally terminates, and it makes no claim of Navier--Stokes global regularity.
