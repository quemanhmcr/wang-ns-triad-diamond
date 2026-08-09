# Material-label sidecar quotient for smooth carrier continuation

## The ambiguity to remove

Several older interfaces used one boolean word, `material_relink`, for physically different operations:

1. an intrinsic material endpoint changes old/new membership;
2. the selected coherent family used by the service/ancestry ledger changes;
3. the smooth PDE role \(Q\) or the registered analysis probe \(\psi\) actually changes.

Only the third item necessarily changes the coefficient object being propagated.  Treating all three as the same carrier stop forces an unnecessary reconstruction and double-counts some nonaffine physics.

## Same carrier means no new label term in Duhamel

Fix the same smooth role \(Q\) and the same registered analysis probe \(\psi\) across a common-slice interval.  The exact coefficient identity is

\[
z(t)=z(s)+I_{HH}[s,t]+I_{interface}[s,t].
\]

There is no factor \(1_O(\zeta)\), no dyadic material address, and no selected coherent-family characteristic function in this equation.

A material address can move because the physical flow is nonaffine.  That does **not** mean the motion is free: the effect of nonaffine dynamics on the carrier coefficient is already in \(I_{interface}\).  What is forbidden is charging the same event twice by adding another independent `material_relink` impulse merely because an intrinsic address or old/new indicator crossed a bookkeeping boundary.

Thus, while \(Q\) and \(\psi\) remain the same, the carrier first-stop faces are still only

\[
|I_{interface}|\ge\frac14|z(t)|,
\qquad
|I_{HH}|\ge\frac12|z(t)|.
\]

If neither is hit,

\[
|z(s)|\ge\frac14|z(t)|
\]

regardless of a simultaneous material-label sidecar event.

## Old/new service ownership can be reread without creating service

After a positive coherent service law has been created, every edge already has two actual intrinsic endpoints.  Let the positive edge weights be \(s_e\).  For any old pool \(O\), exact ownership gives

\[
s_e=s_{OO,e}+s_{ON,e}+s_{NN,e}.
\]

If an endpoint later changes old/new membership, hold the physical edge weights fixed and reread the same partition.  The OO/ON/NN category masses may move, but

\[
S_{OO}^{before}+S_{ON}^{before}+S_{NN}^{before}
=
S_{OO}^{after}+S_{ON}^{after}+S_{NN}^{after}
=S.
\]

Material relabeling therefore creates no service and destroys no service.  It is a sidecar ownership event.

## Selected-family switching keeps its exact Moyal charge

Changing the selected coherent family is also not declared free.  If \(S_{old}\) and \(S_{new}\) are two cell sets with positive cell energies \(E_C\), then

\[
|E(S_{new})-E(S_{old})|
\le
R_{switch}
:=
\sum_{C\in S_{old}\triangle S_{new}}E_C.
\]

The full \(R_{switch}\) remains in the ancestry/service ledger and may force the existing relink branch of the master.  The quotient says only this:

> if \(Q\) and \(\psi\) did not change, paying or recording \(R_{switch}\) does not require discarding the smooth carrier \(Q u\) and constructing another one.

The master may recurse while reusing the same PDE carrier.

## Genuine role/probe change is not quotiented

If \(Q\) changes, or the analysis probe \(\psi\) changes in a way not already represented by the registered interface term, the same coefficient identity cannot simply be reused.  The quotient refuses such an input and delegates to the existing event-role, nonaffine-interface, Xi, or physical relink registration.

So the physically correct split is

\[
\boxed{
\text{material/service sidecar change}
\neq
\text{carrier-role change}.
}
\]

The first can coexist with continued carrier propagation; the second needs an actual interface registration.

## Relation to the older common-slice theorem

The older `registration_first_stop(..., material_relink=True)` theorem is still correct as a conservative superset: stopping early cannot invalidate its lower bounds.  The new theorem refines only the subtype

\[
\text{same }Q + \text{same }\psi + \text{material sidecar only}.
\]

It does not rewrite old certificates or declare all historical `MATERIAL_RELINK` events transparent.

## Master-facing semantics

A single physical time can therefore carry two different outputs:

- a **carrier decision**: continue, interface stop, or HH-regeneration stop;
- zero or more **material sidecars**: membership update and/or selected-family Moyal switch.

The full joint event set is retained without lexicographic priority.  Sidecars remain ancestry/relink currency.  The only quotient is that sidecar currency does not automatically kill a smooth carrier whose PDE role and probe have not changed.

This is the carrier-level statement needed by the shortest critical-shell architecture, where materiality is deliberately assigned only after actual own-scale service exists.

## Scope

This theorem does not make repeated relinks free, does not turn \(R_{switch}\) into work, does not erase ancestry entropy/cycle accounting, and does not cover genuine role changes.  Source/pressure routing and final continuum master assembly remain open.  No 3D Navier--Stokes global-regularity conclusion is asserted.
