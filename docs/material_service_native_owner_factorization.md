# Native material-service owner factorization

Status: **DRAFT THEOREM CANDIDATE — independent PR, not certified, not merge-ready.**

This note attacks the current material/source-service part of mixed genuine-owner recurrence by removing one remaining ambiguity in the causal ontology. It does **not** invent a new material cost. It asks a more basic question:

> once Navier--Stokes has already produced a positive physical service law and that law already has its native physical owner, can OO/ON/NN material classification create another independent recursive owner or another event vertex?

The answer is no. Material classification is positive restriction of an already-existing law. The native owner passes through the restriction exactly, as provenance of the same law and same event.

## 1. Positive service exists before material ownership

Let the already-certified PDE route produce a positive locally finite service measure

\[
 d\sigma \ge 0.
\]

Assume its physical owner has already been resolved as a positive owner decomposition

\[
 d\sigma=\sum_{r\in\mathcal R} d\sigma_r,
 \qquad d\sigma_r\ge0.
\]

The index `r` is a **pre-material physical root**: actual positive HH work after energy reentry, an objective source/SGS service law, an existing strain/deformation law, a heat/dissipation service law, or another independently registered native PDE owner. It is not allowed to be `MATERIAL_RELINK`, `NEW_COHERENT_ANCESTRY`, membership rereading, selected-family `R_switch`, or smooth conservative `K_phys` relink. Using one of those labels as the input owner would be circular, because those labels are read only after the physical law already exists.

For the two intrinsic coherent endpoints `zeta_0,zeta_1` and transported old set `O`, put

\[
\chi_{OO}=1_O(\zeta_0)1_O(\zeta_1),
\]

\[
\chi_{ON}=1_O(\zeta_0)1_{O^c}(\zeta_1)
          +1_{O^c}(\zeta_0)1_O(\zeta_1),
\]

\[
\chi_{NN}=1_{O^c}(\zeta_0)1_{O^c}(\zeta_1).
\]

Pointwise,

\[
\chi_{OO}+\chi_{ON}+\chi_{NN}=1.
\]

No field is decomposed and no new positive part is taken.

## 2. Material restriction commutes with native-owner disintegration

For `C in {OO,ON,NN}`, define

\[
 d\sigma_C:=\chi_C\,d\sigma,
 \qquad
 d\sigma_{r,C}:=\chi_C\,d\sigma_r.
\]

Because multiplication by the same nonnegative measurable indicator distributes over a positive-measure sum,

\[
\boxed{
 d\sigma_C
 =\chi_C\sum_r d\sigma_r
 =\sum_r\chi_C d\sigma_r
 =\sum_r d\sigma_{r,C}.
}
\]

Consequently

\[
\boxed{d\sigma=d\sigma_{OO}+d\sigma_{ON}+d\sigma_{NN}}
\]

and, owner by owner,

\[
\boxed{d\sigma_r=d\sigma_{r,OO}+d\sigma_{r,ON}+d\sigma_{r,NN}}.
\]

This is an identity, not a domination estimate. It has three immediate causal consequences.

1. `ON` is genuine material-interface **provenance**, but its service mass is still the restriction of the physical law which already owned it. A second `MATERIAL_RELINK` charge would clone the same service.
2. `NN` is genuine fresh/new-material **provenance**, but its service mass is still the restriction of the physical law which already owned it. A second `NEW_COHERENT_ANCESTRY` charge would clone the same service.
3. There is no later Hahn split. The input law is already positive; material classification only restricts it.

Thus the correct direction is

\[
\boxed{
\text{native physical owner/event}
\longrightarrow
\text{OO/ON/NN provenance on that same law},
}
\]

not a second causal vertex.

## 3. A material-boundary contact does not by itself modify the carrier equation

The same-carrier theorem already gives, while `Q` and the registered probe `psi` remain fixed,

\[
z(t)=z(s)+I_{HH}[s,t]+I_{interface}[s,t].
\]

No material characteristic function occurs in this identity. Therefore an intrinsic endpoint reaching `partial O` can invalidate an `NN` witness or change the OO/ON/NN sidecar, but **cannot by itself create a new coefficient impulse**. The carrier continues unless one of its actual PDE monitors fires.

This is the dynamic counterpart of the measure identity above: material membership is read from the physical law; it is not an extra term in the PDE.

Selected-family switching is also already separated. Its exact

\[
R_{switch}=\sum_{C\in S_{old}\triangle S_{new}}E_C
\]

remains a Moyal boundary currency, but the same-state anti-theorem proves that `R_switch>0` can coexist with zero state increment and zero positive/negative physical work. It therefore stays a non-event sidecar and cannot be reused here as material work.

## 4. Genuine smooth role change has no independent material remainder

If the smooth role really changes, one must return to the already-certified native interface work, not to a material label. After common observer transport is quotiented, the residual smooth interface has the physical operator split

\[
K_{phys}+S.
\]

The `K_phys` pair matrix is antisymmetric and has exact finite donor closure at the same physical event. Positive `K_phys` relink is therefore conservative same-event provenance with zero recursive depth. The symmetric part is the already-existing strain/deformation owner.

Hence an actual role change has only the already-registered continuations:

\[
\boxed{
\text{same-event }K_{phys}\text{ donor relay}
\quad\text{or}\quad
\text{independently registered native strain/work/source owner}.
}
\]

There is no third primitive `MATERIAL_RELINK` energy source left over. If neither a bound `K_phys` law nor an independently registered native work/source/strain owner is supplied, the candidate theorem fails closed.

## 5. Mixed-owner consequence: material/new ancestry is not an independent letter

After this factorization, a recursive path should be read only after deleting zero-depth material operations:

- membership rereading;
- selected-family Moyal boundary sidecars;
- OO/ON/NN reclassification of already-positive service;
- smooth conservative `K_phys` donor relays.

Fresh `NN` service and interface `ON` service retain their material provenance, but they do not create a second event. The pre-material native owner remains attached to the physical law that already existed. Therefore material/new-ancestry vocabulary cannot by itself break a consecutive high-strain epoch or a consecutive signed-good generated-HH epoch. It breaks such an epoch only when an **independent native owner event** also occurs.

This strictly shrinks the mixed genuine-owner alphabet. The remaining event-facing recurrence is pushed back to the local PDE roots which actually do work, dissipate, deform, or supply source/service. That is exactly the desired direction: fewer primitives, obtained by descending into the PDE rather than by imposing a master convention.

## 6. Finite-atom executable model

`src/material_service_native_owner_factorization.py` implements the exact finite-atom algebra used for regression tests. Every atom carries a nonnegative service weight, a pre-material native owner, and its two old/new endpoint flags.

The implementation verifies both identities

\[
\sigma=\sigma_{OO}+\sigma_{ON}+\sigma_{NN}
\]

and

\[
\sigma_r=\sigma_{r,OO}+\sigma_{r,ON}+\sigma_{r,NN}
\]

for every native owner. It rejects circular material/new-ancestry roots, negative service weights, a later Hahn interpretation, new recursive event creation, and unbacked genuine role changes. A typed smooth-relink donor certificate is accepted only as a same-event relay.

The finite-atom code is not the proof of the Radon statement; it is an exact executable model of the pointwise restriction/distributivity identity and a fail-closed integration guard.

## 7. Scope

If the theorem survives full review and certification, it closes **material/new-ancestry as an independent causal owner class** inside mixed genuine-owner recurrence. It also removes pure material boundary contact as a reason to kill an otherwise unchanged smooth carrier.

It does **not** yet prove termination of the remaining cross-family recurrence among actual HH/high-tail work, strain/dissipation, objective source/SGS/pressure service, and causal reuse. It does not close the degenerate full-signed Young/Christ HH seam, the initial-data interface, or the hypothetical singular-time conclusion.

There is no Navier--Stokes global-regularity claim.

The next theorem after this draft should work on the **compressed native-owner word**, not the historical material labels: prove a cross-owner telescope/transition law using only the local physical ledgers supplied by the surviving HH, strain/dissipation, and objective source roots.
