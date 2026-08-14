# Full-signed Young state-mark factorization

Status: **DRAFT THEOREM CANDIDATE — PR #7, not certified, not merge-ready.**

This note removes a type error from the local HH/Young interface. The canonical positive causal law and the Young/Christ near-extremal theorem answer different questions:

- `dW+` says which positive child-energy work actually occurred;
- Young/Christ says what the three hard physical roles at that same event look like when their full signed trilinear interaction is near sharp equality.

A Christ Gaussian is therefore a **state/role witness**, not causal mass. Once that distinction is kept exact, there is no reason to subtract terminal bad positive work before deciding whether the shared role state is near extremal. The reservation quantity `T_C-b_C` remains a valid stronger robustness certificate from the predecessor theorem, but it is not required to attach a state mark to the surviving canonical good sublaw.

## 1. One physical event carries two exact ledgers

Fix one **event-anchored** deterministic hard Fourier/helicity product cell `C` at a fixed physical interaction time, exactly as in `event_anchored_role_registration.md`. The theorem is applied on that event-state fiber before any aggregation across distinct physical times. A later event may reuse the same hard role labels, but it is a different state and must obtain its own mark. The canonical edge law has already been reconstructed signed and Hahn-split once on edge space. Its hard pushforwards are

\[
g_C=(\pi_\#dW_G^+)(C),\qquad
b_C=(\pi_\#dW_B^+)(C),\qquad
n_C=(\pi_\#dW^-)(C),
\]

and the unchanged signed hard-cell work is

\[
\boxed{T_C=g_C+b_C-n_C.}
\]

These are causal/work data. Geometry-bad positive work `b_C` is already on the existing `TRANSFER_WORK_LOSS` route; geometry-good positive work `g_C` is the only positive causal sublaw eligible to continue.

Independently, the same cell determines the actual hard role triple and its sharp full-cell Young upper

\[
Y_C=m_*A_3\prod_{j=1}^3\|f_{C,j}\|_{3/2}.
\]

Let `xi_C` be the already-existing normalized symbol-freezing discrepancy. On the generated-facing positive branch, the ordinary complex Young/Christ premise is

\[
\delta_C^{state}:=1-\frac{T_C}{Y_C},\qquad
\boxed{\delta_C^{state}+\xi_C\le\delta_{Christ}(\varepsilon_G).}
\]

When `T_C>0`, the existing complex-Young reduction and Christ's theorem mark the actual parent roles by a complex Gaussian near-extremizer. Nothing in this implication is a positive causal measure. The theorem reads the physical state present at the event.

## 2. Intrinsic non-creation law: causal restriction does not restrict state

Let `M_C` denote deterministic/measurable data extracted from the physical role state of this same event-state cell `C` after the full-signed Christ gate: Gaussian parameters, a dual analysis probe, a covariance representative, or another downstream state-only mark.

The surviving good causal law on this cell is still

\[
\mu_{G,C}:=(\pi_\#dW_G^+)|_C.
\]

Adjoin the state mark only as metadata:

\[
\boxed{\widetilde\mu_{G,C}=(\operatorname{id},M_C)_\#\mu_{G,C}.}
\]

If `pr` forgets the mark coordinate, then identically

\[
\boxed{\operatorname{pr}_\#\widetilde\mu_{G,C}=\mu_{G,C}.}
\]



### 2.1 Exact NS role-state invariance under edge-fate restriction

At the fixed physical event time, the hard roles are formed directly from the actual velocity field before causal fate routing:

\[
f_a=P_a u(t),\qquad f_b=P_b u(t),\qquad f_c=P_c u(t).
\]

The canonical good/bad split acts on the **edge work measure** only:

\[
dW_G^+=\mathbf 1_G\,dW^+,\qquad dW_B^+=\mathbf 1_B\,dW^+.
\]

It does not act on `u(t)` and does not define counterfactual fields `u_G`, `u_B`, `P_a u_G`, or `P_a u_B`. Therefore the physical role triple is exactly invariant under the fate restriction:

\[
\boxed{(P_a u(t),P_b u(t),P_c u(t))\ \text{is unchanged by}\ dW^+\mapsto \mathbf1_G dW^+.}
\]

This is the NS-specific core of the theorem. Any state functional or theorem witness evaluated from those same roles is unchanged. What changes is only which already-existing positive work atoms continue in the recursion.

A construction that first synthesizes a new “good-only velocity field” from edge fates would be a different nonlinear representation and is explicitly outside this theorem.

The mark adds no mass, owner, event, recursion depth, scale progress or second Hahn law. It may not be carried to another time merely because the same Fourier/helicity role label reappears. It records a theorem about the shared physical state on already-existing good causal atoms.

The bad sublaw is not reopened:

\[
\mu_{B,C}:=(\pi_\#dW_B^+)|_C\longrightarrow \texttt{TRANSFER_WORK_LOSS}
\]

exactly as before. Canonical negative work remains inside `T_C` as signed physical evidence and keeps its certified same-time cyclic donor provenance; it is neither removed nor used as payment.

Thus there is no statement that bad positive work or negative work “pays” for the good branch. There is no payment: Young/Christ supplied a property of the common roles, not a causal currency.

## 3. Why terminal bad work may coexist with the state witness

The predecessor reserved theorem used

\[
1-\frac{T_C-b_C}{Y_C}=\left(1-\frac{T_C}{Y_C}\right)+\frac{b_C}{Y_C}
\]

to certify that the cell would remain near Young extremal after every unit of terminal bad-positive work were counterfactually withheld. That is a valid and stronger robustness statement.

But the real Navier--Stokes event does not counterfactually delete that interaction before observing its velocity field. All three hard roles are present simultaneously. If their actual full signed trilinear form satisfies Christ, then the corresponding Gaussian/dual mark is a true statement about those actual roles.

Routing one positive submeasure to a terminal fate cannot retroactively change which velocity field existed at the same event. Therefore the additional counterfactual robustness test is unnecessary for **state marking**.

Here “terminal” has exactly the project's existing meaning: terminal in the forward-transfer recursion, **not removed from Navier--Stokes**. Geometry-bad and positive-nonforward work remains real modal energy redistribution in the PDE. Therefore stage-zero routing cannot be interpreted as editing the simultaneous physical interaction before Young reads the role state.


## 4. The Christ-margin seam disappears rather than being estimated

The spare quantity `m_C = delta_Christ - (delta_full + xi_C)` was needed only because the reservation theorem asked whether `b_C/Y_C` fit inside unused Christ margin. The state-mark factorization never asks that question.

For a generated-facing cell there are only two cases:

- `delta_state + xi_C <= delta_Christ`: the state mark is available on the unchanged good `dW+` restriction;
- `delta_state + xi_C > delta_Christ`: there is no state mark, and the existing Young/symbol/transfer failure route applies.

There is no third physical owner at `m_C=0`, and no dynamical quantity needs to be invented to control approach to that analyst-defined boundary. Equality is simply included in the Christ gate.

Consequently the degenerate full-signed Christ-margin seam is a theorem-layer artifact once the Christ output is typed correctly as state information.


## 5. Downstream compatibility

The current hard-role spine already uses the mark in this state-only manner.

1. **Complex Young -> dual Gaussian.** Christ marks the actual complex parent role. The dual-Gaussian theorem uses the mark only to choose an analysis probe and proves a coefficient of that actual role; it does not create work mass.
2. **Hard event -> smooth PDE carrier.** For the event projector `P` and smooth envelope `Q` with `QP=P`, the exact registration `<Pu,phi>=<Qu,Pphi>` attaches that coefficient to the same physical role without replacing a transfer measure.
3. **Bargmann/Moyal parent identity.** The Gaussian coefficient chooses a nearby actual energy anchor. Pushing the surviving positive parent-slot law to that anchor preserves total causal mass; collisions merge as reuse. The nonlinear work atom is explicitly not required to occupy the same coherent cell as the Christ mark.
   This is also the no-double-counting mechanism for shared state: if several surviving causal slots see the same physical parent anchor, they coalesce onto that anchor and are recorded as reuse rather than minting several copies of the parent energy.
4. **Physical pair productivity.** The logarithmic parent-product estimate is integrated under actual retained positive child-work. Parent amplitudes and registered coefficients are properties of the hard roles; adding a state mark does not reweight `dW+`.

Thus no downstream physical currency changes.

## 6. Scope

This theorem candidate closes the **local hard-cell Young handoff** if the type separation above survives review:

- full signed `T_C` decides the role-state Young/Christ mark;
- canonical good `dW+` alone supplies continuing causal mass;
- canonical bad `dW+` remains terminal;
- canonical `dW-` remains signed donor provenance;
- no reservation margin, fresh Hahn law, coherent-POVM positive kernel or failure-payment map is introduced.

It does **not** identify fresh coherent Hahn atoms with canonical causality. A general coherent POVM still needs a separately proved positive kernel if one insists on using its positive atoms as a master causal law. The current hard-event -> dual/Bargmann route does not require that identification.

It also does not prove generic HH recurrence finite, close mixed genuine-owner recurrence, close the initial-data/singular-time interfaces, or prove 3D Navier-Stokes global regularity.

The structural law is

`Causal restriction changes the work measure, not the shared event state.`

It follows by keeping the exact NS state and the exact NS work law as two different physical objects instead of forcing one to serve as the other.
