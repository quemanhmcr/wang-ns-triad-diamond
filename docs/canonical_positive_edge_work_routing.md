# Canonical positive edge-work routing and signed hard-cell compression

Status: **CERTIFIED** — exact implementation/performance SHA `d78a3f4223c16f541548e7628cd73f70e6fdee6c`.

Certification evidence:

- dedicated routing run `31559241272`: **success**, `50,000` stress states;
- independent audit `31559242700`: **success**, `75,000` adversarial states plus an actual dealiased Fourier--Galerkin Navier--Stokes routing audit at resolution/cutoff `24/7`, `64` RK4 steps and `5` physical snapshots;
- full physical-energy causal integration `31559248529`: **success** on the same exact SHA;
- measure-layer compatibility: continuum dedicated `31559244152`, continuum adversarial/actual-NS/multi-child audit `31559245648`, and helical physical-edge run `31559247024`, all **success** on the same SHA.

The dedicated stress had worst `G+B` mass residual `2.2197019035664715e-16`, hard positive pushforward residual `0`, `50,000` retained positive nonforward bad cases, no stage-zero time failures and no geometry-good marking promotions.  The actual-NS audit had signed reconstruction, positive-mass reconstruction and hard pushforward residuals all `0`.  Artifacts are recorded under the matching `recorded-results/<run-id>/` directories.

The physical starting point is already fixed by the continuum helical edge theorem.  On the unordered Fourier/helicity edge space,

\[
 dW=C_F T_e\,d\Lambda_{edge},\qquad
 dA=C_F A_e\,d\Lambda_{edge},\qquad
 dF=J_ec_e\,dA,
\]

with `dW` signed physical child-energy work and `dA` only a positive capacity reference.  The canonical causal law is the Hahn-positive part

\[
 \mu^+=dW^+.
\]

This theorem never replaces it by capacity probability and never takes a second causal Hahn decomposition after an analyst chooses cells.

## 1. Exact physical fate partition

Put

\[
 r_e=(J_e/J_*)c_e,\qquad \eta_0=10^{-4},
\]

and restrict the already-existing positive law by

\[
 G=\{T_e>0:r_e>1-\eta_0\},\qquad
 B=\{T_e>0:r_e\le 1-\eta_0\}.
\]

Because these are complementary measurable restrictions of the same `dW+`,

\[
 \mu^+=\mu_G^++\mu_B^+
\]

exactly.  The source implementation obtains both restrictions by replaying the registered physical fibers; no caller may supply the good/bad masses.

Positive nonforward work has `J_e=0`, hence `r_e=0`.  It remains positive physical child work and therefore belongs to `B`; it is not discarded and receives no invented scale progress.

## 2. Capacity is used only after the bad causal restriction exists

On `B`, the native edge identity gives

\[
 F(B)=\int_B r_eJ_*\,dA\le (1-\eta_0)J_*A(B),
\]

so every nonzero bad positive-work sublaw has

\[
 \epsilon_B=1-\frac{F(B)}{J_*A(B)}\ge \eta_0.
\]

This is a statement about the *same support already selected by `dW+`*.  It is not a capacity-majority argument.

For every certified `0<tau<=0.1`, the existing physical block interface uses

\[
 \delta_\tau=\frac{\tau^2}{1\,036\,800\,000}<10^{-4}.
\]

The source first passes the bad restriction's own `epsilon_B` through the exact transfer-deficit channel factored from `coherent_service_or_flat_gate`.  That typed `FixedTransferLossGate` is the same threshold/cause interface used by the whole physical block gate; only after it returns `triggered=True, cause=physical_transfer_cost` may the source construct the compiler witness.  The computed bad `dW+` mass itself is then passed unchanged to `compile_transfer_measure` and to the joint stop projection.  Thus

\[
 \mu_B^+\to \texttt{fixed\_transfer\_loss}
 \to \texttt{TRANSFER\_WORK\_LOSS}
 \to \texttt{TRANSFER\_COST}.
\]

This is stage-zero block fate.  It has `first_time=None`.  The legacy fine compiler sentinel `0.0` is removed; physical `t=0` remains reserved for the absorbing initial boundary.

## 3. Deterministic hard cells inherit cause and separately compress signed work

The event-anchored role theorem already supplies deterministic disjoint Borel frequency cells together with pointwise orthogonal helical projectors.  Such a hard role assignment is a genuine measurable map `pi` on the physical mode/edge space.  Therefore it has three different, simultaneously valid linear pushforwards:

\[
 P_C=(\pi_#\mu^+)(C),\qquad
 N_C=(\pi_#\mu^-)(C),\qquad
 T_C=(\pi_#\mu)(C)=P_C-N_C.
\]

`P_C` is inherited causal mass. `N_C` is the directly inherited canonical Hahn-negative mass, retained for cancellation provenance; it is not reconstructed by a later hard-cell Hahn split. `T_C` is signed hard-cell trilinear work.  The implementation certifies

\[
 [T_C]_+\le P_C
\]

cellwise and after summation.  The gap is exactly cancellation hidden by the chosen aggregation.  Coarsening or refining the hard labels preserves the total inherited `pi_#mu+`; it may change `sum_C[T_C]_+`.  The latter is therefore retained only as a diagnostic and is never declared a new master causal law.

This is the signed-cell compression needed by Young/Christ.  Young sees `T_C`, not the gross positive mass `P_C`.

## 4. Geometry-good means only Young-eligible; signed-cell binding needs fate purity

The good causal restriction is exported as a `YoungEligiblePositiveWork` handle.  Its fields enforce

- `marking_good=False`;
- `young_certified=False`;
- `registered_generated_continuation=False`.

There is a further physical obstruction that must not be hidden by refinement.  A hard cell can contain good positive causal work, bad positive causal work already sent to transfer loss, and negative signed work.  Young must see the **full signed** trilinear cell work in order to preserve cancellation.  But then a mixed cell would let already-terminal bad positive work help the good branch.  Subtracting only the bad positive part is not allowed either: that nonlinear subtraction is no longer the signed trilinear work of the hard cell.

Therefore the theorem binds a good causal sublaw to signed-cell Young only on a **fate-pure hard cell**

\[
 (\pi_#\mu_G^+)(C)>0,\qquad (\pi_#\mu_B^+)(C)=0.
\]

On such a cell, downstream Young/Christ receives the full

\[
 T_C=(\pi_#dW)(C),
\]

so negative work still cancels before saturation, while no already-terminal bad positive work can assist.  If a hard cell contains both good and bad positive causal mass, the good mass remains exactly accounted as **unresolved Young-eligible mixed-fate work**.  The theorem does not invent a finer cell to remove that obstruction.

Thus geometry-good, fate-pure signed-cell binding, Young-good and generated continuation are distinct physical layers.

## 5. Why the general coherent POVM is not used as the causal adapter

`recursive_coherent_witness_extraction` remains an exact signed identity:

\[
 W_{CDE}=2\Re\langle A_Ew_3,\mathcal N(A_Cw_1,A_Dw_2)\rangle,
 \qquad \sum_{CDE}W_{CDE}=W_{HH}.
\]

But its general positive localization operators are not, in general, a deterministic map on Fourier/helical edges.  Cross-term interference means that Hahn-splitting those signed coherent atoms is not automatically a positive mass-preserving pushforward of canonical `dW+`.  The source and documentation therefore downgrade that positive coherent law to representation-level diagnostics unless a separate positive kernel/disintegration theorem is supplied.

Nothing is lost: the deterministic hard Fourier/helicity event roles already provide the exact master-facing signed compression required for the next Young test.

## 6. Adversarial scope

Dedicated tests/audit are required to cover:

- a genuine helical example with good capacity majority but bad actual-work majority;
- positive nonforward work routed to `B`;
- exact `mu_G^+ + mu_B^+ = mu^+`;
- the native bad deficit, equality with the whole-block transfer gate's threshold/cause channel, and every certified `tau<=0.1` fixed-transfer implication;
- parent swap, helical gauge representation and unit scaling invariance;
- replay rejection of forged ledger summaries;
- hard coarsening/refinement preserving inherited `dW+` while changing visible signed cancellation, and a deliberate coarsening that creates a mixed good/bad cell which remains unresolved rather than being promoted to Young;
- `fixed_transfer_loss.first_time is None`;
- geometry-good never promoting itself to Young-good or `marking_good`;
- an actual dealiased Fourier--Galerkin Navier--Stokes trajectory whose replayed edge law passes the same routing and hard pushforward.

This theorem does not certify Young saturation of the good branch.  In particular, mixed-fate hard cells remain an explicit signed-cell handoff seam rather than being repaired by arbitrary refinement.  It also does not identify a general coherent POVM Hahn law with canonical causality, does not close generic HH recurrence, and makes no claim of global regularity for 3D Navier--Stokes.
