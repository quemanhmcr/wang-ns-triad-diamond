# Adversarial audit — unresolved findings

Audited against main `5e719a4`. Re-audited through latest main `24a7257`; findings 1--2 remain unresolved. Findings 3--4 audit Draft PR #6 at `efbdb06`. Keep each finding in this Draft PR until the corresponding theorem/proof-support bridge is repaired and independently re-audited; findings are cumulative rather than replacement notes.

## 1. Same-carrier stock-only projection can erase simultaneous positive HH work

The inherited-stock theorem correctly added a guard for simultaneous classified residual work, but the typed certificate/master projection still carries no analogous value for actual positive HH child-energy work `W_HH^+`.  Its no-first-stop hypothesis only controls the cumulative coefficient obstruction `|I_HH|`, and the repository already distinguishes that locator from physical positive HH work.

This is not a Navier--Stokes counterexample; it is a countermodel to the proof-support implication needed by the stock-only projection.  In the same selected-role energy algebra, take on `[0,1]`

`c(t)=1+eps sin(2 pi m t)`, `G=R=0`, `F_HH=c_dot`,

with `eps=0.1`, `m=2`.  Then `E0=E1=1`, `K=0`, `W_R^+=0`, while

`I_HH(t)=eps sin(2 pi m t)`, so `sup_t |I_HH(t)|=0.1<1/2`

and the HH coefficient face is never hit.  Nevertheless the actual positive HH child-energy work is

`W_HH^+ = int_0^1 2[c c_dot]_+ dt = 4 m eps = 0.8`,

which is larger than the repository's clean generated-work scale `8E1/15`.

Thus coefficient no-hit does not imply absence of a simultaneous positive HH-work owner.  The sidecar-free central route currently projects a typed inheritance certificate to `owners=()` / `same_carrier_inherited_energy_stock_relay`, while binding only `E0`, the inheritance threshold, and classified residual positive work.  A large simultaneous `W_HH^+` is therefore not preserved by that projection.

Minimal repair: bind the actual same-`Q^2` positive HH work law into the inherited-stock certificate/master projection, just as residual positive work is bound now.  If that work realizes an existing HH owner face, keep the inherited stock component zero-depth but preserve the HH owner event-facing (or fail closed from stock-only projection).  Do not substitute `|I_HH|` for `W_HH^+`, perform temporal deposit matching, or re-Hahn a downstream representation.

## 2. Selected-family service cover loses exact ties by insertion order

The new material-sidecar theorem correctly separates selected-family Moyal boundary energy from Navier--Stokes work and inherited stock.  However its `selected_family_service_no_escape_binding` delegates to `coherent_transfer_cells.service_no_escape`, which collapses the three-way no-escape cover to a single `max(...)` branch.

The exact algebra is only

`P_plus <= E_final + P_minus + R_switch`,

so each candidate is nonnegative and at least one is `>=P_plus/3`.  Exact ties must remain joint; they are not physically ordered.  A minimal exact fixture is

`P_plus=3`, `E_final=P_minus=R_switch=1`.

The balance is saturated and all three candidates lie exactly on the `P_plus/3=1` face.  Current code constructs candidates in the order

`terminal_coherent_energy`, `backflow_or_cancellation`, `relink_symmetric_difference`

and Python `max` returns the first maximum.  Thus the same physical/Moyal data are reported only as `terminal_coherent_energy`; changing dictionary insertion order changes the reported branch.  This conflicts with the repository's established joint-tie semantics and would erase simultaneous provenance if this helper is promoted into owner routing.

This does not invalidate the Moyal inequality or the zero-generation-depth interpretation of `R_switch`; the issue is the proof-support/service projection of the one-third cover.

Minimal repair: return the complete joint set of realized owner candidates (at least every exact maximum tie, preferably every candidate on the certified `>=P_plus/3` face) and let downstream joint-stop/master logic quotient fates without theorem-name or insertion-order priority.  No tie-weight normalization is needed.


## 3. Material-only epoch quotient does not certify the hypotheses of the signed-good HH telescope

Draft PR #6 correctly argues that a material-only observation need not create a second causal event.  The new `native_owner_epoch_quotient`, however, makes the stronger claim that after deleting records with `no_causal_stop=True`, the existing signed-good generated-HH physical-time telescope applies to the resulting event subsequence.

That implication needs more than event kind and event order.  The certified HH telescope explicitly requires, for every consecutive pair,

`N_child,next = N_parent,prev`,

and physical support nesting

`H_next subset [s_prev,b_prev]`,

with one common scaled lifetime parameter.  Its implementation rejects a pair when either the carrier-scale continuation or the support inclusion fails.

The new quotient stores only `(time, kind, projection, witness)` and forms maximal runs solely from the event-kind sequence after material records are deleted.  Its own regression test accepts

`SIGNED_GOOD_HH -> material-only -> SIGNED_GOOD_HH`

without supplying any parent/child frequencies, physical work supports, or a certificate that the same recursive carrier lineage survived the deleted observation.  Thus order-preserving deletion proves that material bookkeeping created no event vertex, but it does not by itself prove the cross-step hypotheses required by the HH telescope.

This is a proof-support/composition gap, not a counterexample to the underlying signed-good HH theorem and not a claim that material rereading is physical work.

Minimal repair: either weaken the quotient conclusion to a purely topological statement, or carry typed `SignedGoodGeneratedHHStep` data plus a zero-depth state-preservation relay and re-check the existing cross-step conditions (`child_next=parent_prev`, common `c`, and support nesting) after quotienting.  A material observation may be erased from event depth without being allowed to manufacture theorem eligibility.

## 4. Descending fresh-SGS certificate can choose an unbound hard-shell candidate

Draft PR #6's new descending fresh-SGS theorem has a sound conditional core: if the **actual selected physical child** satisfies `N_next<=N/2`, the local fresh-service envelope gives a positive parent-frequency floor and repeated half-scale descent is finite.

The executable bridge does not yet bind that actual selection.  The certified upstream fresh-scale theorem proves only

`max(mu_M,mu_2M) >= (2/3) M ||u_j||_2^2`,

so one of the two physical hard shells `(M,2M)` carries the lower bound.  `fresh_service_scale_route` returns the two candidate frequencies and one scalar `hard_shell_mass_lower`, but it does not return which shell actually realizes the maximum (nor the two actual shell masses).

`DescendingFreshSGSRenewalStep` then accepts a caller-supplied `child_frequency` whenever it equals either candidate, and independently accepts a caller-supplied `child_critical_mass` whenever that scalar is above the common lower.  These two facts are not bound to the same physical shell.  This is load-bearing at `j=-1`: the candidates are `(N/2,N)`, so a caller can present the descending `N/2` candidate with the existential lower even when the actual heavy shell is the non-descending `N` candidate.  Exact ties are likewise not represented jointly.

This does not invalidate the theorem restricted to a genuinely certified descending child, and for `j<=-2` both candidates are already `<=N/2`, so the ambiguity is harmless there.  The issue is the typed proof-support bridge used to certify that the selected child is the physical one.

Minimal repair: safest is to restrict the unconditional descending route to `j<=-2`.  To retain the lower `j=-1` case, carry the two actual hard-shell masses (or a typed selected-shell certificate derived from them), select the actual maximum with exact ties kept joint, and enter the descending telescope only on the genuinely realized `N/2` face.  Do not lexicographically choose the lower scale and do not attach the existential max lower to an arbitrary candidate.
