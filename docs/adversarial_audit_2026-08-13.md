# Adversarial audit — unresolved findings

Audited against main `5e719a4`.  Keep each finding in this Draft PR until the corresponding theorem is repaired and independently re-audited; findings are cumulative rather than replacement notes.

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
