# RESEARCH LEDGER — current theorem architecture

Status: current through the certified **hard-tail true-upward supply** theorem, exact theorem SHA `d064bc4d780d9c90d36e64c7b84c3b771b74c896`. This file is deliberately theorem-facing and compact. CI runs, stress counts, residuals, artifact digests, and failure lineage live in `docs/experimental_certification_log.md`; the pre-compaction full ledger is archived under `docs/history/`.

## 1. Non-negotiable physical rules

Start from actual Navier–Stokes/PDE quantities and preserve their type.

- Reconstruct signed nonlinear work `dW` **before Hahn**. The canonical cause is `dW+`; `dW-` may acquire donor provenance but never replaces the positive cause.
- After a canonical Hahn split, downstream labels inherit cause by restriction or positive pushforward. Do not re-Hahn a coarsened signed representation.
- Native capacity is a reference/error envelope, never causal probability, work currency, reset resource, or recurrence budget.
- Persistent energy stock belongs to physical Fourier–helical modes. Hard interaction cells are same-time labels, not wallets.
- Same-time donor provenance is distinct from between-time stock continuity. No FIFO/LIFO/proportional matching of earlier deposits to later withdrawals.
- Observer checkpoints, cover scales, coefficient magnitudes, Duhamel weights, entropy coordinates, and packet labels are not physical clocks unless an independent theorem makes them so.
- Conservative redistribution is real physics but does not automatically create recursive event depth.

## 2. Current physical spine

### Signed work and hard routing
For an unordered helical edge, direct Leray/Fourier work equals the Waleffe registration. The continuum edge law is signed and locally Radon; only afterward is `dW+` formed. Geometry/phase-bad positive work, including positive nonforward work, follows the existing stage-zero `TRANSFER_WORK_LOSS` route. Geometry-good work is merely Young-eligible.

Hard cells inherit three direct laws: signed `pi_#dW`, positive `pi_#dW+`, and negative `pi_#dW-`. In a mixed cell

`T_C = g_C + b_C - n_C`.

Terminal bad-positive assistance is reserved only in the scalar Young/Christ certificate; Young still sees the full signed `T_C`. On the certified low-deficit good core, actual work satisfies `dW/dA > 19/100`. Nondegenerate reservation failure is a domination statement involving existing `n_C` and `b_C`, never a causal payment map.

### Cyclic donor provenance and single charging
For one closed helical triad, the three cyclic modal works satisfy `T0+T1+T2=0` before Hahn. With `P_i=[T_i]+`, `N_i=[-T_i]+`, `Q=sum P=sum N`,

`M(i→j)=N_i P_j/Q`

has canonical `dW-` donor marginal and canonical `dW+` recipient marginal. Generic positive recipients may have two donors. Hard-cell pushforward aggregates overlapping donor provenance into each recipient charge exactly once; coarse self-loops are retained as same-time redistribution with zero added depth and no scale progress.

## 3. Between-time energy and radial transport

For any physical helical-mode set `A`, internal cyclic traffic cancels from divergence:

`W_A+ - W_A- = Phi_in,A - Phi_out,A`.

Navier–Stokes gives

`E_A(t1) + D_A + ∫Phi_out = E_A(t0) + ∫Phi_in`.

This is stock + viscosity + boundary flow, not a bound on gross nonlinear traffic.

For the radial exterior `H_R={|k|>R}`,

`Phi_up(R)=M{|k_d|≤R<|k_r|}`,
`Phi_down(R)=M{|k_r|≤R<|k_d|}`,

and

`E_>R(t1)+D_>R+∫Phi_down = E_>R(t0)+∫Phi_up`.

High→high circulation is real but internal. Integrating crossings in `dR/R` gives exact clipped logarithmic donor→recipient displacement. Equiradial positive transfer has zero radial action, so nonlinear transfer itself supplies neither a minimum dyadic step nor an event count.

## 4. Hard-tail true-upward supply

The old hard-tail inequality using gross positive tail work remains a valid coarse bridge, but the native feeder is now known exactly. Multiplying radial continuity at parent-tail scale `N` gives the owner cover

`N E_>N(t0) + N∫Phi_up ≥ 2 nu D_tail`,

hence at least one of

`N E_>N(t0) ≥ nu D_tail`,
`N∫Phi_up ≥ nu D_tail`.

So high-tail dissipation is supported by **inherited tail stock or actual low→high boundary supply**, never by internal high→high circulation.

Disintegrate each upward atom by the recipient shell `M=2^j N`. The physical energy donor is one recipient interaction parent.

- **Pure-UV HH by support:** both recipient interaction parents are `>M/4`. Since the donor is also `≤N`, this forces `j=1`, hence `M=2N`. Triad closure then puts both parents in the comparable range `(M/4, 3M/2]`. Pure-UV true supply is therefore automatically first-shell and comparable; no locality norm estimate is needed to discover that fact.
- **Deep upward crossing:** if `M≥4N`, the donor satisfies `|k_d|≤N≤M/4`; therefore the atom has resolved-scale parent contact and cannot be pure-UV by support.

These are **atomwise** alternatives. One closed-triad law may carry first-shell pure-UV and deep-contact upward submeasures simultaneously. Resolved-scale contact is only a Fourier-support fact; it is not yet called a smooth-interface owner.

The causal unit remains the parent-block unit `N dW`. Recipient-shell scale `M` is geometry, not a new causal reweighting.

## 5. Other certified recurrence controls

- Smooth physical relink is conservative same-event redistribution; only the symmetric strain/deformation part remains a recursive owner.
- Same-carrier checkpoint segmentation cannot reset cumulative first-hit monitors or create Zeno re-hardening.
- Eventually-pure consecutive high-strain recurrence has a descending physical-scale/dissipation telescope.
- Eventually-pure **signed-good generated-HH** recurrence has a finite physical-time backshift telescope. Generic/non-signed-good HH remains outside that theorem.
- Material labels are sidecars to actual service/work. Shannon/Rényi laws govern breadth/reuse of actual positive ancestry, not time.

## 6. Current frontier

The next work must keep the two true-upward support branches separate.

1. **Resolved-contact branch:** prove a positive, type-correct binding through the actual smooth decomposition `u=V+h` before assigning interface/relink/strain ownership. Frequency contact alone is insufficient.
2. **Pure-UV branch:** exploit the already-rigid `M=2N` comparable geometry and bind the donor-restricted canonical positive submeasure into the existing natural-window/critical-shell continuation without re-Hahn or changing the common `N dW` unit.
3. Local HH still has the degenerate full-signed Young/Christ-margin seam and, only if a future representation truly needs it, a separately proved coherent positive kernel.
4. Globally, mixed genuine-owner recurrence, the initial-data interface, and the hypothetical singular-time interface remain open. There is **no 3D Navier–Stokes global-regularity claim**.

## 7. Working discipline and reading order

Before every new theorem branch: fetch/verify `origin/main`, read this file in full, then branch. Theorem/test/numerical/PDE execution belongs in GitHub Actions; local/cloud shell work is for reading, static inspection, editing, git, and artifact/log handling.

Current reading spine:
`physical_energy_causal_bridge.md` → `high_frequency_dissipation_reentry.md` → `helical_physical_edge_registration.md` → `continuum_helical_edge_measure_registration.md` → `canonical_positive_edge_work_routing.md` → `mixed_fate_reserved_young_handoff.md` → `cyclic_helical_triad_donor_kernel.md` → `cyclic_hard_cell_single_charge_quotient.md` → `helical_mode_set_energy_continuity.md` → `radial_spectral_crossing_layer_cake.md` → `hard_tail_true_upward_supply.md` → `master_no_escape.md`.
