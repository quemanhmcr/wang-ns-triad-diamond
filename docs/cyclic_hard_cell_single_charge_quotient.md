# Cyclic hard-cell single-charge quotient

Status: **CERTIFIED** on exact theorem SHA `964b72b0a614bb44c88a42d41160f11bd550d97a`. Dedicated run `31582837653`, independent audit `31582837613`, and full physical-energy causal integration `31582837580` all completed successfully on that SHA.

The certified cyclic helical-triad donor theorem already identifies canonical Hahn-negative edge work as same-time energy-donor work.  Its recipient marginal is the already-fixed canonical Hahn-positive edge law.  The present theorem asks the next narrower question: what happens when the analyst observes those physical roots only through the deterministic hard cells used by the Young/Christ handoff?

The answer must not be a new cancellation budget.  It must not say that negative work “pays” failed good work.  It must simply push the already-certified physical donor kernel through the same hard representation and preserve exactly who receives the energy.

---

## 1. Start from the physical closed-triad measure, not from a hard-cell scalar

Let `M_triangle` denote the certified positive donor/recipient measure on one regular closed helical triad.  Its donor marginal is canonical `dW-` and its recipient marginal is canonical `dW+`:

\[
(\mathrm{donor})_\#\mathcal M_\triangle=dW^-,
\qquad
(\mathrm{recipient})_\#\mathcal M_\triangle=dW^+.
\]

This identity is already downstream of exact signed three-root Navier--Stokes energy conservation

\[
T_0+T_1+T_2=0
\]

**before Hahn**.  No hard representation appears in the construction of `M_triangle`.

Let

\[
\pi:e\mapsto C
\]

be the deterministic hard Fourier/helicity product-cell map already used by canonical edge routing and the mixed-fate Young handoff.  On positive recipient edges retain the already-certified physical fate map

\[
\phi(e)=
\begin{cases}
G,& (J_e/J_*)c_e>1-10^{-4},\\
B,& (J_e/J_*)c_e\le1-10^{-4}.
\end{cases}
\]

`G` means only Young-eligible.  `B` means the same existing stage-zero `TRANSFER_WORK_LOSS` sublaw.  The fate map is not recomputed from a hard-cell Hahn split.

---

## 2. The hard-cell donor table is only a pushforward

Define

\[
\boxed{
K(C,D,F)
:=
(\pi_{donor},\pi_{recipient},\phi_{recipient})_\#
\mathcal M_\triangle(C,D,F).
}
\]

This is the whole new object.  It is a positive pushforward of the physical same-triad donor measure.

Its donor rows satisfy

\[
\boxed{
\sum_{D,F}K(C,D,F)
=(\pi_\#dW^-)(C)=:n_C.
}
\]

Its recipient columns retain the canonical fate split:

\[
\boxed{
\sum_C K(C,D,G)
=(\pi_\#dW_G^+)(D)=:g_D,
}
\]

\[
\boxed{
\sum_C K(C,D,B)
=(\pi_\#dW_B^+)(D)=:b_D.
}
\]

Therefore

\[
\sum_C K(C,D,G)+\sum_C K(C,D,B)
=(\pi_\#dW^+)(D).
\]

Nothing has been Hahn-split again.  Nothing has been normalized into a probability.  Capacity never enters the transport table.

---

## 3. Why this is an exact single-charge theorem

Suppose the donor hard cells form a disjoint measurable partition `{C_alpha}`.  Restrict the donor side and push to recipients:

\[
\nu_{C_\alpha}^+(D,F)=K(C_\alpha,D,F).
\]

Each restricted law is a positive submeasure of canonical recipient work and has exactly the donor-cell negative mass:

\[
\boxed{
\nu_{C_\alpha}^+(\Omega)=n_{C_\alpha}.
}
\]

Different donor cells may have overlapping recipient support.  That is physical: two energy donors may feed the same recipient mode at the same triad/time.  But the measures recombine exactly as

\[
\boxed{
\sum_\alpha \nu_{C_\alpha}^+(D,F)
=(\pi_\#dW_F^+)(D).
}
\]

Thus an already-canonical recipient charge is **disintegrated among its donor provenances**, not copied once per donor.  The downstream compiler sees the recipient cause once, with exactly its canonical `dW+` mass.

This is distinct from the older `physical_branch_compiler` single-charge rule.  That compiler quotients duplicate theorem manifestations of an already selected positive causal root.  The present theorem proves that cyclic negative-work provenance itself reaches those positive roots without multiplying their measure.

---

## 4. Generic two-donor recipient is not an ambiguity

The previous cyclic theorem proved the anti-theorem that a generic positive child may have two energy donors.  In that sign pattern there is one positive recipient root `j` and two negative donor roots `i_1,i_2`.

At fine hard resolution,

\[
K(C_{i_1},D_j,F_j)>0,
\qquad
K(C_{i_2},D_j,F_j)>0.
\]

But

\[
K(C_{i_1},D_j,F_j)+K(C_{i_2},D_j,F_j)
=(\pi_\#dW^+_{F_j})(D_j).
\]

The recipient is one physical cause with two donor provenances.  Unique donor remains false generically, while single recipient charge remains exact.

---

## 5. Signed-good one-donor/two-recipient triad exposes both downstream fates

On the signed-good forward triad, exactly one interaction parent is the energy donor.  Its loss splits into the forward signed-good child and the positive nonforward side mode:

\[
W_{donor}^-=W_{child}^++W_{side}^+.
\]

The child recipient is geometry-good and remains Young-eligible.  The side recipient has `J=0`, hence is geometry-bad and remains on the existing positive-nonforward `TRANSFER_WORK_LOSS` route.

The hard-cell quotient therefore gives, from one donor row,

\[
\nu_C^+=\nu_{C,G}^++\nu_{C,B}^+,
\]

with both pieces real positive Navier--Stokes work.  The `B` fate terminates that positive sublaw in the **recursive forward-transfer routing**.  It does not mean the side-mode energy disappeared from the PDE.  That modal energy may participate in later ordinary Navier--Stokes interactions.

This distinction is mandatory: recursive route termination is not viscous dissipation.

---

## 6. Coarse self-loops are real and carry zero recursion depth

A deterministic hard map may identify donor and recipient roots:

\[
\pi(e_{donor})=\pi(e_{recipient})=C.
\]

Then

\[
K(C,C,F)>0.
\]

This is not a contradiction and must not be deleted.  The underlying physical roots are distinct; only the coarse analyst labels coincide.  The self-loop is same-time closed-triad energy redistribution viewed below the resolution of the hard map.

It therefore has:

- real positive physical work mass;
- zero additional physical event time;
- zero recursive generation depth;
- no supplied scale progress.

In particular, a hard-cell self-loop defeats any shortcut that tries to read cyclic donor provenance itself as a descending/ascending scale map.

---

## 7. What the theorem does **not** infer from the mixed-fate failure inequality

The mixed-fate theorem proved, under a nondegenerate full-signed Christ margin,

\[
g_C<n_C+(\mu^{-1}-1)b_C
\]

when reservation fails.

The present donor quotient gives physical meaning to `n_C`: it is donor-side manifestation of same-time canonical positive recipient causes.  But the inequality above is still only a scalar domination statement.  It does **not** define a map from failed good work in cell `C` to the donor work `n_C`, and it does not permit the declaration

> “negative work pays the failed good branch.”

The recipient of `n_C` may lie in another hard cell, may be geometry-bad, may be geometry-good, or may even return to the same coarse hard cell.  Its fate is determined by the already-existing recipient edge cause, not by the failure inequality.

---

## 8. Same-time redistribution and between-time modal inventory remain different ledgers

The cyclic kernel answers a same-time question:

> Which simultaneous positive modal work receives the energy withdrawn by this negative modal work?

It does not answer the between-time question:

> Which earlier positive deposit supplied the modal energy later withdrawn at another interaction time?

The latter is governed only at aggregate level by the physical modewise energy balance

\[
E_k(t_1)-E_k(t_0)
=
\int_{t_0}^{t_1}(W_k^+-W_k^-)\,dt
-
2\nu|k|^2\int_{t_0}^{t_1}E_k\,dt.
\]

No canonical FIFO/LIFO pairing of deposits and later withdrawals follows from this identity.  This theorem deliberately creates no such matching.

---

## 9. Floating certification near exact phase cancellation

The exact theorem applies to every mathematically nonzero cyclic work law.  Numerically, the sign of an extremely phase-cancelled work vector is ill-conditioned.  The predecessor theorem already established the correct policy: use native modal capacity only as an immutable floating **error envelope**, never in the transport law, and mint no donor atom when the realized work is below numerical sign resolution.

The hard-cell quotient inherits that policy.  It refuses unresolved kernels rather than manufacturing hard-cell provenance from rounded Hahn signs.  Submeasure-domination residuals are also certified on the native work-mass envelope, never by dividing through a possibly tiny realized recipient Hahn mass.

---

## 10. Actual Navier--Stokes audit

The companion PDE probe evolves the repository's real `2/3`-dealiased incompressible Fourier--Galerkin Navier--Stokes system.  At each physical snapshot it:

1. reads the actual evolving modal coefficients of the selected closed helical triad;
2. reconstructs all three signed root works before Hahn;
3. constructs the certified same-triad donor kernel;
4. pushes donor/recipient roots through both a fine deterministic hard map and a maximal coarse map;
5. checks `dW-` row marginals, `dW+` good/bad column marginals, restricted-donor submeasure conservation, and coarse self-loop conservation;
6. compares the same cutoff-7 Galerkin system on FFT grids `24` and `28` using the native physical work scale.

The audit also runs the negated divergence-free initial field.  Because nonlinear work is cubic, this reverses the initial triad work signs and exposes the opposite physical sign pattern: the ordinary branch begins one donor/two recipients, while the sign-reversed branch begins two donors/one recipient.  After that initial condition is set, both branches are evolved by the same Navier--Stokes equations.

---

## 11. Scope after certification

A successful theorem would close the **representation-level negative-work single-charge seam**:

\[
\text{hard donor }\pi_\#dW^-
\longrightarrow
\text{existing hard recipient }\pi_\#dW^+
\]

with good/bad recipient fate preserved exactly once.

It would **not**:

- prove every geometry-good recipient Young-good;
- resolve a degenerate Christ modulus;
- turn negative work into a budget, reset, clock, or probability;
- infer a causal map from failed good work to negative work;
- create scale progress from a coarse self-loop;
- match earlier modal deposits to later withdrawals;
- terminate generic HH or mixed genuine-owner recurrence;
- prove 3D Navier--Stokes global regularity.

That separate between-time theorem is now certified in `helical_mode_set_energy_continuity.md` on exact SHA `a39d502d9312ac3bd6613a780d60b22e88790863`: persistent stock lives on physical helical modes and satisfies exact mode-set stock/boundary-flow/viscosity continuity without FIFO/LIFO matching or a gross-transfer budget.  The next scale-facing specialization is a radial Fourier boundary.


---

## Certification record

The final certified theorem SHA is `964b72b0a614bb44c88a42d41160f11bd550d97a`.
The dedicated run `31582837653` passed `856` tests, `75,000` hard-cell single-charge stress states, both evolved Navier--Stokes sign patterns on FFT grids `24,28`, and amplitude adversaries.  Independent audit `31582837613` passed `100,000` hard-cell states plus longer positive/negative-phase evolved NS; full integration `31582837580` passed `856` tests and a `50,000`-state theorem run in the causal spine immediately after cyclic donor provenance and before complex Young.

The independent `100,000`-state audit retained `31,762` cases in which distinct donor hard cells overlapped one recipient charge and `44,373` cases with at least one coarse self-loop; worst restricted-donor pushforward native residual was `1.075e-16`.  The longer actual-NS audit retained the one-donor/two-recipient and phase-reversed two-donor/one-recipient patterns with cross-FFT charge residuals at about `1e-16`.

Initial candidate `55a54b02f0a7a77a4b7a73b0403e90311cf0eeb8` is retained in the failure lineage.  It passed the focused theorem and new PDE lanes, but the independent 100k stress exposed a floating guard that compared restricted submeasure domination against an accidentally tiny realized recipient Hahn mass.  Final SHA `964b72b...` uses the certified native work-mass scale only as a floating error envelope for that comparison; it does not place capacity in the kernel or causal law.  The same repair removed one unsupported `--cutoff` argument from an existing dependency workflow.  No theorem threshold, sample count, physical transport, fate map, or causal semantics was weakened.
