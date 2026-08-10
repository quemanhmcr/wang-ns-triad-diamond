# NN-critical heat law to smooth shell carrier seeds

## Status

**EXACT_NN_CRITICAL_HEAT_LAW_TO_SMOOTH_SHELL_CARRIER_SEEDS__NO_CELL_MASS_FLOOR__TEMPORAL_RENEWAL_REMAINS**

The old-incident theorem supplies, on every sufficiently old already-supplied signed-good material epoch, a positive heat sublaw which is simultaneously

- `NN` in coherent material endpoint provenance;
- on the critical resolved shell-time set `G`.

Its mass is at least

\[
\frac14e^{-1/32}S_{heat}.
\]

The next step should not select the largest coherent cell.  The positive law itself is already the physical selector.

## 1. Push the law to shell-time, not to a maximizing packet

Every atom of the `NN intersect G` heat law has an exact deterministic dyadic shell

\[
A_j=\{M/2<|\xi|\le M\}
\]

and a physical time `t`.  Normalize this positive heat sublaw and push it forward by the map which retains `(j,t)` and the NN endpoint provenance.

No atom needs a uniform mass lower bound.  Arbitrarily small atoms remain legitimate members of the probability law.  If many coherent heat edges have the same `(j,t)`, pushforward merely coalesces their weights and preserves total physical heat mass.

Membership in `G` supplies the independent critical fact

\[
M\|P_j u(t)\|_2^2\ge\mu_*,
\qquad
\mu_*=\frac{32\pi^2}{75c^2}.
\]

The heat edge is built from `P_jV`, while this critical mass is for `P_ju`.  Keep these as **two simultaneous exact marks**.  The theorem does not identify the whole `u` shell with NN material.

## 2. The shell itself determines the renewed scale

Set

\[
\boxed{A=\frac34M.}
\]

Then

\[
M/2=2A/3,
\qquad
M=4A/3,
\]

so the complete hard shell lies in

\[
2A/3<|\xi|\le4A/3.
\]

At this scale its physical critical mass is

\[
A\|P_j u(t)\|_2^2
=\frac34M\|P_j u(t)\|_2^2
\ge\frac34\mu_*
=\boxed{\frac{8\pi^2}{25c^2}}.
\]

This is a whole-shell critical coefficient, not a coherent-cell floor.

## 3. Register the whole shell into a smooth PDE carrier

Choose a scalar smooth Fourier envelope `Q_A` which is identically one on the hard shell and whose support has lower edge at least `3A/5` and upper edge at most `3A/2`.  The hard shell sits strictly inside this envelope because

\[
3/5<2/3<4/3<3/2.
\]

Let

\[
f=P_j u(t),
\qquad
\psi=\frac{f}{\|f\|_2}.
\]

Since `Q_AP_j=P_j`,

\[
\boxed{
\langle Q_Au(t),\psi\rangle
=\langle P_ju(t),\psi\rangle
=\|P_ju(t)\|_2.
}
\]

Thus the smooth Fourier carrier has an exact critical terminal coefficient inherited from the entire physical shell.  The dual direction is the shell's own normalized state; no spatial packet, Gaussian argmax, or coherent-cell representative is chosen.

## 4. The seed already has the outer-role support moat

Use the renewed transporter

\[
V_A=S_{A/4}u.
\]

Its low--low output lies below `A/2`.  If the renewed strain action stays in the existing safe corridor `K<=1/30`, the smooth carrier lower edge stays above

\[
\frac35e^{-1/30}A>\frac12A.
\]

So the exact outer moving-role identity applies with strict persistent low--low exclusion.

The scale separation is stronger than required for a near-scale generated parent.  Since `M<=N/4`,

\[
A\le\frac{3N}{16},
\]

hence

\[
\boxed{
\frac{T_A}{T_N}=\left(\frac NA\right)^2\ge\frac{256}{9}.
}
\]

The high-strain event therefore exposes a much longer-lived lower-scale carrier seed.

## 5. Material provenance is retained but not overclaimed

Each seed is selected by an actual NN coherent heat edge.  Retain those two NN intrinsic-zeta endpoint marks as provenance on the seed law.

Do **not** infer that all energy of `Q_Au` or even all energy of `P_ju` lies in NN material.  The strict resolved multiplier has only been used as a contraction in the high-strain ancestor theorem, so no inverse multiplier bound is available.  The correct statement is weaker and exact:

> the same physical heat atom carries an NN material edge mark and a critical `u`-shell mark, and the latter canonically supplies a smooth lower-scale Fourier carrier with a critical coefficient.

This is precisely enough to remove packet selection from the eventwise entrance to renewal while leaving the remaining material attachment honest.

## 6. Scope and downstream status

This theorem supplies the positive law of smooth lower-scale carrier **seeds** and deliberately stops before propagating them.  The temporal seam identified here is now closed downstream by `nn_seed_temporal_first_stop.md`, followed by `critical_annular_carrier_service_reentry.md`; a seed either reaches a named first stop / `t=0` or, on a full no-hit natural corridor, yields an own-scale service lower.

Material attachment was also simplified downstream: `material_label_carrier_quotient.md` shows that pure material-label changes do not create a second carrier impulse, and material ownership is reread from new actual service.  The NN heat edge remains provenance; it is not promoted to ownership of the entire renewed shell.

No per-cell mass floor, full-shell NN assertion, or global-regularity conclusion is made.
