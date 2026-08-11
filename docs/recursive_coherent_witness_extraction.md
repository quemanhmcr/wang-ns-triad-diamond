# Recursive coherent witness extraction: binary causality is already inside quadratic work

The remaining continuum audit should not synthesize a family of packets and then
ask whether the resulting graph resembles Navier--Stokes causality.  Once the
outer transfer-selected divergence-free roles are fixed, the quadratic PDE and
the exact coherent resolution of identity already contain a canonical binary
event space.

This note isolates that exact statement.

## 1. Start with the actual selected quadratic source

For one selected child role let the designated parent roles be `w_1,w_2` and
write the exact high--high source abstractly as

\[
\mathcal N(w_1,w_2).
\]

For Navier--Stokes this is the already defined Leray-projected bilinear source

\[
-Q_c\mathbb P\nabla\cdot
(w_1\otimes w_2+w_2\otimes w_1),
\]

with the strict low-pass transporter separated before this step.  Pressure is
not reintroduced.

Let `{A_C}` and `{A_D}` be coherent localization partitions for the two parent
roles.  Since

\[
\sum_C A_C=I,
\qquad
\sum_D A_D=I,
\]

bilinearity gives exactly

\[
\boxed{
\mathcal N(w_1,w_2)
=
\sum_{C,D}\mathcal N(A_Cw_1,A_Dw_2).
}
\]

There is no Gaussian coefficient choice and no reconstruction error.

## 2. Localize the child only at the work level

Let `{A_E}` be the same canonical coherent type of partition on the child role.
Define

\[
\boxed{
W_{CDE}
=2\Re\left\langle
A_Ew_3,
\mathcal N(A_Cw_1,A_Dw_2)
\right\rangle.
}
\]

Because every `A_E` is self-adjoint and the three partitions resolve the
identity,

\[
\begin{aligned}
\sum_{C,D,E}W_{CDE}
&=2\Re\left\langle
w_3,
\mathcal N(w_1,w_2)
\right\rangle.
\end{aligned}
\]

Thus the triple-indexed atoms reconstruct the **actual child-energy work**
exactly.

This is deliberately a work identity, not a claim that `A_C w` is an invariant
packet or a compactly Fourier-supported solution.  The outer role carries the
strict frequency/helicity information.  Coherent cells carry material
phase-space ancestry.

## 3. Positive atoms are a physical binary causal measure

Take the Hahn split

\[
P=\sum_{C,D,E}[W_{CDE}]_+,
\qquad
N=\sum_{C,D,E}[-W_{CDE}]_+.
\]

Then exactly

\[
\boxed{P-N=W_{HH}},
\]

and therefore

\[
\boxed{P\ge[W_{HH}]_+.}
\]

Every positive atom has, before any graph construction,

- parent material cell `C`;
- parent material cell `D`;
- child material cell `E`;
- one physical interaction time;
- positive actual child-work mass `[W_CDE]_+`.

So after normalization the atoms define a probability law on **binary parent
slots conditioned on a child event**.  This is precisely the structure required
by Shannon/Rényi causal reuse.

The graph is therefore not imposed on the flow.  It is the support of the
positive part of the quadratic Navier--Stokes work measure.

## 4. Backscatter is not hidden

The negative mass `N` is physical backscatter/cancellation.  Taking positive
parts atom by atom can create more total positive mass than the positive part of
the aggregate signed work, but this is not double counting: the exact identity
`P-N=W_HH` records the compensating negative work.

This is the same physical distinction already used by the coherent work ledger.
No absolute forcing norm is introduced.

## 5. Combine with the physical-energy causal gate

The previous theorem gives, on the low-strain generated branch,

\[
W_{HH}^+\ge\frac{8}{15}E_1.
\]

Pointwise coherent atomization followed by time integration gives

\[
\boxed{
\int\sum_{C,D,E}[W_{CDE}(t)]_+dt
\ge W_{HH}^+
\ge\frac{8}{15}E_1.
}
\]

Thus a generated node already carries a quantitatively nontrivial positive
binary coherent work measure.  The amplitude Duhamel measure is not needed to
supply weights.

## 6. The selected physical cross-cell `Xi` is an excision, not another measure

Let a measurable subset of positive triple atoms be omitted by the existing
physical defect moat.  Its positive mass is `Xi_cross`.  Then tautologically

\[
\boxed{
P=P_{ret}+\Xi_{cross}.
}
\]

The moat theorem supplies the quantitative upper bound for this omitted
**physical transfer mass**.  If at one generation it removes a relative fraction
`rho`, then

\[
\boxed{
P_{ret}
\ge(1-\rho)\frac{8}{15}E_1.
}
\]

No packet-count factor and no continuous-to-discrete synthesis term appears.

## 7. Canonical material labels survive recursion

Use the already fixed intrinsic address

\[
\zeta=(L^{-1}X/2,L^Tk).
\]

The three indices `(C,D,E)` are addresses in this same nested dyadic material
hierarchy.  Common affine/Kelvin motion preserves the address, nested refinement
is exactly additive, and representative frequency/covariance changes have their
already summable `Xi_sym`/`Xi_cov` bounds.

At the next synchronized common slice one may reconstruct the exact coherent
analysis from the actual Navier--Stokes field.  One does **not** need to propagate
`A_Cw` as a frozen Gaussian packet through the whole parent lifetime.

A genuine material-cell switch is still a physical relink/backflow/fresh event.
It is not disguised as representation error.

## 8. Why this avoids an important false shortcut

A coherent localization operator is a positive POVM element, not an orthogonal
Fourier projection.  An individual `A_Cw` need not have compact Fourier support.
Therefore we do **not** infer signed-good scale geometry from the synthesized
piece itself.

Scale, helicity and the strict `N/4` transporter belong to the outer selected
Fourier/helical role.  When a representative frequency inside a coherent cell is
used by Hodge/helical geometry, the existing smooth-symbol and covariance
schedules supply the only representation `Xi`.

This separation is essential: coherent cells answer **where/which material
ancestry carries the work**; Fourier/helical roles answer **which interscale
interaction the work represents**.

## 9. What this closes and what remains

This theorem closes the binary parent-pair **witness extraction after the outer
selected PDE roles have been supplied**:

\[
\boxed{
\text{selected physical HH work}
\longrightarrow
\text{exact positive binary coherent event measure}
}
\]

with no packet synthesis, no persistence assumption and no new `Xi`.

At the stage of this theorem, the remaining continuum bridge became more specific.  For every recursively
selected efficient smooth-SGS block one must construct the outer
Fourier/helical roles and prove that their exact moving equation has:

1. the strict low-frequency Kelvin transporter already used by service-or-flat;
2. the designated high--high source above;
3. every non-affine/moving-role remainder classified once into the existing
   source/sideband/transfer/representation ledgers;
4. the same physical transfer normalization used by the block selector.

This note does not assert that this final outer-role extraction has already been
proved for every continuum block, and it makes no global-regularity claim.

## Outer-role update

The outer-role condition in the historical scope is now supplied by `outer_moving_role_extraction.md` and `event_anchored_role_registration.md`.  The nonaffine Heisenberg interface has exact provenance by `nonaffine_role_interface_work.md`: skew work is conservative role flux and symmetric work is off-diagonal strain.  Binary coherent atomization therefore no longer waits on outer-role extraction; it waits only on the recursive first-stop assembly which decides which already-named cause owns each retained physical event.


## Current downstream status

The moving selected-role equation, event hard-role registration, nonaffine interface split, physical pair-work productivity, measurable first-hit extraction, and generic shell/service reentry requested by the historical scope are now supplied by downstream companion theorems.  High-frequency regeneration also has a dedicated physical continuation through common-unit work, Fourier UV locality, and sliding natural time.  This witness theorem therefore remains the exact coherent work-atomization/causal-support layer; the current programme-level frontier is final continuum master assembly rather than construction of another coherent packet state.

## Master-facing caveat after continuum edge-measure registration

The coherent identity in this note remains exact: the signed atoms reconstruct
quadratic Navier--Stokes work and their positive/negative parts expose the
cancellation visible at this coherent resolution.  The continuum edge theorem,
however, fixes a more primitive master-facing causal law before coherent
refinement: the Hahn-positive part `mu^+` of the canonical signed unordered
Fourier/helicity edge measure.

Therefore the statement that positive coherent atoms form a causal binary law
must now be read conditionally.  For the master, a coherent/material binary law
must be proved to be a positive pushforward/disintegration of that already-fixed
`mu^+`; otherwise the coherent Hahn split is an exact work representation but not
a second independent causal law.  This prevents analyst-chosen cell refinement
from changing causal mass while preserving every signed work identity in this
module.

The next theorem is specified in
`docs/canonical_positive_edge_work_routing_frontier.md`.
