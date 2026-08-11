# Canonical positive edge-work routing frontier

Status: **NEXT-THEOREM DESIGN CONTRACT — NOT YET A CERTIFIED THEOREM**.

The continuum helical edge-measure theorem fixes the physical signed law before
any later representation is introduced.  Let

\[
\mu=dW=C_FT_e\,d\Lambda_{edge},
\qquad
C_F=(2\pi)^{-3/2},
\]

on the canonical unordered Fourier/helicity edge space.  Its Hahn decomposition

\[
\mu=\mu^+-\mu^-
\]

is therefore the canonical positive/negative physical child-work law.  The next
problem is **routing this already-existing positive law**, not manufacturing a
new positive law from coherent cells, capacity weights, or theorem labels.

## 1. The Hahn level is physical and must be fixed before analysis refinement

Let `pi` be any later measurable hard/coherent/material label map.  The inherited
causal law is

\[
\boxed{\nu=\pi_\#\mu^+.}
\]

The signed work seen after the same coarsening is

\[
\sigma=\pi_\#\mu.
\]

In general

\[
\boxed{\sigma^+\le \pi_\#\mu^+,}
\]

with strict inequality whenever opposite-sign physical edges cancel inside one
analysis label.  Refining the labels may increase the total positive part of the
coarsened signed representation again.  Neither operation creates or destroys
physical edge work; it only changes how much cancellation is visible after the
chosen representation.

Therefore a later decomposition may do either of two legitimate things:

1. push forward/restrict the already-fixed causal law `mu^+`; or
2. reconstruct signed work for a diagnostic/Young calculation.

It may **not** Hahn-split its own signed cells and silently declare that result a
second master causal law.

For a POVM/coherent localization which is not literally a deterministic map,
the same requirement becomes a positive mass-preserving kernel handoff from
`mu^+`.  If such a kernel identity is not proved, the coherent Hahn atoms remain
an exact signed-work representation but are not yet identified with the master
causal law.

## 2. Partition actual positive work by the edge's own efficiency

For every certified edge define the native signed efficiency

\[
r_e=m_ec_e=\frac{J_e}{J_*}c_e.
\]

Fix the already-certified geometric threshold

\[
\eta_0=10^{-4}.
\]

On the **actual positive edge-work law** define

\[
G=\{T_e>0,\ r_e>1-\eta_0\},
\qquad
B=\{T_e>0,\ r_e\le1-\eta_0\}.
\]

Then this is an exact restriction of one causal law:

\[
\boxed{\mu^+=\mu_G^++\mu_B^+.}
\]

No capacity-majority argument is involved.

Because the one-edge provenance gives `|T_e|<=A_e`, every positive-work subset
has positive capacity whenever it has positive causal mass.  On the bad subset,

\[
\frac{F(B)}{J_*A(B)}
=\frac1{A(B)}\int_B r_e\,dA
\le1-\eta_0,
\]

hence

\[
\boxed{
\epsilon_B=1-\frac{F(B)}{J_*A(B)}\ge\eta_0=10^{-4}.
}
\]

Thus capacity is used only to certify **why the same positive-work subblock is
inefficient**.  The routed causal mass remains `mu_B^+`, not `dA|_B`.

## 3. Geometry-bad positive work has a certified stage-zero transfer-loss route

For the certified coherent range `0<tau<=0.1`, the fixed transfer-loss threshold
is

\[
\delta_\tau=\frac{\tau^2}{1\,036\,800\,000}
\le \frac{1}{103\,680\,000\,000}
<10^{-4}=\eta_0.
\]

Therefore every nonzero bad-geometry causal subblock satisfies

\[
\epsilon_B\ge\eta_0>\delta_\tau.
\]

This scalar inequality is not itself a representation adapter.  The source theorem must bind the measurable restriction `B` of the canonical edge law to the same selected physical block interface used by the transfer compiler; it may not copy `epsilon_B` onto an unrelated coherent block.  Once that same-law block binding is proved, the intended master handoff is

\[
\boxed{
\mu_B^+
\longrightarrow
\texttt{fixed\_transfer\_loss}
\longrightarrow
\texttt{TRANSFER\_WORK\_LOSS}
\longrightarrow
\texttt{TRANSFER\_COST}.
}
\]

This is a **stage-zero fate of the selected block**, not a new first-hit time.
The canonical joint-stop implementation already records it with
`first_time=None` and certificate `stage_zero_fixed_transfer_loss`.

A residual representation bug remains in the fine `physical_branch_compiler`:
its legacy `_first_causal_split()` returns `first_time=0.0` for the same fixed
block fate.  The next source patch should make this `None`.  PDE time `t=0` must
remain reserved for the absorbing initial boundary and must never be reused as a
sentinel for “no event time”.

## 4. Capacity majority is not causal-work majority

The good-core theorem controls capacity mass and, **conditioned on the good
core**, compares normalized `dA` and positive physical work by a certified
Radon--Nikodym bound.  It does not permit a global implication

`good capacity majority => good causal-work majority`.

The obstruction is structural.  Suppose two genuine physical edges have
positive work-to-capacity ratios `rho_g` and `rho_b` with `rho_b>rho_g>0`, the
first geometry-good and the second geometry-bad/nonforward.  Choose capacities
with

\[
1<\frac{A_g}{A_b}<\frac{\rho_b}{\rho_g}.
\]

Then

\[
A_g>A_b
\qquad\text{but}\qquad
T_g=\rho_gA_g<\rho_bA_b=T_b.
\]

So the good edge carries the majority of capacity while the bad edge carries the
majority of actual positive work.  The implementation theorem should encode an
explicit helical NS instance of this inequality as an adversarial regression.

This is why the causal split must be performed directly on `mu^+`.

## 5. Geometry-good is only Young-eligible

The good restriction `mu_G^+` is **not** yet `marking_good=True`.
Geometry/phase efficiency and sharp Young amplitude-profile saturation are
separate questions.

Young/Christ acts on the signed trilinear form of an actual hard product-cell
triple.  If a physical cell `C` contains many canonical edges, its signed work is

\[
T_C=\int_C dW,
\]

not the gross positive edge mass

\[
\int_C dW^+.
\]

Cancellation must therefore be retained before the Young efficiency is computed.
The required next handoff is a **signed cell compression theorem** which proves
that the hard product-cell work used by complex Young is exactly the pushforward
of the same signed continuum edge law with the same unitary normalization and
frozen physical multiplier.

Only after that theorem may one ask whether

\[
\frac{|T_C|}{m_*A_3
\|f_C\|_{3/2}\|g_C\|_{3/2}\|h_C\|_{3/2}}
\]

is close to one.  Low Young saturation is an existing transfer-loss fate; high
Young saturation plus the certified symbol-freezing/phase conditions may then
produce the typed parent marking used by the recursive physical witness.

Hence:

\[
\boxed{
\text{geometry-good}
\neq
\text{Young-good}
\neq
\text{registered generated continuation}.
}
\]

Each arrow needs its own physical theorem.

## 6. Coherent Hahn atomization remains exact, but is not automatically causal

The existing coherent theorem

\[
W_{CDE}=2\Re\langle A_Ew_3,
\mathcal N(A_Cw_1,A_Dw_2)\rangle,
\qquad
\sum_{CDE}W_{CDE}=W_{HH},
\]

is an exact signed Navier--Stokes work identity.  Its positive and negative
atomic parts correctly expose cancellation inside that coherent representation.

After the continuum theorem, however, the master-facing interpretation must be
more precise.  The canonical causal law already exists as `mu^+` on physical
Fourier/helicity edges.  A coherent/material parent label must either:

- inherit `mu^+` through a proved positive pushforward/disintegration; or
- remain a signed diagnostic/Young representation until such an identification
  is proved.

This prevents arbitrary coherent-cell refinement from changing causal mass while
preserving all exact coherent work identities.

## 7. Required source-level theorem shape

A natural next production module is

`src/canonical_positive_edge_work_routing.py`.

Its input should be the replayable physical fibers/atoms of the certified
continuum ledger, not caller-supplied masses or deficits.  A minimal typed output
should record:

- total canonical positive edge-work mass;
- geometry-bad positive-work mass;
- geometry-good positive-work mass;
- exact mass-reconstruction residual;
- bad-subblock native capacity and deficit;
- same-law selected-block binding plus terminal `TRANSFER_WORK_LOSS` projection for the bad sublaw;
- a **Young-eligible but not Young-certified** handle for the good sublaw.

The theorem must fail closed under forged summary fields by replaying the bound
physical edge identities exactly as the continuum ledger already does.

The first implementation should also align the fine compiler with the preferred
joint-stop semantics by replacing the fixed-loss sentinel `first_time=0.0` with
`None`.

## 8. Adversarial regressions required before promotion

The theorem is not ready merely because random states pass.  CI should include:

1. a physical helical counterexample with good-capacity majority but bad-work
   majority;
2. positive nonforward work (`J=0`) routed to bad geometry, never silently
   discarded;
3. exact `mu_G^+ + mu_B^+ = mu^+` mass reconstruction;
4. exact `epsilon_B>=1e-4` on every nonempty bad sublaw;
5. all certified `tau<=0.1` implying the stage-zero transfer-loss gate;
6. parent swap, helical gauge and unit rescaling invariance;
7. provenance tampering rejected by replay;
8. analyst coarsening/refinement unable to replace the inherited causal
   pushforward by a new Hahn law;
9. `fixed_transfer_loss` returning `first_time=None`, while an actual
   `INITIAL_BOUNDARY` hit remains the only absorbing `t=0` case;
10. geometry-good data unable to set `marking_good=True` without a separate
    signed-cell Young certificate.

## 9. Master consequence if this theorem is certified

This theorem would not yet close generic HH recurrence.  It would make the
remaining branch much smaller and more physical:

\[
\boxed{
\text{canonical }dW^+
\to
\begin{cases}
\text{geometry/phase bad} &\to \text{terminal transfer loss},\\
\text{geometry/phase good} &\to \text{signed hard-cell Young test}.
\end{cases}}
\]

The first branch would no longer recurse.  The second would still require the
signed-cell Young/Christ handoff and then the already-existing common-slice /
physical generated-ancestry machinery.

No norm becomes a physical owner.  No capacity mass becomes causal probability.
No analyst refinement creates positive work.  The entire split is read from the
same signed Navier--Stokes interaction law that generated the event.
