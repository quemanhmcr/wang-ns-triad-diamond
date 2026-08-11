# Same-carrier checkpoint segmentation quotient

## Status and proof boundary

Audit repair candidate. Serious execution, stress, and Navier--Stokes numerical
falsification are reserved for GitHub Actions.

This lemma concerns one already-defined smooth event-anchored carrier on one
pre-singular PDE trajectory. It does not construct such a trajectory beyond its
smooth interval, does not terminate genuine owner recurrence, and does not prove
3D Navier--Stokes regularity.

## 1. The object being segmented is a physical path, not an array label

Fix one genuine physical event at terminal time `t`, one smooth carrier `Q_A`, one
terminal dual, and the actual Navier--Stokes trajectory on which they are followed.
The audit representation binds all of

- physical event token;
- carrier token and frequency `A`;
- terminal-dual token and complex terminal coefficient `z(t)`;
- scaled lifetime `c`;
- terminal PDE state token and trajectory token.

These data form one `SameCarrierProvenance`. Equality is exact. A repeated string
such as `"fixed-Q"`, or dimensional quantities that happen to be numerically
close, cannot prove that two rows belong to the same carrier.

For elapsed backward time `ell=t-s`, the cumulative observables are

\[
K_A[s,t],\qquad I_R[s,t],\qquad I_{HH}[s,t].
\]

`K_A` is nonnegative and monotone. The coefficient impulses are stored as their
actual complex cumulative paths. Their magnitudes are derived only when the
closed faces

\[
|I_R|=|z(t)|/4,\qquad |I_{HH}|=|z(t)|/2
\]

are located. This retains phase cancellation and prevents segment magnitudes from
being added or promoted to work.

## 2. Exact finite-segmentation lemma

Let a finite family of segments be restrictions of the preceding path. At every
shared boundary require exact equality of

- the complete carrier provenance;
- the PDE state token;
- native elapsed time;
- cumulative strain action;
- both cumulative complex impulses.

The first row must be the terminal PDE state at elapsed time exactly zero, with
all three cumulative monitors exactly zero. Concatenation then deletes one copy of
each shared boundary and reconstructs the original path identically. Consequently
inserting or deleting finitely many observer cuts leaves unchanged

- the first native elapsed time;
- the complete exact joint first-stop set;
- whether coefficient-energy reentry is required.

This is a restriction/gluing identity. It is not evidence that independently
generated arrays belong to one PDE solution. The provenance and boundary tokens
are the hypothesis that supplies that fact.

## 3. Why complex phase belongs in the certificate

Linear interpolation of endpoint magnitudes is not the magnitude of a linearly
interpolated complex impulse. For example, the complex chord

\[
I(\theta)=0.5-2.5\theta
\]

first reaches `|I|=1` at `theta=0.6`; interpolating the endpoint magnitudes
`0.5 -> 2` would incorrectly report `theta=1/3`.

The regression implementation therefore solves the quadratic circle intersection

\[
|I_0+\theta(I_1-I_0)|^2=\lambda^2
\]

on every complex chord. The analytic theorem remains a statement about the actual
continuous/absolutely-continuous complex impulse, not about this discretization.
An optional numerical tie tolerance is explicitly regression-only; the default
exact first-stop set does not merge distinct native debut times.

## 4. A natural service window is not an arbitrary observer cut

For the fixed carrier, `A` and `c` do not change. Every genuine full-natural
window therefore has the same positive native duration

\[
T_A=cA^{-2}>0.
\]

Starting from a finite event time `t`, at most

\[
\left\lfloor t/T_A\right\rfloor
\]

complete fixed-carrier natural windows fit before `t=0`. Hence genuine
fixed-carrier natural windows cannot have a Zeno accumulation at a positive
physical time.

An analyst may, of course, insert infinitely many plotting or observation cuts
which accumulate in the interior. Those cuts have no service-window status and
cannot each be assigned a fictitious duration `c A_j^-2`. Likewise, shell scales
read as state sidecars at those cuts do not become a carrier-scale lineage.

This distinction removes the apparent contradiction cleanly:

- fixed natural windows use their actual positive PDE durations;
- arbitrary cuts use no physical duration or causal charge at all.

## 5. Interior accumulation is conditional on the actual path

Endpoint limit scalars alone cannot establish a first stop at an accumulation.
The implementation requires a typed prelimit certificate containing the actual
cumulative path from the terminal event to the proposed limit.

There are then three cases.

1. If a closed face was crossed earlier, the premise “no-hit checkpoints
   accumulate here” is false and is rejected.
2. If the first face occurs exactly at the endpoint, that endpoint is the existing
   first stop. Coefficient faces still require physical `Q^2` energy reentry and
   are not work.
3. If no face occurs and the endpoint is interior, continuation is asserted only
   with a matching token certifying that the same PDE trajectory is smooth on an
   open interval around the endpoint. Continuity then preserves the strict margins
   on a smaller open interval.

Values strictly above a threshold cannot be called a first contact at the limit:
continuity would have produced an earlier contact. A `max(1,...)` absolute
tolerance is forbidden because it destroys this statement at small native scales.

At `t=0` the initial boundary absorbs only when cumulative native elapsed time
equals the full event time exactly. Positive remaining native time is not rounded
away by an observer-unit floor.

## 6. Typed checkpoint policy

A few dictionary flags such as `physical_event_created=False` do not constitute a
checkpoint theorem. The continuation policy accepts an actual
`FullNaturalCheckpoint` and verifies exact agreement with the carrier's event time,
frequency `A`, and lifetime `c`.

The checkpoint's hard-shell readings are state sidecars. They cannot

- replace the event-anchored carrier or terminal dual;
- reset the terminal coefficient or cumulative monitors;
- create recursion depth, causal charge, or a physical scale lineage.

Actual high-frequency dissipation/work remains a genuine UV route whenever its
independent physical hypotheses hold. The quotient removes only motion created by
observer rereading.

## 7. Numerical falsification lane

The audit lane evolves the unforced 3D incompressible Fourier--Galerkin
Navier--Stokes system with Leray projection, viscosity, `2/3` dealiasing, and RK4.
On one and the same evolved trajectory it fixes `Q_A` and a terminal dual, forms
the cumulative complex nonlinear pieces and the resolved strain action, then
compares whole-path and multiply segmented first stops.

The experiment must also verify incompressibility, the global energy balance, the
carrier `Q^2` balance, the complex coefficient/Duhamel residual, low--low moat,
nonzero nonlinear activity, partition invariance, and resolution refinement. This
is numerical falsification evidence, not a continuum proof.

## 8. Remaining frontier

After exact same-path segmentation is removed, a free branch may end only through
an actual named physical stop, `t=0`, or continuation of the same smooth PDE path.
The global frontier is still recurrence of genuine physical owners. Those events
must telescope through the native work, source/service, reuse, `Xi`, or true global
resource laws they actually possess.

No artificial norm, synthetic scale tax, or checkpoint count is introduced.
