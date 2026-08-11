# Full-natural service corridor quotient

## Candidate status

Candidate theorem awaiting exact-SHA GitHub Actions certification.

The theorem is deliberately an ontology/assembly correction, not a new scalar
estimate.  It asks where the already-certified own-scale service actually lives
in the Navier--Stokes evolution.

## 1. One full survivor is one physical corridor

For a generic critical shell at frequency `M`, the smooth renewal scale is

\[
A=\frac34M.
\]

If the backward first-stop monitors remain clean for the whole natural interval,
the PDE has traversed the actual physical corridor

\[
I=[t-cA^{-2},t].
\]

The existing critical-shell theorem proves on this same corridor both persistent
carrier energy and a positive bounded heat/increment service law.  Therefore the
classification

`full_natural_own_scale_service`

must not be interpreted as a second event after the corridor.  The service is a
physical observable of the interval which has already consumed the time drop
`cA^-2`.

The corridor object carries the complete native provenance supplied by the PDE
theorem: the incoming hard-shell frequency `M`, renewal-carrier frequency
`A=3M/4`, scaled lifetime `c`, terminal time, endpoint time, actual time drop,
and endpoint carrier-mass lower.  A downstream reader may verify these data but
cannot provide replacements for them.  In particular, a carrier proved near one
`A` cannot be rebound to a foreign hard-shell scale by passing a new number to an
adapter.

Every numerical equality guard is relative to the like-dimensional physical
quantity.  No `max(1,...)` absolute floor is permitted.  This matters on the UV
route because `cA^-2` and the corresponding critical masses can be arbitrarily
small under natural rescaling.

## 2. Materiality is a disintegration of the same law

After the positive service exists, Moyal gives its physical phase-space edge
measure.  OO/ON/NN is then the exact endpoint partition of that nonnegative
measure.  Consequently

\[
S_{OO}+S_{ON}+S_{NN}=S_{service}.
\]

This partition creates no second service mass, causal charge, event time, or
recursion edge.  A downstream theorem may use one positive submeasure to certify
a genuinely new state, but that theorem must supply the actual physical state and
time it creates.

The input weights here are the complete integrated Moyal edge measure of this
corridor, not an arbitrary positive array.  Their total must realize at least the
certified integrated service lower.  A zero or under-mass edge array cannot be
labelled as the same positive service law.

## 3. The endpoint carrier already contains a hard shell

The full-survivor carrier exists at the earlier endpoint itself.  The transported
annular support satisfies

\[
(3/5)e^{-1/30}A<|\xi|<(3/2)e^{1/30}A,
\]

which lies strictly inside `(A/2,2A)`.  Since the smooth multiplier is a
contraction, split into the two exact hard shells

\[
\{A/2<|\xi|\le A\},\qquad \{A<|\xi|\le2A\}.
\]

If

\[
\mu_A=A\|P_Au\|_2^2,\qquad
\mu_{2A}=2A\|P_{2A}u\|_2^2,
\]

then exactly

\[
A\|Q_Au\|_2^2
\le \mu_A+\frac12\mu_{2A}
\le \frac32\max(\mu_A,\mu_{2A}).
\]

Hence

\[
\boxed{\max(\mu_A,\mu_{2A})\ge\frac23A\|Q_Au\|_2^2.}
\]

This hard-shell witness is at the **same physical endpoint**.  Relative to the
incoming shell `M`, the witness frequencies are `3M/4` and `3M/2`.  This is a
comparable-scale statement, not monotone progress.  If the two actual shell
masses tie, both witnesses are retained; no frequency-order or theorem-order
priority is introduced.

If an endpoint hard shell `H` is used to start the next generic-shell corridor,
the next renewal carrier is `3H/4`.  Thus the typed chain is

`parent shell M -> current carrier A=3M/4 -> endpoint shell H in {A,2A} -> next carrier 3H/4`.

Only the corridor step consumes physical time.  The intermediate registrations
are zero-time witness maps, but their distinct frequency roles must not be
identified.

## 4. Master consequence

A chain

`critical shell -> full natural corridor -> own-scale service -> Moyal/material rereading -> endpoint survivor`

contains one physical recursion edge, not several.  The extra arrows are witness
maps on the same corridor.

Therefore an infinite path whose apparent recurrence consists only of these
full-survivor/service layers is simply an infinite chain of full-natural physical
corridors.  The existing physical-time telescope already proves:

- if its frequencies remain bounded, it reaches `t=0` after finitely many
  corridors;
- if it avoids `t=0`, the unresolved route is genuinely UV-unbounded.

This removes **service theorem depth** from the named-owner recurrence problem.
It does not remove independent source/service/reuse events and does not close the
UV-unbounded survivor alternative.

## 5. Anti-theorems encoded by the candidate

- physical service is not denied or discarded; only duplicate event depth is
  removed;
- a service lower is not an additive globally bounded reset;
- OO/ON/NN does not create another causal law;
- endpoint hard-shell existence does not justify a lexicographic tie break;
- a carrier scale cannot be rebound to an independently supplied shell scale;
- an absolute unit-sized tolerance cannot certify a UV corridor or shell mass;
- arbitrary or zero weights cannot be relabelled as the same positive service;
- comparable endpoint scale is not called directional progress;
- no packet persistence or coherent-cell dominance is introduced.

## 6. Physical PDE falsification lane

The companion audit probe evolves the unforced three-dimensional incompressible
Navier--Stokes Fourier--Galerkin system itself.  It uses Leray projection,
viscosity, the quadratic convective term, 2/3 dealiasing, and RK4 time
integration on a fixed multiscale divergence-free field.  Across several spatial
resolutions it reads, from the same evolved physical interval,

- the direct `Q^2` carrier-energy balance, including nonzero nonlinear carrier
  work;
- the global NS energy balance and zero global work of the dealiased quadratic
  nonlinearity;
- the endpoint two-hard-shell cover and its `2/3` witness lower;
- positive bounded increment service integrated on that same corridor.

This is a resolved Fourier--Galerkin approximation to Navier--Stokes, not an
artificial evolution and not a proof of the continuum theorem.  Its purpose is
falsification: a sign, projection, scale, or duplicate-time error in the proposed
physical bookkeeping must also survive contact with an actually evolved
nonlinear incompressible PDE trajectory.

No Navier--Stokes global-regularity claim is made.
