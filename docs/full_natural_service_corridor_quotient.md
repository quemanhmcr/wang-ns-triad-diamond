# Full-natural service corridor quotient

## Status

Certified on exact implementation SHA `e351d0d6bef5a6bd6275083e1d2e706acf717a18` by dedicated GitHub Actions run `31448743219` and full physical-energy causal integration run `31448763557`.

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
the incoming hard-shell critical-mass input, and endpoint carrier-mass lower.  A downstream reader may verify these data but
cannot provide replacements for them.  In particular, a carrier proved near one
`A` cannot be rebound to a foreign hard-shell scale by passing a new number to an
adapter.

That binding is enforced at the producer boundary as well.  The generic-shell
producer receives the actual parent `M` and rejects any requested `A` other than
`3M/4`; its survivor output also carries the exact parent-shell critical-mass
input so a checkpoint successor cannot invent or rebind that resource.  On the high-strain route, the pushed-forward critical seed carries its
own shell time, parent and child frequencies, scaled lifetime, and actual renewed
critical mass.  The corridor must reuse all of them, including the event time and
the terminal coefficient mass; it cannot manufacture a later event so that a
short corridor appears full-natural.

Every numerical equality guard is relative to the like-dimensional physical
quantity.  No `max(1,...)` absolute floor is permitted.  This matters on the UV
route because `cA^-2` and the corresponding critical masses can be arbitrarily
small under natural rescaling.

The corridor therefore carries `cA^-2` in its native local elapsed-time
coordinate.  Its global endpoint `t-cA^-2` is retained only as a diagnostic
reading: at deep UV, binary floating arithmetic may round that subtraction back
to `t` even though the positive duration itself remains representable.  The
physical-time telescope uses the carried duration, not subtraction of two nearly
equal observer-clock readings.  The floating certificate still rejects scales
whose duration itself falls outside its finite representable range; this
numerical lane is evidence, not the continuum proof.

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

## 5. Anti-theorems encoded by the theorem

- physical service is not denied or discarded; only duplicate event depth is
  removed;
- a service lower is not an additive globally bounded reset;
- OO/ON/NN does not create another causal law;
- endpoint hard-shell existence does not justify a lexicographic tie break;
- a carrier scale cannot be rebound to an independently supplied shell scale;
- an absolute unit-sized tolerance cannot certify a UV corridor or shell mass;
- arbitrary or zero weights cannot be relabelled as the same positive service;
- a non-finite edge total cannot be relabelled as a finite positive service law;
- the endpoint-cover reader must replay the completed-corridor certificate rather
  than trust classification flags in an arbitrary dictionary;
- numerical tie slack is capped at the certificate tolerance and cannot promote
  a zero or subcritical shell into the joint witness set;
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

The heat-service row is evaluated from the exact intrinsic Gaussian displacement
law.  In Fourier variables the probe uses

\[
2\left(1-e^{-|\xi|^2/(2A^2)}\right)
\]

and removes the radius-`3/A` tail with the same certified Gaussian tail bound as
the theorem.  It checks both the full annular heat lower and the retained bounded
lower at every evolved time; no hand-picked list of translation directions is
used as a surrogate.

This is a resolved Fourier--Galerkin approximation to Navier--Stokes, not an
artificial evolution and not a proof of the continuum theorem.  Its purpose is
falsification: a sign, projection, scale, or duplicate-time error in the proposed
physical bookkeeping must also survive contact with an actually evolved
nonlinear incompressible PDE trajectory.

No Navier--Stokes global-regularity claim is made.

## 6. Certification

Exact implementation SHA:

`e351d0d6bef5a6bd6275083e1d2e706acf717a18`.

Dedicated GitHub Actions run:

`31448743219` — **success**.

It passed:

- `684` theorem tests;
- `50,000` corridor/service/material/endpoint-shell states;
- worst natural-time identity residual `4.440892098500626e-16`;
- worst OO/ON/NN same-measure partition residual `1.0658141036401503e-14`;
- minimum sampled two-hard-shell cover margin `2.840630292262874e-08`;
- exact joint hard-shell witness tie retained;
- companion continuum-master service-witness barrier failures `0`.

The stored dedicated artifact is under

`recorded-results/31448743219/full-natural-service-corridor-quotient-results/`.

GitHub artifact digest:

`sha256:47a7296117a1dc966c47d325aecaab6754bb73bb9d9f8d59f5fc5b1dd49b0c5e`.

Full physical-energy causal integration:

`31448763557` — **success** on the same exact SHA.

It passed the same `684`-test suite and the complete source, pressure, material,
critical-shell, full-natural-service quotient, high-tail, hard/smooth interface,
first-stop, causal-reuse, physical-branch and master chain.  The final master
stress checked `20,000` episode traces with worst margin `0.0`.

The initial implementation SHA `4d19bf616c88e634698039c098e2c774be01a669`
had only two certificate-wording assertion failures after `682` passing tests;
no theorem stress ran and no theorem equation, physical bound, or routing changed
in the correction.

Certification is regression evidence for the encoded exact identities and guards;
it is not a proof of global no-escape or Navier--Stokes regularity.

## Candidate refinement: the horizon endpoint is a checkpoint, not an event

The certified theorem above already proves that service/material rereading adds no
second event layer.  The next candidate quotient makes the remaining endpoint
semantics explicit.

On a full no-hit branch, the earlier endpoint `t-cA^-2` was selected because the
analysis chose one natural horizon.  No first-stop condition fired there.  Thus
physical time has genuinely elapsed across the corridor, but the endpoint itself
is an **analysis checkpoint** rather than a causal event vertex.

Likewise the two-shell statement in Section 3 must be read only as same-time state
geometry.  The hard shells at `A` and `2A` are actual observables, but the ratios
`3/4` and `3/2` relative to the incoming shell arise from covering one smooth
carrier.  The candidate API therefore records explicitly

- `analysis_checkpoint_reregistration=True`;
- `directional_scale_progress_supplied=False`;
- `high_tail_supplier_admissible=False`;
- `cover_ascent_interpreted_as_dynamics=False`.

This candidate refinement does not alter the certified lower
`max(mu_A,mu_2A)>=(2/3)A||Q_Au||_2^2`; it narrows what may be inferred from that
lower.  Certification of the checkpoint refinement is pending.

The checkpoint refinement also closes a potential selection loophole in using the
two-shell lower.  Downstream re-registration does not accept `A` or `2A` as a
chosen branch.  It receives the actual pair `(mu_A,mu_2A)` and applies the exact
realization lemma internally; unique maxima are physical facts and exact ties stay
joint.

Nor may a downstream ledger concatenate bare checkpoint records.  The successor
must arrive through a typed transition that carries the actual winning
frequency/mass pair into the next critical-shell producer at the same endpoint
state.  This makes the endpoint cover a genuine PDE-state rereading without
turning an observer-chosen candidate shell into dynamics.
