# Full-natural horizon checkpoint quotient

## Candidate status

Candidate theorem awaiting exact-SHA GitHub Actions certification.

This theorem separates three objects which were previously too easy to conflate:

1. a **physical corridor** on which Navier--Stokes has actually evolved;
2. a **physical event** at which a named first-stop/owner law occurs;
3. an **analysis checkpoint** at which a theorem chooses to re-register the same evolving solution.

The distinction removes no physical time, work, energy, service, or shell state.
It removes only event depth created by the analyst's choice of natural horizon.

## 1. The actual corridor clock is `A`, not the incoming shell label

Start from an actual critical hard shell

\[
M\|P_Mu(t)\|_2^2\ge\mu_0>0.
\]

The generic theorem registers the smooth carrier at

\[
A=\frac34M.
\]

Its no-hit natural interval is therefore

\[
I=[t-cA^{-2},t].
\]

If a strain/interface/HH first stop occurs inside `I`, that physical stop is the
next event-facing object.  If the interval reaches `t=0`, the initial boundary is
absorbing.  Only in the remaining branch does the whole interval survive with no
physical stop.

On that branch the elapsed time `cA^-2` is real PDE time.  But the earlier endpoint
was selected because the theorem chose one natural horizon.  No physical condition
fired there.  It is therefore an **analysis checkpoint**, not a new causal/event
vertex.

The checkpoint carries `cA^-2` directly as a native local elapsed coordinate.  Its
absolute endpoint `t-cA^-2` is only a diagnostic: for sufficiently deep UV scales,
the two global clock values can be the same floating-point number even though the
PDE interval is strictly positive.  No certificate recovers the local duration by
subtracting those nearly equal clocks, and no dimensional comparison uses an
observer-unit floor such as `max(1,...)`.

The adapter also replays the complete service-corridor certificate and reads
`M`, its input critical mass, `A`, `c`, and the local duration from the producer.
External arguments may verify these quantities but cannot rebind them.  In
particular, the covariance `(M,c) -> (2M,4c)`, which preserves `c/A^2`, cannot be
used to attach one PDE state to a foreign shell.

## 2. Endpoint hard shells are state witnesses, not automatically scale dynamics

The surviving smooth carrier has transported support inside `(A/2,2A)`.  Hence it
may be reread through the two exact hard shells at `A` and `2A`.  Relative to the
incoming hard shell `M`, these candidate ratios are

\[
\frac{A}{M}=\frac34,
\qquad
\frac{2A}{M}=\frac32.
\]

The upper value is larger than one, but that fact alone is not a dynamical theorem.
It comes from a two-shell cover of one smooth carrier at one checkpoint.

The independently certified physical high-tail route begins only at

\[
\frac{M_{tail}}{N}\ge2.
\]

Therefore

\[
\boxed{\frac32<2}
\]

is encoded as a type boundary: the upper cover witness is not admissible as a
high-tail supplier and supplies no directional scale progress.

This does not deny that the hard shell at `3M/2` is an actual state observable.
It says only that the **reason it was exposed** is checkpoint re-registration, not
an already-proved nonlinear UV-generation mechanism.

## 3. State-certified checkpoint chains carry time but zero event depth

A sequence of no-hit checkpoints

\[
t_0>t_1>\cdots
\]

still obeys the exact physical identity

\[
\sum_j(t_j-t_{j+1})=t_0-t_L.
\]

Time contiguity by itself is insufficient.  Between consecutive corridors the
transition must carry

- the two actual endpoint hard-shell critical masses;
- the unique/joint maximizing witness set determined by those masses;
- the exact endpoint time token;
- one selected witness frequency and the same critical mass as the next producer's
  parent-shell input;
- the same fixed scaled lifetime `c`.

A list of bare checkpoints is therefore not a certified chain, even if its global
clock fields happen to be close.  The typed transition replays endpoint
re-registration and rejects an unrelated shell, a losing cover branch, or a
rebound mass.  For an exact tie both witnesses remain in the transition record;
following either state branch does not declare a causal primary.

The quotient then assigns

- physical time drop: actual corridor time;
- event vertices added: `0`;
- causal charges added: `0`;
- hard-shell rereading: analysis witness only.

Thus the native local time edges telescope, while the absolute endpoint-clock
residual remains only a numerical diagnostic.  Bounded-scale checkpoint
continuation still reaches `t=0` by physical time, but it is not a chain of
recursive physical events.

## 4. The UV obstruction survives, in a cleaner ontology

A hypothetical repeated upper-cover re-registration

\[
M_j=M_0(3/2)^j
\]

uses corridor scales `A_j=3M_j/4`, so

\[
\sum_{j\ge0}cA_j^{-2}
=
\frac{c}{(3M_0/4)^2}\frac1{1-(3/2)^{-2}}
<\infty.
\]

The theorem deliberately **does not** convert this into an infinite event path.
After checkpoint quotient it is an event-free UV continuation seam of the PDE.
A subsequent theorem must decide whether such continuation forces actual high-tail
dissipation/work, a genuine first stop, or another native physical mechanism.

This is stricter than inventing a scale tax: no physical owner has yet been proved,
so none is charged.

## 5. Fail-closed master semantics

The canonical master rejects both

`full_natural_analysis_checkpoint`

and the legacy `FULL_NATURAL_SURVIVOR` disposition from `RecursiveEventState` and
physical-owner bundles.  A full-natural service law remains a same-corridor positive
witness as previously certified.

The resulting architecture is

`physical event -> real no-hit corridor -> analysis checkpoint -> re-registration`

until a new physical first stop/owner event actually occurs.

## Scope

This theorem removes natural-horizon segmentation and two-shell cover ascent from
recursive event depth.  It does **not** prove that event-free UV checkpoint
continuation must enter high-tail work, and it does not prove 3D Navier--Stokes
global regularity.

## 6. Endpoint witness selection is performed by the state, not the analyst

The production checkpoint API does not accept a desired endpoint frequency.  It
accepts exactly the two actual hard-shell critical masses

\[
(\mu_A,\mu_{2A})
\]

at the checkpoint.  The exact two-shell realization theorem then determines the
maximizing witness set.  If one mass is larger, only that physical shell is
returned.  If the masses tie, both witnesses remain joint.  There is no API path
for requesting the upper `3M/2` branch merely because an increasing scale path
would be convenient.

Thus even the checkpoint re-registration geometry is state-driven: the observer
chooses the language used to read the endpoint, but not which physical shell won.
The continuation API additionally requires the next no-hit producer to reuse the
winning pair `(H,mu_H)` at exactly the previous checkpoint state.  Knowing merely
that `H` belongs to the candidate set `{A,2A}` is not enough.

## 7. Numerical falsification lane

The checkpoint anti-tests range across native scales rather than a moderate fixed
unit window.  They include deep-UV clocks for which `A^2` overflows while
`cA^-2` remains a positive representable subnormal, and global endpoint clocks
whose subtraction loses the local interval.

The accompanying physical probe evolves the unforced three-dimensional
incompressible Navier--Stokes Fourier--Galerkin system with Leray projection,
viscosity, the quadratic transport term, `2/3` dealiasing, and RK4 at three
resolutions.  It rereads the carrier balance, heat increment, and hard-shell cover
from the same evolved corridor.  This is a strong attempt to falsify the encoded
identities on actual NS dynamics; it is not a continuum proof and supplies no
global-regularity claim.
