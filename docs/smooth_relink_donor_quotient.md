# Smooth physical relink donor quotient

## Status

**Independently audited theorem**.  The author's candidate was certified on exact
SHA `8f8cdb2f4ad57bd6f70eafc3043a9bb60ee34d03`.  Adversarial audit then exposed
three distinct fail-open seams in native-scale validation, master owner/mass
replay, typed gauge provenance, and donor closure.  The repaired exact audit SHA
is `d40d6c280973ad860378cad8a0cc078fea81ac1a`; the red-to-green record is in
Section 8.

The candidate status string is

`EXACT_SMOOTH_PHYSICAL_RELINK_DONOR_QUOTIENT__GAUGE_QUOTIENTED_KPHYS_PAIR_FLUX__FINITE_SAME_EVENT_NEGATIVE_DONOR_CLOSURE__SMOOTH_RELINK_ZERO_RECURSION_DEPTH`.

This theorem asks a narrower question than the smooth quadratic-carrier interface theorem.  Once arbitrary observer motion has already been quotiented and the residual skew operator `K_phys` has been proved physical, does positive `K_phys` relink create a new recursive Navier--Stokes generation?

The answer encoded here is **no**.  It is genuine physical redistribution, but it is a finite conservative flux among the smooth `Q^2` roles at the same physical event.

## 1. The smooth and hard measures remain different

For a transported smooth square partition

\[
\sum_a A_a^2=I,
\qquad
\eta_a=A_a^2,
\]

the gauge-quotiented smooth theorem gives the signed residual-skew rows

\[
R_a^{K_{phys}}
=-2\operatorname{Re}\langle \eta_a u,K_{phys}u\rangle.
\]

Define the synthesis-pair matrix

\[
\boxed{
T_{ab}^{phys}
=-2\operatorname{Re}\langle \eta_a u,K_{phys}\eta_bu\rangle.
}
\]

Because `K_phys*=-K_phys`,

\[
T_{ab}^{phys}=-T_{ba}^{phys}.
\]

Because `sum_b eta_b=I`,

\[
\boxed{R_a^{K_{phys}}=\sum_bT_{ab}^{phys}.}
\]

These are smooth synthesis-role identities.  They are **not** obtained by identifying `eta_a` with the orthogonal hard event projectors used by the resolved hard-role donor theorem.  The two physical disintegrations remain distinct.  What they share is only the finite antisymmetric-flux lemma once each has independently proved its own pair law.

## 2. Pair-matrix binding is part of the physical certificate

A list of signed relink rows whose sum happens to vanish is not enough.  The production `GaugeQuotientedInterfaceWork` now carries the exact `T_ab^{phys}` matrix produced by the same `K_phys`, state, and transported square partition that produced the relink rows.

The donor quotient accepts only that typed work certificate and verifies

\[
T+T^T=0,
\qquad
T\mathbf 1=R^{K_{phys}},
\qquad
\sum_aR_a^{K_{phys}}=0.
\]

Arbitrary relink arrays or an unbound antisymmetric matrix fail closed.

## 3. Positive relink is incoming donor flux

Set

\[
F[b\to a]=[T_{ab}^{phys}]_+.
\]

Then exactly

\[
\boxed{
R_a^{K_{phys}}
=\operatorname{inflow}_a-\operatorname{outflow}_a.
}
\]

Hence every role with positive relink work is a recipient of same-event flux.  It did not generate that energy.

For any set of roles `C`, antisymmetry cancels every internal pair and gives the subset divergence identity

\[
\sum_{a\in C}R_a^{K_{phys}}
=F(C^c\to C)-F(C\to C^c).
\]

Internal cycles are circulation and disappear from this balance.

## 4. Finite donor closure

Fix **one** role satisfying

\[
R_a^{K_{phys}}>0
\]

and repeatedly add every role which sends positive flux into its current backward
closure.  Since the role set is finite, this closure terminates.

This recipient-specific closure must contain a role with

\[
R_b^{K_{phys}}<0.
\]

Otherwise every role in the closure would have nonnegative net relink work and
the starting recipient would make their total strictly positive, while by
construction the closure receives no positive flux from outside.  That
contradicts the exact subset divergence identity.

The argument is replayed independently for **every** positive recipient.  It is
not enough to close the union of all recipients and find one donor somewhere in
that aggregate: disconnected components may not borrow donor provenance from one
another.

After quotienting internal cycles, a donor path uses at most

\[
\#\{\text{smooth roles}\}-1
\]

edges.  Every edge remains at the same physical event and physical time.

## 5. Master semantics

The label

`conservative_smooth_role_relink`

is therefore not a recursive generation owner.  It is same-event physical provenance.

The canonical energy-reentry route separates two outputs:

- `owner_bundle`: genuine owners which may create a recursive event;
- `same_event_relays`: conservative provenance at the already existing event.

Thus:

- pure positive smooth relink -> finite donor relay, **no child `RecursiveEventState`**;
- smooth relink tied with existing strain -> relink is quotiented at the same event, while strain remains the recursive owner;
- HH generation, high-strain dissipation, inheritance and other genuine branches retain their existing semantics.

The old compatibility API which insists on returning a `PhysicalOwnerBundle` rejects a pure-relink reentry rather than manufacturing a recursive owner.

## 6. Anti-theorems

This theorem does **not** say physical relink is an observer artifact.  Observer motion was already removed before `K_phys` was defined.  The residual relink is physical; its correct ontology is conservative same-event redistribution.

It also does not identify hard and smooth interface measures, choose a primary donor, invent a scale direction, or turn a negative donor row into a new generation.  Donor roles are provenance inside one interaction.

## 7. Remaining frontier

Smooth `K_phys` relink now joins hard skew circulation as a certified structurally closed recursion seam.  Genuine recurrence still includes resolved strain/deformation, actual HH generation, source/dissipation, independent service, material/new ancestry, and causal reuse/Bellman endpoints.

No Navier--Stokes global-regularity claim is made.

## 8. Independent audit certification

The independent audit deliberately began with tests that the candidate was
expected to reject.  Three exact red gates found three different theorem-boundary
failures:

- exact SHA `5911a4289637833e221b339c3ca87fea21a14e7f`, run `31459124052`:
  `6/6` anti-tests failed.  Unit-scale tolerance floors destroyed covariance at
  tiny native work, substantial relative pair/row defects were accepted, and the
  continuum master trusted forged owner labels and arbitrary claimed mass instead
  of replaying the typed native work split;
- exact SHA `b10a823bfce35b310c1234e79871fd3e35885bc6`, run `31460257649`:
  `6` failed and `10` passed.  The typed gauge-quotiented work certificate
  accepted negative, `NaN`, and infinite gauge/skew provenance residuals;
- exact SHA `0d6a5506300b0d04404f8708723c7799cf591ce2`, run `31460546986`:
  `1` failed and `16` passed.  Aggregate backward closure let a recipient in one
  disconnected component borrow the donor found for another component.

The repaired implementation removes every physical `max(1,...)` tolerance floor,
uses native-relative comparisons, replays the actual positive smooth-interface
split inside the master, binds master mass to the replayed native work, validates
gauge provenance as finite and nonnegative, and constructs an independent
negative-donor trace for every positive relink recipient.  In a mixed
relink/strain event only the positive strain component enters the recursive owner
bundle; conservative relink remains same-event provenance.

Exact independently audited SHA:

`d40d6c280973ad860378cad8a0cc078fea81ac1a`.

Adversarial GitHub Actions run `31460849461` completed successfully:

- `17/17` native-scale, provenance, mass, gauge and per-recipient anti-tests;
- `785/785` theorem tests;
- `200,000` bound smooth relink laws over native pair-work scales
  `[1.176e-141,1.577e+141]`;
- worst native-relative pair antisymmetry and row-binding residuals `0.0`;
- worst native-relative total relink residual `1.593e-15`;
- minimum native-relative incoming-minus-recipient-gain margin `0.0`;
- donor-existence and pair-binding failures `0`;
- maximum sampled shortest donor path `3`;
- `100,000` smooth quadratic-carrier dependencies, `100,000` hard-donor
  dependencies and `100,000` continuum-master dependencies, all green.

The same audit run evolved the unforced three-dimensional incompressible
Fourier--Galerkin Navier--Stokes system with Leray projection, viscosity, `2/3`
dealiasing and RK4.  At every evolved state it applied the actual resolved
linearized operator `L_V f=B(V,f)+B(f,V)` and read the `K_phys/S` work split from
physical Hilbert pairings; no random matrix supplied the operator.  For
resolutions `12,16,20`, all `49/49` snapshots had positive relink and mixed
relink/strain work, maximum divergence was `7.223e-17`, global balance residual
was at most `1.173e-11`, `K` antisymmetry residual was `0`, `K` row residual was at
most `6.628e-15`, `S` row residual was at most `8.865e-15`, donor path length was
`1`, and master failures were `0`.  The final positive-relink resolution spread
was `1.234e-6`.

Stored audit artifact:

`recorded-results/31460849461/audit-smooth-relink-donor-results/`.

GitHub artifact digest:

`sha256:cac760c5ae1388d5f3d7cde20e3105272f72a4201aa542e709f4ff7781b51177`.

The exact-head dedicated workflow `31460849463` was also successful with the
same `785` tests, `50,000` native-scale donor laws, the three-resolution physical
PDE probe, and all direct dependencies.  Its artifact digest is
`sha256:ee914926547fb47bbf11f2b4c37f147efd0aa808f35e11f1c247f6e807c1d395`.

Full physical-energy causal integration run `31460849500` was successful on the
same exact SHA with all `785` theorem tests and the complete causal stack.  Its
artifact digest is
`sha256:e58cdeca217a1d36cc8b10f49aea74ebc6b2858f3fc9fbfa5390e3ed67756e3d`.

These numerical trajectories are falsification evidence, not a continuum proof.
The exact theorem is the finite antisymmetric-flux argument conditional on the
gauge-quotiented work certificate being produced by the stated PDE carrier,
state and operator.

### Author pre-audit certification

Exact certified implementation/final fixture SHA:

`8f8cdb2f4ad57bd6f70eafc3043a9bb60ee34d03`.

Dedicated GitHub Actions run:

`31457786141` — **success**.

It passed:

- `716` theorem tests;
- `50,000` bound smooth `K_phys` relink laws;
- worst pair antisymmetry residual `0.0`;
- worst row-binding residual `0.0`;
- worst total relink residual `1.4210854715202004e-14`;
- minimum incoming-minus-recipient-gain margin `0.0`;
- donor-existence failures `0`;
- maximum sampled shortest donor path `3`;
- pair-binding rejection failures `0`;
- the smooth `Q^2` carrier dependency remained green;
- the hard resolved donor lemma remained green;
- companion `50,000`-state continuum-master smooth-relink recursion-barrier failures `0`.

Stored artifact:

`recorded-results/31457786141/smooth-relink-donor-quotient-results/`.

GitHub artifact digest:

`sha256:be54c2dad2836ef6deca51dc69f21b198a81375fad6f9107001536517e071b5c`.

Full physical-energy causal integration:

`31457786115` — **success** on the same exact SHA.

It passed the same `716`-test suite and the complete causal spine with `57`
successful job steps.  The final master episode stress checked `20,000` traces
with worst margin `0.0`.

Integration artifact digest:

`sha256:2bf74da42e08d053fd5348b60892bc511e68795f0c9d016b61751d6af55d3201`.

The companion smooth quadratic-carrier workflow `31457786119` also completed
successfully on the same exact SHA, certifying the pair-matrix binding in both
directions.

Failure/correction provenance:

- initial implementation SHA `bc0bc248356e7dc18235532649c4a3c1ccdd5a8d` had one brittle certificate-wording assertion after `715` passing tests; no theorem stress ran;
- SHA `113133f31bf3852f16bbc80e692dbce26430967a` changed only that fixture and exposed a second brittle substring assertion in the same certificate test, again after `715` passing tests and before theorem stress;
- final SHA `8f8cdb2...` changed only that second fixture assertion to semantic checks.

No theorem identity, pair law, donor closure, master routing, physical bound, or
workflow topology changed in either correction.
