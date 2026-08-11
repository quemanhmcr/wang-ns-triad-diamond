# Smooth physical relink donor quotient

## Status

**Certified theorem** on exact implementation SHA `8f8cdb2f4ad57bd6f70eafc3043a9bb60ee34d03` by dedicated GitHub Actions run `31457786141` and full physical-energy causal integration run `31457786115`.

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

Start with all roles satisfying

\[
R_a^{K_{phys}}>0
\]

and repeatedly add every role which sends positive flux into the current set.  Since the role set is finite, this closure terminates.

It must contain a role with

\[
R_b^{K_{phys}}<0.
\]

Otherwise the closed set would have strictly positive total net relink work while, by construction, receiving no positive flux from outside.  That contradicts the exact subset divergence identity.

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

## 8. Certification

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
