# PDE-native material/source-service normal form

## Status

**DRAFT — independent PR, not certified, not merged.**

Proposed status string:

`DRAFT_PDE_NATIVE_MATERIAL_SOURCE_SERVICE_NORMAL_FORM__NO_PRIMITIVE_MATERIAL_RELINK_GENERATOR__KPHYS_ZERO_DEPTH__SOURCE_STRAIN_HH_RETAIN_NATIVE_OWNERS__FRESH_SERVICE_RELAYS_TO_GENERIC_SHELL`

This draft attacks the current `RESEARCH_LEDGER.md` frontier directly: a recursive material/new-ancestry event is not allowed to survive merely because a material label changed. The theorem descends to the local Navier--Stokes law and asks which physical term actually changed the state or supplied the service.

The conclusion is deliberately narrower than mixed-owner termination and deliberately stronger than another bookkeeping quotient:

> **`MATERIAL_RELINK` and `NEW_COHERENT_ANCESTRY` are not primitive generators in the PDE-native recurrence normal form.**
>
> A material manifestation may remain recursive only through the independently witnessed physical source, strain/deformation, or actual HH generation which produced it. Same-carrier inherited stock, selected-family boundary data, and conservative `K_phys` relink have zero generation depth. Fresh NN service is a supplier to the already-certified hard-shell first-stop law, not a new material clock.

No Navier--Stokes regularity claim is made.

---

## 1. Start from the local NS equation, not from a material label

On a smooth pre-singular interval,

\[
\partial_t u+\mathbb P(u\cdot\nabla u)=\nu\Delta u,
\qquad \nabla\cdot u=0.
\]

For the event-anchored smooth carrier `w=Qu`, the certified common-slice theorem gives, while the same `Q` and registered probe `psi` are used,

\[
z(t)=z(s)+I_{HH}[s,t]+I_{interface}[s,t].
\]

There is no material membership indicator, no old/new label and no selected-family characteristic function in this coefficient identity. Therefore a change of material label cannot be inserted as a third Duhamel impulse.

A coefficient obstruction is only an interval locator. At the first obstructed interval the carrier reenters its actual `Q^2` energy law.

---

## 2. The complete smooth-interface energy law has only `K_phys` and `S`

Complete the smooth carrier by the already-certified quadratic partition

\[
\sum_a A_a^2=I,
\qquad \eta_a=A_a^2,
\]

transported by the common affine/Kelvin gauge `G`. After quotienting that observer motion, the native interface row is

\[
J_a
=-2\operatorname{Re}\langle \eta_a u,K_{phys}u\rangle
-2\operatorname{Re}\langle \eta_a u,Su\rangle,
\]

with

\[
K_{phys}^*=-K_{phys},
\qquad S^*=S.
\]

Thus there is no third local linearized term waiting to be called “material relink”. The two physical pieces already have exact meanings:

1. `K_phys`: conservative role redistribution;
2. `S`: existing symmetric strain/deformation.

For positive native interface work,

\[
[J]_+\le [J^{rel}]_+ + [J^{str}]_+.
\]

This is a positive cover, not a re-Hahn operation on a coarsened material label. Exact ties remain joint.

---

## 3. `K_phys` is physical but has zero recursive generation depth

The gauge-quotiented pair law is

\[
T^{phys}_{ab}
=-2\operatorname{Re}\langle \eta_a u,K_{phys}\eta_b u\rangle,
\qquad
T^{phys}_{ab}=-T^{phys}_{ba}.
\]

Its row sums reconstruct the signed relink work and

\[
\sum_a\sum_b T^{phys}_{ab}=0.
\]

The certified donor-closure theorem traces every positive recipient through finite same-event inflow to negative-net donor roles. Internal cycles cancel. Hence `K_phys` can change which role carries energy while creating neither net energy, scale progress nor a new physical time.

Therefore

\[
\boxed{
\text{pure }K_{phys}\text{ relink}
\Longrightarrow
\text{same-event provenance, zero recursion depth}.
}
\]

Calling this same flux a new `MATERIAL_RELINK` child would double-count one physical interface event.

---

## 4. Material sidecars do not supply the missing generator

The inherited-stock/material-sidecar theorems have already separated two observer/material operations.

Intrinsic membership rereading only repartitions a fixed positive service law. It changes no total service and creates no physical work.

A selected-family change carries the exact Moyal boundary currency

\[
R_{switch}
=\sum_{C\in S_{old}\triangle S_{new}}E_C.
\]

The same-state anti-theorem permits `R_switch>0` while every cell increment, positive/negative physical work and total state change vanish. Hence `R_switch` cannot be the missing generation law.

So neither sidecar can rescue a naked material owner.

---

## 5. Fresh coherent ancestry is a supplier, not a primitive cause

For actual positive fresh NN coherent service, the refinement-invariant scale theorem pushes the same positive service law to the fixed LP band index. If `F_j` is the selected fresh-band service on a scaled interval of length `c`, then an actual hard shell satisfies

\[
\mu_{hard}\ge \frac{F_j}{6c}.
\]

For the full fresh law `F\ge Y/4`, with `p_max` the maximal band fraction,

\[
\mu_{hard}\ge \frac{p_{max}Y}{24c},
\qquad
\mu_{hard}e^{H_\infty^{scale}}\ge \frac{Y}{24c}.
\]

The theorem itself already declares the next owner to be `generic_critical_shell_first_stop`. Materiality is reread only after the shell produces renewed service. Therefore the logical route is

\[
\boxed{
\text{fresh NN service}
\to
\text{hard critical-shell seed}
\to
\text{native shell first stop},
}
\]

not `fresh NN service -> new primitive MATERIAL_RELINK event`.

The fresh label is provenance on the supplier law; it is not another physical clock.

---

## 6. Source-owned material service keeps the source as owner

The objective source compiler already decomposes actual objective variation into physical source classes and routes them through native PDE suppliers:

- local coherent source -> resolved `D_V` -> critical shell;
- viscosity -> resolved `D_V` -> critical shell;
- SGS -> coherent service -> high-tail / old-pool / fresh-scale shell;
- pressure -> actual positive SGS source or resolved pressure-pair shell.

If one of these source events is the first physical cause, downstream material service does not mint a second owner. Single-charge semantics require

\[
\boxed{
\text{source event owns its generated service until the architecture deliberately reenters the next physical shell/work stop}.
}
\]

The strict draft API therefore refuses an abstract `RESOLVED_SOURCE` recurrence unless the caller also records at least one already-certified PDE-facing supplier kind. This is a fail-closed integration guard, not a new estimate.

---

## 7. Actual HH generation remains actual HH generation

If the `Q^2` physical-energy gate selects high--high generation, the causal weight is the positive physical child-energy work

\[
d\mathcal T_{HH}(t,\alpha)
=2[\operatorname{Re}\langle c,F_{HH,\alpha}\rangle]_+dt.
\]

Material/coherent labels may be attached downstream to this already-fixed physical law, but they do not replace it. Hence a material manifestation caused by actual HH work keeps the existing internal `HH_REGENERATION` owner.

---

## 8. Native normal-form theorem

Under the certified hypotheses of event-anchored hard roles and transported smooth carriers, same-carrier coefficient registration, exact `Q^2` carrier energy, common affine/Kelvin gauge quotient, smooth `K_phys` donor closure, material-sidecar stock decomposition, fresh-service scale reentry, objective source routing, and physical HH energy causality, any proposed recursive material/new-coherent-ancestry manifestation must fall into one of the following native cases:

\[
\boxed{
\begin{array}{rcl}
\text{inherited same-carrier stock} &\to& \text{zero generation depth},\\[2mm]
\text{membership / selected-family sidecar} &\to& \text{zero event depth},\\[2mm]
\text{positive }K_{phys}\text{ relink} &\to& \text{same-event donor provenance},\\[2mm]
\text{symmetric deformation} &\to& \text{existing strain/dissipation owner},\\[2mm]
\text{objective source service} &\to& \text{original source owner + native supplier},\\[2mm]
\text{actual HH-generated service} &\to& \text{HH generation owner},\\[2mm]
\text{fresh NN service} &\to& \text{hard-shell supplier relay}.
\end{array}
}
\]

There is no leftover local NS term whose only physical name is “material relink”. Consequently the strict master-facing corollary is

\[
\boxed{
\text{naked MATERIAL\_RELINK / NEW\_COHERENT\_ANCESTRY}
\quad\text{is unresolved provenance, not a recursive PDE cause.}
}
\]

A caller that cannot descend below such a label must fail closed.

---

## 9. What this removes from mixed recurrence

If certified, this normal form would remove the material/new-ancestry label as an **independent primitive recurrence class**. The surviving recursive physical roots are the already-existing native ones: source, strain/dissipation, actual HH/high-tail work, material reuse with its own erosion law, reuse/Bellman endpoints, and the initial boundary.

This is a real reduction of the mixed-owner topology: an infinite path can no longer alternate through arbitrarily many newly named material vertices unless each such vertex is backed by a different native PDE event.

What it does **not** yet prove is termination of every possible mixed sequence of those remaining native owners. In particular, this draft does not close the separate generic-HH/degenerate Young-Christ seam and does not derive a global mixed source/strain/HH telescope.

That non-overclaim is essential: the purpose of this theorem is to make the remaining recurrence problem smaller and more physical, not to rename it as solved.

---

## 10. Methodological content

The theorem adds no material clock, packet lifetime, reset count or new causal probability. It uses only identities already forced by Navier--Stokes: local energy balance, skew/self-adjoint decomposition of resolved transport, antisymmetric conservative flux, actual positive source/HH work, physical service-to-shell capacity, and exact inherited-stock/material-sidecar separation.

The simplification therefore appears only after the full local structure is exposed. Material recurrence is reduced because the local PDE has too few genuine ways to change a carrier once observer freedom and conservative redistribution are removed.
