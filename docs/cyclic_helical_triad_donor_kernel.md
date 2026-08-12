# Cyclic helical-triad donor/recipient kernel

Status: **candidate theorem block; certification requires dedicated, independent adversarial/actual-NS, and full causal-integration GitHub Actions gates.**

The mixed-fate reserved Young theorem moved the sharp local HH obstruction from terminal bad-positive assistance to the physical structure of canonical negative work.  This note follows that negative work down to the smallest closed Navier--Stokes object on which its destination is determined: one closed Fourier/helicity triad.

The theorem does **not** turn negative work into a budget.  It does not invent a new causal law.  It proves that the already-fixed canonical Hahn-negative work is same-time donor work whose recipient marginal is the already-fixed canonical Hahn-positive work.

---

## 1. Closed triad before any positive split

Write three nonzero closed wavevectors

\[
k_0+k_1+k_2=0
\]

with helicities `s_i in {+1,-1}` and modal amplitudes `a_i` in the repository's deterministic helical basis.  The real field has

\[
\widehat u(-k_i)=\overline{a_i}\,h_{s_i}(-k_i),
\qquad
h_s(-k)=\overline{h_s(k)}.
\]

Root the same physical triad at each of its three mode energies.  Root `i` uses

\[
\text{child}=-k_i,\qquad
\text{parents}=k_j,k_\ell,
\qquad \{i,j,\ell\}=\{0,1,2\}.
\]

The repository Waleffe triple product is cyclic:

\[
g(k_1,k_2,k_0)=g(k_2,k_0,k_1)=g(k_0,k_1,k_2)=:g_\triangle.
\]

Put

\[
R_\triangle
=4\Re\!\left(a_0a_1a_2\overline{g_\triangle}\right).
\]

Then the three **actual child-energy works** returned by the existing physical edge theorem are

\[
\boxed{
\begin{aligned}
T_0&=(s_1|k_1|-s_2|k_2|)R_\triangle,\\
T_1&=(s_2|k_2|-s_0|k_0|)R_\triangle,\\
T_2&=(s_0|k_0|-s_1|k_1|)R_\triangle.
\end{aligned}}
\]

Therefore, before Hahn,

\[
\boxed{T_0+T_1+T_2=0.}
\]

This is the nonlinear energy-conservation law of the single physical helical triad.  It is not a norm inequality and not an averaged conservation slogan: each `T_i` is independently reconstructed by the direct Leray/curl physical-edge registration.

---

## 2. The correct continuum quotient is the full `S3` closed-triad quotient

The current canonical edge law quotients only the order of the two parents.  In ordered parent variables this contributes the factor

\[
\frac12.
\]

For a **closed unrooted triad**, all six permutations of `(k0,k1,k2)` describe the same three-mode physical object.  Away from the coincident-wavevector fixed loci, quotient the full `S3` action:

\[
\boxed{d\Lambda_\triangle=\frac16\,d\Lambda_{\rm ordered\ closed}.}
\]

The fixed loci require two wavevectors to coincide and are codimension three in the regular continuum, hence Lebesgue-null.  A finite-group quotient of a Radon base remains Radon.

Now mark one of the three physical energy roots.  Root marking gives

\[
3\cdot\frac16=\frac12,
\]

exactly the existing parent-swap edge quotient.  Thus, a.e.,

\[
\boxed{
\Lambda_{\rm edge}
=
ho_\#\big(\Lambda_\triangle\otimes \#_{\{0,1,2\}}\big),
}
\]

where `rho` sends a closed triad plus root index to the corresponding unordered parent edge.

In the already-certified sum/relative variables

\[
z=x+y,\qquad r=x-y,
\]

the inverse ordered-parent Jacobian is `1/8`.  Therefore the closed-triad root-chart density is

\[
\frac16\cdot\frac18=\boxed{\frac1{48}},
\]

and the three root marks give

\[
3\cdot\frac1{48}=\boxed{\frac1{16}},
\]

which is exactly the current canonical joint edge density.

One cyclic re-root is

\[
\boxed{
z'=-\frac{z+r}{2},
\qquad
r'=\frac{3z-r}{2}.}
\]

Its scalar-coordinate matrix has determinant one; hence the six-dimensional absolute Jacobian is also one.  Cyclic re-rooting does not distort the physical base measure.

Global reality negation `k_i -> -k_i` is **not** quotient here.  The canonical edge law already counts both Fourier signs in the real-field energy convention.  The theorem proves covariance under global reality negation but does not change that measure convention.

---

## 3. Canonical donor and recipient laws on one triad

Define

\[
P_i=[T_i]_+,\qquad
N_i=[-T_i]_+.
\]

Exact signed conservation gives

\[
\boxed{
Q_\triangle:=\sum_iP_i=\sum_iN_i.}
\]

If `Q_triangle=0`, there is no work and no donor atom.  If `Q_triangle>0`, define

\[
\boxed{
M_{i\to j}
=\frac{N_iP_j}{Q_\triangle}.}
\]

Then exactly

\[
\sum_jM_{i\to j}=N_i,
\qquad
\sum_iM_{i\to j}=P_j.
\]

Thus the donor marginal is the original Hahn-negative triad work and the recipient marginal is the original Hahn-positive triad work.

### Why this coupling is not arbitrary

For a generic transport problem, specifying row and column marginals does not determine a unique positive coupling.  Here there are exactly three physical energy slots and their signed works sum to zero.  Every nonzero sign pattern is therefore one of

- one donor and two recipients;
- two donors and one recipient;
- one donor, one recipient, and one zero-work slot.

So one side of the transport table is always a singleton.  The positive transport with the required row/column work marginals is therefore **unique**.  The formula above is not a maximum-entropy or product-probability choice; it is simply the unique finite transport table.

This uniqueness must not be generalized to arbitrary coherent cells or to transport tables with four or more independent work slots.

---

## 4. Measure-level kernel: `dW-` donor marginal, canonical `dW+` recipient marginal

Restore the same unitary Fourier factor `C_F=(2 pi)^(-3/2)`.  On the full closed-triad quotient define

\[
\boxed{
d\mathcal M_\triangle(i,j)
=C_FM_{i\to j}\,d\Lambda_\triangle.}
\]

By the root-marked quotient identity,

\[
\boxed{
(\operatorname{donor})_\#d\mathcal M_\triangle=dW^-,
\qquad
(\operatorname{recipient})_\#d\mathcal M_\triangle=dW^+.
}
\]

The second marginal is not a new positive law.  It is exactly the already-certified canonical Hahn-positive edge law.  The new object is only a same-time, same-closed-triad donor provenance kernel over that law.

No capacity measure appears in the mathematical kernel.
No fresh Hahn split appears.
No physical time is changed.
No new recursive event is created.

For the **floating numerical certificate only**, the sum of the three native modal
capacities is used as an error scale, exactly as the physical edge registrations
already use native source/capacity scales to judge reconstruction residuals.  It
is never used as a probability, kernel weight, causal mass, or recurrence currency.
In particular the code does not divide conservation errors by realized `Q_triangle`,
which can be arbitrarily small under real phase cancellation.

If the independently registered realized work is at or below the certified native
floating resolution, the numerical API fails closed: it reports the triad as
`numerically_unresolved_transport` and mints **no donor atom**.  This is not a
mathematical zero theorem and does not delete exact continuum work; it only refuses
to infer a Hahn sign from floating roundoff.  The exact analytic kernel above still
applies to every nonzero mathematical triad work.

---

## 5. Restricted negative work pushes to a positive submeasure of canonical `dW+`

Let `D` be any measurable restriction of donor/root edge space.  Restrict the donor side of the kernel and push to recipients:

\[
\nu_D^+(E)
:=\mathcal M_\triangle(D\times E).
\]

Then

\[
\boxed{
\nu_D^+(\text{all recipients})=dW^-(D),
\qquad
0\le \nu_D^+\le dW^+.
}
\]

This is the key interface with the mixed-fate theorem.  A hard-cell term

\[
n_C=(\pi_\#dW^-)(C)
\]

is not an unexplained cancellation currency.  It is a restriction of actual donor work, and that restriction has a mass-preserving positive same-time pushforward into existing canonical positive causes.

The theorem does **not** say that this recipient submeasure is terminal.  Its recipient edges must still be routed by their existing physical fate: geometry-bad recipient work terminates through transfer loss; geometry-good recipient work remains on the existing Young/HH route.

So the theorem removes `dW-` as a candidate **new independent owner ontology**, while preserving every unit of negative physical work and its recipient provenance.

---

## 6. Generic anti-theorem: a positive recipient need not have one donor

Unique donor is false generically.

Take three distinct closed-mode magnitudes with all helicities equal.  For one orientation of the common phase, the cyclic coefficient differences have the sign pattern

\[
(-,-,+),
\]

so the positive recipient receives energy from two physical donors at the same triad/time.

The implementation contains an explicit nondegenerate two-donor counterexample.  Therefore no generic HH API may store one distinguished energy donor by theorem fiat.

The **kernel** is generic.
The **unique donor** statement below is a signed-good corollary only.

---

## 7. Signed-good forward core forces opposite parent helicities

Now take one positive forward canonical edge with child magnitude normalized to one and

\[
r_e=(J_e/J_*)c_e>1-10^{-4}.
\]

The existing signed-good theorem gives

\[
\boxed{
\frac35<x,y<\frac58,
\qquad
u:=\log(y/x)\le\frac1{200}
}
\]

after ordering `x<=y`.

The exact single-edge sign exhaustion says the maximal opposite-parent-helicity sign factor is

\[
P_{opp}=(x+y)(1+y-x),
\]

whereas the maximal same-parent-helicity factor is

\[
P_{same}=(y-x)(1+x+y).
\]

Since

\[
\frac{y-x}{x+y}=\tanh(\nu/2)<\frac1{400},
\]

and `x+y<5/4`,

\[
\frac{P_{same}}{P_{opp}}
<\frac1{400}\frac{9}{4}
=\boxed{\frac9{1600}}.
\]

All other geometric factors are identical and the sign-exhausted envelope is at most `J_*`, hence same-helicity parents would force

\[
\frac{J_e}{J_*}<\frac9{1600},
\]

which contradicts signed-good efficiency `>1-10^-4`.

Therefore

\[
\boxed{s_x=-s_y}
\]

on the signed-good positive forward core.

---

## 8. Interaction parents are not the same ontology as energy donors

Because `s_x=-s_y` and `s_z` is either `+1` or `-1`, exactly one parent helicity equals the child helicity.

For positive signed-good child work, the cyclic three-root formula then forces:

- the parent sharing the child helicity has **negative work** and is the unique energy donor;
- the other interaction parent has **positive work** and is a simultaneous side recipient.

This is a useful ontology split:

\[
\boxed{
\text{two quadratic interaction parents}
\quad\neq\quad
\text{two energy donors}.}
\]

Young/Christ and the HH source still use both interaction parents.  The theorem must never be used to linearize the quadratic Navier--Stokes source merely because one parent is the unique energy donor.

---

## 9. The side recipient carries a rigid amount of actual physical work

Let

- `D` = donor-parent frequency / child frequency;
- `S` = side-parent frequency / child frequency.

Both lie in

\[
\frac35<D,S<\frac58.
\]

After simultaneous helicity reversal if needed, the three cyclic work magnitudes have the exact common factor `R>0`:

\[
W_{child}^+=(D+S)R,
\]

\[
W_{donor}^-=(1+S)R,
\]

\[
W_{side}^+=(1-D)R.
\]

Hence

\[
\frac{W_{side}^+}{W_{child}^+}
=\frac{1-D}{D+S}.
\]

The scale window gives the clean strict bounds

\[
\boxed{
\frac3{10}
<\frac{W_{side}^+}{W_{child}^+}
<\frac13.}
\]

Since

\[
W_{donor}^-=W_{child}^++W_{side}^+,
\]

also

\[
\boxed{
\frac34
<\frac{W_{child}^+}{W_{donor}^-}
<\frac{10}{13},}
\]

and

\[
\boxed{
\frac3{13}
<\frac{W_{side}^+}{W_{donor}^-}
<\frac14.}
\]

These are ratios of **actual signed Navier--Stokes energy work** in one closed triad.  They are not Young deficits and not normalized capacity fractions.

---

## 10. The side recipient is positive nonforward work

The side cyclic edge has as one of its parents the original high-frequency child.  Its own child frequency is the lower side-parent frequency.  Therefore

\[
\frac{N_{side\ child}}{N_{side\ parent,top}}<1.
\]

Consequently its upper-scale progress is exactly zero:

\[
\log_+\frac{N_{side\ child}}{N_{side\ parent,top}}=0,
\]

and hence

\[
\boxed{J_{side}=0.}
\]

But Section 9 proves

\[
W_{side}^+>0.
\]

Thus the side recipient is precisely the positive-nonforward physical work which the canonical router already places in the geometry-bad sublaw and terminates as

`TRANSFER_WORK_LOSS` with `first_time=None`.

This is an existing fate theorem applied to a newly identified cyclic recipient.  No new terminal rule is invented here.

### Critical caution

The side mode still gains real physical energy.  `TRANSFER_WORK_LOSS` is the terminal fate of that positive causal sublaw in the **forward-transfer recursion**, not a statement that the energy vanished from Navier--Stokes.  The side mode may interact again later through ordinary PDE dynamics.

Therefore one must **not** infer

> every signed-good step loses thirty percent of energy forever, so the cascade terminates.

That would turn same-time redistribution into a fake additive/multiplicative reset.

---

## 11. Relation to the mixed-fate reserved handoff

The previous theorem showed that a nondegenerate full-signed Christ-margin failure obeys

\[
g_C<n_C+(\mu^{-1}-1)b_C.
\]

The `b_C` term is already terminal positive bad work.  The present theorem gives the native interpretation of

\[
n_C=(\pi_\#dW^-)(C):
\]

it is donor work whose same-time recipient pushforward is a positive submeasure of canonical `dW+`.

Thus the remaining negative-work obstruction can be routed back into **existing positive physical causes** rather than promoted to a new cancellation owner or scalar budget.

What is not yet proved is a global recurrence telescope after repeatedly following these recipient causes.  In particular, recipient work may itself be geometry-good and continue into the HH/Young machinery.  The present theorem fixes provenance; it does not manufacture termination.

---

## 12. Actual Navier--Stokes falsification plan

The dedicated PDE audit evolves the existing real, divergence-free, `2/3`-dealiased Fourier--Galerkin Navier--Stokes system.  It uses the same cutoff-7 near-extremal physical triad already certified in the mixed-fate audit:

\[
z=(7,6,5),
\qquad
x=(5,0,4),
\qquad
y=(2,6,1),
\]

with

\[
|x|=|y|=\sqrt{41},
\qquad
|z|=\sqrt{110},
\]

and helicities

\[
(s_x,s_y,s_z)=(+,-,+).
\]

At every physical snapshot, the audit reads the **actual evolving helical coefficients** from the Galerkin solution, registers the closed modes

\[
(-z,x,y),
\]

and reconstructs all three cyclic physical edge works independently.

It checks:

1. cyclic Waleffe coupling identity;
2. `T0+T1+T2=0` before Hahn;
3. donor and recipient measure marginals;
4. signed-good unique donor and side-recipient ratios at the engineered initial event;
5. side positive nonforward geometry;
6. Navier--Stokes energy balance and incompressibility;
7. representation invariance of the **same cutoff-7 Galerkin system** embedded on different FFT grids, measured against the same native triad work envelope rather than realized positive work;
8. amplitude adversaries, which must preserve dimensionless routing and scale every work cubically.

Different Galerkin truncations are not required to agree.  Only different FFT representations of the same finite physical system are an invariant.

---

## 13. Adversarial obligations

Certification should actively try to break the theorem with:

- all eight helicity assignments;
- one-donor and two-donor generic sign patterns;
- exact zero-work amplitudes;
- arbitrarily small phase-cancelled work, which must not create floating donor provenance merely because a rounded sign is nonzero;
- all six permutations of the closed modes;
- spatial translations of modal phase;
- global reality negation without quotienting it;
- uniform wavevector dilation;
- cubic amplitude scaling;
- nonfinite/provenance rejection;
- restricted donor sublaws whose recipient pushforward must remain below full canonical `dW+` rootwise;
- actual evolved Galerkin NS;
- identical-cutoff FFT-grid representation changes.

No sample count or tolerance should be weakened because the cyclic registration is more expensive than one edge.  If performance becomes a bottleneck, immutable closed-triad geometry may be cached only after differential equivalence with the readable three-root reference path is established.

---

## 14. What this theorem would close, and what it would not

If certified, this theorem closes the **native provenance of canonical negative HH work**:

\[
\boxed{
\text{canonical }dW^-
\xrightarrow{\text{same closed triad, same time}}
\text{positive submeasure of canonical }dW^+.
}
\]

On signed-good forward edges it additionally identifies a unique energy donor and a rigid positive nonforward side recipient carrying between `3/10` and `1/3` of the selected child work.

It would **not**:

- replace `dW+` as the causal law;
- make `dW-` a causal probability;
- turn side work into dissipation or an additive reset;
- reduce two HH interaction parents to one source parent;
- prove every positive recipient has one donor;
- create a general coherent/material positive kernel;
- prove every geometry-good cell Young-good;
- terminate generic HH or mixed genuine-owner recurrence;
- close the initial-data or singular-time interfaces;
- prove 3D Navier--Stokes global regularity.

The next global question after this provenance theorem would be whether the donor/recipient law, the already-terminal nonforward side branch, and the existing native recurrence laws can be assembled into a genuine typed telescope for **mixed owner recurrence**, without turning same-time redistribution into a scalar Bellman currency.
