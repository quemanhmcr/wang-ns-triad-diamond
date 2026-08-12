# Continuum helical edge-measure registration: signed Navier--Stokes work before Hahn

Status: **certified on exact implementation SHA `6bc9099190048a8796b54fb2f4782314b699bd4b` after independent adversarial and actual-Navier--Stokes audit**.

The certified one-edge theorem identifies the physical variables carried by one
unordered helical parent pair. The continuum question addressed and certified
here is not another single-edge estimate. It is a measure-registration question:

> how does the actual Fourier convolution of Navier--Stokes become a signed
> physical edge law without choosing a parent orientation, without losing
> cancellation, and without replacing physical work by capacity?

This note records that bridge.

## 1. Use the repository's unitary Fourier convention

The canonical convention is

\[
\widehat f(k)
=(2\pi)^{-3/2}\int_{\mathbb R^3}e^{-ik\cdot x}f(x)\,dx.
\]

Therefore a Fourier product/convolution carries

\[
\boxed{C_F=(2\pi)^{-3/2}.}
\]

This factor belongs to the **native physical continuum multiplier**.  It is
important not to confuse that exact normalization with the existing clean
productivity bound `4 R A_3`: because `C_F<1`, the latter is a deliberately
larger safe upper bound.  It remains valid, but it is not the exact physical
multiplier normalization used by the present measure theorem.

## 2. Quotient parent order without choosing an orientation

Fix a child wavevector `z` and write `y=z-x`.  The ordered Fourier convolution
uses the `x` variable.  Parent exchange is the measure-preserving involution

\[
\iota_z(x)=z-x.
\]

Let

\[
\pi_z(x)=\{x,z-x\}
\]

be the unordered-pair quotient.  Define

\[
\boxed{
\lambda_z^{\rm unord}
=\frac12(\pi_z)_\# dx.
}
\]

The fixed set `x=z/2` is Lebesgue-null.  On one orbit use the physical source

\[
S_{\{x,y\}}
=\mathbb P_z\big(u_x\times\omega_y+u_y\times\omega_x\big).
\]

Then

\[
\int S_{\{x,y\}}\,d\lambda_z^{\rm unord}
=
\int \mathbb P_z(u_x\times\omega_{z-x})\,dx.
\]

No lexicographic representative such as `x<y` is selected.  The quotient removes
only observer ordering; the physical two-parent orbit remains intact.

### 2.1 The outer child variable is part of the same Radon quotient

The preceding fixed-child quotient can be written jointly in variables intrinsic
to the physical triad:

\[
z=x+y,
\qquad
r=x-y,
\qquad
x=\frac{z+r}{2},
\qquad
y=\frac{z-r}{2}.
\]

The inverse linear map has absolute determinant

\[
\left|\det\frac{\partial(x,y)}{\partial(z,r)}\right|=\frac18.
\]

At the geometric level parent exchange is `r -> -r`.  On a helicity-resolved
physical edge, however, the helicity sign stays attached to its physical parent
mode, so the full involution is

`(r,s_x,s_y,s_z) -> (-r,s_y,s_x,s_z)`.

Let `q:R^3_r -> R^3_r/{+-1}` denote the geometric quotient.  Equivalently, let
`Q` denote the finite-group quotient of `(r,s_x,s_y)` by the combined involution
above, with counting measure on the parent helicity signs.  Because the physical
orbit integrand is the **sum of both ordered parent terms with their signs
transported with the parent modes**, the joint geometric base measure for child
plus unordered parents is

\[
\boxed{
 d\Lambda_{\rm unord}
 =\frac1{16}\,dz\,d(q_\#dr).
}
\]

For the helicity-resolved edge space the same normalization can be written

\[
 d\Lambda_{\rm edge}
 =\frac1{16}\,dz\,d\Big(Q_\#(dr\otimes \#_{s_x,s_y})\Big)\,d\#_{s_z},
\]

with the orbit-summed edge density understood.  This does **not** quotient away
helicity: it only quotients the analyst's ordering of the two complete parent
mode objects `(wavevector,helicity)`.  Mixed assignments `(+,-)` and `(-,+)`
remain distinct physical assignments unless related by swapping the two complete
parent modes.

Equivalently, for every integrable ordered triad density `f`,

\[
\int dz\int dx\,f(z,x,z-x)
=\frac1{16}\int dz\int dr\,
\left[
 f\!\left(z,\frac{z+r}{2},\frac{z-r}{2}\right)
 +
 f\!\left(z,\frac{z-r}{2},\frac{z+r}{2}\right)
\right].
\]

This is an analytic change-of-variables identity, not a finite-grid ansatz.  The
quotient is proper: `|r|` descends to a continuous quotient radius, so the inverse
image of a compact quotient set is closed and bounded in `R^3` times a finite
helicity set and is therefore compact.  Hence the pushforward of Lebesgue times
finite counting measure is locally finite Radon.  The fixed locus
`r=0,s_x=s_y` lies inside the codimension-three null set `r=0`.

The **base measure being Radon is not yet the same statement as the physical
weights having finite variation**.  That second step comes directly from the
Navier--Stokes energy class, without a UV parent cutoff.  Write

\[
E=\int_{\mathbb R^3}|\widehat u(k)|^2\,dk.
\]

For the two orthogonal helical coefficients,

\[
\sum_{s=\pm1}|a_s(k)|\le \sqrt2\,|\widehat u(k)|.
\]

At fixed child `z`, the parent convolution is bounded by Cauchy--Schwarz, while
the unordered-parent measure contributes its exact factor `1/2`.  Thus for every
bounded Borel child block `B`,

\[
\boxed{
A(B)
\le
4\sqrt2\,C_F E^{3/2}
\left(\int_B |z|^2\,dz\right)^{1/2}.
}
\]

The audited one-edge law gives pointwise `|T_e|<=A_e`, including nonforward
edges, and the global geometric envelope gives `|J_e c_e|<=J_*`.  Therefore

\[
|W|(B)\le A(B),
\qquad
|F|(B)\le J_*A(B).
\]

So `dA`, signed `dW`, and signed `dF` are locally finite physical Radon measures
on bounded child blocks.  This bound certifies measure existence only; it is not
a causal budget, does not tax scale, and does not promote `dA` into a probability
law.

Thus aggregating several physical child modes before Hahn splitting is not a new
law: it is a finite probe of the same joint `z`--unordered-parent Radon measure.

## 3. Helicity is resolved exactly at the event

For every nonzero wavevector,

\[
\widehat u(k)
=\sum_{s=\pm1}a_s(k)h_s(k),
\qquad
\widehat\omega(k)
=\sum_{s=\pm1}s|k|a_s(k)h_s(k).
\]

The event-anchored hard-role theorem already makes the two helical sectors an
exact orthogonal fiber resolution.  Expanding the two parents and the child gives
all eight sectors `(s_x,s_y,s_z)`.

For arbitrary divergence-free complex parent/child Fourier vectors on one triad
fiber, the new registration verifies

\[
\boxed{
T_{\rm vector}
=\sum_{s_x,s_y,s_z}T_{s_xs_ys_z}.
}
\]

The same equality holds after multiplying by the common upper-scale progress
`log_+(|z|/max(|x|,|y|))`.

Helicity sign is physical interaction data and is **not erased** by the quotient.
It travels with the physical wavevector under parent exchange.  The typed edge
identity therefore stores the unordered pair of physical parent mode objects
`(wavevector,helicity)` together with the child mode.  Any canonical ordering used
inside the implementation is only an exact storage key for this finite-group
quotient and is never a physical orientation, owner or causal choice.  Only the
phase convention of the helical basis is gauge.

## 4. Construct the signed physical measure before taking positive parts

On every certified helical edge the one-edge theorem supplies

\[
T_e g_e=A_eJ_ec_e,
\qquad
A_e=4|z|\,|a_xa_ya_z|.
\]

The continuum measures are therefore

\[
\boxed{
 dW=C_FT_e\,d\lambda^{\rm unord},
\qquad
 dA=C_FA_e\,d\lambda^{\rm unord},
\qquad
 dF=g_e\,dW=J_ec_e\,dA.
}
\]

Their roles are different:

- `dW` is **signed physical child-energy work**;
- `dA` is a positive interaction-capacity reference measure;
- `dF` is signed upper-progress work.

`dA` is never promoted to a causal probability law.

Only after exact reconstruction of `dW` do we take its Hahn decomposition.  If
`P_edge` and `N_edge` are the positive and negative edge masses, then

\[
P_{edge}-N_{edge}=W,
\qquad
P_{edge}\ge[W]_+.
\]

With several continuum fibers one also has

\[
P_{edge}
\ge\sum_{fibers}[W_{fiber}]_+
\ge\Big[\sum_{fibers}W_{fiber}\Big]_+.
\]

The first inequality may be strict because helicity channels can cancel inside a
fiber; the second may be strict because distinct fibers can cancel.  The theorem
therefore **does not** assert

\[
[\int T]_+=\int[T]_+.
\]

## 5. Native block deficit is an observable of the same edge law

Normalize the geometric multiplier by the global single-edge maximum:

\[
m_e=J_e/J_*,
\qquad -1\le c_e\le1.
\]

For a measurable block with `A(B)>0`, put

\[
R_B=\frac{F(B)}{J_*A(B)},
\qquad
\epsilon_B=1-R_B.
\]

Then exactly

\[
\boxed{
\epsilon_B
=\mathbb E_{dA/A(B)}\big[(1-m_e)+m_e(1-c_e)\big].
}
\]

This is the measure-level version of the existing finite polarization theorem.
The crucial architectural change is that the arrays are no longer free inputs:
`A_e,m_e,c_e` come from the actual continuum NS edge registration.

The adapter to `coherent_service_or_flat_gate` therefore **does not accept** an
`avg_transfer_deficit` argument.  It computes `epsilon_B` from this physical
ledger and supplies that value to the existing gate.

Thus a large transfer deficit is not a theorem label.  It is the measured loss of
signed upper-progress work relative to the actual available trilinear capacity.

### 5.1 This is not the Young norm-saturation deficit

The quantity `epsilon_B=1-F/(J_*A)` measures geometry/phase signed efficiency
relative to the **actual modal amplitude capacity** carried by the continuum
edge measure.  It does not assert that `A(B)` saturates the separate sharp Young
bound `A_3 ||f||_(3/2)||g||_(3/2)||h||_(3/2)`.  Complex Young/Christ therefore
remains a downstream theorem with its own weighted-cell efficiency and symbol-
freezing premise.  The continuum measure theorem supplies real `A,J,c,T`; it
does not collapse shape near-extremality into Hodge/phase flatness.

## 6. Positive nonforward work remains physical

It is possible to have

\[
T_e>0,
\qquad
|z|\le\max(|x|,|y|).
\]

Then

\[
g_e=J_e=0.
\]

Such an edge remains in the positive Hahn child-work law `dW^+`.  It is **not**
called signed-good and it contributes full multiplier defect to the capacity
reference.  In particular positive generation and forward scale progress remain
separate physical predicates.

## 7. The low-deficit branch supplies the physical good-core law

Fix the certified threshold

\[
\eta_0=10^{-4}
\]

and the signed-good core

\[
G_{\eta_0}=\{m_ec_e>1-\eta_0\}.
\]

Markov gives

\[
\frac{A(G_{\eta_0})}{A(B)}
\ge1-\frac{\epsilon_B}{\eta_0}.
\]

Hence if

\[
\epsilon_B<\frac1{20000},
\]

the core carries more than half of the capacity measure.

On the same edge identity,

\[
\frac{dW}{dA}
=\frac{J_*m_ec_e}{g_e}.
\]

Single-edge stability gives, on the certified good core,

\[
|g_e-\gamma_*|
\le\sqrt{\eta_0}+25\eta_0
=\frac1{80}.
\]

Therefore the normalized positive physical child-work law and normalized
capacity law on the same core have Radon--Nikodym condition number strictly less
than `53/50`; the clean bound is

\[
\boxed{
\frac{50}{53}
\le
\frac{d\mathbb P_{phys,G}}{d\mathbb P_{cap,G}}
\le
\frac{53}{50}.
}
\]

This finally supplies the premise previously assumed by the physical
transfer-defect moat.  Combining with the single-edge defect moment gives

\[
\boxed{
\mathbb E_{phys,G}\mathcal D
\le\frac{106}{25}\epsilon_B.
}
\]

No Duhamel weight and no packet-count law appears.

## 8. One deficit, two native consequences

The coherent service-or-flat transfer threshold lies far inside the preceding
`1/20000` good-core regime.  Therefore the same physical observable

\[
\epsilon_B=1-F/(J_*A)
\]

has a clean dichotomy:

1. if it crosses the coherent threshold, the block has the existing named
   `physical_transfer_cost`;
2. if it stays below the threshold, the same edge law automatically lies in the
   physical good-core regime used by Hodge/strain rigidity and the transfer moat.

This is not a global reset currency.  It is a local efficiency observable of one
actual nonlinear block.

## 9. Scope

The theorem closes the **continuum signed edge-measure registration** seam:

`actual hard HH event`
`-> unordered Fourier parent quotient`
`-> exact helical work disintegration`
`-> signed dW / capacity dA / progress dF`
`-> native transfer deficit`
`-> physical transfer-cost OR good-core physical law`.

It does **not** prove that every generic HH block is low deficit.  High-deficit
transfer, positive nonforward work, high-tail physics, source/strain/relink, and
mixed physical recurrence remain real branches.  No global-regularity claim is
made.


## 10. Certification and adversarial audit

Exact implementation SHA: `6bc9099190048a8796b54fb2f4782314b699bd4b`.

Four exact-SHA gates were required:

- independent audit `31494249932`: **success**, including `50,000` arbitrary
  continuum fibers, `1,000` near-extremal core blocks, a same-Galerkin-grid
  representation audit and a three-child same-orbit Navier--Stokes audit;
- continuum dedicated `31494250081`: **success**, `776` tests plus the dedicated
  continuum/dependency stresses;
- helical dependency `31494402810`: **success**, `776` tests and `50,000` physical
  helical edges;
- full causal integration `31494250037`: **success**, `61` successful steps.

The independent audit was not ceremonial.  Earlier green candidates were rejected
after it exposed unit-scale tolerance floors, forgeable summary provenance,
nonfinite arithmetic, an invalid convergence demand across different Galerkin
truncations, missing nonforward `|T_e|<=A_e` provenance, duplicate helicity-sector
identity, and the final gap between a Radon quotient base and locally finite
weighted physical measures.  The certified source closes each issue fail-closed.

Exact artifacts are stored under `recorded-results/31494249932/`,
`recorded-results/31494250081/`, `recorded-results/31494402810/`, and
`recorded-results/31494250037/`.

## 10. Master-facing routing now certified; signed-cell Young remains

This theorem fixes `dW` before positivity.  Canonical positive-edge routing now performs the Hahn decomposition of this **same** Fourier/helicity edge law and splits actual `mu^+=dW^+` by native efficiency `r=(J/J_*)c` at `eta_0=1e-4`.  The bad support is bound to the same physical block and terminates through the fixed transfer-loss interface; positive nonforward work is retained there.  The good support is only Young-eligible.

Every later hard/coherent/material label must still inherit `mu^+` by a proved positive restriction/pushforward.  The deterministic hard Fourier/helicity map now supplies `pi_#mu^+` as causal mass and `pi_#dW` as signed cell work, with no causal re-Hahn.  Full signed cell work may bind good mass to Young only on fate-pure cells with zero inherited bad-positive mass.  Mixed good/bad cells remain explicitly unresolved rather than being repaired by analyst refinement.  See `canonical_positive_edge_work_routing.md` for the certified theorem and `canonical_positive_edge_work_routing_frontier.md` for the remaining handoff seam.
