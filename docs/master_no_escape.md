# Master finite-dimensional no-escape theorem

This theorem glues the previously proved finite-dimensional modules. It is an
abstract closure theorem for the grain-cascade model, not a Navier--Stokes
regularity proof.

## 1. Cost-or-episode architecture

For each block `j=0,...,L-1` let `R_j` be its normalized transfer efficiency.
After transfer-adapted extraction suppose

\[
R_j\le e^{-C_j}+\eta_j,
\qquad C_j\ge0,
\]

and define the exact cross penalty

\[
\xi_j=\log(1+\eta_j e^{C_j}).
\]

Then

\[
-\log R_j\ge C_j-\xi_j.
\]

The Bellman, ancestry, resistance, Hodge, and balanced-spherical modules
supply a fixed threshold `c0>0` so that every **costly** block satisfies

\[
C_j\ge c_0.
\]

The only low-cost blocks are synchronized, flat, hemispherical blocks on a
dominant directional core.

## 2. Barycentric entropy-or-reset lemma

Let

\[
\mu=\sum_iw_i\delta_{x_i},
\qquad b=\int x\,d\mu(x),
\]

and fix `beta in (0,1)`.  The exact barycenter--collision inequality from the
spherical module is

\[
H_2(\mu)\ge\log\frac{2}{1+|b|}.
\]

Hence there are two cases.

If `|b|<=beta`, the block pays the fixed atomic entropy cost

\[
\boxed{H_2(\mu)\ge h_\beta:=\log\frac{2}{1+\beta}>0.}
\]

If `|b|>beta`, define the barycentric potential

\[
P(\mu):=-\log|b(\mu)|.
\]

Then it has the uniform reset bound

\[
\boxed{0\le P(\mu)<P_{max}:=-\log\beta.}
\]

No mass is discarded and no new cross-error is created by this reset.

For exact equal-marginal flat propagation the spherical midpoint identity gives

\[
b_{j+1}=b_j/c_*,
\qquad c_*=\cos(\theta_*/2),
\]

and therefore

\[
P_{j+1}=P_j-\kappa_*,
\qquad \kappa_*=-\log c_*.
\]

The nonsymmetric companion identity

\[
b(\nu)=2c_*b(\lambda)-b(\mu)
\]

is used as the gate into this episode: if old/companion marginals are not
synchronized well enough for the barycentric erosion inequality, then either
the companion barycenter is small (hence it pays collision entropy) or it is
concentrated and is sent to the fresh/reuse ancestry classification.  Thus the
only low-cost transition is precisely the synchronized flat transition for
which the potential erosion below is valid.

A particularly natural threshold is

\[
\beta=c_*.
\]

Then

\[
P_{max}=\kappa_*,
\qquad
h_\beta=\log\frac{2}{1+c_*}\approx0.09516.
\]

In the exact model this means a zero-cost synchronized flat episode contains
at most one transition before another entropy/fresh/reuse event must occur.

## 3. Erosion inside one flat episode

For an exact optimal midpoint step with no fresh direction,

\[
P_{j+1}\le P_j-\kappa_*,
\qquad
\kappa_*=-\log\cos(\theta_*/2)>0.
\]

For a near-extremal block allow a nonnegative perturbation `zeta_j` and a
fixed `kappa0>0`:

\[
\boxed{P_{j+1}\le P_j-\kappa_0+\zeta_j.}
\]

If a flat episode contains `m` low-cost blocks, summing while `P>=0` gives

\[
m\kappa_0\le P_{max}+\sum_{episode}\zeta_j.
\]

## 4. Master no-escape theorem with potential resets

Let `N_K` be the number of costly blocks and `N_F` the number of low-cost flat
blocks.  Costly blocks separate the flat cascade into at most `N_K+1`
episodes. Hence

\[
N_F\kappa_0
\le
(N_K+1)P_{max}+Z,
\qquad Z:=\sum_j\zeta_j.
\]

Since `L=N_K+N_F`,

\[
\boxed{
N_K\ge
\frac{\kappa_0L-P_{max}-Z}{\kappa_0+P_{max}}.
}
\]

Therefore

\[
\sum_jC_j
\ge
\frac{c_0}{\kappa_0+P_{max}}
(\kappa_0L-P_{max}-Z).
\]

Combining with the logarithmic cross penalties gives

\[
\boxed{
-\log\prod_{j<L}R_j
\ge
c_{eff}L
-\frac{c_0(P_{max}+Z)}{\kappa_0+P_{max}}
-\Xi,
}
\]

where

\[
\boxed{
c_{eff}=\frac{c_0\kappa_0}{\kappa_0+P_{max}}>0,
\qquad
\Xi=\sum_j\xi_j.
}
\]

Thus, if the perturbation sum `Z` and cross-error penalty `Xi` remain bounded,

\[
\boxed{
\prod_{j<L}R_j\lesssim e^{-c_{eff}L}.
}
\]

This is the corrected closure theorem.  The potential may reset after every
costly block, but every reset creates a new episode and the number of episodes
is itself controlled by the number of costly blocks.

## 5. Where each costly block comes from

The previously proved finite-dimensional modules supply the alternatives:

1. transfer-weighted component splitting -> Bellman collision entropy;
2. hidden atomic entropy -> component entropy or same-ancestry pair cycles;
3. high resistance -> Poisson Bellman cut or conductance collision entropy;
4. low resistance plus curved reuse -> Hodge cost;
5. low resistance plus low Hodge -> integer gauge synchronization;
6. synchronized flat but barycentrically balanced -> atomic collision entropy;
7. synchronized flat and concentrated -> a dominant cap; nested extraction
   either splits it (Bellman cost) or produces a trackable fresh/reused grain.

A fresh grain gives a fresh ancestry split; a reused grain returns to 2--5.
The only branch not immediately paying a cost is therefore the cap-confined,
synchronized-flat, no-fresh transition treated by the episode erosion above.

The single-edge certificate removes the previous placeholder in branch 4.  In
the abstract capacity weights, if the good-edge Hodge threshold is
`E_H>=h_H`, then

\[
-\log R_{block}\ge\frac12E_H\ge\frac{h_H}{2}.
\]

The smooth-flux cocycle now supplies a PDE-facing version on a packet block.
For signed efficiency deficit `epsilon=1-R`, take the pointwise-good threshold
`eta_0=10^{-4}`.  If `epsilon>=1/20000`, the block already pays that fixed
cost.  Otherwise at least half the capacity mass is on the positive core; Arb
certifies that normalized physical child-transfer weights and capacity weights
differ by at most a factor `53/50`.  Therefore

\[
\boxed{
epsilon\ge\frac{25}{106}\mathcal E_H^{phys}
}
\]

on the low-cost core.  A physical Hodge threshold `E_H^{phys}>=h_H` thus pays

\[
\boxed{
c_{0,H}^{phys}=\min\{1/20000,(25/106)h_H\}>0.
}
\]

Here `E_H^{phys}` is weighted by actual positive child-energy transfer.  A
common spectral moat no longer has to be assumed for a signed-good crossing
core.  The four-bin crossing theorem extracts a subcore carrying at least one
quarter of the transfer-weighted Hodge numerator and Arb certifies its shell
halfwidth `<2/25` with smooth moat margin `>1/25`.  If the PDE-facing argument
insists on using one such physical smooth-cut subcore, the conservative branch
constant is

\[
\boxed{
c_{0,H}^{cross}=\min\{1/20000,(25/424)h_H\}>0.
}
\]

A second physical reuse obstruction is available in **helical phase**.  For a
four-edge reuse diamond let `H_phi` be the gauge-invariant signed phase holonomy
after the six modal phases cancel.  The sharp inequality

\[
\sum_{i=1}^4(1-\cos\delta_i)
\ge4\left(1-\cos\frac{|H_\phi|}{4}\right)
\]

feeds directly into the existing polarization deficit.  At `eta_0=10^-4`, if
`|H_phi|>=1/5` and every diamond edge has capacity weight at least `beta_phi`,
Arb certifies the finite-packet cost

\[
\boxed{c_{0,\phi}=\beta_\phi/250.}
\]

This branch is conditional only on the PDE packet extraction actually producing
such a weighted diamond; the phase identity and coefficient are theorem-level.

If instead a set of outside-neighborhood edges has transfer weight at least
`beta_bad`, it pays at least `beta_bad/100`.  Thus, after fixing the positive
classification thresholds used by the other branches, one may take

\[
\boxed{
c_0=\min\{c_{Bell},c_{anc},c_{res},h_H/2,\beta_{bad}/100,
\beta_\phi/250,h_\beta,c_{fresh},\ldots\}>0.
}
\]

This formula does not claim that the threshold choices have been optimized; it
records that the Hodge/single-edge entry is now a theorem-level positive number
rather than a numerical candidate.

## 6. Remaining analytic bridge

The finite-dimensional single-edge and physical scale-flux inputs are now much
more complete than in the earlier version of this note:

- the smooth transfer-weighted midgap theorem realizes the progress sum as an
  exact smooth SGS tail flux on a common spectral moat, so there is no longer a
  sharp-to-smooth comparison loss;
- phase/backscatter cancellation has an exact positive polarization deficit;
- on the signed-good core, actual child-transfer weights and capacity weights
  differ by less than `53/50` at the certified threshold `eta_0=10^{-4}`;
- pressure can be absorbed into the combined-work ledger after spatial
  localization, with the remaining cancellation branch forcing a critical
  annular charge.

The affine three-dimensional polarization side is also no longer featureless.
A triad-normal helical gauge is exactly `SO(3)` covariant, its Waleffe coupling
phase is constant away from zeros, and reuse between different triad planes is
measured by the spin-one dihedral transition `exp(-is psi)`.  Diamond reuse has
an exact gauge-invariant phase holonomy and the positive cost above.  In an
objective transverse frame, first-order trace-free strain is an opposite-
helicity mixer rather than a Berry-phase drift.  The absolute mixers give an
Arb-certified `13/20` tomography statement, but the signed extremal tensor shows
that a common determinant-one deformation of both parent helicity spinors
preserves their **unnormalized symplectic nonlinear wedge** exactly.  After
quotienting this common `SL(2)` gauge, the transfer-distinguishable strain
observable obeys the sharper physically aligned theorem

\[
\boxed{
\|D_\Pi\|_F^2+\|D_1-D_2\|_F^2+\|D_{child}\|_F^2
\ge\frac12\|S\|_F^2.
}
\]

Thus absolute helicity conversion is not promoted by itself to a phase cost;
relative parent polarization, child polarization, capacity normalization and
scalar shape are the variables the nonlinear transfer actually distinguishes.
Moreover, on the signed-good core Arb gives `sin(theta)>9/10`, and the
triad-normal gauge satisfies the scale-free local estimate

\[
\|\dot h_s\|\le\frac52
(\|\dot{\widehat k}\|+\|\dot{\widehat p}\|)
\]

for a parent role.  Hence the global Chern topology does not introduce a local chart singularity
inside an efficient packet.  The transfer-facing common-parent time ordering is
also closed exactly by the symplectic wedge identity; individual-spinor Magnus
phases are not the quantity charged by the nonlinear transfer.  The remaining
PDE burden is localized **differential** transport/SGS/window forcing and its
registration in ancestry.

The transfer-facing time dependence can now be written without a Magnus
truncation.  After scalar attenuation is factored out,

\[
\frac d{dt}(U^TJV)=U^TJ(D_1-D_2)V,
\]

and the child factor adds only `D_child`.  A common arbitrary time-ordered
`SL(2)` parent history cancels pointwise.  With packet forcing, the additional
term is explicitly multilinear in `F_1,F_2,F_3` and has a direct capacity norm
bound.  On a low-strain moving block, Kelvin direction stability plus the
certified helical frame moat gives

\[
\boxed{E_{pol}\le3h+\frac{15}{2}c\,\kappa M.}
\]

Hence the polarization bridge does not introduce a new localization scale: it
renormalizes the existing curvature coefficient and adds a summable frequency-
cell `O(h)` term.

The residual forcing has now been further quotiented by the affine Gaussian
symmetries.  If `x-X=Lz` and `H=nabla^2 U`, define

\[
B=L^{-1}H[L,L].
\]

After allowing center, carrier, covariance and chirp to osculate the flow, the
first scalar transverse non-affine forcing is exactly

\[
\boxed{
\|F_\perp\|_2^2/\|\psi\|_2^2=(3/8)\|\operatorname{Sym}B\|_F^2.
}
\]

The kernel `Sym B=0` is not declared free.  Under incompressibility it is the
quadratic swirl `V(z)=z cross (M z)` with `M` symmetric trace free; it is routed
to spatial strain/relative-polarization variation.  Common material advection,
including common velocity-Hessian chirp, preserves the signed triad phase lock
exactly, so only differential role velocities/sources contribute phase forcing.
Bulk viscosity is Gaussian tangent as well; viscous boundary flux is not.

The affine Gaussian also carries a physical critical grain without imposing an
isotropic cell.  With `r_g=(det Sigma_x)^(1/6)`, the certified radius-two
ellipsoid satisfies

\[
\boxed{r_g^{-1}\int_{E_2}|u|^2\ge3/10.}
\]

A fresh family therefore obeys the physical energy budget

\[
\boxed{\sum_j r_{g,j}\le P\|u(t)\|_2^2/\eta}
\]

at overlap multiplicity `P`.  In the affine incompressible Gaussian model,
`r_g` is exactly preserved by inviscid strain and satisfies
`d(r_g^2)/dt>=nu` with viscosity.  Static affine anisotropy is an exact Young
symmetry and is **not** a Bellman cost.  Reused elongated grains remain a
spacetime ancestry/dynamic-curvature problem.

The spatial window can be made affine as well.  For
`chi_(L,M)(x)=chi_0(L^-1(x-X)/M)` transported by `Xdot=U(X)`, `Ldot=A(X)L`,

\[
D_t[L^{-1}(x-X)]
=L^{-1}[U(X+Lz)-U(X)-A(X)Lz].
\]

Hence affine motion cancels exactly, while the material remainder is
`O(kappa_aff M)` with `kappa_aff=||L^-1(nabla^2U)[L,L]||`.  The shell lower axis
gives `N^-1||grad chi||<=3C_chi/(2M)`, and the physical convolution commutator
satisfies

\[
\|[\chi,G_N*]f\|_2
\le (3/2)m_1(G)C_\chi M^{-1}\|f\|_2.
\]

Thus the ellipsoidal geometry preserves the same `1/M` versus
`kappa_aff M` balance without an aspect penalty.  The remaining PDE work is to
insert pressure, `RU`, viscous boundary and partition-overlap terms into this
affine window rather than to invent a new localization scale.

The actual Navier--Stokes gap is therefore now the **construction** of the
transfer roles and their spacetime packetization.  Common parent/child log
shells are no longer a separate hypothesis: crossing plus the certified
single-edge gap produces them by four-bin extraction.  Nor is a full iterative
within-block profile decomposition required: on an efficient frozen block,
Christ's inverse Young theorem gives one Gaussian profile and the exact
trilinear replacement loss is `3 eps_G+3 eps_G^2+eps_G^3`.  Frequency-block
support then gives the critical-mass bound

\[
N\|G\|_2^2\ge C_\Omega^{-1/3}(1-\epsilon_G)^2.
\]

The frequency multiplier part of that realization is now explicit.  On the
certified compact crossing shells, freeze the normalized smooth SGS symbol on
relative cells of diameter `h`.  If `L_M` is its dimensionless Lipschitz
constant, sharp Young gives

\[
|T_M-T_{M_h}|\le A_3L_Mh\prod_{j=1}^3\|f_j\|_{3/2}.
\]

Choosing the relative cell sizes so that `sum_j L_{M,j}h_j<infinity` puts this
loss directly into the summable cross-error ledger.  Moreover at one-percent
Christ profile distance, the certified `2/25` shell gives

\[
\boxed{N\|G\|_2^2>3/4.}
\]

What remains is therefore spatial/time rather than multiplier-theoretic: attach
the frozen frequency roles to moving physical packets, propagate the affine
Gaussian grain through its packet lifetime, control partition overlap and
non-pressure leakage, and synchronize the moving windows across generations.

## 7. PDE-facing spacetime grain dynamics

The frozen-time packet bridge now has a quantitative affine/Kelvin continuation.
For the trace-free symmetric strain `D` seen by the extremal triad plane,
Arb certifies the infinitesimal coercivity

\[
\frac12\dot u^2+2\dot v^2
\ge\frac{43}{100}\|D\|_F^2.
\]

A frozen principal strain with eigenvalues `+/-d` and `dT<=1/25` pays

\[
\boxed{
\frac1T\int_0^T\operatorname{Def}(t)dt
\ge\frac1{10}(dT)^2.
}
\]

After removing local rigid rotation, the same conclusion is stable under five-
percent strain variation: if `sup||D(t)-D(0)||<=d/20` and `dT<=1/30`,

\[
\boxed{
\frac1T\int_0^T\operatorname{Def}(t)dt
\ge\frac1{24}(dT)^2.
}
\]

Failure of that coherence forces objective strain variation `>d/20`.  The exact
Navier--Stokes gradient equation routes the variation into

\[
\mathring S
=Q(A)-\nabla^2p+\nu\Delta S,
\qquad
Q(A)=-S^2-\Omega^2+[S,\Omega].
\]

For `d=sigma N^2` and a packet lifetime `T=cN^{-2}`, at least one source channel
then has normalized average size `>=sigma/(60c)`.  In the band-limited packet
model the quadratic and viscous channels force critical mass by Bernstein; the
far pressure-Hessian kernel has the stronger locality gain `5-3=2` after 3D
packing.

The moving window should deform with the affine flow.  Its affine advective
leakage vanishes exactly; velocity curvature competes with the smooth-filter
commutator.  The combined spatial error has the exact form

\[
\boxed{
E(M)=\frac aM+b\kappa M,
\qquad
\kappa=N^{-3}\|\nabla^2u\|_\infty,
}
\]

with optimum

\[
\boxed{
M_*=\sqrt{a/(b\kappa)},
\qquad
E_*=2\sqrt{ab\kappa}.
}
\]

This corrects the earlier heuristic that one fixed expanding spatial moat could
make every localization term summable.  The quadratic schedule remains valid in
the defect-space nested-grain theorem and for the commutator alone, but the
physical moving window must be curvature-adaptive.

These packet-lifetime inequalities are not yet inserted as new universal master
`c_0` branches.  The intrinsic Gram identity removes common triad-plane tilt as
a scalar multiplier defect.  What remains is the full 3D packet-frame
realization with helical polarization/phase transport, adaptive-window
compatibility, and registration of the forced critical mass as fresh or reused.

## 8. PDE-facing localized pressure alternative

For a smooth physical SGS packet block, spatial localization should use the
combined work

\[
G=\Pi+\nabla\cdot(PU).
\]

The resolved local-energy identity then cancels pressure exactly from the
window leakage. Along a chain with nonincreasing weights,

\[
\sum_j a_j W_{j,+}+\sum_j a_jD_j
\le a_0E_0+\sum_j a_j(L_j)_++\sum_j a_jW_{j,-}.
\]

If the raw SGS transfer `S>=0` at a packet window is not retained as combined
work at level `S/2`, pressure cancellation forces

\[
\int_{\operatorname{supp}\nabla\chi}
(|U|^3+|P|^{3/2})
\ge\frac{S}{2\|\nabla\chi\|_\infty}.
\]

Thus the PDE-facing costly/escape classification acquires a concrete new
channel:

\[
\boxed{
\text{forward physical Hodge/SGS transfer}
\to
\text{combined-work depletion}
\quad\text{or}\quad
\text{critical annular charge}.
}
\]

The pressure channel has a further packet-level closure.  Since
`int U dot grad chi=0`, the constant pressure mode cancels.  The far pressure
kernel is therefore a dipole difference with `|grad K|=O(|x-y|^{-4})`.
Three-dimensional packet packing grows as `2^{3n}`, leaving the summable spare
power `4-3=1`.  The local pressure source is controlled by
Calderon--Zygmund and band-limited Bernstein, so a definite local cancellation
forces definite critical `L^2` mass.  In the stated packet model this yields

\[
W_{cancel}\ge\rho
\Longrightarrow
\mu_{fresh}\ge
\min\left\{(\rho/(2C_{near}))^{2/3},\rho/(2C_{far})\right\}>0.
\]

The remaining conditional step is to construct the actual PDE packet frame so
that its near/far pressure split, packing constants and cross errors satisfy
these hypotheses, and then to register the resulting mass in the spacetime
fresh/reuse ancestry ledger.
## Actual packet forcing split

The PDE bridge now uses a strict two-level bookkeeping rule.  A selected
full-velocity role lies above `3N/5`; a strict transporter below `N/4` cannot
produce it by a low--low interaction.  The affine low--high equation is exactly
the Kelvin equation, and a shell-localized divergence-free weak packet removes
microscopic pressure while paying only the existing `O(1/M)` commutator.

After quotienting common affine/Gaussian/helical motion, the field residual is
resolved by Hermite degree.  Only the degree-zero projection is the `F_i` of the
forced symplectic identity.  Parent/profile degree-zero remainder work is the
existing trilinear cross-error `Xi`; `H_1` polarization and `H_3` envelope
curvature are orthogonal sidebands.  Their raw curvature has an affine
connection

\[
\dot B+2A_{aff}B=L^{-1}\operatorname{Hess}(D_t^V V)[L,L],
\]

so sideband dephasing is sourced by pressure-third, differentiated-SGS or
viscous-fourth-derivative channels rather than a new abstract phase error.

At the macroscopic resolved-energy level, differential advection plus `RU`
forces a critical cubic annular charge; `R` is itself bounded by cubic velocity
increments at the filter scale, resolved viscous boundary transport is `O(1/M)`,
and a complete quadratic spatial partition cancels globally.  These terms are
not inserted again into microscopic spinor forcing.

What remains conditional for the master theorem is now the **spacetime
registration**: the `H_3` sideband and the mild-aspect physical `H_1/swirl` branch both have local no-escape theorems and pair-rescue
ancestry routing; persistent
cubic increment charge must enter fresh/reuse, dissipation or Bellman
bookkeeping; and selected-subfamily interfaces must stay summable in `Xi`.


## Hermite sideband no-escape insertion

The actual packet equation now supplies a transfer-facing local alternative for
the scalar `H3` curvature channel.  In the affine interaction frame define

\[
I_3=\int_0^T\|\operatorname{Sym}\widetilde B\|dt.
\]

The local theorem says that a block with non-negligible `I_3` must enter one of
five branches:

\[
\boxed{
\text{source/dephasing},\quad
\text{nonlinear sideband feedback},\quad
\text{large daughter},\quad
\operatorname{Def}\ge {3\over4096}I_3^2,\quad
R_{pair}\ge {3\over4096}I_3^2.
}
\]

The pair branch is not left as a perturbation.  Hermite parity makes pair rescue
an ordinary weighted graph on odd daughter atoms.  A transfer-weighted ancestry
split routes it into

\[
\boxed{
\Xi,\quad
\text{dominant reused daughter},\quad
\text{component Bellman entropy},\quad
\text{same-ancestry cycle attachments}.
}
\]

Thus, once a PDE block supplies a uniform lower threshold `I_3>=iota_3`, the
base-edge branch has the explicit local cost

\[
\boxed{c_{H3}(\iota_3)={3\over4096}\iota_3^2.}
\]

The other branches already have named destinations in the master architecture;
what is not yet proved is a continuum statement guaranteeing a uniform
`iota_3` or pricing the physical `H1/swirl` alternative on every efficient
block.  In particular the intrinsic affine curvature bound

\[
\|\operatorname{Sym}B\|^2+\|C_{hook}\|^2\ge {1\over6}\|B\|^2
\]

must not be converted into an aspect-independent Euclidean polarization cost.


## Physical mild-aspect H1/swirl insertion

The hook sector has a transfer-facing pointwise bridge with no scalar `D_Pi`:

\[
Q_{pol}=\sum_c(\|D_1-D_2\|_F^2+\|D_3\|_F^2),
\qquad
Q_{pol}\ge\frac1{25}\|B_{hook}\|^2
\]

on `cond(L)<=21/20`.  The relative-coordinate action energy is `Q_pol/2`, but
that coordinate evolves by nonunitary `SL(2)` and is not itself one Young role.
For the three physical packet roles,

\[
\boxed{
\sum_i\|F_i^{H1}\|^2
\ge\frac14Q_{pol}
\ge\frac1{100}\|B_{hook}\|^2.
}
\]

On the existing low-strain lifetime branch, each base role has

\[
K=\int_0^T\|G_i\|_{op}dt\le1/30.
\]

Conditioning of both the interaction pullback and physical pushforward gives

\[
\boxed{
J_1\ge I_1/(11T)
\quad\text{or}\quad
\delta_1^2\ge I_1^2/480.
}
\]

After physical feedback, one-of-three-role selection, odd-Hermite convexity and
pair-rescue splitting, the preferred physical local H1 alternative is

\[
\boxed{
\operatorname{Def}\ge I_1^2/184320
\quad\text{or}\quad
R_{pair}\ge I_1^2/184320.
}
\]

Together with `I_B<=sqrt(6)I_3+I_1` and the stronger H3 branch, a physical
low-strain mild-aspect full-curvature impulse has common clean cost

\[
\boxed{c_{mild,curv}^{phys}(I_B)=I_B^2/737280}
\]

outside source/dephasing, nonlinear-feedback, large-daughter and pair-rescue
branches.  Pair rescue routes to `Xi`, dominant reuse, Bellman entropy or
ancestry cycles.

These constants supersede the idealized `1/25600` and `1/102400` calculation
which treated the relative-coordinate pullback as an isometry and did not split
the daughter among three physical roles.  The pointwise Arb theorem
`Q_pol>=1/25||B_hook||^2` is unchanged.

The H1 source branch is no longer unnamed.  On the original mild strip,

\[
J_1\le2\int\|S_{src}\|+54\int\|A\|\|B\|,
\]

so dephasing reaches pressure-third, differentiated SGS, viscous-fourth, or a
conservative velocity-gradient/frame-coupling branch.  The latter is not called
a pure strain defect.

## Filtered SGS / Onsager source insertion

For the strict transporter `V=S_(N/4)u`,

\[
\operatorname{supp}\widehat R\subset B_{N/2}.
\]

Consequently

\[
\boxed{N^{-4}\|S_R\|\le(3/2000)s\|R\|_{3/2}},
\qquad
\boxed{N^{-4}\|S_\nu\|\le(\nu s/5000)\sqrt{d_V}}.
\]

A differentiated-SGS source therefore becomes the same cubic increment currency
already generated by the macroscopic `RU` boundary ledger; it must not be
charged twice.  Standard 3D LP/Bernstein gives

\[
\boxed{
Q_N\le g_1(C_{LP}C_B)^3
[(4/3)(\beta R_G)^2\mu_{low}^{max}+2d_{high}]^{3/2}.
}
\]

Thus an increment event yields a low/base critical-mass reservoir or
high-frequency normalized enstrophy.  Packetizing a winning mass reservoir with
`theta=1/4`, `alpha=1/2` gives a dominant one-quarter atom, ancestry entropy at
least `log 2`, or same-ancestry pair/cycle mass at least `1/4`.

The source-to-increment and increment-to-mass powers cancel, so SGS source
**weight** converts linearly to mass/enstrophy; no temporal-persistence
hypothesis is needed.  The viscous branch has `d_V>=b rho_nu^2`, so temporal
concentration only raises its dissipation cost by Cauchy.

## No-gap affine-aspect insertion

The transfer-facing hook bridge extends to

\[
\boxed{\operatorname{cond}L\le567/500}
\]

with

\[
Q_{pol}\ge\frac1{4000}\|B_{hook}\|^2,
\qquad
\sum_i\|F_i^{H1}\|^2\ge\frac1{16000}\|B_{hook}\|^2.
\]

The conservative physical local costs are

\[
\boxed{\operatorname{Def}\ \text{or}\ R_{pair}\ge I_1^2/28\,800\,000}
\]

and

\[
\boxed{\operatorname{Def}\ \text{or}\ R_{pair}\ge I_B^2/115\,200\,000}.
\]

Beyond the strip, aspect is not charged.  Physical covariance satisfies

\[
\frac d{dt}\log\operatorname{cond}L\le2\|\operatorname{sym}A\|_{op},
\]

with a nonpositive viscous contribution.  Hence on `int||sym A||<=1/30`,

\[
\operatorname{cond}L(t_1)>567/500
\Longrightarrow
\operatorname{cond}L(t_0)>21/20
\]

for the same ancestry.  If no predecessor exists, shell concentration plus the
affine critical-mass theorem gives the fresh-radius price

\[
\boxed{N\int_{E_2}|u|^2>(1/5)\operatorname{cond}(L)^{1/3}}.
\]

Thus aspect itself remains an exact Young symmetry: the master ledger sees
bounded-strip H1 deformation, inherited ancestry, or fresh physical radius --
never a free-standing aspect deficit.

## Reservoir synchronization and source sampling insertion

On the signed-good core, consecutive sticky lineage scales satisfy

\[
N_{q+1}/N_q>8/5.
\]

A materially reused low-frequency reservoir on the low-strain branch can grow in
frequency by less than `21/20`.  Hence its low-band increment coefficient per
unit critical mass contracts by

\[
\boxed{(21/32)^2=441/1024<1/2},
\]

and per unit physical energy by

\[
\boxed{(21/20)^3(5/8)^2=231525/512000<1/2}.
\]

Even with the full global energy cap available at each time, one fixed reservoir
has finite total future service capacity.  Material identity is tracked by the
exact affine covector `q=L^T k`; a new incompatible carrier is relinking, not
free reuse.

For `N^-1`-separated H1 source centers, standard Plancherel--Polya sampling gives
scale-invariant additive source charges:

\[
\sum_a\rho_{R,a}^{3/2}\lesssim Q_{inc},
\qquad
\sum_a\rho_{\nu,a}^2\lesssim d_V.
\]

The strict filtered pressure satisfies

\[
-\Delta P=\partial_i\partial_j(V_iV_j+R_{ij}),
\]

so its pressure-third samples route to low-pass critical mass plus the same
increment charge.  One fixed low-low pressure pair has coefficient ratio

\[
\boxed{(21/20)^4(5/8)^3=194481/655360<1/3}
\]

per generation, hence total future pair capacity below `3/2` times its initial
energy-capacity coefficient.

A maximal separated source subfamily also provides the resolvable-cluster count;
sub-bandwidth multiplicity is merged rather than charged repeatedly.

The whole-old-pool continuation sharpens this further.  If the old reservoir
family has a frame budget `sum E_a<=P E_global`, then its total low-band service
obeys

\[
\boxed{
C_{old}(q)
\le
\alpha^3N_0P E_{global}
(231525/512000)^q
<2^{-q}C_{old}(0).
}
\]

Thus a uniform service threshold forces newly relinked spectral capacity after
finitely many generations.  In the active finite quadratic graph, relinking
components satisfy

\[
\boxed{(n-1)+\beta=2m},
\]

so they are fresh-rich or cycle-rich.

For equal-covariance affine Gaussian synthesis packets, the frame budget itself
is explicit.  With

\[
\zeta=(L^{-1}X/2,L^Tk),
\qquad
|\langle g_a,g_b\rangle|=e^{-|\zeta_a-\zeta_b|^2/2},
\]

a 5-separated family has the Arb-certified Riesz bound

\[
\boxed{
47I/50\le G\le53I/50,
\qquad
\sum|c_a|^2\le(50/47)\|\sum c_ag_a\|_2^2.
}
\]

There is also an exact analysis-side formulation with no discrete frame loss.
For a normalized affine Gaussian window,

\[
\boxed{
\int|\mathcal V_Lu(X,k)|^2\frac{dXdk}{(2\pi)^3}
=\|u\|_2^2.
}
\]

Thus coherent phase-space cells have positive energies summing exactly to the
physical band energy: `P=1`.  The Riesz constant `50/47` is only needed when the
analysis cells are replaced by an actual separated synthesis family.

## Coherent nonlinear work, increment edges and stopping epochs

The transfer-to-coherent-cell step is now exact at the physical-work level.  For
a normalized coherent window,

\[
\boxed{
\mathcal W_C
=2\Re\int_C\mathcal V_gf\,\overline{\mathcal V_gF}\,d\mu,
\qquad
\sum_C\mathcal W_C=2\Re\langle f,F\rangle.
}
\]

Thus a Leray-projected nonlinear band forcing has a canonical coherent-cell work
law.  Material affine cells use `zeta=(L^-1X/2,L^Tk)` and therefore acquire no
common-affine interface forcing.  For a piecewise selected material family,

\[
\boxed{P_+\le E_{final}+P_-+R_{switch}},
\]

so at least one of terminal coherent energy, physical backflow/cancellation or
relinking symmetric-difference mass is `>=P_+/3`.

The actual SGS increment gives an even more transfer-facing positive edge law.
With

\[
Y={(Q/g_1)^{2/3}\over(C_{LP}C_B)^2},
\]

if `d_high<Y/4` then

\[
S_{low}(r)=\sum_{j\le0}M_j\|\delta_ru_j\|_2^2\ge Y/2
\]

for one physical filter displacement `r`.  Moyal decomposes this exactly into

\[
\boxed{
s_{j,C}=M_j\int_C|\mathcal V_{g_j}\delta_ru_j|^2d\mu,
}
\]

and translation covariance gives

\[
\boxed{s_{j,C}\le2M_j\{E_j(C)+E_j(C-r)\}.}
\]

Once old-pool capacity is `<=Y/8`, either

\[
\boxed{\Xi_{cell}\ge Y/8},
\]

or a new coherent cluster has critical mass `>=Y/32`, or fragmentation pays
`log 2` ancestry entropy / `1/4` same-ancestry cycle mass.  This is a genuine
old/interface/new service graph, not a relabeling of aggregate band mass.

Changing only the Gaussian covariance representative is also controlled.  If
`d_log` is the affine-invariant SPD log-covariance distance, Gaussian fidelity
and window-slot Moyal give

\[
\boxed{
\int\big||\mathcal V_{g_\Sigma}f|^2-|\mathcal V_{g_\Theta}f|^2\big|d\mu
\le {d_{log}\over\sqrt2}\|f\|_2^2.
}
\]

Hence small covariance-cell updates enter `Xi_cov`; they are not fresh energy.

The old-pool ratio

\[
r=231525/512000<1/2
\]

now becomes a stopping theorem.  If a sticky epoch has `Y>=Y_0`, it must incur a
named cost by the first `q_*` with `C_0r^{q_*}<=Y_0/8`.  Relinking can restart the
clock only by paying its own interface/fresh/cycle currency.

On the differentiated-SGS H1 branch no persistence hypothesis is needed.  The
source law `Q>=c_Q rho_R^(3/2)` followed by the coherent `2/3` collision gives
`Y>=c_Y rho_R`.  For source weight `Sigma_R`, the scale-matched alternatives
include

\[
\boxed{D_{high}\ge c_Y\Sigma_R/16},
\qquad
\boxed{\int\Xi_{cell}\ge c_Y\Sigma_R/32},
\]

or

\[
\boxed{\int\mu_{coh,new}\ge c_Y\Sigma_R/128}
\]

with the same entropy/cycle alternatives.

## Uniform physical flat gate and explicit spherical erosion

For `0<tau<=1/10`, define

\[
\boxed{\delta_\tau=\tau^2/1\,036\,800\,000.}
\]

On the signed-good, transition-aspect, low-strain physical packet branch, if
transfer/pair costs, H1/H3 source impulses and objective-strain variation all
stay below their certified `tau` thresholds, then

\[
\boxed{\sqrt{E_H^{phys}}\le\tau/3},
\qquad
\boxed{(dT)_{nonconf}\le\tau/3},
\qquad
\boxed{I_B\le\tau/3}.
\]

Thus

\[
\boxed{\mathfrak F_K\le\tau.}
\]

At `tau=1/100`, the common transfer/pair gate is
`1/10,368,000,000,000`.

Signed-good triad geometry also replaces the abstract spherical perturbation.
For physical Hodge energy `H` and parent barycenter mismatch `Delta_b`,

\[
\boxed{
\left|b_c-b_1/c_*\right|
\le2\sqrt H+H/2+(5/8)\Delta_b.
}
\]

Hence on the concentrated branch

\[
\boxed{
\zeta\le4\sqrt H+H+(5/4)\Delta_b.
}
\]

If a `1%` Kelvin-flat block has `Delta_b<=1/100`, Arb gives the uniform physical
master erosion

\[
\boxed{\kappa_0>17/100.}
\]

Equal parent marginals are not required.  For a `1%` flat block, either one
parent barycenter is at most `0.99`, giving collision entropy

\[
\boxed{H_2>1/200},
\]

or both parent marginals contain directional cores of mass at least `7/9` whose
barycenter directions differ by more than one radian and whose `3/10`-chord caps
have gap `>1/3`.  Relative to a distinguished old lineage the second core is a
trackable companion that must be fresh or reused.

## Causal completion: binary branching forces sticky reuse

Fresh companion energy by itself is not a uniform cost.  Instead use quadratic
causality.  In a synchronized layered backward ancestry with one terminal packet,
let `n_j` be the number of distinct packet ancestors and put

\[
r_j=2n_{j+1}-n_j,
\qquad
\rho_j=r_j/(2n_{j+1}).
\]

The connected 3-uniform incidence graph satisfies exactly

\[
\boxed{\beta_{cycle}=\sum_j r_j},
\]

and

\[
\boxed{
{n_0\over2^L}=\prod_j(1-\rho_j),
\qquad
\sum_j-\log(1-\rho_j)=L\log2-\log n_0.
}
\]

Arb strengthens the signed-good scale window to

\[
\boxed{3/5<N_{parent}/N_{child}<5/8.}
\]

With a coherent root budget `P E_global` and root critical mass `NE>=eta`,

\[
\boxed{
\mathcal A_{reuse}
\ge L\log(48/25)
-\log(P E_{global}N_{base}/\eta).
}
\]

For Moyal `P=1` and clean `eta=1/5`,

\[
\boxed{
L\log(36/25)>\log(5E_{global}N_{base})
\Longrightarrow
\max_j\rho_j\ge1/4.
}
\]

Thus a sufficiently deep causally complete flat ancestry must contain a
reuse-rich layer.  This is a finite synchronized theorem only: raw cycle rank is
not yet a transfer-weighted cost, and a full PDE Duhamel extraction is not yet a
synchronized layer graph.

## Transfer-weighted causal reuse

The count-level reuse theorem has a physical transfer-weighted refinement.  For
child causal-transfer law `w`, duplicate each event into the two parent-role
slots with weight `w/2` and push forward through the physical parent label.
Shannon reuse information is

\[
\boxed{
\mathcal R_j
=H(child,role\mid parent)
=H(w_{j+1})+\log2-H(w_j).
}
\]

The `log 2` role entropy is baseline causality, not a defect.  For a one-terminal
ancestry,

\[
\boxed{
\sum_j\mathcal R_j=L\log2-H(w_0).
}
\]

The root coherent-energy bound implies the same `L log(48/25)` lower growth, so
sufficient causal depth forces one layer with `R_j>log(4/3)`.

The Renyi version reaches the existing master currencies directly.  If
`Q_child=sum w^2`, the two slots have `Q_slot=Q_child/2`, and parent pushforward
satisfies

\[
\boxed{
Q_{parent}=Q_{child}/2+R_{hidden}
=\frac12Q_{child}(1+\theta).
}
\]

Across depth `L`,

\[
\boxed{
\sum_j\log(1+\theta_j)=L\log2+\log Q_0.
}
\]

The same causal-depth threshold forces `theta_j>1/3` at one layer.  If
`H2_child<log2`, then

\[
\boxed{R_{hidden}>1/12.}
\]

Otherwise the atomic-to-ancestry collision theorem gives

\[
\boxed{H_{2,ancestry}\ge\tfrac12\log2}
\]

or

\[
\boxed{R_{same\ ancestry}>1/5.}
\]

Thus the finite synchronized causal layer is already expressed in component
entropy / hidden pair-cycle currencies; raw cycle count is no longer the
weighted bottleneck.

## Adjoint Kelvin--Duhamel causal insertion

For one selected packet coefficient in the interaction picture of the same
low-frequency affine/Kelvin transport, choose the backward adjoint dual
`dot psi=-G^* psi`.  Then exactly

\[
\boxed{\frac d{dt}\langle\psi,c\rangle
=\langle\psi,F_{HH}+R_{class}\rangle,}
\qquad
\boxed{z_1=z_0+I_{HH}+I_R}.
\]

With `A=|z_1|`, one has the one-step causal gate

\[
\boxed{|z_0|\ge A/4}
\quad\text{or}\quad
\boxed{|I_R|\ge A/4}
\quad\text{or}\quad
\boxed{|I_{HH}|\ge A/2}.
\]

The first branch is material inheritance, the second already belongs to the
classified interface/source ledger, and only the third is new high--high causal
generation.  A single phase aligned with `I_HH` turns its quadratic parent-pair
Duhamel atoms into a positive generation law with total positive mass at least
`|I_HH|`; the exact two-parent role baseline remains free.  Common affine strain,
Leray pressure and bulk viscosity are therefore not charged as new generation.

Signed-good geometry gives `64/25<T_parent/T_child<25/9`.  A half-child slab
carrying at least half the positive generation mass has parent **natural**
backward windows with common overlap longer than `103/50 T_child`.  This does not
assert coherent parent persistence on that full overlap; persistence/slab
matching is precisely the remaining asynchronous PDE interface problem.

The unresolved master step is now **actual Duhamel-to-causal-layer
synchronization and single-ledger telescoping**.  The asynchronous Navier--Stokes
packet extraction must produce coherent causal layers with the same material
parent/ancestry labels used by Hodge, resistance and holonomy, while every
skipped interaction or relabeling is charged once to `Xi`.  An ancestry reaching
`t=0` must terminate against the smooth initial tail.  These events must then be
telescoped with the physical `kappa_0>0.17` flat erosion and coherent/source
stopping costs without double charging.
