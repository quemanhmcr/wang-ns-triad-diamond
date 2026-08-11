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

This root-count statement remains valid **when an absolute per-root critical-mass floor is supplied independently**.  It is not inferred from Young shape rigidity alone.  The preferred amplitude--entropy theorem below removes this extra hypothesis and retains the same linear reuse slope.

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

The first branch is coefficient inheritance.  The second and third are only
role-interface and HH **coefficient obstructions**: they locate the first interval
on which clean continuation fails, but neither `|I_R|` nor `|I_HH|` is physical
work.  The same smooth carrier must reenter its exact `Q^2` energy identity.  Only
actual positive native interface or HH work selected by that gate can enter a
causal law.  Common affine strain, Leray pressure and bulk viscosity are therefore
neither hidden nor charged as synthetic Duhamel generation.

Signed-good geometry gives `64/25<T_parent/T_child<25/9`.  The companion
asynchronous theorem upgrades this natural-window overlap to an actual causal
**time-layer geometry** without assuming packet persistence.  After the first
half-slab, `alpha_1<=25/128`; recursively

\[
\boxed{\alpha_{j+1}\le(25/64)(\alpha_j+2/5)},
\qquad \boxed{\alpha_j\le10/39\quad(j\ge1)}.
\]

Every generated layer therefore has a common reference slice with sharp margin
`67/195`.  Allowing the next generated support to begin anywhere in the previous
common interval gives the corrected progress

\[
\boxed{s_j-s_{j+1}\ge(1792/4875)T_j},
\qquad
\boxed{\Delta s_L\ge(1792/7605)T_0[(64/25)^L-1]}.
\]

Thus finite-time ancestries either stop or hit `t=0` after finite depth.  The
initial surface is an explicit boundary root, not an interior fresh-grain event.

The unresolved master step is now **canonical material-label registration and
single-ledger telescoping**.  Exact common Kelvin transport is free; coherent-cell
switches and covariance-window changes already have one-boundary Moyal/`Xi_cov`
charges.  These labels must be chosen identically for Duhamel, Shannon/Renyi,
Hodge/resistance and holonomy, and all remaining frequency/profile/covariance
representation errors must be summed once in `Xi` before telescoping flat erosion,
coherent/source stopping and causal reuse costs.


## Adjoint-causal integration regression

On integration SHA `f1f0e62`, coherent transfer cells, coherent increment service,
service-or-flat rigidity, the physical flat episode, causal binary ancestry,
transfer-weighted Shannon and Renyi reuse, the adjoint Kelvin--Duhamel gate and
the master theorem were all green.  The preferred master run is `31249675942`
with `20,000` episode traces and worst margin `0`; the integrated adjoint run
`31249674254` had `339` tests plus `50,000` causal/Duhamel checks.

## Asynchronous Duhamel synchronization regression

The certified parabolic recursion uses only the signed-good scale window and the
one-step adjoint Duhamel gate.  It does not run a frozen packet backwards.  A
common slice exists for every generated layer, with sharp normalized span `<=10/39`, margin `67/195`, and fixed point `10/39`.  Initial-boundary termination is
finite by the geometric lifetime growth.  Preferred theorem artifact:
`31250408864` on `b9c17dc` (`347` tests plus `50,000` synchronization checks).

## Canonical material labels and multi-currency physical telescope

After asynchronous synchronization, use one nested dyadic address of the
intrinsic coherent coordinate `zeta=(L^-1 X/2,L^T k)` for the Duhamel causal
pushforward, Shannon/Renyi reuse and Hodge/resistance/holonomy ancestry.  Common
Kelvin transport preserves the address exactly and nested Moyal refinement has
zero switch cost.  With geometric representative schedules,

\[
\boxed{\Xi_{sym}\le2A_3L_*B_*h_0},
\qquad
\boxed{\Xi_{cov}\le\sqrt2\,\delta_0E_{global}}.
\]

Thus representation-level label registration is summable and tunably small;
actual material cell switches remain physical relink/backflow/fresh events.
Preferred certificate: `31250736899` (`352` tests + `50,000` checks).

The master can also admit globally finite additive physical-resource resets.
If every such reset is assigned one primary class `r`, consumes at least `b_r`,
and the total resource is at most `B_r`, then

\[
N_A\le\sum_r B_r/b_r.
\]

Together with flat erosion,

\[
\boxed{
N_T\ge
\frac{\kappa_0L-P_{max}-Z}{\kappa_0+P_{max}}
-\sum_r\frac{B_r}{b_r}.
}
\]

Hence any reset class with a proved **scale-independent threshold in a globally
bounded resource** modifies only the finite offset; the asymptotic coefficient
remains `c_eff=c0 kappa0/(kappa0+Pmax)`.  Critical fresh mass `N E` and normalized
dissipation `D_V` do not qualify from the energy inequality alone, because their
physical costs decay like `1/N` along geometric scale chains.  Preferred algebra
certificate: `31250867463` (`355` tests + `50,000` synthetic episodes); the
high-strain correction is `31251119223`.

The remaining PDE audit is branch-by-branch: every service-or-flat exit must be
primary-charged exactly once as either multiplicative transfer cost or an
already-proved globally bounded additive resource, while the remaining selected
spatial/profile interfaces are combined into the one total `Xi` ledger.


## High strain and resolved objective-strain correction

The high-strain exit is a definite critical dissipation event for the strict
transporter:

\[
\boxed{
\int_0^{cN^{-2}}\|S\|_{op}dt>1/30
\Longrightarrow
N\int\|\nabla V\|_2^2dt>\frac{32\pi^2}{75c}.
}
\]

This does **not** provide a uniform global reset count: on `N_j=N_0q^j`, a
constant normalized dissipation `D_V` costs only `nu D_V/N_j` in the physical
energy inequality, and the geometric sum is finite.  The same scale-discount
obstruction holds for fresh critical mass `N_j E_j=mu`.  Preferred certificate:
`31251119223`.

The objective-strain variation in the affine packet must use the resolved field
`V`, whose exact corotational source is

\[
\mathring S_V
=-S^2-\Omega^2+[S,\Omega]
-\nabla^2P-\operatorname{sym}\nabla\nabla\!\cdot R+\nu\Delta S.
\]

The filtered support gives

\[
\boxed{
\rho_P\le\mu_V/5700+\|R\|_{3/2}/380,
}
\]

so pressure-Hessian dephasing routes to resolved coherent mass or SGS
increment/service.  Preferred certificate: `31251421182`; the earlier
`31251337503` failure was a stress-diagnostic direction error, not a theorem
counterexample.

## Global divergence-free coherent frame removes the compact-window interface

For divergence-free `u`,

\[
\langle u,\mathbb P(g_ze_a)\rangle=\langle u,g_ze_a\rangle.
\]

Consequently the Leray-projected coherent probes form the same continuous
Parseval/Moyal analysis frame on `L^2_sigma`, while

\[
\langle\mathbb P(g_ze_a),\nabla p\rangle=0
\]

exactly.  The canonical global coherent ancestry therefore needs no compact
moving spatial cutoff merely to define packet coefficients or cancel pressure.
The compact-window commutator, curvature-balanced moat and localized pressure
boundary work remain relevant to optional local/CKN diagnostics, but they are
absent from the global master `Xi`.  Leray/helical variation across a narrow
frequency cell is already contained in the summable symbol-freezing ledger.
Preferred certificate: `31251605854`.


## Physical transfer defect moat closes the canonical cross-cell `Xi`

The divergence-free coherent localization operators give an exact positive
resolution of the identity.  Hence for every continuous trilinear work form,

\[
\boxed{
T(f,g,h)=\sum_{C,D,E}T(A_Cf,A_Dg,A_Eh),
}
\]

with no continuous-to-discrete reconstruction error.  The only global coherent
interface is therefore the **actual physical transfer** carried by cell triples
omitted when the ancestry components are separated.

On the `eta_0=10^{-4}` signed-good core, the physical child-transfer law and the
capacity law differ by at most `53/50`.  Combining this with the certified
single-edge stability defect gives

\[
\boxed{
\mathbb E_{phys}\mathcal D\le\frac{106}{25}\,\epsilon.
}
\]

At recursive depth `j`, put a moat directly in the scalar transfer defect:
choose `R_j>0`, divide `[R_j/2,R_j]` into `M_j` bins, delete the least-transfer
bin, and connect all active triad vertices using edges below its lower boundary.
Every cross-component edge is then in that moat or has
`mathcal D>=R_j/2`, so

\[
\boxed{
\eta_{cross,j}
\le \frac1{M_j}+\frac{2\overline{\mathcal D}}{R_j}.
}
\]

For

\[
M_j=M_0(j+2)^2,
\qquad
R_j=R_0(j+2)^2,
\]

\[
\boxed{
\sum_j\eta_{cross,j}
\le\frac{13}{20}
\left(\frac1{M_0}+\frac{2\overline{\mathcal D}}{R_0}\right).
}
\]

On the low-cost branch this makes the logarithmic cross penalty `Xi` summable and
tunably small, with no packet-count factor and no Gaussian synthesis tail.  If a
retained low-defect component percolates instead of splitting, the exact
incidence law `(n-1)+beta=2m` routes it to fresh-rich or cycle-rich ancestry.
Preferred certificate: `31252438256` on `f7d1d37` (`362` tests, `9` optional
skips, `50,000` weighted triad graphs).


## Canonical cross-cell integration regression

On integration SHA `3d1381f`, coherent localization operators, the physical
transfer defect moat, canonical material labels, Renyi causal reuse, the physical
multi-currency telescope, smooth symbol freezing and the master theorem were all
green.  The preferred master run is `31258603563` with `20,000` episode traces
and worst margin `0`.  The coherent localization run `31258595663` had exact
trilinear reconstruction residual `5.288e-14`; the defect-moat run `31258597148`
had minimum cross-bound margin `6.917e-4` over `50,000` weighted triad graphs.

## Single-charge causal quotient and physical branch compiler

The remaining branch audit now has a theorem-level **single-charge algebra on the
positive physical transfer measure**.  The key correction is that raw theorem
predicates are not assumed mutually exclusive.  They are first quotiented by
physical causal provenance: source-caused H1/H3 dephasing is one source event;
high strain and its forced `D_V` lower bound are one critical-dissipation event;
pair rescue delegates to its reuse/cycle endpoint; a Duhamel `R_class` term
delegates to the source/interface term it contains; and the initial surface is an
absorbing boundary rather than a fresh interior grain.

The exact coherent synthesis and transfer-defect moat give first

\[
\boxed{d\mathcal T=d\mathcal T_\Xi+d\mathcal T^{ret}.}
\]

A fixed multiplicative transfer loss is an absorbing gate on
`dT^ret`.  Otherwise let `tau_r` be the first hitting time of the certified
physical causal root `r`.  The first root owns the retained charge.  If truly
independent roots tie exactly, no lexicographic theorem priority is imposed:
with positive stopping weights `a_r`, use

\[
\boxed{\lambda_r={a_r\over\sum_s a_s},\qquad
 d\mathcal T_r=\lambda_r d\mathcal T^{ret}.}
\]

Duplicate manifestations of the same root are combined before this
normalization.  Consequently

\[
\boxed{
d\mathcal T
=d\mathcal T_\Xi
+d\mathcal T_{mult}
+d\mathcal T_{flat}
+d\mathcal T_{sticky}
+d\mathcal T_{reuse}
+d\mathcal T_{side}
+d\mathcal T_{src}
+d\mathcal T_{diss}
+d\mathcal T_{reset}
+d\mathcal T_0
}
\]

with nonnegative pieces and exact total mass.  Away from exact independent tie
sets this is a disjoint stopping partition; on a tie set it is an exact
Radon--Nikodym partition rather than duplicated measure.

The master projection is deliberately asymmetric.  Transfer/Renyi/direct
sideband endpoints enter `N_T`; Kelvin-flat steps enter `N_F`; a reset enters
`N_A` only with a scale-independent threshold in the same genuinely globally
bounded resource; source, coherent new mass and `D_V` remain **recursive
scale-critical currencies**; `t=0` terminates ancestry; and the selected
cross/representation measure enters `Xi` once.  Thus the compiler cannot be used
to smuggle critical `N E` or `D_V` into the finite reset count.

The stronger proposed identity between the adjoint Duhamel probability law and
physical child-energy transfer is false and is no longer required.  In the exact
flat scalar model `G=R=0`, `c(0)=0`, `c_dot=1`, normalized `dGamma=dt`, while
normalized physical transfer is `dT=2t dt`; the half-time masses are `1/2` and
`1/4`.  This is the unavoidable linear-amplitude versus quadratic-energy
difference.

The correct bridge is the selected-child energy balance.  For

\[
\dot c=Gc+F_{HH}+R_{class},\qquad G=-S_\perp-\nu|k|^2I,
\]

let `E_i=||c(t_i)||^2`, `K=int||S_perp||`, and let `W_HH^+`, `W_R^+` be
positive physical works.  Then

\[
\boxed{E_1\le e^{2K}(E_0+W_{HH}^++W_R^+)}.
\]

On `K<=1/30`, if `E0<E1/5` and `W_R^+<E1/5`,

\[
\boxed{W_{HH}^+\ge8E_1/15}.
\]

Therefore the master-facing causal law is taken directly to be

\[
\boxed{d\mathcal T_{HH}(t,\alpha)
=2[\Re\langle c(t),F_{HH,\alpha}(t)\rangle]_+dt}.
\]

These atoms have exactly the same same-time parent-pair support exposed by the
Duhamel source.  The asynchronous half-slab/parabolic synchronization proof uses
only positivity and support, so it applies to `dT_HH` without modification.
Shannon/Renyi hence receives actual physical child-transfer weights directly.
Raw `dGamma` remains an amplitude/interaction-picture diagnostic; it is not a
second conserved measure and does not need to equal `dT`.

At the response level the exact identity is

\[
2\Re\langle c,F\rangle
=\frac{2}{\|\psi\|^2}\Re\left(\overline{\langle\psi,c\rangle}\,\langle\psi,F\rangle\right)
+2\Re\langle c_\perp,F\rangle,
\]

which makes explicit the state-dependent quadratic energy lift and the genuine
orthogonal child-profile cross work.

Preferred physical-energy bridge certificate: `31262755985` on `84c8652`
(`387` tests plus `50,000` algebra/synchronization states; projection identity
residual `1.066e-14`).

The binary parent-pair witness itself no longer requires a synthesized Gaussian
packet.  Given the exact outer selected roles and coherent partitions
`sum_C A_C=I`,

\[
\mathcal N(w_1,w_2)=\sum_{C,D}\mathcal N(A_Cw_1,A_Dw_2)
\]

and

\[
\boxed{W_{CDE}=2\Re\langle A_Ew_3,\mathcal N(A_Cw_1,A_Dw_2)\rangle}
\]

reconstruct the actual selected HH work exactly.  Positive atoms are therefore an
actual binary causal work measure; negative atoms are backscatter.  Combining
with `W_HH^+>=8E1/15`, the generated branch has at least `8E1/15` atomic positive
mass before cross-cell `Xi`, and at least `(1-rho)8E1/15` after a relative physical
moat excision `rho`.  No continuous-to-discrete reconstruction term appears.
Preferred certificate: `31263195439` on `079f2f4` (`392` tests plus `20,000`
finite POVM/bilinear-work states; work residual `3.678e-15`).

This does not promote an individual `A_Cw` to a compact Fourier packet.  The
outer Fourier/helical role still carries the scale geometry; coherent cells carry
material ancestry.

The outer moving-role commutator now has an exact affine-subtracted Egorov form:

\[
(\partial_tQ+[V\cdot\nabla,Q])f
=\int K_N(y)[V(x)-V(x-y)-Ay]\cdot\nabla f(x-y)dy,
\qquad A=\nabla V(X).
\]

Hence common affine transport is exactly removed and the remainder is controlled
by resolved velocity curvature over the coherent region.  But the point-sampled
center Hessian cannot close this estimate.  The strict-lowpass divergence-free
shear

\[
V_2=a[\sin(rx_1)-\tfrac12\sin(2rx_1)],\qquad r=N/8,
\]

has `grad V(0)=Hess V(0)=0` while `d_1^3V_2(0)=3ar^3`.  The correct replacement
observable is therefore

\[
\mathcal K_C^2
=\operatorname{Var}_\gamma[L^{-1}\nabla V(X+Lz)L],
\]

which is common-affine invariant, vanishes for affine flow, and obeys

\[
\mathcal K_C^2
\le\mathbb E_\gamma
\|L^{-1}(\nabla^2V)(X+Lz)[L,L]\|_F^2
\]

by Gaussian Poincare.  Preferred certificate: `31263818795` on `ba3d345`
(`398` tests plus `50,000` Egorov/curvature states; commutator residual
`5.848e-16`).

The whole higher-Hermite spectrum can in fact be controlled without creating new
currencies by choosing the affine gauge from the coherent eddy itself.  With
`W=L^-1V(X+Lz)`, `vbar=E W`, `Abar=E grad W`, the residual
`R=W-vbar-Abar z` has Gaussian Hermite degree at least two.  Ornstein--Uhlenbeck
and creation/annihilation estimates give

\[
E|R|^2\le\tfrac12K_{coh}^2,
\qquad E|z|^2|R|^2\le7K_{coh}^2,
\]

so the complete Gaussian-core non-affine low--high forcing obeys

\[
\|F_{nonaff}\|_2/\|\psi\|_2
\le(1+|q|/\sqrt2+\sqrt7/2)K_{coh}.
\]

Radius/aspect controls `q` on the scale-matched branch, while

\[
I_K^2\le0.275568824559\,cD_V
\]

on `cond(L)<=567/500`.  Hence full coherent deformation is perturbative when
small and critical-dissipative when large; it does not require `H4/H5/...`
master currencies.  Preferred certificate: `31264233454` on `f8164ce` (`403`
tests plus `50,000` states).

The coherent averaged affine jet now has an exact source calculus.  With
`Xdot=barV`, `Ldot=barA L` and `r=V-barV-barA(x-X)`,

\[
\dot{\bar A}
=-\bar A^2-\langle(A-\bar A)^2\rangle
-\langle\nabla^2P\rangle
-\langle\nabla\nabla\cdot R_{SGS}\rangle
+\nu\langle\Delta A\rangle
-\langle r\cdot\nabla A\rangle.
\]

Gaussian regression/IBP gives

\[
\langle r\cdot\nabla A\rangle
=\langle(z\cdot R)(A-\bar A)\rangle,
\]

and therefore the two averaging corrections are bounded by
`kappa^2 K_coh^2` and `sqrt(7) kappa K_coh^2`.  Their normalized integrated
source weight is `<=1.18115356379 D_V`; with the averaged quadratic local source
the coefficient is `2.0383176489`.  Averaged filtered pressure/SGS/viscosity
inherit the existing clean collisions because probability averaging does not
increase the global supremum.  Preferred certificate: `31264579580` on
`d4841d4` (`408` tests + `50,000` states; identity residual `9.172e-16`).  Run
`31264531046` failed only on exact floating fixture equality before theorem
stress and is recorded separately.

The coherent service-or-flat gate now removes frozen Gaussian-profile persistence
from that assembly.  Define

\[
\mathfrak F_{coh}
=\sqrt{E_H^{phys}}+(dT)_{nonconf}+C(q_{max})I_K.
\]

For the clean default radius/aspect/shell branch,
`qmax=4.71207563594`, `C(qmax)=5.65481629117`; at `tau=1/100`,
`I_K<5.89468014821e-4` makes the complete non-affine Gaussian-core connection
piece `<tau/3`, while the existing Hodge and common-strain gates supply the other
two thirds.  Larger `I_K` is critical `D_V`.  The theorem reports all simultaneous
physical root candidates and delegates the primary choice to the first-causal
compiler.  Preferred certificate: `31264981294` on `184c64e` (`414` tests plus
`50,000` branch states).

## Physical-transfer-weighted amplitude--entropy reuse

Young near-extremality is homogeneous in the parent roles, so the preferred
master does not infer an absolute `N E_root>=eta` from shape.  It also no longer
needs the stronger pointwise Duhamel-parent law.

At an event define the scale-critical analysis coefficient

\[
\boxed{\alpha=\sqrt N\,|\langle u,\phi\rangle|.}
\]

Hard event roles carry actual physical transfer.  On a retained parent-pair cell
the positive child-work density obeys the sharp-Young capacity bound

\[
r_e(t)\le C_YN\,a_c(t)a_{1,e}(t)a_{2,e}(t).
\]

Normalize the actual positive work by its total `W` and compare with normalized
physical time times the finite hard pair-cell reference.  KL positivity gives,
**under the same physical law used by Shannon/Renyi**,

\[
\boxed{
\mathbb E_{d\mathcal T}
\log(\alpha_{p_1}\alpha_{p_2})
\ge
\mathbb E_{d\mathcal T}\log\alpha_c
+\log\Lambda_j.
}
\]

No `dGamma -> dT` parent-pair identification occurs.  Duhamel remains only an
exact support/adjoint identity.

If `M_j` retained hard pair cells are used, `Lambda_j` loses only the factor
`M_j^{-1}`.  With the exact two-parent `1/2` baseline,

\[
\ell_j\ge\frac12\log\Lambda_j+rac12\ell_{j+1},
\]

so

\[
\boxed{
\ell_0
\ge
\sum_{j=0}^{L-1}2^{-(j+1)}\log\Lambda_j
+2^{-L}\log\alpha_L.
}
\]

For polynomial symbol refinement `M_j<=M_0(j+3)^p`, the weighted sum of
`log M_j` is finite.  Thus cell refinement changes only a finite offset.

The Moyal/Bargmann energy anchor remains

\[
\boxed{N_rE_r\ge\beta\alpha_r^2},
\qquad
\beta=1.43386756899\times10^{-5}.
\]

Therefore

\[
H(w_0)+2\mathbb E_{w_0}\log\alpha_r
\le\log\sum_r\alpha_r^2
\]

and the exact Shannon telescope gives

\[
\boxed{
\sum_j\mathcal R_j
\ge
L\log\frac{48}{25}
-\log\frac{P E_{global}N_{base}}{\beta}
+\sum_{j=0}^{L-1}2^{-j}\log\Lambda_j
+2^{1-L}\log\alpha_L.
}
\]

Since `H_2<=H_1`, the same lower bound controls the exact Renyi action.  The
linear reuse coefficient is still `log(48/25)`; arbitrary parent-amplitude
imbalance and polynomial pair-cell refinement only alter the finite logarithmic
offset.  The older `N E_root>=eta` root-count theorem remains a valid special
case when such a floor is independently supplied, but it is not the preferred
closure.

## Exact outer role and event registration

Let

\[
V=S_{N/4}u,
\qquad h=u-V,
\qquad
\mathcal L_Vf=\mathbb P\nabla\!\cdot(V\otimes f+f\otimes V).
\]

For a scalar moving Fourier role `w=Q(t,D)u`, direct Navier--Stokes algebra gives

\[
(\partial_t+\mathcal L_V-\nu\Delta)w
=Q\mathcal B(V,V)-Q\mathcal B(h,h)
+(\partial_tQ+[\mathcal L_V,Q])u.
\]

Use an exact hard event projector `P` for physical transfer/Young/Moyal identity
and a smooth PDE envelope `Q` with

\[
\boxed{QP=P.}
\]

The hard signed-good role begins above `3N/5`; the smooth envelope begins above
`11N/20`.  On `K<=1/30`,

\[
\frac{11}{20}e^{-1/30}N>\frac12N,
\]

while `V tensor V` is supported below `N/2`.  Hence throughout the slab

\[
\boxed{Q\mathcal B(V,V)=0}
\]

and

\[
\boxed{
(\partial_t+\mathcal L_V-\nu\Delta)w
=-Q\mathbb P\nabla\!\cdot(h\otimes h)+R_Q.
}
\]

The event/PDE coefficient registration is exact:

\[
\boxed{\langle Pu,\phi\rangle=\langle Qu,P\phi\rangle.}
\]

Hard frequency/helicity projection is a pointwise vector contraction, so it does
not worsen the `L^3` or `L^2` dual-probe constants.  Helicity is an eventwise
terminal fiber mark transported by the adjoint Kelvin equation; no persistent
helical packet is assumed.

The interface must be read in the representation in which the PDE actually
propagates it.  The envelope `Q` is smooth and generally non-idempotent, so its
carrier energy is

\[
E_Q=\|Qu\|_2^2=\langle u,Q^2u\rangle.
\]

Complete smooth roles by `sum_a A_a^2=I`, preferably with the angle pair
`Q=cos(theta)`, `R=sin(theta)`.  The native smooth-interface work is

\[
J_Q=\langle u,\partial_t(Q^2)u\rangle
-2\Re\langle Q^2u,L_Vu\rangle.
\]

It is not the commutator work alone.  The exact handoff from the outer equation
is

\[
\boxed{
J_Q=2\Re\langle Qu,R_Q\rangle
-2\Re\langle Qu,L_VQu\rangle.
}
\]

Write `L_V=K+S`, with `K^*=-K` and `S^*=S`.  The square partition alone is
not yet a physical relink law.  First identify the common skew generator `G` of
the same affine/Kelvin transport used to propagate the smooth roles and require

\[
\boxed{\dot A_a+[G,A_a]=0\quad\text{for every }a.}
\]

This is the observer-gauge quotient.  Its channel energy work cancels exactly:

\[
\langle u,\partial_t(A_a^2)u\rangle
-2\Re\langle A_a^2u,Gu\rangle=0.
\]

Only after that identity is verified do we write

\[
K=G+K_{phys}.
\]

The residual `K_phys` pair law is antisymmetric and has zero total work, so it
is conservative **physical** relink.  The `S` pair law is symmetric and
reconstructs the already existing resolved strain/deformation work.  Arbitrary
motion of an analysis partition which is not generated by the certified common
transport is rejected before Hahn splitting; conservation of total channel
energy is not enough to make observer motion physical.

At an actual event, the complete orthogonal hard-role partition has its separate
exact `K/S` work law.  The resolved donor quotient traces positive hard-role
skew gain to a simultaneous negative-net donor and removes circulation from
recursive depth.  The smooth and hard measures share operator provenance but
are not identified without an explicit physical-work pushforward.

The guard against importing hard algebra into `Q` is permanent:

\[
\mathcal I_Q(K)+\mathcal I_{I-Q}(K)
=4\Re\langle Q(I-Q)u,Ku\rangle,
\]

which equals `-1` in the encoded two-dimensional counterexample.  The quadratic
complement restores zero.  A large interface or HH coefficient impulse is only
an interval locator.  Actual `Q^2` carrier energy and gauge-quotiented native
work must reenter the physical-energy gate before relink, strain, inheritance or
HH generation is named.

The earlier draft-PR runs `31401197668` / `31401197364` certified the underlying
`Q^2` algebra before this observer-gauge correction.  They are regression
history only and are not used as certification evidence for the present
`dot A+[G,A]=0` quotient or the canonical master obstruction barrier.  New
exact-SHA certification is required before promotion.

## Generated survival on actual physical work

After one support-level cross-cell `Xi` excision, the generated physical HH law
has only three local destinations:

1. a first physical cause already named by the compiler;
2. an earlier genuine HH regeneration, which is recursion rather than a new
   currency;
3. a Young/phase-good, common-slice-registered continuation.

A bad Young/phase event cannot disappear; it carries transfer-loss provenance.
A failed common-slice mark cannot be renamed packet decoherence.  A coefficient
failure is only a typed obstruction locator; it must reenter actual physical
energy/work before source/relink or HH-generation ownership is named.  The exact
initial boundary remains absorbing.

Let the registered survivor set at depth `j` carry fraction

\[
q_j={\mathcal T_j(C_j)\over\mathcal T_j(\Omega_j)}.
\]

Repeating the same physical KL proof on the restricted law gives exactly

\[
\boxed{\Lambda_{j,\rm surv}=q_j\Lambda_j.}
\]

If `q_j>=1/2`, the layer loses at most `log 2`.  Because the Shannon/Renyi lower
uses the geometric coefficient `2^{-j}`, a fully continuing ancestry with all
`q_j>=1/2` pays at most

\[
\boxed{
\sum_{j\ge0}2^{-j}\log q_j\ge-2\log2.
}
\]

This is another finite offset, not a change of linear slope.  If some
`q_j<1/2`, a majority of that layer's actual physical HH work has already left
free continuation through a named stop or earlier regeneration.

## Exact simultaneous causes are joint stops, not fractions

The preferred master does **not** require a Radon--Nikodym split of an exact
first-time tie.  Keep the whole event mass and the complete set `J` of
simultaneous first causes.  Project the joint state by terminal semantics:

- `t=0` is absorbing;
- if any cause already certifies a fixed multiplicative transfer/reuse/sideband
  cost, the joint state is terminal transfer-cost;
- otherwise a valid scale-independent threshold in a genuinely globally bounded
  resource may terminate as the existing additive reset;
- otherwise source/SGS, critical `D_V`, material/new ancestry and earlier HH
  regeneration all remain one `recurse` state.

The last line removes the artificial tie-weight problem: all those simultaneous
causes have the same master fate.  Arbitrary positive dummy tie weights can vary
by orders of magnitude without changing the preferred projection.  The older
RN split remains an optional fine-currency subledger only when genuinely
physical common-unit densities already exist.

Critical `NE` and `D_V` are still scale-critical `O(N^{-1})` physical costs and
are **not** admitted as finite additive resets.

Preferred integrated certificate: `31288063518` on `7e8a99e` (`476` tests plus
physical energy, coherent work, service/flat, complex Young, dual/Bargmann,
common-slice, outer role, nonaffine interface, event registration, physical pair
productivity, recursive witness, joint stop, Shannon/Renyi and compiler stress;
master `20,000` traces with worst margin `0`).  Stored artifact:
`recorded-results/31288063518/`.

## Local first hit and generated-survivor relay are exact; downstream supplier reentry is now explicit

The continuum first-hit problem is now expressed directly by the smooth-SGS physics rather than an artificial vector clock.  Event/support facts appear first, then continuous/absolutely-continuous slab observables, then actual positive HH work, backward common-slice registration, and only afterwards ancestry information.  Native-unit threshold debuts are Borel; exact ties remain an unsplit finite cause set.  Material-cell Moyal content is absolutely continuous and helical phase is monitored branch-free, so optimizer chatter and angle branches create no physical stops.

On a generated no-hit event the exact adjoint identity gives `|z(s)|>=|z(t)|/4`.  The key further simplification is that this coefficient already belongs to the **same smooth material carrier** `w=Q u`:

\[
\|w(s)\|_2^2\ge \frac{|z(s)|^2}{\|\psi(s)\|_2^2}.
\]

There is therefore no reason to manufacture a fresh hard packet at the common slice.  Re-anchoring the coherent chart is composition of common affine/Kelvin transport and leaves intrinsic `zeta=(L^{-1}X/2,L^T k)` exactly invariant.  It creates zero relink mass and zero new `Xi`.

The smooth carrier also does not become a replacement causal law.  If `r` is the signed physical HH work density and `0<=q<=1` is its scalar Fourier envelope, then

\[
\boxed{
\left[\int q^2r\right]_+
\le\int q^2[r]_+
\le\int[r]_+.
}
\]

Thus any positive HH generation needed to build carrier energy is dominated by actual positive physical HH work on the same event support.  Hard Fourier/helical roles are read only **at the actual nonlinear interaction event**, where orthogonal work disintegration is exact.  Between interactions the object being transported is smooth.

If coefficient continuation fails before that next event, neither the HH nor
interface Duhamel impulse is promoted to work.  The same carrier reenters the
direct `Q^2` energy identity with its actual initial/terminal energy, strain
action and positive gauge-quotiented native interface work.  Only an interface branch selected by
that gate is split into conservative smooth relink and existing strain.

A second apparent re-entry interface also disappears.  For any resolved field `V`, with `h=u-V`,

\[
\boxed{
-L_V(Q u)+Q B(V,V)-Q B(h,h)+(L_VQ-QL_V)u=-Q B(u,u).
}
\]

Hence changing from `S_(N/4)u` to the parent-scale `S_(N_p/4)u` is an exact **repartition gauge** of the same Navier--Stokes nonlinearity.  It moves interaction between the transport, HH and Heisenberg descriptions but creates no cutoff-switch forcing or currency.  The support geometry renews with a wide margin: for `3/5<N_p/N<5/8`, two low-strain transports give

\[
\frac{|\xi|}{N_p}\ge\frac{11/20}{5/8}e^{-1/15}
=\frac{22}{25}e^{-1/15}>\frac12,
\]

so `S_(N_p/4)u` still has low-low output below the relayed carrier.  The parabolic lifetime remains `64/25<T_p/T<25/9`.

Dedicated relay run `31291697710` passed `489` tests and `50,000` carrier/affine/Hahn/event states; worst two-step intrinsic-zeta residual was `2.707e-15`.  Corrected cutoff run `31291932513` passed `492` tests and `50,000` arbitrary bilinear/cutoff/scale states; worst old/new cutoff identity residual was `1.099e-14`, worst residual against full `-Q B(u,u)` was `8.563e-15`, and the minimum renewed low-low gap was `3.242566e-01`.  Stored artifacts are `recorded-results/31291697710/` and `recorded-results/31291932513/`.  Precursor run `31291900197` was fixture-only: all `492` tests passed and the certificate stopped because Arb object equality was used for exact rational lifetime endpoints; `4fcdc8c` changed only that check to `Fraction`.

The **generated-survivor re-entry algebra and support geometry are therefore supplied**.

The critical-dissipation route is now physically seeded as well, without promoting `D_V` to a reset.  Write the resolved ball as dyadic annuli with upper radii `M_j=(N/4)2^{-j}` and `mu_j(t)=M_j||P_j u(t)||_2^2`.  Since `sum_j M_j=N/2`, shell-time atoms with `mu_j<mu_*` carry at most `c mu_*/2` of normalized resolved dissipation on a lifetime `cN^-2`.  At high strain `D_V>=D_*=32pi^2/(75c)`; choosing `mu_*=D_*/c` gives

`D_V(mu_j>=mu_*) >= D_V/2`.

Thus at least half of the **actual dissipation law** already carries a lower-frequency critical ancestor, with `M_j<=N/4` and natural lifetime at least `16` child lifetimes.  No shell argmax or packet count appears.

The same law has a spatially coherent form from the PDE's own heat geometry.  If `H_N` is the heat kernel at time `1/(2N^2)`, then on `supp Vhat subset B_(N/4)`,

`e^(-1/32)||grad V||_2^2 <= N^2 int H_N(r)||delta_rV||_2^2dr <= ||grad V||_2^2`.

Hence `S_heat=N^3 int dt int H_N||delta_rV||_2^2dr` lies between `e^(-1/32)D_V` and `D_V`, and at least `e^(-1/32)/2>0.48` of the entire heat-service law retains the same critical shell-time mark.  Moyal then disintegrates that positive law into coherent spatial edges exactly, with the existing translation covariance.

Dedicated `31292418910` passed `496` tests plus `50,000` resolved-shell laws; dedicated `31292625136` passed `501` tests plus `50,000` heat-defect/ancestor laws, with relative Moyal residual `6.278e-16` and increment-covariance residual `8.016e-15`.  Artifacts are stored under `recorded-results/31292418910/` and `recorded-results/31292625136/`.

Critical dissipation is therefore no longer an anonymous recursive scalar.  Its coherent heat-edge law also has an exact material ownership partition.  If `chi_0,chi_1` are old-pool membership of the two intrinsic endpoints, then OO, mixed-interface and NN positive measures sum pointwise to the whole edge law.  Because the partition is performed after Moyal, no old/new field cross term appears.  Both endpoint labels are exactly invariant under common affine/Kelvin transport, and each class inherits the local endpoint-energy bound `|e^(-ik.r)A_1-A_0|^2<=2(|A_0|^2+|A_1|^2)`.

Dedicated run `31293279918` on `0cd4d89` passed `506` tests and `50,000` edge states, with worst affine endpoint residual `2.681e-15`, zero orientation failures and zero affine membership failures.  The artifact is stored at `recorded-results/31293279918/`.

The apparent OO-history conflict is now resolved by first-hit causality.  At the first high-strain contact the whole prior history still satisfies `K(t)<=1/30`, so a reused material frequency grows by at most `exp(1/30)`.  For a deterministic band, the heat service on `T(N)=cN^-2` is bounded by `cM^2E/N`.  Applying this bound shellwise before material ownership, and then using OO only as a positive submeasure, avoids any old/new field decomposition.  On a supplied signed-good material epoch,

`rho_OO <= (5/8)exp(1/15)=0.668086941092 < 441/640 < 7/10`,

so

`C_OO(q)<=C_OO(0)(441/640)^q`,

`sum_q C_OO(q)<=(640/199)C_OO(0)`.

Since every first high-strain contact has `S_heat>=S_*(c)`, once `C_OO(q)<=(1-f)S_*` the exact ownership partition forces `S_ON+S_NN>=fS_*`.  This is a material-capacity theorem, not a `D_V` reset count.  Dedicated run `31298235481` on `a0cf825` passed `513` tests plus `50,000` new states; full integration `31298235456` on the same SHA passed the complete causal stack and master `20,000` traces with worst margin `0`.  The dedicated artifact is stored at `recorded-results/31298235481/`.

The heat law has one more exact structure which removes ON as a separate persistent regime.  Since `delta_r P_j=P_j delta_r`, every coherent heat edge retains the deterministic Fourier shell in which the heat law was disintegrated before Moyal.  Capacity routing therefore uses the canonical old pool with its exact orthogonal-band event mark.  If either endpoint is old, that shell belongs to the current old-frequency envelope, so pointwise `old-incident=OO+ON` and

`S_OO+S_ON <= S_old-shell <= cM_old^2PE_global/N`.

Hence all old-incident heat service contracts by the same `<441/640` factor.  Let `g=e^(-1/32)/2`.  Once old-incident capacity falls below `(g/2)S_*`, NN carries at least `(1-g/2)S_heat`, while the critical shell-time set carries at least `gS_heat`.  Inclusion--exclusion on the same positive heat measure gives

`S_(NN intersect critical) >= (g/2)S_heat = (1/4)e^(-1/32)S_heat`.

Dedicated run `31299237508` on `8a98d43` passed `521` tests plus `50,000` new states; full integration `31299237482` passed the complete stack and master `20,000` traces with worst margin `0`.  Stored artifact: `recorded-results/31299237508/`.

The remaining packet-selection concern also disappears at the eventwise entrance to renewal.  Normalize the positive `NN intersect critical` heat sublaw and push it forward by its deterministic shell-time mark.  On each selected atom, `M||P_ju||_2^2>=mu_*`.  With `A=3M/4`, the hard shell lies in `{2A/3<|xi|<=4A/3}` and

`A||P_ju||_2^2 >= 8pi^2/(25c^2)`.

Take a scalar smooth envelope `Q_A=1` on the whole hard shell with lower support `3A/5`.  The shell's own normalized state `psi=P_ju/||P_ju||` gives `<Q_Au,psi>=||P_ju||` exactly, so every heat-law atom yields a smooth whole-shell carrier seed with a critical coefficient and no coherent-cell mass floor.  Under renewed low strain, `(3/5)e^(-1/30)A>A/2`, giving strict low--low exclusion; and `M<=N/4` yields `T_A/T_N>=256/9`.

The NN mark belongs to the `P_jV` heat edge and the critical coefficient to the simultaneous `P_ju` shell.  They remain distinct exact provenance marks; no inverse low-pass estimate or full-shell NN assertion is inserted.  Dedicated run `31300227437` on `2f64de3` passed `528` tests plus `50,000` seed states; full integration `31300227438` passed the complete stack and master `20,000` traces with worst margin `0`.  Stored artifact: `recorded-results/31300227437/`.

Quantitative old-incident heat capacity and eventwise carrier-seed extraction are therefore supplied on any already-provided signed-good band-addressed material epoch.  Neither OO nor ON can be an indefinitely reusable third regime, and the NN-critical remainder already generates a positive law of lower-scale smooth Fourier carrier seeds without a packet/cell floor.

The temporal extension of each such seed is now also exact.  On the backward natural interval `[max(0,t-cA^-2),t]`, monitor renewed strain, the role-interface coefficient obstruction, the HH coefficient obstruction, and the intrinsic boundary distance of the two retained NN heat-edge endpoints in their native units.  The two coefficient monitors locate physical-energy reentry and do not supply work weights.  Until a first hit, low--low remains excluded and

`z(t)=z(s)+I_HH[s,t]+I_interface[s,t]`.

Thus a branch either reaches a named strain/interface/HH/material-boundary first stop, reaches the absorbing initial boundary, or survives the full natural interval with

`|z(s)|>=|z(t)|/4`,  hence  `A|z(s)|^2>=pi^2/(50c^2)`.

This is a critical **smooth-carrier** survivor with an NN endpoint witness, not a theorem that the whole carrier energy is NN material.  Large interface or HH coefficient impulses are not reinterpreted as work.  Dedicated run `31301650158` and full integration `31301745046` on exact SHA `a2f2a3b` both passed the `534`-test theorem suite; the dedicated `50,000`-corridor stress had worst Duhamel residual `5.184e-15` and zero monitor-order/unit failures.  The artifact is stored under `recorded-results/31301650158/`.

The apparent remaining material-attachment problem is in fact unnecessary.  A full no-hit critical annular carrier creates **its own service at the renewed scale**.  Every prefix still has `|z(s)|>|z(t)|/4`; the registered affine/Kelvin/viscous analysis dual obeys the existing scale-independent natural-window bound

`||psi(s)||<=J||psi(t)||`,  `J=exp(1/30+nu c[(3/2)exp(1/30)]^2)`,

so throughout the full slab

`A||Q_Au(s)||_2^2 >= pi^2/(50c^2J^2)`.

Because the carrier stays in the transported annulus above `(3/5)exp(-1/30)A`, its intrinsic `A`-scale heat defect is uniformly positive.  The Gaussian heat displacement is unbounded, so it is not fed directly into the reservoir theorem.  Using `||delta_r w||_2^2<=4||w||_2^2`, Arb certifies that restricting to `|r|<=3/A` still retains more than half of the annular heat lower.  Hence at every slab time some actual bounded displacement satisfies

`A||delta_r Q_Au||_2^2 >= Y_0(c,nu)>0`,  `|r|<=3/A`,

and the whole natural slab carries normalized bounded heat service at least `cY_0`.  For the default certificate `c=nu=1`, `Y_0=2.33125086914e-4`.  `Y_0` is scale-independent in `A` for fixed block parameters, not a globally uniform reset resource.

Only after this renewed positive service exists is materiality read by exact Moyal from its actual two endpoints.  Therefore the old NN witness is not promoted to whole-carrier ownership, and service is not promoted to near-extremal HH efficiency.  Dedicated run `31303081994` and full integration `31303081937` on `f0ac683` certified this annular-service re-entry; the dedicated radius-3 truncation margin was `3.781e-02` and the worst renewed OO/ON/NN residual `1.137e-13`.

This also reveals a shorter high-strain route.  The positive resolved-dissipation restriction `D_V|_G` already carries at least half of `D_V` on shell-time marks with `M||P_Mu||_2^2>=mu_*`.  Normalize those weights only as a **diagnostic dissipation sampling law**, never as causal HH probabilities, and put `A=3M/4`.  Each atom is immediately a smooth critical whole-shell carrier seed.  Before renewed service exists, inspect only the three native first stops: strain, role-interface coefficient obstruction, and HH-regeneration coefficient obstruction.  Coefficient hits locate reentry of the same carrier into actual `Q^2` energy work; they do not name an owner by amplitude.  A physical hit keeps its routed owner, `t=0` absorbs, and a full no-hit corridor feeds the annular-service theorem above.  Material OO/ON/NN is assigned only afterwards from that new service law.

Thus the shortest high-strain recursion is now

`high strain -> D_V|_G -> critical smooth carrier -> named stop / t=0 / completed full-natural corridor [own-scale service witness]`.

Child-scale heat ownership, old-incident erosion and `NN intersect critical` remain valid material-capacity refinements but are no longer prerequisites for renewal entrance.  The corrected dedicated run `31303385148` and full integration `31303385157` on exact `7ef566a` passed `547` tests and the full causal stack; the first `51840e6` run was fixture-only because a unit-rescaling stress used an absolute rather than corridor-scale timestamp tolerance.  On the corrected SHA the worst Duhamel residual was `4.974e-14`, monitor-order/unit failures were zero, physical log-product margin remained `1.665`, and master `20,000` traces had worst margin `0`.

Accordingly **critical high-strain dissipation is no longer part of the universal-renewal gap**: it already enters an existing recursive stop, the absorbing boundary, or a completed full-natural corridor carrying its same-corridor coherent-service witness without making `D_V` a reset.

The remaining master-facing continuum problem is narrower again because the carrier theorem is now **generic in the critical shell mass**, not specific to high strain.  At any actual shell-time event with

`M||P_Mu(t)||_2^2>=mu_0>0`,

set `A=3M/4`.  Exact whole-shell registration gives `A|z(t)|^2>=(3/4)mu_0`.  Before materiality is assigned there are only three native first stops: renewed strain, role-interface coefficient obstruction and HH-regeneration coefficient obstruction.  A coefficient hit only locates physical-energy reentry; the observed backward horizon remains certified, so a shorter monitor path cannot masquerade as a full-natural corridor or `t=0` root.  On a full no-hit corridor,

`A|z(s)|^2 >= 3mu_0/64`,

`A||Q_Au(s)||_2^2 >= 3mu_0/(64J^2)`,

and the existing radius-three annular heat theorem gives actual bounded own-scale service

`A||delta_r Q_Au||_2^2 >= q_b 3mu_0/(64J^2)`,  `|r|<=3/A`.

This service is a positive witness on the corridor already traversed; exact Moyal OO/ON/NN rereading adds no second recursion edge.  At the same earlier endpoint the surviving carrier support lies in `(A/2,2A)`, so the exact hard shells at `A` and `2A` satisfy `max(mu_A,mu_2A)>=(2/3)A||Q_Au||_2^2`.  With `A=3M/4` these are comparable witnesses at ratios `3/4` or `3/2`, with exact ties retained jointly; no monotone progress is asserted.

Material OO/ON/NN is assigned only afterwards from the renewed positive service law.

Any certified resolved `D_V>=D_0>0` is now a supplier: `mu_0=D_0/c` makes the low-mass dissipation at most `D_0/2`, so at least half the actual `D_V` lies on qualifying shells.  These normalized dissipation weights remain diagnostic and never replace positive HH causal probabilities.  The existing dominant fresh coherent-service branch is another supplier: `M(E_C+E_{C-r})>=theta Y/8` implies `M E_shell>=theta Y/16`, hence `Y/64` at `theta=1/4`.  This uses only whole-shell domination of each cell and creates no cell mass floor.

The specialization `mu_0=32pi^2/(75c^2)` reproduces the high-strain terminal, survivor, carrier and service constants exactly.  Dedicated run `31304867746` and full integration `31304886428` on exact `ccdf9f3` certified the generic theorem: `555` tests, `50,000` stress states, minimum `D_0` retained fraction `0.500135189`, exact zero specialization/fresh-identity margins, worst survivor Duhamel residual `4.885e-15`, and zero order/unit/horizon failures.  The earlier `ac1b162` failure was fixture-only: one scope test expected `pressure/source` instead of the certificate's `source/pressure`.

Thus high strain, any unit-matched resolved-dissipation source, and dominant fresh coherent mass all enter the same local recursion

`critical shell -> named stop / t=0 / completed full-natural corridor [own-scale coherent-service witness]`.

The theorem is shell-local and does not invent signed-good scale progress relative to the supplier block.  Pressure mass occupation and unrelated high-frequency enstrophy remain on their existing reservoir/service/entropy routes rather than being relabeled `D_V`.

The material-label/physical-role quotient proposed at this stage is now canonical in `material_label_carrier_quotient.md`: with `Q_A` and the analysis probe fixed, a pure intrinsic material-label or selected-family bookkeeping change does not create a second carrier impulse.  A genuine role/probe change remains physical interface/relink currency, ancestry charges are not erased, and material ownership is reread from subsequent actual service.

Downstream work has also supplied direct objective-source routing, pressure-pair hard-shell reentry, refinement-invariant fresh SGS scale reentry, physical high-tail regeneration ownership, Fourier UV locality, and sliding natural-window high-tail shell reentry.  The remaining master-facing task is final **continuum assembly** of these certified supplier routes without double counting and without introducing observer-dependent clocks or artificial scale-independent resets.  The master remains conditional at that programme level, and there is no global-regularity proof for 3D Navier--Stokes.
