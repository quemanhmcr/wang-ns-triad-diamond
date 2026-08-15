# Native material vorticity--heat law

Status: **draft primitive theorem layer; no global-regularity claim**.

This note deliberately discards the recent temptation to build the core proof out of
operator thresholds, convex spectral entropies, owner classes, or event clocks.  It
starts from the objects Navier--Stokes itself evolves.

Let

\[
\varpi=d u^\flat
\]

be the physical vorticity two-form, let \(\Phi_t\) be the incompressible material
flow, and set

\[
\beta=\Phi_t^*\varpi,
\qquad
g=\Phi_t^*g_0.
\]

Then, exactly,

\[
\boxed{
\partial_t\beta=\nu\Delta_g\beta,
\qquad d\beta=0,
}
\]

and

\[
\boxed{
\partial_tg=2\Phi_t^*S,
\qquad \det g=1.
}
\]

There is no explicit vortex-stretching source in the material vorticity equation.
Euler changes only the metric seen by heat.  Viscosity is the only term that changes
\(\beta\) itself.

More importantly, these are not two equations waiting for external closure.  They
already form a closed autonomous NS system.  Let

\[
L_g=d\delta_g+\delta_gd=-\Delta_g
\]

be the nonnegative Hodge Laplacian on mean-zero forms.  Since `beta` is exact, the
pulled-back velocity one-form is reconstructed intrinsically by

\[
\boxed{
\alpha=\delta_gL_g^{-1}\beta,
\qquad d\alpha=\beta,
\qquad\delta_g\alpha=0.
}
\]

Put `v=alpha^sharp_g`.  Then the material metric equation is simply

\[
\boxed{\partial_tg=\mathcal L_vg.}
\]

Thus the smallest closed primitive system is

\[
\boxed{
\begin{aligned}
\partial_t\beta&=-\nu L_g\beta,\\
\partial_tg&=\mathcal L_{(\delta_gL_g^{-1}\beta)^\sharp_g}g,\\
d\beta&=0,\qquad\det g=1,\qquad\operatorname{Riem}(g)=0.
\end{aligned}
}
\]

Flatness is not an optional regularity assumption: `g=Phi^*g0`, so every metric slice
is the Euclidean metric in material coordinates.  The nonlinearity is entirely the
self-generated motion of this flat metric; the vorticity two-form itself only sees
heat.

The physical kinetic energy is also native to the same pair:

\[
\boxed{
E=\|\alpha\|_{L^2_g}^2
=\langle\beta,L_g^{-1}\beta\rangle_g.
}
\]

Hence energy, enstrophy and material heat are not separate currencies.  They are the
`H^-1_g`, `H^0_g` and gradient-flow levels of the same closed two-form.

## 1. Vorticity amplification is exactly transverse heat determinant

Use the fixed material volume \(da\) and write a nonzero two-form as

\[
\beta=i_q da.
\]

Because \(d\beta=0\), its principal-symbol characteristic covectors satisfy

\[
\boxed{\xi\cdot q=0.}
\]

Thus the admissible oscillations of the vorticity two-form live in the transverse
plane \(q^\perp\).

The principal symbol of material Hodge heat is \(g^{-1}\).  Since \(\det g=1\),
the cofactor identity in three dimensions gives

\[
\boxed{
\det\!\left(g^{-1}\big|_{q^\perp}\right)
=
\frac{q^Tgq}{|q|^2}.
}
\]

But \(g=F^TF\), so

\[
q^Tgq=|Fq|^2=|\omega|^2
\]

for the Euler-frozen material vector.  Hence

\[
\boxed{
\det\!\left(g^{-1}\big|_{q^\perp}\right)
=
\frac{|Fq|^2}{|q|^2}.
}
\]

The geometric mean of the two transverse diffusion coefficients is therefore
**exactly the vorticity amplification factor** \(|Fq|/|q|\).

This corrects the naive statement that an anisotropic material metric can simply
create a weak diffusion direction and escape.  One transverse coefficient can be
small, but the product of the two coefficients is fixed by the amount of vortex
stretching itself.

The differential version is equally exact.  Holding the material two-form direction
fixed, as Euler does,

\[
\boxed{
\partial_t\log |Fq|^2
=
\partial_t\log\det\!\left(g^{-1}|_{q^\perp}\right).
}
\]

So stretching and transverse-heat-area growth are not delayed effects.  They are the
same deformation read in vector and two-form geometry.

There is a matching primitive magnitude equation in Eulerian coordinates.  Writing

\[
\omega=m\xi,\qquad |\xi|=1,
\]

one gets directly from the vorticity equation

\[
\boxed{
D_tm
=m\,\xi\cdot S\xi
+\nu\big(\Delta m-m|\nabla\xi|^2\big).
}
\]

Thus viscosity does not merely smooth the magnitude.  Directional disorder of the
vorticity carries the exact negative term `-nu m |grad xi|^2`.  On `R^3`, the
Biot--Savart representation gives the complementary cancellation

\[
\boxed{
\xi(x)\cdot S(x)\xi(x)
=-\frac{3}{4\pi}\,PV\!\int
\frac{m(y)}{|r|^3}
(\xi(x)\cdot\hat r)
\big[(\xi(x)\times\xi(y))\cdot\hat r\big]\,dy.
}
\]

The nonlinear source therefore vanishes when the vorticity direction is parallel,
while the same loss of directional coherence is seen directly by heat.  This is an
exact primitive cancellation, but no universal pointwise inequality between the two
terms is claimed.

## 2. Two transverse directions satisfy an exact no-free-area identity

For \(\xi,\eta\in q^\perp\), put \(A=g^{-1}\).  Restriction to the two-dimensional
characteristic plane gives

\[
\boxed{
\det
\begin{pmatrix}
\xi^TA\xi & \xi^TA\eta\\
\eta^TA\xi & \eta^TA\eta
\end{pmatrix}
=
\frac{|Fq|^2}{|q|^2}
|\xi\wedge\eta|^2.
}
\]

Consequently

\[
\boxed{
\frac{|Fq|}{|q|}|\xi\wedge\eta|
\le
|\xi|_{g^{-1}}|\eta|_{g^{-1}}.
}
\]

Any genuine two-direction transverse frequency area is therefore paid by the same
heat symbol whose determinant is created by stretching.

There is a matching physical-space null law.  If an incompressible velocity gradient
has only one derivative direction,

\[
A=a\otimes\xi,
\qquad a\cdot\xi=0,
\]

then its vorticity is proportional to \(\xi\times a\), and

\[
\boxed{S\omega=0.}
\]

A one-direction incompressible structure cannot self-stretch its own vorticity.
This is the physical rank-one null behind the collinear Fourier degeneracy.  Strong
self-stretching must use geometry beyond the only configuration that could try to
hide forever in one weak transverse diffusion direction.

The last sentence is structural, not yet a global coercive estimate: nonlocal strain
from other structures can stretch an almost one-direction vortex.  The theorem is the
exact null and transverse heat geometry, not a claim that every rank-two interaction
has already been closed.

## 3. Stretching writes a permanent heat-area memory

The pointwise determinant law is stronger when integrated in time.  On the Euler
part of the material equation, \(q\) is frozen.  Let

\[
A_\perp(t)=g^{-1}(t)|_{q^\perp}.
\]

The accumulated transverse diffusion matrix is

\[
\mathcal H_\perp(t)=\int_0^t A_\perp(s)\,ds.
\]

For positive \(2\times2\) matrices, \(\sqrt{\det}\) is concave and homogeneous of
degree one.  The Minkowski determinant inequality therefore gives

\[
\boxed{
\sqrt{\det\mathcal H_\perp(t)}
\ge
\int_0^t\sqrt{\det A_\perp(s)}\,ds
=
\int_0^t\frac{|F(s)q|}{|q|}\,ds.
}
\]

This is a genuine persistence law generated by NS itself.  The weak transverse axis
may rotate arbitrarily and the anisotropy may become arbitrarily large; neither can
reset the determinant of the **integrated** heat covariance.  Every unit of geometric
vorticity amplification is written into a nonnegative heat-area memory.

For a spatially homogeneous material metric and a polarized two-form
\(\beta_0=i_q(f_0da)\) with \(f_0\) independent of the material coordinate along
\(q\), the material equation is an exact two-dimensional Gaussian heat equation.
Its covariance is \(2\nu\mathcal H_\perp\), so

\[
\|f(t)\|_\infty
\le
\frac{\|f_0\|_1}
{4\pi\nu\sqrt{\det\mathcal H_\perp(t)}}
\le
\frac{\|f_0\|_1}
{4\pi\nu\int_0^t |F(s)q|/|q|\,ds}.
\]

This affine formula is only an interpretation of the exact memory law, not the full
inhomogeneous NS theorem.  Its significance is that time-dependent rotation of the
weak diffusion direction does **not** defeat accumulated diffusion even in the most
adversarial affine history.

## 4. Full NS has one exact memory--reset inequality

The Euler history law used a fixed material polarization.  In full Navier--Stokes the
material two-form is rewritten by heat, so its vector representative `q(t)` can move.
This does not create a new mechanism.  Anchor the history at the final polarization
`q_T=q(T)`.  Applying the fixed-plane Minkowski law to `q_T` and then only the triangle
inequality gives

\[
\boxed{
\int_0^T|F_tq_t|\,dt
\le
|q_T|
\sqrt{\det\!\left(\int_0^Tg_t^{-1}|_{q_T^\perp}dt\right)}
+\int_0^T|F_t(q_t-q_T)|\,dt.
}
\]

The first term is the transverse heat-area memory of a fixed material plane.  The
second is the **only** remainder by which the actual vorticity history can leave that
plane.  But

\[
q_t-q_T=-\int_t^T\partial_sq_s\,ds,
\]

and `partial_s q_s` is precisely the vector representative of
`partial_s beta=-nu L_g beta`.  Thus the same heat operator appears on both sides of
the attempted escape: it creates the accumulated covariance and it is the only agent
allowed to move the polarization out of that memory.

This is the first direct coupling of persistence and reset in the primitive grammar.
It is still not a global coercive estimate because the reset remainder is pointwise
and the global energy law controls its natural `H^-2_g` action rather than this
strong material norm.  The remaining gap is therefore no longer qualitative: one must
convert the heat-only reset path into enough of this remainder at critical
concentration, or prove that it cannot do so.

## 5. Euler cannot reset the material vorticity; heat is the only reset mechanism

The material equation already says this in its shortest form:

\[
\boxed{\partial_t\beta=\nu\Delta_g\beta.}
\]

There is no Euler term.  In particular, the nonlinear flow can distort \(g\), rotate
physical vorticity through \(F\), and create arbitrarily complicated Eulerian
geometry, but it cannot independently rewrite the material two-form in order to erase
its own heat history.

Any actual change of the material vorticity comes from the same heat operator that
must dissipate it.

This statement has an exact global action identity.  Put

\[
L_g=-\Delta_g\ge0
\]

on exact two-forms.  Since

\[
\beta_t=-\nu L_g\beta,
\]

the instantaneous homogeneous \(H^{-2}\) norm gives

\[
\boxed{
\|\beta_t\|_{\dot H^{-2}_g}^2
=
\nu^2\|\beta\|_{L^2_g}^2.
}
\]

But \(\|\beta\|_{L^2_g}^2=\|\omega\|_2^2\), and the velocity energy law is

\[
E'(t)=-2\nu\|\omega\|_2^2.
\]

Therefore

\[
\boxed{
\frac1\nu\int_0^T
\|\partial_t\beta\|_{\dot H^{-2}_{g_t}}^2dt
=
\frac{E(0)-E(T)}2.
}
\]

This is not an analyst-imposed reset budget.  It is exactly the NS energy identity
read in the only coordinate system where Euler has disappeared from the vorticity
equation.

The metric \(H^{-2}\) is weak and this identity by itself does not close the critical
regularity problem.  Its role is conceptual and structural: **nonlinearity can only
accumulate deformation memory; every material reset is viscous and globally billed by
physical energy loss.**

## 6. The two primitive state velocities are exactly locked

The material formulation contains two moving fields, but Navier--Stokes does not let
them move independently.  Equip the determinant-one metric manifold with its native
affine-invariant speed

\[
\|\dot g\|_{\mathcal M,g}^2
:=\int\operatorname{tr}[(g^{-1}\dot g)^2]\,da.
\]

Since `g_t=2 F^T S F`, similarity invariance of trace gives

\[
\|\dot g\|_{\mathcal M,g}^2=4\|S\|_2^2.
\]

For an incompressible field on `R^3` or the periodic box,

\[
\|\omega\|_2^2=2\|S\|_2^2=\|\beta\|_{L^2_g}^2.
\]

Therefore

\[
\boxed{
\|\dot g\|_{\mathcal M,g}^2
=2\|\beta\|_{L^2_g}^2.
}
\]

Together with the heat-only reset identity,

\[
\boxed{
\|\dot\beta\|_{\dot H^{-2}_g}^2
=\frac{\nu^2}{2}\|\dot g\|_{\mathcal M,g}^2.
}
\]

And the physical energy law becomes the single speed identity

\[
\boxed{
-\dot E
=\nu\|\dot g\|_{\mathcal M,g}^2
=\frac2\nu\|\dot\beta\|_{\dot H^{-2}_g}^2.
}
\]

Thus the two apparent ways of changing the primitive state are not independent
controls.  The same enstrophy fixes both the Euler-generated metric speed and the
viscous material-vorticity reset speed.  In particular, using heat to move the
polarization out of its accumulated transverse memory is inseparable from the same
state amplitude that is deforming the metric and writing new memory.

This speed lock is global in space; it does not yet convert the weak `H^-2` reset norm
into the pointwise reset remainder of the previous section.  That norm conversion is
one precise form of the remaining critical gap.

## 7. Physical directional mismatch has exactly two material sources

The Biot--Savart stretching kernel needs a cross product between vorticity
directions.  In material variables that cross product has an exact decomposition,
not a phenomenological coherence model.

At two labels `a,b`, write

\[
\omega_a=F_aq_a,\qquad\omega_b=F_bq_b.
\]

Then

\[
\boxed{
\omega_a\times\omega_b
=F_a^{-T}(q_a\times q_b)
+F_aq_a\times[(F_b-F_a)q_b].
}
\]

The first term is the material-two-form mismatch seen through the **same inverse
deformation that defines the local heat symbol**.  Indeed

\[
|F_a^{-T}(q_a\times q_b)|^2
=(q_a\times q_b)^Tg_a^{-1}(q_a\times q_b).
\]

The second term is precisely the failure of the deformation to be locally affine
across the two labels.

The two terms in this coordinate formula must **not** be promoted to two mechanisms.
Infinitesimally they are simply the product rule

\[
d(Fq)=(dF)q+F\,dq,
\]

namely two coordinate pieces of the **single covariant variation of the vorticity
two-form**.  Naturality gives

\[
\boxed{
\|\nabla^g\beta\|_{L^2_g}
=\|\nabla\varpi\|_{L^2}.
}
\]

If `F_a=F_b`, the entire physical mismatch is the local heat-covector image of
`q_a cross q_b`; if `q_a=q_b`, the same covariant variation is carried by deformation
non-affinity.  These are gauges of one object, not separate currencies.

This identity is especially relevant because the exact Biot--Savart stretching rate
already contains the factor `xi(x) cross xi(y)`.  The covariant variation required
for nonlinear directional stretching is therefore the **same covariant derivative
which Hodge heat squares in its Dirichlet form**.

## 8. Spatial turnover is already inside the same heat operator

There is one important loophole to a same-label memory statement: Euler could move
the region of strongest deformation to fresh material labels without rewriting the
two-form on the old label.  That is not a new physical mechanism either.

Every material metric here is flat because `g=Phi^*g0`.  Hence the Weitzenboeck
curvature term vanishes and, on the exact vorticity two-form,

\[
\boxed{
\langle\beta,L_g\beta\rangle_g
=\|\nabla^g\beta\|_2^2.
}
\]

Thus spatial non-affinity of the deformation -- the only way an amplification
hotspot can vary strongly from label to label -- is seen directly by the same Hodge
heat Dirichlet form.  In physical variables, writing `omega=m xi`, this is the exact
splitting

\[
\boxed{
|\nabla\omega|^2
=|\nabla m|^2+m^2|\nabla\xi|^2.
}
\]

So heat sees both ways of making a dangerous concentrated region nonuniform:
magnitude concentration and direction/polarization variation.  Kato's inequality

\[
|\nabla^g|\beta|_g|\le|\nabla^g\beta|
\]

shows in particular that moving a large amplification factor through material space
cannot be invisible to Hodge dissipation.

This still is not the missing global coercive estimate.  Palinstrophy is not bounded
by the velocity-energy law alone, so one cannot simply integrate this Dirichlet form
and declare victory.  The structural point is narrower and important: **same-label
persistence and cross-label turnover exhaust the apparent ways to evade heat memory,
and both are already read by the same primitive material heat equation.**

## 5. Extreme metric distortion is globally sparse in material labels

One further consequence uses only the same primitive strain.  Along a material
trajectory,

\[
\frac d{dt}\log\sigma_{\max}(F)
\le \|S\|_{op}.
\]

Cauchy in time, material volume preservation, the identity
\(\|\omega\|_2^2=2\|S\|_2^2\), and the velocity energy law give

\[
\boxed{
\int
[\log\sigma_{\max}F(a,t)]_+^2\,da
\le
\frac{t\,[E(0)-E(t)]}{4\nu}.
}
\]

Hence

\[
|\{a:\sigma_{\max}F(a,t)\ge R\}|
\le
\frac{t[E(0)-E(t)]}{4\nu(\log R)^2}.
\]

This does **not** control the supremum: a hypothetical singular deformation may still
concentrate on a vanishing set of labels.  The value of the theorem is to identify the
remaining loophole without inventing a branch.  After the transverse-memory and
heat-only-reset laws, the only possible escape is a jointly concentrated process in
which the material two-form uses viscous rewriting on ever smaller sets quickly enough
to evade a heat covariance that remembers all Euler stretching.

## 9. Critical scaling separates memory from reset counting

The exact `H^-2` reset identity must not be misread as a finite reset-count theorem.
NS scaling itself rules that out.  A critical concentration at frequency `N` may have
velocity energy of order `N^-1`; on a parabolic interval of length `N^-2`, the
corresponding contribution to

\[
\frac1\nu\int\|\beta_t\|_{\dot H^{-2}_g}^2dt
\]

is again only order `N^-1`.  Along a dyadic sequence this is summable.  Thus the
energy-billed reset action by itself cannot prohibit infinitely many ever-finer
reconfigurations.

The transverse memory has the opposite critical scaling for a **persistent material
polarization**.  If the same material vorticity element is geometrically amplified by
order `N^2` on a natural interval of length `N^-2`, then

\[
\int_{I_N}rac{|F(t)q|}{|q|}dt\gtrsim1.
\]

Hence an infinite same-lineage critical amplification sequence forces divergent
accumulated transverse heat area.  There is no summable per-octave loophole in the
Minkowski memory.

This sharpens the remaining question without introducing a branch architecture.  A
hypothetical escape must continually move the dangerous amplification through
material space or alter the material two-form quickly enough to avoid a memory which
is non-summable on a persistent critical trajectory.  The first motion is read by the
flat Hodge Dirichlet form `||nabla^g beta||^2`; the second is heat itself.

This scale observation is not a proof that the two costs cannot cooperate on a
concentrating set.  It explains why a correct final theorem, if it exists in this
grammar, has to couple **memory and covariant turnover**, not count reset events.

## 10. Falsification guard: the primitive critical vorticity norm is not monotone

It is tempting, after reaching the material two-form grammar, to hope that the
scale-critical primitive norm `int |omega|^(3/2)` is the missing scalar Lyapunov
quantity.  A direct dealiased Galerkin referee rejects that shortcut.  In a mixed
Euler/viscous ensemble of smooth divergence-free states, the instantaneous derivative
of this quantity occurred with both signs (`36` positive and `122` negative samples in
the referee run).

This is numerical falsification evidence, not a theorem about distributions of
states.  Its methodological use is clear: the remaining mechanism is genuinely
historical/covariant.  The proof core should not return to searching for a scalar
critical monotone unless it is forced directly by the material equations.

## 11. What remains

The new primitive reduction is

\[
\boxed{
\text{closed vorticity 2-form}
+\text{volume-preserving material metric}
+\text{Hodge heat}.
}
\]

All three are generated directly by Navier--Stokes.

The strongest law found here is not a pointwise `nonlinearity <= viscosity` estimate.
Such pointwise barriers were explicitly falsified in the previous operator study.  It
is a **history law**:

\[
\boxed{
\text{Euler stretching cannot reset itself; it only accumulates transverse heat area.}
}
\]

and

\[
\boxed{
\text{the only mechanism capable of rewriting that material vorticity is heat itself.}
}
\]

A full no-escape theorem would now have to prove that the finite physical energy-loss
action available to heat cannot rewrite the material two-form on a concentrating
sequence rapidly enough to evade the accumulated transverse heat determinant.  That
last implication is **not proved here**.  The important reduction is that it no longer
involves spectral owner classes, event ancestry, or an analyst-defined persistence
functional: it is a direct question about the two equations NS already evolves.



## 14. Chern--Simons, criticality and viscosity are one curl operator

The primitive velocity one-form `alpha` already carries the metric-independent Abelian
Chern--Simons functional

\[
H=\int_M \alpha\wedge d\alpha.
\]

On co-closed one-forms let

\[
C_g=*_{g}d,\qquad \Lambda_g=|C_g|,\qquad L_g=C_g^2.
\]

Then the three objects which previously looked unrelated are powers of the same NS
operator:

\[
\boxed{
H=\langle\alpha,C_g\alpha\rangle_g,
\qquad
K=\langle\alpha,\Lambda_g\alpha\rangle_g,
\qquad
\text{viscous generator}=L_g=C_g^2.
}
\]

Equivalently, for the vorticity two-form `beta=d alpha`,

\[
E=\langle\beta,L_g^{-1}\beta\rangle_g,
\qquad
K=\langle\beta,L_g^{-1/2}\beta\rangle_g,
\qquad
Z=\|\beta\|_g^2.
\]

Thus the critical half derivative is not inserted by scaling folklore: it is the
positive modulus of the same signed curl operator whose bilinear form is helicity and
whose square is heat.

There is a canonical local realization.  On the half-cylinder `X=M x (0,infinity)`
with product metric `dy^2+g`, put

\[
A(y)=e^{-y\Lambda_g}\alpha,\qquad \mathcal F=d_4A.
\]

Then `d_4 F=delta_4 F=0`: this is a source-free Abelian Maxwell field generated by the
primitive curl operator.  Mode by mode, and therefore exactly,

\[
\boxed{
\int_X|\mathcal F|^2=K,
\qquad
\int_X\mathcal F\wedge\mathcal F=\pm H.
}
\]

The orientation sign in the second formula is conventional.  The positive/negative
signed-curl stocks are precisely the self-dual and anti-self-dual Maxwell energies,

\[
\boxed{
K_+=\frac{K+H}{2},\qquad K_-=\frac{K-H}{2}.
}
\]

For the boundary energy profile

\[
Q(y)=\int_M|\mathcal F(x,y)|^2dx
=2\sum_a a^2E_a e^{-2|a|y},
\]

one has

\[
\boxed{
\int_0^\infty Q(y)dy=K,\qquad
Q(0)=2Z,\qquad
-\frac14Q'(0)=M_3.
}
\]

Moreover

\[
\int_X|\partial_y\mathcal F|^2=M_3,
\qquad
\int_X|\nabla_x\mathcal F|^2=M_3,
\]

so the full four-dimensional Dirichlet energy of the Maxwell curvature is `2 M3`.
Critical viscosity is therefore bulk Maxwell-gradient dissipation of the same canonical
extension.

## 15. Pure helicity is stressless; critical growth is cross-duality stress

For a four-dimensional two-form define the Maxwell stress

\[
T_{AB}
=\mathcal F_{AC}\mathcal F_B{}^C
-\frac14G_{AB}\mathcal F_{CD}\mathcal F^{CD}.
\]

The source-free extension gives

\[
\operatorname{div}_4T=0,\qquad \operatorname{tr}_4T=0.
\]

Writing `F=F_++F_-` for the Hodge self-dual split, the algebra is exact:

\[
\boxed{T(\mathcal F_+)=T(\mathcal F_-)=0,}
\]

and

\[
\boxed{|T(\mathcal F)|^2=4|\mathcal F_+|^2|\mathcal F_-|^2.}
\]

Thus a pure helicity sector is not merely weakly nonlinear in the critical norm: its
canonical four-dimensional stress is identically zero.  Metric deformation can change
critical Maxwell energy only through overlap of opposite dualities.  Globally,

\[
\int_X|T|
\le2\sqrt{K_+K_-}
=\sqrt{K^2-H^2}.
\]

This is the primitive BPS form of the earlier near-Beltrami impedance.

The Euler material metric velocity extends to `X` with no `y` component.  Since the
Maxwell field is the energy minimizer for its boundary gauge class, metric variation
and the divergence-free stress give

\[
\boxed{
K'_{Euler}
=-\int_XT:\dot G
=2\int_M u\cdot T_{y\cdot}(0)
=2\kappa(0).
}
\]

At the boundary

\[
T_{y\cdot}(0)=-\Lambda u\times\omega,
\]

so the last equality is simply the scalar triple product already present in NS.  The
critical source is therefore a boundary Poynting flux of a source-free Maxwell field,
not an independently introduced spectral current.

Combining the bulk gradient identity with viscosity gives the full critical law in one
four-dimensional line:

\[
\boxed{
K'
=-\int_XT:\dot G
-\nu\int_X|\nabla_4\mathcal F|^2.
}
\]

No regularity conclusion is inserted here.  Large opposite-duality Maxwell stress can
still concentrate near the boundary.

## 16. The `sech` scale filter is literal harmonic-depth overlap

Normalize one curl radius by its critical Maxwell energy:

\[
\phi_r(y)=\sqrt{2r}\,e^{-ry},\qquad \|\phi_r\|_{L^2_y}=1.
\]

Then two radii have exact overlap

\[
\boxed{
\int_0^\infty\phi_r(y)\phi_s(y)dy
=\frac{2\sqrt{rs}}{r+s}
=\operatorname{sech}\!\left(\frac12\log\frac rs\right).
}
\]

Thus the earlier Sylvester/Poisson `sech` factor has no independent ontology.  It is
simply the overlap of the canonical Maxwell fields in the scale-depth coordinate.
Comparable curl radii occupy the same harmonic depth; widely separated radii do not.
No shell threshold, packet lineage or scale case is needed to state this locality.


## 17. Full Navier--Stokes is a closed spacetime curvature law

The material Hodge equation is not the only gauge in which the primitive current is
visible.  Return to physical spacetime and put

\[
\alpha=u^\flat,\qquad \beta=d\alpha,\qquad c=\delta\beta,
\qquad B=p+\frac12|u|^2.
\]

Rotational Navier--Stokes is exactly

\[
\alpha_t+dB=-\iota_u\beta-\nu c.
\]

Consequently the real Abelian spacetime connection

\[
\boxed{\mathbb A=\alpha-B\,dt}
\]

has curvature

\[
\boxed{
\mathbb F=d_4\mathbb A
=\beta-dt\wedge e,
\qquad e:=\iota_u\beta+\nu c.
}
\]

Its Bianchi identity is not an analogy:

\[
\boxed{d_4\mathbb F=0
\quad\Longleftrightarrow\quad
\beta_t+d e=0.}
\]

Together with

\[
\boxed{c=\delta\beta,\qquad
u e_{\rm visc}=\nu c,\qquad
u e_{\rm Euler}=\iota_u\beta,\qquad
u u=(\delta L^{-1}\beta)^\sharp}
\]

(with the obvious correction that only the viscous term carries the coefficient
`nu`), this is a closed vorticity formulation of NS; the Bernoulli/pressure scalar is
the temporal gauge potential needed to reconstruct the velocity one-form equation.
No pressure source appears in the curvature equation.

The key pointwise topological null is

\[
\boxed{(\iota_u\beta)\wedge\beta=0.}
\]

In vector language this is only `(u cross omega).omega=0`, but in the spacetime
curvature it has a stronger consequence.  Chern--Simons transgression gives

\[
d_4(\mathbb A\wedge\mathbb F)=\mathbb F\wedge\mathbb F,
\]

while

\[
\boxed{
\mathbb F\wedge\mathbb F
=-2\nu\,dt\wedge c\wedge\beta.
}
\]

Thus Euler has **zero bulk topological source pointwise**.  Helicity is changed only
by the Hodge current used by viscosity.  This is the local gauge form of the usual
Euler helicity conservation and viscous helicity decay law.

A direct dealiased Fourier referee on three independent smooth states gave Bianchi
residuals `4.0e-16`, `3.0e-16`, `3.0e-16`; the Euler wedge null was at
`1e-24--1e-26`; the transgression-density residual was `1e-20` or smaller.

## 18. Enstrophy obeys an exact Poynting--Joule law

The same curvature gives a local balance without introducing a new energy currency.
In vector representatives let `e=-u cross omega+nu curl omega` and
`c=curl omega`.  From `omega_t+curl e=0`,

\[
\boxed{
\partial_t\frac{|\omega|^2}{2}
+\operatorname{div}(e\times\omega)
=-e\cdot c.
}
\]

But

\[
\boxed{
-e\cdot c
=(u\times\omega)\cdot c-\nu|c|^2.
}
\]

Therefore Euler stretching and viscous palinstrophy are not two independent source
channels.  They are the reversible and Ohmic pieces of **one electromotive--current
work law** for the same closed spacetime curvature.  After spatial integration the
flux term disappears and this is exactly the enstrophy balance.

The pointwise Cauchy alignment in this work law can be arbitrarily close to saturation
at the jet level; no local angle gap is claimed.  The missing regularity mechanism is
therefore historical/geometric, not a stronger pointwise Joule inequality.

## 19. The hidden deformation group is `SO(3,3)` on two-forms

On an oriented four-dimensional vector space the wedge product gives `Lambda^2` a
nondegenerate bilinear form of signature `(3,3)`.  If `A in sl(4)` is an infinitesimal
volume-preserving deformation, its natural action on covariant two-forms is

\[
\rho(A)F=A^TF+FA.
\]

Wedge invariance is exactly

\[
\boxed{\rho(A)^TJ+J\rho(A)=0,}
\]

so the exterior-square action lands in `so(3,3)`.  The Hodge star is the involution
whose `+1` and `-1` eigenspaces are the two three-dimensional duality sectors.

The Cartan split of the physical deformation now has no extra ontology:

* if `A` is skew, `rho(A)` is Euclidean-skew and commutes with `J`; it is compact
  rotation and cannot mix dualities;
* if `A=S` is symmetric trace-free, `rho(S)` is Euclidean-symmetric and
  anticommutes with `J`; it is the noncompact boost which **only** mixes opposite
  dualities.

Moreover

\[
\boxed{\|\rho(S)\|_{HS(\Lambda^2)}^2=2|S|^2.}
\]

Equivalently, for a trace-free metric variation `h`,

\[
\boxed{
\{\dot *,*\}=0,
\qquad
\|\dot *\|_{HS}^2=2|h|^2.
}
\]

For the material NS metric this gives

\[
\boxed{
\int_M\|\partial_t *_4\|_{HS}^2=4Z,
\qquad
-E'(t)=\frac\nu2\int_M\|\partial_t *_4\|_{HS}^2.
}
\]

So the speed of the self-dual/anti-self-dual splitting itself is exactly billed by
the physical energy law.  The earlier Krein compact/boost split is therefore a shadow
of the elementary exterior-square representation of volume-preserving deformation.

The Maxwell stress has matching rigidity.  Besides being trace-free,

\[
\boxed{
T^2=\frac{|T|^2}{4}I_4,
}
\]

so its principal values are always `(+sigma,+sigma,-sigma,-sigma)`.  The stress is not
an arbitrary symmetric tensor that an adversarial proof architecture may orient at
will; it is the rank-one cross-duality moment of the same `(3,3)` geometry.

The numerical exterior-square referee over `20,000` random trace-free generators gave
`so(3,3)` residual `2.54e-15`, compact-rotation residual `0`, boost residual
`1.54e-15`, and the exact boost-speed constant to `5.65e-16` relative error.

## 20. Guard: Maxwell depth is canonical, but naive depth transport is false

The harmonic coordinate `y` localizes `|curl|`, but it must not be promoted to a fake
cascade time.  A direct Galerkin referee tested the tempting equation

\[
Q_t+\partial_y J_{\rm Poynting}=0
\]

for the depth energy profile `Q(y)=int |F(y)|^2`.  Even after optimizing one global
normalization of the proposed Poynting flux, the local profile residual was about
`0.115`.  Thus the **Poynting-only local depth continuity law is false**: Poisson
extension does not commute with the nonlinear product.

What is exact is the integrated Maxwell stress work giving `K'_Euler` at the physical
boundary and the harmonic-overlap/sech identities.  Any future five-dimensional or
depth-local law must include the actual field-variation current; it may not be inferred
from stress alone.

## Numerical PDE referee

A smooth nonvanishing 3D vorticity state was checked directly in a `24^3`
2/3-dealiased Fourier--Galerkin solver.  Reconstructing the vorticity equation
independently gave

`omega_t + (u.grad)omega - (omega.grad)u = nu Delta omega`

with relative residual `8.36e-17`.  Thus the statement that the Euler part disappears
from the material two-form derivative is an actual PDE identity in the repository
normalization, not a coordinate slogan.  On the same state the magnitude/direction
balance had relative residual `4.06e-6`; that larger residual is consistent with the
nonlinear product/truncation error introduced by forming `|omega|` and `omega/|omega|`
on a finite Galerkin grid.


## 12. The two primitive state fields obey a local Hodge current law

The weak reset identity of Section 5 is only the lowest rung of a much more rigid
operator relation.  Put

\[
h:=\partial_t g=\mathcal L_v g,
\qquad
\beta=dv^\flat_g.
\]

These are the symmetric and antisymmetric parts of the **same** covariant velocity
gradient:

\[
\boxed{
\nabla^g v^\flat=\frac12(h+\beta).
}
\]

Flatness and incompressibility therefore give the first-order conjugacy

\[
\boxed{
\delta_g\beta=-\operatorname{div}_g h.
}
\]

Since `d beta=0`,

\[
L_g\beta=d\delta_g\beta=-d\,\operatorname{div}_g h.
\]

Consequently the full material Navier--Stokes heat law is exactly

\[
\boxed{
\partial_t\beta
=\nu d\,\operatorname{div}_g(\partial_tg).
}
\]

Define the physical material viscous current

\[
\boxed{
j:=\nu\delta_g\beta=-\nu\operatorname{div}_g g_t.
}
\]

Then the same equation becomes a literal local conservation law

\[
\boxed{
\partial_t\beta+d j=0.
}
\]

Thus heat cannot teleport material polarization from one label to another.  For every
fixed material two-surface `Sigma`, Stokes gives

\[
\boxed{
\frac d{dt}\int_\Sigma\beta
=-\int_{\partial\Sigma}j.
}
\]

This is stronger ontologically than calling `beta_t` a reset: all rewriting is an
actual local flux through material boundaries, and that flux is exactly the divergence
of the metric-deformation velocity generated by the same state.

There is also a local mixed-derivative compatibility.  In flat material coordinates,

\[
\boxed{
\nabla_k h_{ij}-\nabla_i h_{kj}=\nabla_j\beta_{ki}.
}
\]

So spatial variation of the vorticity two-form and non-affinity of metric velocity are
not merely norm-comparable; they are components of the same differentiated velocity
jet.

## 13. The speed lock holds on the entire Hodge scale

The identity `||S||_2^2=||omega||_2^2/2` is modewise, not merely an integrated
coincidence.  On a flat metric, Hodge/rough functional calculus commutes with the
strain--vorticity Riesz map.  Therefore for every real `s` in the common domain,

\[
\boxed{
\|g_t\|_{\dot H^s_g}^2
=2\|\beta\|_{\dot H^s_g}^2.
}
\]

Using `j=nu delta beta` and `beta_t=-d j`, exactness/coexactness gives the single
all-scale ladder

\[
\boxed{
\|\beta_t\|_{\dot H^{s-2}_g}^2
=
\|j\|_{\dot H^{s-1}_g}^2
=
\frac{\nu^2}{2}\|g_t\|_{\dot H^s_g}^2.
}
\]

Three rungs are especially informative:

\[
\begin{array}{c|c|c}
s&\text{material rewrite}&\text{same metric motion}\\ \hline
0&H^{-2}\ \beta_t&H^0\ g_t\\
1&H^{-1}\ \beta_t&H^1\ g_t\\
2&L^2\ \beta_t&H^2\ g_t.
\end{array}
\]

Thus the weak/strong reset gap is not freedom between unrelated mechanisms.  Each
extra derivative demanded from material rewriting is exactly an extra derivative of
the same metric velocity.  In operator language, if `R_g beta=g_t` denotes the Hodge
strain reconstruction, then on exact two-forms

\[
\mathcal R_g^*\mathcal R_g=2I,
\qquad
L_g^{\rm sym}\mathcal R_g=\mathcal R_gL_g^{(2)},
\]

but these are consequences of the local symmetric/antisymmetric velocity-gradient
identity above, not a new proof object.

## 14. The Hodge current is also the exact source of vorticity stress work

The same current appears on the nonlinear side.  The stress-energy tensor of the
vorticity two-form is

\[
T_\beta(X,Y)
=\langle\iota_X\beta,\iota_Y\beta\rangle_g
-\frac12|\beta|_g^2g(X,Y).
\]

In ordinary vector notation this is

\[
T_\beta=\frac12|\omega|^2I-\omega\otimes\omega.
\]

The closed-form Noether identity is

\[
\boxed{
\operatorname{div}_gT_\beta
=\iota_{(\delta_g\beta)^\sharp}\beta.
}
\]

Hence, with `c=delta_g beta`,

\[
\boxed{
\int\omega\cdot S\omega
=-\frac12\int T_\beta:g_t
=\int v\cdot\operatorname{div}T_\beta
=\langle c,\,v\times\omega\rangle.
}
\]

This is the same `c` for which viscosity has

\[
\beta_t=-\nu d c.
\]

So nonlinear enstrophy production and viscous material rewriting do not merely happen
to involve comparable derivatives.  **They use the same Hodge codifferential current.**
The tensor identity is the Noether statement behind this: the enstrophy functional

\[
\mathcal Z(g,\beta)=\frac12\|\beta\|_g^2
\]

is invariant under simultaneous diffeomorphism pullback.  Euler freezes `beta` in the
material gauge and spends the corresponding group derivative entirely through `g_t`;
viscosity moves `beta` down its Hodge gradient.  The two terms are coordinate faces of
the differential of one primitive functional.

For a periodic or decaying incompressible field, Betchov's null-Lagrangian identity
makes the same metric statement even shorter.  With

\[
B=g^{-1}g_t,
\qquad \operatorname{tr}B=0,
\]

`B` is similar to `2S`, and

\[
\boxed{
\int\omega\cdot S\omega
=-\frac16\int\operatorname{tr}(B^3)
=-\frac12\int\det B.
}
\]

Thus global Euler enstrophy production is the oriented cubic volume of the
trace-free logarithmic material-metric velocity.  This is an identity, not a sign law.

## 15. Cross-product skewness removes the apparent second-derivative source

The current formulation also exposes a useful null structure without creating a new
observable.  Since `c=curl omega=L u`,

\[
P_{\rm stretch}=\langle Lu,u\times\omega\rangle.
\]

Integrate by parts once:

\[
P_{\rm stretch}
=\sum_j\langle\partial_j u,\partial_j(u\times\omega)\rangle.
\]

The derivative which lands on the first velocity factor contributes

\[
\partial_j u\cdot(\partial_j u\times\omega)=0
\]

pointwise.  Therefore

\[
\boxed{
P_{\rm stretch}
=\sum_j\int
\partial_j\omega\cdot(\partial_j u\times u)\,dx.
}
\]

The nonlinearity is forced to use a first derivative of vorticity -- exactly the
Hodge derivative appearing in `||delta beta||_2^2=||nabla^g beta||_2^2`.  Equivalently,
if `K_beta a=i_a beta` is the pointwise skew contraction operator, then

\[
P_{\rm stretch}
=\frac12\langle c,[L^{-1},K_\beta]c\rangle.
\]

The commutator form says the same thing: if the Hodge inverse commuted with the local
vorticity rotation, stretching would vanish.  Expanding `[L,K_beta]` produces a
second-derivative `L beta` term and a first-derivative term; the former vanishes in the
quadratic form because `K_{L beta}` is skew.  Only covariant variation of `beta`
survives.

One can complete the enstrophy balance as

\[
\frac12Z'
=-\nu\sum_j\left\|
\partial_j\omega-\frac{\partial_j u\times u}{2\nu}
\right\|_2^2
+\frac1{4\nu}\sum_j\|\partial_j u\times u\|_2^2.
\]

This square is only a readout of the null identity.  The companion term must **not** be
promoted to an independent causal current.  The unresolved theorem remains historical:
large-data NS might still arrange the surviving first-derivative current coherently on
shrinking sets.  No uniform alignment gap or global-regularity conclusion is asserted.


## 21. Viscosity has an intrinsic vortex-line gauge: slip versus Frobenius twist

The physical spacetime curvature law admits a sharper factorization away from
`omega=0` which does not introduce an analysis packet or a selected scale.  Write

\[
\omega=m\xi,\qquad |\xi|=1,\qquad c=\operatorname{curl}\omega,
\]

and split the actual Hodge current orthogonally,

\[
c=c_\perp+c_\parallel,
\qquad
c_\parallel=(c\cdot\xi)\xi.
\]

Because contraction by the vorticity two-form maps transverse velocities onto
transverse one-forms, the whole perpendicular current can be absorbed into a
change of transport velocity.  Define

\[
\boxed{
 v_{\rm slip}=-\nu\frac{\omega\times c}{|\omega|^2},
 \qquad w=u+v_{\rm slip}.
}
\]

With the convention `i_v beta=-v cross omega`, the exact electromotive field is

\[
\boxed{
e=\iota_w\beta+\nu c_\parallel^\flat.
}
\]

Hence Bianchi/Faraday becomes

\[
\boxed{
\partial_t\beta+\mathcal L_w\beta
=-\nu\,d(c_\parallel^\flat).
}
\]

This corrects an over-strong reading of the earlier material-reset language.  The
perpendicular viscous current is not a genuine destruction of vortex-line
transport: it is exactly vortex-line **slip relative to the fluid**.  Only a
parallel one-form remains outside Lie transport.  Even that statement must not be
over-read as a pointwise reconnection theorem: an exact-gradient part is still
pure gauge, and only `d(c_parallel^flat)` actually rewrites the two-form in the
vortex-line frame.

The parallel current has a primitive differential-geometric meaning:

\[
\boxed{
 c_\parallel=m\tau\xi,
 \qquad
 \tau:=\xi\cdot\operatorname{curl}\xi,
}
\]

where

\[
\xi^\flat\wedge d\xi^\flat=\tau\,dV.
\]

Thus `tau` is exactly the Frobenius obstruction of the transverse plane field.
Moreover

\[
\boxed{
\xi\cdot\operatorname{curl}(m\tau\xi)=m\tau^2,
}
\]

so the parallel residual contributes a sign-definite magnitude sink
`-nu m tau^2` in the vortex-line gauge.  Positive amplification is therefore
kinematic deformation of the transported two-form; the genuine parallel residual
cannot directly amplify its magnitude.

The perpendicular part is equally concrete.  If
`kappa=(xi.grad)xi` is vortex-line curvature,

\[
\boxed{
 v_{\rm slip}
 =\nu\left(\kappa-\nabla_\perp\log m\right).
}
\]

Consequently

\[
\boxed{
\nu|c|^2
=\frac{m^2}{\nu}|v_{\rm slip}|^2
 +\nu m^2\tau^2,
}
\]

and the Poynting--Joule density completes the same geometry:

\[
\boxed{
-e\cdot c
=\frac{m^2}{4\nu}|u_\perp|^2
 -\frac{m^2}{\nu}
  \left|v_{\rm slip}+\frac12u_\perp\right|^2
 -\nu m^2\tau^2.
}
\]

There is a revealing limiting geometry.  If `tau=0`, Frobenius gives locally

\[
\omega=\mu\nabla\phi.
\]

Then

\[
c=\nabla\mu\times\nabla\phi,
\qquad
(\xi\cdot\nabla)\xi
=\nabla_\perp\log|\nabla\phi|,
\]

and line curvature cancels from the viscous slip:

\[
\boxed{
 v_{\rm slip}=-\nu\nabla_\perp\log|\mu|,
 \qquad
 \nu|c|^2
 =\nu|\nabla\phi|^2|\nabla_\perp\mu|^2.
}
\]

So the low-twist limit forced by efficient growth is not arbitrary 3D alignment:
it is an integrable transverse foliation on which perpendicular viscosity is
surface-gradient slip of one scalar amplitude.

Finally, the transverse determinant is a transported-two-form identity and does
not require an incompressible transport velocity.  For any orientation-preserving
`F`, with `J=det F` and `g=F^T F`,

\[
\boxed{
\det(g^{-1}|_{q^\perp})
=\frac{|Fq|^2}{J^2|q|^2}.
}
\]

Since `Fq/J` is exactly the vector representative of a transported two-form, the
same inverse-metric area law remains available in the vortex-line gauge even when
`w` is not volume preserving.  What remains open is a history theorem coupling
this vortex-line transport memory to the residual Frobenius-twist term without
replacing either by an analyst-defined event taxonomy.


## 22. Euler moves a fixed Hodge spectrum by a Lax conjugation

The material heat operator is not an arbitrary time-dependent elliptic operator.
Every metric in the material formulation is a pullback of the same Euclidean
metric.  Naturality therefore gives the operator identities

\[
\boxed{
\partial_tL_g=[\mathcal L_v,L_g],
\qquad
\partial_t\delta_g=[\mathcal L_v,\delta_g].
}
\]

Thus `L_g` is a genuine Lax conjugation of the fixed Euclidean Hodge Laplacian.
Its intrinsic spectrum and all intrinsic Poincare/Sobolev constants are unchanged
by Euler.  Fixed material coordinates can display extreme coefficient anisotropy,
but Euler does not manufacture weak heat eigenvalues; it moves the Hodge frame.

For the primitive current `c=delta_g beta`,

\[
\boxed{
\partial_t c
=[\mathcal L_v,\delta_g]\beta
-\nu L_g^{(1)}c.
}
\]

So viscosity heat-damps the current and Euler can regenerate it only through the
motion of the same Hodge frame generated by the state.  A direct flat-coordinate
variation referee showed that the local regeneration contraction can approach its
ordinary Cauchy envelope to better than `99.99%`; no uniform instantaneous gap is
claimed.  This is why the remaining theorem must be historical rather than a
pointwise source-versus-heat inequality.

## 23. Deformation-frame turnover is a zero-curvature `SL(3)` current

Let

\[
F=D_a\Phi,
\qquad
\Gamma=F^{-1}dF,
\qquad
B=F^{-1}F_t.
\]

These are not extra state variables: they are the Maurer--Cartan forms of the
actual deformation gradient.  They satisfy identically

\[
\boxed{
d\Gamma+\Gamma\wedge\Gamma=0,
}
\]

and the spacetime compatibility law

\[
\boxed{
\partial_t\Gamma=D_\Gamma B.
}
\]

Hence spatial deformation turnover and temporal deformation are not independently
orientable.  The connection is pure gauge.

Flatness and incompressibility lock its speed to the same Hodge current:

\[
\boxed{
\|D_\Gamma B\|_{L^2_g}^2
=\|\delta_g\beta\|_{L^2_g}^2,
\qquad
\|\nabla^g g_t\|_{L^2_g}^2
=2\|\delta_g\beta\|_{L^2_g}^2.
}
\]

Thus “moving the dangerous deformation to fresh labels” and “creating viscous
current” are the same derivative-order activity in two gauges.  The identity does
not make palinstrophy globally finite; it removes an apparent independent escape
mechanism.

## 24. The material metric has finite path length at almost every label

The determinant-one metrics form the canonical nonpositively curved symmetric
space

\[
\operatorname{SPD}_1(3)\simeq SL(3)/SO(3)
\]

with affine speed

\[
|g_t|_g^2=\operatorname{tr}[(g^{-1}g_t)^2].
\]

The primitive speed lock and energy identity give

\[
\int_0^T\!\int |g_t|_g^2\,da\,dt
=\frac{E(0)-E(T)}{\nu}.
\]

If

\[
\ell_T(a)=\int_0^T|g_t(a,t)|_g\,dt,
\]

then Cauchy in time yields the branch-free history bound

\[
\boxed{
\int\ell_T(a)^2\,da
\le
\frac{T[E(0)-E(T)]}{\nu}.
}
\]

The same right side controls the squared affine distance
`dist(I,g(a,T))^2`.  Consequently the material metric path has finite total
length for almost every label on every finite smooth interval and has an affine
limit there as the endpoint is approached.  A hypothetical singularity can still
concentrate on a null set of labels; this statement is not a supremum bound.

## 25. Convective self-stretching is connection, not a material force

There is an even shorter Lagrangian identity.  Since

\[
X_t=u(X,t),
\]

full Navier--Stokes gives

\[
\boxed{
X_{tt}=(-\nabla p+\nu\Delta u)(X,t).
}
\]

Differentiating with respect to the material label,

\[
\boxed{
F_{tt}=
(-\nabla^2p+\nu\nabla\Delta u)(X,t)F.
}
\]

Equivalently, in the velocity-gradient equation

\[
D_tA+A^2=-\nabla^2p+\nu\nabla\Delta u,
\]

the explicit `A^2` term cancels exactly against differentiating `F_t=AF`.
The quadratic convective self-stretch is therefore a connection coefficient of
material coordinates, not an independent second-order force.

For Euler, the material vorticity vector is `omega=Fq` with fixed `q`, and the
same cancellation gives

\[
\boxed{
D_t^2\omega=-(\nabla^2p)\omega.
}
\]

Pressure is the volume-preserving constraint reaction of the geodesic motion; it
should not be promoted to an independent source taxonomy.  Viscosity adds only
the actual gradient of its Hodge current.

The degree-one parent of the current law says the same thing in circulation form.
For the pulled-back velocity one-form,

\[
\boxed{
\partial_t\widetilde\alpha+d\pi
=-\nu\,\delta_g d\widetilde\alpha,
}
\]

so every material loop obeys

\[
\boxed{
\frac d{dt}\oint_\gamma\widetilde\alpha
=-\nu\oint_\gamma\delta_g\beta.
}
\]

Taking one exterior derivative recovers `beta_t+nu d delta_g beta=0`.  The
vorticity-flux continuity law is therefore the derivative of one primitive
Kelvin-current law, not a separate reset mechanism.


## 26. The Klein quadric is the primitive vortex-worldsheet manifold

The physical spacetime curvature

\[
\mathbb F=\beta-dt\wedge e
\]

lives in the six-dimensional space `Lambda^2 R^4`, whose wedge quadratic form has
signature `(3,3)`.  A nonzero four-dimensional two-form is decomposable/rank two
exactly when

\[
\boxed{\mathbb F\wedge\mathbb F=0.}
\]

Projectively this null cone is the Klein quadric `Gr(2,4)`.  For Euler,
`e=i_u beta`, so the algebraic null `(i_u beta) wedge beta=0` puts the physical
curvature identically on the Klein quadric.  Since `d_4 mathbb F=0`, the kernel of
this rank-two form is involutive.  It contains exactly the primitive directions

\[
\boxed{\partial_t+u,\qquad \omega.}
\]

Thus Helmholtz/Kelvin freezing is a literal integrable **vortex-worldsheet**
distribution in spacetime.  No material lineage is selected by an observer.

For viscous NS, at fixed spatial `beta` the Klein condition is simply
`e.omega=0`.  Write `c=c_perp+c_parallel`.  The ideal term and `nu c_perp` are
transverse to `omega`, so

\[
e_T=\iota_u\beta+\nu c_\perp
\]

remains on the fixed-`beta` Klein slice, whereas

\[
e_N=\nu c_\parallel
\]

is its orthogonal departure.  Consequently

\[
\boxed{
\operatorname{dist}_{\beta\ {
m fixed}}(\mathbb F,\mathrm{Klein})^2
=\nu^2|c_\parallel|^2,
}
\]

and the normal cost is exactly

\[
\boxed{
\nu|c_\parallel|^2
=\frac1\nu
\operatorname{dist}_{\beta\ {
m fixed}}(\mathbb F,\mathrm{Klein})^2.
}
\]

The nonlinear Poynting--Joule production sees only `c_perp`.  Hence all positive
enstrophy work is tangent to the ideal vortex-worldsheet slice; the fiber-normal
current is purely dissipative.  This is an orthogonal decomposition of the one
Hodge current, not a case taxonomy.

## 27. Pointwise off-Klein curvature is not the true flux-freezing obstruction

There is an important gauge correction.  Faraday only contains `d e`.  To freeze
the same spatial two-form with another local velocity `w`, it is enough to find a
scalar `psi` such that

\[
\boxed{e-\iota_w\beta=d\psi.}
\]

The contraction term annihilates `omega`, so `psi` is determined along vortex
lines by

\[
\boxed{
\omega\cdot\nabla\psi=e\cdot\omega=\nu c\cdot\omega.
}
\]

Away from a vorticity zero this is a first-order ODE along a one-dimensional
characteristic and is locally solvable.  After solving it, `e-d psi` is transverse
to `omega` and therefore equals `i_w beta` for a suitable local flux velocity.
Thus pointwise `c_parallel`, the Pfaffian of the physical spacetime curvature, or
`F wedge F` is **not** by itself a local reconnection theorem.

The genuine obstruction is global/leafwise.  On a closed vortex line `gamma`, a
single-valued `psi` requires

\[
\boxed{
\mathcal E_\gamma:=\oint_\gamma e
=\nu\oint_\gamma c^\flat=0.
}
\]

The ideal one-form contributes nothing on the vortex tangent.  A nonzero period is
therefore a gauge-invariant obstruction to a global flux velocity along that
closed leaf.  More generally the obstruction belongs to the leafwise cohomology
of the vortex-line foliation; closed periods are the cleanest directly visible
part of it.

This global obstruction is not free.  With arclength `s`,
`c.tangent=m tau`, so Cauchy gives

\[
\boxed{
|\mathcal E_\gamma|^2
\le
\nu^2L_\gamma\oint_\gamma m^2\tau^2ds
=
\nu L_\gamma D_{{\rm twist},\gamma}.
}
\]

Thus true failure of a single-valued flux-transport gauge on a closed vortex line
has a direct cost in the sign-definite twist dissipation already present in the
Poynting--Joule law.  What is *not* yet proved is a conversion of these linewise
costs into a global volume no-escape budget; closed vortex lines can be sparse,
long, or organized in concentrating tubes.


## 28. Curl itself has one primitive line-geometry polar law

The slip/twist formulas are not a special construction for vorticity.  They are
the polar geometry of the primitive `curl` operator itself.  Let any nonzero
vector field be

\[
b=m n,\qquad |n|=1.
\]

Define its line curvature, Frobenius twist and transverse amplitude gradient by

\[
\kappa_b=(n\cdot\nabla)n,
\qquad
\tau_b=n\cdot\operatorname{curl}n,
\qquad
A_b=\kappa_b-\nabla_\perp\log m.
\]

Direct differentiation gives the exact four-in-one identity

\[
\boxed{
\operatorname{curl}b
=m\tau_b n+m\,n\times A_b,
}
\]

\[
\boxed{
b\cdot\operatorname{curl}b=m^2\tau_b,
\qquad
b\times\operatorname{curl}b=-m^2A_b,
}
\]

and

\[
\boxed{
|\operatorname{curl}b|^2
=m^2(\tau_b^2+|A_b|^2).
}
\]

Thus `tau_b` and `A_b` are not two analyst-created mechanisms.  They are simply
the longitudinal and transverse components of one application of `curl` relative
to the field's own direction.  `tau_b` is the Frobenius twist of the orthogonal
plane field; `A_b` is its curvature-minus-transverse-concentration defect.

Navier--Stokes applies this same geometry twice in succession.  For `b=u`,

\[
\boxed{
u\times\omega=-|u|^2A_u,
\qquad
u\cdot\omega=|u|^2\tau_u.
}
\]

So the entire Euler/Lamb geometry is the transverse part `A_u`, while the
longitudinal part is the local helicity density.  For `b=omega`,

\[
\boxed{
 c=|\omega|\tau_\omega\xi
   +|\omega|\,\xi\times A_\omega,
\qquad
v_{\rm slip}=\nu A_\omega.
}
\]

Hence the exact nonlinear enstrophy work is only the cross-coupling of two
successive transverse curl defects:

\[
\boxed{
(u\times\omega)\cdot c
=|u|^2|\omega|\,
(\xi\times A_u)\cdot A_\omega.
}
\]

The vorticity twist `tau_omega` cannot enter this production at all, while
viscosity sees both pieces:

\[
\boxed{
|c|^2
=|\omega|^2(\tau_\omega^2+|A_\omega|^2).
}
\]

Thus the local Poynting--Joule law can be read without introducing any independent
coherence, phase, slip or topology currency: NS repeatedly applies one curl
geometry, uses the transverse part of the first application to drive the
transverse part of the second, and Hodge heat squares the complete second curl.
The hard global question remains persistence/concentration of this iterated
transverse geometry, not a missing local algebraic branch.


## 29. The whole pressure-free PDE is one Cartan--Hodge constitutive current law

The preceding material, spacetime, Klein, Maxwell and line-geometry laws can be
placed under a still smaller equation.  Return to physical space and keep only
the closed vorticity two-form `beta`.  In Coulomb/mean-zero gauge the velocity is
not an independent state variable:

\[
\boxed{
u^\flat=\delta L^{-1}\beta.}
\]

Cartan's formula and `d beta=0` give

\[
\mathcal L_u\beta=d\iota_u\beta,
\]

while on a closed two-form Hodge heat is

\[
L\beta=d\delta\beta.
\]

Therefore the entire three-dimensional Navier--Stokes vorticity equation is the
single autonomous constitutive law

\[
\boxed{
\partial_t\beta
=-d e_\beta,
\qquad
e_\beta
=(\iota_u+\nu\delta)\beta,
\qquad
u^\flat=\delta L^{-1}\beta,
\qquad
d\beta=0.
}
\]

Pressure has disappeared because exact one-forms have been quotiented before the
curvature equation is written.  No shell, packet, event clock, spectral owner or
selected persistence variable occurs.  The two physical pieces of the one
one-form current are simply

\[
\boxed{
e_E=\iota_u\beta,
\qquad
e_\nu=\nu\delta\beta.
}
\]

Both lower form degree by one; the same exterior derivative `d` returns them to
the curvature evolution.

The elementary exterior algebra already explains why Euler has two exact nulls.
Pointwise,

\[
\boxed{
\langle u^\flat,\iota_u\beta\rangle=0,
}
\]

and in three spatial dimensions

\[
\boxed{
(\iota_u\beta)\wedge\beta
=\frac12\iota_u(\beta\wedge\beta)=0.
}
\]

Thus energy conservation, helicity conservation and tangency of the ideal
spacetime curvature to the Klein quadric are not separate depletion theorems.
They are different readings of the same Cartan contraction identity.

The same total current `e=e_E+e_nu` generates the first three physical quadratic
balances.  With

\[
E=\langle\beta,L^{-1}\beta\rangle,
\qquad
H=\int\alpha\wedge\beta,
\qquad
Z=\|\beta\|_2^2,
\qquad
c=\delta\beta,
\]

one gets

\[
\boxed{
E'=-2\langle e,u^\flat\rangle=-2\nu Z,
}
\]

\[
\boxed{
H'=-2\int e\wedge\beta
=-2\nu\langle c,\omega\rangle,
}
\]

and

\[
\boxed{
Z'=-2\langle e,c\rangle.
}
\]

The ideal current is identically invisible to the first two pairings and can do
work only in the third.  The Poynting--Joule law, material heat equation, Kelvin
current, Klein worldsheet, slip/twist decomposition and canonical Maxwell
extension are therefore not competing mechanisms.  They are representations of
one closed-curvature constitutive equation.

A direct `16^3`/`20^3` dealiased Fourier referee reconstructed `u` from `beta`,
formed the one current `e`, and recovered the full vorticity PDE plus all three
balances at roundoff scale.  The remaining regularity question can now be stated
without any auxiliary ontology: can the Cartan part of this one current sustain
unbounded negative Joule work against its own Hodge-adjoint part on a finite
energy-loss history, despite Klein tangency, Hodge isospectrality and the global
self-return identities forced by `d^*=delta`?


## 30. The material metric speed obeys the cubic Cartan law of `SL(3)/SO(3)`

The material metric has a still more intrinsic global balance.  Put

\[
K:=g^{-1}g_t.
\]

It is self-adjoint with respect to `g`, trace free, and similar to `2S`.  Hence it
is precisely a tangent vector of the symmetric space

\[
SL(3)/SO(3).
\]

The primitive speed lock and Betchov/Piola identities become

\[
\boxed{
Z=\frac12\int\operatorname{tr}K^2,
\qquad
Z'_E=-\int\det K,
\qquad
\int|\nabla^gK|^2=2\int|\operatorname{curl}\omega|^2.
}
\]

Therefore full Navier--Stokes satisfies the one-line metric-speed law

\[
\boxed{
\frac12\frac d{dt}\int\operatorname{tr}K^2
+\nu\int|\nabla^gK|^2
=-\int\det K.
}
\]

The quadratic form is the canonical metric on the symmetric-space tangent, while
`det K=(1/3)tr K^3` is its unique independent cubic invariant.  Pointwise

\[
|\det K|\le\frac{|K|^3}{3\sqrt6},
\]

with equality at the Weyl-wall eigenvalue ratio `(1,1,-2)` up to scale/sign.
There is no improved pointwise constant from flatness: an affine incompressible
strain realizes every trace-free symmetric tangent at one point.  The missing
regularity mechanism must therefore be spatial/historical.

## 31. Both fundamental `sl(3)` Casimirs are globally null for a genuine gradient

Let the full deformation generator be

\[
B=F^{-1}F_t=P+Q\in\mathfrak{sl}(3),
\]

with Cartan decomposition `P^*=P`, `Q^*=-Q`; `K=2P`.  Pointwise,

\[
\boxed{
\operatorname{tr}B^2
=\operatorname{tr}P^2+\operatorname{tr}Q^2,
}
\]

and

\[
\boxed{
\operatorname{tr}B^3
=\operatorname{tr}P^3+3\operatorname{tr}(PQ^2).
}
\]

For `sl(3)` these degree-two and degree-three polynomials generate the invariant
polynomial algebra.  But `B` is similar to the actual gradient `A=grad u`, and
periodic/decaying exactness gives

\[
\boxed{
\int\operatorname{tr}B^2=0,
\qquad
\int\operatorname{tr}B^3=0.
}
\]

Thus the familiar quadratic strain--rotation equality and cubic Betchov
strain--vorticity identity are not two unrelated statistical facts.  They exhaust
the fundamental invariant-polynomial content of the local `SL(3)` generator, and
both are globally null because the generator is an exact velocity gradient.

## 32. The two Casimir nulls are coefficients of one finite degree-one chord law

The previous two identities collapse again.  At fixed time define

\[
\Psi_s(x)=x+s u(x).
\]

On the periodic box this map is homotopic to the identity.  Its signed Jacobian
therefore has degree one:

\[
\boxed{
\int\det(I+s\nabla u)dx=|\mathbb T^3|
\qquad\text{for every real }s.
}
\]

Since `tr A=0`,

\[
\det(I+sA)
=1-\frac{s^2}{2}\operatorname{tr}A^2+s^3\det A.
\]

The quadratic and cubic Betchov laws are exactly the two nonconstant
coefficients of this single finite identity.

There is a local transgression form as well:

\[
\boxed{
\det(I+sA)-1
=\nabla\cdot\left[
-\frac{s^2}{2}(u\cdot\nabla)u
+\frac{s^3}{3}(\operatorname{cof}A)^Tu
\right].
}
\]

For the actual material flow the coordinate-free finite version is obtained from

\[
Y_s(a,t)=X(a,t)+sX_t(a,t).
\]

Writing `F=D_aX`, differentiation of the determinant along the chord gives

\[
\boxed{
\det D_aY_s-\det F
=\operatorname{div}_a
\int_0^s
(\operatorname{cof}D_aY_r)^T X_t\,dr.
}
\]

Since `det F=1`, local tangent-chord volume defect is always a boundary flux; its
spatial integral is zero.  Strong local Cartan cubic production cannot be a
closed source inside the full `SL(3)` deformation geometry.

A tempting stronger collapse is false.  Direct periodic referees show that the
Leray projection of `div(cof K)-3 omega cross curl omega` is not zero; the
relative residual was about `0.65--0.68`.  The metric cubic variation and
vorticity Noether stress agree in the correct global pairings, not as the same
local force modulo pressure.

## 33. Guard: local `SL(3)` algebra alone permits exact blowup

The need for the global degree/finite-energy structure can be falsified sharply.
On `R^3`, take

\[
u(x,t)=A(t)x,
\]

with

\[
S=\operatorname{diag}(a,a,-2a),
\qquad
\omega=b e_1.
\]

The skew part can be chosen so that

\[
\Omega_t+S\Omega+\Omega S=0
\iff b'=ab.
\]

For

\[
a=(T-t)^{-1},
\qquad b=b_0(T-t)^{-1},
\]

`A_t+A^2` is symmetric and is exactly absorbed by a quadratic pressure.  Since
`Delta u=0`, this is simultaneously an Euler and Navier--Stokes solution for any
viscosity, with local gradient/vorticity blowup.  It has unbounded velocity and
infinite energy, so it is **not** a counterexample to the finite-energy regularity
problem.  It proves something methodological instead: no pointwise or purely
finite-dimensional `SL(3)`, `SO(3,3)` or curl-line algebra can be the missing
regularizer by itself.  The final theorem, if it exists in this grammar, must use
the global exactness/degree/finite-energy history that this affine field violates.


## 34. The same Cartan--Hodge current transports the whole local deformation spectrum

The finite tangent-chord law has an exact time-local form.  Let

\[
A=\nabla u,
\qquad
\chi_A(\lambda)=\det(\lambda I-A).
\]

The derivative of a determinant and the Piola identity for the gradient map
`lambda x-u(x,t)` give

\[
\boxed{
\partial_t\chi_A(\lambda)
=-\nabla\cdot\left[
\operatorname{cof}(\lambda I-A)^T u_t
\right].
}
\]

This is valid for every spectral parameter `lambda`.  Navier--Stokes itself says

\[
\boxed{u_t=-P e,
\qquad e=(\iota_u+\nu\delta)\beta,}
\]

so

\[
\boxed{
\partial_t\chi_A(\lambda)
=\nabla\cdot\left[
\operatorname{cof}(\lambda I-A)^T P e
\right].
}
\]

Thus the same one-form current which evolves vorticity curvature and performs
Poynting--Joule work also transports the **entire local characteristic polynomial
of the deformation generator**.  In trace-free three dimensions,

\[
\chi_A(\lambda)
=\lambda^3-\frac12\operatorname{tr}A^2\,\lambda-\det A,
\]

so the two coefficient conservation laws are

\[
\boxed{
\partial_t\operatorname{tr}A^2
=2\nabla\cdot(Au_t),
}
\]

and

\[
\boxed{
\partial_t\det A
=\nabla\cdot\left[
\left(A^2-\frac12\operatorname{tr}A^2 I\right)u_t
\right].
}
\]

There is no independent local eigenvalue charge left out of this grammar.  If
`tr A=tr A^2=det A=0`, Cayley--Hamilton gives `A^3=0`: the instantaneous generator
is nilpotent and has no nonzero exponential eigenvalue rate.  This is not a
regularity theorem, because time-dependent/nonnormal nilpotent generators can
still produce transient growth.  The structural point is that hyperbolic
characteristic-polynomial content cannot be minted by a bulk source; it is
redistributed by the same NS current.

## 35. Away from vorticity zeros, full viscous vorticity is locally frozen into an incompressible flux velocity

The local gauge correction can be sharpened one more step.  Seek a velocity `w`
and scalar `psi` such that

\[
e-\iota_w\beta=d\psi.
\]

Then Faraday becomes exact Lie transport,

\[
\partial_t\beta+\mathcal L_w\beta=0.
\]

At a point where `omega` is nonzero, the first scalar characteristic equation is

\[
\boxed{
\omega\cdot\nabla\psi=e\cdot\omega.
}
\]

It makes `e-d psi` transverse to the vorticity and therefore determines the
transverse part `w_perp` uniquely.  The parallel component of `w` remains free
because `i_omega beta=0`.  Write it as `lambda omega`.  Since `div omega=0`, the
condition that the flux velocity itself preserve volume reduces to a second
scalar ODE on the same vortex line:

\[
\boxed{
\omega\cdot\nabla\lambda
=-\nabla\cdot w_\perp.
}
\]

Both equations are locally solvable by characteristics wherever `omega != 0`.
Hence

\[
\boxed{
\partial_t\beta+\mathcal L_w\beta=0,
\qquad
\nabla\cdot w=0
}
\]

has a local solution around every nonzero-vorticity point.  In that local
vortex-flux gauge the full viscous two-form has an ordinary incompressible Cauchy
transport formula; viscosity changes the vortex transport velocity instead of
locally erasing the two-form.

This statement is deliberately local.  The two characteristic ODEs can fail to
produce single-valued global potentials/parallel corrections on closed or
recurrent vortex lines, and the construction degenerates at vorticity zeros.  The
first obstruction is the leafwise electromotive period already identified above;
the global incompressibility correction has its own compatibility along recurrent
leaves.  Thus the true remaining seams are global topology and concentration of
the flux velocity gradient, not an arbitrary local heat-reset degree of freedom.


## 38. Enstrophy is flux-weighted conformal vortex-line length

The universal curl-polar vector has a direct intrinsic line geometry.  For any
nonzero field `b=m n`, introduce the conformal metric

\[
\widetilde g_b=m^2g.
\]

A field line has conformal length

\[
\widetilde L_b(\gamma)=\int_\gamma m\,ds.
\]

The conformal Levi--Civita formula gives

\[
\boxed{
\widetilde\nabla_{\widetilde n}\widetilde n
=\frac1{m^2}
\left[(n\cdot\nabla)n-\nabla_\perp\log m\right]
=\frac{A_b}{m^2}.
}
\]

Equivalently, for an instantaneous normal shape variation `V_perp` with `m`
held as the ambient scalar field,

\[
\boxed{
\delta\int_\gamma m\,ds
=-\int_\gamma m A_b\cdot V_\perp\,ds.
}
\]

Thus `A_b` is not an auxiliary slip defect; it is the normal gradient of the
field-generated conformal line length.  For vorticity,

\[
v_{\rm slip}=\nu A_\omega
\]

points in the negative-gradient direction of this instantaneous line-shape
functional.

There is a local flux disintegration which makes the global PDE balance equally
transparent.  On a nonvanishing vortex flow-box let `dPhi` denote the vorticity
flux measure across a transverse section and `s` arclength along the field line.
Since `beta=i_omega dV`,

\[
\boxed{dV=\frac{1}{m}\,d\Phi\,ds.}
\]

Therefore

\[
\boxed{
Z=\int|\omega|^2dV
=\int d\Phi\int_\gamma m\,ds,
}
\]

so enstrophy is exactly **vorticity-flux-weighted conformal line length**.  The
curl Pythagoras gives

\[
\boxed{
P=\int|c|^2dV
=\int d\Phi\int_\gamma
m\left(|A_\omega|^2+\tau_\omega^2\right)ds.
}
\]

The transverse piece is the squared line-length shape gradient; the longitudinal
piece is Frobenius twist.  Finally

\[
\boxed{
\int(u\times\omega)\cdot c\,dV
=-\int d\Phi\int_\gamma m\,u_\perp\cdot A_\omega\,ds.
}
\]

Thus the entire enstrophy law is the flux-weighted geometric work identity

\[
\boxed{
\frac12Z'
=-\int d\Phi\int_\gamma m
\left[
 u_\perp\cdot A_\omega
 +\nu|A_\omega|^2
 +\nu\tau_\omega^2
\right]ds.
}
\]

The dangerous equality geometry is now very concrete:

\[
A_\omega=-\frac{u_\perp}{2\nu},
\qquad
v_{\rm slip}=-\frac12u_\perp,
\qquad
\tau_\omega=0.
\]

Euler must move the vortex-line geometry uphill against the same conformal-length
gradient which viscosity squares.

This is **not** an autonomous curve-shortening theorem.  The first-variation
formula freezes the instantaneous background magnitude while varying line shape;
full Navier--Stokes simultaneously transports the line, changes `m`, and can have
global flux-gauge obstruction.  Static scaling also preserves the relevant
dimensionless geometry.  The value of the law is to identify the exact geometric
quantity whose *persistence* is required for growth, not to claim a one-time
coercive gap.


## 39. Poynting equality is spatially frustrated by compulsory vortex twist

All static alignment gaps fail, but the vortex-line geometry contains a stronger
**differential compatibility**.  Write

\[
\omega=m\xi,
\qquad
A_\omega=(\xi\cdot\nabla)\xi-\nabla_\perp\log m,
\qquad
\tau=\xi\cdot\operatorname{curl}\xi.
\]

The curl-polar decomposition of `xi`, the identity `div curl xi=0`, and
`div(m xi)=0` give the exact cancellation

\[
\boxed{
\xi\cdot\operatorname{curl}A_\omega
=\xi\cdot\nabla\tau.
}
\]

Now use the actual defect from the Poynting--Joule square,

\[
\boxed{
R:=u_\perp+2\nu A_\omega.
}
\]

Because `curl u=omega`, its normal curl is

\[
\boxed{
\xi\cdot\operatorname{curl}R
=m-u_\parallel\tau+2\nu\,\xi\cdot\nabla\tau.
}
\]

Equivalently,

\[
\boxed{
m
=u_\parallel\tau
-2\nu\,\partial_s\tau
+\xi\cdot\operatorname{curl}R,
}
\]

where `partial_s=xi.grad` is arclength differentiation along the vortex line.
This is the spatial persistence law missing from every pointwise square.

The instantaneous defect from the positive Poynting envelope is

\[
\mathcal D_{\rm eq}
=\frac{m^2}{4\nu}|R|^2+\nu m^2\tau^2.
\]

At a single point, `R=tau=0` is compatible with nonzero `m` because derivatives
of `R` and `tau` can be nonzero there.  This is why the previous pointwise
Cauchy-gap searches correctly failed.  But if equality persisted on an open
region, then `R=tau=0` throughout that region, their spatial derivatives vanish,
and the boxed compatibility forces `m=0`.  Therefore

\[
\boxed{
\text{perfect productive Poynting geometry cannot persist on an open region of nonzero vorticity.}
}
\]

A branch-free quantitative consequence is

\[
\boxed{
 m^2
 \le
 3u_\parallel^2\tau^2
 +12\nu^2|\partial_s\tau|^2
 +3|\xi\cdot\operatorname{curl}R|^2.
}
\]

Thus trying to suppress the instantaneous defect cannot annihilate the
obstruction: it must reappear as twist, longitudinal twist variation, or spatial
curl of the equality residual.  The identity is an **escalation law**, not a
finished global estimate; the last two terms live one derivative higher and must
be coupled to the existing Hodge/current hierarchy rather than declared bounded.

An explicit nonvanishing periodic field

\[
\omega=(\sin z,\cos z,0.55\sin x)
\]

was used as an independent spectral referee.  The first identity held at
`9e-15--1e-13` relative error from resolutions `24` through `96`.  The full
residual-curl identity converged from the finite-product error `2.35e-6` at
`24^3` to `6.05e-15` at `64^3` and `1.83e-14` at `96^3`.

## 40. Twist-free leaves turn the same frustration into an exact Stokes flux law

When `tau=0`, Frobenius gives locally

\[
\omega=\mu\nabla\phi
\]

with leaves `Sigma={phi=const}` transverse to the vorticity.  The curl-line
geometry becomes

\[
A_\omega=-\nabla_\Sigma\log|\mu|.
\]

Hence the viscous part of `R` is an exact surface gradient.  For every leaf loop
`gamma=partial D`, Stokes gives

\[
\boxed{
\oint_\gamma R\cdot dl
=\oint_\gamma u\cdot dl
=\int_D\omega\cdot n\,dA
=\Phi_D.
}
\]

Therefore a twist-free leaf carrying nonzero vorticity flux cannot support
`R=0` on a neighborhood.  Quantitatively,

\[
\boxed{
\int_\gamma|R|^2ds
\ge\frac{\Phi_D^2}{L_\gamma}.
}
\]

This is the integrable-plane counterpart of the differential twist-frustration
identity.  If the vortex plane field is nonintegrable, `tau` itself is paid by
the sign-definite twist term.  If it becomes integrable to remove that cost, the
residual is constrained by circulation/flux on the resulting leaves.  These are
two coordinate readings of one continuous geometry; no event or owner split is
introduced.


## 41. The Poynting equality defect is one sourced covariant Gauss field

The twist/residual frustration law itself is the polar decomposition of an even
smaller vector identity.  Define the exact completed-square residual

\[
\boxed{
G:=u\times\omega-2\nu c,
\qquad c=\operatorname{curl}\omega.
}
\]

The enstrophy square becomes

\[
\boxed{
\frac12Z'
=\frac1{4\nu}
\left(\|u\times\omega\|_2^2-\|G\|_2^2\right).
}
\]

Now use only the primitive self-return identities

\[
\nabla\cdot(u\times\omega)=|\omega|^2-u\cdot c,
\qquad
\nabla\cdot c=0,
\qquad
u\cdot(u\times\omega)=0.
\]

They collapse to the local Gauss law

\[
\boxed{
|\omega|^2
=\nabla\cdot G-\frac{u}{2\nu}\cdot G.
}
\]

Thus the very field whose `L2` size measures failure of maximal enstrophy
productivity is forced by the **positive enstrophy density itself**.  It is not a
free alignment variable.

The vortex-line quantities of the previous sections are merely its orthogonal
polar components.  With

\[
R=u_\perp+2\nu A_\omega,
\]

one has

\[
\boxed{
G=-m\,\xi\times R-2\nu m\tau\xi,
}
\]

and therefore

\[
\boxed{
|G|^2=m^2|R|^2+4\nu^2m^2\tau^2.
}
\]

So `R` and `tau` should not be promoted to two mechanisms.  They are the
transverse and longitudinal readings of one sourced equality-residual field.
The normal-curl twist-frustration identity is the corresponding polar differential
reading of this Gauss law.

## 42. The Gauss law carries its own canonical Schrödinger coercivity

Let

\[
D_uG:=\nabla\cdot G-\frac{u}{2\nu}\cdot G.
\]

The exact source equation is `D_u G=|omega|^2`.  Its adjoint on scalars is

\[
D_u^*\phi=-\nabla\phi-\frac{u}{2\nu}\phi.
\]

For every smooth scalar test,

\[
\boxed{
\int|\omega|^2\phi
=-\int G\cdot
\left(\nabla\phi+\frac{u}{2\nu}\phi\right).
}
\]

Incompressibility removes the cross term exactly:

\[
\boxed{
\left\|\nabla\phi+\frac{u}{2\nu}\phi\right\|_2^2
=\|\nabla\phi\|_2^2
+\frac1{4\nu^2}\|u\phi\|_2^2.
}
\]

Hence

\[
\boxed{
\|G\|_2^2
\ge
\sup_{\phi\ne0}
\frac{\left(\int|\omega|^2\phi\right)^2}
{\|\nabla\phi\|_2^2+rac1{4\nu^2}\|u\phi\|_2^2}.
}
\]

The Riesz operator of this denominator is not selected by an analyst; it is the
normal operator of the primitive Gauss law:

\[
\boxed{
D_uD_u^*
=-\Delta+\frac{|u|^2}{4\nu^2}.
}
\]

When this operator is invertible in the relevant mean-zero/energy space, the
right side is

\[
\left\langle
|\omega|^2,
\left(-\Delta+\frac{|u|^2}{4\nu^2}\right)^{-1}
|\omega|^2
\right\rangle.
\]

The constant test `phi=1` recovers

\[
\boxed{
\|G\|_2^2\ge\frac{4\nu^2Z^2}{E},
}
\]

so the earlier invariant-normal tax is only the zero-mode shadow of this full
covariant-divergence coercivity.

A dealiased Fourier referee verified the local Gauss equation to
`2.3e-15--4.4e-15`.  On grids `12,16,20`, direct conjugate-gradient solution of
the Schrödinger problem gave lower-bound/actual-residual fractions about
`0.884,0.813,0.856`; the constant-test bound accounted for roughly
`98.7--99.6%` of that Schrödinger lower in those generic states.  These numbers
are referee diagnostics, not universal constants.  In particular no strict gap
between the actual residual and the minimum covariant-divergence field is
claimed.

The structural gain is different: **near Poynting equality is a concentration
problem for a sourced Gauss field**, not an unconstrained phase-alignment problem.
To make `G` small while `|omega|^2` is large, NS must create spatial variation on
the exact drift-divergence scale dictated by `u/(2nu)`.  This is precisely the
kind of spatial persistence cost invisible to every pointwise alignment test.


## 43. The Gauss residual is one universal identity of the primitive curl operator

The preceding residual equation is not a special accident of the NS velocity.
For every smooth divergence-free vector field `v`, define

\[
\boxed{
\mathfrak G_\nu[v]
:=v\times Cv-2\nu C^2v.
}
\]

The elementary Green/self-return identity for `curl` gives

\[
\boxed{
\left(\nabla\cdot-\frac{v}{2\nu}\cdot\right)
\mathfrak G_\nu[v]
=|Cv|^2.
}
\]

Indeed `div C^2v=0`,
`div(v cross Cv)=|Cv|^2-v.C^2v`, and
`v.G_nu[v]=-2nu v.C^2v`.  Nothing else enters.  The Poynting equality residual is
simply the member `G=G_nu[u]` selected by the actual NS state.

For this member, the original rotational equation

\[
u_t=P(u\times Cu)-\nu C^2u
\]

can be reflected exactly as

\[
\boxed{
u_t=\nu C^2u+P G,}
\]

while simultaneously

\[
\boxed{
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)G=|Cu|^2.
}
\]

The first formula should not be misread as a second physical dynamics or a new
anti-diffusive mechanism.  It is the exact reflection of the same equation about
the viscous current.  Its value is structural: the field which must oppose the
formal backward-heat direction is not freely orientable or arbitrarily small;
it is constrained at that same instant by a positive-source first-order Gauss
law.

There is a corresponding historical consequence of the local vortex-flux gauge.
Where the two characteristic gauge equations produce a smooth divergence-free
`w`,

\[
\partial_t\beta+\mathcal L_w\beta=0,
\qquad \det D\Phi_w=1.
\]

Thus the *full viscous* vorticity two-form has exactly the same transported-two-form
cofactor law and Minkowski transverse-memory inequality along `w` trajectories as
the Euler-frozen material form.  Local viscosity supplies no independent reset of
that memory; the remaining failures are globalization of the flux gauge, zeros of
vorticity, or concentration of the vortex-flux deformation itself.
## 44. The actual NS current is a curved graded differential

The universal curl Gauss identity of Section 43 is the vector shadow of a smaller
operator law which is already present in the Cartan--Hodge formulation of Section
29.  Write

\[
\alpha=u^\flat,
\qquad
\beta=d\alpha,
\qquad
d\beta=0,
\]

and for a real parameter `theta` define the degree-lowering current operator

\[
\boxed{
q_\theta:=\nu\delta+\theta\,\iota_u,
\qquad
q_\theta^*=\nu d+\theta\,\alpha\wedge .
}
\]

Nothing has been added to the PDE.  The two summands are exactly the viscous
codifferential and Euler contraction already present in the physical
electromotive current.  The only algebra needed is the exterior Leibniz rule

\[
\boxed{
\{d,\alpha\wedge\}=d\alpha\wedge=\beta\wedge .
}
\]

Taking adjoints gives the primitive self-return law

\[
\boxed{
\{\delta,\iota_u\}=(\beta\wedge)^*.
}
\]

Because `delta^2=i_u^2=0`, the whole one-parameter current family closes by
curvature:

\[
\boxed{
q_\theta^2
=\theta\nu(\beta\wedge)^*,
}
\]

and more generally

\[
\boxed{
\{q_a,q_b\}
=\nu(a+b)(\beta\wedge)^*.
}
\]

Thus vorticity is literally the failure of the NS current differential to be
nilpotent.  The full physical current is the member `theta=1`,

\[
\boxed{
e=q_1\beta=(\iota_u+\nu\delta)\beta,
}
\]

and the rotational momentum and vorticity equations are simply

\[
\boxed{
\alpha_t+dB=-q_1\beta,
\qquad
\beta_t=-d q_1\beta .
}
\]

On the state itself the curvature return becomes the scalar positive source

\[
\boxed{
q_1^2\beta=\nu|\beta|^2.
}
\]

Since `q_1` annihilates scalars, the state current chain terminates after one
further step:

\[
\boxed{
\beta\xrightarrow{q_1}e
\xrightarrow{q_1}\nu|\beta|^2
\xrightarrow{q_1}0.
}
\]

There is an equivalent covariant-differential notation, but it is only a
repackaging of these same primitive operations.  Put

\[
\mathcal D_\theta
:=d+\frac{\theta}{\nu}\alpha\wedge .
\]

Then

\[
\boxed{
\mathcal D_\theta^2
=\frac{\theta}{\nu}\beta\wedge,
\qquad
q_\theta=\nu\mathcal D_\theta^*.
}
\]

In particular `alpha wedge alpha=0`, so `D_1 alpha=beta` and the full momentum
law modulo its exact Bernoulli gauge is

\[
\boxed{
\alpha_t+dB
=-\nu\mathcal D_1^*\mathcal D_1\alpha .
}
\]

This should not be read as importing an external gauge theory.  The
``connection'' is the velocity one-form itself and its curvature is the physical
vorticity already generated by NS.

The same self-return law gives a Gauss equation for the **actual acceleration**,
not only for the midpoint Poynting residual.  Since

\[
e^\sharp=-u\times\omega+\nu\operatorname{curl}\omega,
\qquad
u_t+\nabla B=-e^\sharp,
\]

one has

\[
\boxed{
-\nabla\cdot e^\sharp
+\frac{u}{\nu}\cdot e^\sharp
=|\omega|^2,
}
\]

or equivalently

\[
\boxed{
\left(\nabla\cdot-\frac{u}{\nu}\cdot\right)
(u_t+\nabla B)
=|\omega|^2.
}
\]

Thus the physical NS acceleration itself is a curvature-sourced current.  The
midpoint residual of Sections 41--43 is not the origin of the law; it is the
special member selected by the exact heat/Euler square.

Finally the same operator generates the vorticity PDE through

\[
\boxed{
\{d,q_1\}=\mathcal L_u+\nu L.
}
\]

Associativity therefore yields the exact current/generator commutator

\[
[\{d,q_1\},q_1]=d q_1^2-q_1^2d.
\]

On the closed state `d beta=0`,

\[
\boxed{
[\mathcal L_u+\nu L,q_1]\beta
=\nu\,d|\beta|^2.
}
\]

Spatial concentration is therefore a compulsory noncommutation of the current
with the generator which that same current creates.  This is an algebraic
self-frustration law, not a pointwise amplitude gap.


## 45. The midpoint current Hodge square contains Schrödinger, strain, Poynting and Gauss

The exact Hilbert squares of the previous sections select the canonical midpoint
between pure Hodge heat and the actual NS current:

\[
\boxed{
q:=q_{1/2}=\nu\delta+\frac12\iota_u,
\qquad
q^*=\nu d+\frac12u^\flat\wedge .
}
\]

Its curvature square is

\[
\boxed{
q^2=\frac{\nu}{2}(\beta\wedge)^*.
}
\]

Let

\[
\boxed{
L_u:=-\nu^2\Delta+\frac{|u|^2}{4}.
}
\]

The Hodge--Dirac square `q^*q+qq^*` has three complementary readings in three
spatial dimensions.  On scalars,

\[
\boxed{
qq^*\big|_{\Omega^0}=L_u.
}
\]

The drift cross term vanishes exactly because `div u=0`.  Thus the Schrödinger
operator of Section 42 is simply `nu^-2 qq^*` in degree zero; it was not selected
by an analyst after the Gauss identity.

On one-forms,

\[
\boxed{
H^{(1)}:=q^*q+qq^*
=L_u I+\nu S,
}
\]

whereas, after Hodge identification of two-forms with vectors,

\[
\boxed{
H^{(2)}:=q^*q+qq^*
=L_u I-\nu S.
}
\]

Hence physical strain is exactly the difference between two complementary
degrees of the same current Laplacian:

\[
\boxed{
2\nu S
=H^{(1)}-\star^{-1}H^{(2)}\star,
}
\]

while their sum is the scalar Schrödinger geometry,

\[
\boxed{
2L_uI
=H^{(1)}+\star^{-1}H^{(2)}\star,.
}
\]

Both Hodge squares are positive, so the same identity gives the two-sided
quadratic-form bound

\[
\boxed{
-L_u I\le\nu S\le L_u I.
}
\]

This is not a Sobolev estimate imposed from outside.  The symmetric velocity
gradient is the degree defect of one positive graded current square.

Now apply the degree-two identity to the actual vorticity two-form.  Closedness
gives

\[
q^*\beta=\frac12u^\flat\wedge\beta,
\]

whereas

\[
\boxed{
(q\beta)^\sharp
=\nu\operatorname{curl}\omega-\frac12u\times\omega
=-\frac12G.
}
\]

Thus the Section 41 completed square is literally the degree-two norm identity
of this same midpoint operator.  The degree-zero block is the Section 42
Schrödinger normal operator, and

\[
\boxed{
q^2\beta=\frac{\nu}{2}|\beta|^2
}
\]

is the Section 43 positive-source Gauss law.  Poynting, Gauss, Schrödinger and
strain are therefore not four independent mechanisms; they are graded
components of one midpoint current algebra.

The general parameter makes the same structure transparent:

\[
\boxed{
H_\theta^{(1)}
=\left(-\nu^2\Delta+\theta^2|u|^2\right)I
+2\theta\nu S,
}
\]

\[
\boxed{
H_\theta^{(2)}
=\left(-\nu^2\Delta+\theta^2|u|^2\right)I
-2\theta\nu S.
}
\]

The midpoint is distinguished because it is the member which simultaneously
appears in the exact Euler/heat reflection square.  No claim is made that the
positivity of these Hodge squares by itself supplies a supercritical coercive
estimate; attempts to extract a free negative `L4` vorticity term cancel inside
the exact current algebra.


## 46. The Gauss law polarizes to every strain direction

Section 43 used `b=omega`.  The primitive current square gives a stronger mixed
identity with the NS state `u` fixed and an arbitrary smooth vector field `b`.
Define

\[
\boxed{
\mathfrak G_\nu[u;b]
:=u\times b-2\nu\operatorname{curl}b.
}
\]

Then elementary curl calculus gives, with **no divergence-free hypothesis on
`b`**,

\[
\boxed{
\left(\nabla\cdot-\frac{u}{2\nu}\cdot\right)
\mathfrak G_\nu[u;b]
=\omega\cdot b.
}
\]

Indeed

\[
\nabla\cdot(u\times b)
=b\cdot\omega-u\cdot\operatorname{curl}b
\]

and the drift term restores the second summand exactly.  The Section 43 law is
only the diagonal member `b=omega`.

For divergence-free `b`, integration by parts also gives

\[
\langle b,Sb\rangle
=\langle u\times b,\operatorname{curl}b\rangle.
\]

Therefore the entire strain quadratic form has the exact Poynting factorization

\[
\boxed{
\nu\langle b,Sb\rangle
=\nu^2\|\operatorname{curl}b\|_2^2
+\frac14\|u\times b\|_2^2
-\frac14\|\mathfrak G_\nu[u;b]\|_2^2.
}
\]

After polarization, on the divergence-free form domain,

\[
\boxed{
4\nu A_u
=P\left(4\nu^2C^2+U_u^*U_u-\mathfrak G_u^*\mathfrak G_u\right)P,
\qquad A_u:=P S,
\qquad U_ub=u\times b.
}
\]

This is an operator-level self-frustration law: strain is heat geometry plus
local Lamb capacity minus the Gram energy of the same sourced current.

The mixed Gauss equation also lifts the scalar Schrödinger bound to an operator
quadratic-form inequality.  Put

\[
D_u:=\nabla\cdot-\frac{u}{2\nu}\cdot,
\qquad
M_\omega b:=\omega\cdot b.
\]

Then

\[
D_u\mathfrak G_u=M_\omega,
\qquad
D_uD_u^*
=-\Delta+\frac{|u|^2}{4\nu^2}.
\]

Consequently, whenever the inverse is defined in the relevant energy space,

\[
\boxed{
\mathfrak G_u^*\mathfrak G_u
\ge
M_\omega^*
\left(-\Delta+\frac{|u|^2}{4\nu^2}\right)^{-1}
M_\omega
}
\]

in quadratic-form order.  Every candidate strain direction with nontrivial
`omega.b` overlap therefore carries a canonical curvature floor.

There is a crucial sign guard.  If `b_+` and `b_-` lie in opposite curl-helicity
sectors, the even heat term vanishes from their cross pairing and

\[
\boxed{
4\nu (A_u)_{\rm odd}
=\left[P\left(U_u^*U_u-\mathfrak G_u^*\mathfrak G_u\right)P\right]_{\rm odd}.
}
\]

The helicity-odd projection of a positive Gram operator need not remain positive.
Hence the critical Reynolds operator cannot be bounded by simply transporting a
negative square through the helicity projection.  The exact result is a
**difference of two Gram currents** in the dangerous sector, not a free
contraction theorem.


## 47. Every Sobolev Hilbert square is an angle between two physical currents

The shifted Hilbert square of the curl-centered law can be written without a new
residual.  Define the two endpoint currents

\[
\boxed{
j_0:=\nu C^2u
}
\]

and

\[
\boxed{
j_1:=\nu C^2u-F_E=-u_t,
\qquad F_E=P(u\times Cu).
}
\]

Thus `j_0` is pure heat, `j_1` is the actual Leray-projected NS current, and their
midpoint is

\[
\boxed{
j_{1/2}:=\frac{j_0+j_1}{2}
=\nu C^2u-\frac12F_E
=-\frac12PG.
}
\]

For

\[
K_s=\|\Lambda^su\|_2^2,
\]

one has for every admissible real `s`

\[
\boxed{
\frac12K_s'
=-\frac1\nu
\langle j_0,j_1\rangle_{H^{s-1}}.
}
\]

The parallelogram law therefore gives

\[
\boxed{
K_s'
=\frac1{2\nu}
\left(
\|\Lambda^{s-1}F_E\|_2^2
-\|\Lambda^{s-1}PG\|_2^2
\right).
}
\]

Equivalently, the residual in the Section 54 Hilbert square is simply

\[
\boxed{
R_s
=\Lambda^{s+1}u-\frac1{2\nu}\Lambda^{s-1}F_E
=-\frac1{2\nu}\Lambda^{s-1}PG
=\frac1\nu\Lambda^{s-1}j_{1/2}.
}
\]

Thus Sobolev growth means exactly that the actual NS current turns obtuse to the
pure heat current in the shifted metric `H^{s-1}`.  The field measuring approach
to perfect backward-heat reflection is not an invented error term; it is the
same midpoint current whose unprojected representative satisfies the positive
Gauss law.

At the critical exponent `s=1/2`,

\[
\boxed{
K'
=\frac1{2\nu}
\left(
\|\Lambda^{-1/2}F_E\|_2^2
-\|\Lambda^{-1/2}PG\|_2^2
\right).
}
\]

The Euler energy null gives

\[
\langle u,F_E\rangle=0,
\qquad
\langle u,PG\rangle=-2\nu Z.
\]

Reading this in the critical dual pair yields the pressure-free compulsory tax

\[
\boxed{
\|\Lambda^{-1/2}PG\|_2^2
\ge\frac{4\nu^2Z^2}{K}.
}
\]

Hence

\[
\boxed{
K'
\le
\frac1{2\nu}\|\Lambda^{-1/2}F_E\|_2^2
-2\nu\frac{Z^2}{K}.
}
\]

The helicity null ` <Cu,F_E>=0 ` supplies a second compulsory projection of the
same residual.  More invariantly, for every `s` the two vectors

\[
\Lambda^{1-s}u,
\qquad
\Lambda^{1-s}Cu
\]

have fixed pairings `Z` and `<Cu,C^2u>` with `R_s`.  Orthogonal projection onto
their Gram plane therefore gives a canonical two-null lower bound on
`||R_s||^2`.  This is not a new stock: it is the normal component which energy
and helicity conservation forbid Euler from cancelling.

Write locally

\[
K_a:=\|\Lambda^a u\|_2^2,
\qquad
H_a:=\langle u,\mathsf J\Lambda^a u\rangle,
\qquad
H_3=\langle Cu,C^2u\rangle .
\]

The Gram-plane projection is explicit:

\[
\boxed{
\|R_s\|_2^2
\ge
\frac{
K_{2-s}Z^2-2H_{3-2s}ZH_3+K_{1-s}H_3^2
}{
K_{1-s}K_{2-s}-H_{3-2s}^2
},
}
\]

with the usual continuous/pseudoinverse reading when the Gram determinant
vanishes.  Since

\[
R_s=-\frac1{2\nu}\Lambda^{s-1}PG,
\]

the same formula is a two-null lower bound on the actual projected midpoint
current in the exact Sobolev metric.  At `s=1/2` it strictly sharpens the
single energy-null tax whenever the helicity-normal direction is independent.
No new invariant is introduced: this is only the orthogonal projection forced by
the two Euler null laws already present in NS.


## 48. Static Gauss coercivity cannot by itself close the Leray-projected gap

The current algebra also exposes a sharp no-go.  The scalar adjoint current is

\[
q_\theta^*\phi
=\nu d\phi+\theta u^\flat\phi.
\]

For this one-form to be divergence-free it must satisfy

\[
\nu\Delta\phi+\theta u\cdot\nabla\phi=0.
\]

On a periodic domain, or under sufficient decay, multiply by `phi` and use
`div u=0`.  One obtains

\[
\nu\|\nabla\phi\|_2^2=0.
\]

Therefore

\[
\boxed{
\phi=\text{constant}
}
\]

is the only scalar Gauss test living completely in the Leray-horizontal sector.
The constant test is precisely the invariant-normal lower bound.  Every
nonconstant Schrödinger refinement necessarily sees the exact/pressure sector as
well.  It is therefore invalid to turn the full static Schrödinger lower into a
coercive estimate for `PG` without additional dynamics.

That dynamics can also be written with no new object.  For the actual current
`e=q_1 beta`, differentiating and using `beta_t=-d e` gives

\[
\boxed{
e_t+(\mathcal L_u+\nu L)e
=\iota_{u_t}\beta+\nu d|\beta|^2.
}
\]

The curvature-gradient term is the same Jacobi defect found in Section 44.  Now
Hodge-decompose

\[
e=h+d\phi,
\qquad
h=Pe=-u_t.
\]

Because `L_u` and `L` commute with `d`, the exact sector remains exact under the
linear current generator.  Leray projection removes the explicit positive
source and leaves

\[
\boxed{
h_t+P(\mathcal L_u+\nu L)h
=-P\,\iota_h\beta.
}
\]

Thus the positive Gauss source is real but it does not inject a free sign into
the horizontal acceleration equation.  It constrains pressure/exact geometry
and feeds the projected current only through the same evolving `u,beta`
coefficients.  This is why persistence, rather than an instantaneous
Schrödinger gap, remains the theorem frontier.


## 49. The critical Reynolds operator is the helicity-odd reading of the same graded current

The critical curl law already shows that the dangerous generator is the
Poisson/Sylvester transform of physical odd strain.  Section 45 now identifies
that strain itself as

\[
\boxed{
S=\frac1{2\nu}
\left(H^{(1)}-\star^{-1}H^{(2)}\star\right),
\qquad A_u=P S\big|_{H_{div}}.
}
\]

Hence the critical midpoint Reynolds operator is not an independent spectral
mechanism.  It is the helicity-odd, log-curl-smoothed reading of the degree
imbalance of the same midpoint current complex.

Combining Section 46 with helicity parity gives the exact dangerous block
identity

\[
\boxed{
4\nu (A_u)_{\rm odd}
=\left[P\left(U_u^*U_u-\mathfrak G_u^*\mathfrak G_u\right)P\right]_{\rm odd}.
}
\]

The first Gram term is local Lamb capacity.  The second Gram term is a current
which, in every test direction `b`, must solve

\[
\boxed{
D_u\mathfrak G_u b=\omega\cdot b.
}
\]

Therefore a sustained supercritical Reynolds direction must do more than make
strain large.  Across the intrinsic Poisson depths of the critical transform it
must maintain a helicity-odd Gram imbalance while its reflected-current
component continuously avoids the curvature floor `omega.b`.

If it does not avoid that floor, `mathfrak G_u^* mathfrak G_u` returns a
compulsory current cost.  If it avoids the floor by reorienting, that turning is
produced by the same Cartan/Jacobi dynamics and the same strain which deforms the
material heat metric.  The material transverse-memory law then records the
history of that deformation.  This identifies a single remaining dynamical
question rather than a new case split:

> **Can a critical helicity-odd Gram imbalance persist for the infinite action
> required by escape while the same curvature/current algebra forces continual
> current turning and the same strain writes irreversible transverse heat
> memory?**

The present section does not answer that question.  It narrows the missing
regularity theorem to a persistence statement connecting three exact readings
of the same state-generated current law.


## 50. Minimal ontology after the graded-current collapse

The new identities can be summarized without promoting any representation to an
independent mechanism.  The physical state creates the one current operator

\[
\boxed{
q_1=\nu\delta+\iota_u,
}
\]

and the two primitive relations

\[
\boxed{
\beta_t=-d q_1\beta,
\qquad
q_1^2=\nu(\beta\wedge)^*.
}
\]

Everything in the new layer is generated from these and the canonical midpoint
between actual current and pure heat:

\[
\boxed{
\begin{aligned}
q_1^2
&\longleftrightarrow \text{vorticity curvature},\\
\{d,q_1\}
&\longleftrightarrow \text{vorticity evolution generator},\\
[\{d,q_1\},q_1]\beta
&\longleftrightarrow \nu d|\beta|^2\text{ current frustration},\\
q_{1/2}^*q_{1/2}+q_{1/2}q_{1/2}^*
&\longleftrightarrow \text{Schrödinger}\pm\text{strain},\\
q_{1/2}\beta
&\longleftrightarrow \text{Poynting midpoint current},\\
q_{1/2}^2\beta
&\longleftrightarrow \text{Gauss positive source},\\
S_{\rm odd}
&\xrightarrow{\text{intrinsic curl Poisson calculus}}
\text{critical Reynolds operator}.
\end{aligned}
}
\]

The interpretation is deliberately modest.  This is a **curved graded current
algebra**, not a claim that an external Lie group has solved Navier--Stokes and
not a global-regularity theorem.  Its value is ontological: structures which
previously appeared as separate Poynting, Gauss, Schrödinger, strain and critical
operator laws are now exact degree/projection readings of one current already
present in NS.

The remaining frontier is correspondingly strict.  A proof still has to show
that the state cannot sustain the required critical backward-heat alignment for
infinite productive action by concentrating/reorienting through the seams left
open by the static identities.  Combined with the productive-Fisher theorem on
the curl side, the exact closure target is

\[
\boxed{
\int^T
\frac{\kappa(0,t)^2}
{N(t)^2\,[E(t)Z(t)-K(t)^2]}
\,dt < \infty
}
\]

on every finite positive-energy smooth interval; critical mean-curl escape would
require the same integral to diverge.  The graded current algebra does not yet
prove this finiteness.  It identifies the only intrinsic dynamical loop from
which such a persistence theorem should come: current turning, helicity-odd Gram
imbalance, and material transverse heat memory are three readings of the same
state-generated law.  No owner, shell, cutoff, entropy currency or external
closure principle is introduced to bridge that gap.
