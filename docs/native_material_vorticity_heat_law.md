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
