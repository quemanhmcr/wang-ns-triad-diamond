# Research Ledger — Wang–Zahl-style rigidity programme for 3D Navier–Stokes

**Status date:** 2026-08-08
**Repository:** `quemanhmcr/wang-ns-triad-diamond`  
**Scope:** finite-dimensional / Gaussian-grain research programme; **not** a proof of global regularity.

This document is the central ledger for the current state of the programme. It records the mathematical architecture, theorem-level finite-dimensional statements, computationally supported constants, countermodels that changed the strategy, and the remaining PDE bridge. Individual derivations live in `docs/`; reproducible runs live in `recorded-results/`.

## 1. Guiding method

The programme deliberately follows an extremizer/rigidity workflow rather than trying to estimate the full Navier–Stokes equation directly:

\[
\boxed{\text{near-counterexample}
\to \text{near equality across scales}
\to \text{rigid coherent structure}
\to \text{incompatibility / quantitative loss}.}
\]

For 3D incompressible Navier–Stokes

\[
\partial_tu+(u\cdot\nabla)u=-\nabla p+\Delta u,
\qquad \nabla\cdot u=0,
\]

with scaling

\[
u_\lambda(x,t)=\lambda u(\lambda x,\lambda^2t),
\]

the intended contradiction mechanism is not “energy alone prevents infinite scale travel”. In fact a critical packet at frequency `N` has amplitude `N`, spatial size `N^{-1}`, lifetime `N^{-2}`, energy cost `N^{-1}`, and dissipation cost over its lifetime `N^{-1}`. Thus a single chain `N_q=2^qN_0` has summable total energy/dissipation cost:

\[
\boxed{\text{finite energy permits an infinite single-channel critical cascade}.}
\]

The programme therefore searches for a second ledger: an efficient cascade must either branch, reuse old structure, or remain in a highly rigid flat configuration, and every one of these behaviours should pay a quantitative cost.

## 2. Packet language and the sticky-cascade principle

For a Littlewood–Paley packet at frequency `N`, a model scale-invariant local mass is

\[
\mathsf M(N,x,t)
=N\int \chi(N(y-x))^2|u_N(y,t)|^2\,dy.
\]

Persistent packets also carry a localized dissipation charge. At a fixed time, or in spacetime for persistent packets, the elementary energy ledger gives schematically

\[
\sum_{P\in\mathcal P}\frac1{N_P}
\lesssim \frac{\|u_0\|_2^2}{\eta^2}.
\]

A single geometric chain remains summable, whereas large branching is expensive. This motivates the first structural dichotomy:

\[
\boxed{\text{attenuation or replication; otherwise rigidity}.}
\]

The word **sticky** is used for an efficient cascade that repeatedly reuses a nested ancestry branch rather than paying for many fresh independent packets.

## 3. Helical single-triad module

Using the helical basis

\[
\widehat u(k)=u^+(k)h^+(k)+u^-(k)h^-(k),
\qquad ik\times h^\pm(k)=\pm |k|h^\pm(k),
\]

a single triad admits an explicit coupling magnitude. With the triad-plane convention used in the code,

\[
\boxed{
|g_{kpq}|=
\frac{\Delta(k,p,q)}{2\sqrt2\,kpq}
|s_kk+s_pp+s_qq|,
}
\]

where `Δ` is the triangle area. The inviscid triad energy derivatives have the common phase factor structure

\[
\dot E_k=(s_pp-s_qq)R,
\quad
\dot E_p=(s_qq-s_kk)R,
\quad
\dot E_q=(s_kk-s_pp)R.
\]

Near equality forces amplitude balance and phase locking.

### Progress-weighted edge functional

The research diagnostic used throughout the finite model is

\[
\boxed{
\mathcal J(k,p,q;s_k,s_p,s_q)=
\log\frac q p\,
\frac{|s_kk-s_pp|}{q}|g_{kpq}|,
\qquad k\le p\le q.
}
\]

It measures transfer strength times logarithmic scale progress.  The functional
`mathcal J` itself is still a normalized finite-edge coefficient rather than the
full Navier--Stokes flux, but the logarithmic progress factor now has an exact
PDE-facing origin.  For one conservative ordered Fourier triad `k<=p<=q`, the
sharp-cutoff outward nonlinear energy flux satisfies

\[
\boxed{
\int_0^\infty \Pi_K^{(e)}\,\frac{dK}{K}
=-\dot E_k\log\frac pk
+\dot E_q\log\frac qp.
}
\]

Thus at the equal-parent extremizer `k=p` the lower segment vanishes exactly and
the surviving term is `dot E_q log(q/p)`.  Near the extremizer the omitted lower
segment is turned on by the same parent-imbalance variable that has cusp
stability.  The certified local box gives, even for the adverse maximizing child
helicity,

\[
\left|F_{low}\right|/F_{upper}<0.08296357712<0.1,
\qquad
F_{Mellin}\ge0.9F_{upper}.
\]

So the remaining PDE bridge is no longer the origin of the logarithm; it is the
global packet/triad synthesis, cancellation, localization, and construction of
the positive transfer weights used by the finite-dimensional ledger.

For symmetric opposite-helicity parents, `k=p=xq`,

\[
\mathcal J_{\rm sym}(x)
=
\frac{\sqrt{4x^2-1}}{4\sqrt2\,x}
\log\frac1x.
\]

The critical point solves

\[
\boxed{-\log x=4x^2-1,}
\]

giving

\[
r_*=x_*\approx0.6109041015867660,
\qquad
R_*=r_*^{-1}\approx1.63691813068957,
\]

\[
\gamma_*:=\log R_*\approx0.492815285342135,
\qquad
J_*\approx0.100110175856189.
\]

The local common-scale curvature is negative,

\[
J''(r_*)\approx-4.40210953306117.
\]

The decimal for `r_*` recorded in the first version of this ledger was stale:
Action `31156944917` rejected the old bracket by rigorous Arb arithmetic.  The
corrected value above is consistent with the exact equation and with the
certified rational bracket

\[
0.61090410158<r_*<0.61090410160.
\]

The full sign/domain optimizer is now certified rather than merely numerical.
For child magnitude one and ordered parents `0<x<=y<1`, exact sign elimination
gives

\[
|s_xx-s_yy|\,|s_xx+s_yy+s_z|
\le (x+y)(1+y-x),
\]

with equality for the opposite-parent-helicity orbit.  Arb branch-and-bound on
the resulting two-variable envelope proves that the symmetric critical point is
the unique global maximum in the normalized forward-triad domain.

## 4. The cusp: the extremizer is anisotropically rigid

Write the normalized parent magnitudes as

\[
x=r+d,\qquad y=r-d.
\]

Because `\mathcal J` contains `\max(x,y)`, it is nonsmooth at the maximizing line `x=y`. The transverse scale-imbalance deficit is therefore linear, not quadratic:

\[
\boxed{
1-\frac{J(r+d,r-d)}{J_*}
\gtrsim c|d|+c'(r-r_*)^2.
}
\]

The one-sided normalized slopes at the optimum are

\[
a_+=2+\frac1{r_*\log(1/r_*)}\approx5.322,
\]

\[
a_-=\frac1{r_*\log(1/r_*)}-2\approx1.322.
\]

Thus a near-extremizer has thickness `O(ε)` in the parent-scale imbalance direction and `O(√ε)` in smooth tangent directions. This motivates the term **triad grain**: the near-extremal set is an anisotropic thin object rather than a point triad.


### Certified mixed single-edge stability

Use the log-scale coordinates already seen by the Hodge module,

\[
u=|\ell_p-\ell_q|,
\qquad
v=\ell_c-\frac{\ell_p+\ell_q}{2}-\gamma_*.
\]

Action `31157463384` gives an inclusion-arithmetic certificate at 160-bit Arb
precision.  On the whole near-extremal rectangle

\[
0\le u\le\frac{2}{25},
\qquad |v|\le\frac{2}{25},
\]

the normalized deficit satisfies

\[
\boxed{
\operatorname{Def}_e
:=1-\frac{J_e}{J_*}
\ge \frac1{50}u+v^2.
}
\]

The certified lower bounds used by the proof are

\[
\partial_u\operatorname{Def}_e>0.0468113935>\frac1{50},
\qquad
\partial_v^2\operatorname{Def}_e>6.805681556>2
\]

on the relevant intervals.  Since `u<=2/25`, the linear cusp dominates
`u^2/4`; together with the exact residual identity

\[
r_p^2+r_q^2=\frac{u^2}{2}+2v^2,
\]

this yields the theorem-level edge-to-Hodge conversion

\[
\boxed{
\operatorname{Def}_e\ge\frac12(r_p^2+r_q^2).
}
\]

Outside the local rectangle, the same Action certifies the uniform global gap

\[
\boxed{
\operatorname{Def}_e\ge\frac1{100}.
}
\]

The compact global certificate used `637` gap boxes plus `93` boxes absorbed by
the local theorem, with maximum subdivision depth `15`.  The accompanying
`100,000` local plus `100,000` global random stress samples are regression tests
only and are not part of the proof.

For normalized transfer weights `w_e`, this is already in the form required by
the Hodge ledger.  If the retained good-edge network has arc conductance `w_e`
on each of the two parent-to-child arcs of triad `e`, then

\[
1-R_{\rm block}
=\sum_e w_e\operatorname{Def}_e
\ge\frac12\sum_e w_e(r_{e,p}^2+r_{e,q}^2)
\ge\frac12\mathcal E_H.
\]

Therefore

\[
\boxed{-\log R_{\rm block}\ge\frac12\mathcal E_H,}
\]

and any Hodge-cost branch with `\mathcal E_H\ge h_H>0` now has the certified
finite-dimensional block cost `c_{0,H}=h_H/2`.  A bad-edge set of total transfer
weight `\beta` instead pays at least `\beta/100` before any Hodge routing is
invoked.

For the exact nonflat three-triad reuse motif, unit arc weights have
`\mathcal E_H=\gamma_*^2/5`.  With equal normalized triad transfer weights
`w_e=1/3`, there is now a completely closed local/global dichotomy.  If all
three triads are in the local rectangle, the inherited arc weights give

\[
-\log R_{\rm block}\ge\frac{\gamma_*^2}{30}
\approx0.00809556352.
\]

If any one triad is outside, its weight `1/3` and the global gap give
`-\log R_{\rm block}\ge1/300`.  Hence the unconditional motif cost is

\[
\boxed{
c_{0,\rm motif}
\ge\min\!\left\{\frac{\gamma_*^2}{30},\frac1{300}\right\}
=\frac1{300}>0.
}
\]

This is a concrete theorem-level positive `c_0` for that normalized reuse block,
not yet a claim that the same number is the universal master constant across
all costly branches.

### Smooth log-scale cocycle and physical SGS transfer weights

The sharp all-scale Mellin identity is **not** itself the correct block edge
functional. There is now a certified countermodel. At

\[
x=\frac{13}{40},\qquad y=\frac{17}{20},\qquad
(s_x,s_y,s_q)=(-,+,-),
\]

the full all-scale Mellin coefficient is Arb-certified to be

\[
0.157176722626149\ldots>\frac32J_*.
\]

The extra gain is dominated by transfer across the lower scale segment. Thus
replacing `mathcal J` by `int Pi_K dK/K` would destroy the equal-parent grain.
This is a genuine countermodel, not a numerical optimizer artifact.

The correct replacement is a **block-local weak scale potential**. Write
`t=log K`. For any conservative finite interaction with modal log frequencies
`ell_i` and rates `dot E_i`, let

\[
\Pi^0(t)=-\sum_i\dot E_iH(t-\ell_i).
\]

If `rho_delta` is an even probability kernel on log scale, with CDF
`Psi_delta`, the graded spectral flux satisfies the exact convolution identity

\[
\boxed{\Pi^\delta=\rho_\delta*\Pi^0.}
\]

Hence smoothing only redistributes flux in `log K` and preserves the all-log
moment:

\[
\boxed{
\int\Pi^\delta(t)dt=\sum_i\dot E_i\ell_i.
}
\]

More importantly, the weak scale-potential identity

\[
\boxed{
\int\Pi^0(t)\Phi'(t)dt=\sum_i\dot E_i\Phi(\ell_i)
}
\]

allows one to isolate a generation boundary. For an even compact log-kernel,
define `Phi' = 2 Psi_delta(t-tau)`. Then `Phi` is exactly flat below
`tau-delta` and exactly `2(t-tau)` above `tau+delta`.

For a forward core of edges with positive child transfer `T_e`, top-parent and
child logs `(p_e,q_e)`, define the **transfer-weighted midgap**

\[
\boxed{
\tau=\frac{\sum_eT_e(p_e+q_e)/2}{\sum_eT_e}.
}
\]

If a common smooth moat separates all parents and children,

\[
p_e\le\tau-\delta,\qquad q_e\ge\tau+\delta,
\]

then all centering errors cancel exactly and

\[
\boxed{
2\int_\tau^\infty\Pi^\delta_{\rm core}(t)dt
=\sum_eT_e(q_e-p_e)
=\sum_eT_e\log\frac{Q_e}{P_e}.
}
\]

Thus one **single smooth physical scale cut** reproduces the complete upper
progress ledger of a packet block exactly. There is no sharp-to-smooth loss.

The common moat itself now has a certified numerical margin. If parent and
child log shells each have halfwidth `2/25` and the smooth transition halfwidth
is `1/20`, Arb gives

\[
\boxed{
\frac{\gamma_*}{2}-2\frac{2}{25}-\frac1{20}
>\frac9{250},
}
\]

with enclosure `0.0364076427...`.

Cancellation is also no longer an unstructured PDE ambiguity. If `A_e` are
unsigned transfer capacities, `m_e=J_e/J_*` and `c_e in [-1,1]` is the
phase/orientation factor, then for

\[
R=\sum_ew_em_ec_e,\qquad w_e=A_e/\sum A,
\]

there is the exact polarization identity

\[
\boxed{
1-R=\sum_ew_e(1-m_e)+\sum_ew_em_e(1-c_e).
}
\]

Near saturation therefore forces both multiplier rigidity and phase/backscatter
rigidity in the same positive capacity measure. Taking the pointwise good
threshold `eta_0=10^{-4}`, the good core has capacity mass at least
`1-10^4(1-R)`. On this core the mixed single-edge theorem gives

\[
\left|\log(q/p)-\gamma_*\right|\le\frac1{80}.
\]

Since actual upper progress satisfies `A_e J_* r_e=T_e log(q/p)`, Arb certifies
that normalized positive child-transfer weights and normalized capacity weights
have pointwise Radon--Nikodym distortion less than

\[
\boxed{\frac{53}{50}},
\]

with enclosure `1.0521543606...`.

Consequently the Hodge branch can now be stated in **physical transfer weights**.
If `epsilon=1-R`, then either

\[
\epsilon\ge\frac1{20000},
\]

or the positive core has at least half the capacity mass and

\[
\boxed{
\epsilon\ge\frac{25}{106}\,\mathcal E_H^{\rm phys}.
}
\]

Thus a physical-transfer Hodge threshold `E_H^{phys}>=h_H` has the explicit
finite-packet cost

\[
\boxed{
c_{0,H}^{\rm phys}
=\min\left\{\frac1{20000},\frac{25}{106}h_H\right\}>0.
}
\]

Finally this graded spectral flux is a physical coarse-grained quantity. Choose
a self-adjoint convolution filter `G_t` with
`|Ghat_t|^2=Psi_delta(t-log|xi|)`. With

\[
\tau_t(u,u)=\overline{u\otimes u}-\bar u\otimes\bar u,
\qquad
\Pi_t^{\rm SGS}(x)=-\nabla\bar u:\tau_t(u,u),
\]

periodicity or spatial decay gives

\[
\boxed{
\int\Pi_t^{\rm SGS}(x)dx=\Pi^\delta(t).
}
\]

The Leray projector changes no global work because the filtered velocity is
divergence free. Pressure is therefore spatial transport, not global
interscale work; after packet-window localization it reappears only as boundary
work in the transition moat. The dealiased pseudospectral regression confirms
these identities to roughly `5e-16` relative error. A naive aliased regression
failed first, correctly exposing fake numerical energy transfer if Galerkin
conservation is not respected.

The localized pressure branch can also be made exact. For a smooth moving
window `chi`, the resolved energy density `e=|U|^2/2` satisfies

\[
\partial_t e+\nabla\cdot J=-\Pi-\nu|\nabla U|^2,
\qquad
J=(e+P)U+RU-\nu\nabla e.
\]

Define the combined pressure--flux work

\[
\boxed{G=\Pi+\nabla\cdot(PU).}
\]

Then pressure cancels **exactly** from the localized leakage ledger:

\[
\boxed{
\int\!\!\int\chi G
+\nu\int\!\!\int\chi|\nabla U|^2
=E_\chi(t_0)-E_\chi(t_1)+\widetilde L_\chi,
}
\]

where

\[
\widetilde L_\chi
=\int\!\!\int e\,\partial_t\chi
+\int\!\!\int\nabla\chi\cdot
\left[eU+RU-\nu\nabla e\right].
\]

For a chain of matching window endpoints and nonincreasing nonnegative weights
`a_j`, summation by parts gives the depletion estimate

\[
\boxed{
\sum_ja_jW_{j,+}+\sum_ja_jD_j
\le a_0E_0+\sum_ja_j(L_j)_++\sum_ja_jW_{j,-}.
}
\]

Thus forward combined work is paid by initial localized energy, positive window
leakage, or negative combined work/backscatter; pressure cannot reset this
budget.

There is also an exact single-window pressure-cancellation dichotomy. Let

\[
S=\int\chi\Pi\ge0,
\qquad W=\int\chi G.
\]

Then

\[
\boxed{
W\ge\frac S2
\quad\text{or}\quad
\left|\int PU\cdot\nabla\chi\right|\ge\frac S2.
}
\]

In the second branch, with `A=supp grad chi`, Holder and Young imply

\[
\boxed{
\int_A\left(|U|^3+|P|^{3/2}\right)
\ge\frac{S}{2\|\nabla\chi\|_\infty}.
}
\]

The right side is scale critical under Navier--Stokes scaling. Hence pressure
can cancel a forward SGS transfer only by creating a definite CKN-type critical
charge in the spatial moat. The pressure problem has therefore been converted
from an uncontrolled nonlocal sign into a concrete **combined-work or critical
annular charge** alternative.

### Crossing extraction, one-shot profiles, and pressure multipoles

The PDE-facing bridge has now been reduced further in three directions.

First, a **common spectral moat is forced by crossing geometry**; it is not an
extra Gaussian hypothesis.  On the signed-good core at `eta_0=10^{-4}`, the
single-edge theorem gives

\[
u\le \frac1{200},\qquad
\left|g-\gamma_*\right|\le\frac1{80}.
\]

If a positive child-transfer edge crosses one physical reference cut `tau_0`,
then its parent/child midgap lies in an interval of length
`gamma_*+1/80`.  Splitting this interval into four equal bins gives one bin
carrying at least one quarter of the positive transfer, and (independently) one
bin carrying at least one quarter of the transfer-weighted Hodge numerator.  In
that bin both parents and the child lie in common shells of halfwidth

\[
\boxed{
\sigma_{cross}
=\frac{\gamma_*}{8}+\frac{5}{8}\frac1{80}+\frac1{200}.
}
\]

Action `31164771160` certifies at 160-bit Arb precision

\[
\sigma_{cross}=0.07441441067\ldots<\frac2{25},
\]

and, with the existing smooth transition halfwidth `1/20`,

\[
\boxed{
\frac{\gamma_*}{2}-2\sigma_{cross}-\frac1{20}
>\frac1{25}.
}
\]

Thus an arbitrary crossing good core contains a one-quarter subcore on which
the exact smooth-midgap SGS identity applies.  A conservative PDE-facing Hodge
coefficient after this scale extraction is therefore

\[
\boxed{
\frac{25}{424},
\qquad
c_{0,H}^{cross}
=\min\left\{\frac1{20000},\frac{25}{424}h_H\right\}>0.
}
\]

Second, the **within-block Gaussian extraction need only be one-shot**.  A
frequency-localized finite-energy block automatically lies in Fourier
`L^{3/2}` because, when `|Omega_N|<=C_Omega N^3`,

\[
\|f\|_{3/2}
\le C_\Omega^{1/6}N^{1/2}\|f\|_2.
\]

Weighted near equality implies ordinary Young near equality, so Michael
Christ's near-extremizer theorem supplies one affine Gaussian extremizing triple
at any prescribed small `L^{3/2}` distance `epsilon_G`, once the block deficit
is sufficiently small.  The new deterministic reduction is exact: replacing
three unit inputs by profiles at distance `epsilon_G` changes normalized
weighted transfer by at most

\[
\boxed{
3\epsilon_G+3\epsilon_G^2+\epsilon_G^3.
}
\]

Moreover, because the original role is supported in the finite frequency block,
the Gaussian profile has at least `1-epsilon_G` of its `L^{3/2}` norm there and
therefore carries critical energy

\[
\boxed{
N\|G\|_2^2
\ge C_\Omega^{-1/3}(1-\epsilon_G)^2.
}
\]

So no global `||uhat||_{3/2}` hypothesis and no infinite list of profiles inside
one efficient block are required.  Christ's inverse modulus is an external
analytic theorem and is not assigned a fake numerical constant by the
repository.

Third, pressure cancellation now has a **multipole collision mechanism**.  The
Newtonian pressure kernel is

\[
K_{ij}(z)=\frac1{4\pi}
\left(3\frac{z_i z_j}{|z|^5}-\frac{\delta_{ij}}{|z|^3}\right).
\]

Since `div U=0`,

\[
\int U\cdot\nabla\chi=0,
\]

so a constant pressure mode performs no boundary work.  For a far source one
may subtract `K_ij(x_0-y)` and use

\[
\sum_{ij}|\nabla K_{ij}(z)|\,|u_i u_j|
\le 10|z|^{-4}|u|^2.
\]

The far pressure therefore enters through a dipole moment

\[
\boxed{
\mathfrak P_{far}
=\sum_{n\ge n_0}2^{-4n}\sum_a\mu_{n,a}.
}
\]

Three-dimensional packet packing grows only as `2^{3n}`.  Hence, if every
packet has critical mass at most `mu_*`, the far pressure tail is bounded by a
constant times

\[
\mu_*\sum_{n\ge n_0}2^{-n},
\]

which is summable.  The spare exponent is exactly `4-3=1`; without the constant
pressure cancellation there would be no spare power.

For the local pressure source, Calderon--Zygmund plus band-limited Bernstein
gives the critical-mass implication

\[
\boxed{
r^{-1}\|V\|_2^2
\ge
\frac1{C_B^2}
\left(\frac{\rho}{(r\|\nabla\chi\|_\infty)C_R}\right)^{2/3}
}
\]

when normalized local pressure work is at least `rho`.  In the abstract finite
packet model, if absence of a fresh packet implies

\[
W_{cancel}\le C_{near}\mu_*^{3/2}+C_{far}\mu_*,
\]

then cancellation at least `rho` forces

\[
\boxed{
\mu_{fresh}
\ge
\min\left\{
\left(\frac{\rho}{2C_{near}}\right)^{2/3},
\frac{\rho}{2C_{far}}
\right\}>0.
}
\]

Action `31164911526` stress-tested this packet collision theorem on `50,000`
adversarial shell configurations after the exact kernel estimate; the numerical
threshold printed by that run belongs only to its chosen test constants, not to
a universal PDE theorem.

### Certified narrow-shell mass and smooth-SGS symbol freezing

The crossing shell also makes the frozen-time critical-mass bridge quantitative.
For the certified radial log shell

\[
\left|\log\frac{|\xi|}{N}\right|\le\frac2{25},
\]

the full spherical shell obeys

\[
|\Omega_N|
\le
\frac{4\pi}{3}
\left(e^{6/25}-e^{-6/25}\right)N^3.
\]

If Christ's theorem has supplied a Gaussian profile at `L^{3/2}` distance
`epsilon_G=1/100`, then at least `99/100` of the profile's normalized
`L^{3/2}` mass remains in this shell.  Action `31165654509` certifies at
160-bit Arb precision

\[
C_\Omega
=2.0299769094616733\ldots
\]

and consequently

\[
\boxed{
N\|G\|_2^2
\ge
0.7740577380943306\ldots
>\frac34.
}
\]

This constant is scale invariant.  The statement is conditional only on the
external inverse-Young conclusion `epsilon_G<=1/100`; the shell-volume and
critical-mass implication themselves are interval-certified.

The frequency-multiplier part of the smooth-SGS packetization can also be made
summable without invoking a black-box decomposition theorem.  On the certified
compact crossing shells, after factoring the dimensional power of `N`, the
smooth filter/Leray/helical SGS symbol is a smooth scale-covariant function of
dimensionless frequencies.  Let `L_M` be its Lipschitz constant on that compact
set.  Freeze the normalized symbol on relative frequency cells of diameter
`h`.  Then

\[
\|M-M_h\|_\infty\le L_Mh.
\]

Sharp Young at the critical exponents gives the exact deterministic estimate

\[
\boxed{
|T_M(f,g,h)-T_{M_h}(f,g,h)|
\le
A_3L_Mh
\|f\|_{3/2}\|g\|_{3/2}\|h\|_{3/2},
\qquad
A_3=\left(\frac{\sqrt3}{2}\right)^3.
}
\]

Thus frozen frequency cells provide the finite/countable interaction edges and
the multiplier-freezing loss enters directly into the existing summable error
ledger.  For a uniform scale-covariant `L_M`, a schedule

\[
h_j=(j+3)^{-2}
\]

is summable; more generally choose `h_j` so that
`A_3 L_{M,j}h_j<=epsilon_j` with `sum epsilon_j<infinity`.
Action `31165838379` passed the full test suite plus `50,000` compact-cell
stress checks.  The random checks are regression evidence only; the theorem is
the sup-symbol bound followed by sharp Young.

This removes the **frequency multiplier freezing** part of the PDE
construction.  What remains is spatial/time packet coherence: localized
partition synthesis, moving windows, overlap/leakage, and preservation of the
frozen interaction ledger through a packet lifetime.

### Affine spacetime grain dynamics: strain, curvature, and objective-source collisions

The frozen-time profile theorem can now be propagated through a genuine
**affine/Kelvin Gaussian model** over one parabolic packet lifetime.  The first
lesson is that absolute strain is the wrong defect: sharp Young has affine
Gaussian symmetries, and a common rigid rotation or common planar scalar strain
can transport the whole extremal triad without changing its shape.  The charged
quantity is the **trace-free symmetric strain restricted to the triad plane**.

For an affine resolved velocity `U(x)=U_0+A(x-X)` a Kelvin carrier satisfies

\[
\dot k=-A^Tk,
\qquad
\frac d{dt}\log|k|=-\widehat k^TS\widehat k,
\qquad S=\frac{A+A^T}{2}.
\]

At the symmetric extremizer, write `phi=theta_*/2`, so
`cos(phi)=1/(2r_*)`.  If the trace-free planar strain is

\[
D=\begin{pmatrix}\delta&\beta\\ \beta&-\delta\end{pmatrix},
\]

then the signed Hodge shape coordinates obey exactly

\[
\dot u=4\beta\cos\phi\sin\phi,
\qquad
\dot v=-2\delta\sin^2\phi.
\]

Consequently

\[
\boxed{
\frac12\dot u^2+2\dot v^2
\ge
4\sin^4\phi\,\|D\|_F^2.
}
\]

Action `31168393589` certifies at 160-bit Arb precision

\[
\boxed{
4\sin^4\phi
=0.4359294243\ldots>\frac{43}{100}.
}
\]

Thus the infinitesimal kernel is exactly the planar-conformal direction; planar
shear and anisotropic stretching cannot be invisible to an extremal triad.

The fixed orientation of that plane is in fact unnecessary for scalar triad
geometry.  Put the two parent carriers into the columns of `K in R^{3x2}` and
set

\[
G=K^TK.
\]

For any common three-dimensional carrier law `Kdot=-B K`,

\[
\boxed{
\dot G=-K^T(B+B^T)K.
}
\]

All three side lengths, and therefore the magnitude of the helical single-edge
multiplier, are functions of this intrinsic `2x2` Gram matrix.  If `E` is any
orthonormal frame of the instantaneous (possibly tilting) triad plane, the
charged driver is simply

\[
D=\left(E^T\operatorname{sym}(B)E\right)^0.
\]

Orthogonal invariance gives the same extremal coercivity

\[
\boxed{
\frac12\dot u^2+2\dot v^2
>\frac{43}{100}\|D\|_F^2
}
\]

for an arbitrarily oriented and tilting 3D plane.  Action `31170015795` checked
`50,000` random 3D plane/driver/Gaussian configurations; the exact Gram identity
is the theorem, while the random run is regression evidence.  Extrinsic common
plane tilt is therefore a gauge for scalar side-length geometry.  The genuinely
3D remainder is coherent helical-polarization/phase transport and the spatial
packet frame.

For a frozen trace-free planar strain with eigenvalues `+/-d`, starting from the
exact extremal triad, if

\[
dT\le\frac1{25},
\]

the whole path remains in the already-certified local single-edge box and

\[
\mathcal H(t):=\frac{u(t)^2}{2}+2v(t)^2
\ge\frac35d^2t^2.
\]

Using `Def>=H/2` gives the packet-lifetime theorem

\[
\boxed{
\frac1T\int_0^T\operatorname{Def}(t)\,dt
\ge\frac1{10}(dT)^2.
}
\]

The Gaussian envelope itself has an exact affine advection--diffusion evolution.
For Fourier precision `P` and spectral peak `kappa`,

\[
\boxed{
\dot P=AP+PA^T+2\nu I,
\qquad
\dot\kappa=-A^T\kappa-2\nu P^{-1}\kappa.
}
\]

The dual center `b=P kappa` satisfies

\[
\boxed{\dot b=Ab,}
\]

so a common Gaussian precision preserves the three-wave resonance relation.
For incompressible `tr A=0`,

\[
\boxed{
\frac d{dt}\log\det P=2\nu\operatorname{tr}(P^{-1})\ge0.
}
\]

At an isotropic Gaussian instant `P=pI`, viscosity contributes only a scalar to
the carrier-length driver, so it is first-order neutral for the triad shape.
For the vector Kelvin amplitude, pressure/Leray projection preserves
transversality and does no direct transverse work:

\[
\dot a=-Aa+2k\frac{k\cdot Aa}{|k|^2}-\nu|k|^2a,
\]

\[
\boxed{
\frac d{dt}|a|^2=-2a\cdot Sa-2\nu|k|^2|a|^2.
}
\]

The natural physical-space window must deform affinely as well.  If

\[
\dot X=U(X,t),
\qquad
\dot F=(\nabla U)(X,t)F,
\qquad
\chi(x,t)=\chi_0(F^{-1}(x-X)/R),
\]

then the affine Taylor velocity cancels exactly:

\[
\boxed{(\partial_t+U_{aff}\cdot\nabla)\chi=0.}
\]

The first uncancelled advective error is curvature.  If
`||nabla^2 U||<=H` on the window,

\[
\boxed{
|(\partial_t+U\cdot\nabla)\chi|
\le
\frac H2 R\|F\|^2\|F^{-1}\|\|\nabla\chi_0\|_\infty.
}
\]

This is the correct moving-frame replacement for a merely translated packet.

#### Coherent strain or objective-strain action

A genuine Navier--Stokes grain sees a time-dependent gradient.  Remove the local
rigid rotation by working in the frame generated by the antisymmetric part of
`A`.  Let `D(t)` be the trace-free symmetric strain on an invariant extremal
triad plane, with `D(0)` eigenvalues `+/-d`.  If

\[
\sup_{0\le t\le T}\|D(t)-D(0)\|_{op}\le\frac d{20},
\qquad
 dT\le\frac1{30},
\]

Action `31168671304` certifies that the entire path stays in the local box with
`|u|,|v|<=7/100`, and

\[
\boxed{
\frac1T\int_0^T\operatorname{Def}(t)dt
\ge\frac1{24}(dT)^2.
}
\]

If this five-percent coherence fails, continuity gives a definite objective
strain action

\[
\boxed{
\int_0^T\|\dot D(t)\|_{op}dt>\frac d{20}.
}
\]

This failure is not a free time-dependence error.  For a smooth incompressible
Navier--Stokes solution,

\[
\boxed{
D_tA+A^2=-\nabla^2p+\nu\Delta A,
\qquad A=\nabla u.
}
\]

Writing `A=S+Omega`, the objective/corotational strain derivative is

\[
\boxed{
\mathring S
=-S^2-\Omega^2+[S,\Omega]-\nabla^2p+\nu\Delta S.
}
\]

Therefore strain dephasing must be generated by self-stretching/vorticity,
pressure Hessian, or viscous smoothing.

If `d=sigma N^2` and `T=cN^{-2}`, five-percent coherence failure forces at least
one of these three normalized source channels to have average size

\[
\boxed{
\rho_{src}\ge\frac{\sigma}{60c}.
}
\]

In the stated band-limited packet model, the quadratic and viscous channels
force critical `L^2` mass by Bernstein.  The far pressure-Hessian kernel is even
more local than the pressure-work kernel: `\nabla^2 K=O(|x|^{-5})`, so three-
dimensional packet packing leaves the spare exponent

\[
\boxed{5-3=2,}
\]

hence a summable `2^{-2n}` far tail.  Action `31169128097` stress-tested these
source-collision formulas on `50,000` packet configurations.  The local
pressure-Hessian coefficient remains a hypothesis of the continuum packet
realization.

### Helical spin connection, phase holonomy, and 3D strain tomography

The genuinely three-dimensional polarization remainder can now be separated into
**gauge transport**, **helicity conversion**, and **reuse holonomy**.  These are
not the same phenomenon.

For the standard local spherical helical gauge

\[
h_s=(e_\theta+i s e_\phi)/\sqrt2,
\qquad i k\times h_s=s|k|h_s,
\]

the Berry connection and curvature are exactly

\[
\mathcal A_s=i h_s^*dh_s=s\cos\theta\,d\phi,
\qquad
\mathcal F_s=-s\sin\theta\,d\theta\wedge d\phi,
\]

so the helicity line bundle has Chern number

\[
\boxed{c_1=-2s.}
\]

Thus no globally smooth phase convention exists.  It would nevertheless be
incorrect to charge the Berry phase of one rotating triad as a cascade defect.
A nondegenerate triad supplies its own physical gauge.  If `n` is its oriented
normal, set

\[
\boxed{
h_s(k;n)=\frac{n\times\widehat k+i s n}{\sqrt2}.
}
\]

This gauge is exactly `SO(3)` covariant:

\[
h_s(Rk;Rn)=Rh_s(k;n).
\]

In the reality-compatible forward-triad convention, the Waleffe coupling in
this gauge is purely imaginary.  Hence away from a coupling zero its phase is
locally constant (`+/- pi/2`).  Action `31171018560` passed the full test suite
plus `50,000` random three-dimensional rotation/transition/coupling checks.
This is a theorem-level **countermodel to a naive single-triad Berry tax**:
rigid rotation may accumulate chart-dependent Berry phase, but it is an exact
symmetry and does not itself reduce transfer.

The geometric phase becomes physical when one Fourier mode is reused by two
triads with different planes.  If their normals differ by signed dihedral angle
`psi` about the common carrier `k`, then

\[
\boxed{
h_s(k;n_2)=e^{-is\psi}h_s(k;n_1).
}
\]

This spin-one transition function is the local incidence connection of the
reuse network.

For the four-edge reuse diamond

\[
a+b=m,\qquad m+c=d,
\]
\[
b+c=n,\qquad a+n=d,
\]

write the forward phase residual of an edge `x+y=z` as

\[
\delta_e
=\arg g_e-\theta_x-\theta_y+\theta_z-\tau_e,
\]

where `tau_e` is the sign-dependent forward target.  All six modal phases cancel
exactly:

\[
\boxed{
\delta_{abm}+\delta_{mcd}-\delta_{bcn}-\delta_{and}
\equiv H_\phi\pmod{2\pi}.
}
\]

The same `H_phi` is reconstructed from the four constant triad-normal coupling
phases plus six spin-weighted dihedral transition phases.  It is invariant under
an arbitrary rigid rotation of the whole diamond.

There is a sharp four-edge phase-lock inequality.  For principal
`|H_phi|<=pi`,

\[
\boxed{
\sum_{i=1}^4(1-\cos\delta_i)
\ge
4\left(1-\cos\frac{|H_\phi|}{4}\right).
}
\]

The equal lift `(H/4,H/4,-H/4,-H/4)` is the sharp minimizer.  Combining this with
the existing positive polarization identity, if each diamond edge has capacity
weight at least `beta_phi`, each multiplier satisfies `m_e>=1-10^-4`, and

\[
|H_\phi|\ge\frac15,
\]

then Arb certifies

\[
4(1-10^{-4})\left(1-\cos\frac1{20}\right)
=0.0049984585\ldots>\frac1{250},
\]

hence the finite-packet phase branch pays

\[
\boxed{
D_{phase}\ge\frac{\beta_\phi}{250}.
}
\]

Action `31171127537` certified this coefficient and stress-tested `50,000`
nondegenerate diamonds.  A later implementation exposes the same holonomy as an
explicit sum of the six spin-dihedral transitions; that identity is rerun in the
final integration below.

The local Kelvin polarization dynamics has a complementary structure.  Let `E`
be a transverse frame of `k^perp`.  Its spin about `k` is free.  Choosing

\[
E^T\dot E=-\operatorname{skew}(E^T\nabla u\,E)
\]

removes rigid transverse rotation.  In this **objective polarization frame** the
linear generator is exactly

\[
\boxed{
-\operatorname{sym}(E^T\nabla u\,E)-\nu|k|^2I.
}
\]

For trace-free transverse strain

\[
D=\begin{pmatrix}\delta&\beta\\ \beta&-\delta\end{pmatrix},
\]

the circular/helical representation is

\[
\boxed{
D_{hel}=
\begin{pmatrix}
0&\delta-i\beta\\
\delta+i\beta&0
\end{pmatrix}.
}
\]

Thus first-order nonconformal strain is **opposite-helicity conversion**, not an
unavoidable diagonal Berry-phase drift.  Frozen `D=diag(d,-d)` sends a pure
positive-helicity amplitude to

\[
u_+(t)=\cosh(dt),\qquad u_-(t)=-\sinh(dt).
\]

Time-varying strain orientations generate local geometric phase through
noncommutativity:

\[
\boxed{
[D_1,D_2]
=2(\delta_1\beta_2-\beta_1\delta_2)
\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
}
\]

Under the five-percent coherence hypothesis, the second Magnus generator obeys

\[
\|\Omega_2\|_{op}
\le(\varepsilon+\varepsilon^2/2)(dT)^2.
\]

At `epsilon=1/20`, `dT=1/30` this is below `10^-4`.  This statement is
**deliberately only a second-Magnus bound**; no full time-ordered exponential
estimate is claimed.  If strain orientation ceases to be coherent, the existing
objective-strain/source collision theorem already charges that variation.
Action `31171484247` passed `50,000` objective-frame/commutator checks.

Finally, scalar triad deformation and helicity conversion together observe the
**entire three-dimensional incompressible symmetric strain**.  At the symmetric
extremizer let `Pi` be the triad plane, `k_1,k_2` the two parent directions, and

\[
D_\Pi=(S|_\Pi)^0,
\qquad
D_i=(S|_{k_i^\perp})^0.
\]

Arb interval arithmetic on the certified `r_*` bracket proves

\[
\boxed{
\|D_\Pi\|_F^2+\|D_1\|_F^2+\|D_2\|_F^2
\ge\frac{13}{20}\|S\|_F^2.
}
\]

The actual dangerous eigenmode is an off-plane shear with generalized ratio

\[
0.6602495167\ldots>\frac{13}{20}.
\]

Since the scalar shape-speed theorem gives

\[
\frac12\dot u^2+2\dot v^2
\ge\frac{43}{100}\|D_\Pi\|_F^2
\]

and each helical mixer satisfies `|zeta_i|^2=||D_i||_F^2/2`, one obtains

\[
\boxed{
\left(\frac12\dot u^2+2\dot v^2\right)
+\frac{43}{50}(|\zeta_1|^2+|\zeta_2|^2)
\ge\frac{559}{2000}\|S\|_F^2.
}
\]

Action `31171360107` Arb-certified the `13/20` tomography constant and checked
`100,000` random symmetric trace-free strains.  This is an **absolute tomography
statement**: it says the strain is visible if one records the three individual
mixers.  It must not by itself be read as a normalized transfer deficit.

The signed complex extremal tensor reveals an additional symmetry.  On every
isosceles forward triad `|x|=|y|=r`, `|x+y|=1`, the child-energy helicity tensor
factorizes exactly as

\[
\boxed{
\mathcal C_{s_1s_2s_3}
=-i C(r)s_3\varepsilon_{s_1s_2},
\qquad C(r)=\frac{\Delta(r,r,1)}{\sqrt2 r},
}
\]

with

\[
J=(\varepsilon_{ij})=
\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]

Thus the parent helicity sector is the symplectic wedge `U^T J V`.  Since

\[
M^TJM=(\det M)J,
\]

there is an exact common-parent gauge:

\[
\boxed{
(MU)^TJ(MV)=U^TJV,
\qquad M\in SL(2).
}
\]

A common determinant-one helicity deformation can therefore mix the two parent
helicities without changing the **unnormalized nonlinear parent numerator**.
This is a theorem-level countermodel to charging absolute helicity conversion as
an independent phase cost.  Capacity normalization can still register the
associated amplitude deformation; the statement here concerns the exact signed
nonlinear numerator.

If the two parents see different determinant-one maps, then exactly

\[
\boxed{
(M_1U)^TJ(M_2V)
=U^TJ(M_1^{-1}M_2)V.
}
\]

Hence the transfer-distinguishable parent variable is the **relative**
polarization deformation.  The child helicity factor remains separate.

This quotient by the common `SL(2)` gauge leads to a sharper physically aligned
strain observable.  At the symmetric extremizer define

\[
Q_{rel}
=\|D_\Pi\|_F^2
+\|D_1-D_2\|_F^2
+\|D_3\|_F^2,
\]

where `D_3` is the trace-free strain on the child polarization plane.  In
child-aligned coordinates

\[
S=\begin{pmatrix}a&b&x\\b&d&y\\x&y&-a-d\end{pmatrix},
\qquad C=\cos^2(\theta_*/2)=\frac1{4r_*^2},
\]

direct algebra gives the exact positive remainder

\[
\boxed{
Q_{rel}-\frac12\|S\|_F^2
=\frac32d^2
+(1+8C-8C^2)b^2
+(7-8C)x^2+y^2.
}
\]

Action `31172607991` Arb-certified the two nontrivial coefficients on the
certified `r_*` bracket, hence

\[
\boxed{
Q_{rel}\ge\frac12\|S\|_F^2.
}
\]

The workflow passed `156` tests and `50,000` tensor/`SL(2)`/strain checks.  The
worst tensor-factorization residual was `5.551e-16`, common-`SL(2)` invariance
residual `7.553e-14`, and the worst sampled transfer-relevant ratio was
`0.510336893`.  The earlier `13/20` theorem remains useful absolute tomography;
this `1/2` theorem uses variables that respect the exact common-parent
symplectic neutrality of the extremal nonlinear interaction.

#### Good-core helical gauge regularity: the Chern obstruction is not local

The nonzero Chern number of the helicity line bundle is a genuine global
obstruction, but the signed-good extremal core stays uniformly away from the
collinear set where the triad-normal gauge can become singular.  At
`eta_0=10^-4`, the single-edge mixed theorem gives

\[
0\le u\le\frac1{200},
\qquad |v|\le\frac1{100}.
\]

With child length normalized to one,

\[
x=r_*e^{-v-u/2},
\qquad
y=r_*e^{-v+u/2},
\]

and the parent angle satisfies the exact identity

\[
\boxed{
\cos\theta
=\frac{e^{2v}}{2r_*^2}-\cosh u.
}
\]

Monotonicity in `u,v,r_*` plus the certified `r_*` bracket gives by Arb

\[
\boxed{
\frac14<\cos\theta<\frac25,
\qquad
\sin\theta>\frac9{10}.
}
\]

If `a,b` are the unit parent directions and

\[
n=\frac{a\times b}{|a\times b|},
\]

then for every differentiable variation through the good core

\[
\boxed{
\|\dot n\|
\le\frac{10}{9}(\|\dot a\|+\|\dot b\|).
}
\]

For the triad-normal helical gauge

\[
h_s(a;n)=\frac{n\times a+i s n}{\sqrt2},
\]

one obtains the clean scale-free bound

\[
\boxed{
\|\dot h_s\|
\le\frac52(\|\dot a\|+\|\dot b\|)
}
\]

for either parent; the child obeys the same `5/2` constant on the sum of the
three carrier-direction rates.  Action `31173321354` Arb-certified the angle
margin and passed `158` tests plus `50,000` derivative checks.  The worst
normal-bound ratio was `0.962540964`, the worst helical-bound ratio
`0.427795976`, and the minimum sampled parent sine `0.930294406`.

Therefore global helicity topology does not force a local chart defect on an
efficient packet block.  Once carrier-direction variation is controlled by the
frequency-cell, affine-transport and curvature ledgers, triad-normal helical
frame variation is linearly subordinate with a universal dimensionless
constant.  Constructing that carrier-direction control in a genuine localized
PDE packet remains part of the continuum bridge.

#### Exact full time-ordered relative polarization and localized packet freezing

The common-parent `SL(2)` quotient can be used directly at the level of the
nonlinear transfer observable, avoiding a nonphysical matrix distance on the
propagators.  After factoring scalar amplification/viscous damping, let the
objective helicity spinors satisfy

\[
\dot U=-D_1U,\qquad \dot V=-D_2V,\qquad \dot Z=-D_3Z,
\]

where each `D_i` is trace free.  Since every trace-free `2x2` matrix satisfies
`D^T J+J D=0`, the parent symplectic wedge obeys the exact arbitrary-time law

\[
\boxed{
\frac d{dt}(U^TJV)=U^TJ(D_1-D_2)V.
}
\]

No commutativity, frozen strain, or Magnus truncation is needed.  A common
noncommuting time-ordered `SL(2)` history cancels pointwise.  For the full
extremal polarization factor

\[
\mathcal P=(U^TJV)(\lambda^TZ),\qquad \lambda=(1,-1)^T,
\]

one has

\[
\boxed{
\dot{\mathcal P}
=[U^TJ(D_1-D_2)V](\lambda^TZ)
-(U^TJV)(\lambda^TD_3Z).
}
\]

Consequently

\[
\boxed{
|\dot{\mathcal P}|
\le
2\sqrt{\|D_1-D_2\|_F^2+\|D_3\|_F^2}
\,\|U\|\|V\|\|Z\|.
}
\]

Action `31174812731` passed `171` tests plus `50,000` arbitrary time-ordered
histories.  The worst common-history wedge residual was `5.661e-15`.  A useful
countermodel has a common hyperbolic propagator with `||M-I||≈2.98e3` and
condition number `≈8.89e6` while the parent wedge remains invariant to numerical
precision.  Therefore Euclidean propagator distance is explicitly recorded as a
**wrong defect variable**.

The identity also survives a genuine packet residual.  If

\[
\dot U=-D_1U+F_1,\quad
\dot V=-D_2V+F_2,\quad
\dot Z=-D_3Z+F_3,
\]

the only new term is

\[
\boxed{
\mathcal R_F
=(F_1^TJV+U^TJF_2)(\lambda^TZ)
+(U^TJV)(\lambda^TF_3),
}
\]

with

\[
|\mathcal R_F|
\le\sqrt2(
\|F_1\|\|V\|\|Z\|+
\|U\|\|F_2\|\|Z\|+
\|U\|\|V\|\|F_3\|).
\]

Thus nonlinear packet forcing is an additive PDE ledger term rather than a
failure of the symplectic cancellation.

The remaining localization of the generators can also be reduced to existing
currencies.  For a Kelvin direction `a`,

\[
\dot a=-(I-aa^T)A^Ta,
\]

and direct algebra gives

\[
\|f_A(a)-f_B(b)\|
\le4L\|a-b\|+\|A-B\|_{op},
\qquad L=\max(\|A\|_{op},\|B\|_{op}).
\]

If `T=cN^-2`, `||A||<=sigma_0 N^2`, the initial directional cell has diameter
`h`, and the spatial gradient variation obeys `||A_x-A_0||<=kappa M N^2`, then

\[
\boxed{
\delta_{dir}(T)
\le e^{4c\sigma_0}(h+c\kappa M).
}
\]

Using the certified good-core angle moat and direct endpoint comparison of the
triad-normal frames gives the transfer-generator freezing bound

\[
\sqrt{
\|\Delta(D_1-D_2)\|_F^2+
\|\Delta D_3\|_F^2}
\le
\sqrt5(\Delta S_F+16\sigma N^2\delta_{dir}).
\]

On the low-strain packet branch `c sigma_0<=1/30`, Arb certifies the clean
integrated polarization remainder

\[
\boxed{
E_{pol}\le3h+\frac{15}{2}c\,\kappa M.
}
\]

Action `31174612248` passed `170` tests plus `50,000` localized packet
configurations; the worst Kelvin-Lipschitz ratio was `0.939777827`, the worst
generator-bound ratio `0.514395148`, and the simplified bound retained positive
margin.

Therefore helical localization creates **no new spatial scale** on this branch.
Adding it to the existing commutator/curvature ledger gives

\[
\boxed{
E_{loc}(M)
\le
\frac aM+\left(b+\frac{15}{2}c\right)\kappa M+3h,
}
\]

with optimum

\[
\boxed{
M_*=\sqrt{\frac{a}{(b+15c/2)\kappa}},
\qquad
E_{loc,*}
\le3h+2\sqrt{a(b+15c/2)\kappa}.
}
\]

The `h` term belongs to the already-summable smooth-symbol cell schedule, and the
helical spatial term is absorbed into the same adaptive `kappa M` moat required
by the moving-window theorem.  What remains is to estimate the actual packet
forcing `F_i` and the other localized PDE residuals, and to synchronize the
result across packet ancestry.

#### Affine-covariant residual forcing and critical-grain energy

The one-shot Young profile is affine Gaussian, so the packet forcing must be
measured in the grain's own physical metric.  Write

\[
x=X+Lz,
\qquad
\Sigma_x=LL^T,
\]

and for the resolved velocity Hessian `H_ijk=partial_j partial_k U_i(X)` define

\[
\boxed{
B_{abc}=(L^{-1})_{ai}H_{ijk}L_{jb}L_{kc},
\qquad q=L^Tk.
}
\]

`B` and `q` are invariant under a common affine change of physical coordinates
(and transform only by the orthogonal gauge in the non-unique factorization of
`Sigma_x`).  For the quadratic Taylor remainder, Wick's formula gives the exact
unprojected Gaussian norm

\[
\frac{\|R_2\cdot\nabla\psi\|_2^2}{\|\psi\|_2^2}
=
\frac14\left(2\|q\cdot B\|_F^2+\operatorname{tr}(q\cdot B)^2\right)
+
\frac1{16}\left(6\|T\|_F^2+9\|\operatorname{tr}T\|^2\right),
\]

where `T=Sym B`.  The quadratic `q dot B` term is a chirp/covariance tangent.
The cubic trace is a center/carrier tangent.  After projecting to the complement
of the full center--carrier--covariance--chirp Gaussian tangent space, the first
true scalar shape forcing is exactly third Hermite chaos:

\[
\boxed{
\frac{\|F_\perp\|_2^2}{\|\psi\|_2^2}
=
\frac38\|\operatorname{Sym}B\|_F^2,
\qquad
\frac{\|F_\perp\|_2}{\|\psi\|_2}
\le\frac{\sqrt6}{4}\|B\|_F.
}
\]

Action `31179739773` passed `174` tests and `50,000` affine/Hermite checks.  An
explicit transformed grain with condition number `1.091e10` retained affine
invariance to residual `5.608e-16`.  Thus Euclidean aspect is not itself a
forcing cost.

The kernel `Sym B=0` is also physical rather than missing coercivity.  Under
differentiated incompressibility it has the exact five-dimensional form

\[
\boxed{
B_{abc}=\varepsilon_{abd}M_{dc}+\varepsilon_{acd}M_{db},
\qquad M=M^T,\quad \operatorname{tr}M=0,
}
\]

hence

\[
\boxed{V(z)=z\times(Mz).}
\]

It is tangent to normalized Gaussian level spheres (`V dot z=0`) and a carrier
sees only a quadratic chirp.  Action `31180257124` passed `180` tests plus
`50,000` swirl checks with representation residual `1.068e-14`.  This mode must
be routed to vector/helical spatial variation rather than charged as scalar
envelope forcing.

That vector channel is explicit.  Define the physical strain variation per
grain coordinate

\[
C_{ijc}=\frac12(H_{ijk}+H_{jik})L_{kc}.
\]

Each `C_c` is symmetric trace free.  Therefore the existing Arb-certified
transfer-relevant tomography applies to each slice.  For `z` standard Gaussian,

\[
\boxed{
\mathbb E_z Q_{rel}(C_cz_c)
=
\sum_cQ_{rel}(C_c)
\ge\frac12\|C\|_F^2.
}
\]

Action `31181076691` passed `187` tests; its `50,000` Hessian/grain checks had
worst RMS ratio `0.506683513`.  Tested swirl kernels had vanishing scalar
third-Hermite signal to `8.639e-16` but positive sampled polarization signal.
No aspect-independent lower bound comparing this physical `C` channel directly
to affine-normalized `||B||` is claimed.

The quadratic phase/chirp tangent has an exact triad explanation.  If

\[
(\partial_t+U\cdot\nabla)\phi_i=\rho_i,
\]

then

\[
\boxed{
(\partial_t+U\cdot\nabla)(\phi_1+\phi_2-\phi_3)
=\rho_1+\rho_2-\rho_3.
}
\]

Carrier resonance and Hessian/chirp lock are consequently preserved by an
arbitrary **common nonlinear** resolved flow; at exact carrier resonance the
common velocity-Hessian source cancels from the chirp equation.  For
role-dependent velocities only the differences `U-U_i` remain.  Action
`31180506627` passed `182` tests plus `50,000` phase-lock checks, with common
Hessian-source residual exactly zero in the tested algebra.  Common non-affine
phase is therefore a gauge, not `F_i`.

Bulk viscosity is tangent as well.  For a complex Gaussian

\[
\psi(y)=C\exp[-\tfrac12y^TGy+i k\cdot y],
\]

\[
\boxed{
\Delta\psi/\psi
=y^TG^2y-2i k^TGy-|k|^2-\operatorname{tr}G,
}
\]

a polynomial of degree at most two.  Thus bulk viscosity changes Gaussian
parameters but creates no third-Hermite transverse forcing.  Viscous window
boundary flux remains an unresolved localized PDE term.

The frozen inverse-Young profile also carries a genuinely affine physical mass.
For `||f-F||_(3/2)<=1/100` and the certified shell, the `|F|^(3/2)` covariance
`Gamma_xi` and physical `L^2` covariance obey the exact uncertainty matrix

\[
\boxed{\Sigma_x\Gamma_\xi=I/3.}
\]

Action `31179827015` Arb-certified

\[
\boxed{
\lambda_{\min}(\Sigma_x)^{1/2}>\frac{2}{3N}
}
\]

and, on the radius-two covariance ellipsoid `E_2`, the actual-role mass

\[
\boxed{
\int_{E_2}|u|^2
\ge\frac3{10}(\det\Sigma_x)^{1/6}.
}
\]

Define the geometric-mean physical radius and affine critical mass

\[
\boxed{
r_g=(\det\Sigma_x)^{1/6},
\qquad
\mathsf M_{aff}(E)=r_g^{-1}\int_E|u|^2.
}
\]

Then

\[
\boxed{\mathsf M_{aff}(E_2)\ge3/10.}
\]

This quantity is exactly invariant under the isotropic Navier--Stokes scaling.
For a fresh family with `M_aff>=eta` and overlap multiplicity `P`, physical
energy gives the exact ancestry budget

\[
\boxed{
\sum_jr_{g,j}\le\frac{P\|u(t)\|_2^2}{\eta}.
}
\]

Moreover if `s=Nr_g` and `A=N lambda_max(Sigma_x)^(1/2)`, the axis lower bound
gives

\[
\boxed{A\le\frac94s^3.}
\]

The same radius is dynamically natural.  Since `Sigma_x=P/2`, the exact affine
Gaussian covariance equation is

\[
\dot\Sigma_x=A\Sigma_x+\Sigma_xA^T+\nu I.
\]

For incompressible `tr A=0`,

\[
\boxed{
\frac d{dt}\log r_g
=\frac\nu6\operatorname{tr}\Sigma_x^{-1}\ge0,
\qquad
\frac d{dt}r_g^2\ge\nu.
}
\]

Thus inviscid affine strain changes aspect but preserves geometric radius
exactly; viscosity increases it.  This prevents shear from manufacturing or
erasing the fresh-radius currency inside the affine packet model.

Thus natural geometric scale implies bounded aspect; an extremely elongated
**fresh** grain consumes a larger physical radius budget.  Crucially, arbitrary
common affine anisotropy remains an exact Young symmetry (tested to aspect
`1e8`), so aspect alone is not a Bellman/replication deficit.  Reused elongated
grains must be handled by spacetime ancestry and the dynamic curvature/
polarization ledgers.  Action `31180882083` passed `186` tests and `50,000`
affine-grain budget checks; the downstream rational provenance was cleaned in
commit `21d8976` so the `2/3` and `3/10` inputs are taken directly from the Arb
shell certificate rather than reconstructed from floating approximations.

#### Ellipsoidal moving windows: affine curvature balance

The affine grain can now carry its own moving physical window.  Let

\[
\chi_{L,M}(x,t)
=\chi_0\!\left(\frac{L(t)^{-1}(x-X(t))}{M}\right),
\]

and transport

\[
\dot X=U(X),
\qquad
\dot L=A(X)L,
\qquad A=\nabla U.
\]

With `z=L^{-1}(x-X)`, direct differentiation gives the exact identity

\[
\boxed{
(\partial_t+U\cdot\nabla)z
=L^{-1}\big(U(X+Lz)-U(X)-A(X)Lz\big).
}
\]

Thus translation and the full affine velocity jet cancel from the material
window derivative.  If

\[
\kappa_{aff}
:=\sup\|L^{-1}(\nabla^2U)[L,L]\|
\]

on the transition region, `|grad chi_0|<=C_chi`, and the base transition is in
`|z/M|<=R_chi`, Taylor's integral remainder gives

\[
\boxed{
|(\partial_t+U\cdot\nabla)\chi_{L,M}|
\le
\frac{C_\chi R_\chi^2}{2}\,\kappa_{aff}M.
}
\]

The opposite localization cost remains `1/M`.  From the Arb-certified shell
axis lower bound `l_min>2/(3N)`,

\[
\boxed{
N^{-1}\|\nabla_x\chi_{L,M}\|_\infty
\le\frac{3C_\chi}{2M}.
}
\]

For a physical coarse-graining kernel `G_N=N^3G(N\cdot)`, there is also an exact
L2 commutator theorem.  Writing

\[
[\chi,G_N*]f(x)
=\int G_N(y)(\chi(x)-\chi(x-y))f(x-y)\,dy
\]

and using the mean-value theorem plus Young gives

\[
\boxed{
\|[\chi,G_N*]f\|_2
\le
\frac{m_1(G)}{N}\|\nabla\chi\|_\infty\|f\|_2,
\qquad
m_1(G)=\int |y||G(y)|\,dy.
}
\]

Therefore on the certified affine grain

\[
\boxed{
\|[\chi_{L,M},G_N*]f\|_2
\le
\frac{3m_1(G)C_\chi}{2M}\|f\|_2.
}
\]

The geometric window ledger consequently has the same exact balance form

\[
\boxed{
E_{aff}(M)
\le
\frac aM+b\kappa_{aff}M,
\qquad
M_*=\sqrt{\frac{a}{b\kappa_{aff}}},
\qquad
E_{aff,*}=2\sqrt{ab\kappa_{aff}}.
}
\]

Action `31182421821` (before the final commutator unit-test addition) passed
`194` tests and `50,000` ellipsoidal-window checks.  It saw worst Taylor
leakage/bound ratio `0.992529785`, shell-gradient ratio `0.999859538`, affine
curvature invariance residual `1.933e-08`, and an extreme transformed condition
number `1.444e9`.  The final Actions rerun records the commutator addition as
well.

This closes the **geometry and generic filter-commutator scaling** of an affine
moving moat.  It does not yet provide the complete localized SGS/pressure
packet identity: pressure boundary work, `RU`, viscous boundary flux and
partition-overlap coefficients must still be inserted for the actual PDE packet
construction.

#### Actual smooth-SGS affine packet equation: two levels and Hermite forcing

The packet bridge must distinguish two descriptions of the same Navier--Stokes
nonlinearity.  The near-extremal Young/helical roles are roles of the **full
velocity** `u`, while `R`, `Pi` and `RU` belong to the coarse-grained resolved
energy equation.  Mixing these two levels would double count SGS transport.

On the signed-good core the smaller parent, after normalizing the child to
frequency `N`, satisfies

\[
|k_{parent}|/N\ge r_*e^{-1/80}>3/5.
\]

Action `31186597005` Arb-certified the clean `3/5` bound.  Choose a strict
low-frequency transporter

\[
V=S_{N/4}u,
\qquad \operatorname{supp}\widehat V\subset\{|\xi|\le N/4\}.
\]

Then

\[
\operatorname{supp}\widehat{V\otimes V}
\subset\{|\xi|\le N/2\},
\]

so the low--low product cannot feed a selected near-extremal role.  The role
nonlinearity therefore splits exactly into **low--high transport plus high--high
transfer**.  The Kelvin background is thus the actual low-frequency part of the
same Navier--Stokes velocity, not an externally imposed model flow.

For the affine jet `V=A(x-X)`, `tr A=0`, the divergence-free low--high
linearization has the exact Fourier equation

\[
\boxed{
(\partial_t-(A^T\xi)\cdot\nabla_\xi)\widehat w
=\left[-A+2\widehat\xi\widehat\xi^TA\right]\widehat w
-\nu|\xi|^2\widehat w.
}
\]

The pressure is determined by transversality,

\[
i\widehat\pi=-2\frac{\xi\cdot A\widehat w}{|\xi|^2},
\]

and is exactly the `2 khat khat^T A` correction rather than an independent
packet source.  Along characteristics this gives

\[
\boxed{
\dot k=-A^Tk,
\qquad
\dot a=(-A+2\widehat k\widehat k^TA)a-\nu|k|^2a.
}
\]

In the objective transverse frame,

\[
\boxed{
\dot c=-\operatorname{sym}(E^TAE)c-\nu|k|^2c,
}
\]

so the exact PDE low--high linearization, Kelvin carrier law and the previously
used two-component polarization generator are one transport identity.  Action
`31187637446` passed `211` tests plus `50,000` checks; the pressure/Kelvin and
objective-coordinate residuals were `8.702e-15` and `2.920e-14` respectively.

The spectral packet cell has an equally exact transport law.  For affine
transport, the correct Heisenberg residual of a scalar multiplier is

\[
\boxed{
\partial_tm-(A^T\xi)\cdot\nabla_\xi m.
}
\]

A cell transported by the dual affine flow makes this zero.  Action
`31186597005` saw zero algebraic residual in `50,000` random affine tests.  Thus
common carrier/covariance motion is not spectral leakage.

A spatially localized role must also remain divergence free.  Rather than
multiplying a helical mode by a window and silently dropping pressure, use the
weak test packet

\[
\boxed{
\Psi_i=M_{i,N}(D)(\chi_{L,M}\phi_i),
}
\]

where `M_(i,N)` is the smooth shell-localized Leray/helical/cell multiplier.
Then `div Psi_i=0` exactly and full Navier--Stokes gives

\[
\boxed{
\frac d{dt}\langle u,\Psi_i\rangle
=\langle u,\partial_t\Psi_i\rangle
+\langle u\otimes u,\nabla\Psi_i\rangle
-\nu\langle\nabla u,\nabla\Psi_i\rangle.
}
\]

Pressure vanishes exactly from this microscopic weak equation.  Since the
complete shell multiplier is smooth away from `xi=0`, its kernel has finite
first moment and

\[
\boxed{
\|[\chi,M_{i,N}]f\|_2
\le \frac{m_1(K)}N\|\nabla\chi\|_\infty\|f\|_2
\le \frac{3m_1(K)C_\chi}{2M}\|f\|_2.
}
\]

Action `31188118452` passed `215` tests plus `50,000` checks; the Leray pressure
residual was `1.791e-15` and the discrete commutator/bound ratio was
`0.996773131`.  Thus enforcing divergence-free/helical localization only
renormalizes the existing `1/M` moat coefficient.

The remaining field residual must **not** be inserted wholesale as the spinor
forcing `F_i`.  In the affine Gaussian/Hermite frame it decomposes as

\[
\boxed{
f_i=f_i^{(0)}+f_i^{H_1,pol}+f_i^{tan,\le2}
+f_i^{H_3,env}+f_i^{higher}.
}
\]

Only the degree-zero projection `f_i^(0)` enters the forced symplectic spinor
identity.  The exact Gaussian-chaos identities are

\[
\boxed{
\mathbb E[(T:H_3)P_2(Z)]=0\quad(\deg P_2\le2),
\qquad
\mathbb E[(T:H_3)^2]=6\|T\|_F^2,
}
\]

and, for the spatial polarization sideband,

\[
\boxed{
\mathbb E\left[\sum_c C_cZ_c\right]=0,
\qquad
\mathbb E\left|\sum_cC_cZ_c\right|^2=\sum_c\|C_c\|_F^2.
}
\]

Hence the `Sym B` envelope residual is an orthogonal `H_3` daughter mode, and
spatial polarization curvature is an `H_1` sideband; both have zero direct
projection into the base Gaussian spinor.  Scalar degree `<=2` contains the
center/carrier/covariance/chirp and bulk-viscosity tangent motion.  Action
`31187104071` passed `206` tests plus `20,000` Gaussian-chaos checks, with
`H_3`/degree-`<=2` projection `1.889e-14` and `H_1`/base projection
`4.521e-16`.

The degree-zero forcing produced by the **non-Gaussian parent remainders** is
already the transfer cross-error.  If actual unit roles `f,g,h` and Gaussian
roles `F,G,H` are each within `epsilon`, then

\[
\boxed{
|T(f,g,H)-T(F,G,H)|
\le 2\epsilon+3\epsilon^2+\epsilon^3,
}
\]

while

\[
|T(f,g,h)-T(f,g,H)|\le\epsilon.
\]

Their sum is exactly the old replacement polynomial

\[
3\epsilon+3\epsilon^2+\epsilon^3.
\]

At one percent it is `0.020301+0.010000=0.030301`.  Action `31187443147`
passed `208` tests plus `50,000` complex trilinear checks.  Thus, at the
**work level**, degree-zero parent/cross forcing belongs to the existing
`eta_j` / `Xi` ledger rather than a new `L^2` forcing currency.  No bound on the
entire nonlinear residual norm is inferred from this statement.

The coarse-grained resolved SGS energy remains a separate macroscopic ledger.
Writing the affine-window transporter as `V_aff` and `W=U-V_aff`, the
pressure-free combined-work leakage is

\[
\widetilde L_\chi
=\int\!\!\int\nabla\chi\cdot(eW+RU-\nu\nabla e).
\]

The nonviscous boundary terms satisfy pointwise

\[
\frac12|U|^2|W|+|R|_F|U|
\le\frac23|U|^3+\frac16|W|^3+\frac23|R|_F^{3/2},
\]

so a leakage `L_cub` forces the scale-critical annular charge

\[
\boxed{
\int_A(|U|^3+|W|^3+|R|_F^{3/2})
\ge\frac{|L_{cub}|}{\|\nabla\chi\|_\infty}
\ge\frac{2M}{3NC_\chi}|L_{cub}|.
}
\]

The resolved viscous boundary term on a lifetime `T=cN^-2` obeys, for the
certified smooth support,

\[
\boxed{
|L_\nu|
\le\frac{33}{20}\frac{c\nu C_\chi}{M}
\sup_t\|U(t)\|_2^2,
}
\]

because Arb certifies `exp(1/20)<11/10`.  Thus viscosity at the resolved
boundary adds to the same `1/M` coefficient; this is distinct from bulk
Gaussian viscosity, which is tangent microscopically.

The SGS stress itself has the exact increment identity

\[
\boxed{
R(x)=\int G_\ell(r)\,\delta_ru\otimes\delta_ru\,dr
-\left(\int G_\ell(r)\delta_ru\,dr\right)^{\otimes2}.
}
\]

For `g_1=||G||_1`,

\[
\boxed{
|R(x)|_F^{3/2}
\le(1+g_1)^{3/2}g_1^{1/2}
\int|G_\ell(r)|\,|\delta_ru(x)|^3dr.
}
\]

Thus a large `RU` branch is an Onsager-type cubic **velocity-increment charge at
the actual filter scale**, not an independent stress-tensor mystery.  Action
`31187768828` passed `212` tests plus `50,000` checks; the stress/increment-cubic
ratio was `0.638619634` and the exact increment-identity residual
`1.349e-13`.

Finally, for a complete quadratic spatial partition `sum eta_alpha=1`,

\[
\boxed{
\sum_\alpha(e\,\partial_t\eta_\alpha+\nabla\eta_\alpha\cdot F)=0.
}
\]

There is no global packet-count boundary loss from the partition itself.
Interfaces become real cross-errors only after selecting a lineage/subfamily.
The first affine-SGS boundary Action `31186600162` failed solely because an
exact rational equality `3/2*11/10=33/20` was tested using Arb ball equality;
this is recorded as certificate-semantics provenance, not a mathematical
countermodel.  The corrected run `31186719319` and the increment-upgraded run
`31187768828` are green.

This gives a concrete no-double-counting rule:

\[
\boxed{
\begin{array}{c|c}
\text{microscopic full-velocity role}&
\text{macroscopic resolved SGS energy}\\ \hline
\text{Kelvin low--high transport}&\Pi\text{ / combined work}\\
\text{base-spinor cross work}\to\Xi&RU\text{ boundary transport}\\
H_1,H_3\text{ sidebands}&\text{pressure annular charge}\\
\text{Gaussian tangents/gauges}&\text{resolved viscous/window leakage}
\end{array}
}
\]

Pressure, `RU`, viscosity and curvature must never be charged in both columns.

#### Affine curvature connection: sideband dephasing sources

The `H_1/H_3` sideband geometry itself has a common affine transport gauge.  For
any smooth transporter `V`, along

\[
\dot X=V(X),\qquad \dot L=A L,
\qquad A=\nabla V(X),
\]

let `H=nabla^2 V(X)`, `F=D_t^V V`, and

\[
B=L^{-1}H[L,L],\qquad A_{aff}=L^{-1}AL.
\]

The material Hessian commutator gives the exact covariant identity

\[
\boxed{
\dot B+2A_{aff}B=L^{-1}(\nabla^2F)[L,L].
}
\]

Thus raw `dot B` is not itself a defect: common affine deformation acts as a
connection on curvature.  For the resolved Navier--Stokes transporter

\[
D_t^V V=-\nabla P-\nabla\cdot R+\nu\Delta V,
\]

so

\[
\boxed{
\dot B+2A_{aff}B
=L^{-1}\left[-\nabla^3P-\nabla^2\nabla\cdot R
+\nu\nabla^2\Delta V\right][L,L].
}
\]

This identifies the only source channels for curvature-sideband dephasing after
the affine connection is removed: pressure third derivatives, differentiated
SGS stress, or viscous fourth velocity derivatives.  It is a source theorem,
not yet a positive sideband-coherence cost.

At this derivative level pressure is even more local.  The 3D pressure kernel
has homogeneity `-3`; three derivatives have homogeneity `-6`.  Dyadic
three-dimensional packet packing therefore leaves

\[
\boxed{6-3=3}
\]

summable powers in the far field.  Action `31189036483` passed `217` tests plus
`50,000` curvature-connection checks; the worst normalized connection residual
was `1.361e-13`, while the resolved-source split and `6-3=3` homogeneity
residuals were zero in the stress suite.

#### Covariant Hermite sideband no-escape and daughter ancestry

The orthogonal `H_1/H_3` channels can now be controlled without introducing a
raw residual norm or assuming an almost-constant complex phase.  The first
ingredient is a Banach-space interval-variation theorem.  For every absolutely
continuous forcing `f:[0,T]->X`,

\[
\boxed{
\int_0^T\|f\|dt
\le
\left\|\int_0^T fdt\right\|
+\frac T2\int_0^T\|\dot f\|dt.
}
\]

Thus, if `A=int ||f||`, either the first Duhamel impulse has norm at least
`A/2`, or the total variation is at least `A/T`.  There is no phase-lock
hypothesis in this statement.

For the affine curvature connection

\[
\dot B+2A_{aff}B=S,
\qquad
S=L^{-1}\nabla^2(D_t^V V)[L,L],
\]

pull back the output index by `Pdot=2 A_aff P`.  Then

\[
\boxed{\widetilde B=PB,\qquad \dot{\widetilde B}=PS.}
\]

Common affine curvature transport has therefore disappeared before coherence is
measured.  The `H_3` envelope forcing in this interaction frame satisfies

\[
\boxed{
\frac{\|f_3\|_2^2}{\|g\|_2^2}
=\frac38\|\operatorname{Sym}\widetilde B\|_F^2,
}
\]

and its dephasing derivative is sourced by `Sym(P S)`, hence ultimately by the
already identified acceleration-Hessian channels.  Action `31192448279` passed
`228` tests; the exact pullback and pulled-H3-source residuals were
`3.252e-15` and `3.860e-15`.

A coherent daughter is also visible in the **critical Young norm**.  Let
`g(z)=exp(-|z|^2/4)`.  The `L^2` Gaussian is standard normal whereas the
`|g|^(3/2)` measure has variance `4/3`.  Standard Gaussian
hypercontractivity, Paley--Zygmund, and interval-certified arithmetic give the
clean implications

\[
\boxed{
\frac{\|H_3\text{-daughter}\|_{3/2}}{\|g\|_{3/2}}
>\frac1{160}
\frac{\|H_3\text{-daughter}\|_2}{\|g\|_2},
}
\]

and, for a two-component `H_1` polarization sideband,

\[
\boxed{
\frac{\|(Cz)g\|_{3/2}}{\|g\|_{3/2}}>\frac1{16}\|C\|_F.
}
\]

The workflow certifies only the arithmetic constants; Gaussian
hypercontractivity is the analytic theorem used as input.

There is an exact **odd-Hermite triad selection rule**.  For a centered resonant
affine Gaussian triad, simultaneous inversion of all centered frequency
deviations preserves the Gaussian trilinear weight.  Hence

\[
\boxed{
T(P_{n_1}G_1,P_{n_2}G_2,P_{n_3}G_3)=0
\quad\text{if }n_1+n_2+n_3\text{ is odd}.
}
\]

In particular

\[
\boxed{T(H_1,G,G)=T(H_3,G,G)=0.}
\]

A single odd daughter cannot rescue the base Gaussian transfer at first order.
If the three roles have base critical norms `b_i` and odd-sideband norms
`rho_i`, every rescue term contains at least two odd factors and sharp Young
gives

\[
\boxed{
|T_{rescue}|
\le A_3(\rho_1\rho_2b_3+\rho_1\rho_3b_2
+\rho_2\rho_3b_1+\rho_1\rho_2\rho_3).
}
\]

For one odd degree-`<=3` perturbation `(1+P)G`, put

\[
\sigma^2=\mathbb E_{|G|^{3/2}}|P|^2.
\]

Degree-three hypercontractivity gives `E|P|^4<=729 sigma^4`.  If
`sigma<=1/80`, parity and uniform convexity imply

\[
\boxed{
\frac{\|(1+P)G\|_{3/2}^{3/2}}{\|G\|_{3/2}^{3/2}}
\ge1+\frac3{16}\sigma^2,
\qquad
\operatorname{Def}_{transfer}\ge\frac1{16}\sigma^2.
}
\]

Action `31192451039` passed `227` tests (`1` optional Arb test skipped in the
non-Arb workflow) plus `50,000` parity/rescue checks; the worst `H_1` and `H_3`
parity residuals were `1.019e-15` and `4.527e-14`.

The affine-normalized differentiated-incompressible curvature itself has no
unclassified static sector.  For

\[
B_{abc}=B_{acb},\qquad B_{aac}=0,
\]

let `T=Sym B`, `t_c=T_aac`, and define

\[
\boxed{
B^E_{abc}=T_{abc}
-\frac12(\delta_{ab}t_c+\delta_{ac}t_b)
+\delta_{bc}t_a.
}
\]

Then `B=B^E+B^H` orthogonally, `Sym B^H=0`, and `B^H` remains divergence
free.  The hook sector is exactly the five-dimensional swirl representation

\[
\boxed{
B^H_{abc}=\varepsilon_{abd}M_{dc}+\varepsilon_{acd}M_{db},
\qquad M=M^T,\quad\operatorname{tr}M=0.
}
\]

The exact norm identities are

\[
\|B^E\|^2=\|T\|^2+3\|t\|^2,
\qquad
\|t\|^2\le\frac53\|T\|^2,
\]

\[
\|B^H\|^2=6\|M\|^2,
\qquad
\|C^H\|^2=\frac14\|B^H\|^2,
\]

where `C^H` is the intrinsic grain-coordinate `H_1` strain sideband.  Hence

\[
\boxed{
\|\operatorname{Sym}B\|_F^2+\|C^H\|_F^2
\ge\frac16\|B\|_F^2.
}
\]

Representation-theoretically, the fifteen-dimensional curvature space splits
as

\[
\boxed{15=(7\oplus3)_{envelope}\oplus5_{swirl}.}
\]

Action `31192445017` passed `227` tests (`1` optional Arb skip) plus `50,000`
checks; all tensor identities were at `1e-14` or better.  This `1/6`
observability is intrinsic to the affine grain.  It does **not** remove the
existing aspect-sensitive caveat when comparing the hook norm to physical
Euclidean relative-helicity curvature.

These ingredients close a local `H_3` no-escape theorem.  Define

\[
I_3=\int_0^T\|\operatorname{Sym}\widetilde B(t)\|_Fdt,
\qquad
J_3=\int_0^T\|\operatorname{Sym}\widetilde S(t)\|_Fdt.
\]

Then either

\[
\boxed{J_3\ge I_3/T,}
\]

or the first-Duhamel `H_3` impulse has

\[
\boxed{\delta_1^2\ge\frac3{32}I_3^2.}
\]

If nonlinear sideband feedback is at least `delta_1/2`, this is already a
high--high sideband-interaction branch.  Otherwise the surviving projected
`H_3` daughter obeys

\[
\boxed{\delta^2\ge\frac3{128}I_3^2.}
\]

The change from the `L^2` Gaussian to the critical Young Gaussian only increases
the `H_3` second moment.  For variance `s=4/3`,

\[
\boxed{
\mathbb E_{\gamma_s}|T:H_3|^2
=6s^3\|T\|^2+9s(s-1)^2\|t\|^2
\ge6\|T\|^2.
}
\]

Thus its critical second-moment parameter satisfies `sigma>=delta`.  If
`sigma>=1/80`, there is a definite daughter-capacity event.  Otherwise the
single-role transfer cost is at least `sigma^2/16`; splitting it against
possible pair-sideband rescue gives the clean constant

\[
\boxed{c_{H3}=\frac3{4096}.}
\]

Every `H_3` curvature impulse therefore has one of the five exits

\[
\boxed{
\begin{array}{ll}
J_3\ge I_3/T,&\text{acceleration-Hessian dephasing source},\\
\|R_{Duhamel}\|\ge\delta_1/2,&\text{nonlinear sideband feedback},\\
\sigma\ge1/80,&\text{definite daughter capacity},\\
\operatorname{Def}_{transfer}\ge(3/4096)I_3^2,&\text{base-edge cost},\\
R_{pair}\ge(3/4096)I_3^2,&\text{pair-sideband rescue}.
\end{array}}
\]

Action `31193125969` passed `229` tests (`1` optional Arb skip) and `50,000`
five-branch checks; the minimum transfer-cost and pair-cost margins were
`3.321e-11` and `4.293e-11`.  The preceding run `31192748440` failed only
because its synthetic small-sideband sampler called `uniform(low,high)` in the
narrow case `low>0.999/80`; it is implementation provenance, not a theorem
countermodel.

Finally, pair rescue itself now enters the existing ancestry machinery without
assuming spatial disjointness of Hermite modes.  Odd parity implies every
nonzero rescue triad contains exactly **two odd daughter endpoints** and one
even/base role, so rescue projects canonically to an ordinary weighted graph.
Split its mass as

\[
W=W_{cross}+W_{same}
\]

according to the ancestry labels of the two odd endpoints.  A large cross
fraction is an existing `Xi` interaction.  On same-ancestry edges define

\[
w_v=\frac{d_v}{2W_{same}}.
\]

For an ancestry class `A`, `W_A=sum_(v in A)w_v` equals exactly the fraction of
same-edge rescue mass lying in `A`.  Hence

\[
\boxed{
Q_{anc}-Q_{at}
=\sum_A\sum_{u\ne v\in A}w_uw_v.
}
\]

For any `h>0`, `0<alpha<1`, the same-edge graph has the exact routing:

- `H_at<h` -> one odd daughter has endpoint weight `>exp(-h)`, a dominant-reuse
  branch;
- `H_at>=h` and `H_anc>=alpha h` -> component/Bellman collision entropy;
- `H_at>=h` and `H_anc<alpha h` -> distinct same-ancestry endpoint-pair mass
  at least
  \[
  \boxed{e^{-\alpha h}-e^{-h}},
  \]
  which enters the existing ancestry-cycle contraction theorem.

Together with the cross branch, quadratic sideband rescue is therefore routed
into `Xi`, dominant reuse, Bellman entropy, or ancestry cycles.  Action
`31193531201` passed `231` tests (`1` optional Arb skip) plus `50,000` weighted
graph checks; the class-mass and collision-chain identities were accurate to
`5.551e-16` and `9.992e-16`.  The preceding run `31193443926` failed because a
test fixture contained only one cross edge out of three while incorrectly
expecting a majority-cross branch; the theorem routed it correctly to Bellman.

#### Physical `H_1` / swirl bridge and physical mild-aspect no-escape

The five-dimensional hook sector is connected to the **transfer-facing
relative-parent/child polarization variables** on a controlled mild-aspect
branch, with no scalar-shape `D_Pi` term.

For hook curvature `B^H`, define the physical strain-gradient slices

\[
\boxed{C_c^{phys}=\operatorname{Sym}(L B_c^H L^{-1})}
\]

and

\[
\boxed{
Q_{pol}(C)=\sum_c\left(
\|D_1(C_c)-D_2(C_c)\|_F^2+
\|D_3(C_c)\|_F^2
\right).
}
\]

At `L=I`, direct expansion in the five symmetric trace-free swirl variables
makes `Q_pol-(1/10)||B^H||^2` block diagonal; Arb certifies positivity throughout
the `r_*` bracket:

\[
\boxed{Q_{pol}(I)\ge\frac1{10}\|B^H\|_F^2.}
\]

After polar/global-scalar normalization, if `kappa=cond(L)`,

\[
\|C^{phys}(L)-C^{phys}(I)\|_F
\le(\kappa-1)\|B^H\|_F,
\]

while the observable has operator norm at most `sqrt(5)`.  Hence

\[
\sqrt{Q_{pol}(L)}
\ge
\left(\frac1{\sqrt{10}}-\sqrt5(\kappa-1)\right)\|B^H\|_F.
\]

On

\[
\boxed{\operatorname{cond}L\le21/20,}
\]

Arb certifies

\[
\boxed{Q_{pol}(L)\ge\frac1{25}\|B^H\|_F^2.}
\]

For a real symmetric trace-free `2x2` generator,

\[
D^2=\frac12\|D\|_F^2I.
\]

Thus the **auxiliary relative-parent/child coordinate** has forcing energy
`Q_pol/2`, hence at least `||B^H||^2/50`.  That auxiliary number is useful for
polarization algebra but is **not yet one physical Young role**.

For the three physical packet roles separately,

\[
\begin{aligned}
\sum_{i=1}^3\|F_i^{H1}\|^2
&=\frac12(\|D_1\|_F^2+\|D_2\|_F^2+\|D_3\|_F^2)\\
&\ge\frac14(\|D_1-D_2\|_F^2+\|D_3\|_F^2)
=\frac14Q_{pol},
\end{aligned}
\]

so the physical role-level bound used by the odd-Hermite transfer theorem is

\[
\boxed{
\sum_i\|F_i^{H1}\|^2
\ge\frac1{100}\|B^H\|_F^2.
}
\]

Action `31196456570` on commit `4a5729c` passed `237` tests plus `50,000`
configurations and was Arb-certified.  It retained the previous worst isotropic
and mild-aspect ratios `0.121618531` and `0.104335578`, while explicitly
recording that the relative coordinate is non-unitary and not itself a Young
role.  The earlier run `31194787885` remains valid provenance for the pointwise
`Q_pol` theorem but predates this role-level distinction.

The spacetime theorem must also respect the fact that homogeneous polarization
transport is `SL(2)`, not unitary.  On the existing low-strain lifetime branch,
each physical role has base action budget

\[
\boxed{
K:=\int_0^T\|G_i(t)\|_{op}dt\le1/30.
}
\]

Therefore every base propagator and inverse has singular values in
`[e^{-K},e^K]`.  Let

\[
I_1=\int_0^T\|B^H(t)\|_Fdt.
\]

Pulling the physical forcing into the interaction picture costs one factor
`e^{-K}`; pushing the Duhamel impulse back to physical roles costs a second.
Arb certifies

\[
\frac{e^{-1/30}}{10}>\frac1{11},
\qquad
 e^{-2/15}>\frac56.
\]

Thus the Banach variation theorem gives either

\[
\boxed{J_1\ge\frac{I_1}{11T}}
\]

for the interaction-picture covariant forcing variation, or a **physical
three-role** first-Duhamel daughter with

\[
\boxed{\delta_1^2\ge\frac{I_1^2}{480}.}
\]

If physical nonlinear sideband feedback cancels at least half of `delta_1`, it
is already a high--high sideband-interaction branch.  Otherwise surviving total
role-daughter energy is at least

\[
\delta^2\ge I_1^2/1920.
\]

At least one of the three physical roles therefore has daughter energy
`>=I_1^2/5760`.  Below critical sideband size `1/80`, the one-role odd-Hermite
uniform-convexity theorem gives deficit at least `I_1^2/92160`; splitting this
against possible pair rescue yields the final physical local alternative

\[
\boxed{
\operatorname{Def}_{transfer}\ge\frac{I_1^2}{184320}
\quad\text{or}\quad
R_{pair}\ge\frac{I_1^2}{184320}.
}
\]

Hence

\[
\boxed{c_{H1}^{phys}=1/184320.}
\]

The irreducible curvature split still gives

\[
I_B\le\sqrt6 I_3+I_1.
\]

Therefore either `I_3>=I_B/(2 sqrt(6))` or `I_1>=I_B/2`.  The `H_3` branch
remains stronger; the `H_1` branch gives the common physical mild-aspect
full-curvature constant

\[
\boxed{c_{mild,curv}^{phys}=1/737280.}
\]

Action `31196458909` on `4a5729c` passed `237` tests plus `50,000` physical
five-branch checks and was Arb-certified for the nonunitary conditioning.  Its
minimum pair-cost, transfer-deficit and full-channel margins were `1.546e-11`,
`5.585e-12`, and `5.238e-09`.

This **supersedes** the idealized constants `1/25600` and `1/102400` printed by
run `31195130386`, which treated the relative-coordinate pullback as if it were
an isometry and did not split the final daughter across the three physical Young
roles.  That run is retained as algebra/provenance, not as the preferred
physical H1 no-escape theorem.  The pointwise `Q_pol>=1/25||B^H||^2` theorem is
unchanged.

The earlier open `H_1` source branch is now expanded quantitatively.  Let
`P_H` be the fixed hook projector, `B^H=P_HB`, and

\[
G_c=L B_c^H L^{-1}.
\]

The exact projector/product rules are

\[
\boxed{
\dot B^H=P_H(-2A_{aff}B+S),
\qquad
\dot G_c=A G_c+L\dot B_c^H L^{-1}-G_cA.
}
\]

No commutation of `P_H` with affine transport is assumed.  On the original
`cond(L)<=21/20`, `K<=1/30` branch, the three-role interaction forcing
variation obeys

\[
\boxed{
J_1
\le
2\int_0^T\|S_{src}\|dt
+54\int_0^T\|A\|\,\|B\|dt.
}
\]

Hence `J_1>=I_1/(11T)` implies

\[
\boxed{
\int\|S_{src}\|dt\ge I_1/(44T)
}
\]

or

\[
\boxed{
\int\|A\|\,\|B\|dt\ge I_1/(1188T).
}
\]

Here `A=grad V`; the second quantity is deliberately called a
**velocity-gradient/frame coupling**, not a pure strain defect.  On the
H1-dominant branch, if the base-gradient action is below `1/2376`, one of
pressure-third, differentiated-SGS, or viscous-fourth source channels has
integrated norm at least `I_1/(132T)`.  Action `31198595416` passed `240` tests
plus `50,000` exact-algebra/source checks; hook-projector and physical product-rule
residuals were `3.472e-15` and `2.278e-15`.  The earlier failed source run used a
strict Arb comparison at an exact rational touching equality and is certificate
semantics, not a mathematical countermodel.

#### Filtered SGS / viscous source collision

The differentiated-SGS source is now reduced to the primitive LES/Onsager
observable without differentiating raw turbulence.  The transporter is

\[
V=S_{N/4}u,
\qquad \operatorname{supp}\widehat V\subset B_{N/4},
\]

so its SGS stress has

\[
\boxed{\operatorname{supp}\widehat R\subset B_{N/2}.}
\]

Vector-valued Hausdorff--Young/Bernstein, the affine factor
`N||L^-1||||L||^2<=cond(L)^2 Nr_g`, and the mild-aspect bound give

\[
\boxed{
N^{-4}\|S_R\|
\le \frac3{2000}s\|R\|_{3/2},
\qquad s=Nr_g.
}
\]

Thus on `s<=s_0`, a normalized differentiated-SGS source `rho_R` forces

\[
\boxed{
\|R\|_{3/2}
\ge\frac{2000}{3s_0}\rho_R,
}
\]

and the exact Germano increment identity gives a cubic increment charge
proportional to `rho_R^(3/2)` at the actual filter scale.  If instead `s>s_0`,
the already-certified affine grain itself has

\[
\boxed{
N\int_{E_2}|u|^2\ge \frac3{10}s_0,
}
\]

which is a radius-energy/ancestry event, not an aspect penalty.

For viscosity, band support of `V` gives

\[
\boxed{
N^{-4}\|S_\nu\|
\le\frac{\nu s}{5000}\sqrt{\mathfrak d_V},
\qquad
\mathfrak d_V=N^{-1}\|\nabla V\|_2^2.
}
\]

Hence on `s<=s_0`, source level `rho_nu` forces

\[
\boxed{
\mathfrak d_V
\ge
\left(\frac{5000\rho_\nu}{\nu s_0}\right)^2.
}
\]

Action `31226401666` passed `246` tests plus `50,000` checks; the worst affine
factor ratio was `0.999999860` and the minimum SGS/viscous collision margins
were `4.789e-08` and `2.377e-08`.

#### Onsager cubic increments: grain, entropy, or dissipation

Let `u=sum_j u_j`, `N_j=2^jN`, and

\[
\mu_j=N_j\|u_j\|_2^2.
\]

Standard 3D `L^3` Littlewood--Paley square-function plus band Bernstein gives,
for filter support `|r|<=R_G/N`,

\[
\boxed{
Q_N
\le
g_1(C_{LP}C_B)^3
\left[
\frac43(\beta R_G)^2\mu_{\le N}^{max}
+2\mathfrak d_{>N}
\right]^{3/2},
}
\]

where

\[
\mathfrak d_{>N}
=\sum_{j\ge1}2^j\mu_j
\simeq N^{-1}\|\nabla P_{>N}u\|_2^2.
\]

For

\[
X=\left[Q_N/(g_1(C_{LP}C_B)^3)\right]^{2/3},
\]

one has the exact low/high collision

\[
\boxed{
\mu_{\le N}^{max}
\ge\frac{3X}{8(\beta R_G)^2}
\quad\text{or}\quad
\mathfrak d_{>N}\ge X/4.
}
\]

There is no logarithmic loss in the number of dyadic bands.  If the winning
aggregate band is packetized with normalized masses `w_a`, choosing
`theta=1/4`, `alpha=1/2` gives

\[
\boxed{
\mu_{atom}\ge\frac14\mu_{band}
}
\]

or

\[
\boxed{H_{anc}\ge\log2}
\]

or

\[
\boxed{Q_{anc}-Q_{at}\ge\frac14.}
\]

Thus cubic increments route into a dominant fresh/reused reservoir, uniform
Bellman entropy / ancestry cycles, or high-frequency enstrophy.  Action
`31226526751` passed `246` tests (`5` optional dependency skips) plus `50,000`
checks.  Its status is exact sequence algebra given the standard LP/Bernstein
analytic input.

The source-to-increment power `rho_R^(3/2)` and the increment collision power
`2/3` cancel exactly.  Therefore source episodes are best partitioned by
**source weight**, not time measure.  If one H1 physical source channel carries

\[
\Sigma_*=\int_0^c\rho_*(\tau)d\tau,
\]

then `Sigma_*>=I_1/(132c)` on the original mild strip and
`Sigma_*>=I_1/(1800c)` on the extended strip below.  Outside a large-radius
branch carrying half the source weight, either the mass or enstrophy branch
carries at least one quarter.  Consequently

\[
\boxed{D_{>N}\ge(c_d/4)\Sigma_*}
\]

on the high-frequency branch.  For viscosity, `d_V>=b rho_nu^2` and Cauchy on
the whole scaled lifetime gives

\[
\boxed{D_V\ge b\Sigma_*^2/(4c).}
\]

Thus temporal concentration is not a free SGS/viscous escape; for viscosity it
only raises the quadratic dissipation price.  Preferred Action `31227535026`
passed `264` tests plus `50,000` checks.  Run `31226868712` failed only because
an algebraically identical viscous expression was checked with an absolute
floating tolerance after a very large coefficient; runs `31227448616` and
`31227450217` failed only because an extended-source parameter was wired into a
function body before its signature.  These are implementation provenance, not
countermodels.

#### Extended aspect bridge and sticky anisotropy ancestry

The physical H1 bridge does not stop abruptly at `21/20`.  The perturbative
formula

\[
\sqrt{Q_{pol}}
\ge
\left[\frac1{\sqrt{10}}-\sqrt5(\kappa-1)\right]\|B^H\|
\]

remains positive on the exact transition strip

\[
\boxed{\kappa=\operatorname{cond}L\le567/500.}
\]

Arb certifies

\[
\boxed{Q_{pol}\ge\frac1{4000}\|B^H\|^2,}
\qquad
\boxed{\sum_i\|F_i^{H1}\|^2\ge\frac1{16000}\|B^H\|^2.}
\]

After the same physical nonunitary Duhamel bookkeeping,

\[
\boxed{
\operatorname{Def}\ \text{or}\ R_{pair}
\ge I_1^2/28\,800\,000,
}
\]

and the full-curvature transition-strip clean constant is

\[
\boxed{
\operatorname{Def}\ \text{or}\ R_{pair}
\ge I_B^2/115\,200\,000
}
\]

outside source/feedback/large-daughter/rescue branches.

High aspect beyond that strip is not charged.  Covariance dynamics gives, with
`S=sym A`,

\[
\boxed{
\frac d{dt}\log\operatorname{cond}L
\le2\|S\|_{op};
}
\]

viscosity contributes nonpositively.  On `int||S||dt<=1/30`,
`exp(1/15)<27/25`, so

\[
\boxed{
\kappa(t_1)>567/500
\Longrightarrow
\kappa(t_0)>21/20
}
\]

for the same transported covariance ancestry.  Strong anisotropy is therefore
sticky/inherited; if the predecessor identification breaks, the occurrence is
fresh/relinked.

For a fresh grain the shell lower axes imply

\[
\boxed{s=Nr_g>\frac23\kappa^{1/3},}
\qquad
\boxed{N\int_{E_2}|u|^2>\frac15\kappa^{1/3}.}
\]

Again this is **not** an aspect deficit: shell concentration forbids arbitrarily
small squeeze, so fresh high aspect can only be produced by elongation and
larger physical radius occupancy.  The extended source calculus gives
`J_1>=I_1/(132T)` -> `int||S_src||>=I_1/(600T)` or
`int||A||||B||>=I_1/(16000T)`; the latter is a conservative full
velocity-gradient/frame coupling, not a pure strain defect.  If that base
coupling action is below `1/32000` on the H1-dominant branch, one physical source
channel is at least `I_1/(1800T)` and feeds the same source-weighted collision
machinery.  Preferred Action `31227536542` passed `264` tests plus `50,000`
checks; viscosity's sampled contribution to anisotropy growth was always
negative.

#### Ancestor-reservoir synchronization: spectral half-life

The low/base reservoir produced by the Onsager collision can no longer be reused
indefinitely at no price.  On the signed-good core,

\[
\boxed{
{N_{q+1}\over N_q}
>e^{\gamma_*-1/80}>\frac85.
}
\]

For one materially reused reservoir covector under the same affine/Kelvin flow,

\[
\dot k=-A^Tk,
\qquad
{d\over dt}\log|k|\le\|\operatorname{sym}A\|_{op}.
\]

Thus on the existing low-strain branch `int||sym A||dt<=1/30`, Arb certifies

\[
\boxed{{M_{q+1}\over M_q}<e^{1/30}<\frac{21}{20}.}
\]

The low-band contribution to the LP increment square is proportional to
`(M/N)^2 mu_M`.  Therefore its coefficient per unit critical mass contracts by

\[
\boxed{
\left({21/20\over8/5}\right)^2
={441\over1024}<\frac12.
}
\]

Since `mu_M=M E_M`, the coefficient per unit physical energy contracts by

\[
\boxed{
\left({21\over20}\right)^3\left({5\over8}\right)^2
={231525\over512000}<\frac12.
}
\]

Even if the same old reservoir is adversarially allowed the entire conserved
energy at every future service time, its service capacity is geometrically
summable and totals less than twice its generation-zero energy-capacity
coefficient.  Material reuse is not merely a label: for `Ldot=A L` and
`kdot=-A^T k`,

\[
\boxed{q=L^Tk=\text{constant}.}
\]

A new carrier/band not obtainable from the transported material state is a
relinking/source event, not free reuse.  Action `31228679822` on `d3f8a51`
passed `268` tests plus `50,000` checks; the worst Kelvin-covector residual was
`1.490e-15`, and the largest sampled critical-mass / physical-energy service
ratios were `0.425474357` and `0.444631270`.

#### Band-limited source sampling and resolvable source clusters

Let a same-scale source family have centers separated by `delta/N`.  For a fixed
band limit, the standard Plancherel--Polya theorem gives

\[
\sum_a|f(x_a)|^p
\le C_{PP}N^3\|f\|_p^p.
\]

Combining this with the affine factor
`||L^-1||||L||^2<=kappa^2 r_g` makes the Navier--Stokes powers cancel exactly.
For differentiated SGS source,

\[
\boxed{
\sum_a\rho_{R,a}^{3/2}
\le
C_{PP}(\kappa_0^2s_0C_{D3})^{3/2}
\|R\|_{3/2}^{3/2}
\lesssim Q_{inc}.
}
\]

The exponent identity is

\[
-5\cdot\frac32+3+3\cdot\frac32=0.
\]

For viscous-fourth source,

\[
\boxed{
\sum_a\rho_{\nu,a}^2
\le
C_{PP,2}(\kappa_0^2s_0\nu C_{41})^2
\mathfrak d_V,
}
\]

with the exact scaling identity `-10+3+6+1=0`.  Hence separated source
replication pays additive increment charge or dissipation.

A maximal `delta/N`-separated subfamily covers every remaining source center by
one bandwidth-resolution cluster.  If all selected centers have `rho>=rho_0`,
then

\[
\boxed{
\#\mathcal C_R\le Q_R/\rho_0^{3/2},
\qquad
\#\mathcal C_\nu\le Q_\nu/\rho_0^2.
}
\]

Sub-resolution multiplicity is therefore merged rather than counted as many
independent source events.  Preferred Action `31228947363` on `c4c87f6` passed
`267` tests (`7` optional skips) plus `50,000` checks; the scaling residual was
`6.262e-16`.  The theorem status is exact scaling/routing given the standard
Plancherel--Polya, Bernstein and Riesz analytic inputs.

#### Pressure-third near field and quadratic reservoir synchronization

For the strict filtered pressure,

\[
\boxed{
-\Delta P=\partial_i\partial_j(V_iV_j+R_{ij}),
}
\]

so `P` remains band limited.  Riesz boundedness plus low-pass Bernstein gives

\[
\boxed{
\|P\|_{3/2}^{3/2}
\le
\sqrt2 C_R^{3/2}
\left[C_B^3\mu_V^{3/2}+\|R\|_{3/2}^{3/2}\right].
}
\]

The preceding sampling theorem therefore routes separated pressure-third source
clusters to a resolved low-frequency mass reservoir or to the same cubic
increment charge.  The SGS part already enters the Onsager ledger.

For one fixed low-low material reservoir pair `(a,b)`, the pre-`3/2` normalized
pressure service coefficient per unit `sqrt(E_aE_b)` is

\[
\mathcal P_{ab}(N)
\sim
{M_{max}^3\sqrt{M_aM_b}\over N^3}.
\]

On a low-strain signed-good lineage it contracts by

\[
\boxed{
\left({21\over20}\right)^4
\left({5\over8}\right)^3
={194481\over655360}
<\frac13.
}
\]

Thus one old pressure pair has total future service less than `3/2` times its
initial pair energy-capacity coefficient, even when both endpoints are allowed
the full global energy cap at every time.  Persistent pressure service must
relink to new pairs, fragment into a pair-entropy/cycle problem, or leave the
low-strain material-reuse branch.  Action `31229118774` on `6e21fe6` passed
`270` tests (`7` optional skips) plus `50,000` checks; the maximum sampled pair
ratio was `0.288962935`.

#### Whole-old-pool erosion and forced relinking

The single-reservoir half-life extends to an **entire old reservoir pool** once
its packet family has a Bessel/frame energy budget.  Suppose at generation zero
all old atoms satisfy `M_(a,0)<=alpha N_0` and, at every later service time,

\[
\boxed{
\sum_aE_{a,q}\le P E_{global}.
}
\]

No orthogonality is assumed; `P` is an explicit frame constant.  Low-strain
material transport and signed-good progress give

\[
M_{a,q}\le(21/20)^qM_{a,0},
\qquad
N_q\ge(8/5)^qN_0.
\]

Since low-band increment service is `M^3E/N^2`, summing over all old atoms yields

\[
\boxed{
\mathsf C_{old}(q)
\le
\alpha^3N_0P E_{global}
\left({231525\over512000}\right)^q
<2^{-q}\alpha^3N_0P E_{global}.
}
\]

Hence

\[
\boxed{
\sum_{q\ge0}\mathsf C_{old}(q)
<2\alpha^3N_0P E_{global}.
}
\]

The energy distribution may be changed adversarially among old atoms at every
generation; only the total frame budget is used.  Thus a uniform positive service
threshold eventually forces **newly relinked spectral capacity**.

Once non-negligible relinking events are represented by the existing active
quadratic triad graph, the exact incidence Euler law for a connected component
with `m` triads, `n` packet vertices and cycle rank `beta` is

\[
\boxed{(n-1)+\beta=2m.}
\]

Therefore

\[
\boxed{n-1\ge m\quad\text{or}\quad\beta\ge m.}
\]

Relinking is fresh-rich or cycle-rich; there is no neutral third regime in the
finite atomic model.  Action `31229543390` on `8c50370` passed `274` tests
(`7` optional skips) plus `50,000` pool checks and `5,000` incidence checks.  Its
status is exact given the frame budget and the existing incidence theorem.

#### Affine coherent-state Riesz budget inside one covariance cell

The frame-budget input above is itself theorem-level for separated coherent
Gaussian synthesis packets in one covariance cell.  For normalized packets

\[
g_{X,k,L}(x)
=C_L e^{-|L^{-1}(x-X)|^2/4}e^{ik\cdot x},
\]

direct Gaussian integration gives

\[
\boxed{
|\langle g_{X,k,L},g_{Y,\ell,L}\rangle|
=
\exp\!\left[
-\frac18|L^{-1}(X-Y)|^2
-\frac12|L^T(k-\ell)|^2
\right].
}
\]

Define the affine-invariant coherent coordinate

\[
\boxed{
\zeta=(L^{-1}X/2,L^Tk)\in\mathbb R^6.
}
\]

Then the overlap is exactly `exp(-|zeta_a-zeta_b|^2/2)`.  For a `4`-separated
family, six-dimensional packing plus Schur gives the Arb-certified Bessel bound

\[
\sum_a|\langle f,g_a\rangle|^2
\le\frac{25}{4}\|f\|_2^2.
\]

More importantly, at separation `5`, Arb certifies

\[
\boxed{
\sup_a\sum_{b\ne a}|G_{ab}|<\frac3{50},
}
\]

hence

\[
\boxed{
\frac{47}{50}I\le G\le\frac{53}{50}I.
}
\]

Therefore **arbitrary synthesis coefficients** in this separated equal-covariance
family obey

\[
\boxed{
\sum_a|c_a|^2
\le\frac{50}{47}
\left\|\sum_ac_ag_a\right\|_2^2.
}
\]

Thus the old-pool erosion theorem has the explicit coefficient budget
`P=50/47` inside one 5-separated covariance cell.  Action `31229863366` on
`84ae98c` passed `286` tests plus `50,000` affine/pair checks and `3,000` finite
Gram probes.  Arb certified the `[47/50,53/50]` spectrum; sampled 5-separated
Gram eigenvalues lay between `0.999999172` and `1.000000828`.

This does not yet prove that an arbitrary transfer-selected decomposition can be
reduced to such separated coherent synthesis cells with summable loss.  Close
phase-space clusters and covariance-cell changes must be registered as
relinking/interface structure rather than silently re-expanded with new
coefficients.

#### Exact coherent Moyal reservoir-energy ledger

The frame-budget problem has an even more canonical analysis-side formulation.
For any `L^2`-normalized affine Gaussian window `g_L`, define the coherent/STFT
transform

\[
\mathcal V_L f(X,k)
=
\int f(x)\overline{g_L(x-X)}e^{-ik\cdot x}\,dx.
\]

Plancherel in `k`, followed by Fubini and translation invariance in `X`, gives the
exact Moyal identity

\[
\boxed{
\int_{\mathbb R^6}|\mathcal V_Lf(X,k)|^2
\frac{dX\,dk}{(2\pi)^3}
=
\|f\|_2^2.
}
\]

Thus any measurable coherent phase-space partition has positive cell energies

\[
E_C=\int_C|\mathcal V_Lf|^2d\mu,
\qquad
\boxed{E_C\ge0,\quad\sum_CE_C=\|f\|_2^2.}
\]

At the continuous analysis level the reservoir frame constant is therefore
exactly

\[
\boxed{P=1.}
\]

This avoids coefficient cancellation entirely: close coherent probes are merged
inside phase-space cells rather than counted as independent synthesis atoms.  For
an orthogonal dyadic old-band decomposition, Moyal supplies the whole-old-pool
energy budget directly; a smooth LP partition pays only its standard fixed
square-function overlap constant.  The 5-separated Riesz theorem remains the
complementary discrete synthesis realization.

Action `31230201131` on `975f680` passed `282` tests (`8` optional skips) plus
`5,000` periodic discrete-Moyal regressions.  The worst relative Moyal and cell
partition residuals were `6.061e-16` and `5.787e-16`.  The continuum identity is
exact by Plancherel/Fubini, not a numerical approximation.

The remaining bridge is now specifically transfer-facing: prove that the
low/base **service selected by the nonlinear Navier--Stokes block** can be
assigned to these positive coherent cells through recursive generations, with
cell changes and rejected interfaces entered into `Xi`/fresh/relink bookkeeping.

#### Coherent nonlinear work and material phase-space cells

The nonlinear transfer-to-cell assignment is now exact at the work level.  For a
normalized coherent window `g`, polarized Moyal gives

\[
\boxed{
\int \mathcal V_g f(z)\overline{\mathcal V_gF(z)}\,d\mu(z)
=\langle f,F\rangle.
}
\]

For any measurable phase-space partition `C_alpha`, define

\[
\boxed{
\mathcal W_\alpha
=2\Re\int_{C_\alpha}\mathcal V_gf\,
\overline{\mathcal V_gF}\,d\mu.
}
\]

Then

\[
\boxed{\sum_\alpha\mathcal W_\alpha=2\Re\langle f,F\rangle.}
\]

Thus a Leray-projected nonlinear band forcing, source term, or other physical
work term has a canonical coherent-cell decomposition without choosing
redundant synthesis coefficients.  Positive cells are actual positive service;
negative cells are actual backflow/cancellation.

Common affine transport is again a gauge, not an interface.  Under

\[
L\mapsto ML,\qquad X\mapsto MX,\qquad k\mapsto M^{-T}k,
\]

the intrinsic coherent coordinate

\[
\boxed{\zeta=(L^{-1}X/2,L^Tk)}
\]

is unchanged and the phase-space Jacobian is one.  Therefore material coherent
cells create no common-affine Heisenberg/interface forcing.

If a piecewise material selected cell family is changed, positivity of the Moyal
cell energies gives

\[
|E(S_+)-E(S_-)|\le E(S_+\triangle S_-).
\]

Integrating the exact cellwise work balance yields

\[
\boxed{
P_+\le E_{final}+P_-+R_{switch},
}
\]

where `P_+` and `P_-` are integrated positive/negative coherent work and
`R_switch` is the symmetric-difference Moyal mass of selection changes.  Hence
one of terminal coherent energy, physical backflow/cancellation, or relinking
mass is at least `P_+/3`.

Preferred Action `31245848657` passed `293` tests (`8` optional skips) plus
`20,000` work/affine/relink checks; polarized-Moyal and cell-work residuals were
`3.265e-16` and `6.209e-16`.  The first rerun on the same code SHA failed only
because new coherent-increment unit fixtures elsewhere in the full pytest suite
did not satisfy their own service hypothesis; this is test-fixture provenance,
not a work-cell countermodel.

#### Actual SGS increments generate coherent service edges

The low/base reservoir selected by the Onsager collision can now be assigned to
phase space using the **increment itself**, rather than aggregate band mass.
Let

\[
Q=\int |G_N(r)|\|\delta_ru\|_3^3dr,
\qquad
Y={ (Q/g_1)^{2/3}\over(C_{LP}C_B)^2}.
\]

Some filter displacement `r` satisfies the corresponding `L^3` square lower
bound.  Standard LP/Bernstein and the existing high-frequency enstrophy currency
give

\[
\boxed{
S_{low}(r):=\sum_{j\le0}M_j\|\delta_ru_j\|_2^2
\ge Y-2d_{high}.
}
\]

In particular, if `d_high<Y/4`, then

\[
\boxed{S_{low}(r)\ge Y/2.}
\]

For every low band and coherent cell define the positive service edge weight

\[
\boxed{
s_{j,C}(r)
=M_j\int_C|\mathcal V_{g_j}\delta_ru_j|^2d\mu.
}
\]

Moyal gives `sum_(j,C) s_(j,C)=S_low`.  Translation covariance is exact:

\[
\boxed{
\mathcal V_g(\delta_ru)(X,k)
=e^{-ik\cdot r}\mathcal V_gu(X-r,k)-\mathcal V_gu(X,k).
}
\]

Hence

\[
\boxed{
s_{j,C}(r)
\le2M_j\{E_j(C)+E_j(C-r)\}.
}
\]

Each service atom is therefore a physical coherent edge joining neighborhoods
separated by the actual SGS displacement `r`.

Classify these edges relative to the transported old pool.  Once whole-old-pool
erosion gives `C_old<=Y/8`, the clean routes are

\[
\boxed{d_{high}\ge Y/4}
\]

or

\[
\boxed{\Xi_{cell}\ge Y/8}
\]

or at least `Y/4` new--new coherent service remains.  With the quarter-dominance
threshold, the last branch gives

\[
\boxed{\mu_{coh,new}\ge Y/32}
\]

for one actual coherent cluster, or else the existing collision chain yields

\[
\boxed{H_{ancestry}\ge\log2}
\qquad\text{or}\qquad
\boxed{Q_{ancestry}-Q_{atomic}\ge1/4.}
\]

Action `31245850938` passed `293` tests (`8` optional skips) plus `50,000`
routing checks.  The worst STFT translation-covariance residual was `6.829e-14`;
the largest sampled local increment-capacity ratio was `0.855742519<1`.  Runs
`31245808950` and `31245810803` failed only because three hand-written fixtures
supplied less total cell service than their theorem hypothesis
`S_low>=Y-2d_high`; the router correctly rejected them.

#### Covariance-cell changes have a deterministic `Xi` bound

Changing the Gaussian analysis covariance is also controlled explicitly.  For
centered normalized Gaussian windows with physical covariances `Sigma,Theta`,

\[
\boxed{
\langle g_\Sigma,g_\Theta\rangle
=2^{3/2}{(\det\Sigma\det\Theta)^{1/4}\over
\det(\Sigma+\Theta)^{1/2}}.
}
\]

If `exp(a_i)` are the generalized covariance eigenvalues, this equals
`prod_i cosh(a_i/2)^(-1/2)`.  The inequality `log cosh x<=x^2/2` gives

\[
\boxed{
\|g_\Sigma-g_\Theta\|_2
\le {1\over2\sqrt2}
\|\log(\Sigma^{-1/2}\Theta\Sigma^{-1/2})\|_F.
}
\]

Moyal in the **window slot** is exact:

\[
\|\mathcal V_{g_\Sigma}f-\mathcal V_{g_\Theta}f\|_2
=\|f\|_2\|g_\Sigma-g_\Theta\|_2.
\]

Therefore the total variation of the positive coherent energy law obeys

\[
\boxed{
\int\left||\mathcal V_{g_\Sigma}f|^2-|\mathcal V_{g_\Theta}f|^2\right|d\mu
\le {1\over\sqrt2}d_{log}(\Sigma,\Theta)\|f\|_2^2.
}
\]

The same bound holds after any common phase-space partition.  Small changes of
the covariance representative therefore enter a deterministic
`Xi_cov`; large covariance jumps remain a genuine covariance/strain/source or
fresh-relink event rather than being hidden as representation gauge.

Preferred Action `31246026975` passed `297` tests (`8` optional skips) plus
`50,000` noncommuting-SPD/window checks.  The overlap formula residual was
`1.438e-14`, the clean distance ratio reached `0.999320817`, and the window-Moyal
residual was `7.042e-16`.  Run `31245962226` failed only because its scalar sanity
fixture used an impossible Gaussian overlap `>1`; the matrix formula itself was
correct.

#### Coherent service stopping and source-weighted SGS ancestry

The preceding cell theorem and the old-pool half-life give a genuine stopping
statement.  Let

\[
r={231525\over512000}<\frac12,
\qquad
C_{old}(q)\le C_0r^q.
\]

If every block in one sticky ancestry epoch has `Y>=Y_0`, then a block which
avoids dissipation, interface, fresh-coherent, entropy and cycle exits requires
`C_old(q)>Y_0/8`.  Hence the cost-free epoch ends by the first `q_*` satisfying

\[
\boxed{C_0r^{q_*}\le Y_0/8.}
\]

A relink may restart the clock, but the relink itself is already a charged exit.
Action `31246194276` passed `301` tests (`8` optional skips) plus `50,000`
stopping checks; the largest sampled cost-free epoch length was `15`.

For the differentiated-SGS H1 source, even the pointwise threshold assumption is
unnecessary.  If

\[
\rho_R(\tau)=N^{-4}\|S_R(\tau)\|,
\qquad
\Sigma_R=\int\rho_Rd\tau,
\]

then the already-certified source collision gives
`Q_inc>=c_Q rho_R^(3/2)`.  The coherent increment threshold takes the `2/3`
power, so

\[
\boxed{Y_{coh}(\tau)\ge c_Y\rho_R(\tau)}
\]

with an explicit `c_Y`.  After the existing large-radius/scale-matched split,
one has

\[
\boxed{D_{high}\ge c_Y\Sigma_R/16}
\]

or

\[
\boxed{\int S_{low}d\tau\ge c_Y\Sigma_R/8.}
\]

If the old material pool has integrated capacity

\[
cC_{old}\le c_Y\Sigma_R/32,
\]

then either

\[
\boxed{\int\Xi_{cell}d\tau\ge c_Y\Sigma_R/32}
\]

or

\[
\boxed{\int S_{new-new}d\tau\ge c_Y\Sigma_R/16.}
\]

A quarter-dominant new service edge gives

\[
\boxed{\int\mu_{coh,new}d\tau\ge c_Y\Sigma_R/128,}
\qquad
\boxed{\sup_\tau\mu_{coh,new}\ge c_Y\Sigma_R/(128c),}
\]

while fragmentation pays `log 2` ancestry entropy or `1/4` same-ancestry
pair/cycle mass.  Thus temporal concentration is not a coherent-cell escape.
Action `31246242892` passed `303` tests (`8` optional skips) plus `50,000`
source-weight checks; the source-homogeneity residual was `8.938e-16`.

#### Uniform service/action-or-Kelvin-flat rigidity

There is now a single positive packet-level gate joining the physical currency
ledgers to a gauge-quotiented flatness statement.  For

\[
0<\tau\le1/10
\]

define

\[
\boxed{
\delta_\tau={\tau^2\over1\,036\,800\,000}.
}
\]

The denominator is `36*28,800,000`; the limiting branch is the extended-aspect
physical H1 pair/deficit theorem, not H3 or Hodge.

On the signed-good physical low-cost core, inside
`cond(L)<=567/500`, on the low-strain packet branch, and outside the already
named nonlinear-feedback / large-daughter exits, suppose

\[
\operatorname{Def}_{block},\quad
\operatorname{Def}_{sideband},\quad
R_{pair}<\delta_\tau,
\]

\[
TJ_3<\frac{\tau}{6\sqrt6},
\qquad
TJ_1<\frac{\tau}{792},
\]

and

\[
T\int\|\dot D_{obj}\|dt<\frac{\tau}{60}.
\]

Then the three independent rigidity sectors satisfy

\[
\boxed{\sqrt{E_H^{phys}}\le\tau/3},
\]

\[
\boxed{(dT)_{nonconf}\le\tau/3},
\]

and

\[
\boxed{I_B\le\tau/3}.
\]

Thus the gauge-quotiented Kelvin connection flatness

\[
\boxed{
\mathfrak F_K
=\sqrt{E_H^{phys}}+(dT)_{nonconf}+I_B
\le\tau.
}
\]

Common translation, rigid rotation, common affine/Kelvin transport,
common-parent `SL(2)` polarization motion and material coherent-cell transport
are not charged.

The full-curvature step uses the irreducible channel dichotomy rather than
adding H1 and H3 constants: if full curvature is non-flat, either
`I_3>=I_B/(2 sqrt(6))` or `I_1>=I_B/2`.  The H1 alternative fixes the clean
constant above.

Contrapositively, failure of `tau`-flatness forces a uniformly positive named
branch: transfer/pair cost, H1/H3 source, objective-strain source action,
high-strain lifetime, fresh/inherited high aspect, nonlinear daughter/feedback,
or helical phase holonomy.  Source branches feed the already-certified coherent
SGS/pressure/viscous ancestry or dissipation ledgers.

At the concrete target

\[
\boxed{\tau=1/100},
\]

the common transfer/pair threshold is

\[
\boxed{
\delta_{1/100}=1/10\,368\,000\,000\,000.
}
\]

Preferred Action `31246845056` passed `316` tests plus `50,000` branch checks;
the maximum sampled `F_K/tau` was `0.563215647` and every named branch was hit.
The small constant is deliberately conservative; the theorem's role is a
uniform positive gate, not optimization of the final exponent.

#### Physical flat episode: the master perturbation is explicit

The spherical master previously allowed an abstract near-flat perturbation
`zeta_j`.  Signed-good physical triad geometry now computes it.

With

\[
x=r_*e^{-v-u/2},\qquad y=r_*e^{-v+u/2}
\]

and `c_e=cos(theta_e/2)`, the law of cosines gives exactly

\[
\boxed{
c_e^2=c_*^2e^{2v}-\sinh^2(u/2),
\qquad c_*=(2r_*)^{-1}.
}
\]

On `|u|,|v|<=2/25`, Arb certifies

\[
\boxed{|c_e-c_*|\le\frac65|v|+\frac15u^2.}
\]

For unit parent directions `p,q`, exact unequal-parent geometry gives

\[
\tan\angle(n_c,m)
=\tan(\theta_e/2)|\tanh(u/2)|,
\]

where `m=(p+q)/|p+q|`.  The local half-angle is above `3/4`, hence

\[
\boxed{\|n_c-m\|\le|u|/2.}
\]

For normalized positive physical transfer weights let

\[
H=E_H^{phys}=\mathbb E(2v^2+u^2/2),
\qquad
\Delta_b=|b_2-b_1|.
\]

Then

\[
\boxed{
\left|b_c-\frac{b_1}{c_*}\right|
\le
2\sqrt H+\frac12H+\frac58\Delta_b.
}
\]

On the concentrated branch `|b_1|>=c_*`, the barycentric potential
`P=-log|b|` obeys

\[
P_c\le P_1-\kappa_*-\log(1-e),
\]

with the displayed right-hand error `e`.  For `e<=1/2` the master perturbation
may be taken as

\[
\boxed{
\zeta
\le4\sqrt{E_H^{phys}}+E_H^{phys}+\frac54\Delta_b.
}
\]

At `tau=1/100`, if parent marginals are synchronized at `Delta_b<=tau`, the
service-or-flat theorem gives `sqrt(H)<=tau/3`, and therefore

\[
\boxed{
\zeta\le\frac{31}{12}\tau+\frac{\tau^2}{9}.
}
\]

Arb then certifies the **physical flat erosion rate**

\[
\boxed{
\kappa_0:=\kappa_*-\zeta>17/100.
}
\]

Action `31247033855` passed `319` tests plus `50,000` physical-coupling checks;
the worst actual barycenter-error / clean-bound ratio was `0.956538235`, and the
clean `kappa_0>0.17` margin was `4.487e-03`.

#### Equal parent marginals are not required: entropy or companion cores

The synchronized-marginal assumption itself is no longer a mandatory input.
On a `tau=1/100` Kelvin-flat block,

\[
H=E_H^{phys}\le\tau^2/9,
\]

and, without assuming `b_1=b_2`, the same triad geometry gives

\[
\boxed{
\left|b_c-\frac{b_1+b_2}{2c_*}\right|
\le2\sqrt H+\frac12H.
}
\]

Choose

\[
\boxed{\beta=99/100.}
\]

If either parent barycenter has norm at most `beta`, the atomic
barycenter--collision theorem gives

\[
H_2\ge\log\frac{200}{199}.
\]

Arb certifies the clean floor

\[
\boxed{H_2>1/200.}
\]

If both barycenter norms exceed `0.99`, physical `|b_c|<=1` forces their
barycenter directions to be separated by more than one radian.  Since a unit-
direction law satisfies

\[
\mathbb E|X-\widehat b|^2=2(1-|b|),
\]

Markov at chord radius `3/10` puts at least `7/9` of each parent transfer mass in
one directional cap.  The two caps have chord gap greater than `1/3`.
Consequently

\[
\boxed{
\text{1\% Kelvin-flat}
\Longrightarrow
H_2\ge1/200
\quad\lor\quad
\text{two separated parent cores, each mass }\ge7/9.
}
\]

Relative to a distinguished old lineage the second core is a **trackable
companion**, which must be classified by ancestry rather than silently replaced
by an equal-marginal symmetry assumption.

Preferred Action `31247281011` passed `322` tests plus `50,000` checks.  The
entropy, one-radian angle and `1/3` cap-gap clean constants had margins
`1.461e-05`, `1.754e-01` and `1.756e-01`.  Run `31247206213` failed only because
`7/9` was compared to its slightly lower floating representation with exact
`>=`; no theorem constant failed.

#### Causal binary ancestry: fresh companions force sticky reuse

A fresh companion at every generation is **not** ruled out by energy alone:
critical packet costs `~1/N_q` are summable.  The missing physical input is
quadratic causality.  A genuinely new comparable-frequency companion must
itself have parents in the preceding causal layer unless it is imported from an
already-existing reservoir.

The following finite theorem is now exact under a synchronized causal layering.
On the signed-good core, Arb certifies the two-sided scale window

\[
\boxed{
\frac35<\frac{N_{parent}}{N_{child}}<\frac58,
}
\]

hence

\[
\boxed{
\frac85<\frac{N_{child}}{N_{parent}}<\frac53.
}
\]

Let `n_j` be the number of distinct packet ancestors at level `j`, with one
terminal packet `n_L=1`, and every child at level `j+1` generated by one active
triad with two parent incidences at level `j`.  Define

\[
\boxed{r_j=2n_{j+1}-n_j\ge0},
\qquad
\boxed{\rho_j={r_j\over2n_{j+1}}}.
\]

For the connected layered 3-uniform incidence graph the cycle rank satisfies the
exact identity

\[
\boxed{\beta_{cycle}=\sum_{j=0}^{L-1}r_j.}
\]

Moreover binary branching has the exact multiplicative reuse law

\[
\boxed{
{n_0\over2^L}
=\prod_{j=0}^{L-1}(1-\rho_j),
}
\]

so

\[
\boxed{
\mathcal A_{reuse}
:=\sum_j-\log(1-\rho_j)
=L\log2-\log n_0.
}
\]

Thus each parent-slot identification is simultaneously a binary-branch saving
and a cycle-rank unit.

Let `N_base` be the distinguished lineage scale at the root generation.  The
signed-good forward/backward windows imply every depth-`L` root ancestor has

\[
N_{root}<\left(\frac{25}{24}\right)^L N_{base}.
\]

If distinct coherent root cells obey a positive energy/frame budget `P E_global`
and every selected affine grain has critical mass `N E>=eta`, then

\[
\boxed{
n_0\le{P E_{global}N_{base}\over\eta}
\left(\frac{25}{24}\right)^L.
}
\]

Hence

\[
\boxed{
\mathcal A_{reuse}
\ge
L\log\frac{48}{25}
-\log\frac{P E_{global}N_{base}}{\eta}.
}
\]

For Moyal energy `P=1` and the clean affine-grain mass `eta=1/5`, if every layer
had `rho_j<1/4`, then `A_reuse<L log(4/3)`.  Therefore

\[
\boxed{
L\log\frac{36}{25}>\log(5E_{global}N_{base})
\Longrightarrow
\max_j\rho_j\ge1/4.
}
\]

A sufficiently deep causally complete ancestry must therefore contain a
**reuse-rich layer** in which at least one quarter of formal binary parent slots
are identified with already-used parents.  This is the causal sticky mechanism
that energy alone was missing: binary parent creation outruns the allowable
signed-good scale dilation.

Preferred Action `31247674599` passed `326` tests plus `50,000` layered graphs;
cycle identity residual was exactly `0`, binary-product residual `1.110e-16` and
reuse-action telescope residual `5.329e-15`.  Run `31247611869` failed only
because a hand-written reuse-count fixture expected cycle rank `2` for the count
vector `[4,3,2,1]`; the exact value is `3`.

**Scope:** this causal theorem is deliberately conditional on a synchronized
layered Duhamel/packet ancestry with a coherent root-cell budget.  Constructing
that layering for a genuine asynchronous Navier--Stokes packet extraction, with
summable interfaces, remains an open continuum bridge.  Raw cycle rank alone is
also not automatically physical flux weight; a reuse-rich layer must still be
routed through the transfer-weighted Hodge/Bellman/resistance machinery.

#### Transfer-weighted causal reuse information

The raw cycle-count caveat can itself be reduced to a physical transfer law.
Let `w_(j+1)` be the normalized causal-transfer distribution on child packets of
one synchronized layer.  Duplicate each event into the two structural parent
roles with joint probability

\[
\boxed{p(c,s)=w(c)/2,\qquad s=1,2,}
\]

and push this law through the physical parent-packet / coherent-cell label map.
The two-role `log 2` increase is exact baseline causality and is **not** charged.
The actual transfer-weighted reuse information is

\[
\boxed{
\mathcal R_j
=H(C,S\mid P)
=H(w_{j+1})+\log2-H(w_j)\ge0.
}
\]

For a one-terminal depth-`L` ancestry this telescopes exactly:

\[
\boxed{
\sum_j\mathcal R_j
=L\log2-H(w_0).
}
\]

Since `H(w_0)<=log n_0`, the same Moyal root-energy and signed-good scale bound
gives

\[
\boxed{
\sum_j\mathcal R_j
\ge
L\log\frac{48}{25}
-\log\frac{P E_{global}N_{base}}{\eta}.
}
\]

Therefore

\[
L\log\frac{36}{25}
>
\log\frac{P E_{global}N_{base}}{\eta}
\]

forces one layer with

\[
\boxed{\mathcal R_j>\log(4/3).}
\]

This is already weighted by the physical causal transfer law, so a reuse signal
cannot be hidden entirely in negligible combinatorial edges.  Action
`31248342830` passed `334` tests plus `50,000` layered laws; the information
telescope residual was `2.887e-15`.

#### Rényi causal reuse reaches existing Bellman/cycle currencies

The collision version closes the finite weighted-routing gap.  If

\[
Q_{child}=\sum_cw_c^2,
\]

the two parent-role slots have the exact cost-free collision baseline

\[
Q_{slot}=Q_{child}/2.
\]

After pushforward through physical parent labels,

\[
\boxed{
Q_{parent}
=Q_{child}/2+R_{hidden},
}
\]

where `R_hidden` is the transfer-weighted hidden pair mass of distinct causal
slots sharing the same parent.  Put

\[
\boxed{
\theta_j={2R_{hidden,j}\over Q_{child,j}}.
}
\]

Then

\[
\boxed{
Q_j=\frac12Q_{j+1}(1+\theta_j)
}
\]

and, because the terminal law is one atom,

\[
\boxed{
\sum_j\log(1+\theta_j)
=L\log2+\log Q_0.
}
\]

The root support/energy bound again yields

\[
\boxed{
\sum_j\log(1+\theta_j)
\ge
L\log\frac{48}{25}
-\log\frac{P E_{global}N_{base}}{\eta}.
}
\]

Thus the clean causal-depth condition forces a layer with

\[
\boxed{\theta_j>1/3.}
\]

At such a layer there is an **existing-currency** alternative.  If

\[
H_{2,child}<\log2,
\]

then `Q_child>1/2` and

\[
\boxed{R_{hidden}>1/12.}
\]

If instead

\[
H_{2,child}\ge\log2,
\]

apply the already-certified atomic-to-ancestry collision chain.  Either

\[
\boxed{H_{2,ancestry}\ge\tfrac12\log2}
\]

or

\[
\boxed{
R_{same\ ancestry}
\ge2^{-1/2}-1/2>1/5.
}
\]

Hence a sufficiently deep synchronized causal ancestry has the clean physical
transfer-weighted route

\[
\boxed{
\text{component/Bellman entropy}
\ \lor\
\text{same-ancestry pair/cycle}>1/5
\ \lor\
\text{parent-slot reuse pair}>1/12.
}
\]

No cost is assigned to the mere existence of the two parent roles; only
collision amplification above the exact factor `1/2` is reuse.  Preferred
Action `31248345125` passed `334` tests plus `50,000` causal laws; the layer
recurrence residual was `1.110e-16` and the multiplicative-action telescope
residual `2.665e-15`.

**Remaining scope:** the finite weighted causal algebra is now closed.  The open
bridge is constructing the synchronized causal layers from the actual
asynchronous Navier--Stokes Duhamel/packet extraction and making the parent /
ancestry labels used above coincide with the coherent/nested labels used by the
master, with all layer cuts paid once through `Xi`.

#### Spatial moat width must balance curvature

The earlier commutator estimate `O(1/M)` remains correct, and the quadratic
schedule `M_j=(j+3)^2` remains correct in the **defect-space nested-grain
extraction**.  It is not, however, a universal spatial moving-window schedule.
For a physical transition radius `R=M/N`, a parabolic lifetime `T=cN^{-2}`, and

\[
\kappa:=N^{-3}\|\nabla^2u\|_\infty,
\]

the affine-window Taylor error has the opposite dependence on `M`.  The
geometric localization ledger is

\[
\boxed{
E(M)=\frac aM+b\kappa M,
}
\]

where `a` is the smooth-filter commutator coefficient and `b` packages lifetime,
deformation, and the reference-window gradient.  The exact optimum is

\[
\boxed{
M_*=\sqrt{\frac{a}{b\kappa}},
\qquad
E_*=2\sqrt{ab\kappa}.
}
\]

There is a structural countermodel to using the fixed spatial schedule
`M_j=(j+3)^2` for both effects: if `kappa_j=(j+3)^(-3)`, then

\[
\sum_jM_j^{-1}<\infty,
\qquad
\sum_jM_j\kappa_j
=\sum_j(j+3)^{-1}=\infty.
\]

The curvature-balanced choice `M_j=kappa_j^{-1/2}=(j+3)^{3/2}` makes both terms
`(j+3)^(-3/2)` and hence summable.  Action `31168888413` verifies the exact
optimization and records this countermodel.

For a band-limited packet, curvature itself has a critical-mass meaning.  In the
unitary Fourier convention, for one scalar Hessian component and
`|xi|<=Lambda N`,

\[
N^{-3}\|\partial_{ij}u\|_\infty
\le C_B\sqrt\mu,
\qquad
\mu=N\|u\|_2^2,
\]

where

\[
C_B=(2\pi)^{-3/2}\sqrt{4\pi/7}\,\Lambda^{7/2}.
\]

Thus if even the optimally balanced moat costs at least `eta`, then in this
packet model

\[
\boxed{
\mu\ge\frac{\eta^4}{16a^2b^2C_B^2}.
}
\]

So excessive spatial curvature is itself a fresh/critical-mass candidate rather
than a neutral localization loss.

## 5. Gaussian packet inverse mechanism

The critical Fourier norm

\[
\|\widehat u\|_{L^{3/2}_\xi}
\]

is scale invariant under Navier–Stokes scaling. This matches the symmetric Young exponents `3/2,3/2,3/2`. The sharp scalar constant in `R^3` is

\[
A_3=\left(\frac{\sqrt3}{2}\right)^3\approx0.649519052838.
\]

For a weighted trilinear form with `|M|\le M_*`, near equality implies three independent conclusions:

1. the unweighted magnitudes nearly extremize sharp Young, hence are close to Gaussian extremizers after affine/translation symmetries;
2. multiplier deficit is small in the normalized interaction measure;
3. the combined transfer phase is locked.

For equal-width isotropic Fourier Gaussians, the scalar overlap factor is exact:

\[
\boxed{
R_{\rm scalar}=
\exp\!\left[
-\frac{|\kappa_1+\kappa_2-\kappa_3|^2}{12\sigma^2}
-\frac{\sigma^2}{3}\sum_{i<j}|x_i-x_j|^2
\right].
}
\]

Near equality therefore forces both Fourier resonance and physical-space meeting. Width imbalance has an explicit separate penalty; the nested-grain theorem below is intentionally stated for an equal-width block after width balancing.

The scale-critical Gaussian experiment shows a linear small-width deficit, matching the cusp prediction; the leading analytic coefficient is approximately `1.65`.

**Current inverse-theorem target:** near-maximal weighted transfer produces a coherent Gaussian triad grain with resonant centres, near-balanced widths, near-optimal helical signs/scale ratios, and locked phase.

## 6. Exact scale holonomy

For the reuse motif

- `a+b -> m`,
- `m+c -> d`,
- `b+c -> n`,

write `ℓ_v=log|v|` and residuals

\[
r_1=\ell_a-\ell_b,
\quad
r_2=\ell_m-\frac{\ell_a+\ell_b}{2}-\gamma_*,
\]

\[
r_3=\ell_m-\ell_c,
\quad
r_4=\ell_b-\ell_c.
\]

Then the following identity is exact:

\[
\boxed{
r_2-r_3+\tfrac12r_1+r_4=-\gamma_*.
}
\]

Hence not all reuse constraints can be simultaneously near zero. This is the first explicit Wang-style loop incompatibility: each edge separately wants the same rigid scale relation, but the loop cannot close.

Numerical optimization gives an observed 3-edge reuse min-ratio around `0.911617`, but this number is **not** used as a certified theorem constant.

## 7. Transfer-preserving component decomposition

Rather than decomposing functions first and later asking whether transfer decouples, the programme builds components from the interaction itself.

For a three-sided partition with critical masses `X_C,Y_C,Z_C`, normalized on each side, define

\[
S=\sum_C(X_CY_CZ_C)^{2/3}.
\]

Hölder gives the exact Bellman inequality

\[
\boxed{S\le1.}
\]

Two equal disconnected perfect components achieve only `S=1/2`. More generally `m` equal independent copies achieve `1/m`. Replication is therefore penalized directly by the convex geometry of the critical norm, before using the physical energy budget.

A stability argument shows that if `S\ge1-\delta`, a common component must carry nearly all mass on all three sides. This produces a quantitative notion of **stickiness from near equality**.

## 8. Transfer-weighted collision entropy and the Bellman cocycle

At a parent interaction component `v`, the correct weight is generated by its transfer contribution,

\[
\lambda_v=
\frac{(X_vY_vZ_v)^{2/3}}{
\sum_u(X_uY_uZ_u)^{2/3}}.
\]

If `c_X(v),c_Y(v),c_Z(v)` are child collision factors and `\rho_v` is a local reuse loss, two Hölder steps give an exact refinement estimate. In logarithmic form,

\[
\boxed{
-\log R
\ge
\Delta\mathcal H+\mathcal R,
}
\]

where `ΔH` is transfer-weighted conditional collision entropy and `\mathcal R=-\log\bar\rho` is the averaged reuse cost. Cross-interaction errors enter additively after logarithm and are harmless if summable.

This is the first true multiscale cocycle of the programme: **branch entropy and reuse loss add across scales**.

## 9. Nested Gaussian grain extraction with summable cross-error

For a finite equal-width Gaussian atomic model, each triad edge has

\[
w_e=\beta_e e^{-D(e)^2}
\]

with an exact defect `D`. A fixed phase-space grid is avoided. Instead, at each parent node, an annular **moat is chosen in defect space by transfer pigeonholing**. Deleting the cheapest annulus and connecting only shorter edges gives child components.

The node-level cross loss obeys

\[
\boxed{
\operatorname{Cross}(P)
\le
\left(M^{-1}+\varepsilon_{\rm tail}\right)T_P.
}
\]

With schedules such as `M_j=(j+3)^2` and summable tail tolerances,

\[
\boxed{\sum_j\operatorname{Cross}_j<\infty.}
\]

Children are recursively constructed inside parents and are never rejoined, so the decomposition is an actual nested grain tree.

### Percolation correction

A single connected component need not be geometrically small: many grains can chain through short interactions. This led to the exact incidence-graph identity

\[
\boxed{(n-1)+\beta=2m,}
\]

for a connected 3-uniform interaction component, where `n` is the number of packet vertices, `m` the number of triads, and `β` the incidence cycle rank. A large connected component is therefore fresh-rich or cycle-rich; chaining is not a neutral escape.

## 10. Hodge formulation of reuse curvature

Replace each triad by two directed parent-to-child arcs. Let `D` be the graph incidence matrix, `W` positive conductances, and `a` the desired scale-increment cochain. Define

\[
\boxed{
\mathcal E_H(a)=
\inf_\phi\|W^{1/2}(a-D^T\phi)\|_2^2.
}
\]

This is the squared distance of the desired increment field from exact gradients. Equivalently it has the cycle-space dual

\[
\mathcal E_H
=
\sup_{z\in\ker D,\,z\ne0}
\frac{(z\cdot a)^2}{z^TW^{-1}z}.
\]

The original nonflat reuse motif has strictly positive Hodge energy; in the normalized toy model it is approximately `0.0485733810`.

### Flat-cycle countermodel

A three-triad butterfly

- `a+b -> m`,
- `a+c -> n`,
- `m+n -> d`,

has cycle rank one but **zero scale Hodge energy**. Exact geometry gives `d=a`. Thus the statement “cycle rank implies scale holonomy” is false.

The countermodel revealed the missing branch: reuse cycles are either **curved** (pay Hodge energy) or **flat** (enter a rigid midpoint geometry).

## 11. Spherical erosion of flat networks

For unit directions `p,q` with spherical midpoint

\[
m=\frac{p+q}{|p+q|},
\]

and an open-hemisphere pole `n`, define

\[
\Phi_n(x)=-\log(n\cdot x).
\]

The exact midpoint barrier is

\[
\boxed{
\Phi_n(m)
\le
\frac{\Phi_n(p)+\Phi_n(q)}2
-\kappa(\theta),
\qquad
\kappa(\theta)=-\log\cos(\theta/2).
}
\]

At the optimal angle `θ_*≈70.1383°`,

\[
\kappa_*\approx0.2003318956.
\]

Hence a flat lineage contained in a fixed open hemisphere consumes a positive spherical potential every generation.

If a direction set escapes every open hemisphere, then `0` belongs to its convex hull. Carathéodory reduces this to at most four directions, and such a certificate must contain a pair separated by at least

\[
\boxed{
\theta_{\rm tet}=\arccos(-1/3)\approx109.4712^\circ.
}
\]

The escape is therefore a broad balanced configuration, not a narrow grain.

### Balanced states pay entropy

For an atomic direction measure with barycenter `b`,

\[
\boxed{
H_2(\mu)\ge\log\frac2{1+|b|}.
}
\]

Thus a nearly balanced state has definite collision entropy. In the equal-marginal flat propagation model,

\[
b_{j+1}=b_j/\cos(\theta_*/2),
\]

so a long-lived balanced chain pays asymptotically `log 2` collision entropy per generation up to a bounded correction.

For nonsymmetric parent/companion marginals, the exact local identity

\[
\boxed{
b(\nu_j)=2c_*b(\lambda_j)-b(\mu_j),
\qquad c_*=\cos(\theta_*/2),
}
\]

shows that the companion either has definite entropy or is strongly concentrated and therefore becomes a trackable fresh/reused grain.

## 12. Atomic entropy cannot hide inside ancestry components

Let `w_i` be atomic weights and `A(i)` ancestry labels. Define

\[
Q_{\rm at}=\sum_iw_i^2,
\qquad
Q_{\rm anc}=\sum_A\left(\sum_{i\in A}w_i\right)^2.
\]

Then

\[
\boxed{
Q_{\rm anc}-Q_{\rm at}
=
\mathbb P\{I\ne J,\ A(I)=A(J)\}.
}
\]

So entropy lost under ancestry coarse-graining reappears exactly as same-ancestry pair mass.

With the pair-biased ancestry law `\alpha_A\propto W_A^2`, if

\[
d=H_{\rm at}-H_{\rm anc},
\]

then

\[
\mathbb E_\alpha q_A=e^{-d}.
\]

Consequently, for every `λ>1`, at least `1-1/λ` of pair-biased ancestry mass lies in components containing at least

\[
\boxed{k_A\ge e^d/\lambda}
\]

reused atoms. Hidden entropy therefore creates multiplicity of reused attachments, hence a large ancestry cycle space after contraction.

## 13. Multicommodity Hodge routing and gauge synchronization

For a family of cycle commodities `(z_r,\mu_r)`, the basis-free Rayleigh inequality is

\[
\boxed{
\mathcal E_H
\ge
\frac{\sum_r\mu_r(z_r\cdot a)^2}
{\sum_r\mu_r z_r^TW^{-1}z_r}.
}
\]

For a reused terminal pair, route one unit electrically through the old ancestry network and one through the new interaction network. Their difference is a canonical cycle flow. The denominator becomes a pair-effective-resistance budget `\mathfrak R`.

If both old and new networks are individually flat, let `d_i` be their gauge difference on reused terminals. Then

\[
\boxed{
\operatorname{Var}_p(d)
\le
\frac12\mathcal E_H\mathfrak R.
}
\]

If scale gauges live in `\gamma_*\mathbb Z`, then small `\mathcal E_H\mathfrak R` forces most reused pair mass into a single integer gauge class. Thus low-curvature, low-resistance reuse synchronizes into the spherical flat-network module.

## 14. Resistance is not an escape: Poisson Bellman stopping

On a weighted tree with normalized conductances `c_e`, resistances `r_e=1/c_e`, terminal law `p`, and tree resistance distance `R(i,j)`, cut edges independently with

\[
q_e=1-e^{-r_e/\Lambda}.
\]

Then the probability a pair remains connected is exactly

\[
\boxed{
\mathbb P(i\leftrightarrow j)=e^{-R(i,j)/\Lambda}.
}
\]

Hence the expected collision probability of the random partition is exactly a Laplace transform of the pair-resistance distribution. A deterministic stopping argument yields a partition with positive component entropy and controlled deleted conductance.

The soft-complexity quantity

\[
\rho_\Lambda=\sum_e\min(c_e,\Lambda^{-1})
\]

is itself not an escape: if it is large, the conductance distribution has the exact lower bound

\[
\boxed{
H_2(c)
\ge
-\log(1-\rho_\Lambda+\Lambda^{-1}).
}
\]

Thus high resistance gives a cheap Bellman cut or atomic conductance entropy; low resistance feeds the Hodge synchronization theorem. A quantile version gives the finite-dimensional congestion trichotomy

\[
\boxed{
\text{Hodge cost}
\quad\text{or}\quad
\text{gauge synchronization}
\quad\text{or}\quad
\text{Bellman/atomic entropy},
}
\]

with no logarithm of the packet count.

## 15. Master finite-dimensional no-escape theorem

The cleanest final reset uses the **barycentric potential**

\[
\boxed{P(\mu)=-\log|b(\mu)|.}
\]

Fix `0<\beta<1`. If `|b|\le\beta`, the block pays the definite entropy cost

\[
H_2(\mu)\ge h_\beta:=\log\frac2{1+\beta}.
\]

If `|b|>\beta`, then the potential is automatically reset into

\[
0\le P<P_{\max}:=-\log\beta.
\]

During a synchronized flat episode,

\[
P_{j+1}\le P_j-\kappa_0+\zeta_j.
\]

Every costly block may start a new flat episode, so episode counting rather than a single global potential is essential. If `N_K` is the number of costly blocks, `N_F` the number of synchronized-flat low-cost blocks, and `Z=\sum\zeta_j`, then

\[
\boxed{
N_F\kappa_0
\le
(N_K+1)P_{\max}+Z.
}
\]

If every costly block pays at least `c_0>0`, and logarithmic cross-error penalty sums to `\Xi`, then a depth-`L` cascade satisfies

\[
\boxed{
-\log\prod_{j<L}R_j
\ge
\frac{c_0\kappa_0}{\kappa_0+P_{\max}}L
-
\frac{c_0(P_{\max}+Z)}{\kappa_0+P_{\max}}
-\Xi.
}
\]

Therefore, if perturbation and cross-error ledgers are summable, efficiency decays exponentially in depth.

A particularly natural reset is `\beta=c_*=\cos(\theta_*/2)`, for which `P_{\max}=\kappa_*`. In the exact flat model this gives an effective rate `c_0/2`: asymptotically at least about half the generations must pay a positive cost.

This theorem is the current **finite-dimensional no-escape closure**. It is conditional only through the constants and PDE-to-grain hypotheses supplied by the remaining bridge.

## 16. What is theorem-level, what is numerical, what is still conditional

### Exact / proved within the stated finite-dimensional models

- helical triad coupling formula under the repository convention;
- Bellman/Hölder component inequality and collision-entropy refinements;
- Gaussian equal-width overlap/resonance formulas;
- exact scale-holonomy identity;
- nested defect-moat cross-loss certificate and summability schedule;
- incidence fresh/cycle identity;
- Hodge primal/dual identities;
- flat butterfly certificate and `d=a` rigidity;
- spherical midpoint barrier, barycenter identities, tetrahedral threshold;
- atomic-to-ancestry collision chain rule and cycle-rank gain after contraction;
- multicommodity Hodge Rayleigh inequality and gauge synchronization;
- tree resistance identities and Poisson stopping theorem;
- barycentric episode-counting master no-escape inequality;
- exact helicity-sign reduction for the full normalized single-edge envelope;
- Arb-certified global uniqueness/maximality of the single-edge optimizer;
- Arb-certified mixed single-edge stability
  `Def_e >= (1/50)u+v^2`, the local edge-to-Hodge constant `c_stab=1/2`, and
  the global outside-neighborhood gap `Def_e>=1/100`;
- exact sharp-cutoff single-triad Mellin flux identity
  `int Pi_K dK/K = -dot E_k log(p/k)+dot E_q log(q/p)` and the Arb-certified
  `>=9/10` retention of the upper progress segment on the local stability box;
- Arb-certified rational countermodel showing the full all-scale Mellin moment
  exceeds `3J*/2` and therefore is not the correct edge extremal functional;
- exact graded-filter log convolution and filter-invariant Mellin moment;
- exact weak scale-potential identity and transfer-weighted common-midgap
  smooth-tail theorem recovering `sum T_e log(q_e/p_e)` with no smoothing loss;
- Arb-certified common smooth moat for shell halfwidth `2/25`, filter halfwidth
  `1/20`, with residual margin `>9/250`;
- exact flux polarization identity separating multiplier deficit from
  phase/backscatter deficit;
- Arb-certified good-core change of measure: at `eta_0=1e-4`, gap deviation
  `<=1/80` and physical child-transfer/capacity weight condition `<53/50`;
- physical-transfer Hodge cost
  `min{1/20000,(25/106) h_H}` on a packet block with the common moat;
- exact equality between global graded spectral flux and the space-average
  smooth SGS transfer `-grad ubar : tau(u,u)`; Leray/pressure does no global work.
- exact localized combined-work identity
  `G=Pi+div(P U)` in which pressure cancels from the window leakage ledger;
- exact nonincreasing-weight finite-chain depletion inequality for positive
  combined work, dissipation, positive leakage and backscatter;
- exact pressure-cancellation dichotomy: a positive raw SGS event either retains
  half its size as combined work or forces the scale-critical annular charge
  `int_A(|U|^3+|P|^(3/2)) >= S/(2||grad chi||_inf)`.
- Arb-certified four-bin crossing extraction: every signed-good positive-transfer
  core crossing one reference cut has a transfer subcore and a Hodge-numerator
  subcore of fraction at least `1/4`, shell halfwidth `<2/25`, and smooth moat
  margin `>1/25`; conservative physical Hodge coefficient `25/424`;
- exact one-shot profile replacement bound
  `3 eps_G+3 eps_G^2+eps_G^3` and exact finite-shell bridge
  `N||G||_2^2 >= C_Omega^(-1/3)(1-eps_G)^2`; Gaussian existence here uses
  Christ's published inverse Young theorem rather than a repository numerical
  certificate;
- exact pressure-kernel constant-mode cancellation and dipole decay
  `|grad K(z)|=O(|z|^-4)`, giving the three-dimensional `4-3=1` summable far
  shell gain; exact local pressure-to-critical-mass Bernstein/CZ reduction and
  the resulting positive fresh-mass threshold inside the stated packet model.
- Arb-certified narrow-shell critical-mass bridge: on the log shell
  `|log(|xi|/N)|<=2/25`, a Gaussian profile at `1%` `L^{3/2}` distance carries
  `N||G||_2^2 > 3/4` (certified enclosure `0.7740577380...`); the existence of
  that `1%` profile still uses Christ's external inverse-Young modulus;
- exact smooth-symbol freezing estimate
  `|T_M-T_Mh| <= A_3 L_M h prod ||f_j||_(3/2)` on each certified compact SGS
  crossing block, so relative frequency-cell errors can be chosen summable
  across generations without a packet-count loss.

- Arb-certified affine/Kelvin extremal strain rigidity in the invariant
  triad-plane model: planar trace-free strain satisfies
  `Hodge_speed^2 >= (43/100)||D||_F^2`; frozen `+/-d` strain with
  `dT<=1/25` pays the time-averaged edge deficit `>= (1/10)(dT)^2`;
- exact affine Gaussian envelope laws
  `Pdot=AP+PA^T+2nu I`, `kappadot=-A^T kappa-2nu P^-1 kappa`, the dual-center
  identity `(P kappa)dot=A(P kappa)`, incompressible log-det law, and Kelvin
  polarization/energy identities; exact affine-deforming window cancellation
  with curvature as the first Taylor remainder;
- exact intrinsic 3D Gram reduction `Gdot=-K^T(B+B^T)K`: scalar triad side-length/multiplier geometry depends only on the evolving `2x2` Gram matrix, so arbitrary common plane tilt is a gauge; the same `43/100` extremal coercivity applies to the trace-free symmetric restriction of the carrier driver on the instantaneous plane;
- Arb-certified five-percent co-rotating planar strain-coherence theorem:
  if `sup||D(t)-D(0)||<=d/20` and `dT<=1/30`, then the averaged edge deficit is
  `>= (1/24)(dT)^2`; failure forces objective strain variation `>d/20`;
- exact Navier--Stokes velocity-gradient/objective-strain identities
  `D_t A+A^2=-Hess p+nu Delta A` and
  `S_circ=-S^2-Omega^2+[S,Omega]-Hess p+nu Delta S`;
- exact curvature-balanced moving-moat law `E(M)=a/M+b kappa M`, optimizer
  `M*=sqrt(a/(b kappa))`, optimum `2 sqrt(a b kappa)`, together with the
  countermodel showing the fixed spatial schedule `M_j=(j+3)^2` can have
  summable commutator error but nonsummable curvature error;
- exact band-limited curvature-to-critical-mass implication and the objective-
  strain source collision in the stated packet model; far pressure Hessian has
  kernel order five, so 3D packet packing leaves the summable exponent `5-3=2`.

- exact spin-one helical bundle formulas `A_s=s cos(theta)dphi`,
  `F_s=-s sin(theta)dtheta wedge dphi`, and Chern number `-2s`; exact
  `SO(3)`-covariant triad-normal gauge and spin-dihedral transition
  `h_s(k;n2)=exp(-is psi)h_s(k;n1)`; single-triad normal-gauge coupling is
  quadrature, so rigid rotation is not a phase deficit;
- exact diamond modal-phase cancellation identity and sharp four-edge phase
  holonomy cost `sum(1-cos delta_i)>=4(1-cos(|H_phi|/4))`; Arb-certified clean
  finite-packet branch `|H_phi|>=1/5`, edge weights `>=beta_phi`, good
  multipliers `m_e>=1-10^-4 => D_phase>=beta_phi/250`;
- exact objective transverse polarization generator `-sym(A_perp)-nu|k|^2I`,
  exact off-diagonal helical representation of trace-free strain, and exact
  strain-area commutator law; a second-Magnus bound remains useful for an
  individual spinor, but the transfer-facing common-parent time ordering is
  closed exactly by the symplectic wedge identity below;
- Arb-certified full 3D strain tomography at the symmetric extremizer:
  `||D_Pi||_F^2+||D_1||_F^2+||D_2||_F^2 >= (13/20)||S||_F^2`, hence the combined
  scalar-shape/helicity-mixing observable is `>=559/2000 ||S||_F^2`;
- exact isosceles extremal helicity-tensor factorization
  `C_(s1,s2,s3)=-i C(r) s3 epsilon_(s1,s2)`, common-parent symplectic identity
  `(M U)^T J(M V)=U^TJV` for `M in SL(2)`, and relative formula through
  `M_1^-1 M_2`; Arb-certified transfer-distinguishable strain observable
  `||D_Pi||^2+||D_1-D_2||^2+||D_child||^2 >= (1/2)||S||^2`.  The `13/20`
  result is absolute tomography, not by itself a transfer-deficit theorem;
- Arb-certified signed-good angular separation
  `1/4<cos(theta)<2/5`, `sin(theta)>9/10`, exact triad-normal derivative bound
  `||dn||<=(10/9)(||da||+||db||)`, and helical-frame Lipschitz bound `5/2`;
  thus the global Chern obstruction does not create a local gauge singularity on
  a near-extremal packet core;

- exact arbitrary-time relative-polarization transport:
  `d(U^TJV)/dt=U^TJ(D1-D2)V` and the corresponding child-factor identity; a
  common time-ordered parent `SL(2)` history cancels pointwise, so no full Magnus
  theorem is needed for the transfer-facing parent observable; exact forced
  extension with additive `F_1,F_2,F_3` residual and capacity-weighted bound;
- Arb-certified low-strain localized polarization bridge: Kelvin direction
  stability plus good-core frame regularity gives
  `E_pol<=3h+(15/2)c kappa M` for `c sigma_0<=1/30`, hence the combined
  localization ledger is `a/M+(b+15c/2)kappa M+3h` with the same square-root
  curvature optimizer; no independent polarization moat scale is introduced;
- exact affine-Gaussian non-affine forcing decomposition: with
  `B=L^-1(nabla^2 U)[L,L]`, quadratic wavefront curvature is Gaussian tangent
  and the first transverse scalar forcing is third Hermite chaos
  `||F_perp||_2^2/||psi||_2^2=(3/8)||Sym B||_F^2`; bulk viscosity is also
  tangent to the complex Gaussian manifold;
- exact incompressible quadratic-swirl kernel `Sym B=0 <=> V(z)=z cross(Mz)`
  with `M` symmetric trace free; this scalar-forcing kernel is routed to spatial
  physical-strain / relative-polarization variation, where
  `E_z Q_rel(S(z)-S(0)) >= (1/2)||C||_F^2` for
  `C_ijc=sym_ij((nabla^2U)_ijk L_kc)`;
- exact material triad phase-lock gauge: common nonlinear resolved advection
  preserves `phi_1+phi_2-phi_3`, carrier resonance and common chirp/Hessian lock;
  with role-dependent velocities only differential velocities/sources remain;
- Arb-certified affine shell/aspect bridge and actual ellipsoidal mass:
  every physical Gaussian standard axis is `>2/(3N)` and the actual role on the
  radius-two covariance ellipsoid satisfies
  `integral_E |u|^2 >= (3/10)(det Sigma_x)^(1/6)`;
- exact affine critical-grain reformulation
  `M_aff=r_g^-1 integral_E|u|^2>=3/10`, `r_g=(det Sigma_x)^(1/6)`, and fresh
  energy budget `sum r_g<=P||u(t)||_2^2/eta`; under incompressible affine
  Gaussian dynamics `d log r_g/dt=(nu/6)tr Sigma_x^-1` and
  `d(r_g^2)/dt>=nu`, so inviscid strain preserves this radius exactly; common
  affine anisotropy is an exact Young symmetry and is explicitly recorded as a
  countermodel to charging aspect ratio itself as Bellman/replication cost;
- exact ellipsoidal moving-window geometry
  `D_t z=L^-1[U(X+Lz)-U(X)-A(X)Lz]`; shell-axis regularity gives
  `N^-1||grad chi_(L,M)||<=3 C_chi/(2M)`, Taylor leakage is
  `<= (C_chi R_chi^2/2) kappa_aff M`, and the physical filter commutator obeys
  `||[chi,G_N*]f||_2 <= (3/2)m_1(G)C_chi M^-1||f||_2`; hence the affine window
  retains the square-root balance `a/M+b kappa_aff M` without an aspect penalty;

- Arb-certified strict low-pass transport separation on the signed-good core:
  every selected role has frequency `>3N/5`, while a transporter supported in
  `|xi|<=N/4` has low-low output in `|xi|<=N/2`; exact affine low-high Fourier
  dynamics is the Kelvin equation and objective 2x2 polarization generator;
- exact divergence-free weak packet coefficient equation using a shell-localized
  Leray/helical/cell multiplier; microscopic pressure vanishes and the
  window/multiplier commutator is `O(1/M)` with the affine shell constant;
- exact Hermite/helicity forcing split: scalar Gaussian tangent modes have degree
  `<=2`, the `Sym B` residual is an orthogonal `H_3` sideband with norm square
  `(3/8)||Sym B||^2`, and spatial polarization curvature is an `H_1` sideband
  with zero base-Gaussian projection; only degree-zero forcing enters the forced
  symplectic spinor identity;
- exact work-level identification of parent/profile degree-zero forcing with the
  existing trilinear cross-error: `(2 eps+3 eps^2+eps^3)+eps` equals the full
  one-shot replacement polynomial `3 eps+3 eps^2+eps^3`;
- exact affine resolved-boundary classification: differential advection plus
  `RU` forces a cubic critical annular charge, resolved viscous boundary flux is
  `O(1/M)` over a parabolic lifetime, a complete quadratic spatial partition has
  exact total boundary cancellation, and `R` has the exact velocity-increment
  representation with the cubic increment upper bound recorded above;

- exact affine curvature connection
  `dot B+2 A_aff B=L^-1 Hess(D_t^V V)[L,L]`; for a resolved Navier--Stokes
  transporter this becomes the source split
  `-nabla^3 P-nabla^2 div R+nu nabla^2 Delta V`; pressure-third far-field
  homogeneity leaves the 3D packing exponent `6-3=3`;

- exact Banach interval-variation / first-Duhamel sideband theorem
  `int||f|| <= ||int f||+(T/2)int||f_dot||`; after pulling back the affine
  curvature connection, `H_3` dephasing is sourced by the acceleration-Hessian
  term rather than common affine motion; Gaussian hypercontractivity plus
  certified arithmetic gives critical Young footprints `>1/160` for `H_3` and
  `>1/16` for a two-component `H_1` sideband;
- exact odd-Hermite Gaussian-triad parity rule: a nonzero centered resonant
  interaction has even total Hermite parity, so one `H_1` or `H_3` daughter has
  zero first-order base-triad matrix element; rescue starts quadratically, and
  for one odd degree-`<=3` sideband with `sigma<=1/80`, the normalized transfer
  deficit is `>=sigma^2/16` unless another odd sideband participates;
- exact differentiated-incompressible curvature irreducible split
  `15=(7+3)_envelope+5_swirl`, with orthogonal reconstruction `B=B_E+B_H` and
  intrinsic observability `||Sym B||^2+||C_hook||^2 >= (1/6)||B||^2`;
- exact local `H_3` no-escape theorem: for
  `I_3=int||Sym B_tilde||` either source variation is `>=I_3/T`, nonlinear
  sideband feedback is large, daughter capacity reaches `1/80`, or transfer
  deficit / pair-sideband rescue is at least `(3/4096) I_3^2`;
- exact odd-sideband pair-rescue ancestry graph: cross-ancestry rescue enters
  `Xi`; same-ancestry endpoint weights obey the atomic/component collision chain
  rule and route into dominant reuse, Bellman component entropy, or quantitative
  same-ancestry pair/cycle mass without assuming spatially disjoint Hermite
  packets;

- Arb-certified physical hook/swirl polarization bridge on the
  mild-aspect branch: with no scalar `D_Pi`, `Q_pol(I)>=(1/10)||B_hook||^2`
  and `cond(L)<=21/20` gives `Q_pol(L)>=(1/25)||B_hook||^2`; the auxiliary
  relative-coordinate forcing energy is `>=1/50||B_hook||^2`, while the sum of
  the **three physical Young-role sideband energies** is
  `>=1/100||B_hook||^2`; high aspect is routed to affine ancestry rather than
  charged as a defect;
- Arb-certified physical mild-aspect `H_1` no-escape theorem on the low-strain
  action branch `K<=1/30`: nonunitary pullback/pushforward gives
  `J_1>=I_1/(11T)` or a physical three-role daughter
  `delta_1^2>=I_1^2/480`; after role splitting and odd-Hermite rescue,
  transfer deficit / pair rescue is `>=I_1^2/184320`; with the `H_3` theorem,
  full mild-aspect curvature has clean cost `I_B^2/737280` outside the named
  source/feedback/daughter/rescue branches;

- Arb-certified physical `H_1` covariant source calculus:
  `J_1<=2 int||S_src||+54 int||A||||B||` on the original mild branch; H1
  dephasing forces a physical curvature source or velocity-gradient/frame
  coupling, and below the stated base-action threshold one of pressure-third,
  differentiated-SGS, or viscous-fourth channels has a quantitative source
  integral;
- Arb-certified filtered differentiated-SGS/viscous collision: because
  `V=S_(N/4)u` and `supp Rhat subset B_(N/2)`, `N^-4||S_R||<=(3/2000)s||R||_(3/2)`
  and `N^-4||S_nu||<=(nu s/5000)sqrt(d_V)`; large radius is routed to affine
  critical mass instead of aspect cost;
- exact dyadic Onsager collision given standard `L^3` Littlewood--Paley and
  Bernstein: cubic increment charge satisfies
  `Q_N<=g1(C_LP C_B)^3[(4/3)(beta R_G)^2 mu_low_max+2 d_high]^(3/2)`, hence a
  low/base critical-mass reservoir or high-frequency normalized enstrophy;
  packet fragmentation yields a dominant `1/4` atom, ancestry entropy `>=log 2`,
  or same-ancestry pair/cycle mass `>=1/4`;
- exact source-weighted SGS/viscous episode routing: the `3/2` and `2/3`
  homogeneities cancel, so SGS source weight converts linearly to mass/enstrophy
  without a temporal-persistence assumption; viscous source concentration only
  increases the quadratic dissipation cost;
- Arb-certified no-gap aspect transition: `cond(L)<=567/500` retains
  `Q_pol>=1/4000||B_hook||^2` and a positive physical H1/full-curvature cost;
  above `567/500`, low symmetric-strain action implies the predecessor was
  already above `21/20` unless the grain is fresh/relinked.  Fresh aspect obeys
  `N int_E|u|^2>(1/5)cond(L)^(1/3)` through physical radius, not a Young deficit;

### Computationally supported, not interval-certified

- numerical reuse-gap constants from nonlinear optimization;
- perturbative robustness tables for near-butterfly / finite-width Gaussian models.

### Still conditional / PDE bridge

1. **Classify every physical service-or-flat exit into one master currency.**
   The multi-currency telescope distinguishes multiplicative transfer costs from
   globally bounded additive resources.  The remaining PDE audit must assign each
   source/fresh/dissipation/initial-boundary exit exactly once.  Critical `N E` or
   normalized `D_V` branches must remain in sticky/causal routing unless a genuinely
   scale-independent global resource threshold is separately proved; no raw source
   norm may be counted as a finite resource without its collision theorem.
2. **Finish the one global selected-interface `Xi` ledger.**  Time-slab synchronization
   costs zero, nested material-label refinement costs zero, and frequency/covariance
   representative schedules are summable.  The canonical global route now uses
   Leray-projected coherent Parseval probes, so compact moving-window commutators and
   localized pressure boundary work are absent from this master ledger.  What remains
   is actual selected cross interaction/profile-defect rejection and any continuous-to-
   discrete synthesis approximation, each charged once.
3. **Run the iterative PDE extraction through the full physical depth.**  Every
   efficient Navier--Stokes block must enter the certified service-or-flat gate with
   the same physical transfer law and the same canonical material labels, so the
   multi-currency master telescope applies to the actual recursive extraction.
No statement in this repository currently closes these PDE gaps, and no claim of Navier--Stokes global regularity is made.

## 17. Reproducible provenance

The most useful recorded runs, in chronological order, are:

| Run | Module / finding |
|---:|---|
| `31142036250` | single-edge scale optimization and numerical reuse gap |
| `31142312810` | fast diamond probe |
| `31142882572` | Gaussian packet inverse experiment |
| `31143091088` | scale-critical weighted Young / cusp experiment |
| `31143784774` | corrected transfer-preserving grain profiles |
| `31144432129` | multiscale Bellman cocycle |
| `31145327819` | nested Gaussian grain extraction + percolation diagnostics |
| `31146490082` | Hodge cycles, flat butterfly, planar erosion |
| `31147082764` | spherical erosion and barycentric entropy bridge |
| `31150171727` | atomic-to-component entropy transfer |
| `31151008574` | multicommodity Hodge routing / gauge synchronization |
| `31152386368` | no-log resistance-to-Bellman stopping |
| `31153769553` | corrected cap-reset master theorem regression run |
| `31154025683` | barycentric master no-escape theorem; `72` tests + `20,000` episode traces |
| `31157463384` | first Arb-certified full single-edge sign/global/local stability; `76` tests + `200,000` adversarial samples |
| `31159084424` | exact sharp-cutoff log-scale flux bridge; `84` tests + local retention grid |
| `31159086953` | updated Arb single-edge certificate including `>=9/10` Mellin-flux retention; `84` tests + `200,000` adversarial samples |
| `31160779428` | full all-scale Mellin adversarial search exposing a different asymmetric maximizer |
| `31161316034` | smooth log-scale cocycle / transfer-weighted midgap theorem; `93` tests + `50,000` blocks |
| `31161626056` | dealiased SGS-to-graded-spectral flux regression; `94` tests + `200` fields |
| `31161811900` | master regression after physical-transfer Hodge insertion; `94` tests + `20,000` episode traces |
| `31161914134` | preferred Arb certificate: smooth moat, weight change, and rational full-Mellin countermodel; `94` tests |
| `31162700474` | localized SGS/pressure-work depletion and critical-annulus trichotomy; `97` tests + `50,000` chain traces |
| `31164771160` | Arb-certified four-bin crossing-to-common-moat extraction; `99` tests + `50,000` crossing blocks |
| `31164911526` | pressure multipole `4-3` packet collision theorem; `104` tests + `50,000` shell configurations |
| `31165172325` | one-shot transfer-preserving Gaussian profile algebra; `108` tests + `100,000` parameter checks |
| `31165654509` | preferred one-shot profile certificate: Arb narrow-shell critical mass `>3/4`; `108` tests + `100,000` checks |
| `31165838379` | smooth-SGS symbol freezing / summable relative-cell error; `111` tests + `50,000` cell checks |
| `31166152074`--`31166171000` | final eight-workflow bridge integration on `37e1380`: single-edge, crossing moat, profile mass, symbol freezing, annular pressure, localized pressure, nested grains and master all green; `111` tests per workflow |

| `31168205564` | affine-grain certificate failed at an exact local-box boundary because Arb cannot certify touching-ball equality; no mathematical countermodel |
| `31168303213` | second affine-grain certificate failed at the exact rational average-cost equality `1/10` for the same interval-comparison reason |
| `31168393589` | preferred affine Gaussian/Kelvin grain dynamics certificate; `118` tests + `50,000` affine/viscous/Kelvin checks |
| `31168671304` | five-percent objective strain-coherence certificate; `121` tests + `20,000` variable-strain traces |
| `31168888413` | curvature-balanced moving moat and fixed-schedule countermodel; `127` tests + `100,000` balance checks |
| `31169128097` | objective-strain NS source collision and pressure-Hessian `5-3=2` locality; `133` tests + `50,000` packet checks |

| `31169598408`--`31169621311` | integrated affine-spacetime bridge on `2b795d2`: affine grain, strain coherence, curvature balance, objective-strain collision, localized pressure, crossing moat, nested grains and master all green; `133` tests per workflow |

| `31170015795` | intrinsic full-3D triad-plane Gram dynamics; `137` tests + `50,000` random plane/driver/Gaussian checks |

| `31171018560` | helical spin/Berry transport: triad-normal `SO(3)` gauge, Chern charge `-2s`, spin-dihedral transition, and strain/helicity mixing; `144` tests + `50,000` checks |
| `31171127537` | Arb-certified diamond phase holonomy branch `D_phase>=beta_phi/250`; `148` tests + `50,000` diamonds |
| `31171360107` | Arb-certified full 3D strain observability `>=13/20` and combined shape/helicity constant `559/2000`; `149` tests + `100,000` strains |
| `31171484247` | objective helical polarization / strain-area commutator; `152` tests + `50,000` checks |
| `31172607991` | extremal helicity tensor / common-parent `SL(2)` neutrality and Arb-certified transfer-relevant strain observable `>=1/2`; `156` tests + `50,000` checks |
| `31173321354` | Arb-certified good-core parent-angle separation and scale-free triad-normal helical-frame Lipschitz theorem; `158` tests + `50,000` derivative checks |

| `31174454743` | first exact relative-polarization transport / common-time-ordered `SL(2)` countermodel; `164` tests + `50,000` histories |
| `31174612248` | Arb-certified localized polarization packet bridge `E_pol<=3h+(15/2)c kappa M`; `170` tests + `50,000` configurations |
| `31174812731` | preferred forced relative-polarization transport: exact additive nonlinear forcing residual; `171` tests + `50,000` histories |

| `31175180687`--`31175212421` | integrated relative-polarization/localization stack on `e0a7855`: forced relative transport, localized polarization, symplectic quotient, smooth-symbol freezing, curvature balance, objective-strain collision, localized pressure and master all green; `171` tests per workflow, master `20,000` traces with worst margin `0` |
| `31179739773` | affine-covariant Gaussian forcing / exact third-Hermite transverse residual; `174` tests + `50,000` affine checks |
| `31179827015` | Arb-certified shell-constrained affine aspect and actual radius-two ellipsoidal mass `>=3/10 r_g`; `178` tests + `50,000` covariance checks |
| `31180257124` | exact five-dimensional incompressible quadratic swirl kernel; `180` tests + `50,000` checks |
| `31180506627` | exact common material triad phase/chirp gauge; `182` tests + `50,000` checks |
| `31180882083` | affine critical-grain energy/radius budget (before provenance-only cleanup); `186` tests + `50,000` checks |
| `31181076691` | affine spatial polarization-curvature RMS bridge `>=1/2`; `187` tests + `50,000` Hessian/grain checks |
| `31181950880`--`31182017032` | eleven-workflow affine-residual integration on `1816ee2`: affine forcing, shell/aspect, swirl, material phase, affine critical grain, polarization curvature, forced relative polarization, localized polarization, objective-strain source, localized pressure and master all green; `189` tests per workflow, master `20,000` traces with worst margin `0` |
| `31182421821` | ellipsoidal moving-window material/gradient curvature balance; `194` tests + `50,000` checks |
| `31182742846`, `31182749191`, `31182756039`, `31182762696`, `31182769409`, `31182775605`, `31182783321`, `31182790116`, `31182796814` | final nine-workflow affine-packet integration on `040b5c7`: affine Hermite forcing, material phase, shell/aspect, affine critical grain, affine moving window+commutator, polarization curvature, forced relative polarization, localized pressure and master all green; `195` tests per workflow, master `20,000` traces with worst margin `0` |
| `31183237575`, `31183245218`, `31183252753`, `31183260464`, `31183267939`, `31183274599` | exact-SHA affine validation on `be4520b`: Hermite/bulk-viscosity gauge, material/differential phase, affine critical grain including `r_g` dynamics, ellipsoidal window+commutator, polarization curvature and master all green; `195` tests per workflow, master `20,000` traces with worst margin `0` |

| `31186597005` | Arb-certified smooth-SGS two-level packet equation: `3/5` role separation, strict `N/4` transporter and zero affine moving-cell Heisenberg residual; `203` tests + `50,000` checks |
| `31186600162` | first affine-SGS boundary certificate failed only at exact Arb equality semantics for `3/2*11/10=33/20`; no mathematical countermodel |
| `31186719319` | corrected affine-SGS boundary ledger: critical cubic leakage, `O(1/M)` viscous boundary and exact partition cancellation; `203` tests + `50,000` checks |
| `31187104071` | Hermite-helicity forcing ledger: exact `H_1/H_3` orthogonality and no direct base-spinor curvature forcing; `206` tests + `20,000` checks |
| `31187443147` | base-spinor parent-remainder work identified with the existing trilinear cross-error / replacement loss; `208` tests + `50,000` checks |
| `31187637446` | exact affine low-high Navier--Stokes Fourier equation -> Kelvin carrier/amplitude -> objective 2x2 spinor; `211` tests + `50,000` checks |
| `31187768828` | preferred affine-SGS boundary continuation including exact SGS velocity-increment identity and cubic increment bound; `212` tests + `50,000` checks |
| `31188118452` | divergence-free localized weak packet equation; Leray pressure cancellation and shell multiplier/window commutator; `215` tests + `50,000` checks |

| `31188755726`--`31188831280` | twelve-workflow actual-packet integration on `88f5173`: smooth-SGS role separation, affine Kelvin PDE, divergence-free packet, Hermite sidebands, base-spinor cross forcing, affine SGS boundary/increments, one-shot profile, symbol freezing, affine window, relative polarization, localized pressure and master all green; `215` tests per workflow, master `20,000` traces with worst margin `0` |
| `31189036483` | exact affine curvature connection / resolved acceleration-Hessian source split and pressure-third `6-3=3` far locality; `217` tests + `50,000` checks |

| `31189330225`, `31189336088`, `31189342467`, `31189348870`, `31189354561`, `31189360561`, `31189366109`, `31189372330`, `31189378786` | exact-SHA actual-packet validation on `f27ed04`: role separation, affine Kelvin PDE, divergence-free weak packet, Hermite sidebands, base-spinor cross forcing, affine SGS boundary/increments, curvature connection, localized pressure and master all green; `217` tests per workflow, master `20,000` traces with worst margin `0` |
| `31192445017` | exact divergence-free curvature irreducible split `15=(7+3)+5` and intrinsic `H3/H1` observability `>=1/6`; `227` tests + `1` optional skip + `50,000` checks |
| `31192448279` | covariant first-Duhamel sideband coherence and critical Young footprint (`H3 >1/160`, two-component `H1 >1/16`); `228` tests + `20,000` checks |
| `31192451039` | odd-Hermite triad selection, quadratic rescue, and small-sideband transfer deficit `>=sigma^2/16`; `227` tests + `1` optional skip + `50,000` checks |
| `31192748440` | H3 no-escape stress failed because the synthetic sampler called `uniform(low,high)` with `low>high` in a narrow threshold interval; implementation failure, no theorem countermodel |
| `31193125969` | preferred local H3 no-escape theorem with clean cost `3 I_3^2/4096`; `229` tests + `1` optional skip + `50,000` five-branch checks |
| `31193443926` | pair-ancestry workflow failed because its hand-written test graph had cross fraction `1/3` while expecting the majority-cross branch; fixture failure, theorem routing was correct |
| `31193531201` | preferred odd-sideband pair-rescue ancestry graph: `Xi` / dominant reuse / Bellman entropy / same-ancestry cycles; `231` tests + `1` optional skip + `50,000` graphs |
| `31194159878`--`31194224545` | twelve-workflow Hermite-sideband integration on `e45eaea`: coherence, parity, curvature irrep, `H_3` no-escape, pair ancestry, curvature source, Hermite-helicity, physical polarization curvature, one-shot profile, relative polarization, localized pressure and master all green; master `231` tests + `1` optional skip + `20,000` traces with worst margin `0` |
| `31194787885` | first Arb-certified mild-aspect `H_1` polarization bridge (`Q_pol>=1/10`, mild `>=1/25`, auxiliary relative-coordinate energy `>=1/50`); valid pointwise provenance but predates the explicit three-physical-role distinction |
| `31195130386` | superseded idealized H1 no-escape algebra (`1/25600`, `1/102400`): treated the relative-coordinate pullback as isometric and did not split the daughter among physical Young roles; retained as correction provenance, not preferred theorem |
| `31196456570` | preferred physical mild-aspect H1 bridge on `4a5729c`: `Q_pol>=1/25`, auxiliary relative energy `>=1/50`, three-role physical sideband energy `>=1/100`; `237` tests + `50,000` configurations, Arb-certified |
| `31196458909` | preferred physical H1 no-escape on `4a5729c`: nonunitary low-strain conditioning, `J_1>=I_1/(11T)`, `delta_1^2>=I_1^2/480`, pair/deficit `>=I_1^2/184320`, full mild-curvature `>=I_B^2/737280`; `237` tests + `50,000` checks |

| `31171921187`--`31171950823` | integrated helical/spacetime stack on `6226fd9`: spin transport, explicit spin-dihedral phase holonomy, full-strain tomography, objective polarization, intrinsic 3D plane, affine grain, strain coherence and master all green; `153` tests per workflow, master `20,000` traces with worst margin `0` |

| `31198595416` | preferred physical H1 covariant-source calculus: `J_1<=2 int||S_src||+54 int||A||||B||`, source/coupling routing; `240` tests + `50,000` checks, exact projector/product-rule residuals at `1e-15` |
| `31226401666` | filtered differentiated-SGS / viscous collision: `N^-4||S_R||<=(3/2000)s||R||_(3/2)`, `N^-4||S_nu||<=(nu s/5000)sqrt(d_V)`; `246` tests + `50,000` checks, Arb-certified |
| `31226526751` | Onsager increment collision: cubic increments -> low/base critical mass or high-frequency enstrophy; packet split -> dominant atom / `log 2` ancestry entropy / `1/4` cycle mass; `246` tests + `5` optional skips + `50,000` checks |
| `31226868712` | source-weighted episode stress failed only because a large algebraic viscous identity was compared with an absolute floating tolerance; diagnostic failure, no theorem countermodel |
| `31227448616`, `31227450217` | unified mild/extended reruns failed because `source_divisor` was referenced in a function body before being added to its signature; implementation wiring only |
| `31227535026` | preferred unified source-weighted H1 episode on `2861f9a`: `C_src=132` mild, `1800` extended; SGS source weight -> mass/enstrophy/dissipation, viscous concentration raises dissipation; `264` tests + `50,000` checks |
| `31227536542` | preferred no-gap affine-aspect theorem on `2861f9a`: transition `567/500`, `Q_pol>=1/4000`, extended H1/full costs, condition-growth/inheritance law, fresh aspect-radius mass; `264` tests + `50,000` checks, Arb-certified |

| `31227852933`--`31227894352` | fourteen-workflow source/sideband/ancestry integration on `9c5f465`: H1 source, SGS/viscous collision, Onsager increment, source-weighted episode, sticky aspect, H1/H3, pair ancestry, curvature, affine critical grain, relative polarization, localized pressure and master all green; master `20,000` traces, worst margin `0` |
| `31228679822` | Arb-certified ancestor-reservoir synchronization on `d3f8a51`: signed-good progress `>8/5`, low-strain Kelvin growth `<21/20`, critical-mass and physical-energy service half-life, exact material covector `q=L^T k`; `268` tests + `50,000` checks |
| `31228947363` | preferred band-limited H1 source sampling on `c4c87f6`: separated SGS/pressure `ell^(3/2)` and viscous `ell^2` source charges, resolvable cluster counting; `267` tests + `7` optional skips + `50,000` checks |
| `31229118774` | pressure-reservoir synchronization on `6e21fe6`: one fixed low-low pressure pair has service ratio `194481/655360<1/3`; `270` tests + `7` optional skips + `50,000` checks |
| `31229543390` | whole-old-reservoir-pool erosion on `8c50370`: under frame budget `P`, old-pool service decays `<2^-q`; uniform service forces relinking, then exact incidence Euler routes connected relinking to fresh-rich or cycle-rich; `274` tests + `7` optional skips + `50,000` pool checks + `5,000` incidence checks |
| `31229863366` | preferred affine coherent-state Bessel/Riesz certificate on `84ae98c`: exact intrinsic overlap, 4-separated Bessel `25/4`, 5-separated Gram spectrum `[47/50,53/50]`, synthesis budget `P=50/47`; `286` tests + `50,000` affine checks + finite Gram probes |
| `31230201131` | exact affine coherent Moyal ledger on `975f680`: positive continuous phase-space cell energy with `P=1` by Plancherel/Fubini; `282` tests + `8` optional skips + `5,000` discrete Moyal regressions |
| `31230314319`--`31230338382` | twelve-workflow coherent-reservoir integration on `7b350ef`: physical scale progress, Moyal `P=1`, coherent Riesz `P=50/47`, whole-old-pool erosion, nested fresh/cycle graph, ancestor and pressure synchronization, band-limited source sampling, Onsager/source episodes, sticky aspect and master all green; master `20,000` traces, worst margin `0` |
| `31245808950`, `31245810803` | first coherent-work/increment reruns failed because three new increment-service fixtures supplied less total cell service than the theorem hypothesis `S_low>=Y-2d_high`; fixture failures, no mathematical countermodel |
| `31245848657` | exact polarized-Moyal coherent nonlinear-work cells and finite-time selected-cell no-escape; `293` tests + `8` optional skips + `20,000` checks |
| `31245850938` | actual SGS increment -> coherent phase-space service edges; clean `Y/4`, `Y/8`, `Y/32`, `log 2`, `1/4` routes; `293` tests + `8` optional skips + `50,000` checks |
| `31245962226` | first covariance-interface run failed because the scalar Gaussian-overlap fixture used a value greater than one; test normalization failure, matrix theorem unchanged |
| `31246026975` | preferred coherent covariance-interface theorem: exact Gaussian fidelity, log-covariance window bound and Moyal energy-TV `Xi_cov`; `297` tests + `8` optional skips + `50,000` SPD/window checks |
| `31246194276` | coherent-service stopping theorem: old-pool ratio `231525/512000<1/2` gives finite cost-free sticky epochs under a uniform service threshold; `301` tests + `8` optional skips + `50,000` checks |
| `31246242892` | source-weighted differentiated-SGS -> coherent ancestry with no persistence assumption; exact `3/2` to `2/3` linear homogeneity and coherent `Xi`/new-mass routing; `303` tests + `8` optional skips + `50,000` checks |
| `31246291317`--`31246322667` | fifteen-workflow coherent-service integration on `a40e2c1`: smooth physical flux, H1 source, SGS/Onsager/source episodes, Moyal work, coherent increment/covariance/stopping/source-weight, old-pool erosion, nested grains, band-limited sampling and master all green; master worst margin `0` |
| `31246845056` | preferred uniform service/action-or-Kelvin-flat packet assembly: `delta_tau=tau^2/1,036,800,000`, `F_K<=tau`; `316` tests + `50,000` branch checks, Arb-certified |
| `31247033855` | preferred physical flat-episode perturbation: signed-good barycenter error and clean master erosion `kappa_0>0.17` at `tau=1/100`; `319` tests + `50,000` checks |
| `31247206213` | first flat-companion run failed only at exact floating comparison of `7/9`; fixture semantics, no theorem countermodel |
| `31247281011` | preferred flat companion gate: `H_2>1/200` or two directional parent cores of mass `>=7/9`, angle `>1` rad and chord gap `>1/3`; `322` tests + `50,000` checks, Arb-certified |
| `31247418331`--`31247440588` | twelve-workflow service/flat/master integration on `b48e9ac`: service-or-flat, physical flat episode, flat companion, smooth physical flux, H1 source, coherent work/increment/source, old-pool erosion, nested grains, spherical erosion and master all green |
| `31247611869` | first causal-binary run failed because fixture `[4,3,2,1]` expected cycle rank `2`; exact reuse savings are `(2,1,0)` and rank `3`; fixture error only |
| `31247674599` | preferred synchronized causal binary ancestry: signed-good `3/5<parent/child<5/8`, exact cycle/product/reuse-action identities and quarter-reuse stopping; `326` tests + `50,000` layered graphs, Arb-certified |
| `31248342830` | transfer-weighted Shannon causal reuse information: exact `H(child,role|parent)` telescope and `log(4/3)` rich-layer threshold; `334` tests + `50,000` causal laws |
| `31248345125` | preferred Renyi causal reuse: exact binary-slot collision baseline, deep-layer `theta>1/3`, then ancestry entropy `>=0.5 log2`, same-ancestry pair `>1/5`, or parent-slot hidden pair `>1/12`; `334` tests + `50,000` causal laws |


| `31248783086`, `31248923266` | first two adjoint Kelvin--Duhamel workflow attempts failed only because the test fixture wrote `2+j` instead of the complex literal `2+1j`; `338` other tests passed, so these are implementation provenance, not mathematical countermodels |
| `31249473037` | preferred adjoint Kelvin--Duhamel causal gate on `3272ee4`: exact adjoint pairing removes common Kelvin/bulk-viscous transport from the source ledger; inherit/residual/high--high alternative `A/4,A/4,A/2`; positive aligned parent-pair generation; natural parent-window overlap `>103/50 T_child`; full pytest + `50,000` instantaneous/phase checks + `5,000` exact block-exponential histories, Duhamel residual `1.461e-15` |
| `31250408864` | preferred asynchronous Duhamel synchronization cone on `b9c17dc`: first parent-event span `<=25/128`, recurrence `alpha_next<=(25/64)(alpha+2/5)`, invariant cone `alpha<=3/8`, original loose common-slice certificate; recurrence/fixed point valid, but its cumulative `10/39` backward-horizon coefficient is superseded by the corrected asynchronous geometry and single-boundary coherent registration rule; `347` tests + `50,000` checks |
| `31250736899` | canonical material coherent labels on `f7be08c`: intrinsic `zeta=(L^-1 X/2,L^T k)` is common-affine invariant; nested dyadic refinement has zero Moyal switch cost; geometric frequency/covariance schedules give summable representation `Xi`; `352` tests + `50,000` checks, affine residual `1.877e-15` |
| `31250867463` | physical multi-currency master on `75d0f09`: conditional single-charge algebra `N_A<=sum B_r/b_r` for **uniform globally bounded** resource classes; flat erosion gives `N_T>=[kappa0 L-Pmax-Z]/[kappa0+Pmax]-N_A`; such true resource resets alter only the finite offset, not `c_eff`; `355` tests + `50,000` episodes; critical `N E`/`D_V` applicability is explicitly corrected by `31251119223` |
| `31251119223` | high-strain lifetime collision on `3d41d55`: Arb-certified `K>1/30 => D_V>32 pi^2/(75c)` for `V=S_(N/4)u`; also records the geometric-scale countermodel showing critical `N E` and `D_V` costs are `O(1/N)` physically and cannot be treated as uniform finite-count resources; `359` tests + `50,000` checks |
| `31251337503` | first resolved-objective-strain stress run failed because the diagnostic tried to saturate a deliberately weaker clean inverse bound (`1/380`) instead of the sharper analytic coefficient; theorem constants had already passed Arb, so this is diagnostic-direction provenance, not a countermodel |
| `31251421182` | preferred resolved objective-strain collision on `01a0915`: exact filtered corotational source identity for `A=grad V`; clean order-two constant `<1/380`, `||V||_3^2<mu_V/15`, pressure Hessian `rho_P<=mu_V/5700+||R||_(3/2)/380`, and viscous constant `<1/1500`; `362` tests + `50,000` checks, identity residual `6.874e-16` |
| `31251605854` | divergence-free coherent Parseval frame on `cfe34a5`: Leray-projected coherent probes preserve Moyal coefficients exactly on `L^2_sigma` and cancel pressure globally, removing compact moving-window/pressure-boundary interfaces from the canonical global ancestry route; `356` tests + `9` optional skips + `20,000` checks |
| `31252438256` | physical transfer-weighted defect moat on `f7d1d37`: on the `eta_0=10^-4` signed-good core, `E_phys[D]<=(106/25) epsilon`; deleting the least-transfer bin in `[R/2,R]` gives `eta_cross<=1/M+2 Dbar/R`; the schedule `M_j=M_0(j+2)^2`, `R_j=R_0(j+2)^2` makes the actual omitted cross-cell physical-transfer `Xi` summable with no packet-count or Gaussian-synthesis tail; `362` tests + `9` optional skips + `50,000` weighted triad graphs |
| `31258595663`, `31258597148`, `31258598281`, `31258599570`, `31258600704`, `31258601942`, `31258603563` | seven-workflow canonical cross-cell integration on `3d1381f`: exact coherent localization, physical transfer defect moat, canonical material labels, Renyi causal reuse, physical multi-currency master, smooth symbol freezing and master all green; master `20,000` traces, worst margin `0`; localization trilinear reconstruction residual `5.288e-14`; defect-moat minimum cross margin `6.917e-4` |
| `31261790706` | single-charge physical branch compiler on `f8cf98d`: exact positive-transfer partition `dT=dT_Xi+sum_r dT_r`, theorem manifestations quotiented by physical causal root, first-causal stopping on the retained measure, symmetric Radon--Nikodym splitting of exact independent ties, absorbing `t=0` boundary, and master projection which refuses critical `N E`/`D_V` as finite additive resets; `dGamma` remains a conditional amplitude-generation kernel until the PDE transfer-law identification is proved; `380` tests + `50,000` causal block states, maximum mass residual `1.421e-14` |
| `31262755985` | physical-energy causal bridge on `84c8652`: exact flat scalar countermodel to raw `dGamma=dT` (`Gamma([0,1/2])=1/2`, `T([0,1/2])=1/4`), exact adjoint-response physical-work split, low-strain energy gate `E1<=e^(2K)(E0+W_HH^++W_R^+)` and clean route `K<=1/30`, `E0,W_R^+<E1/5 => W_HH^+>=8E1/15`; causal layers are weighted by actual `dT_HH=2[Re<c,F_HH,alpha>]_+dt`, while Duhamel supplies same-time parent-pair support; asynchronous synchronization is measure-agnostic; `387` tests + `50,000` states, projection residual `1.066e-14` |
| `31263195439` | recursive coherent witness extraction on `079f2f4`: exact quadratic source disintegration `N(w1,w2)=sum_CD N(A_Cw1,A_Dw2)` and exact actual-work atoms `W_CDE=2 Re<A_Ew3,N(A_Cw1,A_Dw2)>`, with `sum W_CDE` equal to selected HH work; positive/negative Hahn split gives physical binary parent-pair/child event measure, atomic positive mass dominates aggregate positive work, and selected positive cross-cell atoms are excised once into `Xi`; combined with the energy gate, generated mass is `>=8E1/15` before `Xi` and `>=(1-rho)8E1/15` after relative excision; `392` tests + `20,000` POVM/bilinear states, work residual `3.678e-15` |
| `31263818795` | resolved-role Egorov/countermodel on `ba3d345`: exact affine-subtracted identity `(dot Q+[V.grad,Q])f=int K_N(y)[V(x)-V(x-y)-Ay].grad f(x-y)dy`, curvature-only integral-Taylor bound with no higher derivative tower, and explicit strict-lowpass divergence-free shear `V_2=a[sin(rx_1)-sin(2rx_1)/2]`, `r=N/8`, with `grad V(0)=Hess V(0)=0` but `d_1^3V_2(0)=3ar^3`, ruling out center-Hessian-only closure; natural replacement is affine-invariant coherent deformation variance `K_coh^2=Var_gamma[L^-1 grad V(X+Lz)L]`, controlled by coherent curvature energy through Gaussian Poincare; `398` tests + `50,000` states, commutator residual `5.848e-16` |
| `31264233454` | coherent affine projection on `f8164ce`: choose the whole-eddy Gaussian regression `vbar=E_gamma W`, `Abar=E_gamma grad W`; the residual has Hermite degree `>=2`, so OU gap gives `E|R|^2<=K_coh^2/2` and creation/annihilation gives `E|z|^2|R|^2<=7K_coh^2`; full scalar-advection plus vector-amplitude Gaussian-core residual is `<= [1+|q|/sqrt2+sqrt7/2]K_coh`; radius/aspect bounds `q`, and `I_K^2<=0.275568824559 c D_V` on the transition strip, collapsing all higher spatial Hermite deformation into perturbative forcing or critical dissipation rather than new currencies; `403` tests + `50,000` states |
| `31264531046` | first coherent-averaged-strain-source run on `3d508a3` failed in pytest before theorem stress because a fixture compared `190*0.01=1.9000000000000001` to literal `1.9` with exact `==`; fixture/floating semantics only, no theorem formula or constant changed |
| `31264579580` | coherent averaged resolved-strain source on `d4841d4`: exact moving-average identity for `Xdot=barV`, `Ldot=barA L`; exact `dot barA=-barA^2-<a^2>-<HessP>-<grad div R_sgs>+nu<Delta A>-<r.grad A>` and corotational strain version; Gaussian IBP gives `<r.grad A>=<(z.R)(A-barA)>`, so new averaging terms are `<=kappa^2 K_coh^2` and `<=sqrt7 kappa K_coh^2`; normalized Reynolds source weight `<=1.18115356379 D_V`, local quadratic+averaging weight `<=2.0383176489 D_V`, while averaged pressure/SGS/viscous terms inherit existing filtered routes; `408` tests + `50,000` states, identity residual `9.172e-16` |
| `31264981294` | coherent service-or-flat assembly gate on `184c64e`: whole-eddy flatness `F_coh=sqrt(E_H^phys)+(dT)_nonconf+C(qmax)I_K`; for default `Nr_g<=4`, `condL<=567/500`, shell carrier `|k|/N<=e^(2/25)`, `qmax=4.71207563594`, `C=5.65481629117`, so `I_K<5.89468014821e-4` makes full non-affine Gaussian-core connection action `<tau/3` at `tau=1/100`, while larger `I_K` forces critical `D_V`; no frozen Christ-profile persistence is required; the gate returns all threshold-crossing physical roots and delegates primary selection to the first-causal compiler rather than using lexicographic order; `414` tests + `50,000` branch states |
| `31266108674` | complex Young parent marking on `e29ba49`: after one selected frozen scalar/helical cell has weighted efficiency `R` and normalized symbol-freezing error `xi`, exact algebra gives `|T_1(f,g,h)|/A_3>=R-xi` and hence `||f*g||_3/A_3>=R-xi`; Christ's external complex-valued Young stability may therefore mark the actual parent pair by complex Gaussian extremizers once `(1-R)+xi<=delta_Christ(eps)`, with no numerical modulus invented and the same symbol discrepancy charged only once as `Xi`; `425` tests + `50,000` algebra/modulus states |
| `31266550206` | corrected dual-Gaussian analysis registration on `7a7306b`: for a unit-`L^(3/2)` near-Gaussian role, covariance quantization `d_log<=0.4` and shell uncertainty give normalized `N|<f,phi>|^2>=0.207...>1/5`, Leray/outer orthogonal projection preserves the coefficient, and scale coloring plus finite covariance bins gives a depth-independent analysis budget; crucial homogeneity correction: for the actual role the quantum is multiplied by `||f||_(3/2)^2`, so no absolute root mass follows from shape alone; `433` tests + `50,000` states |
| `31266550193` | Bargmann/Moyal root-cell registration on `7a7306b`: a large Gaussian coherent coefficient forces positive Moyal energy in a radius-`sqrt3` intrinsic ball; deterministic maximum-energy cell anchoring preserves the physical parent-slot pushforward and merges reused slots naturally; the normalized cell quantum scales by the squared parent `L^(3/2)` amplitude in the physical role, rather than supplying a false absolute mass floor; `433` tests + `50,000` states |
| `31284828229`, `31284828236` | preferred amplitude-entropy and common-slice registration on `e7181bd`: exact asynchronous margin `67/195`; first-stop adjoint registration gives `|z(s)|>=|z(t)|/4` whenever no earlier HH-generation/classified-residual/material-relink obstruction fires, so two continuing parents retain at least `1/16` of event coefficient product; the scale-critical law `alpha_p1 alpha_p2>=Lambda alpha_c` with `alpha=sqrt(N)|<u,phi>|` implies a transfer-weighted log-amplitude recursion, while Bargmann/Moyal gives `N E_anchor>=beta alpha^2`; log-sum then bounds root entropy and yields the same linear Shannon/Renyi reuse slope `log(48/25)` **without** assuming `N E_root>=eta` for every root; displayed defaults `Lambda=5.97111612285e-4`, `beta=1.43386756899e-5`; `445` tests, `50,000` exact adjoint decompositions and `20,000` adversarial causal DAGs with amplitude ratios up to `7.922e14`; simultaneous registration obstructions are returned as a set and delegated to the single-charge compiler |
| `31285033969` | full physical-energy causal integration on `4108125`: `445` tests plus physical-energy bridge `50k`, recursive coherent work `20k`, service/flat `50k`, complex Young `50k`, dual Gaussian `50k`, Bargmann anchor `50k`, common-slice registration `50k`, amplitude-entropy `20k`, Shannon `50k`, Renyi `50k`, branch compiler `50k`, and master `20,000` traces with worst margin `0`; all modules share actual positive child-energy causal weights and the final conditional status is outer moving-role PDE extraction |
| `31285849905` | exact outer moving-role extraction on `40dcf9d`: for `V=S_(N/4)u`, `h=u-V` and scalar Fourier role `w=Q(t,D)u`, `(dt+L_V-nu Delta)w=Q B(V,V)-Q B(h,h)+(dtQ+[L_V,Q])u`; dual-affine transport and `K<=1/30` keep the selected role above `N/2`, so `Q B(V,V)=0` throughout the slab; `449` tests + `50,000` algebra/support/Egorov states, outer identity residual `5.245e-15` |
| `31286032592` | nonaffine role-interface work on `744463a`: the Heisenberg interface splits as `L_na=K+S`, `K*=-K`, `S*=S`; skew pair work is antisymmetric conservative role flux with zero total, while symmetric interface work is exactly the off-diagonal part of the same resolved strain work; no third source currency and no representation `Xi`; `453` tests + `50,000` complex Hilbert-space partitions |
| `31286361090` | event-anchored hard role / smooth PDE envelope on `1827645`: exact hard Fourier/helical projectors carry physical transfer/Young/Moyal identity, a smooth scalar envelope `Q` satisfies `QP=P`, `(11/20)e^(-1/30)>1/2` preserves low-low exclusion, and `<Pu,phi>=<Qu,Pphi>` exactly; `457` tests + `50,000` states, hard-partition residual `7.573e-16`, coefficient residual `1.472e-16` |
| `31286887248` | physical pair-weighted productivity on `082bbdb`: sharp Young is applied to the **actual positive pair-work density** and KL positivity under normalized `dT` gives `E_dT log(alpha_p1 alpha_p2)>=E_dT log alpha_c+log Lambda_j`; no Duhamel-to-physical pair identification; polynomial pair-cell refinement contributes only a finite geometrically discounted offset; dedicated theorem green |
| `31287088171` | first expanded physical-weighted integration on `47eef1c`: `462` tests and the new outer-role/interface/event-role/physical-pair stack all green before amplitude entropy, Shannon/Renyi, compiler and master; physical-pair KL/productivity margins positive and master `20,000` traces retain worst margin `0` |
| `31287463160` | generated physical first-stop/survival constructor on `6152c5d`: exact `Xi + currencies + HH-regeneration + survivors` mass algebra on actual positive HH work, common-slice/Young failures cannot disappear, and survivor fraction `q` gives the exact conditioned productivity factor `Lambda_surv=q Lambda`; `468` tests + `50,000` physical-work measures, mass residual `4.341e-16` |
| `31287703290` | joint causal-stop projection on `424dfb5`: exact simultaneous first causes are kept as one **unsplit** physical stop; boundary / fixed transfer / true global reset / recurse terminal semantics are invariant under arbitrary positive dummy tie weights, so RN tie splitting is no longer a master hypothesis; `474` tests + `50,000` joint stops, zero weight-invariance failures |
| `31287976953`, `31287976935` | unsplit-master candidate on `0634866` failed in pytest only because a new fixture hard-coded the nonexistent key `recurse_critical` instead of the canonical enum value `sticky_or_scale_critical_recurse`; `475` other tests passed and theorem stress had not started; fixture/wiring provenance only, no theorem formula or constant changed |
| `31288015383` | corrected preferred unsplit generated master on `7e8a99e`: fixture key changed only to `MasterDisposition.RECURSE_CRITICAL.value`; `476` tests + `50,000` recursive physical measures + `50,000` physical-productivity states + `50,000` branch-compiler states all green; minimum positive survivor productivity `3.294e-08` |
| `31288063518` | preferred full physical-energy causal integration on `7e8a99e`: `476` tests; physical energy `50k`, coherent witness `20k`, service/flat `50k`, complex Young `50k`, dual/Bargmann/common-slice `50k` each, outer role/interface/event role/physical pair/recursive witness/joint stop `50k` each, amplitude entropy `20k`, Shannon/Renyi/compiler `50k` each and master `20,000`; physical-pair KL margin `1.161e-4`, log-product margin `1.665`, joint tie-weight invariance failures `0`, master worst margin `0`; artifact stored under `recorded-results/31288063518/` |

| `31249662333`, `31249664073`, `31249665749`, `31249667523`, `31249669387`, `31249671051`, `31249672640`, `31249674254`, `31249675942` | nine-workflow adjoint-causal integration on `f1f0e62`: coherent transfer cells, coherent increment service, service-or-flat rigidity, physical flat episode, causal binary ancestry, transfer-weighted Shannon reuse, Renyi causal reuse, adjoint Kelvin--Duhamel and master all green; adjoint workflow `339` tests + `50,000` checks, master `20,000` traces with worst margin `0` |

The current preferred integrated master regression is run `31288063518` on
`7e8a99e` (`476` tests; the physical-energy/coherent/outer-role/physical-pair/
recursive-witness/joint-stop/Shannon/Renyi/compiler stack is green; master
`20,000` traces, worst margin `0`).  Its artifact is stored under
`recorded-results/31288063518/`.  The canonical cross-cell regression
`31258603563` and accompanying six workflows remain green historical provenance,
as does the earlier twelve-workflow coherent-reservoir stack
`31230314319`--`31230338382`.
The full source/sideband/ancestry integration on that exact SHA is
`31227852933`, `31227856265`, `31227859724`, `31227862808`, `31227865715`,
`31227868626`, `31227871490`, `31227874557`, `31227877505`, `31227881162`,
`31227884346`, `31227887757`, `31227891065`, `31227894352`: physical H1
source, filtered SGS/viscous collision, Onsager increments, source-weighted
episode, sticky aspect, H1 bridge/no-escape, H3 no-escape, pair ancestry,
curvature source, affine critical grain, relative polarization, localized
pressure and master were all green.  The earlier nine-workflow affine master `31182796814`,
eleven-workflow affine
master `31182017032`, relative-polarization master
`31175212421`, integrated helical run
`31171950823` and frequency/pressure bridge run `31166171000`
remains useful provenance.  The earlier
recorded master artifact `31154025683/` remains the canonical stored
episode-trace artifact.  The preferred smooth physical-flux artifacts are
`recorded-results/31161316034/` and `recorded-results/31161626056/`.  The
preferred single-edge/physical-weight certificate is
`recorded-results/31161914134/`.  The localized pressure-work artifact is
`recorded-results/31162700474/`; its new multipole/fresh-mass continuation is
`recorded-results/31164911526/`.  The crossing-scale extraction certificate is
`recorded-results/31164771160/`, and the exploratory one-shot transfer-profile artifact is
`recorded-results/31165172325/`; its preferred Arb narrow-shell continuation is
`recorded-results/31165654509/`.  The preferred smooth-symbol freezing artifact
is `recorded-results/31165838379/`.  The preferred affine-spacetime dynamics artifacts are
`recorded-results/31168393589/`, `recorded-results/31168671304/`,
`recorded-results/31168888413/`, and `recorded-results/31169128097/`, plus the intrinsic 3D continuation `recorded-results/31170015795/`.  The preferred actual-packet PDE artifacts are `recorded-results/31186597005/`, `recorded-results/31187104071/`, `recorded-results/31187443147/`, `recorded-results/31187637446/`, `recorded-results/31187768828/`, `recorded-results/31188118452/`, and `recorded-results/31189036483/`; the certificate-semantics failure is `recorded-results/31186600162/`.  The preferred helical artifacts are
`recorded-results/31171018560/`, `recorded-results/31171127537/`,
`recorded-results/31171360107/`, `recorded-results/31171484247/`, and the
symplectic correction `recorded-results/31172607991/`, good-core frame
regularity artifact `recorded-results/31173321354/`, preferred forced relative-
polarization artifact `recorded-results/31174812731/`, and localized packet
certificate `recorded-results/31174612248/`.  The preferred affine-residual
artifacts are `recorded-results/31179739773/`, `recorded-results/31179827015/`,
`recorded-results/31180257124/`, `recorded-results/31180506627/`,
`recorded-results/31180882083/`, `recorded-results/31181076691/`, and
`recorded-results/31182421821/`; the two failed
interval-boundary provenance runs are `recorded-results/31168205564/` and
`recorded-results/31168303213/`.  Final integration runs are
`31166152074`, `31166155045`, `31166158218`, `31166160711`, `31166163414`,
`31166165985`, `31166168659`, and `31166171000`.  The numerical full-Mellin search is recorded
in `recorded-results/31160779428/`; its qualitative counterexample is separately
Arb-certified in run `31161914134`.

## 18. Current research frontier

The H1 sideband/source problem, cubic SGS increment problem, and the previous
mild/high-aspect gap have all changed form.

On the full physical transition strip

\[
\boxed{\operatorname{cond}L\le567/500,}
\]

hook curvature remains transfer-visible:

\[
\boxed{Q_{pol}\ge\frac1{4000}\|B^H\|^2,}
\qquad
\boxed{\sum_i\|F_i^{H1}\|^2\ge\frac1{16000}\|B^H\|^2.}
\]

Outside source/feedback/large-daughter/pair exits this gives

\[
\boxed{
\operatorname{Def}\ \text{or}\ R_{pair}
\ge I_1^2/28\,800\,000,
}
\]

and, after the irreducible H3/H1 split,

\[
\boxed{
\operatorname{Def}\ \text{or}\ R_{pair}
\ge I_B^2/115\,200\,000.
}
\]

For stronger aspect, covariance dynamics supplies the sticky alternative

\[
\boxed{
\kappa(t_1)>567/500
\Longrightarrow
\kappa(t_0)>21/20
}
\]

on a low symmetric-strain lifetime, unless ancestry breaks and the grain is
fresh/relinked.  Fresh high aspect pays physical radius occupancy

\[
\boxed{
N\int_{E_2}|u|^2>\frac15\kappa^{1/3},
}
\]

not an aspect/Young deficit.

H1 dephasing now reaches actual Navier--Stokes sources.  Differentiated SGS and
viscous-fourth sources no longer remain high-derivative mysteries:

\[
\boxed{
S_R\to R\to |\delta_\ell u|^3
\to
\{\text{low/base mass},\text{entropy/cycles},\text{dissipation}\},
}
\]

\[
\boxed{
S_\nu\to N^{-1}\|\nabla V\|_2^2
\to\text{dissipation}.
}
\]

The key Onsager/LP estimate is

\[
\boxed{
Q_N
\le
g_1(C_{LP}C_B)^3
\left[
\frac43(\beta R_G)^2\mu_{low}^{max}
+2\mathfrak d_{>N}
\right]^{3/2}.
}
\]

Because the SGS source-to-increment power is `3/2` and this collision takes the
`2/3` power, the final mass/enstrophy currency is linear in instantaneous source
density.  Source weight, rather than the duration of a superlevel set, can be
pigeonholed; temporal concentration is therefore not a free SGS/viscous exit.

Ancestor-reservoir synchronization is now theorem-level.  On low-strain material
reuse, one low-band reservoir has a service half-life `<1/2` per generation and
one fixed low-low pressure pair has a one-third-life `<1/3`.  The phase-space
label cannot be reset freely because `q=L^T k` is a material invariant.  Spatial
replication at one generation is also band-limited: a maximal `N^-1`-separated
source family has additive `ell^(3/2)` SGS/pressure charge or `ell^2` viscous
charge, while sub-resolution multiplicity is merged into one source cluster.

Whole-old-pool relabeling is now finite-atomic theorem-level.  Given a frame
budget `P`, the entire material old pool has service capacity `<2^{-q}` up to its
explicit generation-zero coefficient, so a uniform service threshold forces new
spectral capacity after finitely many generations.  Active relinking triads then
obey `(n-1)+beta=2m`, hence every connected relinking component is fresh-rich or
cycle-rich.

For equal-covariance coherent Gaussian pools, the required coefficient budget is
also explicit rather than assumed.  The affine phase coordinate

\[
\zeta=(L^{-1}X/2,L^Tk)
\]

has overlap kernel `exp(-|Delta zeta|^2/2)`, and a 5-separated synthesis family
has

\[
\boxed{47I/50\le G\le53I/50},
\qquad
\boxed{\sum|c_a|^2\le(50/47)\|\sum c_ag_a\|_2^2}.
\]

At the analysis level even the discrete Riesz budget can be bypassed: the affine
coherent transform has the exact Moyal identity, so any measurable phase-space
partition carries positive cell masses with total exactly `||u||_2^2`.  Thus the
analytic energy ledger is `P=1`; the Riesz theorem is needed only when an actual
discrete Gaussian synthesis is chosen.

Nonlinear transfer-to-coherent-cell assignment is now theorem-level at the
physical-work and SGS-service levels.  Polarized Moyal decomposes actual band
work exactly; SGS increments generate positive coherent edges between `C` and
`C-r`; small covariance-window changes have a deterministic `Xi_cov`; and a
sticky old-pool coherent epoch has finite length.  On the differentiated-SGS H1
branch the `rho^(3/2)` source law followed by the `2/3` coherent collision is
linear in source weight, so temporal concentration is not an escape.

Uniform service-or-flat rigidity is now theorem-level in the physical packet
model.  At `tau=1/100`, absence of every named positive branch forces
`F_K<=1/100`; a synchronized flat step has an explicit barycentric erosion rate
`kappa_0>0.17`.  Equal parent marginals are no longer assumed: a flat block pays
collision entropy `>1/200` or exposes two separated parent cores, each carrying
at least `7/9` of its marginal transfer mass.

The new decisive mechanism is **causal binary ancestry versus sticky reuse**.
Under a synchronized causal completion, signed-good parents satisfy

\[
3/5<N_{parent}/N_{child}<5/8,
\]

and the exact layered law

\[
{n_0\over2^L}=\prod_j(1-\rho_j)
\]

shows that binary fresh-parent creation outruns the allowable scale dilation.
Finite coherent root energy forces a reuse-rich layer `rho_j>=1/4` once the
causal depth passes the explicit `log(36/25)` threshold.  This resolves the old
energy-only loophole at finite synchronized level: fresh companions cannot be
inserted forever without either exponentially expanding their causal ancestry or
merging it into reuse cycles.

The transfer-weighted causal algebra is now theorem-level at synchronized
finite depth.  Shannon reuse information removes the negligible-edge caveat, and
the Renyi collision telescope gives a layer which already routes into existing
component-entropy / hidden-pair currencies.  The exact `1/2` two-parent-role
baseline remains free throughout.

The one-step **adjoint Kelvin--Duhamel causal gate is now theorem-level**.  For a
selected divergence-free packet coefficient `c` and the backward dual solving
`dot psi=-G^* psi`, where `G` contains the same common affine/Kelvin transport
and bulk viscosity already quotiented in the packet equation,

\[
\boxed{\frac d{dt}\langle\psi,c\rangle
=\langle\psi,F_{HH}+R_{class}\rangle.}
\]

Hence on one child slab

\[
\boxed{z_1=z_0+I_{HH}+I_R.}
\]

Choosing the terminal dual along the terminal coefficient, `A=|z_1|`, gives the
clean causal alternative

\[
\boxed{|z_0|\ge A/4}
\quad\text{or}\quad
\boxed{|I_R|\ge A/4}
\quad\text{or}\quad
\boxed{|I_{HH}|\ge A/2}.
\]

Thus common affine strain, pressure already removed by Leray, and bulk viscosity
are not charged again as new generation.  On the generated branch a single phase
aligned with `I_HH` makes the positive quadratic parent-pair masses dominate
`|I_HH|`, so no pointwise persistence hypothesis is needed merely to define a
positive one-step causal-generation law.  Signed-good scale geometry also gives
`64/25<T_parent/T_child<25/9`; a half-child slab carrying half the generation
mass has parent **natural** backward windows with common overlap longer than
`103/50 T_child`.  This is only geometric window synchronization, not a claim
that the extracted parent packets remain coherent on that whole interval.

The **asynchronous Duhamel-slab synchronization geometry is now theorem-level**.
After the first half-slab pigeonhole, parent events satisfy `alpha_1<=25/128`.
Choosing the common slice `s_j=a_j-(2/5)T_j^min` gives

\[
\boxed{\alpha_{j+1}\le\frac{25}{64}(\alpha_j+2/5)},
\qquad
\boxed{\alpha_j\le10/39\quad(j\ge1)}.
\]

The sharp common-window margin is `67/195`; `10/39` is the recurrence fixed point.
Only generated nodes continue; inherited and classified-residual nodes stop in
their existing ledgers.  Thus no frozen-packet persistence hypothesis is required
to put generated events on common physical slices.  Because the next generated
support may begin anywhere inside the previous common interval, the corrected
reference-slice progress is

\[
\boxed{s_j-s_{j+1}\ge\frac{1792}{4875}T_j},
\qquad
\boxed{\Delta s_L\ge\frac{1792}{7605}T_0[(64/25)^L-1]}.
\]

so finite-time interior ancestry reaches `t=0` after finite depth unless it stops
earlier.  At the initial surface the exact adjoint gate is truncated; this is an
initial-boundary root, not a fresh interior grain.  For `u_0 in Hdot^m`, band-`M`
roots satisfy `#R_M<=eta^-1 M^(1-2m)||u_0||_Hdot^m^2`.

Canonical material-label registration is now theorem-level at representation
level.  Use one nested dyadic address of the intrinsic coherent coordinate `zeta`
for Duhamel, Shannon/Renyi, Hodge/resistance, holonomy and old-pool ancestry.
Common Kelvin transport and nested refinement have zero charge.  With
`h_j=h_0 2^-j` and `delta_j=delta_0 2^-j`,

\[
\boxed{\Xi_{sym}\le2A_3L_*B_*h_0},
\qquad
\boxed{\Xi_{cov}\le\sqrt2\,\delta_0E_{global}}.
\]

Actual material cell switching remains a physical Moyal relink/backflow/fresh
event; it is not hidden as representation error.

The **physical multi-currency master telescope is also theorem-level**.  If
additive resource class `r` has total budget `B_r` and every primary-charged reset
uses at least `b_r`, then `N_A<=sum_r B_r/b_r`.  Flat erosion and
`L=N_F+N_T+N_A` give

\[
\boxed{
N_T\ge
\frac{\kappa_0L-P_{max}-Z}{\kappa_0+P_{max}}
-\sum_r\frac{B_r}{b_r}.
}
\]

Hence

\[
\boxed{
-\log\prod_{j<L}R_j
\ge c_{eff}L-O(1)-\Xi,
\qquad
c_{eff}=\frac{c_0\kappa_0}{\kappa_0+P_{max}}>0.
}
\]

Only reset classes with a proved **scale-independent threshold in a globally
bounded resource** change merely the finite offset.  Critical packet mass `N E`
and normalized dissipation `D_V` do not qualify from energy alone: their physical
costs scale like `1/N` and are summable along one geometric cascade.  They remain
sticky/causal ancestry currencies rather than automatic finite-count resets.

The high-strain branch is now quantitatively identified rather than left as raw
action.  For `V=S_(N/4)u`,

\[
\boxed{
K:=\int\|S\|_{op}dt>1/30
\Longrightarrow
D_V:=N\int\|\nabla V\|_2^2dt>\frac{32\pi^2}{75c}.
}
\]

This is a critical dissipation currency, not a uniform global reset count.

The objective-strain source branch has also been corrected to the **resolved
Kelvin transporter**.  With `A=grad V`,

\[
\boxed{
\mathring S_V
=-S^2-\Omega^2+[S,\Omega]
-\nabla^2P-\operatorname{sym}\nabla\nabla\cdot R+\nu\Delta S.
}
\]

Because `V` is supported in `B_(N/4)` and `P,R` in `B_(N/2)`, Arb gives the
clean filtered pressure route

\[
\boxed{
\rho_P\le\frac{\mu_V}{5700}+\frac{\|R\|_{3/2}}{380}.
}
\]

Thus pressure-Hessian strain dephasing forces `mu_V>=2850 rho_P` or
`||R||_(3/2)>=190 rho_P`; the latter returns to Germano/Onsager coherent service.
The old unspecified raw-pressure near coefficient is not needed for the actual
resolved packet transporter.

Finally, the canonical global coherent ancestry no longer requires a compact
physical window.  For divergence-free `u`, self-adjoint Leray projection gives

\[
\boxed{
\langle u,\mathbb P(g_ze_a)\rangle=\langle u,g_ze_a\rangle,
}
\]

so projected coherent probes remain a Parseval/Moyal frame on `L^2_sigma`, while

\[
\boxed{\langle\mathbb P(g_ze_a),\nabla p\rangle=0}
\]

exactly.  Hence compact moving-window commutators, curvature-balanced cutoff
errors and localized pressure boundary work are not part of the **canonical global
master `Xi`**; Leray/helical distortion is already the summable frequency-symbol
representation error.  Compact windows remain available for genuinely local CKN
diagnostics.

The canonical **physical cross-cell `Xi` summability mechanism is now theorem-level**.
The positive coherent localization operators reconstruct each role and the
trilinear work exactly, so there is no continuous-to-discrete synthesis error.
On the signed-good physical-transfer core the single-edge defect moment obeys

\[
\boxed{\mathbb E_{phys}\mathcal D\le(106/25)\,\epsilon.}
\]

A moat placed in the scalar transfer defect, rather than in a coordinate grid,
gives at depth `j`

\[
\boxed{\eta_{cross,j}\le M_j^{-1}+2\overline{\mathcal D}/R_j.}
\]

With `M_j=M_0(j+2)^2` and `R_j=R_0(j+2)^2`,

\[
\boxed{\sum_j\eta_{cross,j}\le\frac{13}{20}
\left(M_0^{-1}+2\overline{\mathcal D}/R_0\right).}
\]

Thus the actual omitted cross-component **physical child-transfer mass** is
summable and tunably small, with no packet-count factor or Gaussian-synthesis
tail.  If the low-defect interaction hypergraph remains connected, the exact
`(n-1)+beta=2m` incidence identity routes it to fresh-rich or cycle-rich ancestry
instead of creating an interface escape.

The **single-charge physical branch compiler is now theorem-level at the
positive-transfer-measure / causal-quotient level**.  Raw block predicates are not
assumed mutually exclusive: H1/H3 dephasing is quotiented with the physical
source which caused it, high strain with the critical `D_V` certificate, pair
rescue with its reuse/cycle destination, Duhamel `R_class` with its underlying
source/interface term, and `t=0` with an absorbing boundary state.  The defect
moat first excises `dT_Xi`; a fixed transfer loss then stops the retained measure;
otherwise the first physical causal root owns the charge.  Exact simultaneous
hits of genuinely independent roots are not broken by theorem-name priority: the
positive stopping weights partition the tied transfer measure by normalized
Radon--Nikodym weights.  Hence `dT=dT_Xi+sum_r dT_r` exactly, with no duplicated
physical transfer mass.  Critical source/mass/dissipation destinations remain
marked **recurse**, not finite resets.

The proposed raw Duhamel-to-energy kernel identification has now been **ruled out
as unnecessarily strong**.  Even in the perfectly flat scalar coefficient model
`G=R=0`, `c(0)=0`, `c_dot=1`, normalized adjoint amplitude generation is
`dGamma=dt` while normalized physical child-energy transfer is `dT=2t dt`.  The
correct exact relation is an energy lift: for any adjoint response `psi`,

\[
2\Re\langle c,F\rangle
=\frac{2}{\|\psi\|^2}\Re\left(\overline{\langle\psi,c\rangle}\,\langle\psi,F\rangle\right)
+2\Re\langle c_\perp,F\rangle.
\]

More importantly, the selected-child energy balance directly supplies the
master-facing causal weights.  On the low-strain branch,

\[
E_1\le e^{2K}(E_0+W_{HH}^++W_R^+).
\]

Since `K<=1/30` implies `e^-2K>=14/15`, if material inherited energy and
classified residual positive work are each below `E1/5`, then

\[
\boxed{W_{HH}^+\ge8E_1/15}.
\]

The exact quadratic source atoms therefore carry the positive physical causal
measure `dT_HH=2[Re<c,F_HH,alpha>]_+dt`.  Duhamel is kept for the same-time
parent-pair support/interaction picture, but its raw amplitude weights are not
used by the master.  The half-slab and parabolic synchronization proof is
measure-agnostic, so Shannon/Renyi now receives the same physical transfer law
used by Hodge/flat bookkeeping **by construction**.  Preferred certificate:
`31262755985` on `84c8652` (`387` tests + `50,000` states; projection residual
`1.066e-14`, exact half-time CDF gap `1/4`).

The **binary coherent witness extraction is now theorem-level once the outer
selected PDE roles exist**.  For positive coherent localization operators
`sum_C A_C=I`, bilinearity gives

\[
\mathcal N(w_1,w_2)=\sum_{C,D}\mathcal N(A_Cw_1,A_Dw_2)
\]

and the exact child-work atoms

\[
\boxed{W_{CDE}=2\Re\langle A_Ew_3,\mathcal N(A_Cw_1,A_Dw_2)\rangle}
\]

satisf `sum W_CDE=2 Re<w3,N(w1,w2)>`.  Their positive/negative Hahn split is
physical transfer/backscatter, not an arbitrary packet graph; each positive atom
already has exactly two parent material labels and one child label at one physical
time.  Atomic positive mass dominates the positive part of aggregate HH work, so
the physical-energy gate yields at least `8E1/15` positive binary mass before
cross-cell excision.  If the physical defect moat removes relative positive mass
`rho`, the retained generated law still has at least `(1-rho)8E1/15`.  Coherent
pieces are **not** treated as compact Fourier projections: the outer role carries
scale/helicity and representative errors remain in the existing symbol/covariance
`Xi`.  Preferred certificate: `31263195439` on `079f2f4` (`392` tests + `20,000`
states; work reconstruction residual `3.678e-15`).

The **outer moving-role Egorov defect is now identified exactly**, but it forces a
redesign of the curvature observable.  If `Q_t` is transported by the affine jet
`A=grad V(X)`, then

\[
\boxed{(\partial_tQ+[V\cdot\nabla,Q])f(x)
=\int K_N(y)[V(x)-V(x-y)-Ay]\cdot\nabla f(x-y)dy.}
\]

The integral Taylor remainder depends only on velocity curvature over the
coherent region, so there is no uncontrolled higher-derivative tower.  However a
single center Hessian is **not** sufficient: the strict-lowpass divergence-free
shear `V_2=a[sin(rx_1)-sin(2rx_1)/2]`, `r=N/8`, has `grad V(0)=Hess V(0)=0` but
`d_1^3V_2(0)=3ar^3`.  Thus any bridge that prices the full moving-role residual
only through the osculating center tensor `B(X)` is false.

The physically natural replacement is the affine-invariant coherent deformation
variance

\[
\boxed{\mathcal K_C^2
=\operatorname{Var}_\gammaig(L^{-1}\nabla V(X+Lz)Lig).}
\]

It vanishes for genuinely affine resolved flow and, by Gaussian Poincare,

\[
\boxed{\mathcal K_C^2
\le\mathbb E_\gamma
\|L^{-1}(\nabla^2V)(X+Lz)[L,L]\|_F^2.}
\]

This measures resolved deformation variation across the actual coherent eddy,
not at one privileged point.  Preferred certificate: `31263818795` on
`ba3d345` (`398` tests + `50,000` states; exact commutator residual
`5.848e-16`).

The **full coherent spatial deformation spectrum now collapses to one physical
observable** if the common affine gauge is chosen by Gaussian regression over the
whole eddy.  In intrinsic coordinates `W=L^-1V(X+Lz)`, set

\[
\bar v=\mathbb E_\gamma W,
\qquad
\bar A=\mathbb E_\gamma\nabla_zW,
\qquad
R=W-\bar v-\bar A z.
\]

Gaussian integration by parts makes this the orthogonal projection onto Hermite
degrees `0+1`; hence `R` has degrees `>=2` and

\[
\boxed{\mathbb E|R|^2\le\tfrac12\mathcal K_C^2},
\qquad
\boxed{\mathbb E|z|^2|R|^2\le7\mathcal K_C^2}.
\]

For the Gaussian core `psi=g exp(iq.z)`, scalar residual advection together with
vector rapid-distortion mismatch obeys

\[
\boxed{
\|F_{nonaff}\|_2/\|\psi\|_2
\le(1+|q|/\sqrt2+\sqrt7/2)\mathcal K_C.}
\]

On the scale-matched transition-aspect branch `q` is uniformly bounded by radius
and aspect.  Large accumulated coherent deformation is already critical
dissipation:

\[
\boxed{I_K^2\le0.275568824559\,cD_V.}
\]

Thus higher spatial Hermite degrees do **not** require new master currencies just
to control the moving Gaussian core: small `I_K` is perturbative, large `I_K`
routes to critical `D_V`.  Preferred certificate: `31264233454` on `f8164ce`
(`403` tests + `50,000` states).

The **coherent averaged affine jet now has an exact resolved Navier--Stokes source
calculus**.  Move the analysis eddy by

\[
\dot X=\bar V,
\qquad
\dot L=\bar A L,
\qquad
\bar A=\langle\nabla V\rangle_\gamma.
\]

For `r=V-barV-barA(x-X)` the moving-average identity is

\[
\frac d{dt}\langle f\rangle
=\langle D_t^Vf-r\cdot\nabla f\rangle,
\]

and therefore

\[
\boxed{
\dot{\bar A}
=-\bar A^2-\langle a^2\rangle
-\langle\nabla^2P\rangle
-\langle\nabla\nabla\cdot R_{SGS}\rangle
+\nu\langle\Delta A\rangle
-\langle r\cdot\nabla A\rangle.}
\]

The two terms created solely by averaging are already coherent-deformation
currency.  Gaussian integration by parts gives

\[
\langle r\cdot\nabla A\rangle
=\langle(z\cdot R)(A-\bar A)\rangle,
\]

so

\[
\|\langle a^2\rangle\|\le\kappa^2K_{coh}^2,
\qquad
\|\langle r\cdot\nabla A\rangle\|
\le\sqrt7\,\kappa K_{coh}^2.
\]

On the transition strip their integrated normalized source weight is
`<=1.18115356379 D_V`; including the averaged quadratic local term gives
`<=2.0383176489 D_V`.  These are critical dissipation, not new additive reset
resources.  Probability averaging does not enlarge the existing global `Linf`
bounds for filtered pressure, differentiated SGS or viscosity, so their
mass/SGS-service/`D_V` routes remain unchanged and are charged only once.
Preferred certificate: `31264579580` on `d4841d4` (`408` tests + `50,000`
states; corotational identity residual `9.172e-16`).  The preceding run
`31264531046` was only an exact-float fixture failure and is retained as
provenance.

A **coherent service-or-flat assembly gate is now theorem-level without frozen
Gaussian-profile persistence**.  On a scale-matched transition-aspect block use

\[
\boxed{\mathfrak F_{coh}
=\sqrt{E_H^{phys}}+(dT)_{nonconf}+C(q_{max})I_K.}
\]

For the default clean split `Nr_g<=4`, `condL<=567/500` and outer shell carrier
`|k|/N<=e^(2/25)`,

\[
q_{max}=4.71207563594,
\qquad
C(q_{max})=5.65481629117.
\]

At `tau=1/100`, the existing transfer/objective-strain thresholds make the Hodge
and common-strain pieces each `<tau/3`, while

\[
\boxed{I_K<5.89468014821\times10^{-4}}
\]

makes the **entire** spatial non-affine Gaussian-core connection action `<tau/3`.
If this threshold fails, coherent deformation forces positive critical `D_V`;
large averaged source variation, high strain, high aspect, large radius, phase
holonomy and direct transfer loss keep their already proved physical routes.
The gate explicitly returns every threshold-crossing causal root and chooses no
primary by theorem-name order; the first-causal single-charge compiler owns that
decision.  Existing H1/H3 theorems remain sharper optional diagnostics, but are
no longer a profile-persistence hypothesis required by the assembly.  Preferred
certificate: `31264981294` on `184c64e` (`414` tests + `50,000` branch states).

The preferred causal frontier has now changed again.  Young near-extremality is
separately homogeneous in the parent roles, so an absolute root critical-mass
floor is still **not** inferred from shape.  But the stronger pointwise
Duhamel-parent product law is no longer needed either.

At a physical transfer event use an exact hard Fourier/helical role `P` and a
smooth PDE envelope `Q` with

\[
\boxed{QP=P.}
\]

The event role carries actual transfer, Young and Moyal identity; the envelope
only carries the PDE between causal slices.  On the low-strain branch its radial
lower edge remains

\[
\frac{11}{20}e^{-1/30}N>\frac12N,
\]

so the exact outer moving-role equation is

\[
\boxed{
(\partial_t+\mathcal L_V-\nu\Delta)Qu
=-Q\mathbb P\nabla\!\cdot(h\otimes h)+R_Q,
}
\qquad h=u-S_{N/4}u.
\]

There is no low--low forcing throughout the slab.  The nonaffine Heisenberg term
`R_Q` has exact energy provenance: its skew-adjoint part is antisymmetric
conservative role transfer, while its self-adjoint part is the off-diagonal
piece of the same physical strain work.  It is neither a new generation source
nor representation `Xi`.

For the generated positive physical HH work itself, let `M_j` be the retained
hard pair-cell count at depth `j`.  Sharp Young applied to the actual positive
pair-work density and KL positivity under normalized `dT` give the preferred
productivity law directly in the physical measure:

\[
\boxed{
\mathbb E_{d\mathcal T_j}
\log(\alpha_{p_1}\alpha_{p_2})
\ge
\mathbb E_{d\mathcal T_j}\log\alpha_c
+\log\Lambda_j.
}
\]

No Duhamel parent-pair reweighting occurs.  If `M_j<=M_0(j+3)^p`, then
`Lambda_j` loses only the factor `M_j^{-1}`, and the binary recursion sees it
with geometric weights.  Hence

\[
\boxed{
\ell_0
\ge
\sum_{j=0}^{L-1}2^{-(j+1)}\log\Lambda_j
+2^{-L}\log\alpha_L.
}
\]

The Bargmann/Moyal anchor still gives

\[
N_rE_r\ge\beta\alpha_r^2,
\qquad \beta=1.43386756899\times10^{-5},
\]

so the Shannon/Renyi reuse lower bound becomes

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

Polynomial symbol-cell refinement therefore changes only a finite offset, not
the linear `log(48/25)` reuse slope.

The recursive generated branch is now also measure-level single-charge.  After
one support-level `Xi` excision, every physical HH event either has a first
physical stop, an earlier HH-regeneration stop, or survives Young/phase and
common-slice registration.  If the survivor set carries fraction `q_j` of the
retained generated physical work, the same KL proof restricted to that set gives
exactly

\[
\boxed{\Lambda_{j,\mathrm{surv}}=q_j\Lambda_j.}
\]

Thus `q_j>=1/2` costs at most `log 2` at that layer; over a continuing ancestry

\[
\sum_{j\ge0}2^{-j}\log q_j\ge-2\log2
\]

whenever every `q_j>=1/2`.  If some `q_j<1/2`, a majority of that layer's actual
HH work has already hit a named causal stop or earlier regeneration and is not a
free continuing branch.

Exact simultaneous first causes no longer create a common-unit weighting
problem.  The preferred master keeps the full physical event mass **unsplit**
and retains the entire joint cause set as provenance.  `t=0` is absorbing; any
already-certified fixed transfer/reuse/sideband cost is terminal; otherwise a
true globally bounded uniform resource may terminate; all remaining
source/critical-dissipation/material-new-ancestry/HH-regeneration combinations
stay in one recursive-critical state.  Therefore arbitrary tie weights cannot
change the master.  The older Radon--Nikodym split remains optional fine
bookkeeping only.

Preferred integrated certificate: `31288063518` on `7e8a99e` (`476` tests plus
the full physical-energy/coherent/outer-role/physical-pair/recursive-witness/
joint-stop/Shannon/Renyi/compiler stack; master `20,000` traces, worst margin
`0`).  The stored artifact is `recorded-results/31288063518/`.

A **smooth-SGS measurable first-exit theorem now closes the local first-hit
problem without manufacturing a common vector clock**.  The physical information
appears in its native causal order: event/support selection, smooth slab
observables, actual positive HH work, backward common-slice registration, then
ancestry information.  For each native-unit continuous observable `f_r`, the
closed hit `inf{t:f_r>=theta_r}` or strict-superlevel debut
`inf{t:f_r>theta_r}` is Borel on path space.  For a finite physical family,

\[
\boxed{\tau=\min_r\tau_r,\qquad J=\{r:\tau_r=\tau\}}
\]

is a measurable first-hit time with an unsplit joint cause set.  Independent
strictly increasing changes of units leave the stopping set unchanged; source,
strain, Moyal energy and phase are never scalarized against one another.

Material coherence is anchored in the intrinsic coordinate
`zeta=(L^{-1}X/2,L^T k)` rather than repeatedly re-selected by a winning-cell
index.  If `F(t)=V_{g(t)}u(t)` and `C` is the anchored cell, then

\[
E_C(t)=\|1_CF(t)\|_2^2,\qquad
|E_C'|\le2\sqrt{E_C}
(\|g\|_2\|u_t\|_2+\|g_t\|_2\|u\|_2).
\]

Hence its physical Moyal content is absolutely continuous; dyadic cell
boundaries are Moyal-null, so optimizer/cell-name chatter cannot create material
relink.  Helical phase is represented by the branch-free circle distance
`acos(Re h)`, not by a principal-angle branch.

The existing common-slice identity then gives local no-hit exhaustion.  If no
residual/source, earlier HH regeneration, material relink or `t=0` boundary stop
hits, `|z(s)|>=|z(t)|/4`, so the generated no-hit event is exactly a registered
survivor.  If none of the whole-block service exits hits, the existing theorem
returns Kelvin-flat continuation.  Therefore every already-selected smooth-SGS
block has locally only the two uncharged alternatives already present in the
physics: registered generated survivor or certified Kelvin-flat continuation.

Dedicated certificate `31290886164` on `ed7773d` passed `484` tests plus `50,000`
first-exit/Moyal/phase/registration regressions, with zero monitor-order failures,
zero independent-unit-rescaling failures, minimum Moyal-rate margin `4.039e-09`
and zero branch-cut phase gap.  Full integration `31290886140` on the same SHA
passed `484` tests and the entire causal stack; physical-pair KL margin remained
`1.161e-04`, log-product margin `1.665`, Shannon/Renyi amplitude-entropy margin
`1.958e-01`, and master `20,000` traces had worst margin `0`.  Stored artifacts:
`recorded-results/31290886164/` and `recorded-results/31290886140/`.

Two further exact continuum structures remove the apparent hard-event re-entry machinery.

First, a registered common-slice parent stays in the **same smooth material carrier** `w=Q u`.  Since `z=<w,psi>`,

\[
\|w\|_2^2\ge |z|^2/\|\psi\|_2^2.
\]

Re-anchoring at the common slice is only composition of common affine/Kelvin transport and leaves intrinsic `zeta` exactly fixed.  It therefore creates no relink, no fresh ancestry and no new representation `Xi`.  If `r` denotes the signed physical HH work density and `0<=q<=1` the scalar carrier symbol,

\[
\left[\int q^2r\right]_+\le\int q^2[r]_+\le\int[r]_+.
\]

Thus smooth-carrier generation is dominated by the existing actual positive physical HH law.  Hard Fourier/helical roles are read only at real nonlinear interaction events; they need not persist or be reselected at an intermediate slice.

Second, parent-scale renewal does not create a cutoff interface.  For every resolved field `V`,

\[
-L_V(Q u)+Q B(V,V)-Q B(u-V,u-V)+(L_VQ-QL_V)u=-Q B(u,u),
\]

so replacing `S_(N/4)u` by `S_(N_p/4)u` exactly repartitions the same NS interaction between transport, HH source and Heisenberg interface.  With `3/5<N_p/N<5/8`, the relayed envelope satisfies

\[
|\xi|/N_p\ge(22/25)e^{-1/15}>1/2,
\]

and parent-scale low-low exclusion remains strict; the parabolic lifetime remains `64/25<T_p/T<25/9`.

Dedicated run `31291697710` passed `489` tests plus `50,000` smooth-carrier relay states, with worst two-step zeta residual `2.707e-15`.  Corrected cutoff run `31291932513` passed `492` tests plus `50,000` cutoff/scale states, with worst old/new identity residual `1.099e-14`, worst full-NS residual `8.563e-15`, and minimum renewed low-low gap `3.242566e-01`.  Stored artifacts: `recorded-results/31291697710/` and `recorded-results/31291932513/`.  Precursor `31291900197` was certificate-fixture only: all `492` tests passed before an inappropriate Arb equality check on exact rational lifetime endpoints; `4fcdc8c` changed only that check to exact `Fraction` arithmetic.

The **generated-survivor return is therefore no longer an unresolved hard-packet or cutoff problem**.

High strain now also has a direct physical seed.  For `V=S_(N/4)u`, decompose the resolved ball into `A_j={M_j/2<|xi|<=M_j}`, `M_j=(N/4)2^-j`, and set `mu_j(t)=M_j||P_j u(t)||_2^2`.  Since `sum_jM_j=N/2`, shell-time atoms with `mu_j<mu_*` carry at most `c mu_*/2` of normalized dissipation on a lifetime `cN^-2`.  Combining this with `D_V>=D_*=32pi^2/(75c)` and `mu_*=D_*/c` gives

`D_V(mu_j>=mu_*)>=D_V/2`.

Thus at least half of the actual resolved-dissipation law carries a lower-frequency critical ancestor, with `M_j<=N/4` and natural lifetime at least `16` child lifetimes.  This is a dissipation disintegration, not a shell argmax.

The same law admits an exact physical-space representation using the NS heat kernel `H_N` at time `1/(2N^2)`.  On the strict resolved support,

`e^(-1/32)||grad V||_2^2 <= N^2 int H_N(r)||delta_rV||_2^2dr <= ||grad V||_2^2`.

Hence the heat-increment service `S_heat=N^3 int dt int H_N||delta_rV||_2^2dr` lies between `e^(-1/32)D_V` and `D_V`, and at least `e^(-1/32)/2>0.48` of the full service remains on the same critical shell-time set.  Exact Moyal disintegration and translation covariance then give positive coherent spatial edges without selecting a cell or dividing shell energy among packets.

Dedicated `31292418910` on `bfda899` passed `496` tests plus `50,000` shell/time laws; dedicated `31292625136` on `a08471b` passed `501` tests plus `50,000` heat-defect/ancestor laws, with relative Moyal residual `6.278e-16` and covariance residual `8.016e-15`.  Stored artifacts: `recorded-results/31292418910/` and `recorded-results/31292625136/`.

Material ownership of that positive heat-edge law is now exact as well.  For the two coherent endpoints `(X,k)` and `(X-r,k)`, let `chi_0,chi_1` be membership in a fixed transported old set `O` in intrinsic `zeta` space.  Pointwise on every positive edge density `s`,

`s_OO=chi_0chi_1s`,
`s_ON=[chi_0(1-chi_1)+(1-chi_0)chi_1]s`,
`s_NN=(1-chi_0)(1-chi_1)s`,

with `s_OO+s_ON+s_NN=s` exactly.  This partition is performed after Moyal, so it has no field-decomposition cross term.  Under `L->ML`, `X->MX`, `r->Mr`, `k->M^-T k`, both endpoint `zeta` labels are individually invariant.  Common affine motion therefore cannot create or erase material-interface ownership.  Classwise edge service also obeys the local endpoint-energy bound `|e^(-ik.r)A_1-A_0|^2<=2(|A_0|^2+|A_1|^2)`.

Dedicated ownership run `31293279918` on `0cd4d89` passed `506` tests plus `50,000` material-edge states: worst affine two-endpoint residual `2.681e-15`, worst positive partition residual `1.705e-13`, zero orientation failures and zero affine membership failures.  Stored artifact: `recorded-results/31293279918/`.

First-hit causality now supplies the missing heat-valid old--old history on every **already supplied signed-good material epoch**.  If `K(t)=int||S_V||_op` and `tau` is the first high-strain contact, then `K(t)<=1/30` on the whole closed pre-hit interval and therefore a reused material frequency grows by at most `exp(1/30)`.  For one deterministic band `|xi|<=M`, the NS heat law gives on `T(N)=cN^-2`

`S_heat(f)=N^3 int dt E_H||delta_r f||_2^2 <= c M^2 E/N`.

This estimate is applied shellwise **before** material summation.  Since OO ownership is imposed only after the positive heat/Moyal law exists, `S_OO,j<=S_heat,j` by positivity; no field `u_old` and no old/new quadratic cross term are introduced.  Along a supplied signed-good material epoch `N_(q+1)/N_q>8/5`, so the physical coefficient contracts by

`rho_phys <= (5/8) exp(1/15) = 0.668086941092 < 441/640 < 7/10`.

Thus the clean whole-pool capacity obeys

`C_OO(q) <= c alpha^2 N_0 P E_global (441/640)^q`,

`sum_(q>=0) C_OO(q) <= (640/199) C_OO(0) < (13/4) C_OO(0)`.

Every first high-strain contact already has `S_heat>=S_*(c)>0`.  Hence after the first material age `q_*` with `C_OO(q_*)<=(1-f)S_*`, exact positive ownership forces

`S_ON+S_NN >= f S_*`.

For `f=1/2`, at least half of every sufficiently old high-strain heat event on that supplied epoch is genuine interface-or-new service.  The theorem status is `EXACT_FIRST_HIT_HEAT_RESERVOIR_EROSION__CLEAN_OO_RATIO_441_640__FINITE_NON_OO_FORCE_ON_SUPPLIED_SIGNED_GOOD_EPOCH`.

Dedicated run `31298235481` on `a0cf825` passed `513` tests plus `50,000` first-hit heat-reservoir states; maximum sampled physical one-step ratio was `0.667120651`, the minimum margin to the clean `441/640` envelope was `2.194e-02`, the minimum first-hit Kelvin margin to `21/20` was `1.611e-02`, and the minimum forced ON+NN margin was `1.131e-03`.  Full physical-energy integration `31298235456` on the same SHA passed `513` tests and the complete causal stack; high-strain half-law margin remained `2.721e-01`, heat retained-fraction margin `6.353e-06`, first-hit monitor/unit failures stayed zero, physical log-product margin remained `1.665`, and master `20,000` traces had worst margin `0`.  Stored dedicated artifact: `recorded-results/31298235481/`.

The same heat geometry actually removes the mixed ON class as a separate long-lived routing problem.  The exact heat law is disintegrated by deterministic Fourier shell **before** Moyal, and translations preserve that event-fiber mark:

`delta_r P_j=P_j delta_r`.

For capacity routing, specialize the arbitrary-set ownership theorem to the canonical old reservoir pool already used by Moyal/reservoir erosion: old intrinsic-zeta cells are addressed inside the exact orthogonal shell carrying their current reservoir energy.  The shell mark is eventwise only; intrinsic `zeta` remains the transported material identity.  Then pointwise on the positive heat law

`old-incident = OO+ON`,

and any edge touching old material lies in the current transported old-shell envelope, hence

`S_OO+S_ON <= S_old-shell <= c M_old^2 P E_global/N`.

Thus OO and ON **together** inherit the same first-hit/signed-good contraction

`rho_phys <= (5/8)exp(1/15)=0.668086941092 < 441/640 < 7/10`.

There is no need to manufacture a temporal-relink charge merely because a same-time heat edge is mixed.  Genuine temporal ownership changes remain relink events, but spatial ON heat service is already old-frequency incident and erodes with the old pool.

This combines sharply with the critical-ancestor restriction of the **same positive heat law**.  Put

`g=(1/2)e^(-1/32)=0.484616617238`,

so the critical shell-time set `G` carries at least `g S_heat`.  Choose `epsilon=g/2=0.242308308619`.  After finite material age, geometric old-incident erosion gives `C_old-incident<=epsilon S_*(c)`.  Since `S_heat>=S_*`,

`S_NN >= (1-epsilon)S_heat >=0.757691691381 S_heat`.

Positive-measure inclusion--exclusion, with no independence assumption, then gives

`S_(NN intersect G) >= (g-epsilon)S_heat = (1/4)e^(-1/32) S_heat =0.242308308619 S_heat`.

So every sufficiently old high-strain event on a supplied signed-good band-addressed material epoch carries a fixed positive **epoch-new NN coherent heat sublaw simultaneously marked by a critical resolved-shell ancestor**.  No packet argmax and no per-cell mass floor enter.  The theorem status is `EXACT_OLD_INCIDENT_HEAT_EROSION__OO_PLUS_ON_DECAY__POSITIVE_NN_CRITICAL_HEAT_SEED_ON_SUPPLIED_SIGNED_GOOD_EPOCH`.

Dedicated run `31299237508` on `8a98d43` passed `521` tests plus `50,000` old-incident shell/ownership/epoch states: worst shell-increment commutator residual `5.587e-16`, worst OO/ON/NN partition residual `8.039e-16`, worst old-incident identity residual `4.321e-16`, maximum sampled physical epoch ratio `0.666333530`, minimum margin to `441/640` `2.273e-02`, minimum forced NN-fraction margin `2.727e-03`, and minimum NN/critical-intersection margin `3.075e-02`.  Full integration `31299237482` on the same SHA passed `521` tests and the complete causal stack; high-strain half-law margin remained `2.721e-01`, retained heat fraction margin `6.353e-06`, physical log-product margin `1.665`, first-hit monitor/unit failures stayed zero, and master `20,000` traces had worst margin `0`.  Stored artifact: `recorded-results/31299237508/`.

The positive `NN intersect G` heat sublaw now supplies an **eventwise smooth carrier entrance without selecting a packet**.  Normalize that physical heat sublaw itself and push it forward by its already exact shell-time mark `(j,t)`.  No heat atom or coherent cell needs a uniform mass lower bound.  On `G`,

`M_j ||P_j u(t)||_2^2 >= mu_* = 32 pi^2/(75 c^2)`.

Set the renewed frequency

`A=3M_j/4`.

Then the whole hard shell `{M_j/2<|xi|<=M_j}` becomes `{2A/3<|xi|<=4A/3}` and carries

`A||P_j u(t)||_2^2 >= (3/4)mu_* = 8pi^2/(25c^2)`.

Choose a smooth scalar Fourier envelope `Q_A` equal to one on that entire shell, with support lower edge `>=3A/5` and upper edge `<=3A/2`.  The canonical terminal dual direction is the shell's own normalized state

`psi=P_j u(t)/||P_j u(t)||_2`,

so exactly

`<Q_Au(t),psi>=||P_j u(t)||_2`.

Thus every atom of the NN-critical heat law canonically determines a lower-scale **smooth whole-shell carrier seed with a scale-independent critical coefficient**.  There is no Gaussian/cell argmax and no packet synthesis.  If the renewed `A`-scale strain stays below `1/30`, the smooth envelope lower edge remains `(3/5)e^(-1/30)A>A/2`, while `S_(A/4)u tensor S_(A/4)u` lies below `A/2`, so the outer-role low--low moat is already strict.  Since `M_j<=N/4`, `A<=3N/16` and the natural-lifetime separation is

`T_A/T_N >= 256/9`.

The material statement is deliberately not overstated.  The NN heat edge is built from `P_jV`, whereas the critical mass is for `P_ju`; the strict low-pass was only used as a contraction and has no inverse lower bound.  Therefore each seed retains **two simultaneous distinct marks**: NN coherent heat provenance and critical `u`-shell energy.  The theorem does not declare the whole `u` shell NN material.  Its status is `EXACT_NN_CRITICAL_HEAT_LAW_TO_SMOOTH_SHELL_CARRIER_SEEDS__NO_CELL_MASS_FLOOR__TEMPORAL_RENEWAL_REMAINS`.

Dedicated run `31300227437` on `2f64de3` passed `528` tests plus `50,000` carrier-seed states: minimum renewed critical-mass margin `6.821e-07`, minimum renewed/child natural-lifetime ratio `28.444444444=256/9`, minimum low--low support gap `8.033e-02`, worst shell-probe registration residual `6.636e-16`, and worst seed-law probability residual `6.661e-16`.  Full integration `31300227438` on the same SHA passed `528` tests and the complete causal stack; high-strain half-law margin remained `2.721e-01`, heat retained-fraction margin `6.353e-06`, NN/critical intersection margin `3.075e-02`, physical log-product margin `1.665`, first-hit monitor/unit failures stayed zero, and master `20,000` traces had worst margin `0`.  Stored artifact: `recorded-results/31300227437/`.

This still does **not** make `D_V` a finite reset or a new-new heat edge into a selected transfer parent.  The high-strain route now has simultaneous frequency ancestry, spatial coherent service, exact material endpoint ownership, first-hit erosion of the entire old-incident `OO+ON` service, a fixed positive `NN intersect critical` heat seed, and a canonical **positive law of lower-scale smooth whole-shell carrier seeds with critical coefficients and a strict support moat**, conditional on a supplied signed-good band-addressed material epoch.  Packet/cell selection is no longer the high-strain renewal obstruction.  The carrier-seed theorem itself stops at the event time and deliberately keeps the NN heat witness distinct from the whole `u`-shell critical coefficient; the following temporal theorem supplies the natural-window continuation-or-first-stop step without collapsing those two marks.

The temporal half of that renewed high-strain entrance is now exact.  Start from any atom of the certified `NN intersect critical` heat law and its smooth renewed carrier `w=Q_Au` at event time `t`, with

`A|z(t)|^2 >= 8pi^2/(25c^2)`.

Inspect backward in elapsed time `ell=t-s` over the natural interval `[max(0,t-cA^-2),t]`.  Four native monitors are used without scalarizing their units: renewed strain action, classified nonaffine role-interface coefficient impulse, HH coefficient impulse, and the minimum intrinsic-material boundary distance of the two retained NN heat-edge endpoint witnesses.  Exact simultaneous debuts remain one joint cause set.  While strain stays below `1/30`, the transported smooth envelope remains strictly above the `A/2` low--low output.  The exact outer-role/adjoint identity is

`z(t)=z(s)+I_HH[s,t]+I_interface[s,t]`.

Therefore before the natural endpoint either a named strain/interface/HH/material-boundary stop occurs, `t=0` is reached and absorbs the branch, or the whole natural corridor survives with

`|z(s)|>=|z(t)|/4`,

hence

`A|z(s)|^2 >= pi^2/(50c^2)`.

The material conclusion is deliberately witness-level: the two NN coherent heat endpoints remain outside the old pool until first boundary contact, but the theorem does not declare the whole `u`-shell carrier NN.  A large interface coefficient impulse remains interface/strain provenance rather than physical work, and a large HH impulse remains an earlier-generation obstruction until it passes through the physical-energy gate.

Dedicated GitHub Actions run `31301650158` on `a2f2a3b` passed `534` tests plus `50,000` temporal corridors: minimum clean survivor critical-mass margin `2.049e-02`, minimum quarter-coefficient margin `6.389e-03`, minimum low--low gap `8.032966e-02`, worst exact Duhamel residual `5.184e-15`, and zero monitor-order/unit-rescaling failures.  Full integration `31301745046` on the same exact SHA also passed `534` tests and the complete causal stack.  Stored dedicated artifact: `recorded-results/31301650158/`.  The theorem status is `EXACT_NN_SEED_BACKWARD_NATURAL_FIRST_STOP_CORRIDOR__ONE_QUARTER_CRITICAL_SURVIVOR_OR_NAMED_STOP__MATERIAL_ATTACHMENT_REMAINS`.

The apparent **material/efficiency attachment** obstruction was stronger than the PDE actually needs.  A full no-hit critical annular carrier already manufactures a fresh physical service law at its own renewed scale; material ownership should be read from that law, not inherited as a whole-carrier label from a child-scale NN witness.

Indeed every prefix of a full no-hit natural corridor satisfies the same strict impulse faces, hence `|z(s)|>|z(t)|/4`.  For the registered affine/Kelvin/viscous analysis dual, not the full nonaffine outer-role adjoint, the existing inverse-heat/strain estimate gives

`||psi(s)||_2 <= J ||psi(t)||_2`,

`J=exp(1/30 + nu c [(3/2)exp(1/30)]^2)`.

The full nonaffine transport mismatch is exactly the already monitored classified interface residual; it is not hidden in this norm bound.  Therefore throughout the whole natural slab

`A||Q_Au(s)||_2^2 >= pi^2/(50 c^2 J^2)`.

The transported role remains in the annulus `|xi|>=(3/5)exp(-1/30)A`.  Its intrinsic `A`-scale heat defect therefore obeys

`E_HA ||delta_r Q_Au||_2^2 >= q_* ||Q_Au||_2^2`,

`q_*=2[1-exp(-( (3/5)exp(-1/30) )^2/2)]`.

The Gaussian displacement itself is unbounded and is **not** inserted directly into the bounded-displacement reservoir theorem.  Instead `||delta_r w||_2^2<=4||w||_2^2`; exact 3D Gaussian-tail arithmetic with Arb proves that deleting `|r|>3/A` leaves strictly more than half of the annular heat lower.  Consequently at every point of the full natural slab some actual displacement satisfies

`|r_s|<=3/A`,

`A||delta_(r_s) Q_Au(s)||_2^2 >= Y_0(c,nu)>0`,

with `Y_0=(q_*/2) pi^2/(50 c^2 J^2)`, and the normalized integrated bounded heat service is at least `cY_0`.  For the default certificate `c=nu=1`, `Y_0=2.33125086914e-4`.  This threshold is scale-independent in `A` for fixed canonical block parameters; it is **not** promoted to a universal globally bounded reset constant.

Only after this renewed positive service exists is exact Moyal applied and its two intrinsic endpoints partitioned into OO/ON/NN.  The old NN heat witness remains legitimate seed provenance but is not used to declare the whole smooth carrier new material.  Thus the survivor supplies exactly the uniform bounded-displacement service threshold consumed by the existing old-pool/coherent-service stopping theorem.  This is a coherent-service entrance, not a theorem that the whole renewed carrier is near-extremal HH-transfer efficient.

Dedicated run `31303081994` on `f0ac683` passed `542` tests and `50,000` annular-service states.  Arb certified `ARB_CERTIFIED_RADIUS_3_RETAINS_MORE_THAN_HALF_ANNULAR_HEAT_SERVICE`; the minimum spectral heat margin was `5.527e-03`, the radius-3 truncation margin `3.781e-02`, and the worst renewed OO/ON/NN residual `1.137e-13`.  Full integration `31303081937` on the same SHA also passed.  Stored artifact: `recorded-results/31303081994/`.  Status: `EXACT_FULL_NATURAL_CRITICAL_CARRIER_TO_BOUNDED_ANNULAR_SERVICE_THRESHOLD__MATERIALITY_REREAD_FROM_SERVICE__UNIVERSAL_RENEWAL_REMAINS`.

This observation shortens the high-strain route further.  The child-scale NN heat law is not needed at all for **renewal entrance**.  The resolved-ancestor theorem already supplies the positive diagnostic sampling law `D_V|_G`, where

`G={(j,t): M_j||P_j u(t)||_2^2>=mu_*}`,

`mu_*=32pi^2/(75c^2)`,

and `D_V(G)>=D_V/2`.  Normalizing this retained dissipation mass is only a shell-time sampling device; it is **not** a causal probability and never replaces actual positive HH child-energy work in Shannon/Renyi or transfer causality.

Push each deterministic shell-time mark directly to `A=3M_j/4`.  Then

`A||P_j u(t)||_2^2 >= 8pi^2/(25c^2)`,

with exact whole-shell registration into a smooth `Q_A=1` envelope and natural-lifetime ratio at least `256/9`.  Before any renewed material service is assigned, there are exactly three native backward first-stop monitors: renewed strain, classified role-interface coefficient impulse, and HH-regeneration coefficient impulse.  No synthetic material-boundary monitor is inserted.  A first hit keeps its already named recursive owner; `t=0` absorbs; otherwise the full natural survivor has the quarter coefficient and, by the annular-service theorem above, manufactures its own bounded `A`-scale service.  Material ownership is then reread from that service by exact Moyal.

Hence the shortest high-strain entrance is now

`high strain -> D_V|_G -> critical shell carrier -> named first stop / t=0 / own-scale coherent service`,

not the longer mandatory route through child-scale heat ownership, old-incident erosion and `NN intersect G`.  The latter chain remains fully valid and useful as a finer material-capacity/provenance theorem, but it is no longer a prerequisite for renewal.  `D_V` remains scale-critical recursive currency and is never promoted to a finite reset.

The first short-route run `31303327206` on `51840e6` had all `547` tests pass and failed only in the stress harness because an absolute timestamp tolerance was used after independent unit rescaling while `T_A=cA^-2` can be arbitrarily large in the sampled range.  Commit `7ef566a` changed only that floating regression tolerance to a corridor-scale tolerance.  Corrected dedicated run `31303385148` on exact `7ef566a` passed `547` tests plus `50,000` states: minimum renewed-seed critical-mass margin `5.808e-05`, minimum full-survivor critical-mass margin `2.909e-01`, worst exact Duhamel residual `4.974e-14`, and zero monitor-order/unit failures.  Full integration `31303385157` on the same exact SHA passed the entire causal stack; physical log-product margin remained `1.665`, Shannon/Renyi lower margin `1.958e-01`, and master `20,000` traces had worst margin `0`.  Stored artifact: `recorded-results/31303385148/`.  Status: `EXACT_HIGH_STRAIN_CRITICAL_DISSIPATION_LAW_TO_OWN_SCALE_SERVICE_REENTRY__NN_NOT_REQUIRED_FOR_RENEWAL_ENTRANCE__UNIVERSAL_SOURCE_RELINK_REMAINS`.

Thus **the high-strain material-attachment seam is closed at the level actually required by the existing coherent-service stopping architecture**.  What remains is no longer “attach the old NN witness to the whole renewed carrier”.

The next simplification is now certified at the **critical-shell** level, which is the smaller and more physical abstraction.  Assume at one actual shell-time event only

`M||P_M u(t)||_2^2 >= mu_0 > 0`.

No material label, packet, coherent-cell selector or normalized sampling law belongs to the core theorem.  Put `A=3M/4` and choose the existing smooth scalar `Q_A=1` on the whole hard shell.  With the shell's normalized state as terminal analysis probe,

`A|z(t)|^2 >= (3/4)mu_0`.

Before materiality is assigned, inspect exactly the same three native monitors as in the short high-strain corridor: renewed strain, classified role-interface coefficient impulse and HH-regeneration coefficient impulse.  The generic interface is stricter than the previous caller convention: it records the maximum backward horizon actually observed and refuses to certify either `t=0` or a full-natural survivor unless the monitor path covers the entire required interval.  A named hit keeps its existing recursive owner; `t=0` absorbs; otherwise every prefix satisfies the quarter coefficient bound and hence

`A|z(s)|^2 >= (3/64)mu_0`.

With the already-certified registered affine/Kelvin/viscous analysis-dual cost `J`,

`A||Q_Au(s)||_2^2 >= 3mu_0/(64J^2)`

throughout the whole no-hit slab.  The same annular heat law and Arb-certified radius-three truncation therefore give an actual bounded displacement

`|r_s|<=3/A`,

`A||delta_(r_s)Q_Au(s)||_2^2 >= Y_shell(mu_0) := q_b 3mu_0/(64J^2)>0`,

and normalized full-slab bounded service at least `c Y_shell(mu_0)`.  Only after this positive service exists is exact Moyal used to read OO/ON/NN from the actual endpoints.  Thus the theorem is material-free before service and material-exact after service.

Two previously distinct branches become suppliers of this same shell law.  First, any already-certified resolved-dissipation lower `D_V>=D_0>0` gives `mu_0=D_0/c`: the dyadic low-mass portion costs exactly `D_0/2`, hence at least `D_V-D_0/2>=D_V/2` of the **actual** resolved-dissipation mass lies on qualifying shells.  Normalizing that restriction remains diagnostic sampling only and never enters Shannon/Renyi or substitutes for positive HH child-energy work.  Second, on the existing dominant fresh coherent-service branch, if fresh service is at least `Y/4` and one edge carries fraction `theta`, then local capacity gives `M(E_C+E_{C-r})>=theta Y/8`; because each cell energy is bounded by the whole shell energy,

`M E_shell >= theta Y/16`,

so the canonical `theta=1/4` supplies the clean whole-shell lower `Y/64`.  This is not a per-cell floor and does not freeze a packet.

The high-strain theorem is an **exact specialization**, not a parallel estimate: setting `mu_0=32pi^2/(75c^2)` gives terminal mass `8pi^2/(25c^2)`, quarter-survivor mass `pi^2/(50c^2)`, and exactly the previously certified carrier/service constants.  The theorem is deliberately shell-local; supplier-specific signed-good progress relative to a previous block scale remains a separate geometric fact and is not fabricated here.

The first branch SHA `ac1b162` produced `554` passing tests and one fixture-only failure because the scope regression expected the literal `pressure/source` while the certificate said `source/pressure`; no theorem code or constant failed.  Commit `ccdf9f3` corrected only that literal.  Dedicated run `31304867746` on exact `ccdf9f3` then passed `555` tests plus `50,000` shell/supplier/corridor states: high-strain specialization signed margin exactly `0`, minimum generic `D_0` retained fraction `0.500135189`, fresh-cluster whole-shell identity margin exactly `0`, minimum full-survivor coefficient-mass margin `1.795e-05`, worst survivor Duhamel residual `4.885e-15`, zero monitor-order/unit failures, zero incomplete-horizon certification failures, and sampled joint ties of size `2`.  Full integration `31304886428` on the same exact SHA passed the complete causal stack; physical log-product margin remained `1.665`.  Stored artifact: `recorded-results/31304867746/`.  Status: `EXACT_CRITICAL_SHELL_MASS_TO_OWN_SCALE_SERVICE_REENTRY__GENERIC_DV_AND_FRESH_CLUSTER_SUPPLIERS__MATERIALITY_DEFERRED`.

This materially shrinks the source/SGS frontier.  Any source route which already supplies a **resolved** `D_V>=D_0` now has a direct continuum entrance to critical shell -> named stop / `t=0` / own-scale service; the viscous source lower is one such supplier, and the quadratic objective-strain route supplies another after integration.  A dominant fresh coherent cluster also re-enters through the same shell theorem.  By contrast, pressure/source mass occupation and generic high-frequency enstrophy are not silently renamed `D_V`; they stay on their certified reservoir/service/entropy routes until a correct unit-matched supplier theorem is available.

The proposed material-label quotient is now theorem-level.  For the **same smooth role** `Q` and the **same registered analysis probe** `psi`, the exact carrier identity remains

`z(t)=z(s)+I_HH+I_interface`.

No old-pool indicator, intrinsic-cell address, or selected coherent-family characteristic occurs in this equation.  Hence a pure material membership/address change creates no second carrier impulse: the carrier first-stop set has only the two physical faces `|I_interface|>=|z(t)|/4` and `|I_HH|>=|z(t)|/2`.  If neither is hit, `|z(s)|>=|z(t)|/4` even when material sidecars change.  Nonaffine motion is not made free: its effect is already exactly the classified interface term and must not be charged a second time merely because a bookkeeping boundary was crossed.

Material/service accounting remains exact.  Re-reading OO/ON/NN with the same positive edge weights preserves total service identically.  If the selected coherent family changes, the full Moyal symmetric-difference charge remains

`|E(S_new)-E(S_old)| <= R_switch = E(S_new symmetric_difference S_old)`,

and continues to belong to ancestry/service routing; paying that sidecar charge simply does not force destruction of an unchanged PDE carrier.  If `Q` or `psi` really changes, the quotient refuses the continuation and delegates to the existing event-role/nonaffine-interface/Xi/relink registration.  Thus the physical split is now explicit:

`material/service sidecar change != carrier-role change`.

The older common-slice theorem remains a correct conservative superset; only its same-`Q`, same-`psi` material-sidecar subtype is refined.  Dedicated run `31306299773` on exact `35846d9` passed the complete pytest dependency stack and `50,000` quotient states: minimum transparent quarter-coefficient margin `8.096e-04`, worst same-role Duhamel residual `8.564e-14`, worst service reclassification residual `1.101e-13`, selected-family switch-jump margin `-1.421e-14` at floating precision, zero label-transparency/switch-charge/role-delegate failures, and maximum carrier-stop count exactly `2`.  Full physical-energy integration `31306305152` on the same SHA also passed.  Stored artifact: `recorded-results/31306299773/`.  Status: `EXACT_MATERIAL_LABEL_SIDECAR_QUOTIENT__PURE_LABEL_AND_SELECTED_FAMILY_SWITCH_ADD_NO_SECOND_CARRIER_IMPULSE__ROLE_CHANGE_REMAINS_INTERFACE`.

This removes another observer-imposed obstruction from the shortest continuum architecture.  The next source audit shows that the remaining coherent objective-source branch is itself an **owner compiler**, not a missing packet theorem.

For the exact coherent corotational strain equation, group source action only by its physical owner:

`local_DV + pressure + SGS + viscosity`.

If `A_obj` is the objective-variation action on scaled lifetime `c`, the positive owner weights satisfy

`Sigma_local + Sigma_P + Sigma_R + Sigma_nu >= A_obj/c`,

so at least one owner carries `A_obj/(4c)`; exact ties remain one joint owner set.  No scalarized source clock or theorem-name priority is introduced.

The local coherent quadratic/Reynolds terms already obey `Sigma_local<=C_local D_V`, while the filtered viscous source satisfies `rho_nu<=nu sqrt(d_V)/1500`.  Thus

`D_V >= Sigma_local/C_local`

or, on the viscous owner,

`D_V >= (1500 Sigma_nu/nu)^2/c`.

Both feed the generic critical-shell theorem as scale-critical recursion, never as additive resets.

The differentiated SGS owner simplifies further.  The resolved source collision, Germano inequality and coherent square-service threshold have powers

`rho_R -> ||R||_(3/2) -> Q^(3/2) -> Y^(2/3)`,

so the exponents cancel exactly and give the direct linear law

`Y >= C_Y rho_R`,

`C_Y = 380/[g1(1+g1)(C_LP C_B)^2]`.

Hence integrated objective-SGS source creates integrated coherent square service with **no persistence hypothesis and no affine-radius packet**.  Its existing positive service law routes to high-frequency dissipation, old-pool capacity, selected-interface `Xi`, a dominant fresh critical shell, ancestry entropy, or same-ancestry cycle.  The dominant quarter-service edge supplies pair critical-mass occupation `>=Y/32`, whole-shell occupation `>=Y/64`, hence a pointwise shell event `>=Y/(64c)` which enters the generic shell re-entry theorem.  High-frequency dissipation is deliberately kept under its own owner and is not renamed resolved `D_V`.

Pressure also retains its native two-way structure.  From

`rho_P <= mu_V/5700 + ||R||_(3/2)/380`,

an integrated pressure owner `Sigma_P` gives

`int mu_V >= 2850 Sigma_P`

or

`int ||R||_(3/2) >= 190 Sigma_P`.

The stress alternative is exactly an effective SGS source weight `Sigma_P/2` and enters the same direct coherent-service route.  The low-pass mass alternative is **not** promoted to a generic critical shell: it remains objective pressure-Hessian reservoir occupation.  A derivative-order audit found that the earlier reuse sentence borrowed the **H1 pressure-third** (`D^3P`) pair coefficient for the **objective pressure-Hessian** (`D^2P`) source.  The numerical conclusion was conservative, but the provenance was wrong and is now corrected.

For one objective pressure-Hessian low--low pair, order-two differentiation plus the band-limited `L^(3/2)->L^infinity` step gives, up to the fixed pair-independent Riesz/Bernstein constant,

`C_ab^(2)(N) = M_max^4 sqrt(M_a M_b)/N^4`

per unit `sqrt(E_a E_b)`.  On the already certified signed-good low-strain lineage, reservoir frequencies grow by `<21/20` while the block scale grows by `>8/5`, so a fixed material pair contracts by

`(21/20)^5(5/8)^4 = 4084101/20971520 < 1/5`.

Consequently even the adversarial future energy-cap sum of one fixed objective-Hessian pair is `<5/4` of its generation-zero capacity.  The separate H1 pressure-third source keeps its own coefficient

`C_ab^(3)(N)=M_max^3 sqrt(M_a M_b)/N^3`

and its previously certified ratio `(21/20)^4(5/8)^3=194481/655360<1/3`, with total future capacity `<3/2`.  Thus **derivative order is pressure provenance**: the material-pair mechanism is common, but the lifetimes are not.  Persistent objective-Hessian service must relink pairs, fragment into pair/component entropy or cycle structure, leave the supplied signed-good low-strain lineage, or use the SGS part.

The first theorem SHA `7f1ff35` already passed all `560` tests and the new `50,000`-state compiler stress, but its dedicated workflow then failed in a **dependency-only fixture** because the subsequent resolved-objective source certificate requires `python-flint` and the new workflow had installed only `requirements.txt`.  No theorem formula or test failed.  Commit `8015fb8` added only `python-flint` to that workflow.  Corrected dedicated run `31318529626` on exact `8015fb8` then passed the full theorem suite and all dependency certificates: worst SGS closed-form homogeneity residual `1.207e-15`, minimum owner-pigeonhole margin `2.066e-03`, pressure split identity margin `-1.110e-16`, local-`D_V` identity margin `-1.110e-16`, viscous Cauchy identity margin `0`, and maximum sampled joint owner count `4`.  Full physical-energy integration `31318546158` on the same SHA passed the complete causal stack and master episode stress.  Stored artifact: `recorded-results/31318529626/`.  Status: `EXACT_COHERENT_OBJECTIVE_SOURCE_OWNER_COMPILER__LOCAL_DV_AND_VISCOSITY_TO_CRITICAL_SHELL__SGS_TO_COHERENT_SERVICE__PRESSURE_TO_SGS_OR_RESERVOIR__NO_PACKET_SYNCHRONIZATION`.

Thus the **coherent objective-source assembly seam is closed without introducing a source packet**.  The high-frequency SGS owner is now theorem-level as well, but only after separating the observer's smooth LP scalar from the PDE's actual dissipative state.

Fix once and for all a smooth square-normalized LP analysis--synthesis frame

`u_j=phi_j(D)u`, `sum_j |phi_j(xi)|^2=1`, `u=sum_j phi_j(D)u_j`,

with high-band support `|xi|>=M_j/2`.  This is an exact Calderon frame, not the false overlapping-band identity `u=sum_j u_j`.  If in general the high analysis support is `|xi|>=aM_j` and `sum|phi_j|^2<=B`, Plancherel gives the exact comparison

`D_tail := N int ||grad P_>N u||_2^2 dt >= c_LP D_high`,

`c_LP=a^2/B`.

For the canonical square-normalized frame, `a=1/2`, `B=1`, hence `c_LP=1/4`.  The smooth-LP `D_high` remains a service observable and is never identified with the orthogonal hard-tail PDE currency `D_tail`.

This distinction is necessary: `D_high` alone has no scale-independent shell-mass floor, since placing all currency at dyadic level `j` with `mu_j=2^-jD_high` preserves `2^jmu_j=D_high` while `mu_j->0`.  The missing information is the exact hard-tail energy law.  For `w=P_>N u` and actual positive nonlinear tail work

`W_>^+=int 2[Re<w,-P_>N P div(u tensor u)>]_+ dt`,

Navier--Stokes gives

`N||w(s)||_2^2 + N W_>^+ >= 2 nu D_tail`.

Therefore, with exact ties retained jointly, either inherited scaled tail energy or positive nonlinear tail work carries at least `nu D_tail` and hence at least `nu c_LP D_high`.

On the inherited branch,

`N||w(s)||_2^2=sum_(j>=1) 2^-j mu_j(s)`

with dyadic weights summing to one, so some actual hard shell `M_j>=2N` has

`M_j||P_j u(s)||_2^2 >= nu D_tail`.

This is a genuine critical-shell seed and enters the generic shell first-stop/service theorem.  On the regeneration branch, orthogonal shell work satisfies `W_>^+<=sum_j W_j^+`, and `M_j/N>=2` gives

`sum_j M_j W_j^+ >= 2N W_>^+`.

At each hard shell, choosing `V=S_(M_j/4)u` excludes low--low work by support.  The signed shell work is exactly HH plus resolved mixed/cross-interface work, hence the positive law is covered by actual positive HH and interface owners; one carries at least `nu D_tail`.  This is **not** yet the generated-energy KL branch: HH work must still pass its child-energy gate, and interface work remains physical interface/strain provenance.

Dedicated run `31320980362` on exact `1dcc5ced307d8fe5db450cabbf0c06f5b146b4ec` passed the full pytest stack and `50,000` high-tail/shell/work states: minimum hard-tail gradient lower/upper margins `7.821e-08` / `7.351e-09`, minimum two-owner energy margin `2.946e-07`, inherited-shell margin `2.020e-08`, shell-work disintegration margin `5.521e-08`, HH/interface owner margin `1.806e-07`, and maximum joint owner count `2` at both energy and regeneration splits.  Stored artifact: `recorded-results/31320980362/`.  Status: `EXACT_SQUARE_LP_HIGH_SERVICE_TO_PHYSICAL_TAIL_ENERGY__INHERITED_CRITICAL_SHELL_OR_PHYSICAL_REGENERATION__NO_RESOLVED_DV_RELABEL`.

The derivative-correct pressure theorem was validated in the same exact stack.  Objective-source run `31320980346` passed `50,000` source states and the pressure reservoir dependency: maximum sampled objective-Hessian pair ratio `0.188320537`, margin below `1/5` `1.168e-02`, maximum sampled pressure-third ratio `0.288962935`, margin below `1/3` `4.437e-02`, with both amortized clean-envelope margins nonnegative at floating precision.  Stored artifact: `recorded-results/31320980346/`.  Full physical-energy integration `31320980342` on the same SHA passed the complete causal stack.  Pressure status: `EXACT_OBJECTIVE_PRESSURE_HESSIAN_AND_PRESSURE_THIRD_PAIR_EROSION_GIVEN_SIGNED_GOOD_LOW_STRAIN_LINEAGE`.

The continuum frontier is therefore narrower again.  The anonymous high-frequency-enstrophy exit is gone: it now becomes inherited critical shell or actual positive regeneration/interface work.  The coarse pressure-reservoir entrance has now disappeared as well.  The resolved objective pressure Hessian has an exact **pair-native** renewal law before any material or coherent sidecar is attached.

Let

`H=-N^-4 <Hess P>_gamma`

be the actual coherent-averaged pressure tensor in the corotational strain equation and choose its measurable Frobenius dual `Z=H/||H||_F`, with `Z=0` at `H=0`.  For the exact filtered split `H=H_V+H_R`, decompose only the strict transporter `V=S_(N/4)u` into hard orthogonal dyadic shells at this physical pressure event.  The resolved bilinear term reconstructs exactly as **unordered** shell pairs `{a,b}`: diagonal once, off diagonal as both orientations together.  Thus observer orientation is quotiented before any entropy is read.  Identical unordered pairs are integrated once, and one hard shell label carries one physical frequency; a proof representation is not allowed to manufacture entropy by duplicating one physical atom.

The elliptic map from tensor stress directly to the pressure Hessian has Frobenius multiplier norm at most `|xi|^2`; no hidden `L^(3/2)` Riesz constant is needed.  Combining the already Arb-certified order-two tensor bound `<1/380` with the resolved shell `L^3` mass bound `<1/15` gives the exact clean ordered coefficient

`256/1425 < 1/5`.

Hence one unordered physical pair obeys

`|p_ab| <= (kappa_ab/5)(M_max/N)^4 sqrt(mu_a mu_b)`,

with `kappa_aa=1` and `kappa_ab=2` off diagonal, where `mu_j=M_j||P_jV||_2^2`.  Gaussian averaging is a probability contraction and adds no extra factor.  The countable canonical dyadic law is absolutely source-summable: for `M_j=(N/4)2^-j`,

`sum_(a<=b) cap_ab <= N||V||_2^2/2560`.

Aggregate resolved energy appears here only as a convergence budget, **not** as the pressure renewal state.

Scalarizing by the same event dual gives the exact positive-source cover

`rho_P <= [r_SGS]_+ + sum_(a<=b)[p_ab]_+`.

For integrated pressure source weight `Sigma_P`, either the positive SGS source or the positive resolved pair law carries at least `Sigma_P/2`; exact ties remain joint.  The SGS owner gives `int||R||_(3/2)>=190 Sigma_P` and enters the existing coherent-service route.  On the resolved owner, normalize the actual positive unordered pair law by `q_ab=R_ab/R_pair` and write

`H_2^P=-log sum q_ab^2`.

A countable positive probability law has an attained maximal atom, and `q_max>=sum q^2=exp(-H_2^P)`.  The maximal actual pair therefore exposes, by its integrated capacity, a resolved hard shell with

`mu_child >= 320 exp(-H_2^P) Sigma_P/c`.

The strict low-pass is a scalar `L^2` contraction commuting with the hard shell projector, so

`M||P_MV||_2^2 <= M||P_Mu||_2^2`;

there is no inverse low-pass estimate and no equality claim.  Consequently the native pressure law is the threshold-free conjugacy

`mu_child exp(H_2^P) >= 320 Sigma_P/c`.

Thus **every resolved pressure-pair owner already enters the generic critical-shell first-stop/service theorem**.  The old quarter split is only a diagnostic corollary: `q_max>=1/4` gives `mu_child>=80 Sigma_P/c`, while `q_max<=1/4` gives `H_2^P>=log 4`; equality gives both, but these are not competing causal fates.  Pressure-source entropy is not a child-energy Shannon/Renyi probability and is not independently declared a transfer-cost stop.  It only measures the weakening of the same shell seed.

Composing with the generic shell theorem makes the physical meaning even clearer.  On a full no-hit renewed shell corridor, if `q_b` is the certified bounded annular heat fraction and `J` the registered affine/Kelvin/viscous dual-growth factor, then

`Y_shell exp(H_2^P) >= 15 q_b Sigma_P/(c J^2)`,

and the full-natural integrated own-scale service satisfies

`S_shell exp(H_2^P) >= 15 q_b Sigma_P/J^2`.

So the pressure owner is naturally a **service-versus-complexity relation**, not a reservoir packet.  Material/coherent labels may be attached only after the hard event as sidecars; on a supplied signed-good low-strain material lineage, the derivative-correct fixed-pair `<1/5` erosion theorem remains an optional reuse refinement rather than a renewal entrance.  The coarse inequality `rho_P<=mu_V/5700+||R||_(3/2)/380` remains mathematically valid as a diagnostic only.

Dedicated GitHub Actions run `31342227267` on exact `8bdfc89bbfe38ba87a50d3c30b9e945567fa2b2c` passed `601` tests plus `50,000` pressure tensor/pair/source states: worst dual scalar residual `1.232e-14`, positive-cover margin `-1.066e-14` at floating precision, worst unordered reconstruction residual `1.540e-14`, minimum sampled pair-capacity margin `8.889e-19`, minimum owner-half-split margin `3.091e-04`, dominant-shell clean margin `3.476e-02`, diffuse quarter-corollary entropy margin `9.525e-02`, entropy-shell tradeoff margin `8.585e-03`, and maximum joint primary owner count `2`.  Full causal integration `31342227257` on the same exact SHA passed the complete stack and master episode stress.  Stored artifact: `recorded-results/31342227267/`.  Status: `EXACT_OBJECTIVE_PRESSURE_HESSIAN_DUAL_PAIR_ATOMIZATION__PAIR_OWNER_ALWAYS_TO_CRITICAL_SHELL_WITH_ENTROPY_TRADEOFF__QUARTER_SPLIT_ONLY_A_COROLLARY__AGGREGATE_MU_V_NOT_CANONICAL`.

This closes the **local objective-pressure renewal entrance** without aggregate `mu_V`, packet synchronization, or a mandatory pressure-reservoir lineage.  The compiler migration is now theorem-level as well.  The old function implementing

`rho_P<=mu_V/5700+||R||_(3/2)/380`

is retained only as an explicit `DIAGNOSTIC_ONLY` compatibility interface.  It no longer appears in the canonical master-facing pressure state.

The canonical symbolic pressure interface now records only the actual positive-source split

`SGS pressure source >= Sigma_P/2  OR  resolved unordered pair law >= Sigma_P/2`,

with exact ties joint.  A realized wrapper consumes the **actual** positive SGS/pair source weights from the Frobenius-dual pressure theorem.  If SGS owns the event with actual integrated positive weight `R_SGS`, it uses the exact source bound

`int||R||_(3/2) d tau >= 380 R_SGS`

and enters the same direct coherent-service map as the objective SGS owner.  If the resolved pair law owns the event, it imports the certified

`mu_child exp(H_2^P) >= 320 Sigma_P/c`

and enters the generic critical-shell first-stop theorem.  No `mu_V` reservoir state is manufactured on the way.

Because the generic shell bounded-service theorem is linear in its shell-mass input, this conjugacy survives composition.  On a **full no-hit natural shell survivor** the compiler records

`exp(H_2^P) Y_shell >= Y_shell^generic(320 Sigma_P/c)`

and

`exp(H_2^P) S_shell >= S_shell^generic(320 Sigma_P/c)`.

These are deliberately conditional full-survivor statements: the compiler does not promote service before the shell first-stop corridor is actually resolved.  Pressure `H_2^P` remains a deterministic source-complexity coordinate, not a child-energy causal probability and not a new stop class.  The fixed-material-pair `<1/5` erosion theorem remains an optional post-sidecar reuse refinement only.

The first compiler SHA `e667f37` failed only because the new wrapper passed keyword `pair_dominance_fraction` to the already-certified pressure theorem whose public API names that argument `dominant_fraction`; `602` tests passed and exactly three wrapper fixtures failed with the same `TypeError`.  Commit `e3b1478cda467ea0546c61465bf49185d7daeacd` changed that single wiring keyword.  Dedicated compiler run `31343542038` on exact `e3b1478` then passed `605` tests, `50,000` compiler states, and direct `50,000`-state dependencies for high-tail reentry, critical-shell service, objective pressure pairs, pressure reuse and material-label quotient.  The compiler stress reported worst SGS closed-form residual `1.215e-15`, minimum owner-pigeonhole margin `2.748e-03`, pressure diagnostic split margin `-1.110e-16`, pressure entropy-shell margin `-4.547e-13`, pressure full-survivor service registration margin `-4.441e-16`, maximum joint source-owner count `4`, and maximum sampled joint pressure-owner count `1`; the exact half-tie is covered separately by the deterministic fixture.  Stored artifact: `recorded-results/31343542038/`.  High-frequency integration workflow `31343542019` and full physical-energy causal integration `31343542021` both passed on the same exact SHA; the latter also passed the branch compiler and master episode stress.

Compiler status: `EXACT_COHERENT_OBJECTIVE_SOURCE_OWNER_COMPILER__LOCAL_DV_AND_VISCOSITY_TO_CRITICAL_SHELL__SGS_TO_COHERENT_SERVICE__PRESSURE_TO_SGS_OR_ENTROPY_WEIGHTED_CRITICAL_SHELL__AGGREGATE_MUV_DIAGNOSTIC_ONLY__NO_PACKET_SYNCHRONIZATION`.

This removes the last **coarse pressure-reservoir state** from the canonical objective-source compiler.  A new simplification is now visible on the fresh SGS side, but it must respect a crucial representation distinction.  Coherent-cell entropy itself depends on the chosen phase-space partition, so the pressure argument cannot simply be copied to cells.  The old material pool, however, is already a measurable intrinsic set `O`, and NN/fresh edge ownership is pointwise on the positive service measure.  Therefore one can first quotient away coherent-cell refinement by summing fresh service **within each canonical LP frequency band**.  If `F_j` is that fresh band service, then partition refinement changes no `F_j`, while

`F_j <= M_j||delta_r u_j||_2^2 <= 4 M_j||u_j||_2^2`.

The canonical LP frame is already constructed from a smooth dyadic annular cover and square-normalized, but the current certificate registers only the lower support `|xi|>=M_j/2`.  The next native theorem should fix the harmless missing upper support once and for all, for example

`supp phi_j subset {M_j/2<|xi|<=2M_j}`.

Then `|phi_j|<=1` and one smooth band touches only the two hard annuli at `M_j` and `2M_j`, giving exactly

`M_j||u_j||_2^2 <= mu_(M_j) + (1/2)mu_(2M_j)`

and hence some actual hard shell with mass at least `(2/3)M_j||u_j||_2^2`.  If the fresh band law is normalized and has scale collision entropy `H_2^scale`, a candidate refinement-invariant supplier law is therefore

`mu_hard exp(H_2^scale) >= Y/(24c)`

whenever fresh integrated service is at least `Y/4`.  This would remove the current quarter **coherent-cell** dominance as a renewal entrance while leaving coherent ancestry entropy/cycle available as a finer material refinement.  The upper-annular-frame registration and the exact measurable fresh-band pushforward are the next theorem target; this last formula is a research direction, not yet canonical.

Beyond that, the remaining structural seams are the own-scale continuation of actual positive high-tail regeneration work, supplier-specific signed-good scale geometry, and final continuum master assembly.  No Navier--Stokes global-regularity proof is claimed.

---

### Repository reading order

For a new reader, the recommended order is:

1. `RESEARCH_LEDGER.md` — this document;
2. `docs/single_edge_stability_certificate.md`;
3. `docs/log_scale_flux_bridge.md` and `docs/smooth_log_flux_cocycle.md`;
4. `docs/crossing_moat_extraction.md` and `docs/smooth_sgs_symbol_freezing.md`;
5. `docs/localized_sgs_pressure_ledger.md` and `docs/annular_pressure_collision.md`;
6. `docs/transfer_preserving_profile_extraction.md`, `docs/gaussian_packet_inverse.md` and `docs/packet_inverse_theorem.md`;
7. `docs/affine_gaussian_grain_dynamics.md`, `docs/intrinsic_3d_triad_plane.md`, and `docs/strain_coherence_objective_gradient.md`;
8. `docs/helical_spin_transport.md`, `docs/helical_phase_holonomy.md`, `docs/full_strain_observability.md`, `docs/objective_helical_polarization.md`, `docs/extremal_helicity_symplectic.md`, `docs/helical_frame_lipschitz.md`, `docs/relative_polarization_transport.md`, and `docs/localized_polarization_packet.md`;
9. `docs/affine_gaussian_forcing.md`, `docs/material_phase_lock.md`, `docs/quadratic_swirl_kernel.md`, `docs/affine_polarization_curvature.md`, `docs/affine_shell_aspect.md`, `docs/affine_critical_grain.md`, `docs/affine_window_balance.md`, `docs/sideband_coherence_daughter.md`, `docs/hermite_triad_selection.md`, `docs/curvature_sideband_irrep.md`, `docs/h3_sideband_no_escape.md`, `docs/sideband_pair_ancestry.md`, `docs/h1_swirl_mild_aspect.md`, `docs/h1_swirl_no_escape.md`, `docs/h1_covariant_source.md`, `docs/sgs_source_collision.md`, `docs/onsager_increment_collision.md`, `docs/source_episode_collision.md`, `docs/affine_aspect_sticky.md`, `docs/ancestor_reservoir_sync.md`, `docs/bandlimited_source_sampling.md`, `docs/pressure_reservoir_sync.md`, `docs/reservoir_pool_erosion.md`, `docs/affine_coherent_bessel.md`, `docs/affine_coherent_moyal.md`, `docs/coherent_transfer_cells.md`, `docs/coherent_increment_service.md`, `docs/high_frequency_dissipation_reentry.md`, `docs/coherent_covariance_interface.md`, `docs/coherent_service_stopping.md`, `docs/coherent_sgs_episode.md`, `docs/service_or_flat_rigidity.md`, `docs/physical_flat_episode.md`, `docs/flat_companion_gate.md`, `docs/causal_binary_ancestry.md`, `docs/weighted_causal_reuse.md`, `docs/renyi_causal_reuse.md`, `docs/adjoint_kelvin_duhamel.md`, `docs/asynchronous_duhamel_sync.md`, `docs/physical_energy_causal_bridge.md`, `docs/recursive_coherent_witness_extraction.md`, `docs/resolved_role_egorov.md`, `docs/coherent_affine_projection.md`, `docs/coherent_averaged_strain_source.md`, `docs/coherent_service_or_flat.md`, `docs/complex_young_parent_marking.md`, `docs/dual_gaussian_root_registration.md`, `docs/bargmann_root_cell_registration.md`, `docs/amplitude_entropy_causal_reuse.md`, `docs/common_slice_coefficient_registration.md`, `docs/outer_moving_role_extraction.md`, `docs/nonaffine_role_interface_work.md`, `docs/event_anchored_role_registration.md`, `docs/physical_pair_weighted_productivity.md`, `docs/recursive_physical_witness_constructor.md`, `docs/joint_causal_stop_projection.md`, `docs/smooth_sgs_first_hit_extraction.md`, `docs/smooth_material_carrier_relay.md`, `docs/resolved_cutoff_repartition_relay.md`, `docs/high_strain_resolved_ancestor.md`, `docs/high_strain_heat_increment_service.md`, `docs/heat_edge_material_ownership.md`, `docs/first_hit_heat_reservoir_erosion.md`, `docs/old_incident_heat_erosion.md`, `docs/nn_critical_heat_carrier_seed.md`, `docs/nn_seed_temporal_first_stop.md`, `docs/critical_annular_carrier_service_reentry.md`, `docs/high_strain_critical_carrier_reentry.md`, `docs/critical_shell_service_reentry.md`, `docs/material_label_carrier_quotient.md`, `docs/objective_source_routing_compiler.md`, `docs/objective_pressure_pair_atomization.md`, `docs/material_coherent_labels.md`, `docs/physical_multicurrency_master.md`, `docs/high_strain_dissipation_collision.md`, `docs/resolved_objective_strain_collision.md`, `docs/divfree_coherent_parseval.md`, `docs/coherent_localization_operators.md`, `docs/physical_transfer_defect_moat.md`, and `docs/physical_branch_compiler.md`;
10. `docs/curvature_balanced_moat.md` and `docs/objective_strain_source_collision.md`;
11. `docs/scale_holonomy.md`;
12. `docs/multiscale_bellman.md`;
13. `docs/nested_grain_extraction.md`;
14. `docs/cycle_hodge_flat_rigidity.md` and `docs/spherical_flat_erosion.md`;
15. `docs/atomic_component_entropy.md`;
16. `docs/multicommodity_hodge_routing.md`;
17. `docs/resistance_bellman_stopping.md`;
18. `docs/master_no_escape.md`.
