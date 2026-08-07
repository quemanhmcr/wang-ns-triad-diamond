# Research Ledger — Wang–Zahl-style rigidity programme for 3D Navier–Stokes

**Status date:** 2026-08-07  
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

### Computationally supported, not interval-certified

- numerical reuse-gap constants from nonlinear optimization;
- perturbative robustness tables for near-butterfly / finite-width Gaussian models.

### Still conditional / PDE bridge

1. **Realize the frozen frequency cells as a full three-dimensional moving
   packet frame.** Frequency multiplier freezing and crossing-scale selection are
   now controlled, but the continuum theorem must attach the cells to physical
   packets/windows while preserving the signed localized SGS measure and keeping
   partition overlap/cross interactions summable.
2. **Remove the invariant-plane reduction.** The affine and five-percent
   coherence theorems rigorously charge non-conformal strain on an invariant
   triad plane.  A genuine 3D grain can tilt out of that plane.  Off-plane strain,
   plane precession, helical-polarization transport, and their coupling to the
   local vorticity frame must either be absorbed as gauge motion or charged to a
   new positive deformation/fresh-mass branch.
3. **Integrate the adaptive curvature-balanced moat into the nested packet
   tree.** The old fixed spatial moat schedule is not sufficient.  The continuum
   construction must choose `M` from the local curvature balance `a/M+b kappa M`,
   prove compatibility of neighboring adaptive windows, and verify the local
   pressure-Hessian coefficient.  Excessive curvature already forces critical
   mass in the band-limited packet model.
4. **Close the remaining localized transport terms.** Pressure has a combined-
   work ledger and far multipole bounds; affine advection is cancelled by the
   deforming window.  The remaining `R U` SGS transport, viscous boundary flux,
   local pressure-Hessian contribution and partition-overlap terms must be
   summable or produce Bellman/fresh critical mass.
5. **Spacetime ancestry and synchronization.** The one-shot Gaussian profile
   carries certified frozen-time critical mass and coherent non-conformal strain
   now pays a packet-lifetime cost.  The remaining theorem must register the
   resulting mass as fresh or reused across successive `N^{-2}` lifetimes and
   feed all deformation/source errors into the master `zeta_j`/cross-error
   ledgers.

No statement in this repository currently closes these PDE gaps, and no claim of Navier–Stokes global regularity is made.

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

The current preferred master regression is run `31169621311` on the integrated
affine-spacetime bridge commit `2b795d2` (`20,000` episode traces, worst margin
`0`).  The earlier integrated frequency/pressure bridge run `31166171000`
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
`recorded-results/31168888413/`, and `recorded-results/31169128097/`; the two
failed interval-boundary provenance runs are `recorded-results/31168205564/`
and `recorded-results/31168303213/`.  Final integration runs are
`31166152074`, `31166155045`, `31166158218`, `31166160711`, `31166163414`,
`31166165985`, `31166168659`, and `31166171000`.  The numerical full-Mellin search is recorded
in `recorded-results/31160779428/`; its qualitative counterexample is separately
Arb-certified in run `31161914134`.

## 18. Current research frontier

The finite-dimensional architecture has reached the point where every identified escape variable has a ledger:

\[
\boxed{
\text{branching}
\to \text{Bellman entropy},
\quad
\text{curved reuse}
\to \text{Hodge cost},
}
\]

\[
\boxed{
\text{flat reuse}
\to \text{spherical/barycentric erosion},
\quad
\text{high resistance}
\to \text{Poisson Bellman/atomic entropy}.
}
\]

The single-edge multiplier, smooth scale-flux observable, cancellation
polarization and physical-transfer Hodge change of measure are closed in the
finite packet model.  The common spectral moat is now also extracted from any
signed-good crossing core rather than assumed.  Conservatively, after the
four-bin physical crossing selection one may use

\[
\boxed{
c_{0,H}^{cross}=\min\{1/20000,(25/424)h_H\}>0.
}
\]

The frozen-frequency bridge has now acquired a genuine packet-lifetime dynamics
ledger.  The low-cost alternatives can be organized physically as

\[
\boxed{
\text{coherent non-conformal planar strain}
\to \text{multiplier deficit},
}
\]

with the certified costs `>= (1/10)(dT)^2` in the frozen case and
`>= (1/24)(dT)^2` under five-percent co-rotating variation, or

\[
\boxed{
\text{strain dephasing}
\to Q(A)\ \text{or}\ \nabla^2p\ \text{or}\ \nu\Delta S,
}
\]

where the Navier--Stokes gradient equation supplies the source identity and the
band-limited packet model converts each source into critical mass or a summable
far pressure-Hessian tail.  Spatial localization also has a corrected physical
law:

\[
\boxed{
E(M)=a/M+b\kappa M,
\qquad M_*\sim\kappa^{-1/2}.
}
\]

The old quadratic `M_j` schedule remains valid for the defect-space extraction
and for the commutator term alone, but a fixed spatial schedule cannot control
velocity curvature in general.

The next mathematically decisive work is therefore:

1. construct a **full 3D transfer-preserving moving packet frame**, including
   off-plane tilt/helical-polarization transport, rather than an invariant-plane
   reduction;
2. make the **adaptive curvature-balanced windows** compatible with the nested
   packet ancestry, and prove that local pressure-Hessian, SGS transport `RU`,
   viscous boundary flux and partition overlap are summable or create a fresh /
   Bellman event;
3. synchronize the resulting critical grains over successive `N^{-2}` packet
   lifetimes and feed objective-strain/curvature events into the master
   no-escape ledger.

This is now the frontier.  Further single-edge optimization, another abstract
cycle theorem, or an iterative within-block profile decomposition would not
address the principal missing PDE mechanism.

---

### Repository reading order

For a new reader, the recommended order is:

1. `RESEARCH_LEDGER.md` — this document;
2. `docs/single_edge_stability_certificate.md`;
3. `docs/log_scale_flux_bridge.md` and `docs/smooth_log_flux_cocycle.md`;
4. `docs/crossing_moat_extraction.md` and `docs/smooth_sgs_symbol_freezing.md`;
5. `docs/localized_sgs_pressure_ledger.md` and `docs/annular_pressure_collision.md`;
6. `docs/transfer_preserving_profile_extraction.md`, `docs/gaussian_packet_inverse.md` and `docs/packet_inverse_theorem.md`;
7. `docs/affine_gaussian_grain_dynamics.md` and `docs/strain_coherence_objective_gradient.md`;
8. `docs/curvature_balanced_moat.md` and `docs/objective_strain_source_collision.md`;
9. `docs/scale_holonomy.md`;
10. `docs/multiscale_bellman.md`;
11. `docs/nested_grain_extraction.md`;
12. `docs/cycle_hodge_flat_rigidity.md` and `docs/spherical_flat_erosion.md`;
13. `docs/atomic_component_entropy.md`;
14. `docs/multicommodity_hodge_routing.md`;
15. `docs/resistance_bellman_stopping.md`;
16. `docs/master_no_escape.md`.
