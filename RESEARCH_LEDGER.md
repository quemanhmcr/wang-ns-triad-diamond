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

### Computationally supported, not interval-certified

- numerical reuse-gap constants from nonlinear optimization;
- perturbative robustness tables for near-butterfly / finite-width Gaussian models.

### Still conditional / PDE bridge

1. **Derive the actual smooth-SGS affine packet equations and their differential
   sources.**  Common affine/non-affine phase transport, bulk Gaussian viscosity,
   common-parent `SL(2)` motion, scalar third-Hermite forcing and RMS spatial
   polarization curvature are now separated.  The remaining `F_i` must be
   derived from role-dependent SGS transport, partition/window errors, pressure
   localization and differential resolved velocities in the actual packetization.
2. **Insert the affine moving-window theorem into the full localized SGS/pressure
   identity.**  Ellipsoidal material transport, the `1/M` filter commutator and
   the `kappa_aff M` Taylor leakage are now exact.  What remains is to track the
   actual pressure boundary work, `RU`, viscous boundary flux and partition
   overlap with the same affine windows and prove their coefficients are
   summable or force a source/fresh event.
3. **Integrate affine critical grains into spacetime ancestry.**  A fresh grain
   now carries theorem-level scale-critical mass `M_aff>=3/10` and obeys the
   physical budget `sum r_g<=P E/eta`.  Reused elongated grains still require a
   synchronization theorem across successive lifetimes and a rule for comparing
   evolving covariance ellipsoids without quotienting away genuine differential
   curvature/polarization.
4. **Close the remaining localized transport/source terms.**  Pressure has a
   combined-work ledger and far multipole bounds; bulk viscosity and common
   material phase are packet tangents.  Viscous boundary flux, `R U`, local
   pressure-Hessian contribution, role-dependent phase forcing and partition
   overlap must be summable or create fresh/Bellman/source events.
5. **Feed the affine forcing channels into the master no-escape ledger.**  The
   third-Hermite scalar channel, RMS relative-polarization curvature, affine fresh
   radius budget and existing objective-strain/source alternatives must be
   synchronized with the master `zeta_j`/cross-error ledgers.  No continuum
   theorem currently proves that every efficient PDE block enters one of those
   branches with uniform constants.

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

| `31171921187`--`31171950823` | integrated helical/spacetime stack on `6226fd9`: spin transport, explicit spin-dihedral phase holonomy, full-strain tomography, objective polarization, intrinsic 3D plane, affine grain, strain coherence and master all green; `153` tests per workflow, master `20,000` traces with worst margin `0` |

The current preferred master regression is run `31182017032` on the eleven-workflow
affine-residual integration commit `1816ee2` (`189` tests plus `20,000` episode
traces, worst margin `0`).  The earlier relative-polarization master
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
`recorded-results/31168888413/`, and `recorded-results/31169128097/`, plus the intrinsic 3D continuation `recorded-results/31170015795/`.  The preferred helical artifacts are
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

The finite-dimensional and affine-packet architecture now distinguishes **true
physical defects** from symmetry motion much more sharply.  In addition to the
Bellman/Hodge/erosion/reuse ledgers, the packet dynamics has the gauge hierarchy

\[
\boxed{
\text{common material phase/chirp},
\quad
SO(3)\text{ helical frame},
\quad
SL(2)\text{ common-parent polarization},
\quad
\text{affine Gaussian tangent motion},
}
\]

all of which must be quotiented rather than charged.

After those quotients, the non-affine resolved-flow curvature splits physically
into

\[
\boxed{
\operatorname{Sym}(L^{-1}(\nabla^2U)[L,L])
\to
\text{third-Hermite envelope forcing},
}
\]

and

\[
\boxed{
C_{ijc}=\operatorname{sym}_{ij}((\nabla^2U)_{ijk}L_{kc})
\to
\text{spatial shape/relative-polarization variation}.
}
\]

The scalar kernel is the exact quadratic swirl `z cross(Mz)` rather than a
mysterious escape.  The common velocity-Hessian chirp is a material triad phase
gauge.  Bulk viscosity is Gaussian tangent.

The inverse-Young grain is now also physically registered without an isotropic
fiction.  With

\[
r_g=(\det\Sigma_x)^{1/6},
\qquad
\mathsf M_{aff}=r_g^{-1}\int_{E_2}|u|^2,
\]

Arb and Hausdorff--Young give

\[
\boxed{\mathsf M_{aff}\ge3/10,}
\]

while fresh grains satisfy

\[
\boxed{\sum r_g\le P\|u(t)\|_2^2/\eta.}
\]

Static affine anisotropy itself is an exact Young symmetry and is **not** a
replication cost.  This corrects the earlier natural-cell heuristic.  The
geometric radius is preserved by inviscid incompressible affine strain and is
monotone under viscosity, making it a natural spacetime ancestry currency.

The next mathematically decisive work is therefore:

1. derive the **actual affine smooth-SGS packet equations** and express their
   role-dependent residuals in the already-identified currencies: third-Hermite
   forcing, differential phase velocity, relative-polarization curvature,
   pressure/source terms and partition/window errors;
2. insert the now-exact **ellipsoidal moving window and filter commutator** into
   the localized SGS/pressure identity, deriving pressure, `RU`, viscous boundary
   and partition-overlap coefficients in the affine metric;
3. construct **affine spacetime ancestry**: use the theorem-level fresh radius
   budget for new grains, while controlling reuse by covariance synchronization,
   scale/spin holonomy and the existing dynamic strain/polarization ledgers;
4. feed these affine PDE branches into the master `zeta_j`, `Xi` and
   fresh/reuse/Bellman bookkeeping with summable errors.

This is now the frontier.  Further single-edge optimization, raw Euclidean
propagator norms, automatic isotropic subdivision of affine Gaussians, or a full
Magnus theorem for common parent motion would attack symmetry variables rather
than the remaining Navier--Stokes mechanism.

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
9. `docs/affine_gaussian_forcing.md`, `docs/material_phase_lock.md`, `docs/quadratic_swirl_kernel.md`, `docs/affine_polarization_curvature.md`, `docs/affine_shell_aspect.md`, `docs/affine_critical_grain.md`, and `docs/affine_window_balance.md`;
10. `docs/curvature_balanced_moat.md` and `docs/objective_strain_source_collision.md`;
11. `docs/scale_holonomy.md`;
12. `docs/multiscale_bellman.md`;
13. `docs/nested_grain_extraction.md`;
14. `docs/cycle_hodge_flat_rigidity.md` and `docs/spherical_flat_erosion.md`;
15. `docs/atomic_component_entropy.md`;
16. `docs/multicommodity_hodge_routing.md`;
17. `docs/resistance_bellman_stopping.md`;
18. `docs/master_no_escape.md`.
