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

It measures transfer strength times logarithmic scale progress. It is **not** claimed to be a literature-standard Navier–Stokes flux functional; the PDE bridge must still identify the corresponding exact scale-flux ledger.

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
`w_e=1/3`, the inherited arc weights give

\[
\boxed{
-\log R_{\rm block}\ge\frac{\gamma_*^2}{30}
\approx0.00809556352.
}
\]

This is a concrete theorem-level positive `c_0` for that normalized reuse block,
not yet a claim that the same number is the universal master constant across
all costly branches.

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
  the global outside-neighborhood gap `Def_e>=1/100`.

### Computationally supported, not interval-certified

- numerical reuse-gap constants from nonlinear optimization;
- perturbative robustness tables for near-butterfly / finite-width Gaussian models.

### Still conditional / PDE bridge

1. **Exact PDE scale-flux connection.** Relate the progress-weighted diagnostic `\mathcal J` to a rigorous Navier–Stokes scale-flux ledger with controlled symmetrization/localization errors.
2. **Gaussian atomic extraction from an arbitrary near-extremal PDE block.** Need an iterative profile extraction with controlled synthesis constants and transfer-small remainder, without assuming a global `\|\widehat u\|_{3/2}` bound.
3. **Summable perturbation ledger.** Convert near-extremal PDE errors into the `\zeta_j` and logarithmic cross-error terms required by the master theorem.
4. **PDE weighting of the certified edge deficit.** The finite-dimensional multiplier gap is now certified.  What remains is to prove that a genuine near-extremal Navier–Stokes scale-flux block produces the normalized transfer-weighted atomic edge measure to which
   `sum_e w_e Def_e >= (1/2) E_H` applies, with bad-edge removal and localization errors charged to the existing Bellman/cross-error ledgers.
5. **Critical-mass bridge.** Show that a concentrated companion grain carries enough local `L^2` / local-energy charge to be classified as fresh or reused in the physical packet ledger.
6. **Time synchronization and pressure/localization.** Lift the static/packet graph architecture to a spacetime argument compatible with suitable weak solutions and local energy estimates.

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
| `31157463384` | Arb-certified full single-edge sign/global/local stability; `76` tests + `200,000` adversarial samples |

The current preferred master regression artifact is `recorded-results/31154025683/`.
The preferred single-edge theorem certificate is `recorded-results/31157463384/`.

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

The finite-dimensional single-edge multiplier gap is now closed with a rigorous
positive constant.  The next mathematically decisive work should therefore move
**back toward the PDE**, not add more graph combinatorics.  The two
highest-priority targets are:

1. identify an exact/symmetrized Navier–Stokes scale-flux ledger whose
   transfer weights produce the certified edge deficit and hence the Hodge
   block cost `c_{0,H}=h_H/2`;
2. prove a transfer-adapted Gaussian atomic extraction theorem for genuine
   near-extremal Navier–Stokes blocks, with summable synthesis/localization
   errors.

If these PDE modules supply the transfer weights, critical mass, and summable
error ledgers required by the already-closed finite-dimensional architecture,
the master inequality would turn the no-escape mechanism into a quantitative
obstruction to an indefinitely efficient critical cascade.

---

### Repository reading order

For a new reader, the recommended order is:

1. `RESEARCH_LEDGER.md` — this document;
2. `docs/single_edge_stability_certificate.md`;
3. `docs/gaussian_packet_inverse.md` and `docs/packet_inverse_theorem.md`;
4. `docs/scale_holonomy.md`;
5. `docs/multiscale_bellman.md`;
6. `docs/nested_grain_extraction.md`;
7. `docs/cycle_hodge_flat_rigidity.md` and `docs/spherical_flat_erosion.md`;
8. `docs/atomic_component_entropy.md`;
9. `docs/multicommodity_hodge_routing.md`;
10. `docs/resistance_bellman_stopping.md`;
11. `docs/master_no_escape.md`.
