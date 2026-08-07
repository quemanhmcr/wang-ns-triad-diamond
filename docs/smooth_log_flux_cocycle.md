# Smooth log-scale flux cocycle and transfer-weighted midgap theorem

This module is a finite Fourier/Galerkin theorem plus a PDE-facing design
principle.  It does not by itself construct the required packet decomposition of
a Navier--Stokes solution.

## 1. Conservation turns a graded filter into a log-scale convolution

Write `t=log K` and let a finite conservative interaction have modal log
frequencies `ell_i` and nonlinear energy rates `dot E_i` with

\[
\sum_i \dot E_i=0.
\]

Its sharp cumulative outward spectral flux is

\[
\Pi^0(t)=-\sum_i \dot E_i\,H(t-\ell_i).
\]

Let `rho_delta` be any nonnegative probability density on log scale and

\[
\Psi_\delta(z)=\int_{-\infty}^z\rho_\delta(s)\,ds.
\]

The graded low-pass quadratic energy with spectral weight
`Psi_delta(t-log|xi|)` has nonlinear flux

\[
\Pi^\delta(t)
=-\sum_i\dot E_i\Psi_\delta(t-\ell_i).
\]

Because `Psi=rho*H`, exactly

\[
\boxed{\Pi^\delta=\rho_\delta*\Pi^0}
\]

in the log-cutoff variable.  Thus smooth spectral coarse graining is literally
a renormalization convolution of the sharp scale flux.  No turbulence closure is
used.

In particular, whenever the finite interaction is conservative,

\[
\boxed{
\int_{\mathbb R}\Pi^\delta(t)\,dt
=\int_{\mathbb R}\Pi^0(t)\,dt
=\sum_i\dot E_i\ell_i.
}
\]

The all-log Mellin moment is therefore **filter independent**.  The sharp
identity in `docs/log_scale_flux_bridge.md` was not an artifact of the sharp
cutoff.

This is a conservation/RG identity: changing the graded filter only redistributes
flux in the variable `log K`; it cannot change its total log-scale moment.

## 2. Weak scale-potential identity

For any Lipschitz scale potential `Phi`, distributional integration by parts
with the cumulative energy balance gives

\[
\boxed{
\int_{\mathbb R}\Pi^0(t)\Phi'(t)\,dt
=\sum_i\dot E_i\Phi(\ell_i).
}
\]

This is the correct weak formulation of energy transport in scale.  It is more
useful than committing to one scalar Mellin moment, because the test potential
can isolate the generation boundary relevant to a packet block.

## 3. Countermodel: the full Mellin moment is not the right edge functional

The tempting replacement of `mathcal J` by the full all-scale moment is false.
For a triad `k<=p<=q`,

\[
\int\Pi_K\frac{dK}{K}
=-\dot E_k\log\frac pk+\dot E_q\log\frac qp.
\]

The first term can be a large *forward* transfer from the lowest mode to the
middle mode.  It measures legitimate spectral migration, but it is not the
specific parent-to-child crossing represented by one generation edge.

The Actions adversarial search in `physical-mellin-edge.yml` found a full-Mellin
coefficient materially larger than `J_*`, on an asymmetric helicity branch.  So
all-scale Mellin saturation does **not** force the old single-edge grain.

The correction is structural rather than cosmetic: use a block-local scale
potential which tests the physical flux only across the parent/child gap.

## 4. Smooth midgap tail flux isolates the upper crossing exactly

Now take `rho_delta` **even**, supported in `[-delta,delta]`.  Then

\[
\Psi_\delta(z)+\Psi_\delta(-z)=1.
\]

Define

\[
\Phi_{\tau,\delta}'(t)=2\Psi_\delta(t-\tau),
\qquad
\Phi_{\tau,\delta}(t)=0\quad(t\le\tau-\delta).
\]

Evenness gives the exact high-side continuation

\[
\boxed{
\Phi_{\tau,\delta}(t)=2(t-\tau)
\quad\text{for }t\ge\tau+\delta.
}
\]

Equivalently, since the graded flux is the log convolution of the sharp flux,

\[
\boxed{
2\int_\tau^\infty\Pi^\delta(t)\,dt
=\sum_i\dot E_i\Phi_{\tau,\delta}(\ell_i).
}
\]

Consider one ordered forward triad and place the transition in the empty
parent/child gap:

\[
\ell_k\le\ell_p\le\tau-\delta,
\qquad
\ell_q\ge\tau+\delta.
\]

Then both parents lie on the flat side and the child lies on the linear side,
so

\[
2\int_\tau^\infty\Pi^\delta(t)dt
=2\dot E_q(\ell_q-\tau).
\]

Choosing the **actual edge midgap**

\[
\tau_e=\frac{\ell_p+\ell_q}{2}
\]

gives, with no smoothing error at all,

\[
\boxed{
2\int_{\tau_e}^\infty\Pi_e^\delta(t)dt
=\dot E_q(\ell_q-\ell_p)
=\dot E_q\log\frac qp.
}
\]

Thus the progress factor in `mathcal J` has a smooth-filter realization.  The
factor `2` is not arbitrary: it is exactly the calibration converting distance
from the midpoint cross-section to the full parent-child log separation.

## 5. Transfer-weighted barycentric midgap: one cutoff for an entire block

The edge-dependent cutoff can be removed on a forward packet core.  Let edges
have positive child transfers `T_e>0`, top-parent/child logs `(p_e,q_e)`, and
midgaps

\[
m_e=\frac{p_e+q_e}{2}.
\]

Set the cutoff at the **transfer-weighted midgap barycenter**

\[
\boxed{
\tau=\frac{\sum_eT_em_e}{\sum_eT_e}.
}
\]

Assume a common smooth spectral moat,

\[
p_e\le\tau-\delta,
\qquad q_e\ge\tau+\delta
\quad\text{for every edge in the core}.
\]

Linearity and the preceding isolation lemma give

\[
2\int_\tau^\infty\Pi_{\rm core}^\delta(t)dt
=\sum_e2T_e(q_e-\tau).
\]

But

\[
2(q_e-\tau)-(q_e-p_e)=2(m_e-\tau),
\]

and the transfer-weighted sum of the right side vanishes by definition of
`tau`.  Hence

\[
\boxed{
2\int_\tau^\infty\Pi_{\rm core}^\delta(t)dt
=\sum_eT_e(q_e-p_e)
=\sum_eT_e\log\frac{Q_e}{P_e}.
}
\]

This is the key new cocycle.  A **single smooth physical scale cut** reproduces
the whole upper-progress ledger of a forward block exactly, not approximately.
The same transfer weights that locate the cutoff are the weights entering the
progress ledger.

A sufficient shell condition is easy to state.  If all top parents lie in a
log shell `|p_e-L|<=sigma` and all children in
`|q_e-(L+gamma_*)|<=sigma`, then the weighted `tau` lies within `sigma` of
`L+gamma_*/2`.  A common moat is guaranteed whenever

\[
\boxed{\delta+2\sigma<\gamma_*/2.}
\]

For example `sigma=0.06`, `delta=0.05` leaves a definite margin at the certified
`gamma_*`.

## 6. Cancellation becomes a positive cost, not an ambiguity

Let `A_e>=0` be unsigned trilinear capacities, normalized by
`w_e=A_e/sum A`.  Write the signed edge efficiency as

\[
r_e=m_ec_e,
\qquad 0\le m_e\le1,
\qquad -1\le c_e\le1,
\]

where `m_e` contains the normalized helical/progress multiplier and `c_e` the
phase/orientation factor.  The normalized signed block flux is

\[
R=\sum_ew_em_ec_e.
\]

There is an exact nonnegative decomposition

\[
\boxed{
1-R
=\sum_ew_e(1-m_e)
+\sum_ew_em_e(1-c_e).
}
\]

Therefore near saturation simultaneously implies small mean multiplier deficit
and small phase/backscatter defect.  In particular

\[
\sum_ew_e\operatorname{Def}_e\le1-R.
\]

For every `eta>0`, Markov gives

\[
w\{m_ec_e\le1-\eta\}\le\frac{1-R}{\eta}.
\]

Thus a near-extremal signed PDE block has a large positive-transfer core; bad
phase/cancellation cannot secretly destroy the positive weights needed by the
midgap theorem.  On the local single-edge core,

\[
\operatorname{Def}_e\ge\frac12\mathcal D_e,
\]

so the same physical capacity measure pays the Hodge residual.

## 7. Physical-space localization: the moat schedule has the right commutator scaling

Let a smooth spectral filter have physical kernel

\[
G_K(x)=K^3G(Kx),
\qquad M_1(G)=\int |z||G(z)|dz<\infty.
\]

For a Lipschitz spatial cutoff `chi`, the exact commutator representation is

\[
[G_K,\chi]f(x)
=\int G_K(y)(\chi(x-y)-\chi(x))f(x-y)dy.
\]

Young's inequality gives, for `1<=p<=infinity`,

\[
\boxed{
\|[G_K,\chi]f\|_p
\le\frac{M_1(G)}{K}\|\nabla\chi\|_\infty\|f\|_p.
}
\]

If the spatial transition has width `M/K`, the cost is `O(1/M)`.  Hence the
same schedules already used by nested grain extraction, e.g.

\[
M_j=(j+3)^2,
\qquad \sum_jM_j^{-1}<\infty,
\]

make smooth-filter/spatial-window commutators summable.  This is a genuine PDE
localization mechanism rather than a decorative analogy.

Pressure is absent from the **global** graded spectral transfer because the
scalar spectral filter preserves divergence-free fields and

\[
\langle A_Ku,\nabla p\rangle=0.
\]

After spatial localization pressure reappears only as boundary work in the
window transition region.  Charging that pressure work to the same moat/local
energy ledger remains a PDE task; it is not claimed solved here.

## 8. What remains

The new exact identities remove two previous ambiguities:

1. smoothing the spectral cutoff does not change the all-log flux moment;
2. more importantly, a transfer-weighted common midgap cut recovers the desired
   upper-progress ledger **exactly** on a spectrally separated forward core.

The remaining hard PDE step is now sharply formulated: extract from a genuine
near-extremal Navier--Stokes block a finite/summable packet interaction core with
(a) positive child-transfer mass after the polarization loss, (b) a common
log-scale moat, and (c) summable spatial-window/pressure leakage.  Once that is
available, the certified single-edge deficit is already measured with the same
physical transfer weights needed by Bellman/Hodge.

## 9. Near saturation aligns capacity weights with physical child-transfer weights

There are two positive measures in the architecture which must not be silently
identified.  The weighted-Young/Bellman side naturally produces an unsigned
trilinear capacity `A_e`, while the midgap location is generated by the actual
positive child-energy transfer `T_e`.  Near extremality makes them quantitatively
equivalent.

Let the signed normalized progress efficiency be

\[
r_e=m_ec_e,
\qquad
F_e=A_eJ_*r_e=T_e g_e,
\qquad
g_e=\log(q_e/p_e).
\]

If the block has deficit `epsilon=1-R`, then for every `eta>0` the signed-good
core

\[
G_\eta=\{e:r_e>1-\eta\}
\]

has capacity mass

\[
\boxed{w(G_\eta)\ge1-\epsilon/\eta.}
\]

Take `eta<1/100`.  The certified global exclusion says every edge in `G_eta`
must lie in the local stability box.  Moreover `r_e<=m_e`, so
`Def_e=1-m_e<=eta`.  The mixed certificate gives

\[
u_e\le50\eta,
\qquad |v_e|\le\sqrt\eta.
\]

Since

\[
g_e=\gamma_*+v_e-u_e/2,
\]

we obtain the explicit physical gap rigidity

\[
\boxed{
|g_e-\gamma_*|\le a_\eta:=\sqrt\eta+25\eta.
}
\]

On this core

\[
\frac{J_*(1-\eta)}{\gamma_*+a_\eta}
\le\frac{T_e}{A_e}
\le
\frac{J_*}{\gamma_*-a_\eta}.
\]

Thus the condition number between child-transfer density and capacity density is
at most

\[
\boxed{
\kappa_\eta
=\frac{\gamma_*+a_\eta}
{(1-\eta)(\gamma_*-a_\eta)}
\longrightarrow1
\quad(\eta\downarrow0).
}
\]

This closes a conceptual loop: near saturation does not merely make phases
positive; it forces the *physical* transfer weights locating the smooth midgap
to converge to the same transfer-capacity weights used by the inverse/Bellman
ledger.  No statistical-mechanics analogy is being inserted—the equivalence is
forced by the certified multiplier deficit and the actual child-energy flux.

## 10. The graded spectral flux is the space-average of the physical SGS flux

The smooth cocycle is not restricted to Fourier bookkeeping.  Choose a
self-adjoint convolution filter `G_t` whose Fourier multiplier satisfies

\[
|\widehat G_t(\xi)|^2
=\Psi_\delta(t-\log|\xi|).
\]

Write

\[
\bar u=G_t*u,
\qquad
\tau_t(u,u)=\overline{u\otimes u}-\bar u\otimes\bar u.
\]

The standard resolved/subfilter transfer density is

\[
\boxed{
\Pi_t^{\rm SGS}(x)
=-\nabla\bar u:\tau_t(u,u).
}
\]

For a periodic box or sufficient spatial decay, incompressibility gives

\[
\boxed{
\int \Pi_t^{\rm SGS}(x)dx
=-\frac d{dt_{\rm physical}}\frac12\|\bar u\|_2^2\Big|_{NL}
=\Pi^\delta(t).
}
\]

The first equality follows from the filtered energy equation: the cubic
resolved transport integrates to zero.  In Fourier variables the resolved
quadratic energy has weight `|Ghat|^2=Psi_delta`, which is exactly the graded
spectral energy used above.

The Leray projector does not alter the global nonlinear work because `bar u` is
divergence free and the projector is self-adjoint.  Equivalently, pressure does
no global interscale work.  This is the physically important distinction:
pressure transports energy in space but is not itself the global scale-transfer
mechanism.

Consequently the transfer-weighted midgap theorem has the physical-space form

\[
\boxed{
2\int_\tau^\infty\int \Pi_t^{\rm SGS}(x)\,dx\,dt
=\sum_{e\in G}T_e\log(q_e/p_e)
}
\]

on a forward packet core satisfying the common spectral moat.

After multiplying by a spatial packet window, the pressure term becomes a
boundary flux through the transition region.  Thus the unresolved pressure
problem is now localized very specifically: prove that this boundary work is
summable under the same expanding-moat schedule.  It no longer contaminates the
global definition of the transfer weights.

## 11. Certified change of measure and a physical Hodge constant

Fix the concrete signed-good threshold

\[
\eta_0=10^{-4}.
\]

On `G_{eta_0}` the previous formula simplifies exactly to

\[
a_{\eta_0}=\sqrt{10^{-4}}+25\cdot10^{-4}
=\frac1{80}.
\]

The Arb certificate now verifies

\[
\boxed{
\kappa_{\eta_0}<\frac{53}{50}.
}
\]

After normalizing both measures on the good core, if `w_e` is the capacity law
and `\widetilde w_e` the actual positive child-transfer law, then pointwise

\[
\boxed{
\frac{50}{53}\,w_e
\le \widetilde w_e
\le \frac{53}{50}\,w_e.
}
\]

Therefore every nonnegative quadratic residual, and in particular the Hodge
energy, changes by at most this factor under the physical change of measure.
This is the quantitative bridge missing from the older phrase "transfer-weighted
Hodge": the weights can now be taken to be actual positive child-energy
transfers on the near-saturated core.

There is also a clean cost-or-physical-Hodge consequence.  Let

\[
\epsilon=1-R.
\]

Markov gives capacity mass

\[
w(G_{\eta_0})\ge1-10^4\epsilon.
\]

If `epsilon>=1/20000`, the block already pays the explicit cost `1/20000`.
Otherwise the good core has capacity mass at least `1/2`.  Combining

\[
\operatorname{Def}_e\ge\frac12\mathcal D_e
\]

with the change of measure yields

\[
\epsilon
\ge \frac12\cdot\frac12\cdot\frac{50}{53}
\,\mathcal E_H^{\rm phys}
=\frac{25}{106}\mathcal E_H^{\rm phys}.
\]

Hence a physical-transfer Hodge branch with threshold `h_H` satisfies the
explicit finite-packet cost

\[
\boxed{
\epsilon\ge
\min\left\{\frac1{20000},\frac{25}{106}h_H\right\}>0.
}
\]

This constant is not a global Navier--Stokes theorem constant: the PDE still has
to produce the packet core and common spectral/spatial moat.  But once those
objects are extracted, the Hodge cost is now expressed in the weights of an
actual smooth SGS energy transfer, not in an auxiliary graph measure.
