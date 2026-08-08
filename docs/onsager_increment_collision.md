# Onsager cubic increments: mass, enstrophy, or collision entropy

The SGS source and boundary ledgers both produce a cubic velocity-increment charge at a physical filter scale. This note converts that observable into currencies already present in the no-escape architecture without imposing a global `L^3` hypothesis.

The only analytic input is the standard 3D Littlewood--Paley square-function theorem in `L^3` and band Bernstein. All low/high summation and entropy routing below are exact.

Let

\[
u=\sum_{j\in\mathbb Z}u_j,
\qquad N_j=2^jN,
\qquad \mu_j=N_j\|u_j\|_2^2.
\]

Assume the upper Fourier radius of band `j` is at most `beta N_j`. For increments with `|r|<=R_G/N`, Plancherel gives

\[
\|\delta_ru_j\|_2
\le \min(2,\beta R_G2^j)\|u_j\|_2.
\]

Bernstein therefore has exactly critical scaling:

\[
\boxed{
\|\delta_ru_j\|_3
\le C_B\min(2,\beta R_G2^j)\sqrt{\mu_j}.
}
\]

The square-function estimate then yields

\[
\|\delta_ru\|_3^2
\le (C_{LP}C_B)^2
\sum_j\min(4,(\beta R_G2^j)^2)\mu_j.
\]

For `j<=0`, the increment multiplier kills coarse modes geometrically:

\[
\sum_{j\le0}(\beta R_G)^24^j\mu_j
\le \frac43(\beta R_G)^2\mu_{\le N}^{\max}.
\]

For `j>=1`, put

\[
\mathfrak d_{>N}
=\sum_{j\ge1}2^j\mu_j
\simeq N^{-1}\|\nabla P_{>N}u\|_2^2.
\]

Since `4<=2*2^j` on those bands,

\[
\sum_{j\ge1}4\mu_j\le2\mathfrak d_{>N}.
\]

If `g_1=||G||_1`, the weighted cubic increment charge satisfies

\[
\boxed{
Q_N
\le
g_1(C_{LP}C_B)^3
\left[
\frac43(\beta R_G)^2\mu_{\le N}^{\max}
+2\mathfrak d_{>N}
\right]^{3/2}.
}
\]

This has no logarithmic dependence on the number of frequency bands.

## Collision theorem

Define

\[
X=\left[
\frac{Q_N}{g_1(C_{LP}C_B)^3}
\right]^{2/3}.
\]

Then

\[
X\le \frac43(\beta R_G)^2\mu_{\le N}^{\max}+2\mathfrak d_{>N}.
\]

Consequently

\[
\boxed{
\mu_{\le N}^{\max}
\ge\frac{3X}{8(\beta R_G)^2}
}
\]

or

\[
\boxed{
\mathfrak d_{>N}\ge\frac X4.
}
\]

The first branch is a scale-critical energy reservoir at a base/coarser band. The second is a high-frequency normalized enstrophy event.

## Aggregate band mass cannot hide in infinitely many packets

Suppose a transfer-adapted spatial packetization of the winning band has nonnegative critical masses `mu_a` summing to `mu_band`, and put `w_a=mu_a/mu_band`. Fix `0<theta<1`.

Either

\[
\boxed{
\max_a\mu_a\ge\theta\mu_{band},
}
\]

which is a dominant packet to be classified as fresh or reused, or

\[
Q_{at}=\sum_aw_a^2\le\theta,
\qquad
\boxed{H_{at}\ge-\log\theta.}
\]

With ancestry labels, the existing exact atomic/component chain rule gives for any `0<alpha<1`:

\[
H_{anc}\ge\alpha(-\log\theta)
\]

or

\[
\boxed{
Q_{anc}-Q_{at}
\ge \theta^\alpha-\theta.
}
\]

Thus spatial fragmentation is not free: it becomes Bellman component entropy or repeated same-ancestry attachments/cycles.

## Spacetime interpretation

Use scaled time `tau=N^2t`. If the enstrophy alternative

\[
\mathfrak d_{>N}\ge d_0
\]

holds on a set of scaled-time measure `m`, then the normalized high-frequency dissipation satisfies

\[
\boxed{
N\int\|\nabla P_{>N}u\|_2^2dt
\ge md_0.
}
\]

The later source-weighted theorem sharpens this spacetime interpretation.  For differentiated SGS source, the source-to-increment power `3/2` and the increment-to-mass/enstrophy power `2/3` cancel, so the final currency is linear in instantaneous source density.  No persistence or temporal-superlevel assumption is required.  For the viscous source, Cauchy shows that concentrating source weight in time only increases the quadratic dissipation price.

Thus the preferred physical routing is

\[
\boxed{
\text{cubic SGS increments}
\to
\text{dominant fresh/reused packet}
\ \lor\
\text{Bellman/cycle entropy}
\ \lor\
\text{high-frequency dissipation}.
}
\]
