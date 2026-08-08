# Affine coherent Moyal identity: an exact positive reservoir-energy ledger

Discrete packet coefficients can be redundant. The physical object behind them is the coherent-state analysis energy density. For any L2-normalized affine Gaussian window `g_L`, define

\[
\mathcal V_L f(X,k)
=
\int f(x)\overline{g_L(x-X)}e^{-ik\cdot x}\,dx.
\]

Plancherel in `k` gives, for fixed `X`,

\[
\int |\mathcal V_Lf(X,k)|^2\frac{dk}{(2\pi)^3}
=
\int |f(x)|^2|g_L(x-X)|^2dx.
\]

Integrating in `X` and using `||g_L||_2=1` yields

\[
\boxed{
\int_{\mathbb R^6}|\mathcal V_Lf(X,k)|^2
\frac{dX\,dk}{(2\pi)^3}
=
\|f\|_2^2.
}
\]

No isotropy is assumed. Affine shape disappears through the window normalization.

## Positive phase-space reservoir cells

For any measurable phase-space partition, define

\[
E_\alpha
=
\int_{C_\alpha}|\mathcal V_Lf|^2d\mu.
\]

Then

\[
\boxed{E_\alpha\ge0,\qquad\sum_\alpha E_\alpha=\|f\|_2^2.}
\]

Thus coherent cells give a canonical positive reservoir-energy decomposition with exact frame constant `P=1`. Close coherent probes are merged by integration over one cell instead of appearing as many cancelling synthesis coefficients.

For an orthogonal dyadic frequency partition, applying Moyal separately to each old band gives exactly the energy input needed by `reservoir_pool_erosion.md`. A smooth Littlewood--Paley partition pays only its standard fixed square-function overlap constant.

## Relation to the discrete Riesz theorem

Moyal is continuous and analysis-side. The 5-separated theorem in `affine_coherent_bessel.md` supplies the complementary synthesis statement

\[
47I/50\le G\le53I/50.
\]

Therefore ancestry/Bellman bookkeeping can use positive Moyal cell energies, while a separated Riesz representative family is available whenever a discrete Gaussian synthesis is needed.

The remaining bridge is not a generic Bessel estimate. It is to prove that the **physical transfer service selected by a Navier--Stokes block can be assigned to these coherent phase-space cells**, and that changing or merging the cells through recursive generations creates only summable `Xi` or an already named fresh/relink/source cost.
