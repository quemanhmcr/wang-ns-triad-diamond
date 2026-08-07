# Odd-Hermite triad selection and quadratic daughter rescue

A coherent sideband is useful only if it can be connected back to the physical transfer ledger.  The key mechanism is an exact parity selection rule.

## 1. Resonant Gaussian parity

Let the three affine Gaussian roles be centered at resonant carriers

\[
\kappa_1+\kappa_2=\kappa_3.
\]

After writing centered deviations `eta_1,eta_2` and `eta_3=eta_1+eta_2`, the trilinear Gaussian weight is invariant under simultaneous inversion

\[
(\eta_1,\eta_2,\eta_3)\mapsto(-\eta_1,-\eta_2,-\eta_3).
\]

Therefore for polynomial/Hermite factors of parities `(-1)^n`,

\[
\boxed{
T(P_{n_1}G_1,P_{n_2}G_2,P_{n_3}G_3)=0
\quad\hbox{if }n_1+n_2+n_3\hbox{ is odd}.
}
\]

This statement is independent of anisotropic covariance because it uses only centered Gaussian inversion symmetry.

In particular,

\[
\boxed{T(H_1,G,G)=T(H_3,G,G)=0.}
\]

So one odd daughter cannot rescue the original Gaussian transfer at first order.

## 2. Rescue begins quadratically

Write each role as a base Gaussian part plus an odd sideband, with critical norms `b_i` and `rho_i`.  Expanding the trilinear form, all one-sideband terms vanish.  Sharp Young bounds the rest by

\[
\boxed{
|T_{rescue}|
\le A_3(
\rho_1\rho_2b_3+
\rho_1\rho_3b_2+
\rho_2\rho_3b_1+
\rho_1\rho_2\rho_3).
}
\]

Thus transfer rescue requires a genuine **pair of daughter sidebands**.  In the component graph this is a new sideband interaction component rather than a hidden correction to the base edge.

## 3. If there is no second daughter, one role pays a quadratic norm cost

Take one role `G+R=(1+P)G`, with `P` odd of degree at most three, and let

\[
\sigma^2=\mathbb E_{|G|^{3/2}}|P|^2.
\]

Because the measure is inversion symmetric,

\[
E|1+P|^{3/2}
={1\over2}E\left(|1+P|^{3/2}+|1-P|^{3/2}\right).
\]

For complex `z` with `|z|<=1/2`, the minimum at fixed `|z|` occurs for real `z`, and one-dimensional convexity gives

\[
{ |1+z|^{3/2}+|1-z|^{3/2}\over2}
\ge1+{3\over8}|z|^2.
\]

Degree-three Gaussian hypercontractivity yields

\[
E|P|^4\le729\sigma^4.
\]

If `sigma<=1/80`, then at least half the second moment lies on `|P|<=1/2`, hence

\[
\boxed{
{\|G+R\|_{3/2}^{3/2}\over\|G\|_{3/2}^{3/2}}
\ge1+{3\over16}\sigma^2.
}
\]

The numerator of the base--base--child transfer is unchanged because the one-odd term is exactly zero.  After normalizing the child role,

\[
\boxed{\operatorname{Def}_{transfer}\ge{1\over16}\sigma^2.}
\]

Therefore a coherent odd daughter has only two efficient possibilities:

\[
\boxed{
\text{second odd daughter / pair interaction}
\quad\text{or}\quad
\text{quadratic transfer deficit}.
}
\]

If `sigma>1/80`, there is already a definite sideband-capacity event.  This is the desired daughter-mode alternative without calling affine anisotropy itself a cost.
