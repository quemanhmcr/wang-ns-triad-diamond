# Model packet inverse theorem

Let three L2-normalized isotropic Fourier Gaussian packets have common width
\(\sigma\), frequency centers \(\kappa_j\), and spatial centers \(x_j\). Divide
their scalar trilinear overlap by the aligned resonant value. Then

\[
R=\exp\!\left[-\frac{|\kappa_1+\kappa_2-\kappa_3|^2}{12\sigma^2}
-\frac{\sigma^2}{3}\sum_{i<j}|x_i-x_j|^2\right].
\]

Hence \(R\ge 1-\varepsilon\) implies

\[
|\kappa_1+\kappa_2-\kappa_3|^2
\le 12\sigma^2[-\log(1-\varepsilon)]
\]

and

\[
\sum_{i<j}|x_i-x_j|^2
\le 3\sigma^{-2}[-\log(1-\varepsilon)].
\]

For unequal widths with fixed \(S=\sum\sigma_j^2\), set
\(a_j=\sigma_j^2/S\). The relative width factor is

\[
R_w=(27a_1a_2a_3)^{3/4}.
\]

Pinsker's inequality applied to the uniform distribution and \((a_1,a_2,a_3)\)
gives

\[
\|a-(1/3,1/3,1/3)\|_1^2
\le \frac{8}{9}[-\log R_w].
\]

This theorem is exact inside the Gaussian packet model. Extending it to
arbitrary Navier--Stokes packets requires an inverse convolution theorem plus
control of the helical multiplier. Michael Christ's theorem that near
extremizers of sharp Young convolution on Euclidean space are close to Gaussian
extremizers identifies the appropriate rigidity mechanism, but it does not by
itself prove the Navier--Stokes packet theorem.

## Weighted Young inverse lemma

For a measurable multiplier \(0\le m\le m_*\), define

\[
T_m(f,g,h)=\iint m(p,q)f(p)g(q)\overline{h(p+q)}\,dp\,dq,
\qquad f,g,h\in L^{3/2}(\mathbb R^3).
\]

Let \(A=(\sqrt3/2)^3\) be the sharp scalar constant. If the three functions
are normalized and

\[
|T_m(f,g,h)|\ge(1-\delta)m_*A,
\]

then the triangle inequality and sharp Young inequality imply simultaneously:

1. the magnitudes form a \(\delta\)-near extremizer for sharp Young;
2. the trilinear mass gives mean multiplier deficit at most \(O(\delta)\);
3. the phase of \(m f g\overline h\) has \(L^2\)-oscillation at most
   \(O(\delta)\) with respect to that trilinear mass.

Michael Christ's inverse theorem then makes the magnitudes close in
\(L^{3/2}\) to an affine Gaussian extremizing triple. Applied to the
dimensionless helical transfer multiplier, local cusp stability further forces
that Gaussian cloud into an anisotropic tube around the maximizing triad
manifold.
