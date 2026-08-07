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
