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


## One-shot transfer consequence for efficient blocks

A full iterative profile decomposition inside one low-cost block is unnecessary.
If Christ supplies Gaussian profiles `F,G,H` with each `L^{3/2}` distance at
most `epsilon_G`, weighted Young and trilinearity give

\[
\frac{|T_m(f,g,h)-T_m(F,G,H)|}{m_*A_3}
\le 3\epsilon_G+3\epsilon_G^2+\epsilon_G^3.
\]

For a frequency role supported in `Omega_N`, `|Omega_N|<=C_Omega N^3`, finite
energy already implies Fourier `L^{3/2}` and the Gaussian profile carries

\[
N\|G\|_2^2\ge C_\Omega^{-1/3}(1-\epsilon_G)^2.
\]

Thus the external inverse Young theorem, plus the repository's exact weighted
multiplier/phase stability, gives one transfer-preserving critical-mass profile
on every sufficiently efficient frozen block.  The remaining issue is its
spacetime packet realization, not an infinite within-block extraction.


For the certified radial log shell `|log(|xi|/N)|<=2/25`, Action
`31165654509` makes the preceding mass bridge quantitative.  At
`epsilon_G=1/100`,

\[
|\Omega_N|/N^3
\le \frac{4\pi}{3}(e^{6/25}-e^{-6/25})
=2.0299769094\ldots
\]

and Arb certifies

\[
\boxed{N\|G\|_2^2>3/4.}
\]

The non-explicit part remains only Christ's implication from sufficiently small
Young deficit to one-percent Gaussian proximity.

## Affine covariance is a symmetry, not a defect

The Gaussian supplied by the inverse theorem need not be isotropic.  Write its
physical `L^2` covariance as `Sigma_x`.  The one-percent shell certificate gives

\[
\boxed{\lambda_{\min}(\Sigma_x)^{1/2}>2/(3N)}
\]

and transfers actual mass to the radius-two covariance ellipsoid `E_2`:

\[
\boxed{
\int_{E_2}|u|^2\ge \frac3{10}(\det\Sigma_x)^{1/6}.
}
\]

Thus with

\[
r_g=(\det\Sigma_x)^{1/6},
\qquad
\mathsf M_{aff}=r_g^{-1}\int_{E_2}|u|^2,
\]

one has the scale-critical affine grain bound

\[
\boxed{\mathsf M_{aff}\ge3/10.}
\]

For fresh affine grains of mass at least `eta` and overlap multiplicity `P`, the
physical energy inequality gives exactly

\[
\boxed{
\sum_j r_{g,j}\le P\|u(t)\|_2^2/\eta.
}
\]

This replaces the old heuristic that a long Gaussian automatically creates a
replication/Bellman cost.  At `p=3/2` arbitrary common affine deformation is an
exact Young symmetry.  Fresh anisotropic grains are charged by their physical
geometric radius; reused anisotropic grains must be controlled by spacetime
ancestry and the affine curvature/polarization dynamics.
