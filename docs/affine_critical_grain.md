# Affine critical grains: use physical energy scaling, not a false anisotropy cost

The inverse Young profile is affine Gaussian.  Forcing it into isotropic
`N^{-1}` cells is not natural because affine anisotropy is an exact frozen Young
symmetry.  The correct local-energy currency is instead determined directly by
Navier--Stokes scaling.

For physical covariance `Sigma_x`, define the geometric-mean radius

\[
\boxed{r_g=(\det\Sigma_x)^{1/6}.}
\]

For an affine ellipsoid `E` define

\[
\boxed{
\mathsf M_{aff}(E)=\frac1{r_g}\int_E|u|^2.
}
\]

Under the isotropic Navier--Stokes scaling `u_lambda(x,t)=lambda u(lambda x,
lambda^2t)`, both the energy in the scaled ellipsoid and `r_g` scale by
`lambda^{-1}`.  Hence `M_aff` is exactly scale critical.

Using the clean constants certified by `affine_shell_aspect` Action `31179827015`, on the radius-two covariance
ellipsoid of the one-percent Gaussian profile,

\[
\int_{E_2}|u|^2\ge\frac3{10}r_g,
\]

so

\[
\boxed{\mathsf M_{aff}(E_2)\ge\frac3{10}.}
\]

This is uniform in carrier frequency and aspect ratio.

If a family of fresh affine grains has `M_aff(E_j)>=eta` and physical overlap
multiplicity at most `P`, then

\[
\eta r_{g,j}\le\int_{E_j}|u|^2.
\]

Summing and using the physical energy budget gives

\[
\boxed{
\sum_j r_{g,j}\le\frac{P\|u(t)\|_2^2}{\eta}.
}
\]

This is the affine version of the sticky-cascade packet budget.  At isotropic
natural scale `r_g~N^{-1}` it reduces to the earlier `1/N` counting law.

The shell theorem also gives every principal physical standard axis
`l_i>2/(3N)`.  If

\[
s=Nr_g,\qquad A=Nl_{max},
\]

then `r_g^3=l_1l_2l_3` yields

\[
\boxed{A\le\frac94s^3.}
\]

Thus there is a physically correct anisotropy dichotomy:

- if `Nr_g=O(1)`, the aspect ratio is uniformly bounded;
- if the grain is extremely elongated, then `r_g` is correspondingly larger
  and each **fresh** occurrence consumes more of the conserved energy budget.

No Young or Bellman defect is assigned merely because the Gaussian is
anisotropic.  Reused elongated grains must instead be handled by the ancestry
and dynamic curvature/polarization ledgers.
