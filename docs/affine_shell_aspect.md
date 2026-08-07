# Shell concentration, uncertainty and the affine-aspect escape

The one-shot inverse Young theorem supplies an affine Gaussian profile, so an
isotropic `N^{-1}` packet radius must not be assumed.  The certified frequency
shell nevertheless removes one half of the affine freedom: physical axes cannot
be squeezed arbitrarily below the natural scale.

Write a normalized Fourier Gaussian as

\[
F(\xi)=C\exp\!\left[-\frac12(\xi-\kappa)^TA(\xi-\kappa)\right],
\qquad \|F\|_{3/2}=1.
\]

The probability measure proportional to `|F|^(3/2)` has covariance

\[
\Gamma_\xi=\frac23A^{-1},
\]

while in the unitary Fourier convention the physical `L^2` density of
`\check F` has covariance

\[
\Sigma_x=\frac12A.
\]

Hence the matrix uncertainty relation is exact:

\[
\boxed{\Sigma_x\Gamma_\xi=\frac13I.}
\]

If `f` is supported in the certified shell

\[
|\log(|\xi|/N)|\le\frac2{25}
\]

and `||f-F||_(3/2)<=1/100`, then outside the shell
`||F||_(3/2)<=1/100`.  Thus the `|F|^(3/2)` probability of the shell is at least
`999/1000`.  Every one-dimensional projection is consequently a Gaussian that
puts at least that mass in an interval of length `2 exp(2/25)N`.  Bounding an
interval probability by interval length times the Gaussian peak density gives

\[
\lambda_{\max}(\Gamma_\xi)^{1/2}
\le
\frac{2e^{2/25}}{(999/1000)\sqrt{2\pi}}N.
\]

Arb then certifies

\[
\boxed{\lambda_{\min}(\Sigma_x)^{1/2}>\frac{2}{3N}.}
\]

So a near-extremal affine profile can be long, but no physical principal axis
can be arbitrarily shorter than `N^{-1}`.

There is also a localized mass consequence that applies to the actual role
`f`, not only to the Gaussian profile.  Let

\[
E_2=\{x:(x-X)^T\Sigma_x^{-1}(x-X)\le4\}.
\]

For a normalized Gaussian,

\[
\|F\|_2^2
=\frac{9\sqrt2}{16\sqrt\pi}(\det\Sigma_x)^{1/6}.
\]

The radius-two ellipsoid contains the Maxwell probability

\[
P_2=\operatorname{erf}(\sqrt2)-2\sqrt{2/\pi}e^{-2}.
\]

Hausdorff--Young in the unitary convention gives
`||\check f-\check F||_3<=1/100`; Holder on `E_2` therefore transfers the
Gaussian `L^2` mass to the actual packet.  Arb certifies the clean bound

\[
\boxed{
N\int_{E_2}|\check f|^2
\ge \frac{3}{10}N(\det\Sigma_x)^{1/6}.
}
\]

If `A=N lambda_max(Sigma_x)^(1/2)` is the longest dimensionless physical axis,
then the lower bounds on the other two axes imply

\[
\boxed{
N\int_{E_2}|\check f|^2>\frac15 A^{1/3}.
}
\]

This is an **ellipsoidal critical-mass** statement.  It is not yet a theorem
that converts the mass into a natural-cell fresh packet.

That distinction is necessary.  At `p=3/2`, for every invertible linear map
`S`,

\[
F_S(\xi)=|\det S|^{2/3}F(S\xi)
\]

preserves both the `L^{3/2}` norm and the normalized Young trilinear form.
Therefore arbitrary common affine anisotropy is an exact frozen transfer
symmetry.  Merely covering a long affine Gaussian by many isotropic cells does
not by itself create a Bellman deficit.  The dynamic obstruction must come from
the grain-normalized curvature/forcing or from a separate ancestry theorem for
ellipsoidal mass.
