# High affine aspect is fresh radius or sticky ancestry

Affine anisotropy is an exact frozen Young symmetry, so aspect ratio itself must never be charged as a Bellman deficit. The correct question is dynamical: can a high-aspect grain appear from nowhere over one efficient packet lifetime?

The answer is no. A transition strip still has a quantitative H1 transfer observable, while stronger aspect must be inherited unless the grain is fresh.

## 1. Extend the physical H1 bridge to the full transition strip

The hook/swirl polarization estimate already gives

\[
\sqrt{Q_{pol}(L)}
\ge
\left[
\frac1{\sqrt{10}}-\sqrt5(\kappa-1)
\right]\|B^H\|,
\qquad \kappa=\operatorname{cond}L.
\]

Choose

\[
\boxed{
\kappa_{ext}=\frac{567}{500}
=\frac{27}{25}\frac{21}{20}.
}
\]

Arb certifies at this endpoint

\[
\boxed{
Q_{pol}\ge\frac1{4000}\|B^H\|^2.
}
\]

Because the three physical Young roles carry at least `Q_pol/4`,

\[
\boxed{
\sum_i\|F_i^{H1}\|^2
\ge\frac1{16000}\|B^H\|^2.
}
\]

Repeating the already-certified physical H1 Duhamel bookkeeping with the same low-strain action `K<=1/30` gives the conservative clean alternative

\[
\boxed{
\operatorname{Def}
\ \text{or}\ R_{pair}
\ge\frac{I_1^2}{28\,800\,000},
}
\]

outside dephasing/source, nonlinear-feedback and large-daughter exits.

Together with `I_B<=sqrt(6)I_3+I_1`, the common full-curvature transition-strip constant is

\[
\boxed{
\operatorname{Def}
\ \text{or}\ R_{pair}
\ge\frac{I_B^2}{115\,200\,000}.
}
\]

The constant is deliberately conservative; the point is positivity on the whole strip, not optimization.

## 2. Condition number can only be created by symmetric strain

The physical covariance satisfies

\[
\dot\Sigma=A\Sigma+\Sigma A^T+\nu I.
\]

Let `lambda_max,lambda_min` be its extremal eigenvalues and `kappa=sqrt(lambda_max/lambda_min)=cond(L)`. At simple eigenvalues,

\[
\frac d{dt}\log\kappa
=
 e_{max}^TSe_{max}-e_{min}^TSe_{min}
+rac\nu2\left(\lambda_{max}^{-1}-\lambda_{min}^{-1}\right).
\]

The viscous contribution is nonpositive. Therefore

\[
\boxed{
\frac d{dt}\log\kappa\le2\|S\|_{op}.
}
\]

This remains valid in the Dini-derivative sense at eigenvalue crossings. Hence

\[
\boxed{
\kappa(t_1)
\le
\kappa(t_0)
\exp\left(2\int_{t_0}^{t_1}\|S\|dt\right).
}
\]

On the existing efficient low-strain branch

\[
\int\|S\|dt\le\frac1{30},
\]

Arb gives

\[
\exp(1/15)<27/25.
\]

Consequently

\[
\boxed{
\kappa(t_1)>567/500
\Longrightarrow
\kappa(t_0)>21/20,
}
\]

provided the same covariance ancestry is being transported. Strong anisotropy is therefore **sticky**: it was already present in the predecessor. If no such predecessor exists, the grain is fresh/relinked.

## 3. Fresh high aspect pays physical radius, not an aspect deficit

The shell theorem gives every physical standard axis

\[
\ell_i>\frac{2}{3N}.
\]

If `kappa=ell_max/ell_min` and `s=Nr_g`, then

\[
r_g^3=\ell_1\ell_2\ell_3
\ge \ell_{min}^2\ell_{max}
=\kappa\ell_{min}^3,
\]

so

\[
\boxed{
s>\frac23\kappa^{1/3}.
}
\]

The affine critical-grain theorem gives

\[
N\int_{E_2}|u|^2\ge\frac3{10}s,
\]

therefore a fresh high-aspect occurrence satisfies

\[
\boxed{
N\int_{E_2}|u|^2
>\frac15\kappa^{1/3}.
}
\]

This is not a Young deficit. It is simply the physical fact that shell localization forbids squeezing below `N^-1`; large aspect can then only be made by elongation, which increases the geometric volume radius and physical energy occupancy.

## 4. Extended source branch

At `kappa<=567/500`, the H1 interaction forcing remains large enough that dephasing obeys the clean threshold

\[
J_1\ge I_1/(132T).
\]

The hook source calculation at the enlarged aspect cap has Arb-certified coefficients

\[
J_1
\le
\frac{11}{5}\int\|S_{src}\|dt
+60\int\|A\|\|B\|dt.
\]

Thus

\[
\boxed{
\int\|S_{src}\|dt\ge I_1/(600T)
}
\]

or

\[
\boxed{
\int\|A\|\|B\|dt\ge I_1/(16000T).
}
\]

On the H1-dominant full-curvature branch, base velocity-gradient/frame action below `1/32000` excludes the second alternative, so one of pressure-third, differentiated SGS, or viscous-fourth source channels has

\[
\boxed{
\int\|S_*\|dt\ge I_1/(1800T).
}
\]

This can be fed into the filtered SGS/viscous collision modules with a weaker but still positive source constant.

## 5. No-gap aspect architecture

One efficient packet lifetime now has no anonymous aspect escape:

\[
\boxed{
\begin{array}{ll}
\kappa\le567/500 &: \text{quantitative H1/curvature no-escape},\\
\kappa>567/500 &: \text{inherited high-aspect ancestry or fresh/relinked grain}.
\end{array}}
\]

Aspect itself is never charged. What is charged is either physical transfer deformation on the bounded transition strip or physical radius/freshness when a new elongated grain appears.
