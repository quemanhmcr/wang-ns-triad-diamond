# Gaussian packet inverse module

For L2-normalized isotropic Fourier Gaussian envelopes

\[
f_j(k)=(2\pi\sigma_j^2)^{-3/4}\exp\!\left(-\frac{|k-\kappa_j|^2}{4\sigma_j^2}\right)e^{-ix_j\cdot k},
\]

the scalar trilinear overlap is explicit. For equal widths \(\sigma\), relative to the aligned resonant case,

\[
R_{\rm scalar}
=
\exp\!\left[-\frac{|\kappa_1+\kappa_2-\kappa_3|^2}{12\sigma^2}
-\frac{\sigma^2}{3}\sum_{i<j}|x_i-x_j|^2\right].
\]

At fixed \(S=\sigma_1^2+\sigma_2^2+\sigma_3^2\), the width factor relative to its maximum is

\[
R_{\rm width}=(27a_1a_2a_3)^{3/4},\qquad a_j=\sigma_j^2/S.
\]

Thus near equality quantitatively forces frequency resonance, spatial co-location, and balanced widths. Combining this exact scalar module with the smooth helical coefficient on narrow nondegenerate caps gives a model inverse theorem for coherent Navier--Stokes packet transfer.

## Cusp stability at equal parent scales

Normalize the child magnitude to one and write the parent magnitudes as

\[
x=r+d,\qquad y=r-d.
\]

For the heterochiral maximizing sign pattern, near \(r=r_*\), the use of
\(\log(1/\max(x,y))\) makes the efficiency nonsmooth at \(d=0\). The one-sided
normalized slopes are

\[
a_+=2+\frac{1}{r_*\log(1/r_*)},\qquad
 a_-=\frac{1}{r_*\log(1/r_*)}-2,
\]

both positive. Numerically, \(a_+\approx5.322\) and \(a_-\approx1.322\).
Thus the local deficit has the anisotropic form

\[
1-\frac{J}{J_*}\gtrsim c|d|+c'(r-r_*)^2.
\]

The normal direction measuring parent-scale imbalance is therefore linear,
whereas the common-scale direction is quadratic. An isotropic Gaussian cap of
width \(\sigma\) must lose \(\Theta(\sigma)\), exactly as observed by the Actions
experiment. Near-extremal packet clouds must consequently be much thinner in
the scale-imbalance direction than in tangent directions; this is the proposed
"triad grain" structure.
