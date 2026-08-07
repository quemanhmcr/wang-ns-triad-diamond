# Affine ellipsoidal moving windows and curvature balance

Let an affine Gaussian grain have center `X(t)` and physical frame `L(t)`.  Use
an ellipsoidal window

\[
\chi_{L,M}(x,t)=\chi_0\!\left(\frac{L(t)^{-1}(x-X(t))}{M}\right),
\]

and transport the affine frame by the resolved flow jet

\[
\dot X=U(X),\qquad \dot L=A(X)L,\qquad A=\nabla U.
\]

With `z=L^{-1}(x-X)`, direct differentiation gives the exact identity

\[
\boxed{
(\partial_t+U\cdot\nabla)z
=L^{-1}\big(U(X+Lz)-U(X)-A(X)Lz\big).
}
\]

Thus affine translation and strain cancel exactly.  Only the non-affine Taylor
remainder moves the window relative to the packet.

Define the affine-normalized curvature over the window by

\[
\kappa_{aff}
=\sup\|L^{-1}(\nabla^2U)[L,L]\|.
\]

If `|grad chi_0|<=C_chi` and the transition of `chi_0` lies in `|z/M|<=R_chi`,
Taylor's integral remainder gives

\[
\boxed{
|(\partial_t+U\cdot\nabla)\chi_{L,M}|
\le \frac{C_\chi R_\chi^2}{2}\,\kappa_{aff}M.
}
\]

The opposite localization effect retains the `1/M` scale.  Since the shell
certificate gives every physical standard axis of the Gaussian profile

\[
\ell_{min}>\frac{2}{3N},
\]

one has

\[
\boxed{
N^{-1}\|\nabla_x\chi_{L,M}\|_\infty
\le\frac{3C_\chi}{2M}.
}
\]

For a physical coarse-graining kernel `G_N(y)=N^3G(Ny)`, this normalized
gradient gives an exact commutator estimate rather than only dimensional
scaling.  Indeed

\[
[\chi,G_N*]f(x)
=\int G_N(y)(\chi(x)-\chi(x-y))f(x-y)\,dy,
\]

so the mean-value theorem and Young inequality yield

\[
\boxed{
\|[\chi,G_N*]f\|_2
\le \frac{m_1(G)}{N}\|\nabla\chi\|_\infty\|f\|_2,
\qquad
m_1(G)=\int |y||G(y)|\,dy.
}
\]

Consequently the affine shell bound gives

\[
\boxed{
\|[\chi_{L,M},G_N*]f\|_2
\le \frac{3m_1(G)C_\chi}{2M}\|f\|_2.
}
\]

Thus the filter commutator and material-window remainder have the affine form

\[
\boxed{
E_{aff}(M)\le\frac{a}{M}+b\kappa_{aff}M,
}
\]

with optimizer

\[
\boxed{
M_*=\sqrt{\frac{a}{b\kappa_{aff}}},
\qquad E_{aff,*}=2\sqrt{ab\kappa_{aff}}.
}
\]

This is the ellipsoidal analogue of the earlier isotropic curvature balance and
requires no aspect-ratio penalty.  What remains is to derive the exact `a,b`
coefficients from the smooth-SGS packet/window construction, pressure work and
partition overlap rather than merely from geometric scaling.
