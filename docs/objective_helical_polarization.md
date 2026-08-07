# Objective helical polarization and strain-area phase

For a Kelvin carrier `k(t)`, let `E(t)` be an orthonormal frame of the transverse
plane `k(t)^perp`.  The projected amplitude equation has real-frame generator

\[
-E^T(\nabla u)E-E^T\dot E.
\]

The frame has one free spin about `k`.  Choose it objectively so that

\[
E^T\dot E=-\operatorname{skew}(E^T\nabla u\,E).
\]

Then the transverse generator is exactly

\[
\boxed{
-\operatorname{sym}(E^T\nabla u\,E)-\nu|k|^2I.
}
\]

Rigid rotation has disappeared as a gauge.  Write the symmetric part as
`sigma I+D`, with

\[
D=\begin{pmatrix}\delta&\beta\\\beta&-\delta\end{pmatrix}.
\]

In the circular/helical basis,

\[
D_{hel}=\begin{pmatrix}0&\delta-i\beta\\\delta+i\beta&0\end{pmatrix}.
\]

Thus first-order nonconformal strain is a **helicity conversion** operator; it
is not an unavoidable diagonal Berry-phase drift.

A local geometric phase appears when strain orientations at different times do
not commute.  For two trace-free strains,

\[
\boxed{
[D_1,D_2]
=2(\delta_1\beta_2-\beta_1\delta_2)
\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
}
\]

The coefficient is the signed area spanned by the two anisotropy vectors
`(delta,beta)`.  It is an antisymmetric real rotation, hence an opposite phase
for the two circular polarizations.  This is the local strain analogue of a
geometric phase: **first order mixes helicity; second-order noncommutativity
rotates phase**.

Under the repository five-percent coherence hypothesis

\[
\|D(t)-D_0\|_{op}\le\varepsilon d,
\]

one has

\[
\|[D(t_1),D(t_2)]\|_{op}
\le (4\varepsilon+2\varepsilon^2)d^2.
\]

Hence the second Magnus generator obeys

\[
\boxed{
\|\Omega_2\|_{op}
\le(\varepsilon+\varepsilon^2/2)(dT)^2.
}
\]

At `epsilon=1/20`, `dT=1/30` this is below `10^-4`.  This is explicitly only a
second-Magnus estimate, not a claimed full time-ordered-exponential theorem.
If strain orientation fails the coherence hypothesis, the existing objective-
strain identity routes that variation into self-stretching/vorticity, pressure
Hessian, or viscous strain diffusion.
