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
is not an unavoidable diagonal Berry-phase drift.  This local conversion must
also not be charged automatically as a transfer deficit: at equal-parent
geometry the signed nonlinear parent tensor is symplectic, so a common
`SL(2)` deformation of both parent spinors preserves the unnormalized parent
wedge exactly.  The transfer-distinguishable variables are relative parent
polarization, the child polarization factor, capacity normalization, and scalar
shape; see `docs/extremal_helicity_symplectic.md`.

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
rotates phase**.  A common time-ordered determinant-one map on both parents is
again symplectically neutral in the unnormalized parent wedge; only relative
incidence transport can turn this individual-spinor phase into a parent-sector
transfer obstruction.

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
\n\n## Transfer-facing time ordering is now exact\n\nThe second-Magnus estimate above remains useful for the phase of an individual
spinor.  It is no longer the relevant unresolved theorem for the extremal parent
transfer.  The signed parent tensor uses `U^T J V`, and the exact symplectic
identity gives

\[
\frac d{dt}(U^TJV)=U^TJ(D_1-D_2)V
\]

for arbitrary time-dependent noncommuting `D_i`.  Thus common parent time
ordering cancels exactly in the physical numerator; see
`docs/relative_polarization_transport.md`.\n