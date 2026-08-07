# Objective-strain source collision

The strain-coherence theorem says that a low-cost extremal grain can avoid a
quadratic multiplier loss only if its co-rotating non-conformal strain changes
by a definite fraction during one packet lifetime.  Navier--Stokes determines
what can cause that change.  This note converts the failure of coherence into
three physical source channels.

## 1. The exact source identity

For `A=grad u=S+Omega`, the objective strain derivative in the frame rotating
with the local vorticity is

\[
\mathring S
=-S^2-\Omega^2+[S,\Omega]
-\nabla^2p+\nu\Delta S.
\]

Write

\[
Q(A)=-S^2-\Omega^2+[S,\Omega].
\]

If 5% coherence fails, the preceding theorem gives

\[
\int_0^T\|\mathring S\|dt\ge\frac d{20}.
\]

Take `d=sigma N^2` and `T=cN^-2`.  After dividing by `TN^4`, at least one of
`Q`, `Hess p`, or `nu Delta S` has average normalized size

\[
\boxed{
\rho_{src}\ge\frac{\sigma}{60c}.
}
\]

Thus temporal dephasing is not a fourth free escape variable: it must be paid by
a concrete NS source.

## 2. Self-stretching/vorticity source forces critical mass

The elementary operator bound

\[
\|Q(A)\|\le4\|A\|^2
\]

is sufficient.  For a packet frequency cutoff `Lambda N`, unitary-Fourier
Cauchy--Schwarz gives

\[
N^{-2}\|\nabla u\|_\infty\le C_1\sqrt\mu,
\qquad
\mu=N\|u\|_2^2,
\]

with the scalar derivative constant

\[
C_1=(2\pi)^{-3/2}\sqrt{4\pi/5}\,\Lambda^{5/2}
\]

(up to the fixed vector/operator norm factor used by the packet frame).  Hence

\[
N^{-4}\|Q(A)\|\le4C_1^2\mu.
\]

A source level `rho` therefore forces

\[
\boxed{
\mu\ge\frac{\rho}{4C_1^2}.
}
\]

## 3. Viscous strain variation also forces critical mass

Likewise

\[
N^{-4}\|\nabla^3u\|_\infty\le C_3\sqrt\mu,
\]

where

\[
C_3=(2\pi)^{-3/2}\sqrt{4\pi/9}\,\Lambda^{9/2}.
\]

Thus a viscous source level `rho` gives

\[
\boxed{
\mu\ge\left(\frac{\rho}{\nu C_3}\right)^2.
}
\]

So rapid viscous reshaping of the grain cannot occur without scale-critical
energy at the relevant frequency.

## 4. Pressure Hessian has a stronger far-field locality exponent

The pressure kernel itself is order `|z|^-3`.  Its Hessian is order
`|z|^-5`.  A direct product-rule estimate for

\[
K_{ij}(z)=\frac1{4\pi}\left(3z_i z_j|z|^{-5}-\delta_{ij}|z|^{-3}\right)
\]

gives the clean bound

\[
\boxed{
\sum_{ij}\|\nabla^2K_{ij}(z)\|\,|v_i v_j|
\le150|z|^{-5}|v|^2.
}
\]

A critical-mass packet in shell `2^n/N` therefore contributes normalized
pressure-Hessian size `O(mu 2^-5n)`.  Three-dimensional packing grows like
`2^{3n}`, leaving

\[
\boxed{5-3=2}
\]

powers:

\[
\sum_n2^{-5n}2^{3n}=\sum_n2^{-2n}<\infty.
\]

This is even more local than the pressure-work branch, where constant-mode
subtraction produced only the spare exponent `4-3=1`.

Under absence of a fresh far packet with critical mass above `mu_*`,

\[
N^{-4}\|\nabla^2p_{far}\|\le C_{far}\mu_*.
\]

If a packet construction also supplies a local coefficient
`H_near<=C_near mu_*`, then a pressure-Hessian source `rho` forces

\[
\boxed{
\mu_{fresh}
\ge
\min\{\rho/(2C_{near}),\rho/(2C_{far})\}.
}
\]

The far coefficient is already explicit from kernel decay and packing; the near
coefficient remains part of the concrete spatial packet realization.

## 5. New low-cost dynamical alternative

A near-extremal packet lineage now has the following dynamical alternatives:

\[
\boxed{
\text{coherent non-conformal strain}
\to\text{multiplier cost},
}
\]

or

\[
\boxed{
\text{strain dephasing}
\to
Q(A)\ \text{or}\ \nabla^2p\ \text{or}\ \nu\Delta S,
}
\]

and in the band-limited packet model every source either forces critical mass
locally or has a summable far-pressure tail.  The continuum task is now to
register that mass as fresh/reused and to verify the local pressure-Hessian
coefficient in the same adaptive packet frame.
