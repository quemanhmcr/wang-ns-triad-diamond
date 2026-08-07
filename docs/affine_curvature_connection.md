# Affine curvature connection: dephasing is material-acceleration Hessian

The `H_1/H_3` sideband ledger exposes spatial curvature of the transporter.  To
understand whether that curvature remains coherent through a packet lifetime,
its **common affine deformation must first be quotiented**.

Let a smooth velocity `V` define the packet center/frame by

\[
\dot X=V(X),\qquad \dot L=A L,\qquad A=\nabla V(X),
\]

and put

\[
H=\nabla^2V(X),\qquad
B=L^{-1}H[L,L].
\]

Let

\[
F=D_t^V V=(\partial_t+V\cdot\nabla)V
\]

be the actual material acceleration of the transporter.  The commutator
`[D_t,partial_j]=-A_{mj}partial_m` gives

\[
(D_tH)_{i,jk}
=(\nabla^2F)_{i,jk}
-A_{im}H_{m,jk}-A_{mj}H_{i,mk}-A_{mk}H_{i,jm}.
\]

Differentiating `B=L^-1 H[L,L]` then cancels both input-index stretching terms
and yields

\[
\boxed{
\dot B+2A_{aff}B
=L^{-1}(\nabla^2F)[L,L],
\qquad A_{aff}=L^{-1}AL.
}
\]

This is an affine **connection** on the curvature tensor.  A source-free
curvature field is not constant in raw coordinates; it is parallel transported
by the common affine jet.  Charging `dot B` directly would therefore repeat the
same mistake as charging raw Berry phase or raw parent `SL(2)` motion.

## Resolved Navier--Stokes source split

For a convolution-resolved transporter,

\[
D_t^V V=-\nabla P-\nabla\cdot R+\nu\Delta V.
\]

Hence

\[
\boxed{
\dot B+2A_{aff}B
=L^{-1}\left[
-\nabla^3P-\nabla^2\nabla\cdot R+\nu\nabla^2\Delta V
\right][L,L].
}
\]

Thus rapid dephasing of the affine curvature sideband must ultimately be paid by
one of three physical source channels:

1. pressure third derivatives;
2. differentiated SGS stress;
3. viscous fourth velocity derivatives.

No positive cost is claimed merely from this identity.  It identifies the
correct source if a future sideband-coherence theorem says that the parallel
curvature cannot drift too much without paying.

## Pressure becomes even more local at this derivative level

In three dimensions the pressure singular-integral kernel is homogeneous of
degree `-3`.  Three spatial derivatives are therefore homogeneous of degree
`-6`.  A packet at dyadic distance `2^n` contributes a `2^{-6n}` kernel factor,
while three-dimensional packing contributes at most `2^{3n}` packets.  The far
sum is consequently

\[
\boxed{
\sum_n2^{-3n}<\infty,
\qquad 6-3=3.
}
\]

This is stronger than the previously used pressure-Hessian locality
`5-3=2`.  A continuum source-collision theorem still needs the appropriate
near-field coefficient and differentiated-stress control; only the far
homogeneity/packing exponent is asserted here.
