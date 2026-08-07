# Intrinsic three-dimensional triad-plane dynamics

The first affine-grain theorem was written in a fixed planar reduction.  That
restriction is unnecessary for the **scalar triad geometry**.  A three-wave
relation always spans a two-dimensional plane, and a common linear carrier flow
maps that plane to another plane.  The correct object is its intrinsic Gram
matrix, not its orientation in the laboratory frame.

Let the two parent carriers be the columns of `K in R^{3x2}` and let

\[
G=K^TK.
\]

For any common carrier law

\[
\dot K=-BK,
\]

one has the exact identity

\[
\boxed{
\dot G=-K^T(B+B^T)K.
}
\]

The parent lengths and child length `|k_a+k_b|` are

\[
|k_a|^2=G_{11},\qquad
|k_b|^2=G_{22},\qquad
|k_a+k_b|^2=(1,1)G(1,1)^T.
\]

Thus every side-length quantity in the helical coupling magnitude and the
single-edge progress multiplier is intrinsic to `G`.  Extrinsic tilt of the
whole plane is a gauge direction for this scalar geometry.

Choose any orthonormal frame `E` of the instantaneous triad plane.  The
trace-free intrinsic shape driver is

\[
D=\left(E^T\operatorname{sym}(B)E\right)^0.
\]

At the symmetric extremal triad, orthogonal invariance reduces the calculation
to the previous two-dimensional formula, hence

\[
\boxed{
\frac12\dot u^2+2\dot v^2
\ge
4\sin^4(\theta_*/2)\,\|D\|_F^2
>
\frac{43}{100}\|D\|_F^2.
}
\]

This statement is valid for an arbitrarily oriented and tilting plane.  The
fixed-plane assumption was only a coordinate convenience for the scalar
multiplier dynamics.

For the common Gaussian carrier law

\[
B=A^T+2\nu P^{-1},
\]

the same statement includes anisotropic viscous spectral contraction.  If
`P=pI`, the viscous term is scalar and disappears from `D` in every plane.

If `n` is the unit normal to the carrier plane, differentiating `n^T K=0` gives

\[
\boxed{
\dot n=(I-nn^T)B^Tn.
}
\]

This tilt may be large without changing scalar side lengths when it is a common
isometric motion.  What remains genuinely three-dimensional is the transport of
helical polarization/phase and the spatial packet frame.  A future coherence
theorem should therefore be formulated in the evolving intrinsic plane rather
than requiring a fixed plane in physical coordinates.
