# Divergence-free curvature as an `H3 + H1` irreducible multiplet

The previous sideband theorem treats envelope `H3` and polarization/vector `H1` excitations separately.  This note shows that, in the affine-normalized curvature tensor itself, those are the two irreducible non-affine sectors.

Let

\[
B_{abc}=B_{acb},\qquad B_{aac}=0,
\]

where the second identity is differentiated incompressibility in grain coordinates.  Put

\[
T=\operatorname{Sym}B,\qquad t_c=T_{aac}.
\]

Define

\[
\boxed{
B^E_{abc}=T_{abc}
-\frac12(\delta_{ab}t_c+\delta_{ac}t_b)
+\delta_{bc}t_a.
}
\]

Then `B^E` is divergence free and `Sym B^E=T`.  The remainder

\[
B^H=B-B^E
\]

is divergence free and has zero full symmetrization.  Moreover the two pieces are orthogonal in Frobenius inner product.

The symmetric trace satisfies the sharp three-dimensional inequality

\[
\boxed{\|t\|^2\le\frac53\|T\|^2.}
\]

One proof is to split `T` orthogonally into its trace-free `l=3` part and the minimum-norm pure-trace tensor

\[
T^{tr}_{abc}=\frac15(\delta_{ab}t_c+\delta_{ac}t_b+\delta_{bc}t_a),
\]

for which `||T^tr||^2=(3/5)||t||^2`.

Direct contraction gives

\[
\boxed{\|B^E\|^2=\|T\|^2+3\|t\|^2\le6\|T\|^2.}
\]

The hook part is exactly the quadratic-swirl representation

\[
\boxed{
B^H_{abc}=\varepsilon_{abd}M_{dc}+\varepsilon_{acd}M_{db},
\quad M=M^T,\quad\operatorname{tr}M=0,
}
\]

with inverse

\[
M_{ec}=\frac13\varepsilon_{abe}B^H_{abc}.
\]

Its norm is

\[
\boxed{\|B^H\|^2=6\|M\|^2.}
\]

The normalized grain-coordinate strain sideband

\[
C^H_{adc}=\frac12(B^H_{adc}+B^H_{dac})
\]

obeys

\[
\boxed{\|C^H\|^2=\frac14\|B^H\|^2.}
\]

Combining the orthogonal decomposition with the preceding bounds gives

\[
\boxed{
\|\operatorname{Sym}B\|^2+\|C^H\|^2
\ge\frac16\|B\|^2.
}
\]

Thus, in the intrinsic affine grain, curvature cannot hide from both the scalar `H3` envelope and the vector `H1` swirl channel.  Representation-theoretically the 15-dimensional differentiated-incompressible curvature space splits as

\[
15=(7\oplus3)_{\rm symmetric/envelope}\oplus5_{\rm swirl}.
\]

The `C^H` norm is intrinsic to the affine grain.  This theorem does **not** erase the existing caveat that comparison with physical Euclidean helical-generator curvature can deteriorate with aspect; that comparison remains an ancestry/polarization issue rather than a static affine cost.
