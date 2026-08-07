# Physical `H1`/swirl daughter forcing on mild-aspect grains

The intrinsic curvature decomposition gives a five-dimensional hook sector, but a transfer theorem must use the **physical Euclidean helicity planes**.  No aspect-independent comparison is available.  This note closes a controlled mild-aspect branch without charging anisotropy itself.

Let `B=B^H` be a hook tensor and let `L` be the physical grain factor.  For fixed normalized grain coordinate `c`, the physical symmetric strain-gradient slice is

\[
\boxed{C_c=\operatorname{Sym}(L B_cL^{-1}).}
\]

Define the polarization-only transfer-facing observable

\[
\boxed{
Q_{pol}(C)=\sum_c\left(
\|D_1(C_c)-D_2(C_c)\|_F^2+
\|D_3(C_c)\|_F^2
\right).
}
\]

The scalar shape term `D_Pi` is deliberately absent.

At `L=I`, write the swirl matrix as

\[
M=\begin{pmatrix}a&b&x\\b&d&y\\x&y&-a-d\end{pmatrix}.
\]

With `C=cos(phi)^2=1/(4r_*^2)`, direct expansion of `Q_pol-(1/10)||B||^2` is block diagonal in `(a,d),b,x,y`.  Arb positivity of the scalar blocks and the `(a,d)` Sylvester determinant throughout the certified `r_*` bracket gives

\[
\boxed{Q_{pol}(I)\ge\frac1{10}\|B\|_F^2.}
\]

For general `L`, remove its physical orthogonal polar factor and a common scalar.  If `kappa=cond(L)`, the normalized positive factor has spectrum in

\[
[\kappa^{-1/2},\kappa^{1/2}].
\]

Consequently

\[
\|\operatorname{Sym}(LB_cL^{-1})-\operatorname{Sym}B_c\|_F
\le(\kappa-1)\|B_c\|_F.
\]

For every symmetric trace-free physical strain slice,

\[
\|D_1-D_2\|_F^2+\|D_3\|_F^2\le5\|S\|_F^2,
\]

so the triangle inequality in the observable Hilbert space gives

\[
\sqrt{Q_{pol}(L)}
\ge
\left(\frac1{\sqrt{10}}-\sqrt5(\kappa-1)\right)\|B\|_F.
\]

At

\[
\boxed{\kappa\le21/20}
\]

interval arithmetic verifies

\[
\boxed{Q_{pol}(L)\ge\frac1{25}\|B\|_F^2.}
\]

Finally, any real symmetric trace-free `2x2` generator satisfies

\[
D^2=\frac12\|D\|_F^2I.
\]

Therefore its action on **every** unit complex helicity spinor has the same norm.  The combined relative-parent/child `H1` sideband forcing obeys

\[
\boxed{
\mathbb E\|F_{H1}^{rel}\|^2
\ge\frac1{50}\|B^H\|_F^2.
}
\]

This theorem does not price high aspect.  If `cond(L)>21/20`, the grain stays in the affine fresh/reuse/ancestry branch.  That separation is essential because common affine anisotropy remains an exact Young symmetry.
