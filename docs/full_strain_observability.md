# Full 3D strain observability from triad shape and helicity conversion

Let `S` be a symmetric trace-free `3x3` strain tensor at the symmetric optimal
triad.  Denote its triad plane by `Pi` and the two parent directions by `k_1`
and `k_2`.  For any two-plane `V`, write `(S|_V)^0` for the trace-free part of
the restriction.

The scalar triad-shape dynamics sees

\[
D_\Pi=(S|_\Pi)^0,
\]

whereas the helical polarization of parent `i` sees the transverse strain

\[
D_i=(S|_{k_i^\perp})^0.
\]

At the extremal half-angle `phi`, `cos(phi)=1/(2r_*)`.  In coordinates

\[
S=\begin{pmatrix}
a&b&x\\ b&d&y\\ x&y&-a-d
\end{pmatrix},
\]

the quadratic form

\[
Q(S)=\|D_\Pi\|_F^2+\|D_1\|_F^2+\|D_2\|_F^2
\]

is block diagonal between `(a,d)`, `b`, `x`, and `y`.  Arb interval arithmetic
on the certified `r_*` bracket proves

\[
\boxed{
Q(S)\ge\frac{13}{20}\|S\|_F^2.
}
\]

The true smallest generalized eigenmode is an off-plane shear; its ratio is
about `0.6602495`, so `13/20` is conservative.

The scalar-shape theorem already gives

\[
\mathcal H_{speed}:=\frac12\dot u^2+2\dot v^2
\ge\frac{43}{100}\|D_\Pi\|_F^2.
\]

For the helical circular basis of parent `i`, the opposite-helicity mixing
coefficient `zeta_i` obeys exactly

\[
|\zeta_i|^2=\frac12\|D_i\|_F^2.
\]

Therefore

\[
\boxed{
\mathcal H_{speed}
+\frac{43}{50}\left(|\zeta_1|^2+|\zeta_2|^2\right)
\ge
\frac{559}{2000}\|S\|_F^2.
}
\]

This is a physical tomography theorem for the local strain: a large
incompressible symmetric strain cannot simultaneously hide from the triad's
side-length geometry and from the helicity conversion of both parents.  Rigid
rotation is absent because it belongs to the antisymmetric velocity gradient,
while isotropic dilation is absent because `tr S=0`.
